#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdint.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
typedef struct Bmpheader {
	char identifier[2]; // 0x0000
	unsigned int filesize; // 0x0002
	unsigned short reserved; // 0x0006
	unsigned short reserved2;
	unsigned int bitmap_dataoffset; // 0x000A
	unsigned int bitmap_headersize; // 0x000E
	unsigned int width; // 0x0012
	unsigned int height; // 0x0016
	unsigned short planes; // 0x001A
	unsigned short bits_perpixel; // 0x001C
	unsigned int compression; // 0x001E
	unsigned int bitmap_datasize; // 0x0022
	unsigned int hresolution; // 0x0026
	unsigned int vresolution; // 0x002A
	unsigned int usedcolors; // 0x002E
	unsigned int importantcolors; // 0x0032
} Bitmap;

/*construct a structure of RGB*/
typedef struct Pixel{
	unsigned char*data;
}Pixel;

typedef struct Img{
	Bitmap header;
	unsigned char*palette;
	Pixel**data;
}Img;

typedef struct RLE_node{
	double val;
	struct RLE_node*next;
}RLE_node;

typedef struct Huffman_node{
	short val;
	int count;
	struct Huffman_node*left;
	struct Huffman_node*right;
	struct Huffman_node*next;
}Huffman_node;

typedef struct Huffman_queue{
	struct Huffman_node*front;//表頭,不是第0個資料 
	struct Huffman_node*rear;
}Huffman_queue;

void readheader(FILE* fp, Bitmap *x);
Pixel** malloc_2D_Pixel(int row, int col,unsigned short bytes_perpixel);
uint8_t**malloc_uint8_2D(int row,int col);
double**malloc_2D_double(int row_size,int col_size);
double***malloc_3D_double(int size1,int size2,int size3);
double****malloc_4D(int size1,int size2,int size3,int size4);
Img ImRead(char input_name[]);
void ImWrite(Img img,char filename[]);
void reverse_row(Pixel**arr,int rowSize,int colSize);
Bitmap make_basis_header();
unsigned char*make_basis_palette(unsigned short palette_size);
void cut_8x8(double**data,double matrix_8x8[8][8],int m,int n);
void DCT_8x8(double matrix_8x8[8][8],double****F,double**beta,double ****basic_vector,int m,int n,int u,int v);
void IDCT_8x8(double**F,double**matrix_8x8,double**beta,double****basic_vector,int r,int c);
double**ones(int row,int col);
void uint8(double x[8][8],uint8_t**bmp_x);
void bgr2ycbcr(Img*img,double**y_data,double**cb_data,double**cr_data);
void ycbcr2bgr(Img*img,double**y_data,double**cb_data,double**cr_data);
RLE_node***malloc_3D_RLE_node(int M,int N);

RLE_node*new_RLE_node(double next_val);
void add_RLE_node(RLE_node*last,RLE_node*newnode);
void RLE_8x8(double*F_zigzag,RLE_node*list);
void un_RLE_8x8(RLE_node*list,double*F_zigzag);
long long int RLE_code_len(RLE_node*RLE_code);

Huffman_queue Huffman_queue_init();
void Huffman_queue_push(Huffman_queue*queue,Huffman_node*node);
void Huffman_queue_pop(Huffman_queue*queue,Huffman_node*target);
Huffman_node*huffmandict(RLE_node*RLE_code,int len);
void print_tree(Huffman_node*node,int*len);
int main(int argc, char **argv){
	//有顏色差是因為quantize的問題,rgb轉ycbcr也會造成資料誤差(用四捨五入減緩) 
	
	int W=8;
	int H=8;
	//generate basis vector for 2D-DCT
	double****basic_vector=malloc_4D(8,8,8,8);
	int u;
	/*
	Because we want to put the basic vector picture in the directory we specified, 
	we need to change the workspace to the directory we specified
	*/
	for(u=0;u<H;u++){
		int v;
		for(v=0;v<W;v++){
			int r;
			//Make a 8x8 basic vector matrix x
			double x[8][8];
			for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
				int c;
				for(c=0;c<W;c++){ 
					basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
					x[r][c]=127.5+127.5*basic_vector[u][v][r][c];
					
				}
				
			}
			
		}
	}
	//read image
	Img img=ImRead("Kimberly.bmp");
	int M=img.header.height/H;//圖片高有幾個8x8 matrix
	int N=img.header.width/W;//圖片寬有幾個8x8 matrix
	
	//bgr to ycbcr
	double**y_data=malloc_2D_double(img.header.height,img.header.width);
	double**cb_data=malloc_2D_double(img.header.height,img.header.width);
	double**cr_data=malloc_2D_double(img.header.height,img.header.width);
	bgr2ycbcr(&(img),y_data,cb_data,cr_data);
	
	
	//2D-DCT
	double****y_F=malloc_4D(M,N,8,8);//frequency domain of the 8x8 matrices
	double****cb_F=malloc_4D(M,N,8,8);
	double****cr_F=malloc_4D(M,N,8,8);
	double**beta=ones(H,1);
	beta[0][0]=1/sqrt(2);
	int m;
	//Cut the matrix into multiple 8x8 matrices
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			
			//Get one of the 8x8 matrices
			double y_matrix_8x8[8][8];
			double cb_matrix_8x8[8][8];
			double cr_matrix_8x8[8][8];
			cut_8x8(y_data,y_matrix_8x8,m,n);
			cut_8x8(cb_data,cb_matrix_8x8,m,n);
			cut_8x8(cr_data,cr_matrix_8x8,m,n);
			//subtract 128
			int v;
			for(u=0;u<8;u++){
				for(v=0;v<8;v++){
					y_matrix_8x8[u][v]-=128;
					cb_matrix_8x8[u][v]-=128;
					cr_matrix_8x8[u][v]-=128;
				}
			}
			//Get the frequency domain of the 8x8 matrix with DCT
			for(u=0;u<H;u++){
				for(v=0;v<W;v++){
					DCT_8x8(y_matrix_8x8,y_F,beta,basic_vector,m,n,u,v);
					DCT_8x8(cb_matrix_8x8,cb_F,beta,basic_vector,m,n,u,v);
					DCT_8x8(cr_matrix_8x8,cr_F,beta,basic_vector,m,n,u,v);
					//if(y_F[m][n][u][v]>=0.5||y_F[m][n][u][v]<=-0.5)printf("%.2lf ",y_F[m][n][u][v]);
				}
			}
			//printf("\n\n");
		}
		//return 0;
	}
	
	//quantization
	double Q[8][8]={{16,11,10,16,24,40,51,61},
     {12,12,14,19,26,58,60,55},
     {14,13,16,24,40,57,69,56},
     {14,17,22,29,51,87,80,62},
     {18,22,37,56,68,109,103,77},
     {24,35,55,64,81,104,113,92},
     {49,64,78,87,103,121,120,101},
     {72,92,95,98,112,100,103,99}};
     
    double Q2[8][8]={{1,1,1,1,1,1,1,1},
     {1,1,1,1,1,1,1,1},
     {1,1,1,1,1,1,1,1},
     {1,1,1,1,1,1,1,1},
     {1,1,1,1,1,1,1,1},
     {1,1,1,1,1,1,1,1},
     {1,1,1,1,1,1,1,1},
     {1,1,1,1,1,1,1,1}};
    
    
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			int u;
			//quantize and round each 8x8 matrix
			for(u=0;u<8;u++){
				int v;
				for(v=0;v<8;v++){
					int y_tmp=(int)(y_F[m][n][u][v]/Q[u][v]+0.5);
					y_F[m][n][u][v]=(double)y_tmp;
					
					int cb_tmp=(int)(cb_F[m][n][u][v]/Q2[u][v]+0.5);
					cb_F[m][n][u][v]=(double)cb_tmp;
				
					
					int cr_tmp=(int)(cr_F[m][n][u][v]/Q2[u][v]+0.5);
					cr_F[m][n][u][v]=(double)cr_tmp;
					
				}
			}
		}
	} 
	
	/*double max,min;
	max=-9999999;
	min=99999999;*/
	//DC DPCM
	/*for(m=0;m<M;m++){
		int n;
		for(n=N-1;n>=1;n--){
			y_F[m][n][0][0]=y_F[m][n][0][0]-y_F[m][n-1][0][0];
			cb_F[m][n][0][0]=cb_F[m][n][0][0]-cb_F[m][n-1][0][0];
			cr_F[m][n][0][0]=cr_F[m][n][0][0]-cr_F[m][n-1][0][0];
		}
	}
	for(m=M-1;m>=1;m--){
		y_F[m][0][0][0]=y_F[m][0][0][0]-y_F[m-1][0][0][0];
		cb_F[m][0][0][0]=cb_F[m][0][0][0]-cb_F[m-1][0][0][0];
		cr_F[m][0][0][0]=cr_F[m][0][0][0]-cr_F[m-1][0][0][0];
	} */
	/*printf("max=%.2lf min=%.2lf\n",max,min);
	return 0;*/
	
	
	//Zigzag 
	/*int zz_matrix[8][8]={
	{1,2,6,7,15,16,28,29},
	{3,5,8,14,17,27,30,43},
	{4,9,13,18,26,31,42,44},
	{10,12,19,25,32,41,45,54},
	{11,20,24,33,40,46,53,55},
	{21,23,34,39,47,52,56,61},
	{22,35,38,48,51,57,60,62},
	{36,37,49,50,58,59,63,64}
	};
	double***y_F_zigzag=malloc_3D_double(M,N,8*8);
	double***cb_F_zigzag=malloc_3D_double(M,N,8*8);
	double***cr_F_zigzag=malloc_3D_double(M,N,8*8);
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			int r,c;
			for(r=0;r<8;r++){
				for(c=0;c<8;c++){
					y_F_zigzag[m][n][zz_matrix[r][c]-1]=y_F[m][n][r][c];
					cb_F_zigzag[m][n][zz_matrix[r][c]-1]=cb_F[m][n][r][c];
					cr_F_zigzag[m][n][zz_matrix[r][c]-1]=cr_F[m][n][r][c];
				}
			}
		}
	} */
	
	//RLE,node紀錄數列有幾個0以及下個不是0的數值
	/*RLE_node***y_RLE_lists=malloc_3D_RLE_node(M,N);
	RLE_node***cb_RLE_lists=malloc_3D_RLE_node(M,N);
	RLE_node***cr_RLE_lists=malloc_3D_RLE_node(M,N);
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			RLE_8x8(y_F_zigzag[m][n],y_RLE_lists[m][n]);
			RLE_8x8(cb_F_zigzag[m][n],cb_RLE_lists[m][n]);
			RLE_8x8(cr_F_zigzag[m][n],cr_RLE_lists[m][n]);	
		}
		
	}*/
	//RLE_code(simple)將所有list串在一起 
	/*RLE_node*RLE_code;
	RLE_node*yptr;
	RLE_node*cbptr;
	RLE_node*crptr;
	int len1,len2,len3;
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			if(m==0&&n==0){
				yptr=y_RLE_lists[0][0];
				cbptr=cb_RLE_lists[0][0];
				crptr=cr_RLE_lists[0][0];
			}
			else{
				while(yptr->next!=NULL){
					yptr=yptr->next;
				}
				yptr->next=y_RLE_lists[m][n];
				yptr=yptr->next;
				
				while(cbptr->next!=NULL){
					cbptr=cbptr->next;
				}
				cbptr->next=cb_RLE_lists[m][n];
				cbptr=cbptr->next;
				
				while(crptr->next!=NULL){
					crptr=crptr->next;
				}
				crptr->next=cr_RLE_lists[m][n];
				crptr=crptr->next;
			}
		}
	}
	yptr->next=cb_RLE_lists[0][0]->next;
	cbptr->next=cr_RLE_lists[0][0]->next;
	RLE_code=y_RLE_lists[0][0];*/
	
	//Huffman
	/*
	long long int len=RLE_code_len(RLE_code); 
	//printf("%ld",len);
	Huffman_node*dict=huffmandict(RLE_code,len);
	int max_branch=0;
	printf("%d %d\n",dict->left->count,dict->right->count);
	//print_tree(dict,&max_branch);
	return 0;
	*/
	
	//un Huffman
	/*for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			AC_un_Huffman();
			AC_un_Huffman();	
			AC_un_Huffman();			
		}
	}*/ 
	
	//DC_un_Huffman();
	
	
	
	
	//un RLE
	/*
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			un_RLE_8x8(y_RLE_lists[m][n],y_F_zigzag[m][n]);
			un_RLE_8x8(cb_RLE_lists[m][n],cb_F_zigzag[m][n]);
			un_RLE_8x8(cr_RLE_lists[m][n],cr_F_zigzag[m][n]);
		}
	}*/
	//printf("un RLE finish\n");
	
	
	
	//unzigzag
	/*for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			int r,c;
			for(r=0;r<8;r++){
				for(c=0;c<8;c++){
					y_F[m][n][r][c]=y_F_zigzag[m][n][zz_matrix[r][c]-1];
					cb_F[m][n][r][c]=cb_F_zigzag[m][n][zz_matrix[r][c]-1];
					cr_F[m][n][r][c]=cr_F_zigzag[m][n][zz_matrix[r][c]-1];
				}
			}
		}
	}*/
	
	//un DPCM
	/*for(m=1;m<M;m++){
		y_F[m][0][0][0]=y_F[m][0][0][0]+y_F[m-1][0][0][0];
		cb_F[m][0][0][0]=cb_F[m][0][0][0]+cb_F[m-1][0][0][0];
		cr_F[m][0][0][0]=cr_F[m][0][0][0]+cr_F[m-1][0][0][0];
	}
	for(m=0;m<M;m++){
		int n;
		for(n=1;n<N;n++){
			y_F[m][n][0][0]=y_F[m][n][0][0]+y_F[m][n-1][0][0];
			cb_F[m][n][0][0]=cb_F[m][n][0][0]+cb_F[m][n-1][0][0];
			cr_F[m][n][0][0]=cr_F[m][n][0][0]+cr_F[m][n-1][0][0];
		}
	}*/
	
	//unquantization
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			int u;
			//unquantize and round each 8x8 matrix
			for(u=0;u<8;u++){ 
				int v;
				for(v=0;v<8;v++){ 
					y_F[m][n][u][v]*=Q[u][v];
					cb_F[m][n][u][v]*=Q2[u][v];
					cr_F[m][n][u][v]*=Q2[u][v];
				}
			}
		}
	} 
	
	//2D-IDCT
	Pixel**res_data=malloc_2D_Pixel(img.header.height,img.header.width,3);
	for(m=0;m<M;m++){
		int n;
		for(n=0;n<N;n++){
			int r,c;
			double**y_tmp=ones(8,8);
			double**cb_tmp=ones(8,8);
			double**cr_tmp=ones(8,8);
			//Get the spatial domain of the 8x8 frequency matrix with IDCT
			for(r=0;r<8;r++){
				for(c=0;c<8;c++){
					IDCT_8x8(y_F[m][n],y_tmp,beta,basic_vector,r,c);
					IDCT_8x8(cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
					IDCT_8x8(cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
					y_tmp[r][c]+=128;
					cb_tmp[r][c]+=128;
					cr_tmp[r][c]+=128;
				}
			}
			
			
			for(r=0;r<8;r++){
				for(c=0;c<8;c++){
					y_data[8*m+r][8*n+c]=y_tmp[r][c];
					cb_data[8*m+r][8*n+c]=cb_tmp[r][c];
					cr_data[8*m+r][8*n+c]=cr_tmp[r][c];
				}
			}
			
			
		}
	}
	//ycbcr to bgr
	ycbcr2bgr(&img,y_data,cb_data,cr_data);
	return 0;
	//output result
	ImWrite(img,"output.bmp");
	//ImWrite(img,argv[2]);
	return 0;
}
void print_tree(Huffman_node*node,int*len){
	if(node==NULL){
		return;
	}
	else{
		printf("%d ",node->count);
		(*(len))++;
		print_tree(node->left,len);
		print_tree(node->right,len);
	}
}


Huffman_node*huffmandict(RLE_node*RLE_code,int len){
	//計算每種RLE_node的出現次數，並放入queue 
	unsigned char*used=(unsigned char*)calloc(len,sizeof(unsigned char));
	int ptr1=0;
	Huffman_queue queue=Huffman_queue_init();
	RLE_node*node1=RLE_code->next;
	int node_num=0;
	//printf("%d",len);
	while(ptr1<len){
		if(used[ptr1]==0){
			int count=0;
			int ptr2=ptr1;
			RLE_node*node2=node1;
			
			while(ptr2<len){
				if(node2->val==node1->val){
					used[ptr2]=1;
					count++;
				}
				ptr2++;
				node2=node2->next;
			}
			Huffman_node*node=(Huffman_node*)malloc(sizeof(Huffman_node));
			(node->count)=count;
			(node->val)=node1->val;
			(node->left)=(node->right)=NULL;
			(node->next)=NULL;
			//printf("%d ",count);
			Huffman_queue_push(&queue,node);
			node_num++;
		}
		ptr1++;
		node1=node1->next;
	}
	free(used);
	
	while(node_num>1){
		//printf("%x ",queue.front->next);
		//printf("%d ",node_num);
		Huffman_node*ptr_most=queue.front->next;//最小和第二小的位置
		Huffman_node*ptr=queue.front->next;
		while(ptr!=NULL){
			if((ptr_most->count)>(ptr->count)){
				ptr_most=ptr;
			}
			
			ptr=ptr->next;
		}
		ptr=queue.front->next;
		Huffman_node*ptr_second=queue.front->next;
		while(ptr!=NULL){
			if(ptr_second==ptr_most||((ptr->count)<(ptr_second->count)&&(ptr->count)>=(ptr_most->count)&&ptr!=ptr_most)){
				ptr_second=ptr;
			}
			ptr=ptr->next;
		}
		Huffman_node*node=(Huffman_node*)malloc(sizeof(Huffman_node));
		node->left=ptr_most;
		node->right=ptr_second;
		node->count=ptr_most->count+ptr_second->count;
		node->next=NULL;
		Huffman_queue_pop(&queue,ptr_most);
		Huffman_queue_pop(&queue,ptr_second);
		Huffman_queue_push(&queue,node);
		node_num--;
	}
	Huffman_node*head=queue.front->next;
	printf("%d %d\n",head->left->count,head->right->count);
	return queue.front->next;
}

Huffman_queue Huffman_queue_init(){
	Huffman_queue res;
	res.front=(Huffman_node*)malloc(sizeof(Huffman_node));
	res.rear=res.front;
	return res;
}

void Huffman_queue_push(Huffman_queue*queue,Huffman_node*node){
	if(queue->rear==NULL){
		printf("error ");
	}
	queue->rear->next=node;
	queue->rear=node;
	/*Huffman_node*ptr=queue->front;
	while(ptr->next!=NULL){
		ptr=ptr->next;
	}
	ptr->next=node;*/
}

void Huffman_queue_pop(Huffman_queue*queue,Huffman_node*target){
	Huffman_node*pre=queue->front;
	Huffman_node*ptr=queue->front->next;
	//printf("pop ");
	while(ptr!=target){
		pre=ptr;
		ptr=ptr->next;
	}
	pre->next=ptr->next;
	if(queue->rear==target){
		queue->rear=pre;
	}
}


long long int RLE_code_len(RLE_node*RLE_code){ 
	RLE_node*ptr=RLE_code->next;
	long long int len=0;
	while(ptr!=NULL){//計算每個RLE list的長度 
		len++;
		ptr=ptr->next;	
	}
	return len;	
}


void un_RLE_8x8(RLE_node*list,double*F_zigzag){
	int i=0;
	RLE_node*zero_num=list->next;
	RLE_node*ptr=zero_num->next;
	while(!(ptr->val==0&&zero_num->val==0)){
		int count;
		for(count=1;count<=zero_num->val;count++){
			F_zigzag[i]=0;
			i++;
		}
		F_zigzag[i]=ptr->val;
		zero_num=ptr->next;
		ptr=zero_num->next;
		i++;
	}
	while(i<64){
		F_zigzag[i]=0;
		i++;
	}
}

void RLE_8x8(double*F_zigzag,RLE_node*list){
	int i=0;
	RLE_node*last=list;
	int count_zero=0;
	while(i<64){
		//printf("%.1lf ",F_zigzag[i]);
		if(F_zigzag[i]!=0){
			RLE_node*zero_num=new_RLE_node((1.0)*count_zero);//前面的非0值到這個非0值有幾個0 
			RLE_node*not_zero=new_RLE_node(F_zigzag[i]);
			add_RLE_node(last,zero_num);
			last=zero_num;
			add_RLE_node(last,not_zero);
			last=not_zero;
			count_zero=0;
			
		}
		else{
			count_zero++;
		}
		i++;
	}
	RLE_node*end_block1=new_RLE_node(0);
	RLE_node*end_block2=new_RLE_node(0);
	last->next=end_block1;
	last->next->next=end_block2;
}


RLE_node*new_RLE_node(double val){
	RLE_node*res=(RLE_node*)malloc(sizeof(RLE_node));
	res->val=val;
	res->next=NULL;
	return res;
}


void add_RLE_node(RLE_node*last,RLE_node*newnode){
	last->next=newnode;
}


RLE_node***malloc_3D_RLE_node(M,N){
	RLE_node***res=(RLE_node***)malloc(M*sizeof(RLE_node**));
	int i,j;
	for(i=0;i<M;i++){
		res[i]=(RLE_node**)malloc(N*sizeof(RLE_node*));
		for(j=0;j<N;j++){
			res[i][j]=(RLE_node*)malloc(sizeof(RLE_node));
		}
	}
	return res;
}
double***malloc_3D_double(int size1,int size2,int size3){
	double***res=(double***)malloc(size1*sizeof(double**));
	int i,j;
	for(i=0;i<size1;i++){
		res[i]=(double**)malloc(size2*sizeof(double*));
		for(j=0;j<size2;j++){
			res[i][j]=(double*)malloc(size3*sizeof(double));
		}
	}
	return res;
}
void ycbcr2bgr(Img*img,double**y_data,double**cb_data,double**cr_data){
	int i,j;
	for(i=0;i<(img->header).height;i++){
		for(j=0;j<(img->header).width;j++){
			//printf("%.2lf %.2lf %.2lf ",y_data[i][j],cb_data[i][j],cr_data[i][j]);
			double Y=y_data[i][j];
			double Cb=cb_data[i][j];
			double Cr=cr_data[i][j];
			(img->data[i][j]).data[0]=(unsigned int)((Y + 1.772*Cb)+0.5);
			(img->data[i][j]).data[1]=(unsigned int)((Y-0.344*Cb-0.714*Cr)+0.5);
			(img->data[i][j]).data[2]=(unsigned int)((Y + 1.402*Cr)+0.5);
		}
	}
}


void bgr2ycbcr(Img*img,double**y_data,double**cb_data,double**cr_data){
	int i,j;
	for(i=0;i<(img->header).height;i++){
		for(j=0;j<(img->header).width;j++){
			int B=(img->data[i][j]).data[0];
			int G=(img->data[i][j]).data[1];
			int R=(img->data[i][j]).data[2];
			y_data[i][j]=0.299*R + 0.587*G + 0.114*B;
			cb_data[i][j]=0.564*(B - y_data[i][j]);
			cr_data[i][j]=0.713*(R - y_data[i][j]);
		}
	}
}

double**malloc_2D_double(int row_size,int col_size){
	int i,j;
	double**res=(double**)malloc(row_size*sizeof(double*));
	for(i=0;i<row_size;i++){
		res[i]=(double*)malloc(col_size*sizeof(double));
	}
	return res;
}

//generate four-dimensions matrix of double
double****malloc_4D(int size1,int size2,int size3,int size4){
	int i1;
	double****res;
	res=(double****)malloc(size1*sizeof(double***));
	for(i1=0;i1<size1;i1++){
		res[i1]=(double***)malloc(size2*sizeof(double**));
		int i2;
		for(i2=0;i2<size2;i2++){
			res[i1][i2]=(double**)malloc(size3*sizeof(double*));
			int i3;
			for(i3=0;i3<size3;i3++){
				res[i1][i2][i3]=(double*)malloc(size4*sizeof(double));
				
			}
		}
	}
	return res;
}

//used for debug,print 2D double matrix 
void print_matrix(double**matrix,int row,int col){
	int i;
	for(i=0;i<row;i++){
		int j;
		for(j=0;j<col;j++){
			printf("%lf ",matrix[i][j]);
		}	
		printf("\n");
	}
}

//make a two-dimensional array of double,and its values are 1
double**ones(int row,int col){
	double**res=(double**)malloc(row*sizeof(double*));
	int i;
	for(i=0;i<row;i++){
		res[i]=(double*)malloc(col*sizeof(double));
		int j;
		for(j=0;j<col;j++){
			res[i][j]=1;
		}
	}
	return res;
}

//DCT 8x8 double matrix to frequency domain
void DCT_8x8(double matrix_8x8[8][8],double****F,double**beta,double ****basic_vector,int m,int n,int u,int v){
	int r;
	F[m][n][u][v]=0;
	for(r=0;r<8;r++){
		int c;
		for(c=0;c<8;c++){
			F[m][n][u][v]+=2.0/8*beta[u][0]*beta[v][0]*matrix_8x8[r][c]*basic_vector[u][v][r][c];
		}
	}
}

//IDCT frequency domain to 2D double matrix 
void IDCT_8x8(double**F,double**matrix_8x8,double**beta,double****basic_vector,int r,int c){
	int u,v;
	matrix_8x8[r][c]=0;
	for(u=0;u<8;u++){
		for(v=0;v<8;v++){
			matrix_8x8[r][c]+=2.0/8*F[u][v]*beta[u][0]*beta[v][0]*basic_vector[u][v][r][c];
		}
	}
}

//Cut out an 8x8 double matrix from the image data matrix
void cut_8x8(double**data,double matrix_8x8[8][8],int m,int n){
	int i;//byte_index為某個像素的第幾個byte 
	for(i=0;i<8;i++){
		int j;
		for(j=0;j<8;j++){
			matrix_8x8[i][j]=data[8*m+i][8*n+j];//(double)data[8*m+i][8*n+j].data[byte_index];
		}
	}
}

/*read header*/
void readheader(FILE* fp, Bitmap *x) {
	fread(&x->identifier, sizeof(x->identifier), 1, fp);
	fread(&x->filesize, sizeof(x->filesize), 1, fp);
	fread(&x->reserved, sizeof(x->reserved), 1, fp);
	fread(&x->reserved2, sizeof(x->reserved2), 1, fp);
	fread(&x->bitmap_dataoffset, sizeof(x->bitmap_dataoffset), 1, fp);
	fread(&x->bitmap_headersize, sizeof(x->bitmap_headersize), 1, fp);
	fread(&x->width, sizeof(x->width), 1, fp);
	fread(&x->height, sizeof(x->height), 1, fp);
	fread(&x->planes, sizeof(x->planes), 1, fp);
	fread(&x->bits_perpixel, sizeof(x->bits_perpixel), 1, fp);
	fread(&x->compression, sizeof(x->compression), 1, fp);
	fread(&x->bitmap_datasize, sizeof(x->bitmap_datasize), 1, fp);
	fread(&x->hresolution, sizeof(x->hresolution), 1, fp);
	fread(&x->vresolution, sizeof(x->vresolution), 1, fp);
	fread(&x->usedcolors, sizeof(x->usedcolors), 1, fp);
	fread(&x->importantcolors, sizeof(x->importantcolors), 1, fp);
}

//generate a header of basis file 
Bitmap make_basis_header(){
	Bitmap header;
	header.identifier[0]='B';header.identifier[1]='M';
	header.filesize=1142;
	header.reserved=0;
	header.reserved2=0;
	header.bitmap_dataoffset=1078;
	header.bitmap_headersize=40;
	header.width=8;
	header.height=8;
	header.planes=1;
	header.bits_perpixel=8;
	header.compression=0;
	header.bitmap_datasize=64;
	header.hresolution=0;
	header.vresolution=0;
	header.usedcolors=256;
	header.importantcolors=0;
	return header;
}

//generate a palette of basis file 
unsigned char*make_basis_palette(unsigned short palette_size){
	unsigned char*palette=(unsigned char*)malloc(palette_size*sizeof(unsigned char));
	int i;
	int val=0;
	for(i=0;i<palette_size;i+=4,val++){
		palette[i]=palette[i+1]=palette[i+2]=val;
		palette[i+3]=0;
	}
	return palette;
}

//generate a uint8_t 2D matrix 
uint8_t**malloc_uint8_2D(int row,int col){
	uint8_t**res=(uint8_t**)malloc(row*sizeof(uint8_t*));
	int i;
	for(i=0;i<row;i++){
		res[i]=(uint8_t*)malloc(col*sizeof(uint8_t));
	}
	return res;
}

//copy data of double matrix x to uint8_t matrix bmp_x
void uint8(double x[8][8],uint8_t**bmp_x){
	int i;
	for(i=0;i<8;i++){
		int j;
		for(j=0;j<8;j++){
			if(x[i][j]>=255){
				bmp_x[i][j]=255;
			}
			else
				bmp_x[i][j]=(uint8_t)x[i][j];
		}
	}	
}

//Turn the row of a pixel matrix upside down
void reverse_row(Pixel**arr,int rowSize,int colSize){
	int i1,i2;
	int j=0;
	i1=0;
	i2=rowSize-1;
	while(i1<i2){
		for(j=0;j<colSize;j++){
			Pixel tmp=arr[i1][j];
			arr[i1][j]=arr[i2][j];
			arr[i2][j]=tmp;
		}
		i1++;
		i2--;
	}
}

/* A function of making two dimensions memory locate array*/
Pixel** malloc_2D_Pixel(int row, int col,unsigned short bytes_perpixel){
	Pixel**arr=(Pixel**)malloc(row*sizeof(Pixel*));
	int i;
	for(i=0;i<row;i++){
		arr[i]=(Pixel*)malloc(col*sizeof(Pixel));
		int j;
		for(j=0;j<col;j++){
			arr[i][j].data=(unsigned char*)malloc(bytes_perpixel*sizeof(unsigned char));
		}
	}
	return arr;
}
//read bmp file header、palette and data to a variable of type Img 
Img ImRead(char input_name[]){
	Img res;
	FILE *fp_in;
	char*fn_in= input_name;
	fp_in = fopen(fn_in, "rb");
	//read header
	Bitmap bmpheader;
	readheader(fp_in, &bmpheader);
	unsigned short bytes_perpixel=bmpheader.bits_perpixel/8;
	int H = bmpheader.height;
	int W = bmpheader.width;
	/*
	The number of bytes in each column of img must be a multiple of 4, 
	if it is insufficient, padding 0.
	*/
	int skip = (4 - (bmpheader.width * bytes_perpixel) % 4); 
	if (skip == 4) skip = 0;
	char skip_buff[3]={0,0,0};
	
	//read palette
	int palette_size=bmpheader.bitmap_dataoffset-54;
	unsigned char*palette=(unsigned char*)malloc((palette_size)*sizeof(unsigned char));
	int i;
	for(i=0;i<palette_size;i++){
		fread(palette+i,1,1,fp_in);
	}
	
	//read data
	Pixel**pix_data=malloc_2D_Pixel(bmpheader.height, bmpheader.width,bytes_perpixel);
	for(i=0;i<bmpheader.height;i++){
		int j;
		for(j=0;j<bmpheader.width;j++){
			fread(pix_data[i][j].data,1,bytes_perpixel,fp_in);
		}
		if(skip!=0)fread(skip_buff,skip,1,fp_in);
	}
	res.header=bmpheader;
	res.palette=palette;
	/*
	Normally, pixels are stored from bottom to top and left to right, 
	so we should turn the row of image data upside down
	*/
	reverse_row(pix_data,bmpheader.height,bmpheader.width);
	res.data=pix_data;
	return res;
}

//write image header、palette and data to a bmp file 
void ImWrite(Img img,char filename[]){
	FILE *outfile;
    outfile= fopen(filename, "wb");
    reverse_row(img.data,img.header.height,img.header.width);
    
	//write header to file
	fwrite(&img.header.identifier, sizeof(short), 1, outfile);
	fwrite(&img.header.filesize, sizeof(int), 1, outfile);
	fwrite(&img.header.reserved, sizeof(short), 1, outfile);
	fwrite(&img.header.reserved2, sizeof(short), 1, outfile);
	fwrite(&img.header.bitmap_dataoffset, sizeof(int), 1, outfile);
	fwrite(&img.header.bitmap_headersize, sizeof(int), 1, outfile);
	fwrite(&img.header.width, sizeof(int), 1, outfile);
	fwrite(&img.header.height, sizeof(int), 1, outfile);
	fwrite(&img.header.planes, sizeof(short), 1, outfile);
	fwrite(&img.header.bits_perpixel, sizeof(short), 1, outfile);
	fwrite(&img.header.compression, sizeof(int), 1, outfile);
	fwrite(&img.header.bitmap_datasize, sizeof(int), 1, outfile);
	fwrite(&img.header.hresolution, sizeof(int), 1, outfile);
	fwrite(&img.header.vresolution, sizeof(int), 1, outfile);
	fwrite(&img.header.usedcolors, sizeof(int), 1, outfile);
	fwrite(&img.header.importantcolors, sizeof(int), 1, outfile);
	
	//write palette to file
	int ptr=0;
	int palette_size=img.header.bitmap_dataoffset-54;
	while(ptr<palette_size){
		fwrite(img.palette+ptr,1,1,outfile);
		fwrite(img.palette+ptr+1,1,1,outfile);
		fwrite(img.palette+ptr+2,1,1,outfile);
		fwrite(img.palette+ptr+3,1,1,outfile);
		ptr+=4;
	}
	
	//write data to file
	char skip_buf[3] = { 0, 0, 0 };
	unsigned short bytes_perpixel=img.header.bits_perpixel/8;
	int skip = (4 - (img.header.width*bytes_perpixel) % 4); 
	if (skip == 4) skip = 0;
	int i;
	for(i=0;i<img.header.height;i++){
		int j;
		for(j=0;j<img.header.width;j++){
			fwrite(img.data[i][j].data,bytes_perpixel,1,outfile);
		}
		if(skip!=0){
			fwrite(skip_buf, skip, 1, outfile);
		}
	}
}


