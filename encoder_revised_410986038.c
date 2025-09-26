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
	Pixel**pixels;
}Img;

typedef struct RLE_node{
	short val;
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

typedef struct codebook{
	short*vals;
	char**bitstream;
	int node_num;
}codebook;



int RLE_code_len(RLE_node*RLE_code);
//檔案讀寫相關 
void saveas(unsigned char**matrix,int row,int col,char*filename);
Img ImRead(char input_name[]);
void ImWrite(Img img,char filename[]);
void readheader(FILE* fp, Bitmap *x);
void savetxt_table_8x8(int arr[8][8],char*filename);
FILE*save_8x8_RLE_txt(FILE*fptr,short*list);
FILE* save_header_txt(Img*img,FILE*fp);
FILE*save_8x8_RLE_bin(FILE*fptr,short*list);
FILE*RLE2Huffman_txt(codebook*dict,short*RLE,FILE*fptr);
FILE*RLE2Huffman_bin(codebook*dict,short*RLE,FILE*fptr,int*bin_filesize);//資料結構建立相關 
unsigned char**malloc_2D_uchar(int row,int col);
void reverse_row(Pixel**arr,int rowSize,int colSize);
unsigned char*make_basis_palette(unsigned short palette_size);
Pixel** malloc_2D_Pixel(int row, int col,unsigned short bytes_perpixel);
double**malloc_2D_double(int row_size,int col_size);
void free_2D_double(double**data,int size1,int size2);
double***malloc_3D_double(int size1,int size2,int size3);
double****malloc_4D_double(int size1,int size2,int size3,int size4);
void free_4D_double(double****data,int size1,int size2,int size3,int size4);
short****malloc_4D_short(int size1,int size2,int size3,int size4);
void free_4D_short(short****data,int size1,int size2,int size3,int size4);
short***malloc_3D_short(int size1,int size2,int size3);
void free_3D_short(short***data,int size1,int size2,int size3);
RLE_node***malloc_3D_RLE_node(int M,int N);
void free_3D_RLE_node(RLE_node***data,int size1,int size2,int size3);
void RLE_8x8(short*F_zigzag,short*RLE);
RLE_node*new_RLE_node(short val);
void add_RLE_node(RLE_node*last,RLE_node*newnode);
double**ones(int row,int col);
void cut_8x8(double**data,double matrix_8x8[8][8],int m,int n);

//公式相關
void bgr2ycbcr(Img*img,double**y_data,double**cb_data,double**cr_data);
void DCT_8x8(double matrix_8x8[8][8],double****F,double**beta,double ****basic_vector,int m,int n,int u,int v);
void SQNR(double****F,double****e,int H,int W);
void show_channel_compress_ratio(int M,int N,short***RLE_lists,short***y_F_zigzag,int*RLE_size,int*zigzag_size);
//Huffman相關 
codebook huffmandict(short***y,short***cb,short***cr,int M,int N,char*word_filename);
void print_queue(Huffman_queue*queue);
Huffman_queue Huffman_queue_init();
void Huffman_queue_push(Huffman_queue*queue,Huffman_node*node);
void Huffman_queue_pop(Huffman_queue*queue,Huffman_node*target);
void print_tree(Huffman_node*node,int*len);
void save_codebook_txt(FILE**word_fptr,short*vals,int*val_ptr,char**bitstream,Huffman_node*node,char*buff,int buff_ptr);

//FILE*debug_fptr=fopen("debug_encoder.txt","w");

int main(int argc, char **argv){
	if(argv[1][0]=='0'){
		//printf("write RGB\n");
		Img img=ImRead(argv[2]);
		unsigned char**red=malloc_2D_uchar(img.header.height,img.header.width);
		unsigned char**green=malloc_2D_uchar(img.header.height,img.header.width);
		unsigned char**blue=malloc_2D_uchar(img.header.height,img.header.width);
		//printf("height=%d\n",img.header.height);
		//printf("width=%d",img.header.width);
		int i,j;
		for(i=0;i<img.header.height;i++){
			for(j=0;j<img.header.width;j++){
				blue[i][j]=img.pixels[i][j].data[0];
				green[i][j]=img.pixels[i][j].data[1];
				red[i][j]=img.pixels[i][j].data[2];
			} 
		}
		
		//printf("test\n");
		saveas(red,img.header.height,img.header.width,argv[3]);
		saveas(green,img.header.height,img.header.width,argv[4]);
		saveas(blue,img.header.height,img.header.width,argv[5]);
		FILE*fp=fopen(argv[6],"w");
		save_header_txt(&img,fp); 
	}
	
	else if(argv[1][0]=='1'){
		//將q table寫入txt檔 
			int y_Q[8][8]={{16,11,10,16,24,40,51,61},
     					{12,12,14,19,26,58,60,55},
     					{14,13,16,24,40,57,69,56},
     					{14,17,22,29,51,87,80,62},
     					{18,22,37,56,68,109,103,77},
     					{24,35,55,64,81,104,113,92},
     					{49,64,78,87,103,121,120,101},
     					{72,92,95,98,112,100,103,99}};
     	
    		int cbcr_Q[8][8]={{ 17, 18, 24, 47, 99, 99, 99, 99 },
						{ 18, 21, 26, 66, 99, 99, 99, 99 },
						{ 24, 26, 56, 99, 99, 99, 99, 99 },
						{ 47, 66, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99 } };
		savetxt_table_8x8(y_Q,argv[3]);
		savetxt_table_8x8(cbcr_Q,argv[4]);
		savetxt_table_8x8(cbcr_Q,argv[5]);
		 
		//encoder 1 Kimberly.bmp 3~5table dim.txt 7~9qF 10~12eF
		Img img=ImRead(argv[2]);

		double**y_data=malloc_2D_double(img.header.height,img.header.width);
		double**cb_data=malloc_2D_double(img.header.height,img.header.width);
		double**cr_data=malloc_2D_double(img.header.height,img.header.width);
		
		//bgr to ycbcr 
		bgr2ycbcr(&(img),y_data,cb_data,cr_data);
		
		
		int H=8;
		int W=8;
		int M=img.header.height/8;
		int N=img.header.width/8;
		//generate basis vector for 2D-DCT
		double****basic_vector=malloc_4D_double(8,8,8,8);
		int u,v,r,c;
		/*
		Because we want to put the basic vector picture in the directory we specified, 
		we need to change the workspace to the directory we specified
		*/
		for(u=0;u<H;u++){
			for(v=0;v<W;v++){
				//Make a 8x8 basic vector matrix x
				for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
					for(c=0;c<W;c++){ 
						basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
					}
					
				}
				
			}
		}
		
		//2D-DCT
		double****y_F=malloc_4D_double(M,N,8,8);//frequency domain of the 8x8 matrices
		double****cb_F=malloc_4D_double(M,N,8,8);
		double****cr_F=malloc_4D_double(M,N,8,8);
		double**beta=ones(H,1);
		beta[0][0]=1/sqrt(2);
		int m,n;
		//Cut the matrix into multiple 8x8 matrices
		for(m=0;m<M;m++){
			for(n=0;n<N;n++){
			
				//Get one of the 8x8 matrices
				double y_matrix_8x8[8][8];
				double cb_matrix_8x8[8][8];
				double cr_matrix_8x8[8][8];
				cut_8x8(y_data,y_matrix_8x8,m,n);
				cut_8x8(cb_data,cb_matrix_8x8,m,n);
				cut_8x8(cr_data,cr_matrix_8x8,m,n);
				//subtract 128
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
					
					}
				}
				
			}
		} 
		
		
		free_4D_double(basic_vector,8,8,8,8);
		
		//計算並存取量化結果 
		FILE*q_y_fptr=fopen(argv[7],"wb");
		FILE*q_cb_fptr=fopen(argv[8],"wb");
		FILE*q_cr_fptr=fopen(argv[9],"wb");
		short****q_y_F=malloc_4D_short(M,N,8,8);
		short****q_cb_F=malloc_4D_short(M,N,8,8);
		short****q_cr_F=malloc_4D_short(M,N,8,8);
		//FILE*debug_fptr=fopen("debug_encoder.txt","w"); //debug
		int debug_count=1;
		for(m=0;m<M;m++){
			int n;
			for(n=0;n<N;n++){
				int u;
				//quantize and round each 8x8 matrix
				for(u=0;u<8;u++){
					int v;
					for(v=0;v<8;v++){
						q_y_F[m][n][u][v]=(short)round(y_F[m][n][u][v]/y_Q[u][v]);
						fwrite(&q_y_F[m][n][u][v],sizeof(short),1,q_y_fptr);
						q_cb_F[m][n][u][v]=(short)round(cb_F[m][n][u][v]/cbcr_Q[u][v]);
						fwrite(&q_cb_F[m][n][u][v],sizeof(short),1,q_cb_fptr);
						q_cr_F[m][n][u][v]=(short)round(cr_F[m][n][u][v]/cbcr_Q[u][v]);
						fwrite(&q_cr_F[m][n][u][v],sizeof(short),1,q_cr_fptr);
						//fprintf(debug_fptr,"%d %d %d\n",q_y_F[m][n][u][v],q_cb_F[m][n][u][v],q_cr_F[m][n][u][v]);
						/*if(debug_count==130){
							
						}*/
					}
				}
			}
		}

		//計算並存取量化誤差 
		FILE*e_y_fptr=fopen(argv[10],"wb");
		FILE*e_cb_fptr=fopen(argv[11],"wb");
		FILE*e_cr_fptr=fopen(argv[12],"wb");
		
		double****e_y=malloc_4D_double(M,N,8,8);
		double****e_cb=malloc_4D_double(M,N,8,8);
		double****e_cr=malloc_4D_double(M,N,8,8);
		for(m=0;m<M;m++){
			int n;
			for(n=0;n<N;n++){
				int u;
				for(u=0;u<8;u++){
					int v;
					for(v=0;v<8;v++){
						float ytmp,cbtmp,crtmp;
						ytmp=y_F[m][n][u][v]-(float)q_y_F[m][n][u][v]*y_Q[u][v];
						e_y[m][n][u][v]=ytmp;
						fwrite(&ytmp,sizeof(float),1,e_y_fptr);
						
						cbtmp=cb_F[m][n][u][v]-(float)q_cb_F[m][n][u][v]*cbcr_Q[u][v];
						e_cb[m][n][u][v]=cbtmp;
						fwrite(&cbtmp,sizeof(float),1,e_cb_fptr);
						crtmp=cr_F[m][n][u][v]-(double)q_cr_F[m][n][u][v]*cbcr_Q[u][v];
						e_cr[m][n][u][v]=crtmp;
						fwrite(&crtmp,sizeof(float),1,e_cr_fptr);
					}
				}
			}
		}
		
		//寫入header進dim.txt
		FILE*fp=fopen(argv[6],"w");
		save_header_txt(&img,fp); 
		fclose(fp);
		
		free_4D_short(q_y_F,M,N,8,8);
		free_4D_short(q_cb_F,M,N,8,8);
		free_4D_short(q_cr_F,M,N,8,8);
		//SQNR
		printf("SQNR of Y:\n");
		SQNR(y_F,e_y,M,N);
		printf("SQNR of CB:\n");
		SQNR(cb_F,e_cb,M,N);
		printf("SQNR of CR:\n");
		SQNR(cr_F,e_cr,M,N);
		free_2D_double(y_data,img.header.height,img.header.width);
		free_2D_double(cb_data,img.header.height,img.header.width);
		free_2D_double(cr_data,img.header.height,img.header.width);
	}
	
	
	else if(argv[1][0]=='2'){
		if(argv[3][0]=='a'){//encoder 2 Kimberly.bmp ascii rle_code.txt 
			//read image
			Img img=ImRead(argv[2]);
			int W=8;
			int H=8;
			int M=img.header.height/H;//圖片高有幾個8x8 matrix
			int N=img.header.width/W;//圖片寬有幾個8x8 matrix
	
			//bgr to ycbcr
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			bgr2ycbcr(&(img),y_data,cb_data,cr_data);
		
	
			//generate basis vector for 2D-DCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int u;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			for(u=0;u<H;u++){
				int v;
				for(v=0;v<W;v++){
					int r;
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						int c;
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
				
					}
			
				}
			}
			//2D-DCT
			double****y_F=malloc_4D_double(M,N,8,8);//frequency domain of the 8x8 matrices
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			double**beta=ones(H,1);
			beta[0][0]=1/sqrt(2);
			int m,n;
			//Cut the matrix into multiple 8x8 matrices
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
				
					//Get one of the 8x8 matrices
					double y_matrix_8x8[8][8];
					double cb_matrix_8x8[8][8];
					double cr_matrix_8x8[8][8];
					cut_8x8(y_data,y_matrix_8x8,m,n);
					cut_8x8(cb_data,cb_matrix_8x8,m,n);
					cut_8x8(cr_data,cr_matrix_8x8,m,n);
					//subtract 128
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_matrix_8x8[u][v]-=128;
							cb_matrix_8x8[u][v]-=128;
							cr_matrix_8x8[u][v]-=128;
						}
					}
					//Get the frequency domain of the 8x8 matrix with DCT
					for(u=0;u<H;u++){
						int v;
						for(v=0;v<W;v++){
							DCT_8x8(y_matrix_8x8,y_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cb_matrix_8x8,cb_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cr_matrix_8x8,cr_F,beta,basic_vector,m,n,u,v);
						}
					}
					
				}
			} 
		
		
			free_4D_double(basic_vector,8,8,8,8);
			free_2D_double(y_data,img.header.height,img.header.width);
			free_2D_double(cb_data,img.header.height,img.header.width);
			free_2D_double(cr_data,img.header.height,img.header.width);
			
			
			//quantization
			double y_Q[8][8]={{16,11,10,16,24,40,51,61},
    				 		{12,12,14,19,26,58,60,55},
     						{14,13,16,24,40,57,69,56},
     						{14,17,22,29,51,87,80,62},
     						{18,22,37,56,68,109,103,77},
     						{24,35,55,64,81,104,113,92},
     						{49,64,78,87,103,121,120,101},
     						{72,92,95,98,112,100,103,99}};
     
    			double cbcr_Q[8][8]={ 	{ 17, 18, 24, 47, 99, 99, 99, 99 },
						{ 18, 21, 26, 66, 99, 99, 99, 99 },
						{ 24, 26, 56, 99, 99, 99, 99, 99 },
						{ 47, 66, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99 } };
    		
			short****q_y_F=malloc_4D_short(M,N,8,8);
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					//quantize and round each 8x8 matrix
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							q_y_F[m][n][u][v]=(short)round(y_F[m][n][u][v]/y_Q[u][v]);
							q_cb_F[m][n][u][v]=(short)round(cb_F[m][n][u][v]/cbcr_Q[u][v]);
							q_cr_F[m][n][u][v]=(short)round(cr_F[m][n][u][v]/cbcr_Q[u][v]);
						}
					}
				}
			}
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			
			//DC DPCM
			for(m=0;m<M;m++){
				int n;
				for(n=N-1;n>=1;n--){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]-q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]-q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]-q_cr_F[m][n-1][0][0];
				}
			}
			for(m=M-1;m>=1;m--){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]-q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]-q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]-q_cr_F[m-1][0][0][0];
			} 
			
			
			//zigzag
			int zz_matrix[8][8]={
				{1,2,6,7,15,16,28,29},
				{3,5,8,14,17,27,30,43},
				{4,9,13,18,26,31,42,44},
				{10,12,19,25,32,41,45,54},
				{11,20,24,33,40,46,53,55},
				{21,23,34,39,47,52,56,61},
				{22,35,38,48,51,57,60,62},
				{36,37,49,50,58,59,63,64}
			};
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F_zigzag[m][n][zz_matrix[r][c]-1]=q_y_F[m][n][r][c];
							q_cb_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cb_F[m][n][r][c];
							q_cr_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cr_F[m][n][r][c];
						}
					}
				}
			}
			
			free_4D_short(q_y_F,M,N,8,8);//frequency domain of the 8x8 matrices
			free_4D_short(q_cb_F,M,N,8,8);
			free_4D_short(q_cr_F,M,N,8,8);
			
			short***y_RLE_lists=malloc_3D_short(M,N,64);
			short***cb_RLE_lists=malloc_3D_short(M,N,64);
			short***cr_RLE_lists=malloc_3D_short(M,N,64);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					RLE_8x8(q_y_F_zigzag[m][n],y_RLE_lists[m][n]);
					
					
					RLE_8x8(q_cb_F_zigzag[m][n],cb_RLE_lists[m][n]);
					
					
					RLE_8x8(q_cr_F_zigzag[m][n],cr_RLE_lists[m][n]);
						
				}
			}
			
			
			FILE*fptr=fopen(argv[4],"w");
			fprintf(fptr,"%c%c ",img.header.identifier[0],img.header.identifier[1]);
			fprintf(fptr,"%d ",img.header.filesize);
			fprintf(fptr,"%d ",img.header.reserved);
			fprintf(fptr,"%d ",img.header.reserved2);
		  	fprintf(fptr,"%d ",img.header.bitmap_dataoffset);
			fprintf(fptr,"%d ",img.header.bitmap_headersize);
			fprintf(fptr,"%d ",img.header.width);
			fprintf(fptr,"%d ",img.header.height);
			fprintf(fptr,"%d ",img.header.planes);
			fprintf(fptr,"%d ",img.header.bits_perpixel);
			fprintf(fptr,"%d ",img.header.compression);
			fprintf(fptr,"%d ",img.header.bitmap_datasize);
			fprintf(fptr,"%d ",img.header.hresolution);
			fprintf(fptr,"%d ",img.header.vresolution);
			fprintf(fptr,"%d ",img.header.usedcolors);
			fprintf(fptr,"%d\n",img.header.importantcolors);
			
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					fprintf(fptr,"(%d,%d,Y) ",m,n);
					fptr=save_8x8_RLE_txt(fptr,y_RLE_lists[m][n]);
					fprintf(fptr,"(%d,%d,cb) ",m,n);
					fptr=save_8x8_RLE_txt(fptr,cb_RLE_lists[m][n]);
					fprintf(fptr,"(%d,%d,cr) ",m,n);
					fptr=save_8x8_RLE_txt(fptr,cr_RLE_lists[m][n]);
				}
			}
			fclose(fptr);
			free_3D_short(y_RLE_lists,M,N,64);
			free_3D_short(cb_RLE_lists,M,N,64);
			free_3D_short(cr_RLE_lists,M,N,64);
		}
		
		if(argv[3][0]=='b'){//encoder 2 Kimberly.bmp binary rle_code.bin 
			//read image
			Img img=ImRead(argv[2]);
			int W=8;
			int H=8;
			int M=img.header.height/H;//圖片高有幾個8x8 matrix
			int N=img.header.width/W;//圖片寬有幾個8x8 matrix
	
			//bgr to ycbcr
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			bgr2ycbcr(&(img),y_data,cb_data,cr_data);
		
	
			//generate basis vector for 2D-DCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int u;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			for(u=0;u<H;u++){
				int v;
				for(v=0;v<W;v++){
					int r;
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						int c;
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
				
					}
			
				}
			}
			//2D-DCT
			double****y_F=malloc_4D_double(M,N,8,8);//frequency domain of the 8x8 matrices
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			double**beta=ones(H,1);
			beta[0][0]=1/sqrt(2);
			int m,n;
			//Cut the matrix into multiple 8x8 matrices
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
				
					//Get one of the 8x8 matrices
					double y_matrix_8x8[8][8];
					double cb_matrix_8x8[8][8];
					double cr_matrix_8x8[8][8];
					cut_8x8(y_data,y_matrix_8x8,m,n);
					cut_8x8(cb_data,cb_matrix_8x8,m,n);
					cut_8x8(cr_data,cr_matrix_8x8,m,n);
					//subtract 128
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_matrix_8x8[u][v]-=128;
							cb_matrix_8x8[u][v]-=128;
							cr_matrix_8x8[u][v]-=128;
						}
					}
					//Get the frequency domain of the 8x8 matrix with DCT
					for(u=0;u<H;u++){
						int v;
						for(v=0;v<W;v++){
							DCT_8x8(y_matrix_8x8,y_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cb_matrix_8x8,cb_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cr_matrix_8x8,cr_F,beta,basic_vector,m,n,u,v);
						}
					}
					
				}
			} 
		
		
			free_4D_double(basic_vector,8,8,8,8);
			free_2D_double(y_data,img.header.height,img.header.width);
			free_2D_double(cb_data,img.header.height,img.header.width);
			free_2D_double(cr_data,img.header.height,img.header.width);
			
			
			//quantization
			double y_Q[8][8]={{16,11,10,16,24,40,51,61},
    				 		{12,12,14,19,26,58,60,55},
     						{14,13,16,24,40,57,69,56},
     						{14,17,22,29,51,87,80,62},
     						{18,22,37,56,68,109,103,77},
     						{24,35,55,64,81,104,113,92},
     						{49,64,78,87,103,121,120,101},
     						{72,92,95,98,112,100,103,99}};
     
    			double cbcr_Q[8][8]={ 	{ 17, 18, 24, 47, 99, 99, 99, 99 },
						{ 18, 21, 26, 66, 99, 99, 99, 99 },
						{ 24, 26, 56, 99, 99, 99, 99, 99 },
						{ 47, 66, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99 } };
    		
			short****q_y_F=malloc_4D_short(M,N,8,8);
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					//quantize and round each 8x8 matrix
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							q_y_F[m][n][u][v]=(short)round(y_F[m][n][u][v]/y_Q[u][v]);
							q_cb_F[m][n][u][v]=(short)round(cb_F[m][n][u][v]/cbcr_Q[u][v]);
							q_cr_F[m][n][u][v]=(short)round(cr_F[m][n][u][v]/cbcr_Q[u][v]);
						}
					}
				}
			}
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			
			//DC DPCM
			for(m=0;m<M;m++){
				int n;
				for(n=N-1;n>=1;n--){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]-q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]-q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]-q_cr_F[m][n-1][0][0];
				}
			}
			for(m=M-1;m>=1;m--){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]-q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]-q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]-q_cr_F[m-1][0][0][0];
			} 
			
			
			//zigzag
			int zz_matrix[8][8]={
				{1,2,6,7,15,16,28,29},
				{3,5,8,14,17,27,30,43},
				{4,9,13,18,26,31,42,44},
				{10,12,19,25,32,41,45,54},
				{11,20,24,33,40,46,53,55},
				{21,23,34,39,47,52,56,61},
				{22,35,38,48,51,57,60,62},
				{36,37,49,50,58,59,63,64}
			};
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F_zigzag[m][n][zz_matrix[r][c]-1]=q_y_F[m][n][r][c];
							q_cb_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cb_F[m][n][r][c];
							q_cr_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cr_F[m][n][r][c];
						}
					}
				}
			}
			
			free_4D_short(q_y_F,M,N,8,8);//frequency domain of the 8x8 matrices
			free_4D_short(q_cb_F,M,N,8,8);
			free_4D_short(q_cr_F,M,N,8,8);
			
			short***y_RLE_lists=malloc_3D_short(M,N,64);
			short***cb_RLE_lists=malloc_3D_short(M,N,64);
			short***cr_RLE_lists=malloc_3D_short(M,N,64);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					RLE_8x8(q_y_F_zigzag[m][n],y_RLE_lists[m][n]);
					
					
					RLE_8x8(q_cb_F_zigzag[m][n],cb_RLE_lists[m][n]);
					
					
					RLE_8x8(q_cr_F_zigzag[m][n],cr_RLE_lists[m][n]);
						
				}
			}

			FILE*fptr=fopen(argv[4],"wb");
			fwrite(&img.header.identifier, sizeof(short), 1, fptr);
			fwrite(&img.header.filesize, sizeof(int), 1, fptr);
			fwrite(&img.header.reserved, sizeof(short), 1, fptr);
			fwrite(&img.header.reserved2, sizeof(short), 1, fptr);
			fwrite(&img.header.bitmap_dataoffset, sizeof(int), 1, fptr);
			fwrite(&img.header.bitmap_headersize, sizeof(int), 1, fptr);
			fwrite(&img.header.width, sizeof(int), 1, fptr);
			fwrite(&img.header.height, sizeof(int), 1, fptr);
			fwrite(&img.header.planes, sizeof(short), 1, fptr);
			fwrite(&img.header.bits_perpixel, sizeof(short), 1, fptr);
			fwrite(&img.header.compression, sizeof(int), 1, fptr);
			fwrite(&img.header.bitmap_datasize, sizeof(int), 1, fptr);
			fwrite(&img.header.hresolution, sizeof(int), 1, fptr);
			fwrite(&img.header.vresolution, sizeof(int), 1, fptr);
			fwrite(&img.header.usedcolors, sizeof(int), 1, fptr);
			fwrite(&img.header.importantcolors, sizeof(int), 1, fptr);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					fptr=save_8x8_RLE_bin(fptr,y_RLE_lists[m][n]);
					fptr=save_8x8_RLE_bin(fptr,cb_RLE_lists[m][n]);
					fptr=save_8x8_RLE_bin(fptr,cr_RLE_lists[m][n]);
				}
			}
			
			fclose(fptr);
			
			//show channel compress ratio	
			int RLE_y_size,RLE_cb_size,RLE_cr_size;
			int zigzag_y_size,zigzag_cb_size,zigzag_cr_size;
			printf("channel Y ");
			show_channel_compress_ratio(M,N,y_RLE_lists,q_y_F_zigzag,&RLE_y_size,&zigzag_y_size);
			printf("channel cb ");
			show_channel_compress_ratio(M,N,cb_RLE_lists,q_cb_F_zigzag,&RLE_cb_size,&zigzag_cb_size);
			printf("channel cr ");
			show_channel_compress_ratio(M,N,cr_RLE_lists,q_cr_F_zigzag,&RLE_cr_size,&zigzag_cr_size);
			
			//show total compress ratio
			printf("total compress ratio: %lf\n",(double)(RLE_y_size+RLE_cb_size+RLE_cr_size+54)/(img.header.filesize));
			
			free_3D_short(q_y_F_zigzag,M,N,64);
			free_3D_short(q_cb_F_zigzag,M,N,64);
			free_3D_short(q_cr_F_zigzag,M,N,64);
			free_3D_short(y_RLE_lists,M,N,64);
			free_3D_short(cb_RLE_lists,M,N,64);
			free_3D_short(cr_RLE_lists,M,N,64);
			
			return 0;
			
		}
	}
	else if(argv[1][0]=='3'){
		if(argv[3][0]=='a'){//encoder 3 Kimberly.bmp ascii codebook.txt huffman_code.txt 
			//read image
			Img img=ImRead(argv[2]);
			int W=8;
			int H=8;
			int M=img.header.height/H;//圖片高有幾個8x8 matrix
			int N=img.header.width/W;//圖片寬有幾個8x8 matrix
	
			//bgr to ycbcr
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			bgr2ycbcr(&(img),y_data,cb_data,cr_data);
		
	
			//generate basis vector for 2D-DCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int u;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			for(u=0;u<H;u++){
				int v;
				for(v=0;v<W;v++){
					int r;
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						int c;
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
				
					}
			
				}
			}
			//2D-DCT
			double****y_F=malloc_4D_double(M,N,8,8);//frequency domain of the 8x8 matrices
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			double**beta=ones(H,1);
			beta[0][0]=1/sqrt(2);
			int m,n;
			//Cut the matrix into multiple 8x8 matrices
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
				
					//Get one of the 8x8 matrices
					double y_matrix_8x8[8][8];
					double cb_matrix_8x8[8][8];
					double cr_matrix_8x8[8][8];
					cut_8x8(y_data,y_matrix_8x8,m,n);
					cut_8x8(cb_data,cb_matrix_8x8,m,n);
					cut_8x8(cr_data,cr_matrix_8x8,m,n);
					//subtract 128
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_matrix_8x8[u][v]-=128;
							cb_matrix_8x8[u][v]-=128;
							cr_matrix_8x8[u][v]-=128;
						}
					}
					//Get the frequency domain of the 8x8 matrix with DCT
					for(u=0;u<H;u++){
						int v;
						for(v=0;v<W;v++){
							DCT_8x8(y_matrix_8x8,y_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cb_matrix_8x8,cb_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cr_matrix_8x8,cr_F,beta,basic_vector,m,n,u,v);
						}
					}
					
				}
			} 
		
		
			free_4D_double(basic_vector,8,8,8,8);
			free_2D_double(y_data,img.header.height,img.header.width);
			free_2D_double(cb_data,img.header.height,img.header.width);
			free_2D_double(cr_data,img.header.height,img.header.width);
			
			
			//quantization
			double y_Q[8][8]={{16,11,10,16,24,40,51,61},
    				 		{12,12,14,19,26,58,60,55},
     						{14,13,16,24,40,57,69,56},
     						{14,17,22,29,51,87,80,62},
     						{18,22,37,56,68,109,103,77},
     						{24,35,55,64,81,104,113,92},
     						{49,64,78,87,103,121,120,101},
     						{72,92,95,98,112,100,103,99}};
     
    			double cbcr_Q[8][8]={ 	{ 17, 18, 24, 47, 99, 99, 99, 99 },
						{ 18, 21, 26, 66, 99, 99, 99, 99 },
						{ 24, 26, 56, 99, 99, 99, 99, 99 },
						{ 47, 66, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99 } };
    		
			short****q_y_F=malloc_4D_short(M,N,8,8);
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					//quantize and round each 8x8 matrix
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							q_y_F[m][n][u][v]=(short)round(y_F[m][n][u][v]/y_Q[u][v]);
							q_cb_F[m][n][u][v]=(short)round(cb_F[m][n][u][v]/cbcr_Q[u][v]);
							q_cr_F[m][n][u][v]=(short)round(cr_F[m][n][u][v]/cbcr_Q[u][v]);
						}
					}
				}
			}
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			
			//DC DPCM
			for(m=0;m<M;m++){
				int n;
				for(n=N-1;n>=1;n--){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]-q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]-q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]-q_cr_F[m][n-1][0][0];
				}
			}
			for(m=M-1;m>=1;m--){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]-q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]-q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]-q_cr_F[m-1][0][0][0];
			} 
			
			
			//zigzag
			int zz_matrix[8][8]={
				{1,2,6,7,15,16,28,29},
				{3,5,8,14,17,27,30,43},
				{4,9,13,18,26,31,42,44},
				{10,12,19,25,32,41,45,54},
				{11,20,24,33,40,46,53,55},
				{21,23,34,39,47,52,56,61},
				{22,35,38,48,51,57,60,62},
				{36,37,49,50,58,59,63,64}
			};
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F_zigzag[m][n][zz_matrix[r][c]-1]=q_y_F[m][n][r][c];
							q_cb_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cb_F[m][n][r][c];
							q_cr_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cr_F[m][n][r][c];
						}
					}
				}
			}
			
			free_4D_short(q_y_F,M,N,8,8);//frequency domain of the 8x8 matrices
			free_4D_short(q_cb_F,M,N,8,8);
			free_4D_short(q_cr_F,M,N,8,8);
			
			short***y_RLE_lists=malloc_3D_short(M,N,64);
			short***cb_RLE_lists=malloc_3D_short(M,N,64);
			short***cr_RLE_lists=malloc_3D_short(M,N,64);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					RLE_8x8(q_y_F_zigzag[m][n],y_RLE_lists[m][n]);
					RLE_8x8(q_cb_F_zigzag[m][n],cb_RLE_lists[m][n]);
					RLE_8x8(q_cr_F_zigzag[m][n],cr_RLE_lists[m][n]);
						
				}
			}
			
			codebook dict=huffmandict(y_RLE_lists,cb_RLE_lists,cr_RLE_lists,M,N,argv[4]);//encoder 3 Kimberly.bmp ascii codebook.txt huffman_code.txt 
			FILE*fptr=fopen(argv[5],"w");
			fptr=save_header_txt(&img,fptr);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					//printf("%d %d\n",m,n);
					fptr=RLE2Huffman_txt(&dict,y_RLE_lists[m][n],fptr);
					fptr=RLE2Huffman_txt(&dict,cb_RLE_lists[m][n],fptr);
					fptr=RLE2Huffman_txt(&dict,cr_RLE_lists[m][n],fptr);
				}
			}
			
			
			
			
			
		}
		
		 
		else if(argv[3][0]=='b'){//encoder 3 Kimberly.bmp binary codebook.txt huffman_code.bin 
			//read image
			Img img=ImRead(argv[2]);
			int W=8;
			int H=8;
			int M=img.header.height/H;//圖片高有幾個8x8 matrix
			int N=img.header.width/W;//圖片寬有幾個8x8 matrix
	
			//bgr to ycbcr
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			bgr2ycbcr(&(img),y_data,cb_data,cr_data);
		
	
			//generate basis vector for 2D-DCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int u;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			for(u=0;u<H;u++){
				int v;
				for(v=0;v<W;v++){
					int r;
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						int c;
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
				
					}
			
				}
			}
			//2D-DCT
			double****y_F=malloc_4D_double(M,N,8,8);//frequency domain of the 8x8 matrices
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			double**beta=ones(H,1);
			beta[0][0]=1/sqrt(2);
			int m,n;
			//Cut the matrix into multiple 8x8 matrices
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
				
					//Get one of the 8x8 matrices
					double y_matrix_8x8[8][8];
					double cb_matrix_8x8[8][8];
					double cr_matrix_8x8[8][8];
					cut_8x8(y_data,y_matrix_8x8,m,n);
					cut_8x8(cb_data,cb_matrix_8x8,m,n);
					cut_8x8(cr_data,cr_matrix_8x8,m,n);
					//subtract 128
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_matrix_8x8[u][v]-=128;
							cb_matrix_8x8[u][v]-=128;
							cr_matrix_8x8[u][v]-=128;
						}
					}
					//Get the frequency domain of the 8x8 matrix with DCT
					for(u=0;u<H;u++){
						int v;
						for(v=0;v<W;v++){
							DCT_8x8(y_matrix_8x8,y_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cb_matrix_8x8,cb_F,beta,basic_vector,m,n,u,v);
							DCT_8x8(cr_matrix_8x8,cr_F,beta,basic_vector,m,n,u,v);
						}
					}
					
				}
			} 
		
		
			free_4D_double(basic_vector,8,8,8,8);
			free_2D_double(y_data,img.header.height,img.header.width);
			free_2D_double(cb_data,img.header.height,img.header.width);
			free_2D_double(cr_data,img.header.height,img.header.width);
			
			
			//quantization
			double y_Q[8][8]={{16,11,10,16,24,40,51,61},
    				 		{12,12,14,19,26,58,60,55},
     						{14,13,16,24,40,57,69,56},
     						{14,17,22,29,51,87,80,62},
     						{18,22,37,56,68,109,103,77},
     						{24,35,55,64,81,104,113,92},
     						{49,64,78,87,103,121,120,101},
     						{72,92,95,98,112,100,103,99}};
     
    			double cbcr_Q[8][8]={ 	{ 17, 18, 24, 47, 99, 99, 99, 99 },
						{ 18, 21, 26, 66, 99, 99, 99, 99 },
						{ 24, 26, 56, 99, 99, 99, 99, 99 },
						{ 47, 66, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99},
						{ 99, 99, 99, 99, 99, 99, 99, 99 },
						{ 99, 99, 99, 99, 99, 99, 99, 99 } };
    		
			short****q_y_F=malloc_4D_short(M,N,8,8);
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					//quantize and round each 8x8 matrix
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							q_y_F[m][n][u][v]=(short)round(y_F[m][n][u][v]/y_Q[u][v]);
							q_cb_F[m][n][u][v]=(short)round(cb_F[m][n][u][v]/cbcr_Q[u][v]);
							q_cr_F[m][n][u][v]=(short)round(cr_F[m][n][u][v]/cbcr_Q[u][v]);
						}
					}
				}
			}
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			
			//DC DPCM
			for(m=0;m<M;m++){
				int n;
				for(n=N-1;n>=1;n--){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]-q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]-q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]-q_cr_F[m][n-1][0][0];
				}
			}
			for(m=M-1;m>=1;m--){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]-q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]-q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]-q_cr_F[m-1][0][0][0];
			} 
			
			
			//zigzag
			int zz_matrix[8][8]={
				{1,2,6,7,15,16,28,29},
				{3,5,8,14,17,27,30,43},
				{4,9,13,18,26,31,42,44},
				{10,12,19,25,32,41,45,54},
				{11,20,24,33,40,46,53,55},
				{21,23,34,39,47,52,56,61},
				{22,35,38,48,51,57,60,62},
				{36,37,49,50,58,59,63,64}
			};
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F_zigzag[m][n][zz_matrix[r][c]-1]=q_y_F[m][n][r][c];
							q_cb_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cb_F[m][n][r][c];
							q_cr_F_zigzag[m][n][zz_matrix[r][c]-1]=q_cr_F[m][n][r][c];
						}
					}
				}
			}
			
			free_4D_short(q_y_F,M,N,8,8);//frequency domain of the 8x8 matrices
			free_4D_short(q_cb_F,M,N,8,8);
			free_4D_short(q_cr_F,M,N,8,8);
			
			short***y_RLE_lists=malloc_3D_short(M,N,64);
			short***cb_RLE_lists=malloc_3D_short(M,N,64);
			short***cr_RLE_lists=malloc_3D_short(M,N,64);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					RLE_8x8(q_y_F_zigzag[m][n],y_RLE_lists[m][n]);
					RLE_8x8(q_cb_F_zigzag[m][n],cb_RLE_lists[m][n]);
					RLE_8x8(q_cr_F_zigzag[m][n],cr_RLE_lists[m][n]);
						
				}
			}
			
			codebook dict=huffmandict(y_RLE_lists,cb_RLE_lists,cr_RLE_lists,M,N,argv[4]);//encoder 3 Kimberly.bmp ascii codebook.txt huffman_code.txt 
			FILE*fptr=fopen(argv[5],"wb");
			fwrite(&img.header.identifier, sizeof(short), 1, fptr);
			fwrite(&img.header.filesize, sizeof(int), 1, fptr);
			fwrite(&img.header.reserved, sizeof(short), 1, fptr);
			fwrite(&img.header.reserved2, sizeof(short), 1, fptr);
			fwrite(&img.header.bitmap_dataoffset, sizeof(int), 1, fptr);
			fwrite(&img.header.bitmap_headersize, sizeof(int), 1, fptr);
			fwrite(&img.header.width, sizeof(int), 1, fptr);
			fwrite(&img.header.height, sizeof(int), 1, fptr);
			fwrite(&img.header.planes, sizeof(short), 1, fptr);
			fwrite(&img.header.bits_perpixel, sizeof(short), 1, fptr);
			fwrite(&img.header.compression, sizeof(int), 1, fptr);
			fwrite(&img.header.bitmap_datasize, sizeof(int), 1, fptr);
			fwrite(&img.header.hresolution, sizeof(int), 1, fptr);
			fwrite(&img.header.vresolution, sizeof(int), 1, fptr);
			fwrite(&img.header.usedcolors, sizeof(int), 1, fptr);
			fwrite(&img.header.importantcolors, sizeof(int), 1, fptr);
			int bin_filesize=0;
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					//printf("%d %d\n",m,n);
					fptr=RLE2Huffman_bin(&dict,y_RLE_lists[m][n],fptr,&bin_filesize);
					fptr=RLE2Huffman_bin(&dict,cb_RLE_lists[m][n],fptr,&bin_filesize);
					fptr=RLE2Huffman_bin(&dict,cr_RLE_lists[m][n],fptr,&bin_filesize);
				}
			}
			//show compress ratio
			printf("compress ratio: %lf",(double)(bin_filesize+54)/img.header.filesize);
			
			
			
			

		}
		
		
		
	}
	return 0;
}

FILE*RLE2Huffman_bin(codebook*dict,short*RLE,FILE*fptr,int*bin_filesize){
	int i1,i2;
	i1=0;
	i2=1;
	while(1){
		short val1=RLE[i1];
		int vali=0;
		for(vali=0;vali<dict->node_num;vali++){
			if(dict->vals[vali]==val1){
				break;
			}
		}
		int c=0;
		unsigned int bin_val1=0;
		while(dict->bitstream[vali][c]!='\0'){
			//fprintf(fptr,"%c",dict->bitstream[vali][c]);
			//fwrite(&dict->bitstream[vali][c],sizeof(char),1,fptr);
			bin_val1=bin_val1*2+(dict->bitstream[vali][c]-'0');
			c++;
		}
		fwrite(&bin_val1,sizeof(unsigned int),1,fptr);
		(*bin_filesize)+=sizeof(unsigned int);
		//char space='\n';
		//fprintf(fptr,"%c",space);
		
		
		
		short val2=RLE[i2];
		for(vali=0;vali<dict->node_num;vali++){
			if(dict->vals[vali]==val2){
				break;
			}
		}
		c=0;
		unsigned int bin_val2=0;
		while(dict->bitstream[vali][c]!='\0'){
			bin_val2=bin_val2*2+(dict->bitstream[vali][c]-'0');
			c++;
		}
		fwrite(&bin_val2,sizeof(unsigned int),1,fptr);
		(*bin_filesize)+=sizeof(unsigned int);
		if(RLE[i1]==0&&RLE[i2]==0){
			break;
		}
		i1=i2+1;
		i2=i1+1;
	}
	return fptr;
}





FILE*RLE2Huffman_txt(codebook*dict,short*RLE,FILE*fptr){
	int i1,i2;
	i1=0;
	i2=1;
	while(1){
		short val1=RLE[i1];
		int vali=0;
		for(vali=0;vali<dict->node_num;vali++){
			if(dict->vals[vali]==val1){
				break;
			}
		}
		int c=0;
		while(dict->bitstream[vali][c]!='\0'){
			fprintf(fptr,"%c",dict->bitstream[vali][c]);
			c++;
		}
		char space='\n';
		fprintf(fptr,"%c",space);
		
		short val2=RLE[i2];
		for(vali=0;vali<dict->node_num;vali++){
			if(dict->vals[vali]==val2){
				break;
			}
		}
		c=0;
		while(dict->bitstream[vali][c]!='\0'){
			fprintf(fptr,"%c",dict->bitstream[vali][c]);
			c++;
		}
		fprintf(fptr,"%c",space);
		if(RLE[i1]==0&&RLE[i2]==0){
			break;
		}
		i1=i2+1;
		i2=i1+1;
	}
	return fptr;
}


codebook huffmandict(short***y,short***cb,short***cr,int M,int N,char*word_filename){
	codebook dict;
	short maxnum,minnum;
	maxnum=y[0][0][0];
	minnum=y[0][0][0];
	//找最大/小值
	int i,j,c;
	for(i=0;i<M;i++){
		for(j=0;j<N;j++){
			for(c=0;c<64;c++){
				if(maxnum<y[i][j][c]){
					maxnum=y[i][j][c];
				}
				if(minnum>y[i][j][c]){
					minnum=y[i][j][c];
				}
				if(c>0&&y[i][j][c-1]==0&&y[i][j][c]==0){
					break;
				}
			}
			for(c=0;c<64;c++){
				if(maxnum<cb[i][j][c]){
					maxnum=cb[i][j][c];
				}
				if(minnum>cb[i][j][c]){
					minnum=cb[i][j][c];
				}
				if(c>0&&cb[i][j][c-1]==0&&cb[i][j][c]==0){
					break;
				}
			}
			for(c=0;c<64;c++){
				if(maxnum<cr[i][j][c]){
					maxnum=cr[i][j][c];
				}
				if(minnum>cr[i][j][c]){
					minnum=cr[i][j][c];
				}
				if(c>0&&cr[i][j][c-1]==0&&cr[i][j][c]==0){
					break;
				}
			}
			
			
		}
	} 
	//printf("max=%d\nmin=%d\n",maxnum,minnum);
	
	//return dict;
	//統計每個數字出現的次數
	int amount=(1+maxnum-minnum);
	//printf("amount=%d\n",amount);
	int*record=(int*)calloc(amount,sizeof(int));
	for(i=0;i<M;i++){
		for(j=0;j<N;j++){
			for(c=0;c<64;c++){
				record[y[i][j][c]-minnum]++;
				
				if(c>0&&y[i][j][c-1]==0&&y[i][j][c]==0){
					break;
				}
			}
			for(c=0;c<64;c++){
				record[cb[i][j][c]-minnum]++;
				
				if(c>0&&cb[i][j][c-1]==0&&cb[i][j][c]==0){
					break;
				}
			}
			for(c=0;c<64;c++){
				record[cr[i][j][c]-minnum]++;
				if(c>0&&cr[i][j][c-1]==0&&cr[i][j][c]==0){
					break;
				}
			}
		}
	}
	//放入queue
	Huffman_queue queue=Huffman_queue_init(); 
	int node_num=0;
	int num;
	for(num=0;num<amount;num++){
		if(record[num]>0){
			Huffman_node*node=(Huffman_node*)malloc(sizeof(Huffman_node));
			node->left=NULL;
			node->right=NULL;
			node->count=record[num];
			node->val=num+minnum;
			Huffman_queue_push(&queue,node);
			node_num++;
		}
		
	}
	//return dict;
	
	int tmp_num=node_num;
	//FILE*debug_fptr=fopen("debug_encoder.txt","w");
	while(tmp_num>1){
		//print_queue(&queue);
		Huffman_node*ptr_most=queue.front->next;//找最小和第二小的位置
		Huffman_node*ptr=queue.front->next;
		while(ptr!=NULL){
			if((ptr_most->count)>(ptr->count)){
				ptr_most=ptr;
			}
			
			ptr=ptr->next;
		}
		ptr=queue.front->next;
		Huffman_node*ptr_second=NULL;
		while(ptr!=NULL){
			if(ptr!=ptr_most&&(ptr_second==NULL||((ptr->count)<(ptr_second->count)))){
				ptr_second=ptr;
			}
			ptr=ptr->next;
		}
		/*if(ptr_most==ptr_second){
			printf("find most and second error\n");
		}*/
		Huffman_node*node=(Huffman_node*)malloc(sizeof(Huffman_node));
		node->right=ptr_most;
		node->left=ptr_second;
		/*if(node->left==NULL||node->right==NULL){
			printf("error in queue\n");
		}*/
		
		node->count=ptr_most->count+ptr_second->count;
		node->next=NULL;
		Huffman_queue_pop(&queue,ptr_most);
		Huffman_queue_pop(&queue,ptr_second);
		Huffman_queue_push(&queue,node);
		tmp_num--;
	}
	Huffman_node*head=queue.front->next;
	
	FILE*word_fptr=fopen(word_filename,"w");
	int val_ptr=0;
	dict.bitstream=(char**)malloc(node_num*sizeof(char*));
	dict.vals=(short*)malloc(node_num*sizeof(short));
	char*buff=(char*)malloc(node_num*sizeof(char));
	save_codebook_txt(&word_fptr,dict.vals,&val_ptr,dict.bitstream,head,buff,-1);
	dict.node_num=node_num;
	return dict;
}

void save_codebook_txt(FILE**word_fptr,short*vals,int*val_ptr,char**bitstream,Huffman_node*node,char*buff,int buff_ptr){
	if(node==NULL){
		printf("error in node\n");
		return;	
	}
	if(node->left==NULL&&node->right==NULL){
		fprintf(*word_fptr,"%d ",node->val);//symbol
		fprintf(*word_fptr,"%d ",node->count);
		int i=0;
		//printf("%d %d ",node->val,node->count);//debug
		int word_count=0;
		char*word_buff=(char*)malloc(100*sizeof(char));
		int bit_amount=0;
		while(i<=buff_ptr){//convert a bitstream which come from Huffman tree to many codewords
			int res=0;
			int count=1;
			while(i<=buff_ptr&&count<=8){
				//fprintf(debug_fptr,"%c",buff[i]);//debug
				res=res*2+(buff[i]-'0');
				bit_amount++;
				i++;
				count++;
			}
			//if(res<0||res>=256)printf("character error\n");
			word_count++;
			char word=(char)res;
			//printf("%d ",res);//debug
			word_buff[word_count-1]=word;
			//fprintf(*word_fptr,"%c",word);
		}
		fprintf(*word_fptr,"%d ",word_count);
		fprintf(*word_fptr,"%d ",bit_amount);
		for(i=0;i<word_count;i++){
			fprintf(*word_fptr,"%c",word_buff[i]);
		}
		
		fprintf(*word_fptr,"%c",'\n');//record bits stream correspond to value of node 
		vals[*(val_ptr)]=node->val;
		bitstream[*(val_ptr)]=(char*)malloc((buff_ptr+1+1)*sizeof(char));
		i=0;
		while(i<=buff_ptr){
			bitstream[*(val_ptr)][i]=buff[i];
			i++;
		}
		bitstream[*(val_ptr)][i]='\0';//represent end of bit stream
		//printf("%s\n",bitstream[*(val_ptr)]);
		(*val_ptr)++;
	}
	else{
		buff[buff_ptr+1]='0';
		save_codebook_txt(word_fptr,vals,val_ptr,bitstream,node->left,buff,buff_ptr+1);
		buff[buff_ptr+1]='1';
		save_codebook_txt(word_fptr,vals,val_ptr,bitstream,node->right,buff,buff_ptr+1);
	}
}





void show_channel_compress_ratio(int M,int N,short***RLE_lists,short***y_F_zigzag,int*RLE_size,int*zigzag_size){//compare size before RLE and size after RLE
	int m,n;
	*(RLE_size)=0;
	*(zigzag_size)=M*N*64*sizeof(short);
	for(m=0;m<M;m++){
		for(n=0;n<N;n++){
			//printf("%d %d\n",m,n);
			int ptr1,ptr2;
			ptr1=0;
			ptr2=1;
			while(1){
				(*RLE_size)+=sizeof(short);
				if(RLE_lists[m][n][ptr1]==0&&RLE_lists[m][n][ptr2]==0){
					break;
				}
				ptr1=ptr2+1;
				ptr2=ptr1+1;
			}
				
		}
	}

	printf("compress ratio: %lf\n",(double)(*RLE_size)/(*zigzag_size));
	
}



void SQNR(double****F,double****e,int H,int W){
	double p=0;
	double p0=0;
	double sqnr=0;
	int i;
	int j;
	for(i=0;i<H;i++){
		for(j=0;j<W;j++){
		int m,n;
			for(m=0;m<8;m++){
				for(n=0;n<8;n++){
					p+=(double)pow(F[i][j][m][n],2);
					p0+=(double)pow(e[i][j][m][n],2);
				}		
			}
			sqnr=10*log10(p/p0);
			printf("sqnr=%lf ",sqnr);
			p=0;
			p0=0;
		}
		printf("\n");
	}
}
void print_queue(Huffman_queue*queue){
	Huffman_node*node=queue->front->next;
	while(node!=NULL){
		printf("%d ",node->count);
		node=node->next;
	}
	printf("\n");
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

FILE*save_8x8_RLE_bin(FILE*fptr,short*list){
	/*RLE_node*ptr=list->next->next;
	RLE_node*skip=list->next;*/
	int skip=0;
	int ptr=1;
	while(1){
		fwrite(&list[skip],sizeof(short),1,fptr);
		fwrite(&list[ptr],sizeof(short),1,fptr);
		if(list[skip]==0&&list[ptr]==0){
			break;
		}
		skip=ptr+1;
		ptr=skip+1;
	}
	//printf("\n");
	return fptr;
}



FILE*save_8x8_RLE_txt(FILE*fptr,short*list){
	int i1=0;
	int i2=1;
	while(1){
		fprintf(fptr,"%d %d ",list[i1],list[i2]);
		if(list[i1]==0&&list[i2]==0){
			break;
		}
		i1=i2+1;
		i2=i1+1;
	}
	fprintf(fptr,"%c",'\n');
	return fptr;
}


FILE*save_header_txt(Img*img,FILE*fp){
	//寫入header進dim.txt
	fprintf(fp,"%c%c ",img->header.identifier[0],img->header.identifier[1]);
	fprintf(fp,"%d ",img->header.filesize);
	fprintf(fp,"%d ",img->header.reserved);
	fprintf(fp,"%d ",img->header.reserved2);
	fprintf(fp,"%d ",img->header.bitmap_dataoffset);
	fprintf(fp,"%d ",img->header.bitmap_headersize);
	fprintf(fp,"%d ",img->header.width);
	fprintf(fp,"%d ",img->header.height);
	fprintf(fp,"%d ",img->header.planes);
	fprintf(fp,"%d ",img->header.bits_perpixel);
	fprintf(fp,"%d ",img->header.compression);
	fprintf(fp,"%d ",img->header.bitmap_datasize);
	fprintf(fp,"%d ",img->header.hresolution);
	fprintf(fp,"%d ",img->header.vresolution);
	fprintf(fp,"%d ",img->header.usedcolors);
	fprintf(fp,"%d\n",img->header.importantcolors);
	return fp;
}



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

unsigned char**malloc_2D_uchar(int row,int col){
	unsigned char**res=(unsigned char**)malloc(row*sizeof(unsigned char*));
	int i;
	for(i=0;i<row;i++){
		res[i]=(unsigned char*)malloc(col*sizeof(unsigned char));
	}	
	return res;
}
double**malloc_2D_double(int row_size,int col_size){
	int i,j;
	double**res=(double**)malloc(row_size*sizeof(double*));
	for(i=0;i<row_size;i++){
		res[i]=(double*)malloc(col_size*sizeof(double));
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
double****malloc_4D_double(int size1,int size2,int size3,int size4){
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
short***malloc_3D_short(int size1,int size2,int size3){
	short***res=(short***)malloc(size1*sizeof(short**));
	int i,j;
	for(i=0;i<size1;i++){
		res[i]=(short**)malloc(size2*sizeof(short*));
		for(j=0;j<size2;j++){
			res[i][j]=(short*)malloc(size3*sizeof(short));
		}
	}
	return res;
}

short****malloc_4D_short(int size1,int size2,int size3,int size4){
	short****res=(short****)malloc(size1*sizeof(short***));
	int i1;
	for(i1=0;i1<size1;i1++){
		res[i1]=(short***)malloc(size2*sizeof(short**));
		int i2;
		for(i2=0;i2<size2;i2++){
			res[i1][i2]=(short**)malloc(size3*sizeof(short*));
			int i3;
			for(i3=0;i3<size3;i3++){
				res[i1][i2][i3]=(short*)malloc(size4*sizeof(short));
			}
		}
	}
	return res;
}
RLE_node***malloc_3D_RLE_node(int M,int N){
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
void RLE_8x8(short*F_zigzag,short*RLE){
	int i=0;
	short count_zero=0;
	int ptr=0;
	while(i<64){
		if(F_zigzag[i]!=0){
			
			RLE[ptr]=count_zero;
			RLE[ptr+1]=F_zigzag[i];//前面的非0值到這個非0值有幾個0 
			count_zero=0;
			ptr+=2;
		}
		else{
			count_zero++;
		}
		i++;
	}
	RLE[ptr]=0;
	RLE[ptr+1]=0;
}

RLE_node*new_RLE_node(short val){
	RLE_node*res=(RLE_node*)malloc(sizeof(RLE_node));
	res->val=val;
	res->next=NULL;
	return res;
}


void add_RLE_node(RLE_node*last,RLE_node*newnode){
	last->next=newnode;
}
int RLE_code_len(RLE_node*RLE_code){ 
	RLE_node*ptr=RLE_code->next;
	long long int len=0;
	while(ptr!=NULL){//計算每個RLE list的長度 
		len++;
		ptr=ptr->next;	
	}
	return len;	
}
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

void saveas(unsigned char**matrix,int row,int col,char*filename){
	int i,j;
	FILE *fptr;
    fptr = fopen(filename,"w+");
    //printf("save %s",filename);
    unsigned int tmp;
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			//printf("%x",matrix[i][j]);
			//fputc(matrix[i][j],fptr);
			tmp=matrix[i][j];
			fprintf(fptr,"%u ",tmp);
		}
		fprintf(fptr,"%c",'\n');
		//fputc('\n',fptr);
		//printf("\n");
	}
	fclose(fptr);
	//printf("\n");
}

void savetxt_table_8x8(int arr[8][8],char*filename){
	FILE *fptr;
    fptr = fopen(filename,"w+");
	int i,j;
	for(i=0;i<8;i++){
		for(j=0;j<8;j++){
			fprintf(fptr,"%d ",arr[i][j]);
		}	
		fprintf(fptr,"%c",'\n');
	}
}

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
	res.pixels=pix_data;
	return res;
}
void ImWrite(Img img,char filename[]){
	FILE *outfile;
    outfile= fopen(filename, "wb");
    reverse_row(img.pixels,img.header.height,img.header.width);
    
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
			fwrite(img.pixels[i][j].data,bytes_perpixel,1,outfile);
		}
		if(skip!=0){
			fwrite(skip_buf, skip, 1, outfile);
		}
	}
}
void bgr2ycbcr(Img*img,double**y_data,double**cb_data,double**cr_data){
	int i,j;
	for(i=0;i<(img->header).height;i++){
		for(j=0;j<(img->header).width;j++){
			unsigned int B=(img->pixels[i][j]).data[0];
			unsigned int G=(img->pixels[i][j]).data[1];
			unsigned int R=(img->pixels[i][j]).data[2];
			//if(B>255||G>255||R>255)printf("overflow\n");
			y_data[i][j]=0.299*R + 0.587*G + 0.114*B;
			cb_data[i][j]=0.5643*(B - y_data[i][j]);
			cr_data[i][j]=0.7133*(R - y_data[i][j]);
		}
	}
}
void cut_8x8(double**data,double matrix_8x8[8][8],int m,int n){
	int i;//byte_index為某個像素的第幾個byte 
	for(i=0;i<8;i++){
		int j;
		for(j=0;j<8;j++){
			matrix_8x8[i][j]=data[8*m+i][8*n+j];//(double)data[8*m+i][8*n+j].data[byte_index];
		}
	}
}
void DCT_8x8(double matrix_8x8[8][8],double****F,double**beta,double ****basic_vector,int m,int n,int u,int v){
	int r;
	F[m][n][u][v]=0;
	for(r=0;r<8;r++){
		int c;
		for(c=0;c<8;c++){
			F[m][n][u][v]+=2.000000/8*beta[u][0]*beta[v][0]*matrix_8x8[r][c]*basic_vector[u][v][r][c];
		}
	}
}


void free_4D_short(short****data,int size1,int size2,int size3,int size4){
	int i1,i2,i3;
	for(i1=0;i1<size1;i1++){
		for(i2=0;i2<size2;i2++){
			for(i3=0;i3<size3;i3++){
				free(data[i1][i2][i3]);
			}
			free(data[i1][i2]);
		}
		free(data[i1]);
	}
	free(data);
}
void free_2D_double(double**data,int size1,int size2){
	int i1=0;
	for(i1=0;i1<size1;i1++){
		free(data[i1]);
	}
	free(data);
}
void free_4D_double(double****data,int size1,int size2,int size3,int size4){
	int i1,i2,i3;
	for(i1=0;i1<size1;i1++){
		for(i2=0;i2<size2;i2++){
			for(i3=0;i3<size3;i3++){
				free(data[i1][i2][i3]);
			}
			free(data[i1][i2]);
		}
		free(data[i1]);
	}
	free(data);
}
void free_3D_RLE_node(RLE_node***data,int size1,int size2,int size3){
	int i1,i2;
	for(i1=0;i1<size1;i1++){
		for(i2=0;i2<size2;i2++){
			free(data[i1][i2]);
		}
		free(data[i1]);
	}
	free(data);
}
void free_3D_short(short***data,int size1,int size2,int size3){
	int i1,i2;
	for(i1=0;i1<size1;i1++){
		for(i2=0;i2<size2;i2++){
			free(data[i1][i2]);
		}
		free(data[i1]);
	}
	free(data);
}
