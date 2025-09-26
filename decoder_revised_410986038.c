#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdint.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
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
typedef struct codebook{
	short*vals;
	char**bitstream;
	int node_num;
}codebook;
unsigned char**malloc_2D_uchar(int row,int col);
void print_header(Img img);
void ImWrite(Img img,char filename[]);
Pixel** malloc_2D_Pixel(int row, int col,unsigned short bytes_perpixel);
double**malloc_2D_double(int row_size,int col_size);
void free_4D_double(double****data,int size1,int size2,int size3,int size4);
short****malloc_4D_short(int size1,int size2,int size3,int size4);
RLE_node***malloc_3D_RLE_node(int M,int N);
short***malloc_3D_short(int size1,int size2,int size3);

void RLE_8x8(short*F_zigzag,RLE_node*list);
RLE_node*new_RLE_node(short val);
void add_RLE_node(RLE_node*last,RLE_node*newnode);

double**ones(int row,int col);
void reverse_row(Pixel**arr,int rowSize,int colSize);
void readheader(FILE* fp, Bitmap *x);
void readHeader_txt(Img*img,char*filename);
void read_table_txt(int table[8][8],char*filename);
double****malloc_4D_double(int size1,int size2,int size3,int size4);
void IDCT_8x8(double**F,double**matrix_8x8,double**beta,double****basic_vector,int r,int c);
void ycbcr2bgr(Img*img,double**y_data,double**cb_data,double**cr_data);
Img ImRead(char input_name[]);
FILE*read_8x8_RLE_txt(FILE*fptr,short*list);
FILE*read_8x8_RLE_bin(FILE*fptr,short*list);
void un_RLE_8x8(short*list,short*F_zigzag);
FILE*print_RLE_lists(short*RLE_list,FILE*fptr);
int s2i(char*s,int*ptr1,int*ptr2);
codebook read_codebook_txt(char*filename);
short decode_huffman(codebook*book,char*data);
short decode_huffman_by_val(codebook*book,unsigned int*bin_val,unsigned int word_val);
void char2bin(unsigned char c,char*bin);
int main(int argc,char**argv){
	if(argv[1][0]=='0'){
		FILE*rptr=fopen(argv[3],"r");
		FILE*gptr=fopen(argv[4],"r");
		FILE*bptr=fopen(argv[5],"r");
		Img img;
		readHeader_txt(&img,argv[6]);
		int i,j;
		img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);
		unsigned char c;
		for(i=0;i<img.header.height;i++){
			for(j=0;j<img.header.width;j++){
				fscanf(bptr,"%u",&img.pixels[i][j].data[0]);
				fscanf(gptr,"%u",&img.pixels[i][j].data[1]);
				fscanf(rptr,"%u",&img.pixels[i][j].data[2]);
			}
		}
		ImWrite(img,argv[2]);
	}
	
	
	if(argv[1][0]=='1'){
		if(argc==11){
			//decoder 1 QResKimberly.bmp Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw


			
			int y_Q[8][8];
			int cbcr_Q[8][8];
			read_table_txt(y_Q,argv[4]);
			read_table_txt(cbcr_Q,argv[5]);
			Img img;
			//get header
			readHeader_txt(&img,argv[7]);
			/*print_header(img);//debug
			return 0;*/
			
			
			img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);
			int M=img.header.height/8;
			int N=img.header.width/8;
			int H=8;
			int W=8;
			//get frequency quantized from txt file
			FILE*q_y_fptr=fopen(argv[8],"rb");
			FILE*q_cb_fptr=fopen(argv[9],"rb");
			FILE*q_cr_fptr=fopen(argv[10],"rb");
			
			short****q_y_F=malloc_4D_short(M,N,8,8);
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			
			//從紀錄量化值的raw檔讀取下來 
			
			int m,n,u,v;
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					for(u=0;u<8;u++){
						for(v=0;v<8;v++){
							fread(&q_y_F[m][n][u][v],sizeof(short),1,q_y_fptr);
							fread(&q_cb_F[m][n][u][v],sizeof(short),1,q_cb_fptr);
							fread(&q_cr_F[m][n][u][v],sizeof(short),1,q_cr_fptr);
							
						}
					}
					
				}
			}
			
			
			double****y_F=malloc_4D_double(M,N,8,8);
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			
			
			
			//unquantize
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					for(u=0;u<8;u++){
						for(v=0;v<8;v++){
							y_F[m][n][u][v]=(float)(q_y_F[m][n][u][v]*y_Q[u][v]);
							cb_F[m][n][u][v]=(float)(q_cb_F[m][n][u][v]*cbcr_Q[u][v]);
							cr_F[m][n][u][v]=(float)(q_cr_F[m][n][u][v]*cbcr_Q[u][v]);
							
						}
					}
				}
			}
			
			
			
			//generate basis vector for 2D-DCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int r,c;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			for(u=0;u<H;u++){
				for(v=0;v<W;v++){
					//Make a 8x8 basic vector matrix x
					double x[8][8];
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
							x[r][c]=127.5+127.5*basic_vector[u][v][r][c];
						}
					}
				}
			}
			
			//做IDCT 
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			double**beta=ones(H,1);
			beta[0][0]=1.000000/sqrt(2);
			Pixel**res_data=malloc_2D_Pixel(img.header.height,img.header.width,3);
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							
							IDCT_8x8(y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
							y_tmp[r][c]+=128.000000;
							cb_tmp[r][c]+=128.000000;
							cr_tmp[r][c]+=128.000000;
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
			
			
			
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			//show R/G/B SQNR 
			Img o_img=ImRead(argv[3]);
			int i,j;
			double B_SQNR,G_SQNR,R_SQNR;
			double power,e_power;
			int color;
			for(color=0;color<3;color++){
				power=e_power=0;
				for(i=0;i<o_img.header.height;i++){
					for(j=0;j<o_img.header.width;j++){
						power+=(double)pow(o_img.pixels[i][j].data[color],2);
						e_power+=pow(o_img.pixels[i][j].data[color]-img.pixels[i][j].data[color],2);
					}
				}
				double sqnr=10*log10(power/e_power);
				if(color==0)printf("Blue ");
				else if(color==1)printf("Green ");
				else printf("Red ");
				printf("SQNR=%lf\n",sqnr);
			}
			
			ImWrite(img,argv[2]);
			
			
			return 0;
		
		}
		
		else if(argc==13){
			//decoder 1 ResKimberly.bmp 3~5Qtable dim.txt 7~9qF 10~12qe	
			
			int y_Q[8][8];
			int cbcr_Q[8][8];
			read_table_txt(y_Q,argv[3]);
			read_table_txt(cbcr_Q,argv[4]);
			Img img;
			//get header
			readHeader_txt(&img,argv[6]);
			img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);
			int M=img.header.height/8;
			int N=img.header.width/8;
			int H=8;
			int W=8;
			//get frequency quantized from txt file
			FILE*q_y_fptr=fopen(argv[7],"rb");
			FILE*q_cb_fptr=fopen(argv[8],"rb");
			FILE*q_cr_fptr=fopen(argv[9],"rb");
			
			short****q_y_F=malloc_4D_short(M,N,8,8);
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			
			//從紀錄量化值的txt檔讀取下來 
			
			int m,n,u,v;
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					for(u=0;u<8;u++){
						for(v=0;v<8;v++){
							fread(&q_y_F[m][n][u][v],sizeof(short),1,q_y_fptr);
							fread(&q_cb_F[m][n][u][v],sizeof(short),1,q_cb_fptr);
							fread(&q_cr_F[m][n][u][v],sizeof(short),1,q_cr_fptr);
							
						}
					}
					
				}
			}
			
			
			double****y_F=malloc_4D_double(M,N,8,8);
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			
			
			
			//unquantize,讀取紀錄誤差的txt檔並和頻率相加 
			FILE*e_y_fptr=fopen(argv[10],"rb");
			FILE*e_cb_fptr=fopen(argv[11],"rb");
			FILE*e_cr_fptr=fopen(argv[12],"rb");
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					for(u=0;u<8;u++){
						for(v=0;v<8;v++){
							float y_e,cb_e,cr_e;
							fread(&y_e,sizeof(float),1,e_y_fptr);
							fread(&cb_e,sizeof(float),1,e_cb_fptr);
							fread(&cr_e,sizeof(float),1,e_cr_fptr);
							
							y_F[m][n][u][v]=y_e+(float)(q_y_F[m][n][u][v]*y_Q[u][v]);
							cb_F[m][n][u][v]=cb_e+(float)(q_cb_F[m][n][u][v]*cbcr_Q[u][v]);
							cr_F[m][n][u][v]=cr_e+(float)(q_cr_F[m][n][u][v]*cbcr_Q[u][v]);
							
						}
					}
				}
			}
			
			
			//generate basis vector for 2D-DCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int r,c;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			for(u=0;u<H;u++){
				for(v=0;v<W;v++){
					//Make a 8x8 basic vector matrix x
					double x[8][8];
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
							x[r][c]=127.5+127.5*basic_vector[u][v][r][c];
						}
					}
				}
			}
			
			//做IDCT 
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			double**beta=ones(H,1);
			beta[0][0]=1.000000/sqrt(2);
			//Pixel**res_data=malloc_2D_Pixel(img.header.height,img.header.width,3);
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							
							IDCT_8x8(y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
							y_tmp[r][c]+=128.000000;
							cb_tmp[r][c]+=128.000000;
							cr_tmp[r][c]+=128.000000;
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
			
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			ImWrite(img,argv[2]);
			
			
			return 0;
		
		}
	}
	if(argv[1][0]=='2'){
		if(argv[3][0]=='a'){//decoder 2 QResKimberly.bmp ascii rle_code.txt 
			FILE*fptr=fopen(argv[4],"r");
			//read header
			Img img;
			fscanf(fptr,"%c%c ",&(img.header.identifier[0]),&(img.header.identifier[1]));
			fscanf(fptr,"%d",&(img.header.filesize));
			fscanf(fptr,"%d",&(img.header.reserved));
			fscanf(fptr,"%d",&(img.header.reserved2));
			fscanf(fptr,"%d",&(img.header.bitmap_dataoffset));
			fscanf(fptr,"%d",&(img.header.bitmap_headersize));
			fscanf(fptr,"%d",&(img.header.width));
			fscanf(fptr,"%d",&(img.header.height));
			fscanf(fptr,"%d",&(img.header.planes));
			fscanf(fptr,"%d",&(img.header.bits_perpixel));
			fscanf(fptr,"%d",&(img.header.compression));
			fscanf(fptr,"%d",&(img.header.bitmap_datasize));
			fscanf(fptr,"%d",&(img.header.hresolution));
			fscanf(fptr,"%d",&(img.header.vresolution));
			fscanf(fptr,"%d",&(img.header.usedcolors));
			fscanf(fptr,"%d",&(img.header.importantcolors));
			img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);
			int M=img.header.height/8;
			int N=img.header.width/8;
			int H=8;
			int W=8;
			//FILE*debug_fptr=fopen("debug_decoder.txt","w");
			short***y_RLE_lists=malloc_3D_short(M,N,64); 
			short***cb_RLE_lists=malloc_3D_short(M,N,64); 
			short***cr_RLE_lists=malloc_3D_short(M,N,64); 
			int m;
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					fptr=read_8x8_RLE_txt(fptr,y_RLE_lists[m][n]);
					//debug_fptr=print_RLE_lists(y_RLE_lists[m][n],debug_fptr);//debug
					fptr=read_8x8_RLE_txt(fptr,cb_RLE_lists[m][n]);
					//debug_fptr=print_RLE_lists(cb_RLE_lists[m][n],debug_fptr);//debug
					fptr=read_8x8_RLE_txt(fptr,cr_RLE_lists[m][n]);
					//debug_fptr=print_RLE_lists(cr_RLE_lists[m][n],debug_fptr);//debug
					
				}
			}
			//un RLE
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					un_RLE_8x8(y_RLE_lists[m][n],q_y_F_zigzag[m][n]);
					un_RLE_8x8(cb_RLE_lists[m][n],q_cb_F_zigzag[m][n]);
					un_RLE_8x8(cr_RLE_lists[m][n],q_cr_F_zigzag[m][n]);
				}
			}

			short****q_y_F=malloc_4D_short(M,N,8,8);//frequency domain of the 8x8 matrices
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			/*free_3D_short(y_RLE_lists);
			free_3D_short(cb_RLE_lists);
			free_3D_short(cr_RLE_lists);*/
			
			//unzigzag
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
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F[m][n][r][c]=q_y_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cb_F[m][n][r][c]=q_cb_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cr_F[m][n][r][c]=q_cr_F_zigzag[m][n][zz_matrix[r][c]-1];
						}
					}
				}
			}
			
			
			//un DPCM
			for(m=1;m<M;m++){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]+q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]+q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]+q_cr_F[m-1][0][0][0];
			}
			for(m=0;m<M;m++){
				int n;
				for(n=1;n<N;n++){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]+q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]+q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]+q_cr_F[m][n-1][0][0];
				}
			}
			
			
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
			
			
			
			
			
			double****y_F=malloc_4D_double(M,N,8,8);
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			
			
			
			//unquantize
			//FILE*debug_fptr=fopen("debug_decoder.txt","w");//debug
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_F[m][n][u][v]=(float)(q_y_F[m][n][u][v]*y_Q[u][v]);
							cb_F[m][n][u][v]=(float)(q_cb_F[m][n][u][v]*cbcr_Q[u][v]);
							cr_F[m][n][u][v]=(float)(q_cr_F[m][n][u][v]*cbcr_Q[u][v]);
							//fprintf(debug_fptr,"%lf %lf %lf\n",y_F[m][n][u][v],cb_F[m][n][u][v],cr_F[m][n][u][v]);//debug
							
						}
					}
				}
			}
			/*fclose(debug_fptr);
			return 0;*/
			
			//generate basis vector for 2D-IDCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int r,c;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			int u,v;
			for(u=0;u<H;u++){
				for(v=0;v<W;v++){
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
					}
				}
			}
			
			//做IDCT 
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			double**beta=ones(H,1);
			beta[0][0]=1.000000/sqrt(2);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							
							IDCT_8x8(y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
							y_tmp[r][c]+=128.000000;
							cb_tmp[r][c]+=128.000000;
							cr_tmp[r][c]+=128.000000;
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
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			ImWrite(img,argv[2]);
			return 0;
		}
		
		else if(argv[3][0]=='b'){
			//decoder 2 QResKimberly.bmp binary rle_code.bin 
			FILE*fptr=fopen(argv[4],"rb");
			Img img;
			fread(img.header.identifier,2, 1, fptr);
			fread(&img.header.filesize, sizeof(int), 1, fptr);
			fread(&img.header.reserved, sizeof(short), 1, fptr);
			fread(&img.header.reserved2, sizeof(short), 1, fptr);
			fread(&img.header.bitmap_dataoffset, sizeof(int), 1, fptr);
			fread(&img.header.bitmap_headersize, sizeof(int), 1, fptr);
			fread(&img.header.width, sizeof(int), 1, fptr);
			fread(&img.header.height, sizeof(int), 1, fptr);
			fread(&img.header.planes, sizeof(short), 1, fptr);
			fread(&img.header.bits_perpixel, sizeof(short), 1, fptr);
			fread(&img.header.compression, sizeof(int), 1, fptr);
			fread(&img.header.bitmap_datasize, sizeof(int), 1, fptr);
			fread(&img.header.hresolution, sizeof(int), 1, fptr);
			fread(&img.header.vresolution, sizeof(int), 1, fptr);
			fread(&img.header.usedcolors, sizeof(int), 1, fptr);
			fread(&img.header.importantcolors, sizeof(int), 1, fptr);
			int M=img.header.height/8;
			int N=img.header.width/8;
			short***y_RLE_lists=malloc_3D_short(M,N,64); 
			short***cb_RLE_lists=malloc_3D_short(M,N,64); 
			short***cr_RLE_lists=malloc_3D_short(M,N,64); 
			int m;
			int H=8;
			int W=8;
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					fptr=read_8x8_RLE_bin(fptr,y_RLE_lists[m][n]);
					fptr=read_8x8_RLE_bin(fptr,cb_RLE_lists[m][n]);
					fptr=read_8x8_RLE_bin(fptr,cr_RLE_lists[m][n]);
					
					
				}
			}
		
			img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);

			//un RLE
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					un_RLE_8x8(y_RLE_lists[m][n],q_y_F_zigzag[m][n]);
					un_RLE_8x8(cb_RLE_lists[m][n],q_cb_F_zigzag[m][n]);
					un_RLE_8x8(cr_RLE_lists[m][n],q_cr_F_zigzag[m][n]);
				}
			}

			short****q_y_F=malloc_4D_short(M,N,8,8);//frequency domain of the 8x8 matrices
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			/*free_3D_short(y_RLE_lists);
			free_3D_short(cb_RLE_lists);
			free_3D_short(cr_RLE_lists);*/
			
			//unzigzag
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
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F[m][n][r][c]=q_y_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cb_F[m][n][r][c]=q_cb_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cr_F[m][n][r][c]=q_cr_F_zigzag[m][n][zz_matrix[r][c]-1];
						}
					}
				}
			}
			
			
			//un DPCM
			for(m=1;m<M;m++){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]+q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]+q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]+q_cr_F[m-1][0][0][0];
			}
			for(m=0;m<M;m++){
				int n;
				for(n=1;n<N;n++){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]+q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]+q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]+q_cr_F[m][n-1][0][0];
				}
			}
			
			
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
			
			
			
			
			
			double****y_F=malloc_4D_double(M,N,8,8);
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			
			
			
			//unquantize
			//FILE*debug_fptr=fopen("debug_decoder.txt","w");//debug
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_F[m][n][u][v]=(float)(q_y_F[m][n][u][v]*y_Q[u][v]);
							cb_F[m][n][u][v]=(float)(q_cb_F[m][n][u][v]*cbcr_Q[u][v]);
							cr_F[m][n][u][v]=(float)(q_cr_F[m][n][u][v]*cbcr_Q[u][v]);
							//fprintf(debug_fptr,"%lf %lf %lf\n",y_F[m][n][u][v],cb_F[m][n][u][v],cr_F[m][n][u][v]);//debug
							
						}
					}
				}
			}
			/*fclose(debug_fptr);
			return 0;*/
			
			//generate basis vector for 2D-IDCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int r,c;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			int u,v;
			for(u=0;u<H;u++){
				for(v=0;v<W;v++){
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
					}
				}
			}
			
			//做IDCT 
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			double**beta=ones(H,1);
			beta[0][0]=1.000000/sqrt(2);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							
							IDCT_8x8(y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
							y_tmp[r][c]+=128.000000;
							cb_tmp[r][c]+=128.000000;
							cr_tmp[r][c]+=128.000000;
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
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			
			ImWrite(img,argv[2]);
			return 0;
			
			
		}
	}
	if(argv[1][0]=='3'){
		
		if(argv[3][0]=='a'){//decoder 3 QResKimberly.bmp ascii codebook.txt huffman_code.txt 
			//read codebook
			codebook book=read_codebook_txt(argv[4]);
			
			//read header
			FILE*fptr=fopen(argv[5],"r");
			Img img;
			fscanf(fptr,"%c%c ",&(img.header.identifier[0]),&(img.header.identifier[1]));
			fscanf(fptr,"%d",&(img.header.filesize));
			fscanf(fptr,"%d",&(img.header.reserved));
			fscanf(fptr,"%d",&(img.header.reserved2));
			fscanf(fptr,"%d",&(img.header.bitmap_dataoffset));
			fscanf(fptr,"%d",&(img.header.bitmap_headersize));
			fscanf(fptr,"%d",&(img.header.width));
			fscanf(fptr,"%d",&(img.header.height));
			fscanf(fptr,"%d",&(img.header.planes));
			fscanf(fptr,"%d",&(img.header.bits_perpixel));
			fscanf(fptr,"%d",&(img.header.compression));
			fscanf(fptr,"%d",&(img.header.bitmap_datasize));
			fscanf(fptr,"%d",&(img.header.hresolution));
			fscanf(fptr,"%d",&(img.header.vresolution));
			fscanf(fptr,"%d",&(img.header.usedcolors));
			fscanf(fptr,"%d\n",&(img.header.importantcolors));
			img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);
			int M=img.header.height/8;
			int N=img.header.width/8;
			int H=8;
			int W=8;
			//print_header(img);
			
			//read bitstream and decode
			int m;
			int n;
			short***y_RLE_lists=malloc_3D_short(M,N,64); 
			short***cb_RLE_lists=malloc_3D_short(M,N,64); 
			short***cr_RLE_lists=malloc_3D_short(M,N,64); 
			char*buff=(char*)malloc(200*sizeof(char));
			short num1,num2;
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					int ptr1,ptr2;
					ptr1=0;
					ptr2=1;
					while(1){
						fscanf(fptr,"%s\n",buff);
						num1=decode_huffman(&book,buff);
						fscanf(fptr,"%s\n",buff);
						num2=decode_huffman(&book,buff);
						//printf("%d %d ",num1,num2);
						y_RLE_lists[m][n][ptr1]=num1;
						y_RLE_lists[m][n][ptr2]=num2;
						if(num1==0&&num2==0){
							break;
						}
						ptr1=ptr2+1;
						ptr2=ptr1+1;
						
					}
					//printf("\n");
					//return 0;
					ptr1=0;
					ptr2=1;
					while(1){
						fscanf(fptr,"%s\n",buff);
						num1=decode_huffman(&book,buff);
						fscanf(fptr,"%s\n",buff);
						num2=decode_huffman(&book,buff);
						cb_RLE_lists[m][n][ptr1]=num1;
						cb_RLE_lists[m][n][ptr2]=num2;
						if(num1==0&&num2==0){
							break;
						}
						ptr1=ptr2+1;
						ptr2=ptr1+1;
						
					}
					ptr1=0;
					ptr2=1;
					while(1){
						fscanf(fptr,"%s\n",buff);
						num1=decode_huffman(&book,buff);
						fscanf(fptr,"%s\n",buff);
						num2=decode_huffman(&book,buff);
						cr_RLE_lists[m][n][ptr1]=num1;
						cr_RLE_lists[m][n][ptr2]=num2;
						if(num1==0&&num2==0){
							break;
						}
						ptr1=ptr2+1;
						ptr2=ptr1+1;
						
					}
				}
			}
			
			//un RLE
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					un_RLE_8x8(y_RLE_lists[m][n],q_y_F_zigzag[m][n]);
					un_RLE_8x8(cb_RLE_lists[m][n],q_cb_F_zigzag[m][n]);
					un_RLE_8x8(cr_RLE_lists[m][n],q_cr_F_zigzag[m][n]);
				}
			}

			short****q_y_F=malloc_4D_short(M,N,8,8);//frequency domain of the 8x8 matrices
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			/*free_3D_short(y_RLE_lists);
			free_3D_short(cb_RLE_lists);
			free_3D_short(cr_RLE_lists);*/
			
			//unzigzag
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
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F[m][n][r][c]=q_y_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cb_F[m][n][r][c]=q_cb_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cr_F[m][n][r][c]=q_cr_F_zigzag[m][n][zz_matrix[r][c]-1];
						}
					}
				}
			}
			
			
			//un DPCM
			for(m=1;m<M;m++){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]+q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]+q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]+q_cr_F[m-1][0][0][0];
			}
			for(m=0;m<M;m++){
				int n;
				for(n=1;n<N;n++){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]+q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]+q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]+q_cr_F[m][n-1][0][0];
				}
			}
			
			
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
			
			
			
			
			
			double****y_F=malloc_4D_double(M,N,8,8);
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			
			
			
			//unquantize
			//FILE*debug_fptr=fopen("debug_decoder.txt","w");//debug
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_F[m][n][u][v]=(float)(q_y_F[m][n][u][v]*y_Q[u][v]);
							cb_F[m][n][u][v]=(float)(q_cb_F[m][n][u][v]*cbcr_Q[u][v]);
							cr_F[m][n][u][v]=(float)(q_cr_F[m][n][u][v]*cbcr_Q[u][v]);
							//fprintf(debug_fptr,"%lf %lf %lf\n",y_F[m][n][u][v],cb_F[m][n][u][v],cr_F[m][n][u][v]);//debug
							
						}
					}
				}
			}
			/*fclose(debug_fptr);
			return 0;*/
			
			//generate basis vector for 2D-IDCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int r,c;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			int u,v;
			for(u=0;u<H;u++){
				for(v=0;v<W;v++){
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
					}
				}
			}
			
			//做IDCT 
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			double**beta=ones(H,1);
			beta[0][0]=1.000000/sqrt(2);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							
							IDCT_8x8(y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
							y_tmp[r][c]+=128.000000;
							cb_tmp[r][c]+=128.000000;
							cr_tmp[r][c]+=128.000000;
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
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			
			ImWrite(img,argv[2]);
			return 0;
			
			
		
		}
		else if(argv[3][0]=='b'){
			//read codebook
			//printf("test\n");
			codebook book=read_codebook_txt(argv[4]);
			unsigned int*bin_val=(unsigned int*)malloc(book.node_num*sizeof(unsigned int));
			int val_i=0;
			for(val_i=0;val_i<book.node_num;val_i++){
				int j=0;
				unsigned int res=0;
				while(book.bitstream[val_i][j]!='\0'){
					res=res*2+(book.bitstream[val_i][j]-'0');
					j++;
				}
				bin_val[val_i]=res;
				//printf("%d %d\n",book.vals[val_i],bin_val[val_i]);
			}
			//return 0;
			//read header
			FILE*fptr=fopen(argv[5],"rb");
			Img img;
			fread(img.header.identifier,2, 1, fptr);
			fread(&img.header.filesize, sizeof(int), 1, fptr);
			fread(&img.header.reserved, sizeof(short), 1, fptr);
			fread(&img.header.reserved2, sizeof(short), 1, fptr);
			fread(&img.header.bitmap_dataoffset, sizeof(int), 1, fptr);
			fread(&img.header.bitmap_headersize, sizeof(int), 1, fptr);
			fread(&img.header.width, sizeof(int), 1, fptr);
			fread(&img.header.height, sizeof(int), 1, fptr);
			fread(&img.header.planes, sizeof(short), 1, fptr);
			fread(&img.header.bits_perpixel, sizeof(short), 1, fptr);
			fread(&img.header.compression, sizeof(int), 1, fptr);
			fread(&img.header.bitmap_datasize, sizeof(int), 1, fptr);
			fread(&img.header.hresolution, sizeof(int), 1, fptr);
			fread(&img.header.vresolution, sizeof(int), 1, fptr);
			fread(&img.header.usedcolors, sizeof(int), 1, fptr);
			fread(&img.header.importantcolors, sizeof(int), 1, fptr);
			img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);
			
			int M=img.header.height/8;
			int N=img.header.width/8;
			int H=8;
			int W=8;
			//read bitstream and decode
			int m;
			int n;
			short***y_RLE_lists=malloc_3D_short(M,N,64); 
			short***cb_RLE_lists=malloc_3D_short(M,N,64); 
			short***cr_RLE_lists=malloc_3D_short(M,N,64); 
			//char*buff=(char*)malloc(200*sizeof(char));
			short num1,num2;
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					int ptr1,ptr2;
					ptr1=0;
					ptr2=1;
					while(1){
						//fscanf(fptr,"%s\n",buff);
						unsigned int word_val;
						fread(&word_val,sizeof(unsigned int),1,fptr);
						num1=decode_huffman_by_val(&book,bin_val,word_val);
						
						fread(&word_val,sizeof(unsigned int),1,fptr);
						num2=decode_huffman_by_val(&book,bin_val,word_val);
						
						//printf("%d %d ",num1,num2);
						y_RLE_lists[m][n][ptr1]=num1;
						y_RLE_lists[m][n][ptr2]=num2;
						if(num1==0&&num2==0){
							break;
						}
						ptr1=ptr2+1;
						ptr2=ptr1+1;
						
					}
					//printf("\n");
					ptr1=0;
					ptr2=1;
					while(1){
						//fscanf(fptr,"%s\n",buff);
						unsigned int word_val;
						fread(&word_val,sizeof(unsigned int),1,fptr);

						num1=decode_huffman_by_val(&book,bin_val,word_val);
						
						fread(&word_val,sizeof(unsigned int),1,fptr);
						num2=decode_huffman_by_val(&book,bin_val,word_val);
						
						//printf("%d %d ",num1,num2);
						cb_RLE_lists[m][n][ptr1]=num1;
						cb_RLE_lists[m][n][ptr2]=num2;
						if(num1==0&&num2==0){
							break;
						}
						ptr1=ptr2+1;
						ptr2=ptr1+1;
						
					}
					//return 0;
					ptr1=0;
					ptr2=1;
					while(1){
						//fscanf(fptr,"%s\n",buff);
						unsigned int word_val;
						fread(&word_val,sizeof(unsigned int),1,fptr);
						num1=decode_huffman_by_val(&book,bin_val,word_val);
						
						fread(&word_val,sizeof(unsigned int),1,fptr);
						num2=decode_huffman_by_val(&book,bin_val,word_val);
						
						//printf("%d %d ",num1,num2);
						cr_RLE_lists[m][n][ptr1]=num1;
						cr_RLE_lists[m][n][ptr2]=num2;
						if(num1==0&&num2==0){
							break;
						}
						ptr1=ptr2+1;
						ptr2=ptr1+1;
						
					}
				}
			}
			//un RLE
			short***q_y_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cb_F_zigzag=malloc_3D_short(M,N,8*8);
			short***q_cr_F_zigzag=malloc_3D_short(M,N,8*8);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					un_RLE_8x8(y_RLE_lists[m][n],q_y_F_zigzag[m][n]);
					un_RLE_8x8(cb_RLE_lists[m][n],q_cb_F_zigzag[m][n]);
					un_RLE_8x8(cr_RLE_lists[m][n],q_cr_F_zigzag[m][n]);
				}
			}

			short****q_y_F=malloc_4D_short(M,N,8,8);//frequency domain of the 8x8 matrices
			short****q_cb_F=malloc_4D_short(M,N,8,8);
			short****q_cr_F=malloc_4D_short(M,N,8,8);
			/*free_3D_short(y_RLE_lists);
			free_3D_short(cb_RLE_lists);
			free_3D_short(cr_RLE_lists);*/
			
			//unzigzag
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
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int r,c;
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							q_y_F[m][n][r][c]=q_y_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cb_F[m][n][r][c]=q_cb_F_zigzag[m][n][zz_matrix[r][c]-1];
							q_cr_F[m][n][r][c]=q_cr_F_zigzag[m][n][zz_matrix[r][c]-1];
						}
					}
				}
			}
			
			
			//un DPCM
			for(m=1;m<M;m++){
				q_y_F[m][0][0][0]=q_y_F[m][0][0][0]+q_y_F[m-1][0][0][0];
				q_cb_F[m][0][0][0]=q_cb_F[m][0][0][0]+q_cb_F[m-1][0][0][0];
				q_cr_F[m][0][0][0]=q_cr_F[m][0][0][0]+q_cr_F[m-1][0][0][0];
			}
			for(m=0;m<M;m++){
				int n;
				for(n=1;n<N;n++){
					q_y_F[m][n][0][0]=q_y_F[m][n][0][0]+q_y_F[m][n-1][0][0];
					q_cb_F[m][n][0][0]=q_cb_F[m][n][0][0]+q_cb_F[m][n-1][0][0];
					q_cr_F[m][n][0][0]=q_cr_F[m][n][0][0]+q_cr_F[m][n-1][0][0];
				}
			}
			
			
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
			
			
			
			
			
			double****y_F=malloc_4D_double(M,N,8,8);
			double****cb_F=malloc_4D_double(M,N,8,8);
			double****cr_F=malloc_4D_double(M,N,8,8);
			
			
			
			//unquantize
			//FILE*debug_fptr=fopen("debug_decoder.txt","w");//debug
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					for(u=0;u<8;u++){
						int v;
						for(v=0;v<8;v++){
							y_F[m][n][u][v]=(float)(q_y_F[m][n][u][v]*y_Q[u][v]);
							cb_F[m][n][u][v]=(float)(q_cb_F[m][n][u][v]*cbcr_Q[u][v]);
							cr_F[m][n][u][v]=(float)(q_cr_F[m][n][u][v]*cbcr_Q[u][v]);
							//fprintf(debug_fptr,"%lf %lf %lf\n",y_F[m][n][u][v],cb_F[m][n][u][v],cr_F[m][n][u][v]);//debug
							
						}
					}
				}
			}
			/*fclose(debug_fptr);
			return 0;*/
			
			//generate basis vector for 2D-IDCT
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int r,c;
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			int u,v;
			for(u=0;u<H;u++){
				for(v=0;v<W;v++){
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
						}
					}
				}
			}
			
			//做IDCT 
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			double**beta=ones(H,1);
			beta[0][0]=1.000000/sqrt(2);
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							
							IDCT_8x8(y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
							y_tmp[r][c]+=128.000000;
							cb_tmp[r][c]+=128.000000;
							cr_tmp[r][c]+=128.000000;
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
			free_4D_double(y_F,M,N,8,8);
			free_4D_double(cb_F,M,N,8,8);
			free_4D_double(cr_F,M,N,8,8);
			
			
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			
			ImWrite(img,argv[2]);
			return 0;
		
		}
		
	}
}

void char2bin(unsigned char c,char*bin){
	unsigned int val=(unsigned int)c;
	printf("%d ",val);
	int bin_ptr=0;
	while(val>1){
		bin[bin_ptr]='0'+val%2;
		val/=2;
		bin_ptr++;
	}
	if(val==0)bin[bin_ptr]='0';
	else bin[bin_ptr]='1';
	bin_ptr++;
	int ptr1,ptr2;
	ptr1=0;
	ptr2=bin_ptr-1;
	while(ptr1<ptr2){
		char tmp=bin[ptr1];
		bin[ptr1]=bin[ptr2];
		bin[ptr2]=tmp;
		ptr1++;
		ptr2--;
	}
	bin[bin_ptr]='\0';
	
}
short decode_huffman_by_val(codebook*book,unsigned int*bin_val,unsigned int word_val){
	int val_i;
	for(val_i=0;val_i<(book->node_num);val_i++){
		if(bin_val[val_i]==word_val){
			break;
		}
	}
	if(val_i==book->node_num){
		printf("search %s codebook error\n",word_val);
		return -1;
	}
	return book->vals[val_i];
}

short decode_huffman(codebook*book,char*data){
	int val_i;
	//printf("%s",data);
	for(val_i=0;val_i<(book->node_num);val_i++){
		
		//printf("%d %s\n",book->vals[val_i],book->bitstream[val_i]);
		if(strcmp(book->bitstream[val_i],data)==0){
			break;
		}
	}
	if(val_i==book->node_num){
		printf("search %s codebook error\n",data);
		return -1;
	}
	return book->vals[val_i];
}


codebook read_codebook_txt(char*filename){
	codebook book;
	FILE*fptr=fopen(filename,"r");
	int val_ptr=0;
	book.vals=(short*)malloc(107*sizeof(short));
	book.bitstream=(char**)malloc(107*sizeof(char*));
	while(!feof(fptr)){
		int symbol;
		int count;
		unsigned char*buff=(unsigned char*)malloc(100*sizeof(unsigned char));
		int word_count=0;
		int bit_amount;
		fscanf(fptr,"%d %d %d %d ",&symbol,&count,&word_count,&bit_amount);
		int word_ptr;
		//printf("%d %d %d %d\n",symbol,count,word_count,bit_amount);
		for(word_ptr=0;word_ptr<word_count;word_ptr++){
			fscanf(fptr,"%c",&buff[word_ptr]);
		}
		
		int i;
		int stream_ptr=0;
		book.node_num=107;
		book.vals[val_ptr]=symbol;
		book.bitstream[val_ptr]=(char*)malloc((bit_amount+1)*sizeof(char));
		//printf("%d %d ",book.vals[val_ptr],count);
		for(i=0;i<word_count;i++){//convert the characters of a word in buff to intergers   
			int wordnum=buff[i];
			//printf("%d ",wordnum);
			char tmp[200];
			int tmp_i=0;//represent bits number of this word 
			while(wordnum>1){//calculate binary
				tmp[tmp_i]=wordnum%2+'0';
				wordnum/=2;
				tmp_i++;
			}
			if(wordnum==1)tmp[tmp_i]='1';
			else tmp[tmp_i]='0';
			tmp_i++;
			int ptr1,ptr2;
			ptr1=0;
			ptr2=tmp_i-1;
			while(ptr1<ptr2){//reverse
				char c=tmp[ptr1];
				tmp[ptr1]=tmp[ptr2];
				tmp[ptr2]=c;
				ptr1++;
				ptr2--;
			
			}
			if((i==word_count-1)&&(tmp_i<(bit_amount-8*(word_count-1)))){//last word
				int zero_count;
				int zero_need=(bit_amount-8*(word_count-1))-tmp_i;
				for(zero_count=1;zero_count<=zero_need;zero_count++){
					book.bitstream[val_ptr][stream_ptr]='0';
					stream_ptr++;
				}
			}
			else if((i<(word_count-1))&&tmp_i<8){
				int zero_count;
				int zero_need=8-tmp_i;
				for(zero_count=1;zero_count<=zero_need;zero_count++){
					book.bitstream[val_ptr][stream_ptr]='0';
					stream_ptr++;
				}
			}
			
			int bit_count;
			for(bit_count=0;bit_count<tmp_i;bit_count++){
				book.bitstream[val_ptr][stream_ptr]=tmp[bit_count];
				//printf("%c",book.bitstream[val_ptr][stream_ptr]);
				stream_ptr++;
			}
		}
		
		
		
		book.bitstream[val_ptr][stream_ptr]='\0';
		/*int j=0;
		while(book.bitstream[val_ptr][j]!='\0'){
			printf("%c",book.bitstream[val_ptr][j]);
			j++;
		}
		printf("\n");*/
		val_ptr++;
		if(symbol==1){
			break;
		}
	}
	//fclose(fptr);
	return book;
}

int s2i(char*s,int*ptr1,int*ptr2){
	//printf("%d %d\n",*ptr1,*ptr2);
	int res=0;
	int sign=1;
	if(s[*ptr1]=='-'){
		sign=-1;
	}
	else{
		res=s[*ptr1]-'0';
	}
	int i;
	for(i=(*ptr1)+1;i<*ptr2;i++){
		res=res*10+(s[i]-'0');
	}
	return res*sign;
}

FILE*print_RLE_lists(short*RLE_list,FILE*fptr){
	int i1,i2;
	i1=0;
	i2=1;
	while(1){
		fprintf(fptr,"%d %d ",RLE_list[i1],RLE_list[i2]);
		if(RLE_list[i1]==0&&RLE_list[i2]==0){
			break;
		}
		
		i1=i2+1;
		i2=i1+1;
	}
	fprintf(fptr,"%c",'\n');
	return fptr;
}


FILE*read_8x8_RLE_bin(FILE*fptr,short*list){
	int i1,i2;
	i1=0;
	i2=1;
	while(1){
		fread(&list[i1],sizeof(short),1,fptr);
		fread(&list[i2],sizeof(short),1,fptr);
		if(list[i1]==0&&list[i2]==0){
			break;
		}
		i1=i2+1;
		i2=i1+1;
	}
	
	return fptr;
}



void un_RLE_8x8(short*list,short*F_zigzag){
	int i=0;
	int zero_num=0;
	int ptr=1;
	while(!(list[zero_num]==0&&list[ptr]==0)){
		int count;
		for(count=1;count<=list[zero_num];count++){
			F_zigzag[i]=0;
			i++;
		}
		F_zigzag[i]=list[ptr];
		zero_num=ptr+1;
		ptr=zero_num+1;
		i++;
	}
	while(i<64){
		F_zigzag[i]=0;
		i++;
	}
}
FILE*read_8x8_RLE_txt(FILE*fptr,short*list){
	
	char c;
	while(c=fgetc(fptr)){
		if(c==')'){
			break;
		}
	}
	short tmp[64];
	int i=0;
	int count=0;
	while(1){
		fscanf(fptr,"%d %d ",&list[i],&list[i+1]);
		if(list[i]==0&&list[i+1]==0)break;
		count+=2;
		i+=2;
	}
	
	return fptr;
}



void read_table_txt(int table[8][8],char*filename){
	int i,j;
	FILE*fptr=fopen(filename,"r");
	for(i=0;i<8;i++){
		for(j=0;j<8;j++){
			int tmp;
			fscanf(fptr,"%d",&tmp);
			table[i][j]=tmp;
		}
		//printf("\n");
		//fgetc(fptr);
	}
}


void readHeader_txt(Img*img,char*filename){
	FILE*dim_fp=fopen(filename,"r");
	fscanf(dim_fp,"%c%c\n",&(img->header.identifier[0]),&(img->header.identifier[1]));
	fscanf(dim_fp,"%d\n",&(img->header.filesize));
	fscanf(dim_fp,"%d\n",&(img->header.reserved));
	fscanf(dim_fp,"%d\n",&(img->header.reserved2));
	fscanf(dim_fp,"%d\n",&(img->header.bitmap_dataoffset));
	fscanf(dim_fp,"%d\n",&(img->header.bitmap_headersize));
	fscanf(dim_fp,"%d\n",&(img->header.width));
	fscanf(dim_fp,"%d\n",&(img->header.height));
	fscanf(dim_fp,"%d\n",&(img->header.planes));
	fscanf(dim_fp,"%d\n",&(img->header.bits_perpixel));
	fscanf(dim_fp,"%d\n",&(img->header.compression));
	fscanf(dim_fp,"%d\n",&(img->header.bitmap_datasize));
	fscanf(dim_fp,"%d\n",&(img->header.hresolution));
	fscanf(dim_fp,"%d\n",&(img->header.vresolution));
	fscanf(dim_fp,"%d\n",&(img->header.usedcolors));
	fscanf(dim_fp,"%d\n",&(img->header.importantcolors));
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

void print_header(Img img){
	printf("\n");
	printf("identifier=%c%c\n",img.header.identifier[0],img.header.identifier[1]);
	printf("filesize=%d\n",img.header.filesize);
	printf("reserved=%d\n",img.header.reserved);
	printf("reserved2=%d\n",img.header.reserved2);
	printf("bitmap_dataoffset=%d\n",img.header.bitmap_dataoffset);
	printf("bitmap_headersize=%d\n",img.header.bitmap_headersize);
	printf("width=%d\n",img.header.width);
	printf("height=%d\n",img.header.height);
	printf("planes=%d\n",img.header.planes);
	printf("bits_perpixel=%d\n",img.header.bits_perpixel);
	printf("compression=%d\n",img.header.compression);
	printf("bitmap_datasize=%d\n",img.header.bitmap_datasize);
	printf("hresolution=%d\n",img.header.hresolution);
	printf("vresolution=%d\n",img.header.vresolution);
	printf("usedcolors=%d\n",img.header.usedcolors);
	printf("importantcolors=%d\n",img.header.importantcolors);
	//printf("palette=%d\n",img.header.palette);
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
void ImWrite(Img img,char filename[]){
	FILE *outfile;
    outfile= fopen(filename, "wb");
    reverse_row(img.pixels,img.header.height,img.header.width);
    //printf("test\n");
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
void IDCT_8x8(double**F,double**matrix_8x8,double**beta,double****basic_vector,int r,int c){
	int u,v;
	matrix_8x8[r][c]=0;
	for(u=0;u<8;u++){
		for(v=0;v<8;v++){
			matrix_8x8[r][c]+=2.000000/8*F[u][v]*beta[u][0]*beta[v][0]*basic_vector[u][v][r][c];
		}
	}
}
void ycbcr2bgr(Img*img,double**y_data,double**cb_data,double**cr_data){
	int i,j;
	for(i=0;i<(img->header).height;i++){
		for(j=0;j<(img->header).width;j++){
			double Y=y_data[i][j];
			double Cb=cb_data[i][j];
			double Cr=cr_data[i][j];
			double btmp=round((Y+1.772*Cb));
			
			if(btmp<0){
				//printf("b overflow %lf\n",btmp);
				btmp=0;
			}
			else if(btmp>255){
				//printf("b overflow %lf\n",btmp);
				btmp=255;
			}
			double gtmp=round(Y-0.34414*Cb-0.71414*Cr);
			if(gtmp<0){
				//printf("g overflow %lf\n",gtmp);
				gtmp=0;
			}
			else if(gtmp>255){
				//printf("g overflow %lf\n",gtmp);
				gtmp=255;
			}
			double rtmp=round((Y + 1.402*Cr));
			if(rtmp<0){
				//printf("r overflow %lf\n",rtmp);
				rtmp=0;
			}
			else if(rtmp>255){
				//printf("r overflow %lf\n",rtmp);
				rtmp=255;
			}
			img->pixels[i][j].data[0]=(unsigned char)btmp;
			img->pixels[i][j].data[1]=(unsigned char)gtmp;
			img->pixels[i][j].data[2]=(unsigned char)rtmp;
			//printf("ycbcr2bgr test\n");//debug
		}
	}
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
	printf("%x",res);
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
void RLE_8x8(short*F_zigzag,RLE_node*list){
	int i=0;
	RLE_node*last=list;
	short count_zero=0;
	while(i<64){
		if(F_zigzag[i]!=0){
			RLE_node*zero_num=new_RLE_node(count_zero);//前面的非0值到這個非0值有幾個0 
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

RLE_node*new_RLE_node(short val){
	RLE_node*res=(RLE_node*)malloc(sizeof(RLE_node));
	res->val=val;
	res->next=NULL;
	return res;
}


void add_RLE_node(RLE_node*last,RLE_node*newnode){
	if(last==NULL)printf("error");
	last->next=newnode;
}
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
