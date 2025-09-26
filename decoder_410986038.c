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
unsigned char**malloc_2D_uchar(int row,int col);
void print_header(Img img);
void ImWrite(Img img,char filename[]);
Pixel** malloc_2D_Pixel(int row, int col,unsigned short bytes_perpixel);
double**malloc_2D_double(int row_size,int col_size);
short****malloc_4D_short(int size1,int size2,int size3,int size4);
RLE_node***malloc_3D_RLE_node(int M,int N);
short***malloc_3D_short(int size1,int size2,int size3);

void RLE_8x8(short*F_zigzag,RLE_node*list);
RLE_node*new_RLE_node(short val);
void add_RLE_node(RLE_node*last,RLE_node*newnode);

double**ones(int row,int col);
void reverse_row(Pixel**arr,int rowSize,int colSize);
void readHeader_txt(Img*img,char*filename);
void read_table_txt(double table[8][8],char*filename);
double****malloc_4D_double(int size1,int size2,int size3,int size4);
void IDCT_8x8(double**F,double**matrix_8x8,double**beta,double****basic_vector,int r,int c);
void ycbcr2bgr(Img*img,double**y_data,double**cb_data,double**cr_data);
FILE*read_8x8_RLE_txt(FILE*fptr,short*list);
FILE*read_8x8_RLE_bin(FILE*fptr,short*list);
void un_RLE_8x8(short*list,short*F_zigzag);
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
				fscanf(bptr,"%hhu",&img.pixels[i][j].data[0]);
				fscanf(gptr,"%hhu",&img.pixels[i][j].data[1]);
				fscanf(rptr,"%hhu",&img.pixels[i][j].data[2]);
			}
		}
		ImWrite(img,argv[2]);
	}
	
	
	if(argv[1][0]=='1'){
		if(argc==11){
			//decoder 1 QResKimberly.bmp Kimberly.bmp 4~6Qtable dim.txt 8~10qF
			
			//quantization table
			double y_Q[8][8];
			double cbcr_Q[8][8];
			read_table_txt(y_Q,argv[4]);
			read_table_txt(cbcr_Q,argv[5]);
			Img img;
			//get header
			readHeader_txt(&img,argv[7]);
			img.pixels=malloc_2D_Pixel(img.header.height,img.header.width,3);
			
			
			//get frequency quantized from txt file
			FILE*q_y_fptr=fopen(argv[8],"rb");
			FILE*q_cb_fptr=fopen(argv[9],"rb");
			FILE*q_cr_fptr=fopen(argv[10],"rb");
			
			
			int M=img.header.height/8;
			int N=img.header.width/8;
			int H=8;
			int W=8;
			double****q_y_F=malloc_4D_double(M,N,8,8);
			double****q_cb_F=malloc_4D_double(M,N,8,8);
			double****q_cr_F=malloc_4D_double(M,N,8,8);
			int m,n,u,v;
			//從紀錄量化值的raw檔讀取下來 
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					for(u=0;u<8;u++){
						for(v=0;v<8;v++){
							short qtmp,cbtmp,crtmp;
							fread(&qtmp,sizeof(short),1,q_y_fptr);
							fread(&cbtmp,sizeof(short),1,q_cb_fptr);
							fread(&crtmp,sizeof(short),1,q_cr_fptr);
							q_y_F[m][n][u][v]=(double)qtmp;
							q_cb_F[m][n][u][v]=(double)cbtmp;
							q_cr_F[m][n][u][v]=(double)crtmp;
							
							
						}
					}
				}
			}
			
			
			
			//unquantization
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					//unquantize and round each 8x8 matrix
					for(u=0;u<8;u++){ 
						for(v=0;v<8;v++){ 
							q_y_F[m][n][u][v]=(q_y_F[m][n][u][v]*y_Q[u][v]);
							q_cb_F[m][n][u][v]=(q_cb_F[m][n][u][v]*cbcr_Q[u][v]);
							q_cr_F[m][n][u][v]=(q_cr_F[m][n][u][v]*cbcr_Q[u][v]);
							
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
			beta[0][0]=1/sqrt(2);
			Pixel**res_data=malloc_2D_Pixel(img.header.height,img.header.width,3);
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							IDCT_8x8(q_y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
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
			
			//ycbcr -> rgb
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			ImWrite(img,argv[2]);
			//sqnr
		
		}
		
		else if(argc==13){
			//decoder 1 ResKimberly.bmp 3~5Qtable dim.txt 7~9qF 10~12qe	
			
			double y_Q[8][8];
			double cbcr_Q[8][8];
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
			double****q_y_F=malloc_4D_double(M,N,8,8);
			double****q_cb_F=malloc_4D_double(M,N,8,8);
			double****q_cr_F=malloc_4D_double(M,N,8,8);
			//從紀錄量化值的txt檔讀取下來 
			int m,n,u,v;
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					for(u=0;u<8;u++){
						for(v=0;v<8;v++){
							short ytmp,cbtmp,crtmp;
							fread(&ytmp,sizeof(short),1,q_y_fptr);
							fread(&cbtmp,sizeof(short),1,q_cb_fptr);
							fread(&crtmp,sizeof(short),1,q_cr_fptr);
							q_y_F[m][n][u][v]=(double)ytmp;
							q_cb_F[m][n][u][v]=(double)cbtmp;
							q_cr_F[m][n][u][v]=(double)crtmp;
							
						}
					}
					
				}
			}
			
			
			
			//unquantize
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					//unquantize and round each 8x8 matrix
					for(u=0;u<8;u++){ 
						for(v=0;v<8;v++){ 
							q_y_F[m][n][u][v]*=y_Q[u][v];
							q_cb_F[m][n][u][v]*=cbcr_Q[u][v];
							q_cr_F[m][n][u][v]*=cbcr_Q[u][v];
						}
					}
				}
			} 
			//讀取紀錄誤差的txt檔並和頻率相加 
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
							q_y_F[m][n][u][v]+=y_e;
							q_cb_F[m][n][u][v]+=cb_e;
							q_cr_F[m][n][u][v]+=cr_e;
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
			beta[0][0]=1/sqrt(2);
			Pixel**res_data=malloc_2D_Pixel(img.header.height,img.header.width,3);
			for(m=0;m<M;m++){
				for(n=0;n<N;n++){
					double**y_tmp=ones(8,8);
					double**cb_tmp=ones(8,8);
					double**cr_tmp=ones(8,8);
					//Get the spatial domain of the 8x8 frequency matrix with IDCT
					for(r=0;r<8;r++){
						for(c=0;c<8;c++){
							IDCT_8x8(q_y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
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
			
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			ImWrite(img,argv[2]);
			
		
		}
	}
	if(argv[1][0]=='2'){
		if(argv[3][0]=='a'){//decoder 2 QResKimberly.bmp ascii rle_code.txt 
			FILE*fptr=fopen(argv[4],"r");
			//read header
			Img img;
			fscanf(fptr,"%c%c ",&(img.header.identifier[0]),&(img.header.identifier[1]));
			fscanf(fptr,"%u ",&(img.header.filesize));
			fscanf(fptr,"%hhu ",&(img.header.reserved));
			fscanf(fptr,"%hhu ",&(img.header.reserved2));
			fscanf(fptr,"%u ",&(img.header.bitmap_dataoffset));
			fscanf(fptr,"%u ",&(img.header.bitmap_headersize));
			fscanf(fptr,"%u ",&(img.header.width));
			fscanf(fptr,"%u ",&(img.header.height));
			fscanf(fptr,"%hhu ",&(img.header.planes));
			fscanf(fptr,"%hhu ",&(img.header.bits_perpixel));
			fscanf(fptr,"%u ",&(img.header.compression));
			fscanf(fptr,"%u ",&(img.header.bitmap_datasize));
			fscanf(fptr,"%u ",&(img.header.hresolution));
			fscanf(fptr,"%u ",&(img.header.vresolution));
			fscanf(fptr,"%u ",&(img.header.usedcolors));
			fscanf(fptr,"%u",&(img.header.importantcolors));
			int M=img.header.height/8;
			int N=img.header.width/8;
			int H=8;
			int W=8;
			
			
			
			short***y_RLE_lists=malloc_3D_short(M,N,64); 
			short***cb_RLE_lists=malloc_3D_short(M,N,64); 
			short***cr_RLE_lists=malloc_3D_short(M,N,64); 
			int m;
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					fptr=read_8x8_RLE_txt(fptr,y_RLE_lists[m][n]);
					fptr=read_8x8_RLE_txt(fptr,cb_RLE_lists[m][n]);
					fptr=read_8x8_RLE_txt(fptr,cr_RLE_lists[m][n]);
					
					
				}
			}
			
			
			
			/*
			Because we want to put the basic vector picture in the directory we specified, 
			we need to change the workspace to the directory we specified
			*/
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int u;
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
			
			double****q_y_F=malloc_4D_double(M,N,8,8);//frequency domain of the 8x8 matrices
			double****q_cb_F=malloc_4D_double(M,N,8,8);
			double****q_cr_F=malloc_4D_double(M,N,8,8);
			free(y_RLE_lists);
			free(cb_RLE_lists);
			free(cr_RLE_lists);
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
			//unzigzag
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
			
			//unquantization
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					//unquantize and round each 8x8 matrix
					for(u=0;u<8;u++){ 
						int v;
						for(v=0;v<8;v++){ 
							q_y_F[m][n][u][v]*=y_Q[u][v];
							q_cb_F[m][n][u][v]*=cbcr_Q[u][v];
							q_cr_F[m][n][u][v]*=cbcr_Q[u][v];
						}
					}
				}
			} 
			
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			
			//2D-IDCT
			Pixel**res_data=malloc_2D_Pixel(img.header.height,img.header.width,3);
			double**beta=ones(H,1);
			beta[0][0]=1/sqrt(2);
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
							IDCT_8x8(q_y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
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
			img.pixels=res_data;
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			ImWrite(img,argv[2]);	
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
		
			double****basic_vector=malloc_4D_double(8,8,8,8);
			int u;
			for(u=0;u<H;u++){
				int v;
				for(v=0;v<W;v++){
					int r;
					//Make a 8x8 basic vector matrix x
					//double x[8][8];
					for(r=0;r<H;r++){//Calculate the value of the basic vector matrix
						int c;
						for(c=0;c<W;c++){ 
							basic_vector[u][v][r][c]=cos(M_PI*u*(2*r+1)/2/H) * cos(M_PI*v*(2*c+1)/2/W);
							//x[r][c]=127.5+127.5*basic_vector[u][v][r][c];
					
						}
				
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
			
			double****q_y_F=malloc_4D_double(M,N,8,8);//frequency domain of the 8x8 matrices
			double****q_cb_F=malloc_4D_double(M,N,8,8);
			double****q_cr_F=malloc_4D_double(M,N,8,8);
			free(y_RLE_lists);
			free(cb_RLE_lists);
			free(cr_RLE_lists);
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
			//unzigzag
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
			
			//unquantization
			for(m=0;m<M;m++){
				int n;
				for(n=0;n<N;n++){
					int u;
					//unquantize and round each 8x8 matrix
					for(u=0;u<8;u++){ 
						int v;
						for(v=0;v<8;v++){ 
							q_y_F[m][n][u][v]*=y_Q[u][v];
							q_cb_F[m][n][u][v]*=cbcr_Q[u][v];
							q_cr_F[m][n][u][v]*=cbcr_Q[u][v];
						}
					}
				}
			} 
			
			double**y_data=malloc_2D_double(img.header.height,img.header.width);
			double**cb_data=malloc_2D_double(img.header.height,img.header.width);
			double**cr_data=malloc_2D_double(img.header.height,img.header.width);
			
			//2D-IDCT
			Pixel**res_data=malloc_2D_Pixel(img.header.height,img.header.width,3);
			double**beta=ones(H,1);
			beta[0][0]=1/sqrt(2);
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
							IDCT_8x8(q_y_F[m][n],y_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cb_F[m][n],cb_tmp,beta,basic_vector,r,c);
							IDCT_8x8(q_cr_F[m][n],cr_tmp,beta,basic_vector,r,c);
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
			img.pixels=res_data;
			ycbcr2bgr(&img,y_data,cb_data,cr_data);
			ImWrite(img,argv[2]);	
			
			
		}
	}
	if(argv[1][0]=='3'){
		
		if(argv[3][0]=='a'){//decoder 3 QResKimberly.bmp ascii codebook.txt huffman_code.txt 
			
		
		}
		
		
		
		
		
	}
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



void read_table_txt(double table[8][8],char*filename){
	int i,j;
	FILE*fptr=fopen(filename,"r");
	for(i=0;i<8;i++){
		for(j=0;j<8;j++){
			int tmp;
			fscanf(fptr,"%d ",&tmp);
			table[i][j]=tmp;
		}
		//printf("\n");
		fgetc(fptr);
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
    printf("test\n");
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
			matrix_8x8[r][c]+=2.0/8*F[u][v]*beta[u][0]*beta[v][0]*basic_vector[u][v][r][c];
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
			(img->pixels[i][j]).data[0]=(unsigned int)((Y + 1.772*Cb)+0.5);
			(img->pixels[i][j]).data[1]=(unsigned int)((Y-0.344*Cb-0.714*Cr)+0.5);
			(img->pixels[i][j]).data[2]=(unsigned int)((Y + 1.402*Cr)+0.5);
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
