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

Pixel** malloc_2D_Pixel(int row, int col,unsigned short bytes_perpixel);
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

int main(){
	Img img;
	img.pixels=malloc_2D_Pixel(8,8,3);
	int i,j;
	for(i=0;i<8;i++){
		for(j=0;j<8;j++){
			img.pixels[i][j].data[0]=1;
			img.pixels[i][j].data[1]=2;
			img.pixels[i][j].data[2]=3;
		}
	}
	
}
