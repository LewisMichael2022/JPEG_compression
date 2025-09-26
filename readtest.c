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
char get_fptr_char(FILE*fptr);
FILE*move_fptr(FILE*fptr);
int main(){

	/*while(!feof(fptr)){
		char*buff=(char*)malloc(100*sizeof(char));
		fgets(buff,100,fptr);
		if(buff[0]=='\n'&&buff[1]=='\0'){
			printf("%c",'\n');
			continue;
		}
		
		//printf("%s",buff);
		int ptr1,ptr2;
		ptr1=ptr2=0;
		while(buff[ptr2]!=' '){
			ptr2++;
		}
		int symbol=s2i(buff,&ptr1,&ptr2);
		
		
		ptr1=ptr2+1;
		ptr2=ptr1;
		while(buff[ptr2]!=' '){
			ptr2++;
		}
		int count=s2i(buff,&ptr1,&ptr2);
		
		ptr1=ptr2+1;
		char*word=buff+ptr1;
		printf("%d %d %s",symbol,count,word);
	
	}*/
	FILE*fptr=fopen("readtest.txt","r");
	/*printf("before:%x\n",fptr);
	char c=fgetc(fptr);
	printf("after:%x",fptr);*/
	
	
	int line=0;
	int change_line=0;
	char*buff=(char*)malloc(100*sizeof(char));
	while(!feof(fptr)){
		int symbol;
		int count;
		fscanf(fptr,"%d %d ",&symbol,&count);
		printf("%d %d ",symbol,count);
		int i=0;
		char c;
		while(1){
			c=get_fptr_char(fptr);
			if(c=='\n'){
				break;
			}
			printf("%x ",fptr);
			//printf("%d ",i);
			buff[i]=c;	
			fptr=move_fptr(fptr);
			i++;
		}
		while(1){
			c=get_fptr_char(fptr);
			if(c!='\n'){
				break;
			}
			//printf("%d ",i);
			buff[i]=c;
			fptr=move_fptr(fptr);
			i++;
		}
		
		
	}

}
char get_fptr_char(FILE*fptr){
	return fgetc(fptr);
}
FILE*move_fptr(FILE*fptr){
	if(fgetc(fptr)==EOF){
		printf("error\n");
	}
	return fptr;
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
