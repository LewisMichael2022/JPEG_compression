#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdint.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
typedef struct Huffman_node{
	short val;
	int count;
	struct Huffman_node*left;
	struct Huffman_node*right;
	struct Huffman_node*next;
}Huffman_node;

typedef struct Huffman_queue{
	struct Huffman_node*front;//���Y,���O��0�Ӹ�� 
	struct Huffman_node*rear;
}Huffman_queue;

int main(){
		
}

Huffman_node*gen_Huffman_tree(short***y_RLE,short***cb_RLE,short***cr_RLE){
	
}






