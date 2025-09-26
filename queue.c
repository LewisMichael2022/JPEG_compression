#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdint.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
//環狀不行是因為無法pop指定index的元素，也不能sort 
/*typedef struct RLE_node{
	int zero_num;
	double next_val;
	struct RLE_node*next_node;
}RLE_node;

typedef struct RLE_Huffman_queue_node{
	RLE_node*node;
	int count;
	struct RLE_Huffman_queue_node*left;
	struct RLE_Huffman_queue_node*right;
	struct RLE_Huffman_queue_node*next;
}RLE_Huffman_queue_node;

typedef struct RLE_Huffman_queue{
	RLE_Huffman_queue_node*front;//表頭,不是第0個資料 
	RLE_Huffman_queue_node*rear;
}RLE_Huffman_queue;
*/


typedef struct Huffman_node{
	int val;
	int count;
	struct Huffman_node*left;
	struct Huffman_node*right;
	struct Huffman_node*next;
}Huffman_node;

typedef struct Huffman_queue{
	Huffman_node*front;//表頭,不是第0個資料 
	Huffman_node*rear;
}Huffman_queue;

/*
RLE_Huffman_queue_node*new_RLE_Huffman_queue_node(RLE_node*rle_node,int count);
RLE_Huffman_queue RLE_Huffman_queue_init();
void RLE_Huffman_queue_push(RLE_Huffman_queue*queue,RLE_Huffman_queue_node*node);
void RLE_Huffman_queue_pop(RLE_Huffman_queue*queue,RLE_Huffman_queue_node*target);
void display(RLE_Huffman_queue*queue);
*/
Huffman_queue Huffman_queue_init();
void Huffman_queue_push(Huffman_queue*queue,Huffman_node*node);
void Huffman_queue_pop(Huffman_queue*queue,Huffman_node*target);
void display(Huffman_queue*queue);
int main(){
	Huffman_queue queue=Huffman_queue_init();
	int i=0;
	for(i=0;i<3;i++){
		Huffman_node*node=(Huffman_node*)malloc(sizeof(Huffman_node));
		node->count=i;
		Huffman_queue_push(&queue,node);
	}
	Huffman_node*ptr=(queue.front)->next;
	while(ptr!=NULL){
		printf("%d ",ptr->count);
		ptr=ptr->next;
	} 
	
	//display(&queue);
	return 0;
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
}

