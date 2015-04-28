#ifndef graph_h_
#define graph_h_
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <vector>
using namespace std;
/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
// POWER METHOD OR WEIGHTED FREQUENCY (Like Map Equation code)
// Powermethod == 1 means using PageRank with TAU as the teleportation parameter
int POWERMETHOD = 0;
// TAU - TELEPORTATION 
double TAU = 0.15;
/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parameters for permanent memory management (START)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const int BLK_NUM2=1024;
const int BLK_SZ2=4194304;	//1024*1024*4

char *memBlkAr2[BLK_NUM2];	//An array of pointers to allocated blocks
char *curMemPos2=NULL;	//The first free pos (startpos) in the first free block for writing
char *curMemEnd2=NULL;	//The end pos of the first free block for writing
int curBlk2=0;	//curBlk: the ID of the current block
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parameters for permanent memory management (END)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//check whether the pemblock can allocate a memory block with length size
void inline allocatepermemory(unsigned int size){
	if(size>=(size_t)(curMemEnd2 - curMemPos2)||(curMemPos2==NULL)||(curMemEnd2 - curMemPos2)>BLK_SZ2){//free mem in cur block is not enough
		//free mem in cur block is not enough
		//allocate a new block
		if(curBlk2<BLK_NUM2-1){
			++curBlk2;
			curMemPos2=(char*)malloc(BLK_SZ2);
			if(curMemPos2==NULL){
				printf("system is unable to allocate more static memory\n");
				printf("number of static block is %d\n",curBlk2);
				exit(0);
			}
			memBlkAr2[curBlk2]=curMemPos2;
			curMemEnd2=curMemPos2+BLK_SZ2;
		}else{
			printf("system is unable to allocate more static memory\n");
			exit(0);
		}
	}
}
/*@author: linhong (linhong.seba.zhu@gmail.com)
*/
struct Entry{
	int key;
	double *weight;
	int vid;
	int *nbv;
};



inline int binarysearch(int *arrays, int l, int h, int v){
	int cur=l;
	int low=l;
	int high=h;
	while(low<high){
		cur=low+((high-low)/2);
		if(arrays[cur]==v)
				return cur;
			else{
				if(arrays[cur]>v)
					high=cur;
				else
					low=cur+1;
			}
		}
	return -low-1;
}
inline int binarysearch(vector<int> arrays, int l, int h, int v){
	int cur=l;
	int low=l;
	int high=h;
	while(low<high){
		cur=low+((high-low)/2);
		if(arrays.at(cur)==v)
				return cur;
			else{
				if(arrays.at(cur)>v)
					high=cur;
				else
					low=cur+1;
			}
		}
	return -low-1;
}

inline void insert(int *&arrays, int pos, int v, int length){
	int k=0;
	for(k=length;k>pos;k--)
		arrays[k]=arrays[k-1];
	arrays[pos]=v;
}


inline void insert(double *&arrays, int pos, double v, int length){
	int k=0;
	for(k=length;k>pos;k--)
		arrays[k]=arrays[k-1];
	arrays[pos]=v;
}
void readgraph(const char *filename, int &N, Entry *&G){
	FILE *refile=fopen(filename,"r");
	if(refile==NULL){
		printf("could not open graphfile\n");
		exit(0);
	}
	fscanf(refile,"%d\n",&N);
	G=(Entry*)malloc(sizeof(Entry)*N);
	if(G==NULL){
		printf("could not allocate memory\n");
		exit(0);
	}
	int i=0;
	int j=0;
	int u;
	int v;
	double w;
	int deg;
	for(i=0;i<N;i++){
		fscanf(refile,"%d,%d",&u,&deg);
		G[u].vid=u;
		G[u].key=deg;
		if(deg>0){
			G[u].nbv=(int*)malloc(sizeof(int)*deg);
			G[u].weight=(double*)malloc(sizeof(double)*deg);
			for(j=0;j<deg;j++){
				fscanf(refile,":%d,%lf",&v,&w);
				G[u].nbv[j]=v;
				G[u].weight[j]=w;
			}
		}else{
			G[u].nbv=NULL;
			G[u].weight=NULL;
		}
	}
	fclose(refile);
	// for(int i=0;i<N;i++){
		// for(int j=0;j<G[i].key;j++){
			// cout<<i<<"\t"<<G[i].nbv[j]<<"\t"<<G[i].weight[j]<<endl;
		// }
	// }
}
#endif