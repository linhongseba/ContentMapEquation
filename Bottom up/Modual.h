#ifndef Modual_h_
#define Modual_h_
#include<stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include "graph.h"
/*
*@author: linhong (linhong.seba.zhu@gmail.com)
@last update time: 02, July, 2013
*/
class Modual{	
	public:
		int *nbv;
		int *src;
		double *weights;
		Modual();
		Modual(int num);
		~Modual();
		int size;
		int maxsize;
		void addedge(int u, int v, double w);
		void clear();
		void reset();
};

Modual::Modual(){
	maxsize=5;
	size=0;
	nbv=(int*)malloc(sizeof(int)*(maxsize+1));
	src=(int*)malloc(sizeof(int)*(maxsize+1));
	weights=(double*)malloc(sizeof(double)*(maxsize+1));
}
Modual:: Modual(int num){
	maxsize=num;
	size=0;
	nbv=(int*)malloc(sizeof(int)*(num+1));
	src=(int*)malloc(sizeof(int)*(num+1));
	weights=(double*)malloc(sizeof(double)*(num+1));
}
Modual::~Modual(){
}

void Modual::addedge(int u, int v, double w){
	if(size==0){
		src[size]=u;
		nbv[size]=v;
		weights[size]=w;
		size++;
	}else{
		if(size>=maxsize){
				//reallocate memory;
			maxsize+=10;
			nbv=(int*)realloc(nbv,sizeof(int)*maxsize);
			src=(int*)realloc(src,sizeof(int)*maxsize);
			weights=(double*)realloc(weights,sizeof(double)*maxsize);
			}
		src[size]=u;
		nbv[size]=v;
		weights[size]=w;
		size++;
	}

	// if(size==0){
		// nbv[size]=v;
		// weights[size]=w;
		// size++;
	// }else{
		// int cur=binarysearch(nbv,0,size,v);
		// if(cur<0){
			//check memory;
			// if(size>=maxsize){
				//reallocate memory;
				// maxsize+=10;
				// nbv=(int*)realloc(nbv,sizeof(int)*maxsize);
				// weights=(double*)realloc(weights,sizeof(double)*maxsize);
			// }
			// insert(nbv,-cur-1,v,size);
			// insert(weights,-cur-1,w,size);
			// size++;
		// }
	// }
}

void Modual::clear(){
	if(nbv!=NULL){
		free(nbv);
		nbv=NULL;
	}
	if(src!=NULL){
		free(src);
		src=NULL;
	}
	if(weights!=NULL){
		free(weights);
		weights=NULL;
	}
	size=0;
	maxsize=0;
}
void Modual::reset(){
	size=0;
}
#endif