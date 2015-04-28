#ifndef function_h_
#define function_h_
#include<stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
using namespace std;
/*
This files declares and defines the basic variables,
struct, functions of all the project
*/
/*@author: linhong (linhong.seba.zhu@gmail.com)
*/

//////////////////////////////////////////////////////////////////////////
/////////////global structure definition(start)//////////////////////////
/////////////////////////////////////////////////////////////////////////
struct Entry{
	int key;
	double *weight;
	int vid;
	int *nbv;
};

//////////////////////////////////////////////////////////////////////////
/////////////global structure definition(end)///////////////////////
/////////////////////////////////////////////////////////////////////

///----------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////
//global variable definition (start)//////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//variables that are related to memory management
//===========================================================================================================//

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parameters for tmp memory management (START)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const int BLK_NUM=1024;
const int BLK_SZ=4194304;	//1024*1024*4

char *memBlkAr[BLK_NUM];	//An array of pointers to allocated blocks
char *curMemPos=NULL;	//The first free pos (startpos) in the first free block for writing
char *curMemEnd=NULL;	//The end pos of the first free block for writing
int curBlk=0;
int endBlk=0;	//curBlk: the ID of the current block
									//endBlk: the ID of the last allocated block
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parameters for tmp memory management (END)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//==============================================================================================================//
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

//==============================================================================================================
//other variables
int *ccid; //connected component id
char *isvisit; //used in dfs traversal as a marker
int ccnum=0; //number of connectedc omponnet
int maxccsize=0; //maximum size of connected component
int maxccid=0; //the id of largerst connected componnet
/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/

/////////////////////////////////////////////////////////////////
//======================variable related to code length================================//
double *pa; //pagerank score pa[i] is the pagerank score for node i
double *qi;// qi[i] is coresponding to 
	///tau*(n-ni)/(n-1)\sum_{u\in part[i]}pa[u]+(1-tau)\sum_{u \in part[i]}{v\not\in part[i]}pa[u]w[u][v]
double *pam; //pam[i] is the pagerank score for modual i
int *Ni; //Ni[i] keeps the number of nodes in modual i
double *HPi; //HPi[i]=H(Pi)=
	//-qi[i]/(qi[i]+pam[i])*log(-qi[i]/(qi[i]+pam[i]))-\sum_{u\in part[i]}pa[u]/(qi[i]+pam[i])*log(pa[u]/(qi[i]+pam[i]))
// POWER METHOD OR WEIGHTED FREQUENCY (Like Map Equation code)
// Powermethod == 1 means using PageRank with TAU as the teleportation parameter
// Make sure TAU = 0 if Powermethod == 0 //content[i] is the content code length for modual i
 
//#define POWERMETHOD             0
int POWERMETHOD = 0;

// TAU - TELEPORTATION 
//#define TAU                     0
float TAU = 0.15;

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////global function declaration (start)///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//====================================================================================//
//basic functions
inline int binarysearch(int *arrays, int l, int h, int v); //search v against array arrays
inline int binarysearch(vector<int> arrays, int l, int h, int v); //search v against vector arrays
inline void insert(int *&arrays, int pos, int v, int length); //insert v into the arrays at pos 
//note that length is not updated here, will update outside;
inline void insert(double *&arrays, int pos, double v, int length);
void inline allocatepermemory(unsigned int size);
void inline allocatetmpmemory(unsigned int size);

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////global function declaration (end)///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////global function definition (start)///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

#define LOG2E 1.44269504088896340736 //log2(e)

inline double log2(double x){
    return  log(x) * LOG2E;
}
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
//========================================================================================================//
//functions that are related to memory management//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//check whether tmpblock can allocate a memory block with length size
void inline allocatetmpmemory(unsigned int size){
		if(size>=(size_t)(curMemEnd - curMemPos)||(curMemPos==NULL)||(curMemEnd - curMemPos)>BLK_SZ){//free mem in cur block is not enough
			if(curBlk < endBlk){
				//we have already allocated free blocks
				curMemPos=memBlkAr[++curBlk];
				curMemEnd=curMemPos+BLK_SZ;
			}else{
				//allocate a new block
				++endBlk;
				if(endBlk>=BLK_NUM){
					printf("system is unable to allocate more temporal memory\n");
					printf("number of block is %d\n",endBlk);
					exit(0);
				}
				curMemPos=(char*)malloc(BLK_SZ);
				if(curMemPos==NULL){
					printf("system is unable to allocate more temporal memory\n");
					printf("number of block is %d\n",curBlk);
					printf("number of static block is %d\n",curBlk2);
					exit(0);
				}
				memBlkAr[++curBlk]=curMemPos;
				curMemEnd=curMemPos+BLK_SZ;
			}
		}//end of if free mem is not enough
	}

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
					printf("number of block is %d\n",curBlk);
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
//void computeresults(const char* groundTruthFilename, const char *outputFilename,int *&partition, int numpart, int N){
//	ifstream fin(groundTruthFilename);
//    int *truelabel=new int[N];
//    int n1;
//	int t1;
//
//    // Input user stance
//    for (int k=0; k<N; k++){
//        fin >> n1;
//        fin >> t1;
//        truelabel[n1] = t1;
//    }
//
//}

// Create a MATLAB file that can be run to plot the communities with ground truth information
//void outputCommunities_MatlabFile(const char* groundTruthFilename, const char *outputMatlabFilename,int *&partition, int numpart, int N,const char *graphfilename){
//    ifstream fin(groundTruthFilename);
//    int *stance=new int[N];
//    int n1;
//	int t1;
//    // Input user stance
//    for (int k=0; k<N; k++){
//        fin >> n1;
//        fin >> t1;
//        stance[n1] = t1;
//    }
//    
//    // Create new module numbers for elements where x[i]>0
//    ofstream fout(outputMatlabFilename);
//    
//    // Group the communities together - nodeNum+1 (for MATLAB) \t community \t stance
//    fout << "Communities =  [\n";
//    for (int c=0; c<numpart; c++){
//        
//        for(int i=0; i<N; i++){
//            if (partition[i]== c){
//                fout << i+1 << ",\t" << c << ",\t" << stance[i] << ";\n";
//            }
//        }
//    }
//    fout << "];\n\n\n";
//
//  //  
//  //  // Output Edge Network to Matlab
//  //  fout << "A =  [\n";
//  //  for(int i=0; i<N; i++){
//		//for(int j=0;j<A[i].key;j++){
//		//	fout << A[i].weight[j] << "\t" ;
//		//}
//  //      // for(int j=0; j<N; j++){
//  //          // fout << A[i][j] << "\t" ;
//  //      // }
//  //      fout << ";\n";
//  //  }
//  //  fout << "];\n\n";
//	fout<<"load "<<graphfilename;
//
//    // Store each of the columns of Communities Matrix
//    fout << "[N,M] = size(Communities);\n";
//    fout << "nodeNums = Communities(:,1);\n";
//    fout << "classes = Communities(:,2);\n";
//    fout << "groundTruth = Communities(:,3);\n\n";
//    
//    // Find the nodes with edges
//    fout << "hasEdges = find(sum(A)>0);\n\n";
//    fout << "A2 = A(hasEdges, hasEdges);\n";
//
//    // Original Adjacency
//    fout << "figure, spy(A2)\n\n";
//    fout << "B = A(nodeNums, nodeNums);\n";
//    
//    // Reordered Adjacency according to communities
//    fout << "figure, spy(B)\n\n";
//    
//    // Nodes separated by community with color indicating a position (-1,0,1,other)
//    fout << "figure, hold on\n";
//    fout << "for i=1:N\n";
//    fout << "\tif groundTruth(i) == 0\n";
//    fout << "\t\tplot( classes(i),nodeNums(i),'k.', 'MarkerSize', 20)\n";
//    fout << "\telseif groundTruth(i) == 1\n";
//    fout << "\t\tplot( classes(i),nodeNums(i),'g.', 'MarkerSize', 20)\n";
//    fout << "\telseif groundTruth(i) == -1\n";
//    fout << "\t\tplot( classes(i),nodeNums(i),'r.', 'MarkerSize', 20)\n";
//    fout << "\telse\n";
//    fout << "\t\tplot( classes(i),nodeNums(i),'b.', 'MarkerSize', 20)\n";
//    fout << "\tend\n";
//    fout << "end\n\n";
//                
//    fout << "xlabel('Class Number')\n";
//    fout << "ylabel('Node Number')\n\n\n";
//	delete []stance;
//   return;
//}
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////global function definition (end)///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
#endif