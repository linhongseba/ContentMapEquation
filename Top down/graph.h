#ifndef graph_h_
#define graph_h_
#include"function.h"
#include <fstream>
using namespace std;
/*
This file defines all the functions that are related to graph manipulations
read graphs
page rank
connected component
and anything else
*/
/*@author: linhong (linhong.seba.zhu@gmail.com)
*/

//read graph from file in disk
void readgraph(const char *filename, int &N, Entry *&G){
	FILE *refile=fopen(filename,"r");
	if(refile==NULL){
		printf("could not open graphfile\n");
		exit(0);
	}
	fscanf(refile,"%d\n",&N);
	G=(Entry*)malloc(sizeof(Entry)*N);
	if(G==NULL){
		printf("could not allocate memory for graph\n");
		exit(0);
	}
	int i=0;
	int j=0;
	for(i=0;i<N;i++){
		G[i].key=0;
		G[i].nbv=NULL;
		G[i].vid=0;
		G[i].weight=0;
	}
	int u;
	int v;
	double w;
	int deg;
	for(i=0;i<N;i++){
		fscanf(refile,"%d,%d",&u,&deg);
		if(u>=N&&u<0)
			continue;
		G[u].vid=u;
		G[u].key=deg;
		if(deg>0){
			if(deg>100000){
				G[u].nbv=(int*)malloc(sizeof(int)*deg);
				G[u].weight=(double*)malloc(sizeof(double)*deg);
			}else{
				allocatepermemory(sizeof(int)*deg);
				G[u].nbv=(int*)curMemPos2;
				curMemPos2+=(sizeof(int)*deg);
				allocatepermemory(sizeof(double)*deg);
				G[i].weight=(double*)curMemPos2;
				curMemPos2+=(sizeof(double)*deg);
			}
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

/*
load the dictionary file into memory
and store as a struct array dict
*/
void readdictionary(const char *filename, int N, Entry *&dict, int &wordnum){
	ifstream fin(filename);
	if (!fin.good())
		 cout << "FAILED\n" << filename;
	dict=(Entry*)malloc(sizeof(Entry)*N);
	if(dict==NULL){
		printf("could not allocate memory for dictionary\n");
		exit(0);
	}
	 int i=0;
	 int j=0;
	 int nodeid;
	 int wordlength;
	 int word;
	 double weight;
	 for(i=0;i<N;i++){
		 dict[i].key=0;
		 dict[i].nbv=NULL;
		 dict[i].weight=NULL;
		 dict[i].vid=0;
	 }
	 for(i=0;i<N;i++){
		 fin>>nodeid;
		 if(nodeid>=N||nodeid<0)
			 continue;
		 if(fin.eof())
                break;
		 fin>>wordlength;
		 dict[nodeid].vid=nodeid;
		 dict[nodeid].key=wordlength;
		 if(wordlength>0){
			 if(wordlength>100000){
				 dict[nodeid].nbv=(int*)malloc(sizeof(int)*wordlength);
				 dict[nodeid].weight=(double*)malloc(sizeof(double)*wordlength);	
			 }else{
				 allocatepermemory(sizeof(int)*wordlength);
				 dict[nodeid].nbv=(int*)curMemPos2;
				 curMemPos2+=(sizeof(int)*wordlength);
				 allocatepermemory(sizeof(double)*wordlength);
				 dict[nodeid].weight=(double*)curMemPos2;
				 curMemPos2+=(sizeof(double)*wordlength);
			 }
			 for(j=0;j<wordlength;j++){
				 fin>>word;
				 fin>>weight;
				 dict[nodeid].nbv[j]=word;
				 dict[nodeid].weight[j]=weight;
				 if(word>wordnum)
					 wordnum=word;
			 }
		 }else{
			 dict[nodeid].nbv=NULL;
			 dict[nodeid].weight=NULL;
		 }
	 }
}

void readpartition(const char *filename, int N, int *&partid, int &numpart){
	ifstream fin(filename);
	if (!fin.good()){
		 cout << "FAILED reading " << filename<<endl;
		 exit(0);
	}
	int i=0;
	int id;
	int nodeid;
	numpart=0;
	for(i=0;i<N;i++){
		fin>>nodeid;
		fin>>id;
		partid[nodeid]=id;
		if(id>numpart)
			numpart=id;
	}
	numpart++;
	fin.close();
}

int DFS2(Entry *&G, int v,int nodenum,int ccnum){
	int *visitstack=(int*)malloc(sizeof(int)*nodenum);
	int maxstacksize=nodenum;
	int used=0;
	visitstack[used]=v;
	used++;
	int count=0;
	int w=0;
	int u=0;
	bool flag=false;
	int i=0;
	while(used>0){
		u=visitstack[used-1];
		if(isvisit[u]=='b'){
			used--;
			continue;
		}
		flag=true;;
		for(i=0;i<G[u].key;i++){
			w=G[u].nbv[i];
			if(isvisit[w]=='v'){
				if(used==maxstacksize){
					visitstack=(int*)realloc(visitstack,sizeof(int)*(maxstacksize+50));
					maxstacksize+=50;
					if(visitstack==NULL){
						printf("system could not allocate more memory\n");
						exit(0);
					}
				}
				visitstack[used]=w;
				used++;
				flag=false;
				isvisit[w]='g';
			}
		}
		if(flag==true&&isvisit[u]!='b'){
			isvisit[u]='b';
			count++;
			ccid[u]=ccnum;
		}
	}
	//curMemPos=pos;
	if(visitstack!=NULL){
		free(visitstack);
		visitstack=NULL;//release memory in local function
	}
	
	return count;
}

//find the connected componnet of graph G
//the results are store in array CC
void Connectedcomponent(Entry *&G, int nodenum){
	isvisit=new char[nodenum];
	memset(isvisit,'v',sizeof(char)*nodenum);
	int i=0;
	ccnum=0;
	int size=0;
	for(i=0;i<nodenum;i++){
		if(isvisit[i]=='v'){
			size=DFS2(G,i,nodenum,ccnum);
			if(size>maxccsize){
				maxccsize=size;
				maxccid=ccnum;
			}
			ccnum++;
		}		
	}
	if(isvisit!=NULL){
		delete []isvisit;
		isvisit=NULL;
	}
}
#endif
