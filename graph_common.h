#ifndef graph_common_h_
#define graph_common_h_

using namespace std;

void normalizegraph(Entry *&G, const int N){
    // compute row sums
    double* RowSum = new double[N];
    for(int i=0;i<N;i++){
        RowSum[i]=0;
		for(int j=0;j<G[i].key;j++){
			RowSum[i]+=G[i].weight[j];
		}
    }
    
    // normalized the row sum to one
    for(int i=0; i<N; i++){
        if (RowSum[i]>0){
			for(int j=0;j<G[i].key;j++){
				G[i].weight[j]=G[i].weight[j]/RowSum[i];
			}
        }
    }
	delete []RowSum;
    return;
}

////////////////////////////////////////////////////////////
//compute the page rank of nodes in graph G///////////////
//////the results are stored in vector x////////////////
void computepagerank(Entry *&G, int N, double *&x) {
	double* c = new double[N];
	int i=0;
	int j=0;
	// Power Method - PageRank to get stationary distributions
	if (POWERMETHOD) {
        double d=0;
        double temp;
        double dangling;
        for(i=0; i<N; i++){
            x[i]=1.0/N;
        }

        do {
            for(i=0;i<N;i++){ // initialize c
                c[i]=0.0;
            }

            for(i=0;i<N;i++){ // distribute PageRank score among neighbors
				for (int j=0;j<G[i].key;j++) {
					int neigh=G[i].nbv[j];
					c[neigh]+=x[i]/G[i].key*G[i].weight[j]; // weighted PageRank
				}
            }

            dangling=0.0;
            for(i=0;i<N;i++){ // compute the score distributed by dangling nodes
                if (G[i].key == 0){
                    dangling+=x[i]/N;
                }
            }

            d=0.0;
            for(i=0;i<N;i++){
                c[i]=TAU/N + (1.0-TAU)*(c[i] + dangling); // add teleportation
                d+=c[i]; // compute the sum
            }

            temp = 0;
            for(i=0;i<N;i++){
                c[i]/=d;
                temp+=fabs(x[i]-c[i]);
                x[i]=c[i];
            }
        } while(temp > 5.0e-5);
    }
    else{  // Use relative weight for undirected weighted networks
        double d = 0;
        for(i=0;i<N;i++){
            c[i]=0;
			for(j=0;j<G[i].key;j++){
				c[i]+=G[i].weight[j];
			}
            d += c[i];
        }
        for(int i=0;i<N;i++){
			if(d>0){
				x[i]=c[i]/d;
            }
			else{
				x[i]=0;
            }
        }
    }
	delete []c;
}
void readpagerank(const char *filename, int N, double *&x){
	FILE * rfile=fopen(filename,"r");
	if(rfile==NULL){
		printf("could not open file to read\n");
		exit(0);
	}
	int uid;
	double value;
	while(feof(rfile)==false){
		fscanf(rfile,"%d\t%lf\n",&uid,&value);
		if(uid>N||uid<1)
			continue;
		x[uid-1]=value;
	}
}
#endif // graph_common_h