/* --------------------------------------------
content_algorithm_withDictionary.cpp
Modified by Laura Smith on 6/11/13.
Updated by Linhong on 07/21/2013
g++ dictionary_algorithm.cpp  inputGraph.h ModuleDictionary.h outputGraph.h -o dictionary
 ./dictionary graphFile.txt
 
 -------------------------------------------- */

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include "graph.h"
#include "../graph_common.h"
#include "ModuleDictionary.h"
#include "Modual.h"
#include "outputGraph.h"
#include "../Runtimecounter.h"
//#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
using namespace std;

//=========================================================================
// Existing Graph Input Files:
//  - ContentToy.txt
//  - ContentToyPlus36.txt
//  - KDDToy_k6.txt
//  - ZKC.txt
//
// *** Note: Use '-o' to indicate step-by-step output
// *** Note: Use '-c' to indicate content column
// *** Note: Use '-d' to indicate dictionary column
// *** Note: Use '-gt' to indicate there is a ground truth file (Format: NodeNumber \t Community Number)
// *** Note: Use '-df' to indicate the dictionary filename
//=========================================================================

double getentropy(double prob){
	double res=0;
	if(prob<=0)
		return 0;
	res=-prob*log2(prob);
	return res;
}

double *pam;
double *HPi;
double *pa;
double *qi;
int *Ni;

void computeEntropyPi(int partid, int N, int *partition, int numpart){ 
	//compute H(Pi) for for all the modual 
	int m=0;
	double prob=0;
	double sumqp=0;
	for(m=0;m<numpart;m++){
		partid=m;
		HPi[partid]=0;
		sumqp=qi[partid]+pam[partid];
		if(sumqp==0)
			prob=0;
		else
			prob=qi[partid]/sumqp;
		HPi[partid]=getentropy(prob);
	}
	int i=0;
	for(i=0;i<N;i++){
		partid=partition[i];
		sumqp=qi[partid]+pam[partid];
		if(sumqp==0)
			prob=0;
		else
			prob=pa[i]/sumqp;
		HPi[partid]+=getentropy(prob);
	}
}

void InitNi(int *partition, int N){
    /*
    Init pam and Ni
    */
	int i=0;
	int pid;
	for(i=0;i<N;i++){
		pam[i]=0;
		Ni[i]=0;
	}
	for(i=0;i<N;i++){
		pid=partition[i];
		if(pid<0)
			continue;
		pam[pid]+=pa[i];
		Ni[pid]++;
	}
}

double structlength(Entry *&G, int *partition, int N, int numpart, double &HP, double &HQ){
	int i=0;
	int j=0;
	double qsum=0;
	int pid;
	int qid;
	int nv;
	HP=0;
	HQ=0;
	//qi[i] is coresponding to 
	///tau*(n-ni)/(n-1)\sum_{u\in part[i]}pa[u]+(1-tau)\sum_{u \in part[i]}{v\not\in part[i]}pa[u]w[u][v]
	for(i=0;i<numpart;i++){
		qi[i]=TAU*(N-Ni[i])*pam[i]/(N-1);
	}
	for(i=0;i<N;i++){ //linear time to graph size
		pid=partition[i];
		for(j=0;j<G[i].key;j++){
			nv=G[i].nbv[j];
			qid=partition[nv];
			if(pid!=qid){
				qi[pid]+=pa[i]*G[i].weight[j]*(1-TAU);
			}
		}
	}
	//pi[i] pi[i]=q[i]+sum_{u\in part[i]}pa[u]
	qsum=0;
	for(i=0;i<numpart;i++){
		qsum+=qi[i];
	}
	double res=0;
	double prob1;
	computeEntropyPi(i,N,partition,numpart); //linear time
	for(i=0;i<numpart;i++){
		if(qsum==0)
			prob1=0;
		else
			prob1=qi[i]/qsum;
		HQ+=getentropy(prob1);
		HP+= (qi[i]+pam[i])*HPi[i];
	}
	HQ*=qsum;
	//cout<<"left equation is "<<HQ<<endl;
	//cout<<"right equation is "<<HP<<endl;
	res=HQ+HP;
	return res;
}

void exit_with_help(char *name){
	 cout<<"usage:"<<endl;
	 cout<<name<<" graphfilename";
	 cout<<"optional:"<<endl;
	 cout<<"Use '-o' to indicate step-by-step output"<<endl;
     cout<<"Use '-of' to indicate a different output folder"<<endl;
	 cout<<"Use '-c' to indicate content column"<<endl;
	 cout<<"Use '-d' to indicate dictionary column"<<endl;
	 cout<<"Use '-gt' to indicate there is a ground truth file (Format: NodeNumber \t Community Number)"<<endl;
	 cout<<"Use '-df' to indicate the dictionary filename"<<endl;
     cout<<"Use '-of' to indicate the output folder"<<endl;
     cout<<"Use '-suffix' to indicate an output file suffix"<<endl;
     exit(1);
}

int main(int argc, char *argv[]){
	Runtimecounter RC;
	RC.start();
	if(argc < 2)
		exit_with_help(argv[0]);
	int OUTPUT_STEPS = 0;
    int ADD_CONTENT = 0;
    int ADD_DICTIONARY = 0;
    int GROUND_TRUTH = 0;

    string OUTFILE_NAME = "result_";
    string original_name = string(argv[1]);
	string pagefile;
    unsigned slashpos = original_name.rfind('/');
    if (slashpos != string::npos) {
        OUTFILE_NAME += original_name.substr(slashpos+1, original_name.length());
    }
    else {
        OUTFILE_NAME += original_name;
    }
    unsigned dotpos = OUTFILE_NAME.rfind('.');
    if (dotpos != string::npos) {
        OUTFILE_NAME = OUTFILE_NAME.substr(0, dotpos);
    }

    // If reweighted graph is to be used, we check to see if the file exists
    //int NoFileFlag = 0;
    // This is the file that contains the dictionary
    std::string dictionaryFile;
    // Commandline input arguments
	int i;
	int j;
	int m1;
	int k;
	int neigh;
	int e1;
	int neigh2;
	int m2;
    if (argc>2){
        for (i=2; i<argc; i++){
            std::string val_i = argv[i];
            if(val_i == "-o")
                OUTPUT_STEPS = 1;
            else if(val_i == "-c")
                ADD_CONTENT = 1;
            else if(val_i == "-d")
                ADD_DICTIONARY = 1;
            else if(val_i == "-dir")
                POWERMETHOD = 1;
			else if(val_i == "-tau"){
				TAU=atof(argv[i+1]);
			    i++;
			}
			else if(val_i == "-of"){
				OUTFILE_NAME = string(argv[i+1]) + "/" + OUTFILE_NAME;
			    i++;
			}
			else if(val_i == "-suffix"){
				OUTFILE_NAME = OUTFILE_NAME + string(argv[i+1]);
			    i++;
			}
            else if(val_i == "-gt")
                GROUND_TRUTH = 1;
			else if(val_i=="-df"){
                dictionaryFile=argv[i+1];
				i++;
			}else if(val_i=="-p"){
				pagefile=string(argv[i+1]);
				i++;
			}
        }
    }
	if (POWERMETHOD == 0)
		TAU=0;
	//Initialize memory blocks (START)
	curBlk2=0;
	memBlkAr2[0]=curMemPos2=(char*)malloc(BLK_SZ2);
	curMemEnd2=curMemPos2+BLK_SZ2;
	//Initialize memory blocks (END)
	
    // This is the file that contains the edges in the format: Node 1, Node 2, Weight
    std::string graphFile = argv[1];
    
    std::string groundTruthFilename;
	// Ground Truth File
        if(GROUND_TRUTH){
            groundTruthFilename = graphFile.substr(0, graphFile.size()-4) + "_GroundTruth.txt";
        }
    // Get the graph information
    int N=0;
	Entry *G;
	//load the graph into memory
	readgraph(graphFile.c_str(), N, G);
	// Output final graphFile Name
    cout << graphFile << endl;
    // Initialize vectors and matrices
    pa = new double[N];
	HPi=new double[N];
	//x is the pa (pagerank: expected frequence distribution)
    //computepagerank(G, N, pa);
	readpagerank(pagefile.c_str(),N,pa);
    double sum_p_logp = 0;
    for(i=0; i<N; i++){
		if(pa[i]>0)
			sum_p_logp -= pa[i]*log(pa[i])/log(2.0);
    }
	Entry *W=G;
    // Row normalizes W to sum to one
    normalizegraph(W, N);
    //----------------------------------------------------------------/
    // Load in the content files
    // Create a dictionary of module dictionaries
    ModuleDictionary * Dictionaries=NULL;
    if (ADD_CONTENT || ADD_DICTIONARY){
		Dictionaries = new ModuleDictionary [N];
        int nodeNum = 0;
        int dictionaryLength1 = 0;
		int windex;
		double wweight;
        std::ifstream fin(dictionaryFile.c_str());
        
        // Check to see if the dictionary file name is good
        if (!fin.good()) {
            std::cout << "FAILED\n" << dictionaryFile;
            exit(1);
        }
        
        // Load in the initial dictionaries for the nodes
        while(true){
            
            // Node Number for the dictionary
            fin >>nodeNum ;
            if(fin.eof())
                break;

            // Number of words in the dictionary
            fin >> dictionaryLength1 ;
            if(dictionaryLength1>0){
				Dictionaries[nodeNum].InitDictionary(dictionaryLength1);
				for(j=0; j<dictionaryLength1; j++){
                // Word
					fin >> windex;
                // Frequency (sums to one for each word)
					fin >> wweight;
					Dictionaries[nodeNum].addwords(windex,wweight,j);
				}
			}else{
				Dictionaries[nodeNum].set_words(NULL);
				Dictionaries[nodeNum].set_wordFrequency(NULL);
				Dictionaries[nodeNum].set_numWords(0);
				Dictionaries[nodeNum].set_size(0);
			}
			Dictionaries[nodeNum].normalize();
			Dictionaries[nodeNum].set_p_module(pa[nodeNum]);
			//if(nodeNum==0)
				//Dictionaries[nodeNum].printDictionary();
		}
		cout<<"finish loading dictionary"<<endl;
    }
    

    //----------------------------------------------------------------/
	
	///////////////////////////////////////////////////////////////////////////////////////
	//////////////initializing partitioning (start)////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // Begin the merging of nodes into modules
    int* partition = new int[N];
    int numPartitions = N;

    // Initialize partition vector - each node in its own community
    for(i=0;i<N;i++){
        partition[i] = i;
    }
    
    // Exit probabilities
    qi = new double[numPartitions];
	//q is qil in the paper
    
    // Number of nodes in module i
    Ni= new int[numPartitions];

    // Initialize to zero
    for(i=0;i<numPartitions;i++){
        qi[i] = 0;
        Ni[i] =0;
    }
    
    // Count the number of elements in each partition
    for(i=0;i<N;i++){
        Ni[partition[i]] +=1;
    }
    
    // Initialize expressions for L
    double Q = 0; //ql
    double Q_log = 0; //qlogq
    double PQ_log = 0;//pil log pil in paper Appendix in Equation 4
    pam = new double[numPartitions];  // pam gives the sum_(alpha in i) p_i
    //the right part of Equation 4 in paper Appendix
	
	
	// Calculate q and pam with every node being in its own module
	for(i=0;i<numPartitions;i++)
		pam[i]=0;
    for(i=0;i<N;i++){
		pam[partition[i]]+=pa[i];
		qi[partition[i]]+=(TAU*(N-Ni[partition[i]])/(N-1))*pa[i];
		for(int j=0;j<W[i].key;j++){
			neigh=W[i].nbv[j];
			if(partition[neigh]!=partition[i])
				qi[partition[i]]+=(1-TAU)*pa[i]*W[i].weight[j];
		}
	}
	for(m1=0;m1<numPartitions;m1++){
		Q += qi[m1];
         if (qi[m1]>0)
             Q_log += qi[m1]*log(qi[m1])/log(2.0);
         if (qi[m1]+pam[m1]>0)
             PQ_log += (qi[m1]+pam[m1])*log(qi[m1]+pam[m1])/log(2.0);
	}
	
	
    // Initialize minL and minL_original
    double minL = Q*log(Q)/log(2.0) -2.0*Q_log + PQ_log+sum_p_logp;
	cout<<minL<<endl;
    double contentTotal = 0;
    double dictionaryTotal = 0;
    
    // Calculate the content term and full dictionary term
    if(ADD_CONTENT){
        for(i=0;i<numPartitions;i++){
            contentTotal +=   Dictionaries[i].calcContentEntropy(); // Individual dictionaries
        }
        minL += contentTotal;
    }
    if(ADD_DICTIONARY){
        for(i=0;i<numPartitions;i++){
            dictionaryTotal += Dictionaries[i].calcDictionaryEntropy(); // Total dictionary
        }
        minL += dictionaryTotal;
    }
    // Merging the two terms since the change in content/dictionaries is performed together
    contentTotal += dictionaryTotal;
    
    
    //------------------------------------------------------------------
    // Output initial network description length  (start)
    //------------------------------------------------------------------
    
    // Calculate the initial network description length
    double left_column0 = 0;
    double right_column0 = 0;
    //double initial_q_sum = 0;
    double content_column0 = 0;
    double dictionary_column0 = 0;
	double initL=structlength(W, partition, N, numPartitions, left_column0, right_column0);
/*     for(m1=0;m1<numPartitions;m1++){
        initial_q_sum += q[m1];
    } */
    
   /*  // Calculate Within Module (RightColumn) Contribution and Between Module (LeftColumn) Contribution
    for(m1=0;m1<numPartitions;m1++){
        if (initial_q_sum > 0 && q[m1]/initial_q_sum > 0){
            left_column0  -= q[m1]*log(q[m1]/initial_q_sum)/log(2.0);
        }
        if ((q[m1]+pam[m1]) > 0){
            right_column0 -= q[m1]*log(q[m1]/(q[m1]+pam[m1]))/log(2.0);
        }
    }
	for(i=0;i<N;i++){
		m1=partition[i];
		if ((q[m1]+pam[m1]) > 0){
			right_column0 -=x[i]*log(x[i]/(q[m1]+pam[m1]))/log(2.0);
		}
	} */
    
    // Calculate Content Contribution
    if(ADD_CONTENT || ADD_DICTIONARY){
        for(m1=0;m1<numPartitions;m1++){
            if(ADD_CONTENT){
                content_column0 +=Dictionaries[m1].calcContentEntropy();
            }
            if(ADD_DICTIONARY){
                dictionary_column0 += Dictionaries[m1].calcDictionaryEntropy();
            }
        }
    }
    std::cout<< "\n";
    
    // Output the initial  values
    std::cout << endl << "\n\n" << setw(13) << "Left Column:" << setw(16)<< "Right Column:";
    std::cout << setw(20)<< "Content Column:"<< setw(20)<< "Dictionary Column:";
    std::cout << setw(13) <<"Combined:\n";
    
    std::cout << endl  << setw(13) << left_column0 << setw(16) << right_column0;
    std::cout << setw(20)<< content_column0 << setw(20)<< dictionary_column0;
    std::cout << setw(13)<< left_column0+right_column0+content_column0 + dictionary_column0;
    std::cout << endl << endl;
    cout<<minL<<"\t"<<initL<<endl;
    
    //------------------------------------------------------------------
    // End initial description length output    (end)
    //------------------------------------------------------------------
    
	///////////////////////////////////////////////////////////////////////////////////////
	//////////////initilizing partitioning (end)//////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // Initialize storage variables
    double q_i_store = 0;
    double Q_store = 0;
    double Q_log_store = 0;
    double PQ_log_store = 0;
    double contentStore = 0;
    int m1_store=-1, m2_store=-1;
	double Lm12;
	double desccodinglength=0;
	double preL; //the previous iteration L
	preL=minL;
	char *isaccessed=new char[N]; //isaccessed is used to mark whether a modual has been visited or not;
	//psumpart[i] is the sum of x[node] for all the nodes that are in partition i;
	int iterations = 1;
	//interparttion is used to stor interpartition edges;
	Modual *interparttion=new Modual[N];
	double q_sum_store;
	tempWords=(int*)malloc(sizeof(int)*tempsize);
	tempFrequencies=(double*)malloc(sizeof(double)*tempsize);
	do{
		//initilize interparttion edges;
        m1_store=-1;
        m2_store=-1;
		for(i=0;i<N;i++){
			for(j=0;j<W[i].key;j++){
				neigh=W[i].nbv[j];
				if(partition[i]!=partition[neigh]){
					//put neigh, G[i].weight[j] into Modual[partition[i]];
					interparttion[partition[i]].addedge(i, neigh,W[i].weight[j]);
				}
			}
		}
		//////////////////////////////////////////////////////////////////////
		//iterate through all possible combination of m1, m2
		//to search for the best combination of (m1, m2)   (start)////////////
		//////////////////////////////////////////////////////////////////////
		double p_sum_i = 0;
		double p_sum_log = 0;    // Local version of PQ_log
		double q_sum = 0;        // Local version of Q
		double q_sum_log = 0;    // Local version of Q_log
		double contentLocal = 0; // Local version of contentTotal + DictionaryTotal
		for(m1=0;m1<numPartitions;m1++){ 
			memset(isaccessed,'f',sizeof(char)*N);
			//if(iterations>=2474&&m1==782)
				//cout<<"m1 value "<<m1<<endl;
			//cout<<"OK bool here"<<endl;
			m2=0;
			for(j=0;j<interparttion[m1].size;j++){
				//if(iterations>=2474&&m1==782&&j==30)
					//cout<<"j value "<<j<<endl;
				neigh=interparttion[m1].nbv[j];
				m2=partition[neigh];
				if(isaccessed[m2]=='f'){
					isaccessed[m2]='t';
                    //cout<<iterations<<"\t"<<m1<<"\t"<<m2<<endl;
					// Initialize to zero            
					p_sum_i = 0;
					p_sum_log = 0;    // Local version of PQ_log
					q_sum = 0;        // Local version of Q
					q_sum_log = 0;    // Local version of Q_log
					contentLocal = 0; // Local version of contentTotal + DictionaryTotal
					p_sum_i+=pam[m1];
					p_sum_i+=pam[m2];
					q_sum+=qi[m1];
					q_sum+=qi[m2];
					for(k=0;k<interparttion[m1].size;k++){
						e1=interparttion[m1].src[k];
						neigh2=interparttion[m1].nbv[k];
						if(partition[neigh2]==m2)
							q_sum-=interparttion[m1].weights[k]*pa[e1];
					}
					for(k=0;k<interparttion[m2].size;k++){
						e1=interparttion[m2].src[k];
						neigh2=interparttion[m2].nbv[k];
						if(partition[neigh2]==m1)
							q_sum-=interparttion[m2].weights[k]*pa[e1];
					}
					// Calculate expressions for L
					q_sum_store = q_sum;
					if(q_sum > 0){
						q_sum_log = q_sum*log(q_sum)/log(2.0);
					}
					if (q_sum + p_sum_i > 0){
						p_sum_log = (q_sum + p_sum_i)*log(q_sum+p_sum_i)/log(2.0);
					}

					// Don't reloop through everything.  Instead, update only what changed.
					q_sum += Q-qi[m1]-qi[m2];
					q_sum_log += Q_log-qi[m1]*log(qi[m1])/log(2.0)-qi[m2]*log(qi[m2])/log(2.0);
					p_sum_log += PQ_log-(qi[m1] + pam[m1])*log(qi[m1] + pam[m1])/log(2.0)-(qi[m2] + pam[m2])*log(qi[m2] + pam[m2])/log(2.0);

					
					// Calculating L given a merge of m1 and m2
					Lm12 = q_sum*log(q_sum)/log(2.0) - 2.0*q_sum_log + p_sum_log+sum_p_logp;
					
					if(ADD_CONTENT && ADD_DICTIONARY){
						contentLocal = contentTotal+deltaDictionaryEntropy(Dictionaries[m1], Dictionaries[m2]);
						Lm12 += contentLocal;
					}
					else if(ADD_CONTENT){
						//if(iterations>=438)
							//cout<<"OK in line 387"<<endl;
						contentLocal = contentTotal+deltaContentOnlyEntropy(Dictionaries[m1], Dictionaries[m2]);
						Lm12 += contentLocal;
						//if(iterations>=438)
							//cout<<"OK in line 387"<<endl;
					}
					else if(ADD_DICTIONARY){
						contentLocal = contentTotal+deltaDictionaryOnlyEntropy(Dictionaries[m1], Dictionaries[m2]);
						Lm12 += contentLocal;
					}
					// Store if it gives a better measure
					if(Lm12-minL<0.00000001 ){
						//cout<<"come here"<<endl;
						minL =Lm12;
						m1_store = m1;
						m2_store = m2;
						q_i_store = q_sum_store;
						Q_store = q_sum;
						Q_log_store = q_sum_log;
						PQ_log_store = p_sum_log;
						if(ADD_CONTENT || ADD_DICTIONARY){
							contentStore = contentLocal;
						}
					}
					
				}//merge m1 and m2;
			}
		}

		//////////////////////////////////////////////////////////////////////
		//iterate through all possible combination of m1, m2
		//to searh for the best combination of (m1, m2)   (end)/////////////
		//////////////////////////////////////////////////////////////////////
		
		
		////////////////////////////////////////////////////////
		//merge m1_store and m2_store  (start)/////////////////
		//////////////////////////////////////////////////////
		  // Output the intial module values
		  if(OUTPUT_STEPS){
			
           //cout << endl << "\n\nModule\tq_i\tpam_i:\n\n\n";
           //for (i=0; i<numPartitions; i++)
           //    cout << endl  << i << "\t" << qi[i] << "\t" << pam[i];
           //cout << endl << endl;
          }
        //cout<<m1_store<<"\t"<<m2_store<<endl;
		  if(m1_store==-1||m2_store==-1)
			  break;
        // Merge the two selected modules
        int m_min = min(m1_store, m2_store);
        int m_max = max(m1_store, m2_store);
		//cout<<"OK before merging dict"<<endl;
        if (ADD_CONTENT || ADD_DICTIONARY){
			/*if(iterations>438)
				cout<<"OK in line 444"<<endl;*/
            mergeDictionaries(Dictionaries[m_min], Dictionaries[m_max]);
			/*if(iterations>438)
				cout<<"OK in line 447"<<endl;*/
        }
		////cout<<"OK merge dictionaries"<<endl;

        for (i=0; i<N; i++){
            if (partition[i] == m1_store || partition[i]==m2_store){
                partition[i] = m_min;
            }
            else if (partition[i]>m_max){
                partition[i] -=1;
            }
        }
		
		////////////////////////////////////////////////////////
		//merge m1_store and m2_store  (end)///////////////////
		//////////////////////////////////////////////////////
		
		
		//remember to release the memory interpartition later;
		for(i=0;i<numPartitions;i++){
				if(i<(numPartitions-1))
					interparttion[i].reset();
				else
					interparttion[i].clear();
		}

		//cout<<"OK release memory\n"<<endl;
		///////////////////////////////////////////////////
		////update information after merge (start)//////////
		////////////////////////////////////////////////////
		// Update vectors n, q, pam, Dictionaries
		//cout<<m_min<<" "<<m1_store<<" "<<m2_store<<endl;
        Ni[m_min] = Ni[m1_store] + Ni[m2_store];
        qi[m_min] = q_i_store;
        pam[m_min] = pam[m1_store] + pam[m2_store]; //check the value here, it is not correct
        for (i=m_max; i<(numPartitions-1); i++){
            Ni[i] = Ni[i+1];
            qi[i] = qi[i+1];
            pam[i] = pam[i+1];
            if (ADD_CONTENT || ADD_DICTIONARY){
				/*if(iterations>438)
					cout<<"OK in line 484"<<endl;*/
                replaceDictionaries(Dictionaries[i], Dictionaries[i+1]);
			    /*if(iterations>438)
					cout<<"OK in line 487"<<endl;*/
            }
        }
		
	/*	for(i=m_max;i<numPartitions;i++){
			for(j=0;j<Dictionaries[i].get_numWords();j++){
				cout<<" "<<Dictionaries[i].get_words()[j];
			}
			cout<<endl;
		}*/
        // Update values to stored from previous iteration
        Q = Q_store;
        Q_log = Q_log_store;
        PQ_log = PQ_log_store;
        if(ADD_CONTENT || ADD_DICTIONARY){
            contentTotal = contentStore;
        }
        // Update minL and store its old value
        minL = Q*log(Q)/log(2.0) -2.0*Q_log + PQ_log+sum_p_logp;
        
        if(ADD_CONTENT || ADD_DICTIONARY){
            minL += contentTotal;  // Content contribution + Full Dictionary contribution
        }
		// Update number of partitions
        numPartitions -=1;
        //cout<<numPartitions<<endl;
		desccodinglength=preL-minL;
		preL=minL;
		iterations++;
		///////////////////////////////////////////////////
		////update information after merge (start)//////////
		////////////////////////////////////////////////////
		//if(iterations>=2474)
			//cout<<"iteration number "<<iterations<<endl;
		//cout<<desccodinglength<<endl;
	}while(desccodinglength>0.00000001&&numPartitions>1);
	cout<<"finish partitioning"<<endl;
    // Output the final partition
	//string outfile=graphFile.substr(0, graphFile.size()-4);
    string outfile=OUTFILE_NAME;
	//if(ADD_CONTENT &&ADD_DICTIONARY)
	//	outfile=outfile+"_CD";
	//else
	//	if(ADD_CONTENT)
	//		outfile=outfile+"_C";
	//else
	//	if(ADD_DICTIONARY)
	//		outfile=outfile+"_D";
    //outfile=outfile+ "_smap.txt";
    outfile=outfile+".txt";
    // cout << endl << "Communities for " << graphFile << ":\n";
    ofstream fout(outfile.c_str());
    for (i=0; i<N; i++){
        //if (x[i]>0)
            fout << i << "\t" << partition[i]<<endl;
    }
	fout.close();
    //cout << endl << endl;
    cout<<numPartitions<<endl;
    
     //Output the final module values
     //std::cout << endl << "Module\tq_i\tpam_i:\n";
     //for (i=0; i<numPartitions; i++){
     //    if (qi[i] > 0)
     //       std::cout << i << "\t" << qi[i] << "\t" << pam[i]<<endl;
     //}
    
    // Calculate the final values
    double left_column = 0;
    double right_column = 0;
	InitNi(partition,N);
    //for (i=0; i<numPartitions; i++){
        //if (qi[i] > 0)
            //std::cout << i << "\t" << qi[i] << "\t" << pam[i]<<endl;
    //}
    structlength(W, partition, N, numPartitions, left_column, right_column);
    //for (i=0; i<numPartitions; i++){
    //    printf
    //    cout << qi[i] << endl;
    //    if (qi[i] > 0)
    //        std::cout << i << "\t" << qi[i] << "\t" << pam[i]<<"\t"<<HPi[i]<<endl;
    //}
    // This calculates the final content and full dictionary contributions
    double content_column = 0;
    double dictionary_column = 0;
    if (ADD_CONTENT){
        for(m1=0;m1<numPartitions;m1++){
            content_column +=Dictionaries[m1].calcContentEntropy();
        }
    }
    if (ADD_DICTIONARY){
        for(m1=0;m1<numPartitions;m1++){
            dictionary_column += Dictionaries[m1].calcDictionaryEntropy();
        }
    }
    
    // Individual Module Dictionary Contributions
    if ( OUTPUT_STEPS && (ADD_CONTENT || ADD_DICTIONARY)){
        std::cout << "\n\nDictionaries:\n";
        for(int m1=0;m1<numPartitions;m1++){
            std::cout << "\n\n   Module:  " << m1;
            std::cout <<"   -------   Dictionary Contribution:  " <<  Dictionaries[m1].calcContentEntropy() << "\n";
            Dictionaries[m1].printDictionary();
        }
        std::cout<< "\n";
    }
    
    // Output the final  values
    std::cout << endl << "\n\n" << setw(13) << "Left Column:" << setw(16)<< "Right Column:"<< setw(20)<< "Content Column:"<< setw(20)<< "Dictionary Column:" << setw(13) <<"Combined:\n";
    std::cout << endl  << setw(13) << left_column << setw(16) << right_column << setw(20)<< content_column<< setw(20)<< dictionary_column << setw(13)<< left_column+right_column+content_column+dictionary_column;
    std::cout << endl << endl;
    
    std::cout << endl << "Iterations:\t" << iterations << endl << endl;
    
    
// --------------- Outputting Results to Files -------------------//
    
    // Filenames
    string outFile, outFileNewNums,outputResultsFilename, outputMatlabFilename;

    if(ADD_CONTENT || ADD_DICTIONARY){

        
        
        // File Names for output / Matlab code
        outFile = OUTFILE_NAME + "_CommunityResults_WithContent.txt";
        outFileNewNums = OUTFILE_NAME + "_CommunityResults_WithContent_NewNums.txt";
        outputResultsFilename = OUTFILE_NAME + "_CommunityResults_Classes_WithContent.txt";
        outputMatlabFilename = OUTFILE_NAME + "_CommunityResults_Classes_WithContent.m";
    }
    else{
        
        // File Names for output / Matlab code
        outFile = OUTFILE_NAME + "_CommunityResults.txt";
        outFileNewNums = OUTFILE_NAME + "_CommunityResults_NewNums.txt";
        outputResultsFilename = OUTFILE_NAME + "_CommunityResults_Classes.txt";
        outputMatlabFilename = OUTFILE_NAME + "_CommunityResults_Classes.m";
    }
    
    // Output the final community partitions
    //outputCommunities(outFile, x, partition, N);
    
    // Output the final community partitions with new module numbers that exclude singleton modules
    //outputCommunities_newModuleNumbers(outFileNewNums, x, partition, N);
    
    // If there is a ground truth file, then 
    if (GROUND_TRUTH){
        // Nodes grouped by communities, placed with ground truth
        compareToGroundTruth(groundTruthFilename, outputResultsFilename, pa, partition,  N);
        
        // Create a matlab file that will plot the communities
        outputCommunities_MatlabFile(groundTruthFilename, outputMatlabFilename, pa, partition,  N, G);
    }
    delete []pa;
	delete []qi;
	delete []pam;
	delete []isaccessed; //isaccessed is used to mark whether a modual has been visited or not;
	delete []Ni;
	delete []partition;
	delete []HPi;
	while(curBlk2 >= 0){
		if(memBlkAr2[curBlk2]!=NULL){
			free(memBlkAr2[curBlk2]);
			curBlk2--;
		}else
			curBlk2--;
	}
	if(Dictionaries!=NULL){
		for(i=0;i<N;i++){
			Dictionaries[i].clear();
		}
		delete []Dictionaries;
		Dictionaries=NULL;
	}
	if(interparttion!=NULL){
		for(i=0;i<N;i++){
			interparttion[i].clear();
		}
		delete []interparttion;
		interparttion=NULL;
	}
	if(G!=NULL){
		for(i=0;i<N;i++){
			if(G[i].nbv!=NULL)
				free(G[i].nbv);
			if(G[i].weight!=NULL)
				free(G[i].weight);
		}
		free(G);
		G=NULL;
	}
	if(tempWords!=NULL){
		free(tempWords);
		tempWords=NULL;
	}
	if(tempFrequencies){
		free(tempFrequencies);
		tempFrequencies=NULL;
	}
	RC.stop();
	//cout<<"total running time is "<<RC.GetRuntime()<<"  seconds"<<endl;
	//_CrtDumpMemoryLeaks();
    return 0;
}