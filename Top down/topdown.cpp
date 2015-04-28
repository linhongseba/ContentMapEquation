#include "partition.h"
#include "../Runtimecounter.h"
#include "../graph_common.h"

void exit_with_help(char *name){
	 cout<<"usage:"<<endl;
	 cout<<name<<" graphfilename";
	 cout<<"optional:"<<endl;
	 //cout<<"Use '-o' to indicate step-by-step output"<<endl;
     cout<<"Use '-of' to indicate the output folder"<<endl;
     cout<<"Use '-suffix' to indicate an output file suffix"<<endl;
	 cout<<"Use '-c' to indicate content column"<<endl;
	 cout<<"Use '-d' to indicate dictionary column"<<endl;
	 cout<<"Use '-dir' to indicate a directed graph (default is undirected)"<<endl;
	 cout<<"Use '-tau <tau>' to indicate the teleportation probability (only used for directed graphs, default is 0.15)"<<endl;
	 cout<<"Use '-g' to indicate there is a ground truth file (Format: NodeNumber \t Community Number)"<<endl;
	 cout<<"Use '-f' to indicate the dictionary filename"<<endl;
	 cout<<"user '-e' to indicate extended content column"<<endl;
     cout<<"user '-trials <trials>' to return the best partitioning of <trials> trials"<<endl;
	 exit(1);
}
int main(int argc, char*argv[]){
 // Commandline input arguments
	if(argc<2)
		exit_with_help(argv[0]);
	Runtimecounter RC;
	RC.start();
	srand((unsigned)time(NULL));

	//int OUTPUT_STEPS = 0;
    int ADD_CONTENT = 0;
    int ADD_DICTIONARY = 0;
    //int GROUND_TRUTH = 0;
    string OUTFILE_NAME = "result_";
    int TRIALS = 1;

    string original_name = string(argv[1]);
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

	int add_e=0;
	int i=2;
	//Initialize memory blocks (START)
	curBlk2=0;
	memBlkAr2[0]=curMemPos2=(char*)malloc(BLK_SZ2);
	curMemEnd2=curMemPos2+BLK_SZ2;
	//Initialize memory blocks (END)

	char *graphfile;
	char *dictfile;
	allocatepermemory(sizeof(char)*300);
	graphfile=curMemPos2;
	curMemPos2+=(sizeof(char)*150);
	dictfile=curMemPos2;
	curMemPos2+=(sizeof(char)*150);
	strcpy(graphfile,argv[1]);
	int len=strlen(argv[1])-4;
	strncpy(dictfile,argv[1],len);
	dictfile[len]='\0';
	strcat(dictfile,"_Dictionary.txt");
	int k=10;
	if (argc>2){
        for(i=2;i<argc;i++){
            std::string val_i = argv[i];
			    if (val_i == "-c")
				    ADD_CONTENT =1;
			    //else if (val_i == "-o")
				//     OUTPUT_STEPS=1;
			    else if (val_i == "-d")
				    ADD_DICTIONARY = 1;
                else if(val_i == "-dir")
                    POWERMETHOD = 1;
			    else if(val_i == "-tau"){
				    TAU=atof(argv[i+1]);
			        i++;
			    }
			    //else if (val_i == "-g")
				//    GROUND_TRUTH = 1;
			    else if (val_i == "-e")
				    add_e=1;
                else if (val_i == "-f"){
				    if(i>=argc-1)
					    exit_with_help(argv[0]);
				    strcpy(dictfile,argv[i+1]);
				    i++;
                }
                else if (val_i == "-k"){
				    if(i>=argc-1)
					    exit_with_help(argv[0]);
				    k=atoi(argv[i+1]);
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
			    else if(val_i == "-trials"){
			    	TRIALS = atoi(argv[i+1]);
			        i++;
			    }
                else
				    exit_with_help(argv[0]);
        }
	}
    if (POWERMETHOD == 0)
        TAU = 0;
	int N;
	Entry *G=NULL;
	Entry *Dict=NULL;
	readgraph(graphfile,N,G);
	pa=new double[N];
	computepagerank(G,N,pa);
	normalizegraph(G,N);
	if(ADD_CONTENT||ADD_DICTIONARY||add_e){
		tfsize=0;
		readdictionary(dictfile,N,Dict,tfsize);
		tfsize++;
		normalizegraph(Dict,N);
	}
	temptf=new double[tfsize+1];
	int *partid = NULL;
	//int numpart;
    double mdl=1000000;
	ccid=new int[N];
	if(ccid==NULL){
		printf("system could not allocate more memory\n");
		exit(0);
	}
	cout<<endl;
	cout<<argv[1]<<endl;
    double cur_mdl;
    for(int i=0;i<TRIALS;i++){
        cout<<"    trial #"<<i<<endl;
        int *cur_partid=(int *)malloc(sizeof(int)*N);
	    int cur_numpart=0;
        cur_mdl=SCpartition(G,Dict,N,cur_partid,cur_numpart,ADD_CONTENT,ADD_DICTIONARY,add_e,k);
        if (cur_mdl<mdl){
            if(partid!=NULL)
		        free(partid);
            partid=cur_partid;
            //numpart=cur_numpart;
            mdl=cur_mdl;
        }
        else{
            free(cur_partid);
        }
    }

    cout<<endl<<"final minimum codelength: "<<mdl<<endl;

    // clean up
    delete []ccid;

	//char *outfile;
	string outfile=OUTFILE_NAME;
	allocatepermemory(sizeof(char)*450);
	curMemPos2+=(sizeof(char)*150);
	curMemPos2+=(sizeof(char)*150);
	curMemPos2+=(sizeof(char)*150);
	//if(ADD_CONTENT&&ADD_DICTIONARY)
	//    outfile=outfile+"_comb";
	//else
	//	if(ADD_CONTENT)
	//        outfile=outfile+"_cont";
	//	else
	//		if(ADD_DICTIONARY)
	//			outfile=outfile+"_dict";
	//		else
	//			if(add_e)
	//				outfile=outfile+"_e";
	//			//else
	//				//outfile=outfile+"_s";
    string groundfile = outfile + "_GroundTruth.txt";
	string wfile = outfile + ".m";
	outfile=outfile+".txt";
	ofstream fout(outfile.c_str());
	for(i=0;i<N;i++){
		fout<<i<<"\t"<<partid[i]<<endl;
	}
	fout.close();
	//outputCommunities_MatlabFile(groundfile,wfile,partid,numpart,N,graphfile);
	delete []qi;
	delete []HPi;
	delete []pam;
	delete []Ni;
	if(partid!=NULL)
		free(partid);
	if(pa!=NULL)
		delete []pa;
	if(G!=NULL){
		for(i=0;i<N;i++){
			if(G[i].key>100000){
				if(G[i].nbv!=NULL){
					free(G[i].nbv);
				}
				if(G[i].weight!=NULL){
					free(G[i].weight);
				}
			}
		}
		free(G);
		G=NULL;
	}
	if(Dict!=NULL){
		for(i=0;i<N;i++){
			if(Dict[i].key>100000){
				if(Dict[i].nbv!=NULL){
					free(Dict[i].nbv);
				}
				if(Dict[i].weight!=NULL){
					free(Dict[i].weight);
				}
			}
		}
		free(Dict);
		Dict=NULL;
	}
	while(curBlk2 >= 0){
		if(memBlkAr2[curBlk2]!=NULL){
			free(memBlkAr2[curBlk2]);
			curBlk2--;
		}else
			curBlk2--;
	}
	curBlk=0;
	while(curBlk <= endBlk){
		if(memBlkAr[curBlk]!=NULL){
			free(memBlkAr[curBlk]);
			curBlk++;
		}else
			curBlk++;
	}
	RC.stop();
	//cout<<"total running time is\t"<<RC.GetRuntime()<<"\tsecond"<<endl;
	//cout<<argv[1]<<"\t"<<RC.GetRuntime()<<endl;
	return 1;
}