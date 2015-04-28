/* --------------------------------------------
 
 outputGraph.h
 Modified by Laura Smith on 6/11/13.
 
 Functions in this file are for outputting the 
 results of content_algorithm_withDictionary.cpp
 
 -------------------------------------------- */

#include<fstream>
#include<string>
#include<iostream>
#ifndef _outputGraph_h
#define _outputGraph_h

//  Output the communities: NodeNum \t CommunityNum \t p_alpha
void outputCommunities(std::string filename, double x[], int partition[], int N){
    std::ofstream fout(filename.c_str());

    for (int i=0; i<N; i++){
        if (x[i]>0)
           fout << i << "\t" << partition[i] << "\t" << x[i] << "\n";
    }
    return;
}

// Output the communities, excluding singleton modules: NodeNum \t CommunityNum \t p_alpha
void outputCommunities_newModuleNumbers(std::string filename, double x[], int partition[], int N){
    std::ofstream fout(filename.c_str());
    int *newPartition=new int[N];
    for (int i=0; i<N; i++){
        newPartition[i]  = -1;
    }
    int counter = 0;
    
    for (int i=0; i<N; i++){
        if (x[i]>0){
            if (newPartition[partition[i]] == -1){
                newPartition[partition[i]] = counter;
                counter ++;
            }
            fout << i << "\t" << newPartition[partition[i]] << "\t" << x[i] << "\n";
        }
    }
	delete []newPartition;
    return;
}

// Output the communities, excluding singleton modules. Compare to the ground truth communities.
void compareToGroundTruth(std::string groundTruthFilename, std::string outputResultsFilename, double x[], int partition[], int N){
    std::ifstream fin(groundTruthFilename.c_str());
    
    int *nodes=new int[N];
    int *stance=new int[N];
     int n1, s1;
    // Input user stance
    for (int k=0; k<N; k++){
        
       
        fin >> n1;
        fin >> s1;
        
        nodes[k] = n1;
        stance[k] = s1;
    }
    
    
    // Create new module numbers for elements where x[i]>0
    int *newPartition=new int[N];
    for (int i=0; i<N; i++){
        newPartition[i]  = -1;
    }
    int counter = 0;
    
    for (int i=0; i<N; i++){
        if (x[i]>0){
            if (newPartition[partition[i]] == -1){
                newPartition[partition[i]] = counter;
                counter ++;
            }
        }
    }
    
    std::ofstream fout(outputResultsFilename.c_str());

    // Group the communities together
    for (int c=0; c<counter; c++){
        fout << "Community:  " << c << "\n";
        for(int i=0; i<N; i++){
            // Output nodeNum \t Stance
            if (newPartition[partition[i]] == c){
                fout << nodes[i] << "\t" << stance[i] << "\n";
            }
        }
        fout << "\n";
    }
	delete []newPartition;
	delete []nodes;
	delete []stance;
    return;
}



// Create a MATLAB file that can be run to plot the communities with ground truth information
void outputCommunities_MatlabFile(std::string groundTruthFilename, std::string outputMatlabFilename, double x[], int partition[], int N,Entry* A){
    
    std::ifstream fin(groundTruthFilename.c_str());
    std::cout<< "\n"<< "\n" << groundTruthFilename << "\n"<< "\n";
    
    int *nodes=new int[N];
    int *stance=new int[N];
    int n1, s1;
    // Input user stance
    for (int k=0; k<N; k++){
        
       
        fin >> n1;
        fin >> s1;
        
        nodes[k] = n1;
        stance[k] = s1;
    }
    
    // Create new module numbers for elements where x[i]>0
    int *newPartition=new int[N];
    for (int i=0; i<N; i++){
        newPartition[i]  = -1;
    }
    int counter = 0;
    
    for (int i=0; i<N; i++){
        if (x[i]>0){
            if (newPartition[partition[i]] == -1){
                newPartition[partition[i]] = counter;
                counter ++;
            }
        }
    }
    
    std::ofstream fout(outputMatlabFilename.c_str());
    
    // Group the communities together - nodeNum+1 (for MATLAB) \t community \t stance
    fout << "Communities =  [\n";
    for (int c=0; c<counter; c++){
        
        for(int i=0; i<N; i++){
            if (newPartition[partition[i]] == c){
                fout << nodes[i]+1 << ",\t" << c << ",\t" << stance[i] << ";\n";
            }
        }
    }
    fout << "];\n\n\n";

    
    // Output Edge Network to Matlab
    fout << "A =  [\n";
    for(int i=0; i<N; i++){
		for(int j=0;j<A[i].key;j++){
			fout << A[i].weight[j] << "\t" ;
		}
        // for(int j=0; j<N; j++){
            // fout << A[i][j] << "\t" ;
        // }
        fout << ";\n";
    }
    fout << "];\n\n";

    // Store each of the columns of Communities Matrix
    fout << "[N,M] = size(Communities);\n";
    fout << "nodeNums = Communities(:,1);\n";
    fout << "classes = Communities(:,2);\n";
    fout << "groundTruth = Communities(:,3);\n\n";
    
    // Find the nodes with edges
    fout << "hasEdges = find(sum(A)>0);\n\n";
    fout << "A2 = A(hasEdges, hasEdges);\n";

    // Original Adjacency
    fout << "figure, spy(A2)\n\n";
    fout << "B = A(nodeNums, nodeNums);\n";
    
    // Reordered Adjacency according to communities
    fout << "figure, spy(B)\n\n";
    
    // Nodes separated by community with color indicating a position (-1,0,1,other)
    fout << "figure, hold on\n";
    fout << "for i=1:N\n";
    fout << "\tif groundTruth(i) == 0\n";
    fout << "\t\tplot( classes(i),nodeNums(i),'k.', 'MarkerSize', 20)\n";
    fout << "\telseif groundTruth(i) == 1\n";
    fout << "\t\tplot( classes(i),nodeNums(i),'g.', 'MarkerSize', 20)\n";
    fout << "\telseif groundTruth(i) == -1\n";
    fout << "\t\tplot( classes(i),nodeNums(i),'r.', 'MarkerSize', 20)\n";
    fout << "\telse\n";
    fout << "\t\tplot( classes(i),nodeNums(i),'b.', 'MarkerSize', 20)\n";
    fout << "\tend\n";
    fout << "end\n\n";
                
    fout << "xlabel('Class Number')\n";
    fout << "ylabel('Node Number')\n\n\n";

   return;
}



#endif
