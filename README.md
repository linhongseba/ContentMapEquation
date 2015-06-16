# ContentMapEquation
The implementation for the paper Laura M Smith, Linhong Zhu, Kristina Lerman, and Allon G Percus. Partitioning Networks with Node Attributes by Compressing Information Flow. In arXiv preprint arXiv:1405.4332.

#Installation
1.Requirement:

g++, gcc compiler

2.platform support:

Mac, Windows, Linux

3.Compiling:

go into either the Top down or Bottom up, type make to automatically compile and generate the executable files

#Usage
1. Bottom up approach

Usage: graphfilename [option]

 "Use '-o' to indicate step-by-step output";
	 
	 "Use '-of' to indicate the output folder"
	 
	 "Use '-suffix' to indicate an output file suffix";
	 
	 "Use '-c' to indicate content column";
	 
	 "Use '-d' to indicate dictionary column";
	 
	 "Use '-dir' to indicate a directed graph (default is undirected)";
	 
	 "Use '-tau <tau>' to indicate the teleportation probability (only used for directed graphs, default is 0.15)";
	 
	 "Use '-g' to indicate there is a ground truth file (Format: NodeNumber \t Community Number)";
	 
	 "Use '-df' to indicate the dictionary filename";
	 
2.  Top down approach

   usage: graphfilename [option]

	 "Use '-o' to indicate step-by-step output";
	 
	 "Use '-of' to indicate the output folder"
	 
	 "Use '-suffix' to indicate an output file suffix";
	 
	 "Use '-c' to indicate content column";
	 
	 "Use '-d' to indicate dictionary column";
	 
	 "Use '-dir' to indicate a directed graph (default is undirected)";
	 
	 "Use '-tau <tau>' to indicate the teleportation probability (only used for directed graphs, default is 0.15)";
	 
	 "Use '-g' to indicate there is a ground truth file (Format: NodeNumber \t Community Number)";
	 
	 "Use '-f' to indicate the dictionary filename";
	 
	 "user '-e' to indicate extended content column";
	 
	 "user '-trials <trials>' to return the best partitioning of <trials> trials";

#Input format 
    + format for graph

The first line is number of nodes, and starting from the second lins is the adjacence list of each node formated as follows:

node_id,degree_d:neighboreid1,weight1:neighborid2,weight2:...neighboridd,weightd

Note that the node_id is within the range [0,n-1], where n is number of nodes, and the list of neighbors are sorted in ascending order too.

An example of input graph format is as follows:

3

0,2:1,1.0:2,1.0

1,2:0,1.0:2,1.0

2,2:0,1.0:1,1.0

where this graph is a triangle with three vertices

    + format for features

Each line is the feature vector representation for a node formatted as follows:

[node_id] tab [number_of_features d] tab [feature_index1]tab[feature_weight1]...[feature_indexd] tab [feature_weightd]

An example of input feature format is as follows:

0	2	0	5.000000e-01	1	5.000000e-01

1	2	0	5.000000e-01	1	5.000000e-01

2	2	0	5.000000e-01	1	5.000000e-01

3	2	2	5.000000e-01	3	5.000000e-01

4	2	2	5.000000e-01	3	5.000000e-01

5	2	2	5.000000e-01	3	5.000000e-01

6	1	3	1

7	1	3	1

8	1	3	1

9	1	3	1

10	1	3	1

#Output format

Each line i gives the partitioning id for the i-th vertex


