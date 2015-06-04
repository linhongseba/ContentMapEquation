# ContentMapEquation
The implementation for the paper Laura M Smith, Linhong Zhu, Kristina Lerman, and Allon G Percus. Partitioning Networks with Node Attributes by Compressing Information Flow. In arXiv preprint arXiv:1405.4332.

#Installation
1.Requirement:
g++, gcc compiler

2.platform support:
Mac, Windows, Linux

3.Compiling
go into either the Top down or Bottom up, type make to automatically compile and generate the executable files

#Usage
1. Bottom up approach

2. Top down approach
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
1. format for graph
2. format for features

#Output format


