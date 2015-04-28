/* --------------------------------------------
 ModuleDictionary.h
 Created by Laura Smith on 6/11/13.
 Updated by Linhong on 07/02/13.
 -------------------------------------------- */
#ifndef _ModuleDictionary_h
#define _ModuleDictionary_h
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "graph.h"
using namespace std;
int tempsize=800;
int*tempWords;
double *tempFrequencies;
double ContentEntropy(int length, double pvalue){
	double entropy = 0;
	int i=0;
    for (i=0; i<length; i++){
        entropy += tempFrequencies[i]*log(tempFrequencies[i]);
    }
    entropy = -entropy*pvalue/log(2.0);
    
    return entropy;
}
double DictionaryEntropy(int length, double pvalue){
	 double entropy = 0;
	 int i=0;
    for (i=0; i<length; i++){
        if(tempFrequencies[i]*pvalue>0)
            entropy += tempFrequencies[i]*pvalue*log(tempFrequencies[i]*pvalue);
    }
    entropy = -entropy/log(2.0);
    return entropy;
}


/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
// ModuleDictionary Class

class ModuleDictionary{
public:
    ModuleDictionary();
    
    // Retrieve private entities
    double get_p_module();
    int get_numWords();
	int get_size();
    int* get_words();
    double* get_wordFrequency();
    void normalize();
    // Set private entities
    void set_p_module(double p_moduleIn);
    void set_numWords(int numWordsIn);
    void set_words(int* wordsIn);
    void set_wordFrequency(double *wordFrequencyIn);
	void set_size(int size);
	void copy_all(int* wordsIn, double *wordFrequencyIn, int length);
    
    // Input and output dictionaries
    void inputDictionary(double p_moduleIn, int numberOfWords, int* wordsIn, double* frequencyIn);
    void printDictionary();
	void InitDictionary(int number);
	void addwords(int wordsin, double frequencyin, int index);
    
    // Determine the individual and all dictionary contributions
    double calcContentEntropy();
    double calcDictionaryEntropy();
	
	void clear();
    
private:
    double p_module;
    int numWords;
    int* words;
    double* wordFrequency;
	int memsize;
};

// Default Constructor
ModuleDictionary::ModuleDictionary(){
    p_module = 0;
    numWords = 0;
	memsize=0;
	this->words=NULL;
	this->wordFrequency=NULL;
    return;
}
// Retrieve private entities
double ModuleDictionary::get_p_module(){ return p_module;}
int ModuleDictionary::get_numWords(){ return numWords;}
int ModuleDictionary::get_size(){return memsize;}
int* ModuleDictionary::get_words(){ return words;}
double* ModuleDictionary::get_wordFrequency(){ return wordFrequency;}


// Set the p_module = probability of being at that module
void ModuleDictionary::set_p_module(double p_moduleIn){
    p_module = p_moduleIn;
    return;
}

// Set the number of words
void ModuleDictionary::set_numWords(int numWordsIn){
    numWords = numWordsIn;
    return;
}
void ModuleDictionary::set_size(int size){
	memsize=size;
}
void ModuleDictionary::normalize(){
	double sum = 0;
	for (int i = 0; i < numWords; i++){
		sum += wordFrequency[i];
	}
	for (int i = 0; i < numWords; i++){
		if (sum > 0)
			wordFrequency[i] /= sum;
	}
};
void ModuleDictionary::InitDictionary(int number){
	memsize=number;
	if(memsize>=100000){
		words=(int*)malloc(sizeof(int)*memsize);
		wordFrequency=(double*)malloc(sizeof(double)*memsize);
	}else{
		allocatepermemory(sizeof(int)*memsize);
		words=(int*)curMemPos2;
		curMemPos2+=(sizeof(int)*memsize);
		allocatepermemory(sizeof(double)*memsize);
		wordFrequency=(double*)curMemPos2;
		curMemPos2+=(sizeof(double)*memsize);
	}
	if(words==NULL||wordFrequency==NULL){
		cout<<"could not allocate memory"<<endl;
		exit(0);
	}
	numWords =number;
}
void ModuleDictionary::addwords(int wordsin, double frequencyin, int index){
	words[index]=wordsin;
	wordFrequency[index]=frequencyin;
}
// Set the words Vector
void ModuleDictionary::set_words(int* wordsIn){
    words = wordsIn;
    return;
}
inline void ModuleDictionary::copy_all(int* wordsIn, double *wordFrequencyIn, int length){
	//cout<<"memsize is "<<memsize<<" length is "<<length<<endl;
	if(length>=memsize){
		if(length>=400)
			memsize=(length+100);
		else
			memsize=400;
		if(memsize>=100000){
			words=(int*)malloc(sizeof(int)*memsize);
			wordFrequency=(double*)malloc(sizeof(double)*memsize);
		}else{
			allocatepermemory(sizeof(int)*memsize);
			words=(int*)curMemPos2;
			curMemPos2+=(sizeof(int)*memsize);
			allocatepermemory(sizeof(double)*memsize);
			wordFrequency=(double*)curMemPos2;
			curMemPos2+=(sizeof(double)*memsize);
		}
		if(words==NULL||wordFrequency==NULL){
			cout<<"could not allocate memory"<<endl;
			exit(0);
		}
	}
	int i=0;
	for(i=0;i<length;i++){
		words[i]=wordsIn[i];
		wordFrequency[i]=wordFrequencyIn[i];
	}
}

// Set the wordFrequency Vector
void ModuleDictionary::set_wordFrequency(double* wordFrequencyIn){
    wordFrequency = wordFrequencyIn;
    return;
}


// This adds in a dictionary to a module
void ModuleDictionary::inputDictionary(double p_moduleIn, int numberOfWords, int* wordsIn, double* frequencyIn){
    p_module =  p_moduleIn;
    numWords = numberOfWords;
    words = wordsIn;
    wordFrequency = frequencyIn;
    return;
}
void ModuleDictionary::clear(){
	if(words!=NULL){
		if(memsize>=100000){
			free(words);
		}
		words=NULL;
	}
	if(wordFrequency!=NULL){
		if(memsize>=100000){
			free(wordFrequency);
		}
		wordFrequency=NULL;
	}
	memsize=0;
}
// This calculates the single module dictionary entropy
inline double ModuleDictionary::calcContentEntropy(){
    double entropy = 0;
	int i=0;
    for (i=0; i<numWords; i++){
        entropy += wordFrequency[i]*log(wordFrequency[i]);
    }
    entropy = -entropy*p_module/log(2.0);
    
    return entropy;
}

// This calculates the total dictionary entropy
inline double ModuleDictionary::calcDictionaryEntropy(){
    double entropy = 0;
	int i=0;
    for (i=0; i<numWords; i++){
        if(wordFrequency[i]*p_module>0)
            entropy += wordFrequency[i]*p_module*log(wordFrequency[i]*p_module);
    }
    entropy = -entropy/log(2.0);
    return entropy;
}


inline void merge(ModuleDictionary& mod1, ModuleDictionary& mod2, int &k, int newsize){
	int i=0;
	int j=0;
	double combinedFrequency;
	if(newsize>=tempsize){
		tempWords=(int*)realloc(tempWords,sizeof(int)*newsize);
		tempFrequencies=(double*)realloc(tempFrequencies,sizeof(double)*newsize);
		tempsize=newsize;
	}
	while(i<mod1.get_numWords()&&j<mod2.get_numWords()){
		if(mod1.get_words()[i] == mod2.get_words()[j]){
			tempWords[k]=mod1.get_words()[i];
			combinedFrequency = mod1.get_wordFrequency()[i]*mod1.get_p_module() + mod2.get_wordFrequency()[j]*mod2.get_p_module();
			tempFrequencies[k]=combinedFrequency/(mod1.get_p_module()+mod2.get_p_module());
			i++;
			j++;
			k++;
		}else
			if(mod1.get_words()[i] < mod2.get_words()[j]){
				tempWords[k]=mod1.get_words()[i];
				combinedFrequency = mod1.get_wordFrequency()[i]*mod1.get_p_module();
				tempFrequencies[k]=combinedFrequency/(mod1.get_p_module()+mod2.get_p_module());
				i++;
				k++;
			}
			else{
				tempWords[k]=mod2.get_words()[j];
				combinedFrequency = mod2.get_wordFrequency()[j]*mod2.get_p_module();
				tempFrequencies[k]=combinedFrequency/(mod1.get_p_module()+mod2.get_p_module());
				j++;
				k++;
			}
	}
	while(i<mod1.get_numWords()){
		tempWords[k]=mod1.get_words()[i];
		combinedFrequency = mod1.get_wordFrequency()[i]*mod1.get_p_module();
		tempFrequencies[k]=combinedFrequency/(mod1.get_p_module()+mod2.get_p_module());
		i++;
		k++;
	}
	while(j<mod2.get_numWords()){
		tempWords[k]=mod2.get_words()[j];
		combinedFrequency = mod2.get_wordFrequency()[j]*mod2.get_p_module();
		tempFrequencies[k]=combinedFrequency/(mod1.get_p_module()+mod2.get_p_module());
		j++;
		k++;
	}
}
// Replaces mod1 with the merged dictionary of mod1 and mod2
void mergeDictionaries(ModuleDictionary& mod1, ModuleDictionary& mod2){
	int newsize=mod1.get_numWords()+mod2.get_numWords();
    if(newsize>0){
		int k=0;
		merge(mod1,mod2,k,newsize);
		//cout<<"OK in merging dictionaries"<<endl;
		// Replace module one with the merged modules
		//cout<<mod1.get_numWords()<<" "<<mod2.get_numWords()<<endl;
		mod1.set_p_module(mod1.get_p_module() + mod2.get_p_module());
		//cout<<"K is "<<k<<"tempsize is "<<tempsize<<endl;
		mod1.copy_all(tempWords,tempFrequencies,k);
		//cout<<"OK in update mod1"<<endl;
		//mod1.set_words(tempWords);
		//mod1.set_wordFrequency(tempFrequencies);
		mod1.set_numWords(k);
	}else{
		mod1.set_p_module(mod1.get_p_module() + mod2.get_p_module());
		mod1.set_numWords(0);
		mod1.set_words(NULL);
		mod1.set_wordFrequency(NULL);
		mod1.set_size(0);
	}
	//mod1.printDictionary();
}

// Print the module dictionary to the screen
void ModuleDictionary::printDictionary(){
    std::cout << "\n" << p_module << "\n" ;
	int i=0;
    for(i=0; i<numWords; i++){
        std::cout <<i <<"  " << words[i] << "  " << wordFrequency[i] << "\n";
    }
    return;
}



// Replaces mod1 with the merged dictionary of mod1 and mod2
inline void replaceDictionaries(ModuleDictionary& mod1, ModuleDictionary& mod2){
    mod1.set_p_module(mod2.get_p_module());
	mod1.set_words(mod2.get_words());
	mod1.set_wordFrequency(mod2.get_wordFrequency());
	mod1.set_size(mod2.get_size());
    mod1.set_numWords(mod2.get_numWords());
    return;
}

// Calculates the difference in dictionary entropy of the merged dictionary of mod1 and mod2
inline double deltaDictionaryEntropy(ModuleDictionary &mod1, ModuleDictionary &mod2){
    // Individual Modules Content
    double priorEntropy = mod1.calcContentEntropy() + mod2.calcContentEntropy();
    
    // Total Dictionary Content
    priorEntropy += mod1.calcDictionaryEntropy() + mod2.calcDictionaryEntropy();

    //compute the merged module
    //mergeDictionaries(mod1, mod2);
	int newsize=mod1.get_numWords()+mod2.get_numWords();
	double pvalue=mod1.get_p_module()+mod2.get_p_module();
	int k=0;
	merge(mod1,mod2,k,newsize);
    // Individual Modules Content
    //double newEntropy = mod1.calcContentEntropy();
	double newEntropy=ContentEntropy(k,pvalue);
    
    // Total Dictionary Content
    newEntropy +=DictionaryEntropy(k,pvalue);

    return newEntropy-priorEntropy;
}


// Calculates the difference in content entropy of the merged dictionary of mod1 and mod2
inline double deltaContentOnlyEntropy(ModuleDictionary &mod1, ModuleDictionary &mod2){
    // Individual Modules Content
    double priorEntropy = mod1.calcContentEntropy() + mod2.calcContentEntropy();
    
    // Replace mod1 with the merged module
    //mergeDictionaries(mod1, mod2);
	//compute the merged module
	int newsize=mod1.get_numWords()+mod2.get_numWords();
	double pvalue=mod1.get_p_module()+mod2.get_p_module();
	int k=0;
	merge(mod1,mod2,k,newsize);
    
    // Individual Modules Content
    double newEntropy = ContentEntropy(k,pvalue);
    return newEntropy-priorEntropy;
}


// Calculates the difference in dictionary only entropy of the merged dictionary of mod1 and mod2
inline double deltaDictionaryOnlyEntropy(ModuleDictionary &mod1, ModuleDictionary &mod2){
 
    // Total Dictionary Content
    double priorEntropy = mod1.calcDictionaryEntropy() + mod2.calcDictionaryEntropy();
    
    // Replace mod1 with the merged module
    //mergeDictionaries(mod1, mod2);
    int newsize=mod1.get_numWords()+mod2.get_numWords();
	double pvalue=mod1.get_p_module()+mod2.get_p_module();
	int k=0;
	merge(mod1,mod2,k,newsize);
    // Total Dictionary Content
    double newEntropy =DictionaryEntropy(k,pvalue);
    return newEntropy-priorEntropy;
}
#endif
