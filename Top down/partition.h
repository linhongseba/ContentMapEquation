#include"graph.h"
#include <time.h>

void Initpart(Entry *&G, int *&partition, int N, int &numpart, int k){
	int i=0;
	for(i=0;i<N;i++){
		partition[i]=ccid[i];
	}
	numpart=ccnum;
	int a;
	int b;
	double prob;
	int newid;
	for(i=0;i<N;i++){
		if(partition[i]==maxccid){
			a=rand();
			b=rand();
			if(a+b>0)
				prob=(double)a/(a+b);
			else
				prob=0;
			if(prob>(1.0/k)){
			    newid=a%k;
				partition[i]=numpart+newid;
				if(partition[i]>=N)
					partition[i]=N-1;
			}
		}
	}
	numpart+=k;
	if(numpart>=N)
		numpart=N;
}

void reorderpartition(int *&partition, int *newpartid, int N,int &numpart){
	int i=0;
	int unum=0;
	for(i=0;i<N;i++){
		if(Ni[i]>0){
			newpartid[i]=unum;
			Ni[unum]=Ni[i];
			pam[unum]=pam[i];
			qi[unum]=qi[i];
			HPi[unum]=HPi[i];
			unum++;
		}
	}
	for(i=unum;i<N;i++){
		Ni[i]=0;
		pam[i]=0;
		qi[i]=0;
		HPi[i]=0;
	}
	for(i=0;i<N;i++){
		partition[i]=newpartid[partition[i]];
	}
	numpart=unum;
}

double getentropy(double prob){
	double res=0;
	if(prob<=0)
		return 0;
	res=-prob*log2(prob);
	return res;
}

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

/*
Init pam and Ni
*/
void InitNi(int *partition, int N){
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
double *temptf;// temptf length is equal to number of words/features +1
int tfsize; //number of words/features

//
double singlemodualcontent(Entry*dict, int v, int partid){
	double res=0;
	int j=0;
	for(j=0;j<dict[v].key;j++){
		res+=getentropy(dict[v].weight[j]);
	}
	res*=pam[partid];
	return res;
}

double singlemodualextend(Entry *dict, int v, int partid){
	double res=0;
	int j=0;
	for(j=0;j<dict[v].key;j++){
		res+=getentropy(dict[v].weight[j]);
	}
	res*=(pam[partid]+qi[partid]);
	return res;
}

double singlemodualdict(Entry*dict, int v, int partid){
	double res=0;
	int j=0;
	for(j=0;j<dict[v].key;j++){
		res+=getentropy(dict[v].weight[j]*pam[partid]);
	}
	return res;
}

double singlemodualcomb(Entry *&dict, int v, int partid){
	double res=0;
	int j=0;
	for(j=0;j<dict[v].key;j++){
		res+=getentropy(dict[v].weight[j]*pam[partid]);
		res+=pam[partid]*getentropy(dict[v].weight[j]);
	}
	return res;
}

double modualcontent(Entry* &dict, int *partition, int N, int m){
	//compute the contentcodelength for modual m (partition m)
	int i=0;
	int j=0;
	int word;
	double weight;
	double sum;
	double res;
	res=0;
	for(i=0;i<tfsize;i++)
		temptf[i]=0;
	for(i=0;i<N;i++){
		if(partition[i]!=m)
			continue;
		for(j=0;j<dict[i].key;j++){
			word=dict[i].nbv[j];
			weight=dict[i].weight[j];
			temptf[word]+=weight*pa[i];
		}
	}
	sum=pam[m];
	for(i=0;i<tfsize;i++){
		if(sum>0)
			res+=getentropy(temptf[i]/sum);
	}
	res*=pam[m];
	return res;
}

/*
compute the sum of content codelength for all the module
*/
double contentlength(Entry* &dict, int *partition, int N, int numpart){
	int m;
	double res;
	res=0;
	for(m=0;m<numpart;m++){
		res+=modualcontent(dict,partition,N,m);
	}
	return res;
}

/*
//compute the dictionary codelength for modual m (partition m)
*/
double modualdict(Entry* &dict, int *partition, int N, int m){
	int i=0;
	int j=0;
	int word;
	double weight;
	double res;
	res=0;
	for(i=0;i<tfsize;i++)
		temptf[i]=0;
	for(i=0;i<N;i++){
		if(partition[i]!=m)
			continue;
		for(j=0;j<dict[i].key;j++){
			word=dict[i].nbv[j];
			weight=dict[i].weight[j];
			temptf[word]+=weight*pa[i];
		}
	}
	for(i=0;i<tfsize;i++){
		res+=getentropy(temptf[i]);
	}
	return res;
}

double dictlength(Entry* &dict, int *partition, int N, int numpart){
	int m;
	double res;
	res=0;
	for(m=0;m<numpart;m++){
		res+=modualdict(dict,partition,N,m);
	}
	return res;
}

double modualextend(Entry* &dict, int *partition,int N, int m){
	int i=0;
	int j=0;
	int word;
	double weight;
	double sum;
	double res;
	res=0;
	for(i=0;i<tfsize;i++)
		temptf[i]=0;
	for(i=0;i<N;i++){
		if(partition[i]!=m)
			continue;
		for(j=0;j<dict[i].key;j++){
			word=dict[i].nbv[j];
			weight=dict[i].weight[j];
			if(word>=0&&word<tfsize)
				temptf[word]+=weight*pa[i];
		}
	}
	sum=pam[m]+qi[m];
	for(i=0;i<tfsize;i++){
		if(sum>0)
			res+=getentropy(temptf[i]/sum);
	}
	res*=sum;
	return res;
}

double extendcontentlength(Entry* &dict, int *partition, int N, int numpart){
	int m;
	double res;
	res=0;
	for(m=0;m<numpart;m++){
		res+=modualextend(dict,partition,N,m);
	}
	return res;
}

double modualcombo(Entry* &dict, int *partition, int N,int m){
	int i=0;
	int j=0;
	int word;
	double weight;
	double sum;
	double res;
	for(i=0;i<tfsize;i++)
		temptf[i]=0;
	for(i=0;i<N;i++){
		if(partition[i]!=m)
			continue;
		for(j=0;j<dict[i].key;j++){
			word=dict[i].nbv[j];
			weight=dict[i].weight[j];
			temptf[word]+=weight*pa[i];
		}
	}
	res=0;
	sum=pam[m];
	for(i=0;i<tfsize;i++){
		if(sum>0)
			res+=getentropy(temptf[i]/sum);
	}
	res=pam[m]*res;
	for(i=0;i<tfsize;i++)
		res+=getentropy(temptf[i]);
	return res;
}

double dictcontentlength(Entry* &dict, int *partition, int N, int numpart){
	int m;
	double res;
	res=0;
	for(m=0;m<numpart;m++){
		res+=modualcombo(dict,partition,N,m);
	}
	return res;
}

double deltastruct(Entry *&G, int *partition, int N, int numpart, int v, int targetid,double &newqisrc,double &newqitarget, double &newHPisrc, double &newHPitarget){
	//the change of code length from structure description
	//when a node v is moving from its original partition to targetpartition;
	//qi[i] is coresponding to 
	///tau*(n-ni)/(n-1)\sum_{u\in part[i]}pa[u]+(1-tau)\sum_{u \in part[i]}{v\not\in part[i]}pa[u]w[u][v]

	double originalleft;
	double originalright;
	double newpisrc;
	double newpitarget;
	double qsum=0;
	int i=0;
	int nbv;
	for(i=0;i<numpart;i++){
		qsum+=qi[i];
	}
	int src=partition[v];
	originalleft=0;
	for(i=0;i<numpart;i++){
		if(qsum>0)
			originalleft+=qsum*getentropy(qi[i]/qsum);
	}
	originalright=(qi[src]+pam[src])*HPi[src]+(qi[targetid]+pam[targetid])*HPi[targetid];
	newqisrc=qi[src]-TAU*pam[src]/(N-1)-TAU*pa[v]/(N-1);
	newqitarget=qi[targetid]+TAU*pam[targetid]/(N-1)+TAU*pa[v]/(N-1);
	for(i=0;i<G[v].key;i++){
		nbv=G[v].nbv[i];
		if(partition[nbv]==src){
			newqisrc+=((1-TAU)*G[v].weight[i])*pa[v];
		}
		if(partition[nbv]!=src){
			newqisrc-=((1-TAU)*G[v].weight[i])*pa[v];
		}
		if(partition[nbv]==targetid){
			newqitarget-=((1-TAU)*G[v].weight[i])*pa[v];
		}
		if(partition[nbv]!=targetid){
			newqitarget+=((1-TAU)*G[v].weight[i])*pa[v];
		}
	}
	newpisrc=newqisrc+pam[src]-pa[v];
	newpitarget=newqitarget+pam[targetid]+pa[v];
	double prob;
	if(newpisrc>0)
		prob=newqisrc/newpisrc;
	else
		prob=0;
	newHPisrc=getentropy(prob);
	if(newpitarget>0)
		prob=newqitarget/newpitarget;
	else
		prob=0;
	newHPitarget=getentropy(prob);
	for(i=0;i<N;i++){
		if(i!=v&&partition[i]==src){
			if(newpisrc>0)
				prob=pa[i]/newpisrc;
			else
				prob=0;
			newHPisrc+=getentropy(prob);
		}
		if(i==v||partition[i]==targetid){
			if(newpitarget>0)
				prob=pa[i]/newpitarget;
			else
				prob=0;
			newHPitarget+=getentropy(prob);
		}
	}
	double newleft=0;
	double newqsum=qsum-qi[src]-qi[targetid]+newqisrc+newqitarget;
	for(i=0;i<numpart;i++){
		if(newqsum>0){
			if(i==src)
				newleft+=newqsum*getentropy(newqisrc/newqsum);
			else
				if(i==targetid){
					newleft+=newqsum*getentropy(newqitarget/newqsum);
				}else
					newleft+=newqsum*getentropy(qi[i]/newqsum);
		}

	}
	if(targetid==numpart)
		newleft+=newqsum*getentropy(newqitarget/newqsum);
	double newright=newpisrc*newHPisrc+newpitarget*newHPitarget;
	double delta=newleft+newright-originalleft-originalright;
	return delta;
}

double deltadist(Entry*&dict,  int *partition, int N, int numpart, int v, int target){
	int src;
	double originalsrc;
	double originaltar;
	double newsrc;
	double newtar;
	src=partition[v];
	originalsrc=modualdict(dict,partition,N,src);
	originaltar=modualdict(dict,partition,N,target);
	//change partition, change pam
	partition[v]=target;
	pam[src]-=pa[v];
	pam[target]+=pa[v];
	newsrc=modualdict(dict,partition,N,src);
	newtar=modualdict(dict,partition,N,target);
	double delta=newsrc+newtar-originalsrc-originaltar;
	//set back changes
	partition[v]=src;
	pam[src]+=pa[v];
	pam[target]-=pa[v];
	return delta;
}

double deltacontent(Entry*&dict,  int *partition, int N, int numpart, int v, int target){
	int src;
	double originalsrc;
	double originaltar;
	double newsrc;
	double newtar;
	src=partition[v];
	originalsrc=modualcontent(dict,partition,N,src);
	originaltar=modualcontent(dict,partition,N,target);
	//change partition, change pam
	partition[v]=target;
	pam[src]-=pa[v];
	pam[target]+=pa[v];
	newsrc=modualcontent(dict,partition,N,src);
	newtar=modualcontent(dict,partition,N,target);
	double delta=newsrc+newtar-originalsrc-originaltar;
	//set back changes
	partition[v]=src;
	pam[src]+=pa[v];
	pam[target]-=pa[v];
	return delta;
}

double deltaextend(Entry *&dict, int *partition, int N, int numpart, int v, int target, double newqisrc, double newqitar){
	int src;
	double originalsrc;
	double originaltar;
	double newsrc;
	double newtar;
	src=partition[v];
	double oldqisrc=qi[src];
	double oldqitar=qi[target];
	originalsrc=modualextend(dict,partition,N,src);
	originaltar=modualextend(dict,partition,N,target);
	//change partition, change pam
	partition[v]=target;
	pam[src]-=pa[v];
	pam[target]+=pa[v];
	//change qi
	qi[src]=newqisrc;
	qi[target]=newqitar;
	newsrc=modualextend(dict,partition,N,src);
	newtar=modualextend(dict,partition,N,target);
	double delta=newsrc+newtar-originalsrc-originaltar;
	//set back changes
	partition[v]=src;
	pam[src]+=pa[v];
	pam[target]-=pa[v];
	qi[src]=oldqisrc;
	qi[target]=oldqitar;
	return delta;
}

double deltacontentdist(Entry*&dict,  int *partition, int N, int numpart, int v, int target){
	int src;
	double originalsrc;
	double originaltar;
	double newsrc;
	double newtar;
	src=partition[v];
	originalsrc=modualcombo(dict,partition,N,src);
	originaltar=modualcombo(dict,partition,N,target);
	//change partition, change pam
	partition[v]=target;
	pam[src]-=pa[v];
	pam[target]+=pa[v];
	newsrc=modualcombo(dict,partition,N,src);
	newtar=modualcombo(dict,partition,N,target);
	double delta=newsrc+newtar-originalsrc-originaltar;
	//set back changes
	partition[v]=src;
	pam[src]+=pa[v];
	pam[target]-=pa[v];
	return delta;
}

double SCpartition(Entry*G, Entry*dict, int N, int *&partition, int numpart, int addcontent, int adddict, int adde, int k){
	int i=0;
	int j=0;
	int *temppart=new int[N];
	int temppartnum=0;
	int maxpart;
	maxpart=(int)(4*sqrt((double)N));
	if(maxpart>N/2)
		maxpart=N/2;
	ccnum=0;
	Connectedcomponent(G,N);
	cout<<"connected component number "<<ccnum<<endl;
	if(maxpart<=ccnum)
		maxpart=ccnum+2*k;
	if(maxpart>N)
		maxpart=N;
	qi=new double[N];
	HPi=new double[N];
	double codelength=4*N;
	double initcodelength;
	int initnum=(int)(sqrt((double)N));
	pam=new double[N];
	Ni=new int[N];
	double hp,hq;
	int *newpartid=new int[N];
	for(i=0;i<initnum;i++){
		Initpart(G,temppart,N,temppartnum,k);
		InitNi(temppart,N);
		reorderpartition(temppart,newpartid,N,temppartnum);
		InitNi(temppart,N);
		for(j=0;j<N;j++){
			qi[j]=0;
			HPi[j]=0;
		}
		initcodelength=structlength(G,temppart,N,temppartnum,hp,hq);
		if(adde)
			initcodelength=extendcontentlength(dict,temppart,N,temppartnum);
		if(addcontent&&adddict)
			initcodelength+=dictcontentlength(dict,temppart,N,temppartnum);
		else
			if(addcontent)
				initcodelength+=contentlength(dict,temppart,N,temppartnum);
			else
				if(adddict)
					initcodelength+=dictlength(dict,temppart,N,temppartnum);
		if(initcodelength<codelength){
			codelength=initcodelength;
			for(j=0;j<N;j++)
				partition[j]=temppart[j];
			numpart=temppartnum;
		}
	}
	delete []temppart;
	InitNi(partition,N);
	for(j=0;j<N;j++){
		qi[j]=0;
		HPi[j]=0;
	}
	initcodelength=structlength(G,partition,N,numpart,hp,hq);
	if(adde)
		initcodelength+=extendcontentlength(dict,partition,N,numpart);
	if(addcontent&&adddict)
		initcodelength+=dictcontentlength(dict,partition,N,numpart);
	else
		if(addcontent)
			initcodelength+=contentlength(dict,partition,N,numpart);
		else
			if(adddict)
				initcodelength+=dictlength(dict,partition,N,numpart);
    cout << "initial code length:" << initcodelength << endl;
	//double desclength=0;
	int v=0;
	double mindelta=0;
	double qisrcstore=0;
	double HPisrcstore=0;
	double qitarstore=0;
	double HPitarstore=0;
	int vstore=-1;
	int targetstore=-1;
	//int maxid;
	int src=0;
	int tar=0;
	double delta=0;
	double newqisrc=0;
	double newqitar=0;
	double newHPisrc=0;
	double newHPitar=0;
	int iternum=0;
	int swapnum=0;
	//double maxmodual;
	//double localmodual;
	//double maxpa;
	//int id;
	int neigh;
	//double localpa;
	//int preid=-1;
	char * ismarker=new char[N];
	memset(ismarker,'n',sizeof(char)*N);
	bool flag=false;
	do{
		//find the maximum module;
		/*maxid=0;
		maxmodual=0;
		for(i=0;i<numpart;i++){
			localmodual=(qi[i]+pa[i])*HPi[i];
			if(adde)
				localmodual+=modualextend(dict,partition,N,i);
			if(addcontent&&adddict)
					localmodual+=modualcombo(dict,partition,N,i);
				else
					if(addcontent)
						localmodual+=modualcontent(dict,partition,N,i);
					else
						if(adddict)
							localmodual+=modualdict(dict,partition,N,i);
			if(localmodual>maxmodual){
				maxmodual=localmodual;
				maxid=i;
			}
		}*/
		mindelta=0;
		////find the v in maxid such that moving v from maxid to targetid reduce the cost
		vstore=-1;
		targetstore=-1;
		flag=false;
		//maxpa=0;
		//for(i=0;i<N;i++){
		//	if(i!=preid&&partition[i]==maxid&&G[i].key>0){
		//		localpa=pa[i];
		//		if(adde)
		//			localpa+=singlemodualextend(dict,i,partition[i]);
		//		if(addcontent&&adddict)
		//			localpa+=singlemodualcomb(dict,i,partition[i]);
		//		else
		//			if(addcontent)
		//				localpa+=singlemodualcontent(dict,i,partition[i]);
		//			else
		//				if(adddict)
		//					localpa+=singlemodualdict(dict,i,partition[i]);
		//		for(j=0;j<G[i].key;j++){
		//			neigh=G[i].nbv[j];
		//			if(partition[neigh]!=partition[i]){
		//				localpa+=pa[i]*G[i].weight[j];
		//			}
		//		}
		//		if(localpa>maxpa){
		//			maxpa=localpa;
		//			id=i;
		//		}
		//	}
		//}
		for(i=0;i<N;i++){
			if(G[i].key<1)
				continue;
			v=i;
			src=partition[v];
			//if(src!=maxid)
				//continue;
			mindelta=0;
			////find the v in maxid such that moving v from maxid to targetid reduce the cost
			vstore=-1;
			targetstore=-1;
			memset(ismarker,'n',sizeof(char)*N);
			for(j=0;j<G[v].key;j++){
				neigh=G[v].nbv[j];
				tar=partition[neigh];
				if(tar!=src&&ismarker[tar]=='n'){
					ismarker[tar]='y';
					delta=deltastruct(G,partition,N,numpart,v,tar,newqisrc,newqitar,newHPisrc,newHPitar);
					if(adde)
						delta+=deltaextend(dict,partition,N,numpart,v,tar,newqisrc,newqitar);
					if(addcontent&&adddict)
						delta+=deltacontentdist(dict,partition,N,numpart,v,tar);
					else
						if(addcontent)
							delta+=deltacontent(dict,partition,N,numpart,v,tar);
						else
							if(adddict)
								delta+=deltadist(dict,partition,N,numpart,v,tar);
					if(delta-mindelta<0){
						mindelta=delta;
						qisrcstore=newqisrc;
						qitarstore=newqitar;
						HPisrcstore=newHPisrc;
						HPitarstore=newHPitar;
						vstore=v;
						targetstore=tar;
					}
				}
			}
			tar=numpart;
			if(numpart<N){
				delta=deltastruct(G,partition,N,numpart,v,tar,newqisrc,newqitar,newHPisrc,newHPitar);
				if(adde)
					delta+=deltaextend(dict,partition,N,numpart+1,v,tar,newqisrc,newqitar);
				if(addcontent&&adddict)
					delta+=deltacontentdist(dict,partition,N,numpart+1,v,tar);
				else
					if(addcontent)
						delta+=deltacontent(dict,partition,N,numpart+1,v,tar);
					else
						if(adddict)
							delta+=deltadist(dict,partition,N,numpart+1,v,tar);
				if(delta-mindelta<0){
					mindelta=delta;
					qisrcstore=newqisrc;
					qitarstore=newqitar;
					HPisrcstore=newHPisrc;
					HPitarstore=newHPitar;
					vstore=v;
					targetstore=tar;
				}
			}
			if(vstore!=-1&&targetstore!=-1){
				//preid=vstore;
				//move vstore to targetstore
				//update pi,qi,hpi,ni, pam
				//update numpart;
				//update partition
				src=partition[vstore];
				pam[src]-=pa[vstore];
				pam[targetstore]+=pa[vstore];
				Ni[src]-=1;
				Ni[targetstore]+=1;
				qi[src]=qisrcstore;
				qi[targetstore]=qitarstore;
				HPi[src]=HPisrcstore;
				HPi[targetstore]=HPitarstore;
				partition[vstore]=targetstore;
				if(targetstore==numpart){
					numpart++;
				}
				flag=true;
				swapnum++;
			}
		}
		//if(vstore!=-1&&targetstore!=-1){
		//	preid=vstore;
		//	//move vstore to targetstore
		//	//update pi,qi,hpi,ni, pam
		//	//update numpart;
		//	//update partition
		//	src=partition[vstore];
		//	pam[src]-=pa[vstore];
		//	pam[targetstore]+=pa[vstore];
		//	Ni[src]-=1;
		//	Ni[targetstore]+=1;
		//	qi[src]=qisrcstore;
		//	qi[targetstore]=qitarstore;
		//	HPi[src]=HPisrcstore;
		//	HPi[targetstore]=HPitarstore;
		//	partition[vstore]=targetstore;
		//	if(targetstore==numpart){
		//		numpart++;
		//	}
		//}else
		//	break;
		iternum++;
		cout<<"iternum is "<<iternum<<endl;
	}while(flag&&iternum<50);
	double finalcodelength;
	double elength=0;
	double clength=0;
	double dlength=0;
	double cdlength=0;
	//rename partitioning
	reorderpartition(partition,newpartid,N,numpart);
	InitNi(partition,N);
	cout<<"swapnum is\t"<<swapnum<<endl;
	cout<<"number of partition is\t"<<numpart<<endl;
	finalcodelength=structlength(G,partition,N,numpart,hp,hq);
	cout<<"right\tleft\tstruct\tcontent\tdict\tCD\textend\tfinal"<<endl;
	cout<<hp<<"\t"<<hq<<"\t"<<finalcodelength;
	if(adde){
		elength=extendcontentlength(dict,partition,N,numpart);
		finalcodelength+=elength;
	}
	if(addcontent&&adddict){
		cdlength=dictcontentlength(dict,partition,N,numpart);
		finalcodelength+=cdlength;
	}
	else
		if(addcontent){
			clength=contentlength(dict,partition,N,numpart);
			finalcodelength+=clength;
		}
	else
		if(adddict){
			dlength=dictlength(dict,partition,N,numpart);
			finalcodelength+=dlength;
		}
		cout<<"\t"<<clength;
		cout<<"\t"<<dlength;
		cout<<"\t"<<cdlength;
		cout<<"\t"<<elength;
		cout<<"\t"<<finalcodelength<<endl;
	//double compression=(initcodelength-finalcodelength)/initcodelength*100;
	//cout<<"compression ratio is\t"<<compression<<"%"<<endl;
	delete []newpartid;
    return finalcodelength;
}