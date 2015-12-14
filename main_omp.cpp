 #include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <fstream>
#include <cctype>
#include <ctime>
#include <vector>
#include <set>
#include <sstream>
#include "pstream.h"
#include <sys/time.h>
using namespace std;
using redi::pstream;
#define PZERO 0.00001f
#define MZERO -0.00001f
int main(int argc, char** argv)
{
	FILE * ifp = fopen(argv[1], "r");
	struct timeval startTime, endTime;
    long totalTime;
    gettimeofday(&startTime, NULL);
	
	
	int n_id, n_sz;
	int n = 0;
	char ch[100];
	fscanf(ifp, "%d", &n);
	int * kn = (int *)malloc(n* sizeof(int));
	int * dimension = (int *)malloc((n+1)* sizeof(int));
	dimension[0]=1;
	int sumall=0;
	int maxdim = 0;
	for(int i=0;i<n;++i){
		fscanf(ifp, "%d", kn+i);
		sumall+=kn[i];
		dimension[i+1]=dimension[i]*kn[i];
		maxdim = (maxdim>kn[i])?maxdim:kn[i];
		
	}
	

	int * util = (int *)malloc(dimension[n]* n * sizeof(int*));
	int * currentIndex = (int *)malloc(n* sizeof(int));
	for(int i = 0;i<dimension[n];++i){
		
		fscanf(ifp, "%[^\[]",&ch);
		fscanf(ifp, "%c\n",&ch);
		fscanf(ifp, "%d", currentIndex);
		for(int j=1; j<n-1; j++){
			fscanf(ifp, "%d", currentIndex+j);
		}
		fscanf(ifp, "%d%c", currentIndex+n-1,ch);
		//sub index to index
		int ind = 0;
		for(int j=0;j<n;++j){
			ind+=(currentIndex[j]-1)*dimension[j];
		}
		fscanf(ifp, "%[^\[]",&ch);
		fscanf(ifp, "%c\n",&ch);
		for(int j=0; j<n; j++){
			fscanf(ifp, "%d", util+ind*n+j);
		}

		
		fscanf(ifp, "]");
	}
	fclose(ifp);

	gettimeofday(&endTime, NULL);
	totalTime =  (endTime.tv_sec - startTime.tv_sec) * 1000000L;
    totalTime += (endTime.tv_usec - startTime.tv_usec);
	printf("load time: %f sec\n",double(totalTime)/1000000);
	gettimeofday(&startTime, NULL);
	//startTime=clock();

	char num[1000];
	//FILE * ofp = fopen("outp.txt", "w");
	//fprintf(ofp, "%d\n",num_eq);
	int combnum=1,shiftcnt=0,shiftmask;
	int* combmask=(int*)malloc(n*sizeof(int));
	int* shiftall=(int*)malloc(n*sizeof(int));
	for(int comb=0;comb<n;++comb)
	{
		shiftmask = (1<<kn[comb])-1;
		combmask[comb]=shiftmask;
		shiftall[comb]=shiftcnt;
		shiftcnt+=kn[comb];
		
	}
	combnum<<=shiftcnt;
	#pragma omp parallel num_threads(1) in(combnum) in(combmask,shiftall,kn,dimension,currentIndex:length(n)) in(util:length(n*dimension[n]))
  {
  //parallelize via OpenMP on MIC
	int count;
	#pragma omp for
	for(count=0;count<combnum;++count){
		vector<string> eqs;

		int num_eq = 0;
		string eq0,eqj;
		set<pair<int,int>> idset;
		set<pair<int,int>>::iterator iter;

		
		//num_eq = 0;
		int valid = 1;
		for(int cnt=0;cnt<n;++cnt){
			int curmask = combmask[cnt]&(count>>shiftall[cnt]);
			if(curmask==0){
				valid = 0;
				break;
			}
		} 
		if(!valid)continue;
		idset.clear();
		eqs.clear();
		for(int i=0;i<n;++i){
			int curnum = (count>>shiftall[i]) & combmask[i];
			for(int j=0;j<kn[i];++j){
				if((curnum&1)){
					idset.insert(pair<int,int>(i,j));
					num_eq++;
				}
				curnum>>=1;
			}
		}
		for(int i=0;i<n;++i){
			

			int done = 0;
			int cntcure = 0;
			for(int j=0;j<kn[i];++j){
				iter = idset.find(pair<int,int>(i,j));
				if(iter==idset.end()){
					continue;
				}
				eqj.clear();
				for(int k=0;k<dimension[n]/kn[i];++k) {
					int smallk =k%dimension[i];
					int curk = (k-smallk)*kn[i] + smallk + dimension[i]*j;
					int idc = 0;
					int validcur=1;
					string tmpstr;
					if(util[curk*n+i]>=0){//
						tmpstr+=" + ";
					}
					
					tmpstr+=to_string(util[curk*n+i]);
					for(int w=n-1;w>=0;--w){
						currentIndex[w] = curk/dimension[w];
						curk = curk%dimension[w];
						if(w==i)continue;
						int curnum = (count>>shiftall[i]) & combmask[i];
						iter = idset.find(pair<int,int>(w,currentIndex[w]));
						if(iter==idset.end()){
							validcur =0;
							break;
						}
									
						tmpstr+=" * ";
						tmpstr+="x";
						
						tmpstr+=to_string(w+1);
						
						tmpstr+=to_string(currentIndex[w]+1);;
						
					}
					if(validcur){
						eqj = tmpstr+eqj;
					}
				}	
				if(cntcure==0){
					eq0=eqj;
				}
				else{
					eqj += " - ( "+ eq0 + " );";
					eqs.push_back(string(eqj.c_str()));
				}
				cntcure++;
			}
		}
		for(int i=0;i<n;++i){
			string eo;
			for(int j=0;j<kn[i];++j){
				iter = idset.find(pair<int,int>(i,j));
				if(iter==idset.end()){
					continue;
				}
				eo+="x";
				eo+=to_string(i+1);
				eo+=to_string(j+1);
				eo+=" + ";
			}
			eo+=" 0 - 1;";
			eqs.push_back(string(eo.c_str()));
		}
		/*printf("num eqs %d\n", num_eq);
		for(int i=0;i<eqs.size();++i){
			printf("%s\n", eqs[i].c_str());
		}*/
		
		pstream phc; // for running phc
		phc.open("./phc -b");
		if(!phc.is_open())
		{
			printf("Could not open the phc pack executable.\n");
		}
		//printf("n=%d\n",eqs.size());
			// start inputting
		phc << "n" << endl; // is system on input file? no
		phc << num_eq << endl; // number of polynomials
		phc << num_eq << endl; // number of unknowns
		for(int i=0;i<eqs.size();++i){
			phc << eqs[i].c_str()<<endl;
		}
		phc << "n" << endl;
		string opstr;
		vector<string> varname;
		double* solution = (double*)malloc(num_eq*sizeof(double));
		double** res = (double**)malloc(n*sizeof(double*)); 
		for(int place=0;place<n;++place){
			res[place] = (double*)malloc(kn[place]*sizeof(double)); 
		}
		while(phc >> opstr)
		{
			#pragma omp critical 
			{
			if(opstr.compare("solution")==0){
				phc >> opstr;
				if(opstr.compare("for")==0){
					varname.clear();
					phc >> opstr;//t
					phc >> opstr;//:
					//phc >> opstr;
					bool done=1;
					for(int i=0;i<num_eq;++i){
						string vname;
						double reald,compd;
						phc >> vname; //variable name
						phc >> opstr; //: skip
						phc >> reald;
						phc >> compd;
						if(reald<MZERO){
							done=0;
							//printf("%f\n",reald);
							break;
						}
						if(compd<MZERO||compd>PZERO){
						
							done=0;
							//printf("%f\n",compd);
							break;
						}
						solution[i] = reald;
						//printf("%d %d\n",vname.c_str()[1]-'0',vname.c_str()[2]-'0' );
						res[vname.c_str()[1]-'1'][vname.c_str()[2]-'1']=reald;
						varname.push_back(vname);
						
					}
					if(done){
						set<pair<int,int>>::iterator loopiter;
						int valid=0;
						for(int vs = 0;vs <num_eq;++vs){	
							if(solution[vs]>PZERO){
								valid=1;
								break;
							}
						}
						if(!valid){
							continue;
						}
						/*for(loopiter=idset.begin();loopiter!=idset.end();++loopiter){
							int ci = (*loopiter).first;
							int cj = (*loopiter).second;
							printf("ci = %d; cj = %d\n",ci,cj);
						}*/
						for(loopiter=idset.begin();loopiter!=idset.end();++loopiter){
							int ci = (*loopiter).first;
							int cj = (*loopiter).second;
							//printf("ci = %d; cj = %d\n",ci,cj);
							
							double probme=0;
							for(int k=0;k<dimension[n]/kn[ci];++k) {
								double curprob=1;
								int smallk =k%dimension[ci];
								int curk = (k-smallk)*kn[ci] + smallk + dimension[ci]*cj;
								int cpcurk = curk;
								for(int w=n-1;w>=0;--w){
									currentIndex[w] = curk/dimension[w];
									curk = curk%dimension[w];
									if(w==ci)continue;
									iter = idset.find(pair<int,int>(w,currentIndex[w]));
									if(iter==idset.end()){
										curprob=0;
										break;
									}
									curprob*=res[w][currentIndex[w]];
								}
								probme += curprob*util[cpcurk*n+ci];
							}
							
							
							for(int oj =0 ;oj< kn[ci];++oj){
								double probother =0;
								iter = idset.find(pair<int,int>(ci,oj));
								if(oj==cj)continue;
								if(iter!=idset.end()){//if in support
									continue;
								}
								for(int k=0;k<dimension[n]/kn[ci];++k) {
									double curprob=1;
									int smallk =k%dimension[ci];
									int curk = (k-smallk)*kn[ci] + smallk + dimension[ci]*oj;
									int cpcurk = curk;
									for(int w=n-1;w>=0;--w){
										currentIndex[w] = curk/dimension[w];
										curk = curk%dimension[w];
										if(w==ci)continue;
										iter = idset.find(pair<int,int>(w,currentIndex[w]));
										if(iter==idset.end()){
											curprob=0;
											break;
										}
										curprob*=res[w][currentIndex[w]];
									}
									probother += curprob*util[cpcurk*n+ci];
								}
								
								if(probother-probme>PZERO){
									//printf("check prob %d %d with %d: %f, %f\n", ci,cj, oj, probme,probother);
									//fflush(stdout);
									valid = 0;
								}
							}
							if(!valid){
								break;
							}
						}
						if(valid){
							
								printf("solution:\n");
								for(int i=0;i<num_eq;++i){
									printf("%s: %f\n",varname[i].c_str(),solution[i]);
								}
							
						}
					}
				}
			}
			}
			
		}
		#pragma omp critical 
		{
			phc.clear();
			phc.close();
		}
		
	}
	}
	gettimeofday(&endTime, NULL);
	totalTime =  (endTime.tv_sec - startTime.tv_sec) * 1000000L;
    totalTime += (endTime.tv_usec - startTime.tv_usec);
	printf("running time: %f sec\n",double(totalTime)/1000000);
	//fclose(ofp);
}
	