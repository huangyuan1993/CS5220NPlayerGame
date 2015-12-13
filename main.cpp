#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <fstream>
#include <cctype>
#include <ctime>
#include <vector>
#include <sstream>
#include "pstream.h"
#include <sys/time.h>
using namespace std;
using redi::pstream;
#define PZERO 0.00000000001f
#define MZERO -0.00000000001f
int main(int argc, char** argv)
{
	FILE * ifp = fopen(argv[1], "r");
	struct timeval startTime, endTime;
    long totalTime;
    gettimeofday(&startTime, NULL);
	int n_id, n_sz; 
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &n_id); 
	MPI_Comm_size(MPI_COMM_WORLD, &n_sz);
	int n = 0;
	char ch[100];
	fscanf(ifp, "%d", &n);
	int * kn = (int *)malloc(n* sizeof(int));
	int * dimension = (int *)malloc((n+1)* sizeof(int));
	dimension[0]=1;
	int sumall=0;

	for(int i=0;i<n;++i){
		fscanf(ifp, "%d", kn+i);
		sumall+=kn[i];
		dimension[i+1]=dimension[i]*kn[i];

	}
	
	int ** util = (int **)malloc(dimension[n] * sizeof(int*));
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
		util[ind] = (int*) malloc(n* sizeof(int));
		fscanf(ifp, "%[^\[]",&ch);
		fscanf(ifp, "%c\n",&ch);
		for(int j=0; j<n; j++){
			fscanf(ifp, "%d", util[ind]+j);
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
	vector<string> eqs;
	//construct equation
	int num_eq = sumall*2;
	string eq0,eqj;
	char num[1000];
	//FILE * ofp = fopen("outp.txt", "w");
	//fprintf(ofp, "%d\n",num_eq);
	int combnum=1,shiftcnt=0,shiftmask;
	int* combmask=(int*)malloc(*sizeof(int));
	int* shiftall=(int*)malloc(*sizeof(int));
	for(int comb=0;comb<n;++comb)
	{
		shiftmask = (1<<kn[comb])-1;
		combmask[comb]=shiftmask<<shiftcnt;
		shiftall[comb]=shiftcnt;
		shiftcnt+=kn[comb];
		
	}
	combnum<<=shiftcnt;
	for(count=0;count<combnum;++count;){
		
		int valid = 1;
		for(int cnt=0;cnt<n;++cnt){
			int curmask = combmask(cnt)&count);
			if(curmask==0)count+=(1<<shiftall[n]);
		}
		if(count%n_sz!=n_id)continue;
		for(int i=0;i<n;++i){
			int done = 0;
			eq0.clear();
			eq0 += "r";
			eq0+=to_string(i+1);
			eq0+="1";
			for(int k=0;k<dimension[n]/kn[i];++k) {
				int smallk =k%dimension[i];
				int curk = (k-smallk)*kn[i]+smallk;
				int idc = 0;
				if(util[curk][i]>=0 ){
					eq0+=" + ";
				}
				
				eq0+=to_string(util[curk][i]);
				for(int w=n-1;w>=0;--w){
					currentIndex[w] = curk/dimension[w];
					curk = curk%dimension[w];
					if(w==i)continue;
					eq0+=" * ";
					eq0+="x";		
					eq0+=to_string(w+1);
					eq0+=to_string(currentIndex[w]+1);
				}
				
			}
			
			for(int j=1;j<kn[i];++j){
				eqj.clear();
				//using r
				eqj += "r";
				
				eqj+=to_string(i+1);
				
				eqj+=to_string(j+1);
				//using r
				for(int k=0;k<dimension[n]/kn[i];++k) {
					int smallk =k%dimension[i];
					int curk = (k-smallk)*kn[i] + smallk + dimension[i]*j;
					int idc = 0;
					if(util[curk][i]>=0){//k=1
						eqj+=" + ";
					}
					
					eqj+=to_string(util[curk][i]);
					for(int w=n-1;w>=0;--w){
						currentIndex[w] = curk/dimension[w];
						curk = curk%dimension[w];
						if(w==i)continue;
						eqj+=" * ";
						eqj+="x";
						
						eqj+=to_string(w+1);
						
						eqj+=to_string(currentIndex[w]+1);;
						
					}
					
				}	
				eqj += " - ( "+ eq0 + " );";
				eqs.push_back(string(eqj.c_str()));
			}
		}
		for(int i=0;i<n;++i){
			string eo;
			for(int j=0;j<kn[i];++j){
				eo+="x";
				eo+=to_string(i+1);
				eo+=to_string(j+1);
				eo+=" + ";
				
				
				sprintf(num, "r%d%d * x%d%d ;",i+1,j+1,i+1,j+1);
				eqs.push_back(string(num));
			}
			eo+=" - 1;";
			eqs.push_back(string(eo.c_str()));
		}
		for(int i=0;i<eqs.size();++i){
			printf("%s\n", eqs[i].c_str());
		}
		
		pstream phc; // for running phc
		phc.open("./phc -b");
		if(!phc.is_open())
		{
			printf("Could not open the phc pack executable.\n");
		}
		printf("n=%d\n",eqs.size());
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
		while(phc >> opstr)
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
						varname.push_back(vname);
						
					}
					
					if(done){
						printf("solution:\n");
						for(int i=0;i<num_eq;++i){
							printf("%s: %f\n",varname[i].c_str(),solution[i]);
						}
					}
				}
			}
			
		}
		phc.clear();
		phc.close();
		
		
	}
	gettimeofday(&endTime, NULL);
	totalTime =  (endTime.tv_sec - startTime.tv_sec) * 1000000L;
    totalTime += (endTime.tv_usec - startTime.tv_usec);
	printf("running time: %f sec\n",double(totalTime)/1000000);
	//fclose(ofp);
}
	