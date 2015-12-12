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
using namespace std;
int main(int argc, char** argv)
{
	FILE * ifp = fopen(argv[1], "r");
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
	
	//construct equation
	int num_eq = sumall*2;
	string eq0,eqj;
	char* num[100];
	FILE * ofp = fopen("outp.txt", "w");
	fprintf(ofp, "%d\n",num_eq);
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
		
		printf("cons1");
		for(int j=1;j<kn[i];++j){
			eqj.clear();
			eqj += "r";
			
			eqj+=to_string(i+1);
			
			eqj+=to_string(j+1);
			for(int k=0;k<dimension[n]/kn[i];++k) {
				int smallk =k%dimension[i];
				int curk = (k-smallk)*kn[i] + smallk + dimension[i]*j;
				int idc = 0;
				if(util[curk][i]>=0){
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
			fprintf(ofp, "%s\n",eqj.c_str());
		}
	}
	for(int i=0;i<n;++i){
		string eo;
		for(int j=0;j<kn[i];++j){
			eo+="x";
			eo+=to_string(i+1);
			eo+=to_string(j+1);
			eo+=" + ";
			fprintf(ofp, "r%d%d * x%d%d ;\n",i+1,j+1,i+1,j+1);
			
		}
		eo+=" - 1;";
		fprintf(ofp, "%s\n",eo.c_str());
	}
	fclose(ofp);
}
	