#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <string.h>
#include <time.h>
#include <algorithm>

int sampleSize=1000000;
const int maxObjectiveNumber=20;
const int maxSize=1000000;
double refPoint;
int objectiveNumber;
double dados[maxSize][maxObjectiveNumber];

int lerArquivos(char* arquivo);
double calculateApproximateHypervolume(int size);
int main(int argc, char* argv[]){
	if(argc < 3){
		printf("\nuse: hv <file> <ref_point>");
	}
	unsigned long seed=clock()+time(NULL);
	srand(seed);
	
	refPoint=atof(argv[2]);
	lerArquivos(argv[1]);
}

//function used in the file reading
void parse( char *record, char *delim, char arr[][1024],int *fldcnt){
	char*p=strtok(record,delim);
	int fld=0;
	
	while(p != NULL)
	{
		strcpy(arr[fld],p);
		fld++;
		p=strtok('\0',delim);
	}		
	*fldcnt=fld;
}
//function to read the files an throw they values to the matrices in the memory
int lerArquivos(char* arquivo){
	char tmp[4096];
	int fldcnt=0;
	char arr[1000][1024];
	int recordcnt=0;
	FILE *in=fopen(arquivo,"r");         // open file on command line 
	
	if(in==NULL)
	{
		perror("Error opening the file");
		// 		printf(" (%s) ", arquivo);
		exit(EXIT_FAILURE);
	}
	while(fgets(tmp,sizeof(tmp),in)!=0) // read a record 
	{
		if(tmp[0] != '#' && tmp[0] != '\n'){
			parse(tmp,(char*)" \t",arr,&fldcnt);   // whack record into fields
			if(fldcnt>1)
				objectiveNumber=fldcnt;
			for(int coluna=0;coluna<fldcnt;coluna++){
				for(int i=0;i<strlen(arr[coluna]);i++){
					if(arr[coluna][i] == ',')
						arr[coluna][i]='.';
				}
				dados[recordcnt][coluna]=(double)atof(arr[coluna]);
			}
			recordcnt++;
		}else
			if(recordcnt>1)
				printf("%f\n",calculateApproximateHypervolume(recordcnt));
	}
	fclose(in);
	return recordcnt;
}
double calculateApproximateHypervolume(int size){
	int countDominated=0;
	double generated[objectiveNumber];
	double totalVolume = pow(refPoint, objectiveNumber);
	
	for (int i = 0; i < sampleSize; i++) {
		for (int j = 0; j < objectiveNumber; j++) {
// 			//Generates each dimension of a point
// 			// 			generated[j] = refPoint * (rand()/(double)RAND_MAX);
			generated[j] = (refPoint * (rand()/(double)RAND_MAX));
		}
		
		//Verify if the point is dominated or not
		for (int k = size-1; k >= 0; k--) {

			bool dominatedTmp = true;
			for (int d = 0; d < objectiveNumber; d++) {
				if (dados[k][d] > generated[d]) {
					dominatedTmp = false;
					break;
				}
			}
			if (dominatedTmp) {
				countDominated++;//counts the number of dominated
				break;
			}
		}
	}
// 	printf("%d %d\n",countDominated, size);
	
	return (double) countDominated / (double) sampleSize * totalVolume;
}