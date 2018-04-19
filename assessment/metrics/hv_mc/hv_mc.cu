#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <string.h>
#include <time.h>

const int maxObjectiveNumber=20;
int objectiveNumber;
//const int tamanhoAmostra=1000000;
const int tamanhoAmostra=1000;


const int maxSize=10000;
double dados[maxSize][maxObjectiveNumber];
__device__ double d_dados[maxSize][maxObjectiveNumber];
__device__ int d_contaAcertos[maxSize];
int contaAcertos[maxSize];

double refPoint;
int tamRef;
curandState* devStates; //devstates of the curand generator

int blockSize;
int nBlocks;

int lerArquivos(char* arquivo);
int main(int argc, char* argv[]){
	if(argc < 3){
		printf("\nuse: hv <file> <ref_point>");
	}
	cudaMalloc ( &devStates, maxSize*sizeof( curandState ) );
	
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

__global__ void calculaHipervolumeAmostragem(curandState* devStates, int objectiveNumber, int tamRef, double refPoint, int cont, int it);
//function to read the files an throw they values to the matrices in the memory
int lerArquivos(char* arquivo){
	int cont=0;
	char tmp[4096];
	int fldcnt=0;
	char arr[1000][1024];
	int recordcnt=0;
	FILE *in=fopen(arquivo,"r");         // open file on command line 
	
	if(in==NULL)
	{
		perror("Error opening the file\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(tmp,sizeof(tmp),in)!=0){ // read a record 
		parse(tmp,(char*)" \t",arr,&fldcnt);   // whack record into fields
		if ((fldcnt != 1 || (!strcmp(arr[0]," \n") && strcmp(arr[0],"#\n") ) ) ){
			
			for(int coluna=0;coluna<fldcnt;coluna++){
// 				for(int i=0;i<strlen(arr[coluna]);i++){
// 					if(arr[coluna][i] == ',')
// 						arr[coluna][i]='.';
// 				}
				dados[recordcnt][coluna]=(double)atof(arr[coluna]);
			}
			recordcnt++;
			objectiveNumber=fldcnt-1;
		}else{
			if(recordcnt > 0){
				tamRef=recordcnt;
				if(tamRef > maxSize){
					printf("\nERROR! Front size is larger than the maximum allowed! (%d)\n", tamRef);
					exit(1);
				}
				
				//blockSize = 8;
				//blockSize=64;
				blockSize=256;
				
				
				nBlocks = tamanhoAmostra/blockSize + (tamanhoAmostra%blockSize == 0?0:1);
				cudaMemcpyToSymbol(d_dados, &dados, sizeof(double)*maxSize*maxObjectiveNumber, 0, cudaMemcpyHostToDevice);
				int tamMul=1000;
				for(int i=0;i<tamMul;i++){
					calculaHipervolumeAmostragem <<< nBlocks, blockSize >>> (devStates, objectiveNumber, tamRef, refPoint, cont, i);
					cudaThreadSynchronize();
				}
				cudaMemcpyFromSymbol(&contaAcertos, d_contaAcertos, sizeof(int)*tamanhoAmostra, 0, cudaMemcpyDeviceToHost);
				int total=0;
				for(int i=0;i<tamanhoAmostra;i++)
					total+=contaAcertos[i];
				
				double volumeTotal = pow(refPoint, objectiveNumber);
				double hv = (double) total / (double) tamanhoAmostra * volumeTotal;
				hv/=tamMul;
				cont++;
				printf("hv(%d) = %.10f\n", cont, hv);
				
				const char* lastError=cudaGetErrorString(cudaGetLastError());
				
				if(strcmp(lastError, "no error")){
					printf("\nCuda error status: %s\n", lastError);
					exit(1);
				}
				
// 				printf("\nac: %d hv: %f",total, hv/tamMul);
// 				printf("\n val: %s ", arr[0]);
// 				printf("rec: %d, fld: %d\n", recordcnt, objectiveNumber);
				recordcnt=0;
				//break;
			}
		}
	}	
	fclose(in);
	//objectiveNumber=fldcnt;
	return recordcnt;
}

__global__ void calculaHipervolumeAmostragem(curandState* devStates, int objectiveNumber, int tamRef, double refPoint, int cont, int it){
// 	for(int i=0;i<tamRef;i++){
// 		for(int j=0;j<objectiveNumber;j++)
// 			printf("%f " , dados[i][j]);
// 		printf("\n");
// 	}
// 	return 0;
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx<tamanhoAmostra){
		if(cont==0 && it ==0)
			curand_init ( clock(), idx, 0, &devStates[idx] );
		
		//curandState localState = devStates[idx];
		//int contaAcertos = 0;
		if(it==0)
			d_contaAcertos[idx]=0;
		__shared__ double gerado[maxObjectiveNumber];
		// Modificar a linha abaixo caso os dados não estejam normalizados
		// ou estejam normalizados em um intervalo de tamanho diferente
		// 1.1 = 2.1 - 1.0, ou seja, limite superior - limite inferior
		
		//no meu caso 1.0 = 1 - 0
		
		// Se os dados tiverem diferentes limites para os diferentes objetivos
		// calcular o volumeTotal pelo produtório das diferenças dos limites....
		//double volumeTotal = pow(refPoint, objectiveNumber);
		bool dominado = true;


		//for (int i = 0; i < tamanhoAmostra; i++) {
			for (int j = 0; j < objectiveNumber; j++) {
				//Gera cada dimensao de um ponto...
				//gerado[j] = 1.0 + 1.1 * random.nextDouble();
				//gerado[j] = refPoint * rand()/RAND_MAX;
				gerado[j] = refPoint * (double)curand_uniform(&devStates[idx]);
				//limites uniformes...
				//gerado[j] = limiteInferior + (limiteSuperior - limiteInferior) * random.nextDouble();
				//diferente limites para cada dimensao...
				//gerado[j] = limiteInferior[j] + (limiteSuperior[j] - limiteInferior[j]) * random.nextDouble();
// 				printf("%f ",gerado[j]);
			}
// 			printf(" (%f)\n", refPoint);

			dominado = false;

			//Verifica se o ponto é dominado ou não...
			for (int k = 0; ((k < tamRef) && (!dominado)); k++) {
				double* temp = d_dados[k];
				
				bool dominadoTemp = true;

				for (int d = 0; d < objectiveNumber; d++) {
					if (temp[d] > gerado[d]) {
						dominadoTemp = false;
					}
// 					printf("%f ",temp[d]);
				}
// 				printf("\n");

				if (dominadoTemp) {
					dominado = true;
				}
			}

			if (dominado) {
				d_contaAcertos[idx]++;//conta o número de pontos dominados...
			}

			//if(i%1000 == 0)
			//System.out.println(i);
		//}
		//devStates[idx] = localState;
		//System.out.println(contaAcertos);

		//O hipervolume é proporcional a quantidade de pontos dominados e ao volume da área amostrada...
		//return (double) contaAcertos / (double) tamanhoAmostra * volumeTotal;
	}

}

