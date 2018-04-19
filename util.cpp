//verify if a vector has any nan
//params - the vector to be verified
//		 - the number of size of the vector
bool vectorNan(double* x, int &size){
	for( int i = 0; i < size; i++ ){
		if( x[i] != x[i] ){
			return true;
		}
	}
	return false;
}
//verify if a vector has any nan
//params - the vector to be verified
//		 - the number of size of the vector
bool vectorZero(double* x, int size){
	for( int i = 0; i < size; i++ )
		if( x[i] != 0.0 )
			return false;
	return true;
}
//print a vector
//params - the vector to be printed
//		 - the number of size of the vector
void printVector(double* vec, int vecSize){
	for(int i=0;i<vecSize;i++){
		printf("%f ", vec[i]); //%e
	}
}
//saves a vector to a file
//params - the vector to be printed
//		 - the number of size of the vector
//		 - the file name
void printVectorToFile(double* vec, int vecSize, char* fileName){
	if(vectorNan(vec, vecSize)){
		fprintf(stderr,"\nERROR! NaN on the output file - %s\n", fileName);
		for(int i=0;i<vecSize;i++)
			fprintf(stderr,"%f ", vec[i]);
		
		fprintf(stderr,"\n");
		exit(1);
	}
	
	FILE * pFile;
	pFile = fopen (fileName,"aw"); //aw append and write
	if (pFile==NULL){
		fprintf(stderr,"FILE ERROR\n");
		exit(1);
	}
	
	for(int i=0;i<vecSize;i++)
		fprintf(pFile, "%5.12f ", vec[i]);
	fprintf(pFile, "\n");
	
	fclose (pFile);
}
//insert a blank line in the output file
void insertBlankLine(const char* fileName){
	FILE * pFile;
	pFile = fopen (fileName,"aw"); //aw append and write
	if (pFile==NULL){
		fprintf(stderr,"FILE ERROR\n");
		exit(1);
	}
	fprintf(pFile, "\n\n");
	fclose (pFile);
}
//print a given line to a file
void printToFile(const char* fileName, const char* line){
	FILE * pFile;
	pFile = fopen (fileName,"aw"); //aw append and write
	if (pFile==NULL){
		fprintf(stderr,"FILE ERROR\n");
		exit(1);
	}
	fprintf(pFile, "%s\n", line);
	fclose (pFile);
}
//clear the file
void clearFile(const char* fileName){
	FILE * pFile;
	pFile = fopen (fileName,"w"); //aw append and write
	if (pFile==NULL){
		fprintf(stderr,"FILE ERROR\n");
		exit(1);
	}
	fclose (pFile);
}

//execute a given command in a terminal and return its result
void exec(const char* cmd, char* res) {
	FILE* pipe = popen(cmd, "r");
	if (!pipe) res= (char*) "ERROR";
	char buffer[128]={0};
	while(!feof(pipe)) {
		if(fgets(buffer, 128, pipe) != NULL)
			//res += buffer;
			sprintf(res, "%s%s", res, buffer);
	}
	pclose(pipe);
}

//function used in the file reading
void parse(char *record, const char *delim, char arr[][1024],int *fldcnt){
	char*p=strtok(record,delim);
	int fld=0;
	
	while(p != NULL && strcmp(p, "\n")){
		strcpy(arr[fld],p);
		fld++;
		p=strtok('\0',delim);
	}		
	*fldcnt=fld;
}

// function to read the files an throw they values to the matrices in the memory
int readFile(double **data, const char* file){
	char tmp[4096];
	int fldcnt=0;
	char arr[1000][1024];
	int recordcnt=0;
	FILE *in=fopen(file,"r");// open file on command line 
	
	if(in==NULL){
		fprintf(stderr, "Error opening the file (%s)\n", file);
		exit(EXIT_FAILURE);
	}
	while(fgets(tmp,sizeof(tmp),in)!=0){ // read a record 
		parse(tmp,(char*)" \t",arr,&fldcnt);   // whack record into fields
		if(strcmp(tmp, "\n")){ //if tmp is not just an empty line
			for(int c=0;c<fldcnt;c++){
				for(int i=0;i<(int)strlen(arr[c]);i++){
					if(arr[c][i] == ',')
						arr[c][i]='.';
				}
				if(!strcmp(arr[c], "\n"))//if the field is not a line break
					break;
				data[recordcnt][c]=(double)atof(arr[c]);
			}
			recordcnt++;
		}
	}
	fclose(in);
	return recordcnt;
}

// function to count the columns of the first line of a given file
int countColumns(const char* file){
	char tmp[4096];
	int fldcnt=0;
	char arr[1000][1024];
	FILE *in=fopen(file,"r");// open file on command line 
	
	if(in==NULL){
		fprintf(stderr, "Error opening the file (%s)\n", file);
		exit(EXIT_FAILURE);
	}
	if(fgets(tmp,sizeof(tmp),in)!=0){
		parse(tmp,(char*)" \t",arr,&fldcnt); // whack record into fields
	}else
		fprintf(stderr, "Error counting columns of the file (%s): \"%s\" \n", file, tmp);
	
	return fldcnt;
}

double round_5(const double &in){
// 	int tmp = round(in*100000);
// 	return tmp/100000.0;
	return in;
}
//verify the dominance relation between two vectors - original from jmetal
//return 1 if sol1 dominates sol2, -1 if sol2 dominates sol1, 0 if they do not dominate each other, 2 if they are equal
//params - the objective vector of the first solution to be compared
//		 - the objective vector of the second solution to be compared
//		 - the number of objectives
// return- the number corresponding to the comparison
int dominance(const double *sol1, const double *sol2, const int &objectiveNumber){
	int dominate1 = 0; // dominate1 indicates if some objective of solution1 dominates the same objective in solution2
	int dominate2 = 0; // dominate2 is the complementary of dominate1.    
	int flag; //stores the result of the comparison
	
	double value1, value2;
	for (int i = 0; i < objectiveNumber; i++) {
		//value1 = sol1[i];
		//value2 = sol2[i];
		value1=round_5(sol1[i]);
		value2=round_5(sol2[i]);
		
		if (value1 < value2) {
			flag = -1;
		} else if (value1 > value2) {
			flag = 1;
		} else {
			flag = 0;
		}
		
		if (flag == -1) {
			dominate1 = 1;
		}
		
		if (flag == 1) {
			dominate2 = 1;
		}
	}
	
	if (dominate1 == dominate2) {
		return 0; //No one dominate the other
	}
	if (dominate1 == 1) {
		return 1; // solution1 dominate
	}
	return -1;    // solution2 dominate   
}

// //verify the dominance relation between two vectors - new version created to speed up comparisons especially for many-objective
// //return 1 if sol1 dominates sol2, -1 if sol2 dominates sol1, 0 if they do not dominate each other, 2 if they are equal
// //params - the objective vector of the first solution to be compared
// //		 - the objective vector of the second solution to be compared
// //		 - the number of objectives
// // return- the number corresponding to the comparison
// int dominance(const double *sol1, const double *sol2, const int &objectiveNumber){
// 	bool smaller=false;
// 	bool larger=false;
// 	double value1, value2;
// 	for (int i = 0; i < objectiveNumber; i++) {
// 		value1=round_5(sol1[i]);
// 		value2=round_5(sol2[i]);
// 		
// 		if(value1 < value2)
// 			smaller=true;
// 		else
// 			if(value1 > value2)
// 				larger=true;
// 		
// 		if(smaller && larger){
// 			return 0; //No one dominate the other
// 		}
// 	}
// 	if(smaller)
// 		return 1;//solution 1 dominate
// 	if(larger)
// 		return -1;//solution 2 dominate
// 		
// // 	if(smaller==0 && larger==0)
// // 		return 0;
// }
	
//class created to use the std::sort being able to use parameters
class compareSolutions{
	int obj;
public:
	compareSolutions(int p) : obj(p) {}
	
	bool operator()(const Solution &i, const Solution &j) { //compare two solutions
		return i.objectiveVector[obj] < j.objectiveVector[obj];
	}
};

// order by crowding distances, the smallest goes first
bool crowdingComparatorSol(const Solution &i,const Solution &j) { return (i.crowdingDistance>j.crowdingDistance); }

// order by weighted distances, the smallest goes first
bool weightedDistanceComparatorSol(const Solution &i,const Solution &j) { return (i.weightedDistance<j.weightedDistance); }

// order by neighbor distance, the smallest goes first
bool neighborsComparator(const Neighbor &i,const Neighbor &j) { return (i.distance<j.distance); }

//update the crowding distance of the given set of solutions
//params - the set of solutions to update the crowding distances
//		 - the number of solutions to update the crowding distances
void updateCrowdingDistances(Solution* solutions, int solutionNumber){
	double objetiveMaxn, objetiveMinn, distance, diff;

	for(int s=0;s<solutionNumber;s++){
		solutions[s].crowdingDistance=0;
	}
	
	for(int obj=0;obj<objectiveNumber;obj++){
		
		//shellSort(solutions, obj, solutionNumber);
		std::sort(solutions, solutions+solutionNumber, compareSolutions(obj));
		
		objetiveMinn=solutions[0].objectiveVector[obj];
		solutions[0].crowdingDistance=MAXDOUBLE;
		objetiveMaxn=solutions[solutionNumber-1].objectiveVector[obj];
		solutions[solutionNumber-1].crowdingDistance=MAXDOUBLE;
		diff=objetiveMaxn-objetiveMinn;
		
		if(diff!=0){ //if there is no difference in this objective, don not calculate: its useless
			for(int s=1;s<solutionNumber-1;s++){
				distance=solutions[s+1].objectiveVector[obj]-solutions[s-1].objectiveVector[obj];
				distance/=(objetiveMaxn-objetiveMinn);
				solutions[s].crowdingDistance+=distance;
			}
		}
	}
}

//return the scalar value according to the scalarizing function in use, given the objective and weight vectors of the size of the number of objectives in current use
double getScalarValue(double* objectiveVector, double* weight){
// 	double normalizedSol[objectiveNumber];//normalize
// 	normalizeObjectives(normalizedSol, objectiveVector);//normalize
	double* normalizedSol=objectiveVector;//no normalization
	
	if(!strcmp(archSubSwarms, "pbi"))
		return PBI(normalizedSol, weight);
	
	if(!strcmp(archSubSwarms, "tch"))
		return TCH(normalizedSol, weight);
	
	if(!strcmp(archSubSwarms, "wcp")){
		double sum=0;
		for(int o=0;o<objectiveNumber;o++)
			sum+=(weight[o]*normalizedSol[o])*(weight[o]*normalizedSol[o])*(weight[o]*normalizedSol[o]); //x^3
			return sum;
	}
	if(!strcmp(archSubSwarms, "wsum")){
		double sum=0;
		for(int o=0;o<objectiveNumber;o++)
			sum+=(weight[o]*normalizedSol[o]);
			return sum;
	}
	if(!strcmp(archSubSwarms, "r-ideal")){
		double normalizedWeights[objectiveNumber];
		normalizeObjectives(normalizedWeights, weight);
		return calculateEuclideanDistance(normalizedSol, normalizedWeights, objectiveNumber);
	}
	
	if(objectiveNumber==2 && weight[0]==1){
		return normalizedSol[0];
	}
	
	if(objectiveNumber==2 && weight[1]==1){
		return normalizedSol[1];
	}
	
	//if arrive here, something wrong happened
	fprintf(stderr,"Error on scalarized value calculation, no return!\n");
	exit(1);
}

//update the weighted distance of a given set of solutions
//params - the set of solutions to update the weighted distances
//		 - the number of solutions to update the weighted distances
void updateWeightedDistances(Solution* solutions, int solutionNumber, double* weight){
	for(int s=0;s<solutionNumber;s++){
		solutions[s].weightedDistance=getScalarValue(solutions[s].objectiveVector, weight);
	}
}

//calculate the NWSum of a leader candidate in relation to a solution
double calculateNWSum(Solution &solution, Solution &GBest){
	double sum=0,resParc=0,NWSum=0;
	double solNorm[objectiveNumber], gbNorm[objectiveNumber];
	normalizeObjectives(solNorm, solution.objectiveVector);
	normalizeObjectives(gbNorm, GBest.objectiveVector);
	
	for(int k=0;k<objectiveNumber;k++)
		sum+=solNorm[k];
	
	for(int o=0;o<objectiveNumber;o++){
		if(sum!=0)
			resParc=solNorm[o]/sum;
		
		NWSum+=resParc*gbNorm[o];
	}
	return NWSum;
}

double factorial(const int &n){
	double fat = 1;
	for(int i = n;i>0;i--){
		fat*=i;
	}
	return fat;
}
double combination(const int &m, const int &n){
	if(n==m)
		return 1;
	else{
		double fatM = factorial(m);
		double fatN = factorial(n);
		double fatNM = factorial(m-n);
		return (fatM)/(fatN*fatNM);
	}
}
double calculateSigma(const double &f1, const double &f2){
	double value = (f1*f1) - (f2*f2);
	double denominator = (f1*f1)+ (f2*f2);
	if(denominator!=0){
		return  value/denominator;
	}else
		return 0;
}

void calculateSigmaVector(double* out, double* in){
	int  cont = 0;
	for(int i = 0; i<objectiveNumber-1; i++){
		for(int j = i+1; j<objectiveNumber;j++){
			double obj1=in[i];// normalize(objectiveVector[i], smallerPositions[i], largerPositions[i]);
			double obj2=in[j];// normalize(objectiveVector[j], smallerPositions[j], largerPositions[j]);
			out[cont++] = calculateSigma(obj1, obj2);
		}
	}
}

//calculate the Euclidean distance between two vectors
double calculateEuclideanDistance(double* vector1, double* vector2, int &size){
	double sum = 0;
	for (int i = 0; i < size; i++) {
		//sum += pow(vector1[i]-vector2[i],2);
		sum+= ((vector1[i]-vector2[i])*(vector1[i]-vector2[i]));
	}
	return sqrt(sum);
}
double calculateTchebycheffDistance(double* vector1, double* vector2, int size){
	double max=0;
	for (int i = 0; i < size; i++)
		max=std::max(max, fabs(vector1[i]-vector2[i]));
	return max;
}
double calculateMinkowskiDistance(double* vector1, double* vector2, double p, int size){
	double sum = 0;
	for (int i = 0; i < size; i++) {
		sum += pow(fabs(vector1[i]-vector2[i]),p);
	}
	return pow(sum, 1.0/p);
}
int compute_b_mga(Solution* solutions, int repActualSize){
	double max_value = 0;
	// 	for(int i=0;i<repActualSize;i++){
	// 		for(int j = 0; j<objectiveNumber;j++){
	// 			double value=normalize(solutions[i].objectiveVector[j], smallerObjs[j], largerObjs[j]);
	// // 			double value=solutions[i].objectiveVector[j];
	// 			if(value> max_value)
	// 				max_value = value;
	// 		}
	// 	}
	//our values are normalized, hence the maximum value possible is 1
	max_value=1;
	
	return (int) floor(log2(max_value)) + 1;
}
void box_mga(const Solution &solution, double* output, const int b, const int objectiveNumber){
	double solNorm[objectiveNumber];
	normalizeObjectives(solNorm, solution.objectiveVector);
	
	for(int i = 0; i<objectiveNumber; i++){
		double z=solNorm[i];
		
		//output[i] = floor(z / pow(2.0, b) + 0.5); //this 0.5 is taken from the code of Andre, I do not know why it is there
		//output[i]=(int) floor(z*pow(2,-b));
		output[i] = floor(z / pow(2.0, b));
	}
}


double normalize(const double &valor,const double &min,const double &max){
	double diff=max-min;
	if(diff > 0){
		double result=((valor-min)/diff);
		// 		printf("antes: %f, depois: %f\n", valor, result);
// 		if(result >=0 && result <= 1)
			return result;
// 		else{
// 			fprintf(stderr,"min: %f, max: %f, valor: %f, result %f\n", min, max, valor, result);
// 			return valor-min;
// 		}
	}else{
		if(valor > 0 && (max != min)){
			fprintf(stderr,"min: %f, max: %f, valor: %f\n", min, max, valor);
			// 			exit(1);
		}
		return (valor-min);//((valor-min)/(max-min));
	}
	//return valor;
	//normalize(h_swarm[0].particles[i].solution.objectiveVector[p], h_swarm[0].repository.smallerObjs[p], h_swarm[0].repository.largerObjs[p])
}

//update the largest and smallest objective vectors from the repositories to normalize the solutions
//the vectors are updated from the population in the evaluation stage
void updateMaxMinObjs(){
	for(int s=0;s<swarmNumber;s++){
		for(int i=0; i< swarms[s].repository.getActualSize(); i++){
			for(int o=0;o<objectiveNumber;o++){
				maxObjectives[o]=std::max(maxObjectives[o], swarms[s].repository.getSolution(i).objectiveVector[o]);
				minObjectives[o]=std::min(minObjectives[o], swarms[s].repository.getSolution(i).objectiveVector[o]);
			}
		}
	}
}

void normalizeObjectives(double* out, const double* in){
	for(int o=0;o<objectiveNumber;o++)
		out[o]=normalize(in[o], minObjectives[o], maxObjectives[o]);
}

// void normalizeDecision(double* out, const double* in){
// 	for(int d=0;d<decisionNumber;d++)
// 		out[d]=normalize(in[d], inferiorPositionLimit[d], superiorPositionLimit[d]);
// }

//update the largest and smallest decision vectors to truncate the solutions in I-Multi
void updateLargestSmallestDecVectors(Swarm &swarm){
	for(int d=0;d<decisionNumber;d++){
		swarm.repository.largestDecision[d]=MAXDOUBLE*-1;
		swarm.repository.smallestDecision[d]=MAXDOUBLE;
	}
	
	//update the largest and smallest positions vectors from the repository
	for(int i=0; i< swarm.repository.getActualSize(); i++){
		for(int d=0;d<decisionNumber;d++){
			swarm.repository.largestDecision[d]=std::max(swarm.repository.largestDecision[d], swarm.repository.getSolution(i).decisionVector[d]);
			swarm.repository.smallestDecision[d]=std::min(swarm.repository.smallestDecision[d], swarm.repository.getSolution(i).decisionVector[d]);
		}
	}
}

void weightClustering(Solution* solutions, int solNumber){
	//calculates the largest and smaller positions found so far in the search to normalize
// 	double max[objectiveNumber];
// 	double min[objectiveNumber];
	bool used[solNumber];
	
// 	for(int i=0;i<objectiveNumber;i++){
// 		max[i]=MAXDOUBLE*-1;
// 		min[i]=MAXDOUBLE;
// 	}
// 	for(int i=0;i<solNumber;i++){
// 		used[i]=false;
// 		for(int o=0;o<objectiveNumber;o++){
// 			max[o]=std::max(max[o],solutions[i].objectiveVector[o]);
// 			min[o]=std::min(min[o],solutions[i].objectiveVector[o]);
// 		}
// 	}
	
// 	printf("\n");
// 	for(int o=0;o<objectiveNumber;o++){
// 		if(max[o] == globalLargerObjs[o])
// 			printf("\nigual max");
// 		else
// 			printf("\ndiferente max %f != %f", max[o], globalLargerObjs[o]);
// 		
// 		if(min[o] == globalSmallerObjs[o])
// 			printf("\nigual min");
// 			else
// 				printf("\ndiferente min %f != %f", min[o], globalSmallerObjs[o]);
// 	}
// 	printf("\n");
	
	if(solNumber < swarmNumber){ //if there is not even one solution per swarm, let them empty
		for(int i=0;i<swarmNumber;i++){
			for(int j=0;j<swarms[i].getSize();j++){
				swarms[i].particles[j].solution.initialize();//in this case, generate a random solution
				swarms[i].particles[j].solution.evaluate();//in this case, evaluate the random solution before putting it in the repository
				swarms[i].repository.add(swarms[i].particles[j].solution);
			}
		}
	}
	
	//first step, each cluster selects a particle and removes it from the pool, but a particle can be selected by more than one cluster
	for(int j=0;j<swarmNumber;j++){
		double minValue=MAXDOUBLE;
		int index=-1;
		for(int i=0;i<solNumber;i++){
// 			if(!used[i]){ //if the solutions was used, do not use it again
				double dist=getScalarValue(solutions[i].objectiveVector, swarms[j].repository.weight);
				//after calculating the distance regarding this particle, test if this is the best
				if(dist < minValue){
					minValue=dist;
					index=i;
				}
			}
// 		}
		if(solNumber > 0){
			if(index==-1){
				fprintf(stderr,"\nWEIGHT CLUSTERING ERROR %f\n", minValue);
				exit(1);
			}
			
			swarms[j].repository.add(solutions[index]);
			used[index]=true;
		}
	}
	
	//grab each solution from the initial set and atribute it to the better cluster
	for(int i=0;i<solNumber;i++){
		if(!used[i]){ //if the solutions was used, do not use it again
		
			double minValue=MAXDOUBLE;
			int index=-1;
			
			//chech if there exists empty repositories
// 			bool existsEmpty=false;
// 			for(int j=0;j<swarmNumber;j++)
// 				if(multiswarm[j].repository.getActualSize()==0){
// 					existsEmpty=true;
// 					break;
// 				}
			
			for(int j=0;j<swarmNumber;j++){
				double dist=getScalarValue(solutions[i].objectiveVector, swarms[j].repository.weight);
			
				//after calculating the distance regarding this cluster, test if this is the better
				if(dist < minValue){
	// 				if(!existsEmpty || multiswarm[j].repository.getActualSize()==0){ //if there are no empty repositories, or if this is an empty repository, allow to add
						minValue=dist;
						index=j;
	// 				}
				}
			}
			
	// 		printf("%d, %f, %d, %d \n", index, minValue, i, multiswarm[index].repository.getActualSize());
			
			if(index==-1){
				fprintf(stderr,"\nWEIGHT CLUSTERING ERROR %f\n", minValue);
				exit(1);
			}
			
			swarms[index].repository.add(solutions[i]);
			used[i]=true;//not needed, but just in case
		}
	}
	// 	//setting the centroid only to ensure retrocompatibility
	for(int i=0;i<swarmNumber;i++){
// 		swarms[i].repository.organize();
		for(int s=0; s<decisionNumber; s++ ){
			double sum=0;
			for(int p=0;p<swarms[i].repository.getActualSize();p++)
				sum+=swarms[i].repository.getSolution(p).decisionVector[s];
			swarms[i].centroid[s]=(sum/swarms[i].repository.getActualSize());
		}
	}
}

//find the extreme solutions and the solution closer to the ideal
//params - solutions of the front
//		 - extreme solutions (output)
//		 - the number of objectives
// void obtainIdealSolution(Solution* solutions, Solution ideal, int objectiveNumber, int repActualSize){
// // 	for(int i=0;i<objectiveNumber;i++)
// // 		ideal.objectiveVector[i]=MAXDOUBLE;
// 	//Obtain the extreme solutions and calculates the ideal
// 	for(int i=0;i<repActualSize;i++){
// 		for(int j = 0; j<objectiveNumber;j++){
// 			if(solutions[i].objectiveVector[j]<=ideal.objectiveVector[j]){
// 				ideal.objectiveVector[j] = solutions[i].objectiveVector[j];
// 			}
// 		}
// 	}
// }

//efficiently calculate a weighted tchebycheff distance
double r2_max(double* V, double* A, int objectiveNumber){
	double max=(MAXDOUBLE*-1);
	double normA[objectiveNumber];
	normalizeObjectives(normA, A);
	
	for(int j=0;j<objectiveNumber;j++){
		//double valor=V[j]*fabs(ponto[j]-A[j]);
		//double valor=V[j]*fabs(0-A[j]); //considering the origin as reference point
		//double valor= normalize(V[j], smallerObjs[j], largerObjs[j])*normalize(A[j], smallerObjs[j], largerObjs[j]);
		double valor= V[j] * normA[j]; //considering the origin as reference point
		if( valor > max)
			max=valor;
	}
	return max;
}

// //efficiently calculate a asf metric from "Improved Metaheuristic Based on the R2 Indicator for Many-Objective Optimization"
// double r2_max(double* V, double* A, int objectiveNumber, double* smallerObjs, double* largerObjs){
// 	double max=(MAXDOUBLE*-1);
// 	double value;
// 	for(int j=0;j<objectiveNumber;j++){
// 		//double valor=V[j]*fabs(ponto[j]-A[j]);
// 		//double valor=V[j]*fabs(0-A[j]); //considering the origin as reference point
// 		//double valor= normalize(V[j], smallerObjs[j], largerObjs[j])*normalize(A[j], smallerObjs[j], largerObjs[j]);
// 		if(V[j] > 0)
// 			value= normalize(A[j], smallerObjs[j], largerObjs[j])/V[j]; //considering the origin as reference point
// 		else
// 			value= normalize(A[j], smallerObjs[j], largerObjs[j])/1e-10; //considering the origin as reference point
// 			
// 		if( value > max)
// 			max=value;
// 	}
// 	return max;
// }

double r2_min(double* V, Repository &rep){
	double min=MAXDOUBLE;
	for(int i=0;i<rep.getActualSize();i++){
		double valor=r2_max(V, rep.getSolution(i).objectiveVector, objectiveNumber);
		if(valor < min)
			min=valor;
	}
	return min;
}
double calculateR2(Repository &rep, double **referencia, int tamRef){
	if(tamRef <= 0){
		fprintf(stderr,"\nNumber of reference points for R2 calculation is %d, are you reading them?\n", tamRef);
		exit(1);
	}
	
	double R2=0;
	for(int i=0;i<tamRef;i++){
		R2+=r2_min(referencia[i], rep);
	}
	R2*=(1.0/tamRef);
	return R2;
}

double calculateApproximateHypervolume(Repository &rep){
	double refPoint=1.01;
	int countDominated=0;
	double generated[objectiveNumber];
	double totalVolume = pow(refPoint, objectiveNumber);
	int sampleSize=1000000;
	double diff[objectiveNumber];
	for (int j = 0; j < objectiveNumber; j++)
		diff[j]=maxObjectives[j]-minObjectives[j];
	
	for (int i = 0; i < sampleSize; i++) {
		for (int j = 0; j < objectiveNumber; j++) {
			//Generates each dimension of a point
// 			generated[j] = refPoint * (rand()/(double)RAND_MAX);
			generated[j] = ((refPoint * (rand()/(double)RAND_MAX))*diff[j])+minObjectives[j];
		}
		
		//Verify if the point is dominated or not
		for (int k = 0; k < rep.getActualSize(); k++) {
// 			double normalized[objectiveNumber];
// 			for(int o=0;o<objectiveNumber;o++)
				//normalized[o]=normalize(rep.getSolution(k).objectiveVector[o], smallerObjs[o], largerObjs[o]);
// 			if(dominance(normalized, generated, objectiveNumber) == 1) //if the generated solution is dominated by one from the repository
// 			if(dominance(rep.getSolution(k).objectiveVector, generated, objectiveNumber) == 1){ //if the generated solution is dominated by one from the repository
// 				countDominated++;//counts the number of dominated
// 				break;
// 			}
			bool dominatedTmp = true;
			for (int d = 0; d < objectiveNumber; d++) {
				if (rep.getSolution(k).objectiveVector[d] > generated[d]) {
					dominatedTmp = false;
				}
			}
			if (dominatedTmp) {
				countDominated++;//counts the number of dominated
				break;
			}
		}
	}
// 	printf("%d\n",countDominated);
	return (double) countDominated / (double) sampleSize * totalVolume;
}

// // double calculateHypervolume2(Repository &rep, double refPoint);
// double calculateHypervolume(Repository &rep, double refPoint){
// 	char tmpFile[1000]={0};
// 	char command[1000]={0};
// 	char resCommand[1000]={0};
// 	char refPoints[1000]={0};
// 	double hv=0;
// 	for(int o=0;o<objectiveNumber;o++)
// 		sprintf(refPoints, "%s%f ", refPoints, refPoint);
// 	
// 	sprintf(tmpFile, "/tmp/hv_tmp_%d.txt",getpid());
// 	sprintf(command, "./assessment/metrics/hv/wfg %s %s | head -1 | cut -d '=' -f 2",tmpFile, refPoints);
// 	clearFile((char*)tmpFile); //create an empty temp file
// 	printToFile((char*)tmpFile, (char*)"#");
// 	for(int i=0;i<rep.getActualSize();i++){
// 		double normalized[objectiveNumber];
// 		normalizeObjectives(normalized, rep.getSolution(i).objectiveVector);
// 		printVectorToFile(normalized, objectiveNumber, (char*)tmpFile);
// 	}
// 	printToFile((char*)tmpFile, (char*)"#");
// 	exec((char*)command, resCommand);
// 	hv=atof(resCommand);
// 	
// // 	printf("%f == %f\n", hv, calculateApproximateHypervolume(rep, smallerObjs, largerObjs));
// 	
// // 	double hv2=calculateHypervolume2(rep, 1.01);
// // 	if((int)(hv*10000000) != (int)(hv2*10000000))
// // 		printf("%f != %f\n", hv, hv2);
// 	
// 	if(hv == 0)
// 		fprintf(stderr,"\nWARNING! HYPERVOLUME OF A SET: %f -- COMMAND: %s\n", hv, resCommand);
// 	
// 	return hv;
// }

////Functions from the WFG hypervolume 
#include "additionalCode/wfg.h"
#include "additionalCode/wfg.c"
double calculateHypervolume2(Repository &rep, double refPoint){
// 	double refPoint=1.01;
	// allocate memory
	int maxdepth = objectiveNumber - 2; 
	fs = (FRONT*)malloc(sizeof(FRONT) * maxdepth);
	for (int i = 0; i < maxdepth; i++){
		fs[i].points = (POINT*)malloc(sizeof(POINT) * rep.getActualSize()); 
			for (int j = 0; j < rep.getActualSize(); j++) 
				fs[i].points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * (objectiveNumber - i - 1));
	}
	FRONT *f = (FRONT*)malloc(sizeof(FRONT));
	f->points = (POINT*)malloc(sizeof(POINT) * rep.getActualSize());
	for (int j = 0; j < rep.getActualSize(); j++) 
		f->points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * (objectiveNumber));
	///////////////
	for(int i=0;i<rep.getActualSize();i++){
		double normSol[objectiveNumber];
		normalizeObjectives(normSol, rep.getSolution(i).objectiveVector);
		for(int o=0;o<objectiveNumber;o++){
// 			normRep[i][o]=fabs(normSol[o]-refPoint);
			f->points[i].objectives[o]=fabs(normSol[o]-refPoint);
		}
	}
	f->nPoints=rep.getActualSize();
	n=objectiveNumber;
	safe = 0;
	fr = 0;   // current depth 

	double HV=hv(*f);

	////////dealocate memory
	for (int i = 0; i < maxdepth; i++){
		for (int j = 0; j < rep.getActualSize(); j++) 
			free(fs[i].points[j].objectives);
		free(fs[i].points);
		
	}
	free(fs);
	
	for (int j = 0; j < rep.getActualSize(); j++) 
		free(f->points[j].objectives);

	free(f->points);
	free(f);
	///////

	return HV;
}

//comparator if two vectors are equal
//param - the vectors to be compared to
bool isEqual(double* vec1, double* vec2, int size){
	for(int i=0;i<size;i++){
		if(vec1[i] != vec2[i])
			return false;
	}
	return true;
}

// //calculate the product between two vectors
// double vectorProduct(double* vec1, double* vec2, int size){
// 	double sum = 0;
// 	for(int i=0; i<size; i++)
// 		sum+= vec1[i]*vec2[i];
// 	return sum;
// }

//calculate the inner product
double innerProduct(const double* vec1, const double* vec2, int size) {
	double sum = 0;
	for (int i = 0; i < size; i++)
		sum += vec1[i] * vec2[i];
	return sum;
}

//calculate the norm of a vector
double vectorNorm(double* x, int &size){
	double sum = 0;
	for(int i=0;i<size;i++)
		sum = sum + x[i]*x[i];
	return sqrt(sum);
}
//calculate the PBI of a vector given a weight vector
double PBI(double* sol, double* weight){
// // // // 	//some things are different from the common equation because our problems are minimization (see the small letters in PBI paper)
// // // // 	double tmp[objectiveNumber];
// // // // 	for(int i=0;i<objectiveNumber;i++)
// // // // 		tmp[i]=(sol[i]-0.0)*weight[i];//considering the origin as ideal point (r* - f) * w
// // // // 	double d1= vectorNorm( tmp, objectiveNumber )/vectorNorm(weight, objectiveNumber);
// // // // 	
// // // // 	for(int i=0;i<objectiveNumber;i++)
// // // // 		tmp[i]=sol[i]-(0.0+d1*weight[i]);
// // // // 	double d2 = vectorNorm(tmp, objectiveNumber);
// // // // 	
// // // // 	return d1 + 5*d2;
	
	
	//from MOEA/DD
	// normalize the weight vector (line segment)
	double lambda[objectiveNumber];
	double nd = vectorNorm(weight, objectiveNumber);
	for (int i = 0; i < objectiveNumber; i++)
		lambda[i] = weight[i] / nd;

	double realA[objectiveNumber];
	double realB[objectiveNumber];

	// difference between current point and reference point
	for (int n = 0; n < objectiveNumber; n++)
// 		realA[n] = (sol[n] - 0.0);
		realA[n] = (sol[n] - minObjectives[n]);

	// distance along the line segment
	double d1 = abs(innerProduct(realA, lambda, objectiveNumber));

	// distance to the line segment
	for (int n = 0; n < objectiveNumber; n++)
// 		realB[n] = (sol[n] - (0.0 + d1 * lambda[n]));
		realB[n] = (sol[n] - (minObjectives[n] + d1 * lambda[n]));
		
	double d2 = vectorNorm(realB, objectiveNumber);

	return d1 + 5.0 * d2;
	
}
//calculate the weighted Tchebycheff of a vector given a weight vector
double TCH(double* sol, double* weight){
	double max=0;
	for (int i = 0; i < objectiveNumber; i++)
// 		max=std::max(max, fabs(sol[i]-0.0)*weight[i] );//considering the origin as ideal point (f-r*)
		max=std::max(max, fabs(sol[i]-minObjectives[i])*weight[i] );//considering the ideal
	return max;
}

void updateContributingHypervolume(Repository &rep){
// 	rep.organize();
	
	double previousHV=0;
	if(objectiveNumber <=8)
		previousHV=calculateHypervolume2(rep, 1.01);
	else
		previousHV=calculateApproximateHypervolume(rep);
		
// 	for(int i = rep.getActualSize()-1; i>=0;i--){
	for(int i = 0; i<rep.getActualSize();i++){
		Solution tmp;
		tmp=rep.getSolution(i);
// 		memcpy(&tmp, &rep.getSolution(i), sizeof(Solution));
		
		rep.exclude(i);
// 		rep.organize();
		double currentHV=0;
		if(rep.getActualSize() > 1){
			if(objectiveNumber <=8){
				currentHV=calculateHypervolume2(rep, 1.01);//calculateHypervolume2 uses the WFG hv code within this algorithm
// 				double currentHV2=calculateHypervolume2(rep, 1.01);
// 				if(abs(currentHV-currentHV2) > 0.0000000001){
// 					printf("different: hv1: %f hv2: %f\n", currentHV, currentHV2);
// 					exit(1);
// 				}else
// 					if(currentHV != currentHV2)
// 						printf("difference: %.30f\n", abs(currentHV-currentHV2));
				
			}
			else
				currentHV=calculateApproximateHypervolume(rep);
		}
		tmp.crowdingDistance=previousHV-currentHV; //stores in the crowding distance field -- if HV decreased after I removed the solution, it is important
		
// 		rep.add(tmp);
		rep.forceAdd(tmp); //force add to allow being used as archiver, otherwise enter in an infinite loop
		
// 		if(tmp.crowdingDistance < 0)
// 			printf("sol: %d - prev hv: %f, currHV: %f, diff: %f updatedHV: %f\n",i, previousHV, currentHV, tmp.crowdingDistance, calculateHypervolume(rep));
	}
}

void updateContributingR2(Repository &rep){
	double previousR2=calculateR2(rep, reference, refSize);
	
// 	for(int i = rep.getActualSize()-1; i>=0;i--){
	for(int i = 0; i<rep.getActualSize();i++){
		Solution tmp;
		tmp=rep.getSolution(i);
// 		memcpy(&tmp, &rep.getSolution(i), sizeof(Solution));
		
		rep.exclude(i);
// 		rep.organize();
		if(rep.getActualSize() > 1)
			tmp.crowdingDistance=calculateR2(rep, reference, refSize)-previousR2; //if the value of the r2 increased after I take this solution, this solution is important
// 		rep.add(tmp);
		rep.forceAdd(tmp); //force add to allow being used as archiver, otherwise enter in an infinite loop
	}
}

// //calculate the determinant of a matrix in a naive way, very costly
// double matrixDeterminant(double **a,int n) { //http://stackoverflow.com/questions/21220504/matrix-determinant-algorithm-c
// 	int p, h, k, i, j;
// 	double temp[decisionNumber][decisionNumber], det=0;
// 	if(n==1) {
// 		return a[0][0];
// 	} else if(n==2) {
// 		det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
// 		return det;
// 	} else {
// 		for(p=0;p<n;p++) {
// 			h = 0;
// 			k = 0;
// 			for(i=1;i<n;i++) {
// 				for( j=0;j<n;j++) {
// 					if(j==p) {
// 						continue;
// 					}
// 					temp[h][k] = a[i][j];
// 					k++;
// 					if(k==n-1) {
// 						h++;
// 						k = 0;
// 					}
// 				}
// 			}
// 			det=det+a[0][p]*pow(-1,p)*matrixDeterminant(temp,n-1);
// 		}
// 		return det;
// 	}
// }
// //multiplies two NxN matrices
// void matrixMultiplication(double **mat1, double **mat2, double **output, int n){
// 	for(int i=0;i<n;i++){
// 		for(int j=0;j<n;j++){
// 			double sum=0;
// 			for(int k=0;k<n;k++){
// 				sum += mat1[i][k]*mat2[k][j];
// 			}
// 			output[i][j] = sum;
// 		}
// 	}
// }
// //calculate the outer product of two vectors (vec1 x vec2^t)
// void outerProduct(double* vec1, double* vec2, double **output, const int &n){
// 	for(int i=0;i<n;i++){
// 		for(int j=0;j<n;j++){
// 			output[i][j]=vec1[i]*vec2[j];
// 		}
// 	}
// }

double kullbackLeiblerDivergence(const CMAESModel model0, const CMAESModel model1){
//https://en.wikipedia.org/wiki/Kullback–Leibler_divergence
	double invCovMat1[decisionNumber][decisionNumber];
	double covMat0[decisionNumber][decisionNumber];
	
	for(int i=0;i<decisionNumber;i++){//calculate the inverse of the covariance matrix as: Σ^−1 = B D^-2 B^t
		for(int j=0;j<decisionNumber;j++){
			invCovMat1[i][j]=0;
			for(int k=0;k<decisionNumber;k++){
				invCovMat1[i][j]+=model1.B[i][k]*(1.0/(model1.D[k]*model1.D[k]))*model1.B[j][k];
			}
		}
	}
	//complete the covariance matrix, since originally it is only a triangular matrix
	for(int i=0;i<decisionNumber;i++){
		for(int j=0;j<decisionNumber;j++){
			if(j<i)
				covMat0[i][j]=model0.C[i][j];
			else
				covMat0[i][j]=model0.C[j][i];
		}
	}
	
	double tmpMatrix[decisionNumber][decisionNumber];
	for(int i=0;i<decisionNumber;i++){//multiplicates Σ_1^−1 * Σ_0
		for(int j=0;j<decisionNumber;j++){
			tmpMatrix[i][j]=0;
			for(int k=0;k<decisionNumber;k++){
				tmpMatrix[i][j]+=invCovMat1[i][k]*covMat0[k][j];
			}
		}
	}
	double part1=0;
	for(int i=0;i<decisionNumber;i++)
		part1+=tmpMatrix[i][i];//tr(Σ_1^−1*Σ_0)
	
	double tmp[decisionNumber], part2=0;
	//tmp=(µ_1 − µ_0)' * Σ_1^−1
	for(int i=0;i<decisionNumber;i++){
		tmp[i]=0;
		for(int j=0;j<decisionNumber;j++)
			tmp[i]+=(model1.mean[j]-model0.mean[j])*invCovMat1[j][i];
		
		//second part  -- tmp * (µ1 − µ0)
		part2+=tmp[i]*(model1.mean[i]-model0.mean[i]);
	}
	part2-=decisionNumber;//(µ_1 − µ_0)' * Σ_1^−1 * (µ_1 − µ_0)-k
	
	double part3=log(model1.det/model0.det);//ln(Σ_1^−1 / Σ_0)

	return (part1+part2+part3)/2.0;
}

//located here due to cross referencing issues on compilation time -- archiver that removes from the repository the solution with the largest log-likelihood
bool Repository::largestLikelihoodArchiver(Solution &candidate){
	insert(candidate);
// 	organize();
	
	double likelihood=MAXDOUBLE*-1;
	int index=-1;
	
	for(int i=0;i<actualSize;i++){
		double lik=0;
		//lik=calculateEuclideanDistance(norm, weight, objectiveNumber);
		lik=outerSwarm->model.logLikelihood(solutions[i].decisionVector);
// 		lik=logLikelihood(solutions[i].decisionVector, *outerSwarm);
		if(lik > likelihood){
			likelihood=lik;
			index=i;
		}
	}
	if(index==-1){
		fprintf(stderr,"\nLARGEST LIKELIHOOD ARCHIVER ERROR %f\n", likelihood);
		exit(1);
	}
	bool ret=true;
	if(candidate.isEqual(solutions[index]))
		ret=false;
	
	exclude(index);
	
	return ret;
}
//located here due to cross referencing issues on compilation time -- archiver that removes from the repository the solution with the smallest log-likelihood
bool Repository::smallestLikelihoodArchiver(Solution &candidate){
	insert(candidate);
// 	organize();
	
	double likelihood=MAXDOUBLE;
	int index=-1;
	
	for(int i=0;i<actualSize;i++){
		double lik=0;
		//lik=calculateEuclideanDistance(norm, weight, objectiveNumber);
		lik=outerSwarm->model.logLikelihood(solutions[i].decisionVector);
// 		lik=logLikelihood(solutions[i].decisionVector, *outerSwarm);
		if(lik < likelihood){
			likelihood=lik;
			index=i;
		}
	}
	if(index==-1){
		fprintf(stderr,"\nLARGEST LIKELIHOOD ARCHIVER ERROR %f\n", likelihood);
		exit(1);
	}
	bool ret=true;
	if(candidate.isEqual(solutions[index]))
		ret=false;
	
	exclude(index);
	
	return ret;
}
//located here due to cross referencing issues on compilation time -- archiver that removes from the repository the solution with the largest euclidean distance to the mean (in decision space)
bool Repository::largestDistanceArchiver(Solution &candidate){
	insert(candidate);
// 	organize();
	
	double likelihood=MAXDOUBLE*-1;
	int index=-1;
	
	for(int i=0;i<actualSize;i++){
		double lik=0;
// 		double decNorm[decisionNumber];
// 		normalizeDecision(decNorm, solutions[i].decisionVector);
		double *decNorm=solutions[i].decisionVector;
		
		
		lik=calculateEuclideanDistance(decNorm, outerSwarm->model.mean, decisionNumber);
		if(lik > likelihood){
			likelihood=lik;
			index=i;
		}
	}
	if(index==-1){
		fprintf(stderr,"\nLARGEST LIKELIHOOD ARCHIVER ERROR %f\n", likelihood);
		exit(1);
	}
	bool ret=true;
	if(candidate.isEqual(solutions[index]))
		ret=false;
	
	exclude(index);
	
	return ret;
}
//located here due to cross referencing issues on compilation time -- archiver that removes from the repository the solution with the smallest euclidean distance to the mean (in decision space)
bool Repository::smallestDistanceArchiver(Solution &candidate){
	insert(candidate);
// 	organize();
	
	double likelihood=MAXDOUBLE;
	int index=-1;
	
	for(int i=0;i<actualSize;i++){
		double lik=0;
// 		double decNorm[decisionNumber];
// 		normalizeDecision(decNorm, solutions[i].decisionVector);
		double *decNorm=solutions[i].decisionVector;
		
		lik=calculateEuclideanDistance(decNorm, outerSwarm->model.mean, decisionNumber);
		if(lik < likelihood){
			likelihood=lik;
			index=i;
		}
	}
	if(index==-1){
		fprintf(stderr,"\nLARGEST LIKELIHOOD ARCHIVER ERROR %f\n", likelihood);
		exit(1);
	}
	bool ret=true;
	if(candidate.isEqual(solutions[index]))
		ret=false;
	
	exclude(index);
	
	return ret;
}
//changes a given string from lower case to upper case characters
void toUpperCase(char *str){
	int i=0;
	while (str[i]){
		char c=str[i];
		str[i]=toupper(c);
		i++;
	}
}

//*************************************************************************BEGIN OF CLUSTERING QUALILY INDICATORS**************************************************//

int clusteringGetVectorSize(){
	if(!strcmp(clusteringType, "decision"))
		return decisionNumber;
	if(!strcmp(clusteringType, "objectives"))
		return objectiveNumber;
	if(!strcmp(clusteringType, "both"))
		return decisionNumber+objectiveNumber;
	
	fprintf(stderr,"\nDEFINING CLUSTERING SIZE ERROR \n");
	exit(1);
}

void clusteringNormalizePoint(const Solution &in, double* out, double* min, double* max){
	int size=clusteringGetVectorSize();
	
	if(!strcmp(clusteringType, "decision"))
		for(int j=0;j<size;j++)
			out[j]=normalize(in.decisionVector[j], min[j], max[j]);
	
	if(!strcmp(clusteringType, "objectives"))
		for(int j=0;j<size;j++)
			out[j]=normalize(in.objectiveVector[j], min[j], max[j]);
	
	if(!strcmp(clusteringType, "both")){
		for(int j=0;j<size;j++){
			if(j<objectiveNumber)
				out[j]=normalize(in.objectiveVector[j], min[j], max[j]);
			else
				out[j]=normalize(in.decisionVector[j-objectiveNumber], min[j], max[j]);
		}
	}
}

//calculate the smaller and larger objective vectors of multiple swarms
void multiSwarmLimits(Swarm* multiswarm, double* min, double* max){
	int size=clusteringGetVectorSize();
	for(int i=0;i<size;i++){
		max[i]=MAXDOUBLE*-1;
		min[i]=MAXDOUBLE;
	}
	
	if(!strcmp(clusteringType, "decision")){
		for(int s=0;s<swarmNumber;s++){
			for(int i=0;i<multiswarm[s].repository.getActualSize();i++){
				for(int o=0;o<size;o++){
					max[o]=std::max(max[o],multiswarm[s].repository.getSolution(i).decisionVector[o]);
					min[o]=std::min(min[o],multiswarm[s].repository.getSolution(i).decisionVector[o]);
				}
			}
		}
	}
	
	if(!strcmp(clusteringType, "objectives")){
		for(int s=0;s<swarmNumber;s++){
			for(int i=0;i<multiswarm[s].repository.getActualSize();i++){
				for(int o=0;o<size;o++){
					max[o]=std::max(max[o],multiswarm[s].repository.getSolution(i).objectiveVector[o]);
					min[o]=std::min(min[o],multiswarm[s].repository.getSolution(i).objectiveVector[o]);
				}
			}
		}
	}
	
	if(!strcmp(clusteringType, "both")){
		for(int s=0;s<swarmNumber;s++){
			for(int i=0;i<multiswarm[s].repository.getActualSize();i++){
				for(int o=0;o<size;o++){
					if(o<objectiveNumber){
						max[o]=std::max(max[o],multiswarm[s].repository.getSolution(i).objectiveVector[o]);
						min[o]=std::min(min[o],multiswarm[s].repository.getSolution(i).objectiveVector[o]);
					}else{
						max[o]=std::max(max[o],multiswarm[s].repository.getSolution(i).decisionVector[o-objectiveNumber]);
						min[o]=std::min(min[o],multiswarm[s].repository.getSolution(i).decisionVector[o-objectiveNumber]);
					}
				}
			}
		}
	}
}

//used in silhouette analysis
double dissimilarity(const Solution &point, Swarm &swarm, bool ownCluster, double* min, double* max){ //dissimilarity measure used in the silhouette analysis
	if(swarm.repository.getActualSize()==1 && ownCluster) //if i am the only member of this cluster, size 0 to myself
		return 0;
	
	double sum=0;
	int size=clusteringGetVectorSize();
	double normalizedPoint[size];
	double normalizedClusterPoint[size];
	
	clusteringNormalizePoint(point, normalizedPoint, min, max);
	
	for(int i=0;i<swarm.repository.getActualSize();i++){
		clusteringNormalizePoint(swarm.repository.getSolution(i), normalizedClusterPoint, min, max);
		
		sum+=calculateEuclideanDistance(normalizedPoint, normalizedClusterPoint, size);
	}
	if(ownCluster)
		return sum/(swarm.repository.getActualSize()-1); //a particle in its own cluster does not contribute for the sum, and shoud not contribute to the average as well
	else
		return sum/swarm.repository.getActualSize();

	
}
//measures the quality of the clustering -- higher values are better
double silhouette(Swarm* multiswarm){
	if(swarmNumber < 2)
		return 1;// if there is only one cluster, the solution is perfecly clustered
		
	double a=0;//dissimilarity of a solution with its own cluster
	double b=MAXDOUBLE;//dissimilarity of a solution with the closest cluster, except the one it belongs to
	double temp=MAXDOUBLE;
	double silGeneral=0;//the silhouette averaged per all swarms
	
	int size=clusteringGetVectorSize();
	double max[size];
	double min[size];
	multiSwarmLimits(multiswarm, min, max);

// 	fprintf(stderr,"\n");
	for(int s=0;s<swarmNumber;s++){
		double sil=0;
		for(int p=0;p<multiswarm[s].repository.getActualSize();p++){
			for(int si=0;si<swarmNumber;si++){
				if(si != s)
					temp=dissimilarity(multiswarm[s].repository.getSolution(p), multiswarm[si], false, min, max);
				else
					a=dissimilarity(multiswarm[s].repository.getSolution(p), multiswarm[si], true, min, max);
				
				b=std::min(temp,b);
			}
			sil+=(b-a)/std::max(a,b);//sum of the silhouettes for each point
			
			b=MAXDOUBLE;
			temp=MAXDOUBLE;
		}
		if(multiswarm[s].repository.getActualSize() > 0)
			silGeneral+=sil/multiswarm[s].repository.getActualSize();//sum of the silhouettes for each swarm
		
// 		fprintf(stderr,"%f ", sil/multiswarm[s].repository.getActualSize());
		// 		printf("\nSilhouette of cluster %d is %f -- %f %f", s, sil, a, b);
	}
// 	fprintf(stderr," (%f) \n", silGeneral/swarmNumber);
	return silGeneral/swarmNumber;
}
//used in Dunn index
double meanPairwiseInnerDistance(Swarm &swarm, double* min, double* max){
	if(swarm.repository.getActualSize() < 2)
		return 0;
	
	int size=clusteringGetVectorSize();
	
	double normalizedP[size];
	double normalizedQ[size];
	double sum=0;

	for(int p=0;p<swarm.repository.getActualSize();p++){
			clusteringNormalizePoint(swarm.repository.getSolution(p), normalizedP, min, max);

		for(int q=p+1;q<swarm.repository.getActualSize();q++){
			clusteringNormalizePoint(swarm.repository.getSolution(q), normalizedQ, min, max);
			
			sum+=calculateEuclideanDistance(normalizedP, normalizedQ, size);
		}
	}
	
	int comparisons=swarm.repository.getActualSize()*(swarm.repository.getActualSize()-1)/2;
	
	return sum/comparisons;
}
//used in Dunn index -- the average of the distances between all elements of both clusters
double interClusterDistance(Swarm &swarm1, Swarm &swarm2, double* min, double* max){
	int size=clusteringGetVectorSize();
	double normalizedP1[size];
	double normalizedP2[size];
	double sum=0;
	
	for(int p1=0;p1<swarm1.repository.getActualSize();p1++){
		clusteringNormalizePoint(swarm1.repository.getSolution(p1), normalizedP1, min, max);

		for(int p2=0;p2<swarm2.repository.getActualSize();p2++){

			clusteringNormalizePoint(swarm2.repository.getSolution(p2), normalizedP2, min, max);
			
			sum+=calculateEuclideanDistance(normalizedP1, normalizedP2, size);
		}
	}
	
	int comparisons=swarm1.repository.getActualSize()*swarm2.repository.getActualSize();

	return sum/comparisons;
}
//measures the quality of the clustering -- higher values are better
double dunnIndex(Swarm* multiswarm){
	if(swarmNumber < 2)
		return 0;
	
	double minInterCluster=MAXDOUBLE;
	double maxInnerCluster=0;
	int size=clusteringGetVectorSize();
	
	double max[size];
	double min[size];
	multiSwarmLimits(multiswarm, min, max);
	
	for(int i=0;i<swarmNumber;i++){
		maxInnerCluster=std::max(maxInnerCluster,meanPairwiseInnerDistance(multiswarm[i], min, max));
		for(int j=i+1;j<swarmNumber;j++){
			minInterCluster=std::min(minInterCluster, interClusterDistance(multiswarm[i], multiswarm[j], min, max));
		}
	}
	
// 	printf("\n %f -- %f -- %f \n", minInterCluster,maxInnerCluster,minInterCluster/maxInnerCluster);
	
	if(maxInnerCluster > 0)
		return minInterCluster/maxInnerCluster;
	else
		return 0;
}

double scatterWithinCluster(Swarm &swarm, double* min, double* max){
	double sum=0;
	int size=clusteringGetVectorSize();
	
	double normalizedCentroid[size];
	double normalizedClusterPoint[size];
	for(int j=0;j<size;j++)
		normalizedCentroid[j]=normalize(swarm.centroid[j], min[j], max[j]);
	
	for(int i=0;i<swarm.repository.getActualSize();i++){
// 		for(int j=0;j<objectiveNumber;j++)
// 			normalizedClusterPoint[j]=normalize(swarm.repository.getSolution(i).objectiveVector[j], min[j], max[j]);
		clusteringNormalizePoint(swarm.repository.getSolution(i), normalizedClusterPoint, min, max);
		
		sum+=calculateEuclideanDistance(normalizedCentroid, normalizedClusterPoint, size);
	}
	
	return sum/swarm.repository.getActualSize();
}
//used in Davies-Bouldin index -- the euclidean distance between two centroids
double interClusterDistance2(Swarm &swarm1, Swarm &swarm2, double* min, double* max){
	int size=clusteringGetVectorSize();
	double normalizedCentroid1[size];
	for(int j=0;j<size;j++)
		normalizedCentroid1[j]=normalize(swarm1.centroid[j], min[j], max[j]);
	
	double normalizedCentroid2[size];
	for(int j=0;j<size;j++)
		normalizedCentroid2[j]=normalize(swarm2.centroid[j], min[j], max[j]);
	
	return calculateEuclideanDistance(normalizedCentroid1, normalizedCentroid2, size);
}
//measures the quality of the clustering -- lower values are better
double daviesBouldinIndex(Swarm* multiswarm){
	int size=clusteringGetVectorSize();
	double max[size];
	double min[size];
	multiSwarmLimits(multiswarm, min, max);
	
	double DB=0;
	for(int i=0;i<swarmNumber;i++){
		double D_i=0;
		for(int j=i+1;j<swarmNumber;j++){
			double R_ij=(scatterWithinCluster(multiswarm[i], min, max)+scatterWithinCluster(multiswarm[j], min, max))/interClusterDistance2(multiswarm[i], multiswarm[j], min, max);
			D_i=std::max(D_i, R_ij);
		}
		DB+=D_i;
	}
	return DB/swarmNumber;
}

//*************************************************************************END OF CLUSTERING QUALILY INDICATORS**************************************************//

//Input: a set of solutions, the size of this set and the K parameter of k-means.
//Output: An array of swarms
void KMeans(Solution* solutions, int solNumber){
	if(solNumber == swarmNumber){//one solution per swarm, nothing to do
		for(int i=0;i<swarmNumber;i++){
			swarms[i].repository.add(solutions[i]);
			// 			memcpy(&multiswarm[i].centroid, &solutions[i].decisionVector, sizeof(double)*decisionNumber);
			for(int d=0;d<decisionNumber;d++)
				swarms[i].centroid[d]=solutions[i].decisionVector[d];
		}
		return;
	}
	
	int maxIter=100; //max iterations to find the right centroids
	int stringlength=clusteringGetVectorSize();
	
	
	// 	double newCentroid[swarmNumber][stringlength], Distances[swarmNumber];
	// 	int wholeClusters[swarmNumber][solNumber], wholeClustersSize[swarmNumber], init_central_pos[swarmNumber];
	bool flag = true;
	// 	printf("\n%d -- %d -- %d\n", swarmNumber, stringlength, solNumber);
	double* Distances = new double[swarmNumber];
	int* wholeClustersSize = new int[swarmNumber];
	int* init_central_pos = new int[swarmNumber];
	double** newCentroid = new double*[swarmNumber];
	int** wholeClusters = new int*[swarmNumber];
	for(int i=0;i<swarmNumber;i++){
		newCentroid[i] = new double[stringlength];
		wholeClusters[i] = new int[solNumber];
	}
	
	//calculates the largest and smaller positions found so far in the search to normalize
	double max[stringlength];
	double min[stringlength];
	for(int i=0;i<stringlength;i++){
		max[i]=MAXDOUBLE*-1;
		min[i]=MAXDOUBLE;
	}
	for(int i=0;i<solNumber;i++){
		for(int d=0;d<stringlength;d++){
			// 			if(!strcmp(clusteringType, "decision")){
			if(stringlength==decisionNumber){
				max[d]=std::max(max[d],solutions[i].decisionVector[d]);
				min[d]=std::min(min[d],solutions[i].decisionVector[d]);
			}
			// 			if(!strcmp(clusteringType, "objectives")){
			if(stringlength==objectiveNumber){
				max[d]=std::max(max[d],solutions[i].objectiveVector[d]);
				min[d]=std::min(min[d],solutions[i].objectiveVector[d]);
			}
			// 			if(!strcmp(clusteringType, "both")){
			if(stringlength==decisionNumber+objectiveNumber){
				if(d<objectiveNumber){
					max[d]=std::max(max[d],solutions[i].objectiveVector[d]);
					min[d]=std::min(min[d],solutions[i].objectiveVector[d]);
				}else{
					max[d]=std::max(max[d],solutions[i].decisionVector[d-objectiveNumber]);
					min[d]=std::min(min[d],solutions[i].decisionVector[d-objectiveNumber]);
				}
			}
		}
	}
	// //end of largest positions calculation	
	
	// Randomly Choose the new central point.
	init_central_pos[0] = rand()%solNumber; //RANDOMNUMBER( solNumber );
	for(int i=1; i<swarmNumber; i++ ) {
		flag = true;
		while(flag) {
			flag = false;
			init_central_pos[i] = rand()%solNumber;//RANDOMNUMBER( solNumber );
			
			for(int j=0; j<i; j++ ) {
				if( init_central_pos[i] == init_central_pos[j] )
					flag = true;
			}
		}
	}
	///////// Periodically Choose ////////////
	for(int i=0; i<swarmNumber; i++ ) {
		clusteringNormalizePoint(solutions[init_central_pos[i]], newCentroid[i], min, max);
	}
	
	/////////////////////////////////////////////////////////////////////////
	////////////////////  Start K-means Clustering  /////////////////////////
	/////////////////////////////////////////////////////////////////////////
	flag = true;
	int Iter = 0;
	while( flag ) {
		
		/*** initialization of wholeClusters & wholeClustersSize ***/
		for(int i=0; i<swarmNumber; i++ ) {
			wholeClustersSize[i] = 0;
			for(int j=0; j<solNumber; j++ )
				wholeClusters[i][j] = 0;
		}
		/**************** End of the Initilization ****************/
		
		for(int i=0; i<solNumber; i++ ) {
			for(int j=0; j<swarmNumber; j++ ) { //distance calculation from each solution to all centroids
				
				double normalizedPoint[stringlength];
				clusteringNormalizePoint(solutions[i], normalizedPoint, min, max);
				
				if(!strcmp(clusteringDistance, "euclidean"))
					Distances[j]=calculateEuclideanDistance(normalizedPoint, newCentroid[j], stringlength);
				if(!strcmp(clusteringDistance, "tchebycheff"))
					Distances[j]=calculateTchebycheffDistance(normalizedPoint, newCentroid[j], stringlength);
				if(!strcmp(clusteringDistance, "minkowski0.5"))
					Distances[j]=calculateMinkowskiDistance(normalizedPoint, newCentroid[j], 0.5, stringlength);
				if(!strcmp(clusteringDistance, "minkowski4"))
					Distances[j]=calculateMinkowskiDistance(normalizedPoint, newCentroid[j], 4, stringlength);
			}
			int cluster_id=0;
			double min = Distances[0];
			for( int x=1; x<swarmNumber; x++ ) {
				if( Distances[x] < min ) {
					min = Distances[x];
					cluster_id = x;
				}
			}
			
			wholeClusters[ cluster_id ][ wholeClustersSize[cluster_id] ] = i; // add sol i to the bottom of wholeClusters[cluster_id] matrix
			wholeClustersSize[cluster_id]++;
		}
		
		// // 		//update the centroids using the mean
		// // 		for(int i=0; i<swarmNumber; i++ ) {
		// // 			for(int j=0; j<stringlength; j++ ) {
		// // 				newCentroid[i][j] = 0;
		// // 				
		// // 				for( int x=0; x<wholeClustersSize[i]; x++ ) {
		// // 					if(!strcmp(clusteringType, "decision"))
		// // 						newCentroid[i][j] += normalize(solutions[wholeClusters[i][x]].decisionVector[j], min[j], max[j]);//solutions[wholeClusters[i][x]][j];// GETREALSELECTED( j, wholeClusters[i][x] );
		// // 					if(!strcmp(clusteringType, "objectives"))
		// // 						newCentroid[i][j] += normalize(solutions[wholeClusters[i][x]].objectiveVector[j], min[j], max[j]);//solutions[wholeClusters[i][x]][j];// GETREALSELECTED( j, wholeClusters[i][x] );
		// // 					if(!strcmp(clusteringType, "both")){
		// // 						if(j<objectiveNumber)
		// // 							newCentroid[i][j] += normalize(solutions[wholeClusters[i][x]].objectiveVector[j], min[j], max[j]);//solutions[wholeClusters[i][x]][j];// GETREALSELECTED( j, wholeClusters[i][x] );
		// // 						else
		// // 							newCentroid[i][j] += normalize(solutions[wholeClusters[i][x]].decisionVector[j-objectiveNumber], min[j], max[j]);//solutions[wholeClusters[i][x]][j];// GETREALSELECTED( j, wholeClusters[i][x] );
		// // 					}
		// // 				}
		// // 				newCentroid[i][j] = newCentroid[i][j]/(double)wholeClustersSize[i];
		// // 			}
		// // 		}//END of updating the centroid using the mean
		
		//update the centroid using the solution that has the smallest distance to the others
		for(int i=0; i<swarmNumber; i++ ) {
			double smallestDistance=MAXDOUBLE;
			for( int j=0; j<wholeClustersSize[i]; j++ ) {
				double normalizedPointJ[stringlength];
				clusteringNormalizePoint(solutions[wholeClusters[i][j]], normalizedPointJ, min, max); //j-th solution of the cluster i
				double distance=0;
				for( int x=0; x<wholeClustersSize[i]; x++ ) {
					double normalizedPointX[stringlength];
					clusteringNormalizePoint(solutions[wholeClusters[i][x]], normalizedPointX, min, max); //x-th solution of the cluster i
					
					if(!strcmp(clusteringDistance, "euclidean"))
						distance+=calculateEuclideanDistance(normalizedPointJ, normalizedPointX, stringlength);
					if(!strcmp(clusteringDistance, "tchebycheff"))
						distance+=calculateTchebycheffDistance(normalizedPointJ, normalizedPointX, stringlength);
					if(!strcmp(clusteringDistance, "minkowski0.5"))
						distance+=calculateMinkowskiDistance(normalizedPointJ, normalizedPointX, 0.5, stringlength);
					if(!strcmp(clusteringDistance, "minkowski4"))
						distance+=calculateMinkowskiDistance(normalizedPointJ, normalizedPointX, 4, stringlength);
				}
				
				if(distance < smallestDistance){
					memcpy(newCentroid[i], normalizedPointJ, sizeof(double)*stringlength);
					smallestDistance=distance;
				}
			}
		}//END of updating the centroid using the solution that has the smallest distance to the others
		
		Iter++;
		if( Iter > maxIter ) flag = false;
	}//END of K-Means Clustering
	
	for(int i=0; i<swarmNumber; i++ ){
		if(swarms[i].repository.getActualSize() != 0){
			fprintf(stderr,"\nERROR! Repository not empty for K-Means clustering\n");
			exit(1);
		}
		for(int c=0;c<wholeClustersSize[i];c++){
			swarms[i].repository.add(solutions[wholeClusters[i][c] ]); //tries to insert the solutions in the repository
		}
// 		swarms[i].repository.organize();
		
		for(int s=0; s<decisionNumber; s++ ){
			if(!strcmp(clusteringType, "decision"))
				swarms[i].centroid[s]= (newCentroid[i][s]*(max[s]-min[s]))+min[s]; //unormalize
				
				if(!strcmp(clusteringType, "objectives")){ //original (used in the experiments)
					double sum=0;
					for(int p=0;p<swarms[i].repository.getActualSize();p++)
						sum+=swarms[i].repository.getSolution(p).decisionVector[s];
					swarms[i].centroid[s]=(sum/swarms[i].repository.getActualSize());
				}
				
				if(!strcmp(clusteringType, "both"))
					swarms[i].centroid[s]= (newCentroid[i][s+objectiveNumber]*(max[s+objectiveNumber]-min[s+objectiveNumber]))+min[s+objectiveNumber]; //unormalize
		}
	}
	
	// // 	for(int i=0; i<swarmNumber; i++ ){
	// // 		for(int s=0; s<stringlength; s++ )
	// // 			printf("%.3f  ", newCentroid[i][s]);
	// // // 		printf("(%d)\n", wholeClustersSize[i]);
	// // 		printf("\n\n\n");
	// // 		for(int c=0;c<multiswarm[i].repository.getActualSize();c++){
	// // 			for(int j=0;j<stringlength;j++)
	// // 				printf("%.3f  ", multiswarm[i].repository.getSolution(c).objectiveVector[j]);  //solutions[wholeClusters[i][c]][j]);
	// // 			printf("\n");
	// // 		}
	// // 		printf("\n\n\n");
	// // 	}
	// // 	//printf("\n------------------------------\n");
	// // 	exit(1);
	
	for(int i=0;i<swarmNumber;i++){
		delete[] newCentroid[i];
		delete[] wholeClusters[i];
	}
	delete[] Distances;
	delete[] wholeClustersSize;
	delete[] init_central_pos;
	delete[] newCentroid;
	delete[] wholeClusters;
	
}

//return the smallest euclidean distance between a vector and a set of solutions
double smallestEuclideanDistance(double* ponto, Solution* solutions, int solNumber) {
	double smallestDistance = MAXDOUBLE;
	//get the smallest distance between a point in the pareto true front and the approximation of the pareto front obtained
	for(int j = 0; j<solNumber;j++){
		smallestDistance = std::min(calculateEuclideanDistance(ponto, solutions[j].objectiveVector, objectiveNumber), smallestDistance);
	}
	return smallestDistance;
}
//return the smallest euclidean distance between a solution and a set of points
double smallestEuclideanDistance(double* sol, double **reference, int tamRef){
	double smallestDistance = MAXDOUBLE;
	//get the smallest distance between a solution in the approximation of the pareto front obtained and a point in the pareto true front
	for(int j = 0; j<tamRef;j++){
		smallestDistance = std::min(calculateEuclideanDistance(sol, reference[j], objectiveNumber), smallestDistance);
	}
	return smallestDistance;
}

double calculateIGDp(Repository &rep){
	char tmp[10], filename[50];
	int fileSize=30000; //30000 should be enough
	double **trueFront = new double*[fileSize];
	for(int i=0; i<fileSize;i++)
		trueFront[i] = new double[objectiveNumber];
	strcpy(tmp, problemName);
	toUpperCase(tmp);
	sprintf(filename, "assessment/metrics/pareto/%s_%d",tmp,objectiveNumber);
	int tamTrue=readFile(trueFront, filename);
	
	double igdp=0;
	double sum = 0;
	//for all points in the true front front
	for(int i = 0; i<tamTrue; i++){
		double smallestDist = smallestEuclideanDistance(trueFront[i], rep.getSolutions(), rep.getActualSize()); 
		sum+=smallestDist*smallestDist;
	}
	igdp = sqrt( (1.0/tamTrue) * sum ); //igdp
	
	for(int i=0; i<fileSize;i++)
		delete[] trueFront[i];
	delete[] trueFront;

	return igdp;
}

//initialize the neighborhood based on the distances in the weight space
void initializeNeighborhood(){
	int neighborhoodSize=globalNeighborhoodSize;
	if(swarmNumber < globalNeighborhoodSize){
		fprintf(stderr,"WARNING! Swarm number smaller than neighborhood size. Neighborhood size set to %d\n", swarmNumber);
		neighborhoodSize=swarmNumber;
	}
// 	int H[21]; H[2]=50;H[3]=12;H[5]=6;H[8]=5;H[10]=5;H[15]=3;H[20]=3;
	int H[21]; H[2]=99;H[3]=12;H[5]=6;H[8]=3;H[10]=3;H[15]=2;H[20]=2;
	
	double h=1;
	if(objectiveNumber <= 5)
		h=H[objectiveNumber]/5.0;
	
	if(objectiveNumber == 2)
		h=H[objectiveNumber]/20.0; //avg of 18.5 with 200 weight vectors
	
// 	double h=1.0;
// 	double h=H[objectiveNumber]/4.0;
	
	double lambda=sqrt(2)*h/(double)H[objectiveNumber];
	
	for(int i=0;i<swarmNumber;i++){
		if(swarms[i].neighborhood != NULL)
			delete[] swarms[i].neighborhood;
		swarms[i].neighborhood = new Neighbor[swarmNumber];
		swarms[i].neighborhoodSize=0;
		for(int j=0;j<swarmNumber;j++){
			swarms[i].neighborhood[j].index=j;
			swarms[i].neighborhood[j].distance=calculateEuclideanDistance(swarms[i].repository.weight, swarms[j].repository.weight, objectiveNumber);
			if(swarms[i].neighborhood[j].distance <= lambda)
				swarms[i].neighborhoodSize++;
			
			if(updateNeighborhoodMetric==7 && i != j)//initialize randomly and keeps this way (except for the first(itself) )
				swarms[i].neighborhood[j].distance=rand();
		}
		std::sort(swarms[i].neighborhood, swarms[i].neighborhood+swarmNumber, neighborsComparator);
		if(neighborhoodSize != 0)//if neighborhood is set statically, overwrite previous setting
			swarms[i].neighborhoodSize=neighborhoodSize;
	}
// 	int ng=1;
// 	for(int j=0;j<swarmNumber;j++){
// 		printf("%f %d -- ", swarms[ng].neighborhood[j].distance, swarms[ng].neighborhood[j].index);
// 		printVector(swarms[j].repository.weight, objectiveNumber);
// 		if(swarms[ng].neighborhood[j].distance <= lambda)
// 			printf(" Y ");
// 		else
// 			printf(" N ");
// 		printf("\n");
// 	}
// 	
// 	int sum=0;
// 	for(int ng=0;ng<swarmNumber;ng++){
// 		printf("lambda: %f optimal size: %d\n", lambda, swarms[ng].neighborhoodSize);
// 		sum+=swarms[ng].neighborhoodSize;
// 	}
// 	printf("avg: %f\n", sum/(double)swarmNumber);
// 	exit(1);
	
}

//average of the log-likelihood given by each model to all the best solutions of each subproblem weighted by the distance between the subproblems, to be equal to the expected matrix
double averageWeightedLogLikelihood(){
	double avgWeightedLogLikelihood=0;
	int evalNumber=0;
	double distMax=swarms[0].neighborhood[swarmNumber-1].distance;
	
	for(int sw=0;sw<swarmNumber;sw++){
// 		if(swarms[sw].repository.getActualSize() != 1){
// 			fprintf(stderr,"ERROR ON WEIGHTED LOG LIKELIHOOD CALCULATION. Repository has %d solutions\n", swarms[sw].repository.getActualSize());
// 			exit(1);
// 		}
		for(int s=0;s<swarmNumber;s++){
// 			double weight=1-(abs(sw-s)/(double)(swarmNumber-1));//first calculation
			double weight=1-(calculateEuclideanDistance(swarms[sw].repository.weight, swarms[s].repository.weight, objectiveNumber)/distMax);
			
// 			printf("%f %f\n", weight, weightDist);
// 			printf("%f ", weight);
			
			for(int so=0;so<swarms[s].repository.getActualSize();so++){
// 				avgWeightedLogLikelihood+=weight*logLikelihood(swarms[s].repository.getSolution(so).decisionVector, swarms[sw]);//sum of the weighted likelihood a model gives to all solutions
				avgWeightedLogLikelihood+=weight*swarms[sw].model.logLikelihood(swarms[s].repository.getSolution(so).decisionVector);//sum of the weighted likelihood a model gives to all solutions
				evalNumber++;
			}
		}
// 		printf("\n");
	}
// 	exit(1);
// 	return avgWeightedLogLikelihood/(swarmNumber*swarmNumber);
	return avgWeightedLogLikelihood/evalNumber;
}

double averageScalarObjFromPopulation(){
	double avgObj=0;
	int evalNumber=0;
	for(int s=0;s<swarmNumber;s++){
		for(int so=0;so<swarms[s].getSize();so++){
			avgObj+=getScalarValue(swarms[s].particles[so].solution.objectiveVector, swarms[s].repository.weight);
			evalNumber++;
		}
		
	}
	return avgObj/evalNumber;
}

double averageScalarObjFromRepository(){
	double avgObj=0;
	int evalNumber=0;
	for(int s=0;s<swarmNumber;s++){
		for(int so=0;so<swarms[s].repository.getActualSize();so++){
			avgObj+=getScalarValue(swarms[s].particles[so].solution.objectiveVector, swarms[s].repository.weight);
			evalNumber++;
		}
	}
	return avgObj/evalNumber;
}

double likelihoodQuality(int which){
	for(int i=0;i<swarmNumber;i++){
		int index=-1;
		double value=MAXDOUBLE*-1;

		if(which == 1){
			//first metric (2.4.1) -- The quality of model Mj for subproblem i is the average likelihood given by model Mj to all solutions in the population weighted by the scalarized fitness of each solution using the weight vector of i (Wi). This can be expressed as:
// 			double scalarSum=0;
			for(int j=0;j<swarmNumber;j++){
				double lik=meanScalarizedFitness(swarms[i].repository.weight, swarms[j]);
				if(lik > value){
					value=lik;
					index=j;
				}
// 				scalarSum+=getScalarValue(swarms[j].repository.getSolution(0).objectiveVector, swarms[i].repository.weight);
			}
// 			printf("%03d %.4f -> %03d %.4f ", i, scalarSum/swarmNumber, index, value);
			return value;
		}
		
		if(which == 2){
			//second metric //just the likelihood
			index=-1;
			value=MAXDOUBLE*-1;
			for(int j=0;j<swarmNumber;j++){
				double lik=swarms[j].model.logLikelihood(swarms[i].repository.getSolution(0).decisionVector);
				if(lik > value){
					value=lik;
					index=j;
				}
			}
// 			printf("%03d %.4f ", index, value);
			return value;
		}
		
		if(which ==3){
			//third metric (2.4.3) -- The quality of model Mj for subproblem i is the correlation between the likelihood given by Mj to all solutions and the scalar value of each solution calculated using the weight vector of i
			index=-1;
			value=MAXDOUBLE;
			for(int j=0;j<swarmNumber;j++){
				double lik=correlationLikelihoodScalarFitness(swarms[i].repository.weight, swarms[j]);
				if(lik < value){
					value=lik;
					index=j;
				}
			}
// 			printf("%03d %.4f ", index, value);
			return value;
		}
		
		if(which ==4){
			//fourth metric (2.4.4) -- The similarity between models Mi and Mj can be measured as the Kullback-Leibler divergence between the probability functions represented by the two models
			index=-1;
			value=MAXDOUBLE;
			for(int j=0;j<swarmNumber;j++){
				double lik=kullbackLeiblerDivergence(swarms[i].model,swarms[j].model);
				if(lik < value && i != j){
					value=lik;
					index=j;
				}
			}
// 			printf("%03d %.4f\n", index, value);
			return value;
		}
	}
}

int localOptimizerObjFunc(int n, int m, double *x, double *f, double *con, void *state_){
	/* Parameter adjustments */
	--con;
	--x;
// 	bool violate=false;
	for(int i=0;i<n;i++){
		con[i+1] = x[i+1]- inferiorPositionLimit[i];
		con[n+i+1] =  superiorPositionLimit[i]- x[i+1];
		
// 		if(x[i] > superiorPositionLimit[i]){
// 			violate=true;
// 			break;
// 		}
// 		if(x[i] < inferiorPositionLimit[i]){
// 			violate=true;
// 			break;
// 		}
	}
	x++;
	
	double obj_f;
// 	if(!violate){
		Problem problem;
		double objectiveVector[objectiveNumber];
		problem.evaluate(x, objectiveVector);
		
		obj_f=getScalarValue(objectiveVector, globalTmpWeight);
// 	}else
// 		obj_f = 4000;

	*f=obj_f ;
	return 0;
}

void mergeNeighboringRepositories(Repository &repOut, Swarm &sw){
	//gathering the solutions to be used to learn and putting them on repOut
	if(decomposition){
		if(!strncmp(solSet, "matrix",6)){
			if(sw.neighborhoodSize != 1){
				fprintf(stderr, "\nERROR! MERGING MATRICES ON NON-INDEPENDENT PROBLEMS (neighborhoodSize != 1)\n");
				exit(1);
			}
			repOut.initialize(archSubSwarms, sw.repository.getActualSize());
			repOut.add(sw.repository);
		}else{
		
			int solsUsed=-1;
			if( (rand()/(double)RAND_MAX) < delta || sw.neighborhoodSize == 1 )//uses only the neighborhood
				solsUsed=sw.neighborhoodSize;
			else
				solsUsed=swarmNumber;
			
			//------------------------------------//
			if(!strcmp(algorithm, "moead")){
				int repSize=0;
				for(int i=0;i<solsUsed;i++){//sum the repository sizes from all neighbors
					int neighbor=sw.neighborhood[i].index;
					repSize+=swarms[neighbor].repository.getActualSize();
				}
				
				repOut.initialize(archSubSwarms, repSize);
				for(int i=0;i<solsUsed;i++){//merge the neighboring solutions
					int neighbor=sw.neighborhood[i].index;
					repOut.add(swarms[neighbor].repository);
				}
				
// 				printf("\n-------------------------------------------------------------------------------\n");
// 				
// 				for(int i=0;i<repOut.getActualSize();i++){
// 					printVector(repOut.getSolution(i).objectiveVector, objectiveNumber);
// 					printf("\n");
// 				}
// 				printf("\n\n\n");
// 				
// 				for(int i=0;i<swarmNumber;i++){
// 					printVector(swarms[i].repository.getSolution(0).objectiveVector, objectiveNumber);
// 					printf("\n");
// 				}
// 				
// 				printf("\n-------------------------------------------------------------------------------\n");
				
			}
			//-------------------------------------//
			
			//---------getting closer to the approach of "injecting cma-es into moead*------"
			
			if(!strcmp(algorithm, "cmaes-mopso")){
		// 		int repSize=solsUsed*10;
// 				int repSize=0;
// 				for(int i=0;i<solsUsed;i++){//sum the repository sizes from all neighbors
// 					int neighbor=sw.neighborhood[i].index;
// 					repSize+=swarms[neighbor].repository.getActualSize();
// 					repSize+=swarms[neighbor].getSize();
// 				}
// 				repSize/=2;
// 	// 			repSize/=10;
				
				repOut.initialize(archSubSwarms, sw.getSize());// mu is \lambda/2, however lambda is set as mu in the paper
				for(int i=0;i<objectiveNumber;i++)
					repOut.weight[i]=sw.repository.weight[i];
				
				//merge the best solution (with respect to my weight) from the population of each subproblem 
				for(int i=0;i<solsUsed;i++){
					int neighbor=sw.neighborhood[i].index;
					int bestIndex=-1;
					double bestQuality=MAXDOUBLE;
					for(int p=0;p<swarms[neighbor].getSize();p++){
						double quality=getScalarValue(swarms[neighbor].particles[p].solution.objectiveVector, sw.repository.weight);
						if(quality < bestQuality || bestQuality == MAXDOUBLE){
							bestQuality=quality;
							bestIndex=p;
						}
					}
					if(bestIndex == -1){
						printf("Error on finding index of the best solution\n");
						exit(1);
					}
					repOut.add(swarms[neighbor].particles[bestIndex].solution);
				}
				
// 				for(int i=0;i<solsUsed;i++){//merge the neighboring solutions from the population
// 					int neighbor=sw.neighborhood[i].index;
// 					for(int j=0;j<swarms[neighbor].getSize();j++)
// 						repOut.add(swarms[neighbor].particles[j].solution);
// 				}
			}
			
			
			//---------END of getting closer to the approach of "injecting cma-es into moead*------"

			
			
			
			
			
			
			
			
	// 		repOut.initialize(archSubSwarms, 25);
	// 		for(int i=0;i<objectiveNumber;i++)
	// 			repOut.weight[i]=sw.repository.weight[i];
	// 		
	// 		for(int i=0;i<swarmNumber;i++)
	// 			for(int j=0;j<swarms[i].repository.getActualSize();j++)
	// 				repOut.add(swarms[i].repository.getSolutions()[j]);
	// // 		repOut.add(sw.repository);
	// 		
	// 		for(int i=0;i<swarmNumber;i++){
	// 			for(int j=0;j<swarms[i].getSize();j++)
	// 				repOut.add(swarms[i].particles[j].solution);
	// 		}
			
			
			// 		for(int i=0;i<sw.getSize();i++)
			// 			repOut.add(sw.particles[i].solution);
			
			// // 	my test, add the solution generated at the iteration
			// 		repOut.initialize(archSubSwarms, sw.neighborhoodSize*sw.getSize()+1);
			// 		repOut.initialize(archSubSwarms, sw.getSize()+2);
			// 		for(int i=0;i<objectiveNumber;i++)
			// 			repOut.weight[i]=sw.repository.weight[i];
			
			// 		for(int n=0;n<sw.neighborhoodSize;n++){
			// 			int neighbor=sw.neighborhood[n].index;
			// 			for(int i=0;i<sw.getSize();i++){
			// 				repOut.add(swarms[neighbor].particles[i].solution);
			// 			}
			// 		}
			
			
	// 		for(int i=0;i<sw.getSize();i++){
	// 			repOut.add(sw.particles[i].solution);
	// 		}
			
			
			// 		int idx=(rand() % sw.neighborhoodSize);
			// 		int neighbor=sw.neighborhood[idx].index;
			// 		repOut.add(swarms[neighbor].repository);
			
			// 		if( rand()/(double)RAND_MAX < 0.01 )
			// 			repOut.add(sw.repository);
			
			// 		for(int i=0;i<repOut.getActualSize();i++){//merge the neighboring solutions
			// 			printVector(repOut.getSolution(i).objectiveVector, objectiveNumber);
			// 			printf("\n");
			// 		}
			// 		printf("%d -- %d == %d\n", repOut.getActualSize(), sw.neighborhoodSize, repSize);
		}
	}else{
		fprintf(stderr, "\nERROR! MERGING NEIGHBORHOOD REPOSITORIES ON NON-DECOMPOSITION SITUATION\n");
		exit(1);
	}
// 	repOut.organize();
}

double correlation(double *x, double *y, const int &size){
	//https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Definition
	double meanX=0, meanY=0;
	for(int i=0;i<size;i++){
		meanX+=x[i];
		meanY+=y[i];
	}
	meanX/=size;
	meanY/=size;
	
	double p1=0,p2=0,p3=0;
	for(int i=0;i<size;i++){
		p1+=x[i]*y[i];
		p2+=(x[i]*x[i]);
		p3+=(y[i]*y[i]);
	}
	p1-=size*meanX*meanY;
	p2-=size*meanX*meanX;
	p3-=size*meanY*meanY;
	
	if(p1!=0 && p2!=0 && p3!=0)
		return p1/( sqrt(p2)*sqrt(p3) );
	else
		return 0;
}

//Quality of a model in relation to a given subproblem i (weight), calculated as the mean weighted (by scalar function i) likelihood given by the model to the best solution of all subproblems
double meanScalarizedFitness(double* weight, Swarm &sw){
	int counter=0;
	double sum=0;
	for(int i=0;i<swarmNumber;i++){
		for(int j=0;j<swarms[i].repository.getActualSize();j++){
			const Solution *sol = &swarms[i].repository.getSolution(j);
			double scalar=getScalarValue(sol->objectiveVector, weight);
// 			double lik=logLikelihood(sol->decisionVector, sw);
			double lik=sw.model.logLikelihood(sol->decisionVector);
			sum+=scalar*lik;
			counter++;
		}
	}
	return sum/counter;
}
//Correlation between the likelihood of all the population and the scalar fitness according to the given weight
double correlationLikelihoodScalarFitness(double* weight, Swarm &sw){
	double scalarFitness[swarmNumber];//Average scalar fitness of all the best solutions for each swarm (usually one)
	double likelihood[swarmNumber];//Average log-likelihood of all the best solutions for each swarm (usually one)
	
	for(int i=0;i<swarmNumber;i++){
		double sumFit=0;
		double sumLik=0;
		for(int j=0;j<swarms[i].repository.getActualSize();j++){
			const Solution *sol = &swarms[i].repository.getSolution(j);
			sumFit+=getScalarValue(sol->objectiveVector, weight);
// 			sumLik+=logLikelihood(sol->decisionVector, sw);
			sumLik+=sw.model.logLikelihood(sol->decisionVector);
		}
		scalarFitness[i]=sumFit/swarms[i].repository.getActualSize();
		likelihood[i]=sumLik/swarms[i].repository.getActualSize();
	}
	return correlation(scalarFitness, likelihood, swarmNumber);
}
//returns the average likelihood the neighborhood gives to a solution, which reflects the "opinion" of the neighborhood about the quality of the solution
double neighborhoodOpinion(Swarm &swarm, double* sol){
	double opinion=0;
	int neighborsConsidered=0;
// 	for(int n=0;n<swarm.neighborhoodSize;n++){//consider all
	for(int n=0;n<1;n++){//consider just myself
		int neighbor=swarm.neighborhood[n].index;
		if(!swarms[neighbor].init){//if it is valid
// 			opinion+=logLikelihood(sol, swarms[neighbor]);
			opinion+=swarms[neighbor].model.logLikelihood(sol);
			neighborsConsidered++;
		}
	}
	return opinion/(double)neighborsConsidered;
}

void updateNeighborhood(){
	if(swarmNumber < 2 && (updateNeighborhoodMetric >= 1 && updateNeighborhoodMetric <= 5) ){
		fprintf(stderr,"ERROR ON UPDATING NEIGHBORHOOD! Less than 2 swarms\n");
		exit(1);
	}
	
	//first metric //weighted likelihood -- smaller is better (because instead of the direct likelihood, it is considered the difference between the current and the maximum)
	if(updateNeighborhoodMetric==1){
		double likelihoods[swarmNumber][swarmNumber], scalarizedFitness[swarmNumber][swarmNumber];
		double maxLik=MAXDOUBLE*-1;
		//first step: calculation of the matrices for speed up the code
		for(int i=0;i<swarmNumber;i++){//swarms
			for(int j=0;j<swarmNumber;j++){//swarms
				likelihoods[i][j]=0;
				scalarizedFitness[i][j]=0;
				for(int k=0;k<swarms[j].repository.getActualSize();k++){//repository of swarm i
					const Solution *sol = &swarms[j].repository.getSolution(k);
// 					likelihoods[i][j]+=logLikelihood(sol->decisionVector, swarms[i]);
					likelihoods[i][j]+=swarms[i].model.logLikelihood(sol->decisionVector);
					scalarizedFitness[i][j]+=getScalarValue(sol->objectiveVector, swarms[i].repository.weight);
				}
				likelihoods[i][j]/=swarms[j].repository.getActualSize();
				maxLik=max(maxLik,likelihoods[i][j]);
				scalarizedFitness[i][j]/=swarms[j].repository.getActualSize();
			}
		}
		
		
		//second step, update the values
		for(int i=0;i<swarmNumber;i++){//swarms
			for(int j=0;j<swarmNumber;j++){//swarms
				swarms[i].neighborhood[j].index=j;
				if(i != j){
					swarms[i].neighborhood[j].distance=0;
					for(int k=0;k<swarmNumber;k++)
// 						swarms[i].neighborhood[j].distance+=scalarizedFitness[i][k]*likelihoods[j][k]*-1;//multiplies by -1 so now smaller is better, like the distance
						swarms[i].neighborhood[j].distance+=scalarizedFitness[i][k]*(maxLik-likelihoods[j][k]);//since the likelihood now is the difference between the current and the maximum, the two indicators are the smaller the better, like the distance
					swarms[i].neighborhood[j].distance/=swarmNumber;
// 					printf("%f == %f\n", (swarms[i].neighborhood[j].distance), (meanScalarizedFitness(swarms[i].repository.weight, swarms[j])*-1));//comparison between the optimized (faster) value and the previous
				}
				else
					swarms[i].neighborhood[j].distance=MAXDOUBLE*-1;//make sure the first neighbor is myself
			}
		}
	}

	//second metric //just the likelihood -- higher is better
	if(updateNeighborhoodMetric==2){
		for(int i=0;i<swarmNumber;i++){
			for(int j=0;j<swarmNumber;j++){
				swarms[i].neighborhood[j].index=j;
				if(i != j)
// 					swarms[i].neighborhood[j].distance=logLikelihood(swarms[i].repository.getSolution(0).decisionVector, swarms[j])*-1;//multiplies by -1 so now smaller is better, like the distance
					swarms[i].neighborhood[j].distance=swarms[j].model.logLikelihood(swarms[i].repository.getSolution(0).decisionVector)*-1;//multiplies by -1 so now smaller is better, like the distance
				else
					swarms[i].neighborhood[j].distance=MAXDOUBLE*-1;//make sure the first neighbor is myself
			}
		}
	}
	
	//third metric -- smaller is better (negative correlation)
	if(updateNeighborhoodMetric==3){
		double likelihoods[swarmNumber][swarmNumber], scalarizedFitness[swarmNumber][swarmNumber];
		//first step: calculation of the matrices for speed up the code
		for(int i=0;i<swarmNumber;i++){//swarms
			for(int j=0;j<swarmNumber;j++){//swarms
				likelihoods[i][j]=0;
				scalarizedFitness[i][j]=0;
				for(int k=0;k<swarms[j].repository.getActualSize();k++){//repository of swarm i
					const Solution *sol = &swarms[j].repository.getSolution(k);
					likelihoods[i][j]+=swarms[i].model.logLikelihood(sol->decisionVector);
					scalarizedFitness[i][j]+=getScalarValue(sol->objectiveVector, swarms[i].repository.weight);
				}
				likelihoods[i][j]/=swarms[j].repository.getActualSize();
				scalarizedFitness[i][j]/=swarms[j].repository.getActualSize();
			}
		}

		//second step, update the values
		for(int i=0;i<swarmNumber;i++){//swarms
			for(int j=0;j<swarmNumber;j++){//swarms
				swarms[i].neighborhood[j].index=j;
				if(i != j){
					swarms[i].neighborhood[j].distance=correlation(scalarizedFitness[i], likelihoods[j], swarmNumber);
// 					printf("%f == %f\n", (swarms[i].neighborhood[j].distance), (correlationLikelihoodScalarFitness(swarms[i].repository.weight, swarms[j])) );//comparison between the optimized (faster) value and the previous
				}
				else
					swarms[i].neighborhood[j].distance=MAXDOUBLE*-1;//make sure the first neighbor is myself
			}
		}
	}
	//fourth metric -- smaller is better
	if(updateNeighborhoodMetric==4){
		for(int i=0;i<swarmNumber;i++){
			for(int j=0;j<swarmNumber;j++){
				swarms[i].neighborhood[j].index=j;
				if(i != j)
					swarms[i].neighborhood[j].distance=kullbackLeiblerDivergence(swarms[i].model,swarms[j].model);
				else
					swarms[i].neighborhood[j].distance=MAXDOUBLE*-1;//make sure the first neighbor is myself
			}
		}
	}
	//fifth metric -- smaller is better
	if(updateNeighborhoodMetric==5){//just the distance in the decision space between the solutions
		for(int i=0;i<swarmNumber;i++){
			for(int j=0;j<swarmNumber;j++){
				swarms[i].neighborhood[j].index=j;
				if(i != j)
					swarms[i].neighborhood[j].distance=calculateEuclideanDistance(swarms[i].repository.getSolution(0).decisionVector,swarms[j].repository.getSolution(0).decisionVector, decisionNumber);
				else
					swarms[i].neighborhood[j].distance=MAXDOUBLE*-1;//make sure the first neighbor is myself
			}
		}
	}
	//sixth metric -- random
	if(updateNeighborhoodMetric==6){//random
		for(int i=0;i<swarmNumber;i++){
			for(int j=0;j<swarmNumber;j++){
				swarms[i].neighborhood[j].index=j;
				if(i != j)
					swarms[i].neighborhood[j].distance=rand();
				else
					swarms[i].neighborhood[j].distance=MAXDOUBLE*-1;//make sure the first neighbor is myself
			}
		}
	}
	if(updateNeighborhoodMetric >=1 && updateNeighborhoodMetric <=6){//if it is valid, reorder the neighbors
		for(int i=0;i<swarmNumber;i++)//just sort everything
			std::sort(swarms[i].neighborhood, swarms[i].neighborhood+swarmNumber, neighborsComparator);
// 		printf("OK, strategy %d\n",updateNeighborhoodMetric);
	}
	
// 		printf("%d\n",i);
	
// 	for(int i=0;i<swarmNumber;i++){
// 		printf("Neighbors of %d: ", i);
// 		for(int j=0;j<swarms[i].neighborhoodSize;j++){
// 			printf("%d ",swarms[i].neighborhood[j].index);
// 		}
// 		printf("\n");
// 	}
// 	exit(1);
}

//returns the likelihood rank of a given solution in relation to a given (ORDERED) repository, where the likelihoods are calculated in relation to a model contained in a given swarm
double likelihoodRanking(const Solution &sol, Repository &rep, const CMAESModel model){
	double likelihood=model.logLikelihood(sol.decisionVector);
	
// 	for(int i=0;i<rep.getActualSize();i++){
// 		rep.getSolutions()[i].crowdingDistance=logLikelihood(rep.getSolution(i).decisionVector, sw);//stores temporarily in the crowding distance field
// 	}
// 	std::sort(rep.getSolutions(), rep.getSolutions()+rep.getActualSize(), crowdingComparatorSol);
	
	if(likelihood > rep.getSolution(0).crowdingDistance){
// 		printf("BEST: lik: %f, sol0: %f\n", likelihood, rep.getSolution(0).crowdingDistance);
		return 0;
	}
	for(int i=1;i<rep.getActualSize();i++){
		if(rep.getSolution(i).crowdingDistance > rep.getSolution(i-1).crowdingDistance){
			fprintf(stderr,"ERROR ON CALCULATING LIKELIHOOD RANKING! Input repository is not ordered.\n");
			exit(1);
		}
		if(likelihood <= rep.getSolution(i-1).crowdingDistance && likelihood >= rep.getSolution(i).crowdingDistance){
// 			printf(" OK i: %d lik before: %f lik: %f lik after: %f\n", i, rep.getSolution(i-1).crowdingDistance,likelihood,rep.getSolution(i).crowdingDistance);
			if(likelihood == rep.getSolution(i-1).crowdingDistance)
				return i-1;
			if(likelihood == rep.getSolution(i).crowdingDistance)
				return i;
			return (i+(i-1.0))/2.0;
		}
	}
	return rep.getActualSize();
}
//update the quality of the models as the correlation between the likelihood and the scalar fitness regarding all solutions in the global repository
void updateModelsQuality(){
	if(repGlobal->getActualSize() > 10){//avoid too small repositories
// 		repGlobal->organize();
		
		double likelihood[repGlobal->getActualSize()], scalarFitness[repGlobal->getActualSize()];
		for(int i=0;i<swarmNumber;i++){
			Swarm* sw=&swarms[i];
			for(int j=0;j<repGlobal->getActualSize();j++){
				likelihood[j]=sw->model.logLikelihood(repGlobal->getSolution(j).decisionVector);
				scalarFitness[j]=getScalarValue(repGlobal->getSolution(j).objectiveVector, sw->repository.weight);
			}
			sw->model.modelQuality=correlation(likelihood, scalarFitness, repGlobal->getActualSize())*-1;//inverse of the correlation because the higher the likelihood the smaller the scalar fitness value
// 			printf("quality (more than 10)%f sz: %d\n", sw->modelQuality*-1, repGlobal->getActualSize());
// 			printf("%.2f ", sw->modelQuality);
		}
	}else{
		double likelihood[swarmNumber], scalarFitness[swarmNumber];
		for(int i=0;i<swarmNumber;i++){
			Swarm* sw=&swarms[i];
			for(int j=0;j<swarmNumber;j++){
				likelihood[j]=sw->model.logLikelihood(swarms[j].repository.getSolution(0).decisionVector);
				scalarFitness[j]=getScalarValue(swarms[j].repository.getSolution(0).objectiveVector, sw->repository.weight);
			}
			sw->model.modelQuality=correlation(likelihood, scalarFitness, swarmNumber)*-1;//inverse of the correlation because the higher the likelihood the smaller the scalar fitness value
// 			printf("quality (less than 10)%f sz: %d\n", sw->modelQuality*-1, repGlobal->getActualSize());
// 			printf("%.2f ", sw->modelQuality);
		}
	}
// 	printf("\n");
}

int calculateDominationRanks(){//used in moead-dd
	int nonCalculated=0;
	//counting and initializing
	for(int s=0;s<swarmNumber;s++){//for all the swarms
		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
			swarms[s].repository.getSolutions()[p].dominanceLevel=-1;
			nonCalculated++;
		}
	}
	
	int currentDomLevel=0;
	while(nonCalculated>0){
		for(int s1=0;s1<swarmNumber;s1++){//for all the swarms
			for(int p1=0;p1<swarms[s1].repository.getActualSize();p1++){//for all the particles in a given swarm repository
				Solution* sol1=&swarms[s1].repository.getSolutions()[p1];
				if(sol1->dominanceLevel == currentDomLevel || sol1->dominanceLevel==-1){//only compares dominance if the solution is not ranked or is in the current level (in case it was disconsidered previously)
					bool dominated=false;//the solution p1 from s1 is nonDominated
					for(int s2=0;s2<swarmNumber;s2++){//for all the swarms
						for(int p2=0;p2<swarms[s2].repository.getActualSize();p2++){//for all the particles in a given swarm repository
							Solution* sol2=&swarms[s2].repository.getSolutions()[p2];
							if(sol2->dominanceLevel==currentDomLevel || sol2->dominanceLevel==-1){
								int dom=dominance(sol1->objectiveVector, sol2->objectiveVector, objectiveNumber); //1 sol1 dominates, -1 sol2 dominates
								if(dom==-1){//if sol 2 dominates sol1
// // 									printf("sol: ");
// // 									printVector(sol1->objectiveVector, objectiveNumber);
// // 									printf(" dominates sol: ");
// // 									printVector(sol2->objectiveVector, objectiveNumber);
// // 									printf("\n");
									
									dominated=true;
									break;
								}else{
									if(dom==1)//if sol2 is dominated, it is not considered in the current level anymore (optimization)
										sol2->dominanceLevel=currentDomLevel+1;
								}
							}
						}
						if(dominated)
							break;
					}
					if(!dominated){
						nonCalculated--;
						sol1->dominanceLevel=currentDomLevel;
					}else//if sol1 is dominated, it is not considered in the current level anymore (optimiZation)
						sol1->dominanceLevel=currentDomLevel+1;
				}
			}
		}
// 		printf("non-calculated %d, domLevel %d\n",nonCalculated, currentDomLevel);
		currentDomLevel++;
	}
	return currentDomLevel-1;//return the last level index
}

//returns the right domination level of the given solution regarding all the swarms
int calculateDominationRank(const Solution &sol){
	int lastLevel=-1;
	bool isHere=false;//mark if the solution being added is in the global front
	for(int s=0;s<swarmNumber;s++){//for all the swarms
		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
			if(swarms[s].repository.getSolution(p).dominanceLevel > lastLevel)
				lastLevel=swarms[s].repository.getSolution(p).dominanceLevel;
			if(isEqual(swarms[s].repository.getSolution(p).objectiveVector, sol.objectiveVector, objectiveNumber))
				isHere=true;
		}
	}
	if(!isHere){
		return -2;
	}
	
	bool dominated=false;
	for(int currentDomLevel=0;currentDomLevel<=lastLevel;currentDomLevel++){//for all the dominance levels
		dominated=false;//the given solution is nonDominated
		for(int s=0;s<swarmNumber;s++){//for all the swarms
			for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
				Solution* sol2=&swarms[s].repository.getSolutions()[p];
				if(sol2->dominanceLevel==currentDomLevel){//only compares with solutions from this dominance level
					int dom=dominance(sol.objectiveVector, sol2->objectiveVector, objectiveNumber); //1 sol dominates, -1 sol2 dominates
					if(dom==-1){//if sol 2 dominates sol						
						dominated=true;
						break;
					}
				}
			}
			if(dominated)
				break;
		}
		if(!dominated){
			return currentDomLevel;
		}
	}
// 	if(sol.dominanceLevel==-1){//if the solution was not classified so far, it is the worst
// 		lastLevel;
// 	}
	return lastLevel+1;
}

int updateDominationRanksAdd(Solution &sol){
	sol.dominanceLevel=-1;
	int lastLevel=-1;
	bool isHere=false;//mark if the solution being added is in the global front
	for(int s=0;s<swarmNumber;s++){//for all the swarms
		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
			if(swarms[s].repository.getSolution(p).dominanceLevel > lastLevel)
				lastLevel=swarms[s].repository.getSolution(p).dominanceLevel;
			if(isEqual(swarms[s].repository.getSolution(p).objectiveVector, sol.objectiveVector, objectiveNumber))
				isHere=true;
		}
	}
	if(!isHere)
// 		return lastLevel;
		exit(1);
	
	bool dominated;
	for(int currentDomLevel=0;currentDomLevel<=lastLevel+1;currentDomLevel++){//for all the dominance levels
		dominated=false;//the given solution is nonDominated
		for(int s=0;s<swarmNumber;s++){//for all the swarms
			for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
				Solution* sol2=&swarms[s].repository.getSolutions()[p];
				if(sol2->dominanceLevel==currentDomLevel || sol2->dominanceLevel==-1){//only compares with solutions from this dominance level or with solutions being updated
					
					int dom=dominance(sol.objectiveVector, sol2->objectiveVector, objectiveNumber); //1 sol dominates, -1 sol2 dominates
					if(dom==-1){//if sol 2 dominates sol
						dominated=true;
						break;
					}else{
						if(dom==1){//if sol2 is dominated, it is not considered in the current level anymore (optimization)
							sol.dominanceLevel=currentDomLevel;//stays here temporarily
							updateDominationRanksAdd(*sol2);//recursively update the position of the dominated solution
						}
					}
				}
			}
			if(dominated)
				break;
		}
		if(!dominated){
			sol.dominanceLevel=currentDomLevel;
// // 			break;
			return lastLevel;
		}
	}
	printf("It is not supposed to be here (%d)\n", lastLevel);
	return lastLevel;
}

void updateDominationRanksRemove(Solution &sol){
	//check if current solution exists, otherwise it was removed
	bool isHere=false;//mark if the solution being added is in the global front
	for(int s=0;s<swarmNumber;s++){//for all the swarms
		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
			if(swarms[s].repository.getSolution(p).dominanceLevel == sol.dominanceLevel)
				if(isEqual(swarms[s].repository.getSolution(p).objectiveVector, sol.objectiveVector, objectiveNumber)){
					isHere=true;
					break;
				}
		}
		if(isHere)
			break;
	}

	int before=sol.dominanceLevel;
	sol.dominanceLevel=calculateDominationRank(sol);

	if(!isHere || before > sol.dominanceLevel){//if the solution is not present or if the rank improved
// 		printf("bef: %d after: %d \n", before, sol.dominanceLevel);
			
		for(int s=0;s<swarmNumber;s++){//for all the swarms
			for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
				Solution* sol2=&swarms[s].repository.getSolutions()[p];
				if(sol2->dominanceLevel==before+1){//only compares with solutions dominated previously
					int dom=dominance(sol.objectiveVector, sol2->objectiveVector, objectiveNumber); //1 sol dominates, -1 sol2 dominates
					if(dom==-1){//if sol 2 dominates sol and the solution is present
						printf("This is very weird, the solution from dominance level l+1 dominates the solution from l\n");
						exit(1);
					}else{
						if(dom==1){//if sol2 is dominated
							updateDominationRanksRemove(*sol2);
						}
					}
				}
			}
		}
	}
	if(before < sol.dominanceLevel)
		printf("It is now worse\n");
}


//auxiliary structure
struct indexNicheCount{
	int index;
	int nicheCount;
	double scalarSum;//tie break criterion of part 1 of algorithm 5
};
// order by niche count, the largest goes first
bool nicheCountComparator(const indexNicheCount &i, const indexNicheCount &j) { return (i.nicheCount>j.nicheCount); }
// order by scalarSum, the largest goes first
bool scalarSumComparator(const indexNicheCount &i, const indexNicheCount &j) { return (i.scalarSum>j.scalarSum); }
// order solutions by domination level, the largest goes first
bool dominationLevelComparator(const Solution &i, const Solution &j) { return (i.dominanceLevel>j.dominanceLevel); }
// order by weighted distances, the largest goes first
bool largestWeightedDistanceComparator(const Solution &i,const Solution &j) { return (i.weightedDistance>j.weightedDistance); }

double moead_ddScalarSum(int index){
	double result=0;
	for(int p=0;p<swarms[index].repository.getActualSize();p++){
		result+=getScalarValue(swarms[index].repository.getSolution(p).objectiveVector,swarms[index].repository.weight);
	}
	return result;
}

void moead_ddRemoveWorst(){
	//as in locate_worst from algorithm 5 of http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6964796&tag=1
	indexNicheCount ind[swarmNumber];
	
	//identify the most crowded sub-region (sub-problem)
	for(int s=0;s<swarmNumber;s++){
		ind[s].index=s;
		ind[s].nicheCount=swarms[s].repository.getActualSize();
	}
	std::sort(ind, ind+swarmNumber, nicheCountComparator);
	
	if(ind[0].nicheCount == ind[1].nicheCount){//if there is at least one tie in niche count, sort again by scalar value sum
		ind[0].scalarSum=moead_ddScalarSum(0);
		int count=1;
		for(int s=1;s<swarmNumber;s++){
			if(ind[s].nicheCount == ind[s-1].nicheCount){
				ind[s].scalarSum=moead_ddScalarSum(s);
				count++;
			}
			else
				break;
		}
		std::sort(ind, ind+count, scalarSumComparator);
	}
	Repository* worstRegion=&swarms[ind[0].index].repository;//stores the index of the worst subproblem identified
	//end of identifying the most crowded sub-region
	
	std::sort(worstRegion->getSolutions(), worstRegion->getSolutions()+worstRegion->getActualSize(), dominationLevelComparator);
	if(worstRegion->getActualSize() > 1 &&  worstRegion->getSolution(0).dominanceLevel == worstRegion->getSolution(1).dominanceLevel){//if there is at least one tie in dominance level, sort again by scalar value
		worstRegion->getSolutions()[0].weightedDistance=getScalarValue(worstRegion->getSolution(0).objectiveVector,worstRegion->weight);
		int count=1;
		for(int s=1;s<worstRegion->getActualSize();s++){
			if(worstRegion->getSolution(s).dominanceLevel == worstRegion->getSolution(s-1).dominanceLevel){
				worstRegion->getSolutions()[s].weightedDistance=getScalarValue(worstRegion->getSolution(s).objectiveVector,worstRegion->weight);
				count++;
			}
			else
				break;
		}
		std::sort(worstRegion->getSolutions(), worstRegion->getSolutions()+count, largestWeightedDistanceComparator);
	}

	Solution sol=worstRegion->getSolution(0);
	worstRegion->exclude(0);//since the solutions were reordered, exclude the first, that is the worst
	updateDominationRanksRemove(sol);
// 	worstRegion->organize();
	if(worstRegion->getActualSize() == 0){
		printf("problem on worst\n");
// 		for(int s=0;s<swarmNumber;s++){
// 			printf("index %d has niche count of %d and scalar sum of %f\n", ind[s].index, ind[s].nicheCount, ind[s].scalarSum);
// 		}
// 		printf("chosen %d\n\n", ind[0].index);
	}
	
}

void moead_ddInsert(Solution solution){
	//associate the solution with the best sub-region (weight vector)
	double minValue=MAXDOUBLE;
	int index=-1;

	for(int j=0;j<swarmNumber;j++){
		double dist=getScalarValue(solution.objectiveVector, swarms[j].repository.weight);
		if(dist > MAXDOUBLE)//if the distance is too big, disregard the solution
			return;
		//after calculating the distance regarding this region, test if this is the better
		if(dist <= minValue){
			minValue=dist;
			index=j;
		}
	}
	if(index==-1){
		fprintf(stderr,"\nMOEA/D-DD INSERT ERROR %f\n", minValue);
		exit(1);
	}
	
// 	int numbers1[swarmNumber];
// 	int numbers2[swarmNumber];
	int lastLevel=-1;
	
	
// // 	//**** comparison with previous scheme ***//
// // 	for(int s=0;s<swarmNumber;s++)//for all the swarms
// // 		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
// // // 			printf("2nd dom: %d\n", swarms[s].repository.getSolution(p).dominanceLevel);
// // 			numbers1[s]=swarms[s].repository.getSolution(p).dominanceLevel;
// // 		}
// // 			
// // 	
// // 	lastLevel=calculateDominationRanks();
// // 	for(int s=0;s<swarmNumber;s++){//for all the swarms
// // 		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
// // // 			printf("2nd dom: %d\n", swarms[s].repository.getSolution(p).dominanceLevel);
// // 			numbers2[s]=swarms[s].repository.getSolution(p).dominanceLevel;
// // 		}
// // 		if(numbers1[s] != numbers2[s]){
// // 			printf("ERROR!! part1 update: %d != scratch: %d calculate: %d\n", numbers1[s], numbers2[s], calculateDominationRank(solution));
// // 			printf("\n\n");
// // 		}
// // 	}
// // 		//**** comparison with previous scheme ***//
	solution.dominanceLevel=-2;
	bool add=swarms[index].repository.forceAdd(solution);//associated
	
	for(int s=0;s<swarmNumber;s++)//for all the swarms
		for(int p=0;p<swarms[s].repository.getActualSize();p++)//for all the particles in a given swarm repository
			if(swarms[s].repository.getSolution(p).dominanceLevel==-2){
				lastLevel=updateDominationRanksAdd(swarms[s].repository.getSolutions()[p]);
			}
		
		
// 	//**** comparison with previous scheme ***//
// 	for(int s=0;s<swarmNumber;s++)//for all the swarms
// 		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
// 			numbers1[s]=swarms[s].repository.getSolution(p).dominanceLevel;
// 		}
// 			
// 	lastLevel=calculateDominationRanks();
// 	for(int s=0;s<swarmNumber;s++){//for all the swarms
// 		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
// // 			printf("2nd dom: %d\n", swarms[s].repository.getSolution(p).dominanceLevel);
// 			numbers2[s]=swarms[s].repository.getSolution(p).dominanceLevel;
// 		}
// 		if(numbers1[s] != numbers2[s]){
// 			printf("ERROR!! part2 update: %d != scratch: %d calculate: %d\n", numbers1[s], numbers2[s], calculateDominationRank(solution));
// 			printf("\n\n");
// 		}
// 	}
// 	//**** comparison with previous scheme ***//

	
	if(add){//only removes the worst solution if the current solution was successfully added - in some cases the solution is not added because it is equal to the previous
		if(lastLevel==0){//all solutions are non-dominated
			moead_ddRemoveWorst();
		}else{
			//count the number of solutions in the last level
			int lastLevelCount=0;
			int lastLevelRegionIndex=-1;
			int lastLevelSolutionIndex=-1;
			for(int s=0;s<swarmNumber;s++){//for all the swarms
				Repository* rep;
				rep=&swarms[s].repository;
				for(int p=0;p<rep->getActualSize();p++){//for all the particles in a given swarm repository
					if(rep->getSolution(p).dominanceLevel==lastLevel){
						lastLevelCount++;
						lastLevelRegionIndex=s;//stores the index of the subswarm with the last solution from the last dominance level, in case this is the only one
						lastLevelSolutionIndex=p;//stores the index of the solution from the last domination level, in case this is the only one
					}
				}
			}//end of counting
			if(lastLevelCount==0){
				printf("Error! No solutions in last non-domination level\n");
				exit(1);
			}
			
			if(lastLevelCount==1){//if Fl has only one solution
				if(swarms[lastLevelRegionIndex].repository.getActualSize()>1){//if the subregion has more than one solution
					Solution sol=swarms[lastLevelRegionIndex].repository.getSolution(lastLevelSolutionIndex);
					swarms[lastLevelRegionIndex].repository.exclude(lastLevelSolutionIndex);//remove the worst solution
					updateDominationRanksRemove(sol);
				}else{//if the bad solution is the only one of the subregion
					moead_ddRemoveWorst();
				}
			}else{//if Fl has more than one solution
				//identify the most crowded sub-region (sub-problem) that has at least one solution in Fl
				indexNicheCount ind[swarmNumber];
				
				for(int s=0;s<swarmNumber;s++){
					ind[s].nicheCount=0;
					ind[s].index=s;
					for(int p=0;p<swarms[s].repository.getActualSize();p++){
						if(swarms[s].repository.getSolution(p).dominanceLevel==lastLevel){//check if at least one solution is in Fl
							ind[s].nicheCount=swarms[s].repository.getActualSize();
							break;
						}
					}
				}
				std::sort(ind, ind+swarmNumber, nicheCountComparator);
				
				if(ind[0].nicheCount == ind[1].nicheCount){//if there is at least one tie in niche count, sort again by scalar value sum
					ind[0].scalarSum=moead_ddScalarSum(0);
					int count=1;
					for(int s=1;s<swarmNumber;s++){
						if(ind[s].nicheCount == ind[s-1].nicheCount){
							ind[s].scalarSum=moead_ddScalarSum(s);
							count++;
						}
						else
							break;
					}
					std::sort(ind, ind+count, scalarSumComparator);
				}
				
				Repository* worstRegion=&swarms[ind[0].index].repository;//stores the index of the worst subproblem identified
				//end of the identification of the most crowded subregion associated with those solutions in Fl

				if(worstRegion->getActualSize() > 1){//if there is more than one solution, finds the worst
					int indexWorst=-1;
					double scalarWorst=MAXDOUBLE*-1;
					for(int p=0;p<worstRegion->getActualSize();p++){
						double dist=getScalarValue(worstRegion->getSolution(p).objectiveVector, worstRegion->weight); 
						if(dist>scalarWorst){
							scalarWorst=dist;
							indexWorst=p;
						}
					}
					
					Solution sol=worstRegion->getSolution(indexWorst);
					worstRegion->exclude(indexWorst);
					updateDominationRanksRemove(sol);
				}else{ //if there is only one, delete from other region
					moead_ddRemoveWorst();
				}
			}
		}
	}
// 	exit(1);
}

void moead_ddUpdate(Swarm &sw){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
		for(int p=0;p<sw.getSize();p++){
			//use local optimizer for tandem problems on decomposition
			if(!strncmp(problemName, "tandem", 6) && decomposition){
				if(rand()/(double)RAND_MAX > 0.5)//reduce risk of local optimum
					sw.particles[p].solution.localOptimizer(sw.repository.weight);
			}
			
			moead_ddInsert(sw.particles[p].solution);
		}
		
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "MOEA/D-DD updated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void printAllForRoberto(){
	char _fAll[300]; //name of file to store the optimal solutions
	char _sAll[300]; //name of file to store the optimal solutions
	
	sprintf(_sAll, "%s_all",_s);
	sprintf(_fAll, "%s_all",_f);

	
	
	for(int s=0;s<swarmNumber;s++)//for all the swarms
		for(int p=0;p<swarms[s].repository.getActualSize();p++){//for all the particles in a given swarm repository
			printVectorToFile(swarms[s].repository.getSolution(p).decisionVector, decisionNumber, _sAll);
			printVectorToFile(swarms[s].repository.getSolution(p).objectiveVector, objectiveNumber, _fAll);
		}
	
	insertBlankLine(_sAll);
	insertBlankLine(_fAll);
}