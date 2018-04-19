//class solution
class Solution{
public:
	Problem problem;
	
	CMAESModel *model=NULL;
	Solution* offspring=NULL;
	
	double *decisionVector=NULL;
	double *objectiveVector=NULL;
	
	double crowdingDistance;
	double weightedDistance;// stores the weightedDistance relative to the swarm
	bool dominated;
	int dominanceLevel;//stores the dominance level for algorithms that use it, like the MOEA/D-DD
	bool evalSolution;//only evaluate this solution if this flag is true
	
	Solution();//constructor
	~Solution();//destructor
	Solution(const Solution &source);//copy
	Solution& operator= (const Solution &source);//assignment
	
	//initialize the decision variables with random values
	void initialize();

	//compare two Solutions
	//param - the solution to be compared to
	bool isEqual(const Solution &sol) const;
	
	double localOptimizer(double* weight);

	//evaluate this solution
	void evaluate();
	
	//prints the decison and objective vectors of a solution
	void print();
	
};

// // // constructor randomly initialize the variables
Solution::Solution(){
	crowdingDistance = -1; //if the crowding distance is negative, it is not set
	dominated=true;
	evalSolution=true;
	dominanceLevel=-1;
	offspring=NULL;
	
	if(!strcmp(algorithm, "pareto-cmaes"))
		model= new CMAESModel;
	decisionVector = new double[decisionNumber];
	objectiveVector = new double[objectiveNumber];
	memset(decisionVector, 0, sizeof(double)*decisionNumber);
}

//Destructor
Solution::~Solution(){
	delete[] objectiveVector;
	delete[] decisionVector;
	
	if(model != NULL)
		delete model;
}

//Copy constructor
Solution::Solution(const Solution &source){
// 	if (this != &source) {
// 		if (objectiveVector != NULL) {
// // 			delete[] objectiveVector;
// 		}
// 	}
	offspring=source.offspring;
	
	if(!strcmp(algorithm, "pareto-cmaes"))
		model= new CMAESModel;
	objectiveVector = new double[objectiveNumber];
	decisionVector = new double[decisionNumber];
	
	crowdingDistance=source.crowdingDistance;
	weightedDistance=source.weightedDistance;
	dominated=source.dominated;
	dominanceLevel=source.dominanceLevel;
	evalSolution=source.evalSolution;
	
	memcpy(objectiveVector, source.objectiveVector, sizeof(double)*objectiveNumber);
	memcpy(decisionVector, source.decisionVector, sizeof(double)*decisionNumber);
// 	memcpy(&model, &source.model, sizeof(CMAESModel));
	if(!strcmp(algorithm, "pareto-cmaes"))
		*model=*source.model;
}

//Assignment operator
Solution& Solution::operator= (const Solution &source){
	offspring=source.offspring;
	
	crowdingDistance=source.crowdingDistance;
	weightedDistance=source.weightedDistance;
	dominated=source.dominated;
	dominanceLevel=source.dominanceLevel;
	evalSolution=source.evalSolution;
	
	memcpy(objectiveVector, source.objectiveVector, sizeof(double)*objectiveNumber);
	memcpy(decisionVector, source.decisionVector, sizeof(double)*decisionNumber);
// 	memcpy(&model, &source.model, sizeof(CMAESModel));
	if(!strcmp(algorithm, "pareto-cmaes"))
		*model=*source.model;
	return *this;
}

void Solution::initialize(){
	if(!strcmp(algorithm, "cmaes-mopso") || !strcmp(algorithm, "pareto-cmaes")){//if it is CMAES-MOPSO, start the solution with a gaussian
		
// 		model->sample(model->mean, decisionVector);
		
		for(int i=0;i<decisionNumber;i++){
			double mean=inferiorPositionLimit[i]+( (superiorPositionLimit[i]-inferiorPositionLimit[i])/2.0 );
			double sigma=0.3;//as in the CMA-ES code
			double normalizedValue= mean+sigma*cmaes_random_Gauss();
// 			while(normalizedValue > 1 || normalizedValue < 0)//resample until is within the bounds
// 				normalizedValue= mean+sigma*cmaes_random_Gauss();
			
			normalizedValue=std::min(normalizedValue,superiorPositionLimit[i]);//limit to the maximum value of 1
			normalizedValue=std::max(normalizedValue,inferiorPositionLimit[i]);//limit to the minimum value of 0
			
// 			decisionVector[i]=inferiorPositionLimit[i]+ (normalizedValue* (superiorPositionLimit[i]-inferiorPositionLimit[i]) );
			decisionVector[i]=normalizedValue;
		}
	}else{
		for(int i=0;i<decisionNumber;i++)
			decisionVector[i]=inferiorPositionLimit[i]+ ((rand()/(double)RAND_MAX)* (superiorPositionLimit[i]-inferiorPositionLimit[i]) );
	}
	
// 	problem.evaluate(decisionVector, objectiveVector);
}

//evaluate the solution according to a predetermined problem
void Solution::evaluate(){
	if(evalSolution){
		problem.evaluate(decisionVector, objectiveVector);
		crowdingDistance= -1;
		weightedDistance=-1;
		dominanceLevel=-1;
		dominated=true;
		for(int o=0;o<objectiveNumber;o++){
			maxObjectives[o]=std::max(maxObjectives[o], objectiveVector[o]);
			minObjectives[o]=std::min(minObjectives[o], objectiveVector[o]);
		}
	}
}

//comparator if two solutions are equal
//param - the solution to be compared to
bool Solution::isEqual(const Solution &sol) const{
	for(int i=0;i<objectiveNumber;i++){
		//if( sol.objectiveVector[i] != objectiveVector[i] ){
		if( round_5(sol.objectiveVector[i]) != round_5(objectiveVector[i]) ){
			return false;
		}
	}
	return true;
}

void Solution::print(){
	printVector(objectiveVector,objectiveNumber);
	printf("-> ");
	printVector(decisionVector,decisionNumber);
	printf("\n");
}

double Solution::localOptimizer(double* weight){
// 	printf("localOptimizer\n\n");
	double  rhobeg = 0.5, rhoend = 0.0001, rc=0;
	int iprint=0, maxfun = 1000;
// 	int iprint=0, maxfun = 10;
	example_state state;
	state.nprob = 0;
	
	if(objectiveVector[0]<1999){
		globalTmpWeight=weight;
		rc = cobyla(decisionNumber, 2*decisionNumber, decisionVector, rhobeg, rhoend, iprint, &maxfun, localOptimizerObjFunc, &state);
	}
	return rc;
}