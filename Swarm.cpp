struct Neighbor{
	int index;
	double distance;
};
// void cma_es(Swarm &sw);
void r_BOA(Swarm &sw, int swarm_edaSize, bool init);
class Swarm{
	public:
		CMAESModel model;
		
		Particle* particles=NULL;
		
		Neighbor* neighborhood=NULL; //neighborhood for MOEA/D like approach
		int neighborhoodSize;
		
		Repository repository; //repository of the swarm
		double *centroid=NULL; //centroid of the swarm in I-Multi like approaches
		
		Swarm();
		~Swarm();
		int getSize();
		void setSize(int size);
		
		//initialize the particles of the population
		void initializeParticles();
		//special method to initialize the particles of the population setting as initial local leader a ramdom solution from a given set
		void initializeParticles(Solution* candidateSet, int candidateNumber);
		//initialize the repository with the first non-dominated solutions
		void initializeRepository();
		//choose the global leaders according to the predefined strategy (PSO)
		void chooseGlobalLeaders(Repository &rep);
		//calculate the velocity of the particles (PSO)
		void calculateVelocity();
		//update the position of the particles based on its velocity, leaders and current positions (PSO)
		void updatePosition();
		//apply the turbulence factor in a predefined percentage of the particles (PSO)
		void turbulence();
		//evaluate the particles according to its current position (PSO)
		void evaluation();
		//update the repository with the best new non-dominated solutions
		void updateRepository();
		//update the local leader of the particles
		void updateParticlesMemory();
		
		void PSO();
		
		void rBOA();
		
		void CMAES();
		
		void pareto_CMAES();
		
		void DifferentialEvolution();
		
		void SBX();
		
		bool init;
	private:
		int swarmSize; //population
		Swarm(const Swarm &source){}//copy
		Swarm& operator= (const Swarm &source){return *this;}//assignment
};

Swarm::Swarm(){
	swarmSize=0;
	centroid = new double[decisionNumber+objectiveNumber];
	repository.outerSwarm=this;
	particles=NULL;
	neighborhood=NULL;
	init=true;	
}

Swarm::~Swarm(){
	delete[] centroid;
	repository.outerSwarm=NULL;
	if(particles != NULL){
		delete[] particles;
		particles=NULL;
	}
	if(neighborhood != NULL){
		delete[] neighborhood;
		neighborhood=NULL;
	}
}

int Swarm::getSize(){
	return swarmSize;
}

void Swarm::setSize(int size){
	if(swarmSize == size)//if both are the same size, skip
		return;

	if(particles != NULL)
		delete[] particles;

	swarmSize=size;
	particles = new Particle[swarmSize];
}

void Swarm::initializeParticles(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	funcEvals=0;

	for(int p=0;p<swarmSize;p++){
		sprintf(particles[p].leaderType, "%s", leader);
		particles[p].initialize();
		particles[p].weight=repository.weight;//particle has a pointer pointing to the weight vector of its repository
	}
	if(repository.getActualSize() > 0)
		repository.clear();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Particles initialized in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::initializeParticles(Solution* candidateSet, int candidateNumber){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		sprintf(particles[p].leaderType, "%s", leader);
		particles[p].initialize(candidateSet, candidateNumber);
		particles[p].weight=repository.weight;//particle has a pointer pointing to the weight vector of its repository
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Particles initialized in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::initializeRepository(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	repository.initialize(archiver, repositorySize);
// 	for(int p=0;p<swarmSize;p++){
// 		repository.add(particles[p].solution); //tries to insert the solutions in the repository
// 	}
// 	repository.organize(); //finishes with an organized reposiory
	updateRepository();
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Repository initialized in ms =\t%03.2f\n\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::chooseGlobalLeaders(Repository &rep){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	if(rep.getActualSize()==0){
		fprintf(stderr,"\nERROR ON CHOOSING GLOBAL LEADER! EMPTY REPOSITORY\n");
		exit(1);
	}
	
	if(!strcmp(particles[0].leaderType, "cd") || !strcmp(algorithm, "hmopso")) //if the selection method is crowding distance or hmopso, which uses several leader methods, then update the crowding distances
		updateCrowdingDistances(rep.getSolutions(), rep.getActualSize());
	
	for(int p=0;p<swarmSize;p++)
		particles[p].chooseGlobalLeader(rep.getSolutions(), rep.getActualSize());
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Global leaders chosen in in ms=\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::calculateVelocity(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++)
		particles[p].computeSpeed();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Velocities computed in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::updatePosition(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		particles[p].updatePosition();
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Positions updated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::turbulence(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
// 		if(p % 6 == 0){ //original of JMetal
		if( (rand()/(double)RAND_MAX) < 0.15 ){ //15% mutation
			particles[p].truncatePositionIMulti(particles[p].solution.decisionVector, range, repository.smallestDecision, repository.largestDecision);//truncate only to problem bounds before turbulence
			particles[p].turbulence();
			if(p==0 && !strcmp(algorithm, "cmaes-mopso"))
				for(int i=0;i<decisionNumber;i++)
					model.prevSol[i]=particles[p].solution.decisionVector[i];//update the previous solution for the mo-cma-es update
		}
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Turbulence applied in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::evaluation(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
// 	for(int o=0;o<objectiveNumber;o++){//reset the largest and smallest values before evaluation
// 		maxObjectives[o]=MAXDOUBLE*-1;
// 		minObjectives[o]=MAXDOUBLE;
// 	}
// 	updateMaxMinObjs();
	if(!strcmp(truncType, "rdmInSwarm"))
		updateLargestSmallestDecVectors(*this);
	for(int p=0;p<swarmSize;p++){
		if(!strcmp(algorithm, "imulti") || !strcmp(algorithm, "cmulti") || !strcmp(algorithm, "cmaes-mopso")){
			//tmp_trunc[s]+=(double)swarm[s].particles[p].truncatePositionIMulti(swarm[s].centroid, range, swarm[s].repository.smallerDecision, swarm[s].repository.largerDecision); //to count how many solutions were truncated
			particles[p].truncatePositionIMulti(centroid, range, repository.smallestDecision, repository.largestDecision);
		}else
			particles[p].truncatePositionIMulti(particles[p].solution.decisionVector, range, repository.smallestDecision, repository.largestDecision);
		
		
// 		double objectives[objectiveNumber];
// 		for(int i=0;i<objectiveNumber;i++)
// 			objectives[i]=particles[p].solution.objectiveVector[i];
		
		particles[p].solution.evaluate();
		
		// 		if(!diversityPhase){
		// 			//truncDims/=swarmSize;
		// 			tmp_trunc[s]/=swarmSize;
		// 			//printf("%f ",truncDims);
		// 		}
	}
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Solutions evaluated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::updateParticlesMemory(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++)
		particles[p].updateLocalLeader();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Local leaders updated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::updateRepository(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	//use local optimizer for tandem problems on decomposition
	if(!strncmp(problemName, "tandem", 6) && decomposition){
		for(int p=0;p<swarmSize;p++)
			particles[p].solution.localOptimizer(repository.weight);
	}
	
	repository.particlesEntered=0; //initialize the counters of numbers of particles that enters in the repository per iteration
	for(int p=0;p<swarmSize;p++){
		int updatedSolutions=0;
		if(repository.add(particles[p].solution)){ //tries to insert the solutions in the repository
			updatedSolutions++;
			totalUpdates++;
		}
		if(decomposition){
			if(useGlobalRepOnDecomposition)
				repGlobal->add(particles[p].solution);//tries to add the particle to the global repository
			
			int solsUsed=-1;
			if( (rand()/(double)RAND_MAX) < delta || neighborhoodSize == 1 )//uses only the neighborhood
				solsUsed=neighborhoodSize;
			else
				solsUsed=swarmNumber;
			
			//shuffles the neighborhood so when updating a limited number of solutions, there is no bias
			int neighborhoodShuffled[solsUsed];
			for(int i=0;i<solsUsed;i++)
				neighborhoodShuffled[i]=neighborhood[i].index;//shuffles the neighborhood
			std::random_shuffle( neighborhoodShuffled+1, neighborhoodShuffled+solsUsed ); //do not include the first neighbor (itself) on the shuffling
			
			for(int i=1;i<solsUsed;i++){ //tries to insert the particle in all swarms, except the first, which is this one and is already inserted
// 				int neighbor=neighborhood[i].index;
				int neighbor=neighborhoodShuffled[i];
				if(updatedSolutions < maxReplacements){
					if(swarms[neighbor].repository.add(particles[p].solution)){
						updatedSolutions++;
						totalUpdates++;
// 						printf("swarm %d updated neighbor %d -- %d\n", this->neighborhood[0].index, neighbor, updatedSolutions);
// 						swarms[neighbor].repository.organize();//organize the repository if the particle entered
					}
				}
			}
		}
	}
// 	repository.organize();
	
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Repository updated in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::PSO(){
	
	updateParticlesMemory();
	
	if(decomposition){
		Repository repTemp;
		mergeNeighboringRepositories(repTemp, *this);
		chooseGlobalLeaders(repTemp);
	}else
		chooseGlobalLeaders(repository);
	
	calculateVelocity();
	
	updatePosition();
	
// 	turbulence();//now it is outside, for all algorithms
}

void Swarm::rBOA(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	for(int p=0;p<swarmSize;p++){
		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
			fprintf(stderr,"\n ERROR! NaN on decision vector before EDA application\n");
			exit(1);
		}
		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
			fprintf(stderr,"\n ERROR! NaN on objective vector before EDA application\n");
			exit(1);
		}
	}
	
	if(repository.getActualSize() > 2 || decomposition){ //if there are enough particles, the solutions are updated normally (not overriten)
		r_BOA(*this, swarmSize, init);
		if(init) init=false;
	}
	
	for(int p=0;p<swarmSize;p++){
		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
			fprintf(stderr,"\n ERROR! NaN on objective vector after EDA application\n");
			exit(1);
		}
		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
			fprintf(stderr,"\n ERROR! NaN on decision vector after EDA application\n");
			exit(1);
		}
	}
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "rBOA processed in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}
	
void Swarm::CMAES(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
// // 	for(int p=0;p<swarmSize;p++){
// // 		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
// // 			fprintf(stderr,"\n ERROR! NaN on decision vector before EDA application\n");
// // 			exit(1);
// // 		}
// // 		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
// // 			fprintf(stderr,"\n ERROR! NaN on objective vector before EDA application\n");
// // 			exit(1);
// // 		}
// // 	}
	
// 	cma_es(*this);
	//******************  treating the input data ********************//
	bool rankMuUpdate=true;
	
	Repository repTemp;
	repTemp.initialize(archSubSwarms, getSize()+repository.getActualSize());
	bool suc=false;//current solution not better than the previous
	bool *injected;
	
	if(rankMuUpdate){//if it is rankMuUpdate
					
		//gathering the solutions to be used to learn and putting them on repTemp
		if(decomposition){
			mergeNeighboringRepositories(repTemp, *this);
		}else{
			if(!strcmp(solSet, "population") || !strcmp(solSet, "both"))
				for(int i=0;i<getSize();i++)
					if(!particles[i].solution.dominated)
						repTemp.add(particles[i].solution);
					
			if(!strcmp(solSet, "repository") || repTemp.getActualSize()==0 || !strcmp(solSet, "both"))
				repTemp.add(repository);
		}
		
		
		if(decomposition){
			updateWeightedDistances(repTemp.getSolutions(), repTemp.getActualSize(), repository.weight);
			std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), weightedDistanceComparatorSol);
		}else{
			if(!strcmp(cmaRanking, "cd")){
				updateCrowdingDistances(repTemp.getSolutions(), repTemp.getActualSize());
				std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), crowdingComparatorSol);
			}
			if(!strcmp(cmaRanking, "hv")){
				updateContributingHypervolume(repTemp); //stores in the crowding distance field, since both are the higher the contribution, the better, there is no problem
				std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), crowdingComparatorSol);
			}
			if(!strcmp(cmaRanking, "r2")){
				if(refSize==-1){
					perror("ERROR ON CMA-ES R2 RANKING! Reference file not set.");
				}
				updateContributingR2(repTemp); //stores in the crowding distance field, stores in the crowding distance field //if the value of the r2 increase after I take this solution, it is important, so the higher is the contribution, the better
				std::sort(repTemp.getSolutions(), repTemp.getSolutions()+repTemp.getActualSize(), crowdingComparatorSol);
			}
		}
		
		//**************************** end of treatment of input data***************//
		
		injected = new bool[repTemp.getActualSize()];
		for(int s=0;s<repTemp.getActualSize();s++){ //for all solutions
			injected[s]=true;
			for(int r=0;r<getSize();r++){//check if the current solution is in the population from the previous generation, otherwise it was injected
				if(isEqual(repTemp.getSolution(s).decisionVector, particles[r].solution.decisionVector, decisionNumber)){//if this solution is in the population, it is not injected
					injected[s]=false;
					break;
				}
			}
		}

		model.rankMuUpdate(repTemp, injected);
		
	}else{//if rankOneUpdate
		//verify if the current solution is in the repository, if it is, it was better than the previous solution, otherwise the previous solution was better and the mutation was not successfull
		for(int i=0;i<getSize();i++){
			for(int j=0;j<repository.getActualSize();j++){
				if(isEqual(particles[i].solution.objectiveVector, repository.getSolution(j).objectiveVector, objectiveNumber)){//compares only the objective vectors to be faster
					suc=true;//the current solution is in the repository
					break;
				}
			}
			if(suc)
				break;
		}
		model.rankOneUpdate(repository.getSolution(0), suc);
	}
	
	for(int i=0;i<getSize();i++){
		while(!model.sample(model.mean, particles[i].solution.decisionVector)){//repeat until is true
			if(rankMuUpdate)
				model.rankMuUpdate(repTemp, injected);//if there is problem on sampling, reinitializes everything and resample
			else
				model.rankOneUpdate(repository.getSolution(0), suc);
// 			fprintf(stderr, " Resample\n");
		}
	}
	
	
	if(injected != NULL)
		delete[] injected;
	
	if(!strncmp(solSet, "matrix",6) && decomposition){//if combining the covariance matrices directly (no exchange of individuals)
		//http://link.springer.com/chapter/10.1007%2F978-3-642-37192-9_52		
		int mode=atoi(&solSet[6]);
		Swarm *SwN=NULL;
		
		if(mode < 1 || mode > 4){
			fprintf(stderr, "ERROR! MODE NOT VALID solSet must be set to matrix[1-4]\n");
			exit(1);
		}
		
		if(mode==1){//only closest neighbor
			SwN=&swarms[neighborhood[1].index];//the closest neighbor
		}
		if(mode==2){//closest neighbor in 90% of the time
			if( (rand()/(double)RAND_MAX) < delta)//uses only the neighborhood
				SwN=&swarms[neighborhood[1].index];//the closest neighbor
			else
				SwN=&swarms[neighborhood[(rand()%(swarmNumber-1))+1].index];//a randowm swarm from the population, except for 0 (itself)
		}
		if(mode==3){//only in the neighborhood
			SwN=&swarms[neighborhood[(rand()%(20-1))+1].index];//a random neighbor, except for 0 (itself)
		}
		if(mode==4){//in the neighborhood in 90% of the time
			if( (rand()/(double)RAND_MAX) < delta)//uses only the neighborhood
				SwN=&swarms[neighborhood[(rand()%(20-1))+1].index];//a random neighbor, except for 0 (itself)
			else
				SwN=&swarms[neighborhood[(rand()%(swarmNumber-1))+1].index];//a randowm swarm from the population, except for 0 (itself)
		}
			
		if(!SwN->model.init && !model.init){//if both matrices are valid
			
			model.combineMatrices(SwN->model);
// 			double outputMatrix[decisionNumber][decisionNumber];
// 			for(int i=0;i<decisionNumber;i++){//Calculates: alpha(1-cCov)C+(1-alpha)(1-cCov)C_n
// 				for(int j=0;j<decisionNumber;j++){
// 	// 				outputMatrix[i][j]=alpha*(1-cCov)*sw.C[i][j] + (1-alpha)*(1-cCov)*SwN->C[i][j];
// 					outputMatrix[i][j]=alpha*model.C[i][j] + (1-alpha)*SwN->model.C[i][j];
// 					model.C[i][j]=outputMatrix[i][j];
// 				}
// 			}
		}
	}
	
// // 	for(int p=0;p<swarmSize;p++){
// // 		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
// // 			fprintf(stderr,"\n ERROR! NaN on objective vector after EDA application\n");
// // 			exit(1);
// // 		}
// // 		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
// // 			fprintf(stderr,"\n ERROR! NaN on decision vector after EDA application\n");
// // 			exit(1);
// // 		}
// // 	}
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "CMA-ES processed in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::pareto_CMAES(){
	if(showTimes) gettimeofday(&tmpTime, NULL);
	
	bool rankMuUpdate=true;
	
// 	printf("sz: %d\n",repository.getActualSize());
	
// 	for(int p=0;p<swarmSize;p++){
// 		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
// 			fprintf(stderr,"\n ERROR! NaN on decision vector before EDA application\n");
// 			exit(1);
// 		}
// 		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
// 			fprintf(stderr,"\n ERROR! NaN on objective vector before EDA application\n");
// 			exit(1);
// 		}
// 	}

// 	cma_es(*this);
	setSize(repository.getActualSize());//resize the repository, if the new size is different than the previous, the solutions are reset, which is no problem here
	
	if(rankMuUpdate){
		//each particle in the repository (non-dominated) learns from the neighborhood
		for(int s=0;s<repository.getActualSize();s++){
			//calculate the distance from the current particle to all others
			Solution *sol=&repository.getSolutions()[s];
			for(int s2=0;s2<repository.getActualSize();s2++){
				Solution *sol2=&repository.getSolutions()[s2];
				sol2->weightedDistance=calculateEuclideanDistance(sol->objectiveVector, sol2->objectiveVector, objectiveNumber);
			}
			std::sort(repository.getSolutions(), repository.getSolutions()+repository.getActualSize(), weightedDistanceComparatorSol);
			
			Repository repTemp;
			int sz=min(repository.getActualSize(), globalNeighborhoodSize);
			repTemp.initialize("pbi", sz);//initialized with pbi only to avoid dominance check. at this point all solutions are non-dominated
			for(int s2=0;s2<sz;s2++)
				repTemp.add(repository.getSolutions()[s2]);
			
			bool injected[sz]={true};
			injected[0]=false;//solution 0 is itself
			sol->model->rankMuUpdate(repTemp, injected);
			
			//generate new solutions with the updated repository
			while(!sol->model->sample(sol->model->mean, particles[s].solution.decisionVector)){//repeat until is true
				sol->model->rankMuUpdate(repTemp, injected);
// 				fprintf(stderr, " Resample (%d)\n",s);
			}
		}
		
		for(int s=0;s<getSize();s++){
			if( (rand()/(double)RAND_MAX) < 0.15 ){ //15% mutation
				Particle newPart;
				newPart.solution=particles[s].solution;
				newPart.turbulence();
				newPart.solution.evaluate();
				repository.add(newPart.solution);
			}
			if( (rand()/(double)RAND_MAX) < 0.15 ){ //15% local update
				Solution newSol;
				newSol=particles[s].solution;
				double weight[2];
				weight[0]=1;
				weight[1]=0;
				newSol.localOptimizer(weight);
				newSol.evaluate();
				repository.add(newSol);
			}
			particles[s].solution.evaluate();
			repository.add(particles[s].solution);
		}
		
	}else{
		//each particle in the repository (non-dominated) generates one solution
		for(int s=0;s<repository.getActualSize();s++){
			repository.getSolutions()[s].offspring=NULL;//clears the pointer
			particles[s].solution=repository.getSolution(s);
			repository.getSolutions()[s].offspring=&particles[s].solution;
			
			while(!particles[s].solution.model->sample(repository.getSolution(s).decisionVector, particles[s].solution.decisionVector)){//repeat until is true
				particles[s].solution.model->rankOneUpdate(repository.getSolution(s), true);//learn from its successfull parent
// 				fprintf(stderr, " Resample (%d)\n",s);
			}
		}
		
		//trying to insert the new solutions in the repository
		for(int s=0;s<getSize();s++){
// 			if( (rand()/(double)RAND_MAX) < 0.15 ) //15% mutation
// 				particles[s].turbulence();
			particles[s].solution.evaluate();
			bool suc=repository.add(particles[s].solution);
			
			if(suc){//if it is not successfull, do not update the offspring, it was discarded anyway
				for(int s2=repository.getActualSize()-1;s2>=0;s2--){
					if(isEqual(repository.getSolution(s2).objectiveVector, particles[s].solution.objectiveVector, objectiveNumber)){//comparing in the objective space to be faster
						repository.getSolutions()[s2].model->rankOneUpdate(repository.getSolution(s2), suc);
						break;
					}
				}
				
			}
		}
		
		//updating the models of the parents
		for(int s=0;s<repository.getActualSize();s++){
			if(repository.getSolutions()[s].offspring != NULL){//updates only the sigma of the parent if it generated an offspring (ignores new solutions) 
				repository.getSolutions()[s].model->rankOneOnlySigma(!repository.getSolutions()[s].offspring->dominated);//if the solution is non-dominated, it is successfull, otherwise it is not
			}
		}
	}
	
	if(!strncmp(solSet, "matrix",6)){//if combining the covariance matrices directly (no exchange of individuals)
		//http://link.springer.com/chapter/10.1007%2F978-3-642-37192-9_52
		double alpha=0.7;//relative weight, set as in: https://www.researchgate.net/publication/221007885_Particle_Swarm_CMA_Evolution_Strategy_for_the_Optimization_of_Multi-Funnel_Landscapes?enrichId=rgreq-cdc85ea6a42fcdc663a5d2f7b4c47f80-XXX&enrichSource=Y292ZXJQYWdlOzIyMTAwNzg4NTtBUzoxMDM2MDQwMDIyMzAyNzJAMTQwMTcxMjUyNzA1Ng%3D%3D&el=1_x_2
		int mode=atoi(&solSet[6]);
		
		for(int s=0;s<getSize();s++){
			Solution *sol=&repository.getSolutions()[s];
			
// 			//find the closest neighbor of the current solution
// 			Solution *closestNeighbor;
// 			double minDist=MAXDOUBLE;
// 			for(int s2=repository.getActualSize()-1;s2>=0;s2--){
// 				Solution *sol2=repository.getSolutions()[s2];
// 				double dist=calculateEuclideanDistance(sol.objectiveVector, sol2.objectiveVector, objectiveNumber);
// 				if(dist < minDisd && s != s2){//do not consider itself
// 					minDisd=dist;
// 					closestNeighbor=sol2;
// 				}
// 			}
			for(int s2=0;s2<repository.getActualSize();s2++){
				Solution *sol2=&repository.getSolutions()[s2];
				sol2->weightedDistance=calculateEuclideanDistance(sol->objectiveVector, sol2->objectiveVector, objectiveNumber);
			}
			std::sort(repository.getSolutions(), repository.getSolutions()+repository.getActualSize(), weightedDistanceComparatorSol);
			
			
		
			const Solution *neighbor=NULL;
			
			if(mode < 1 || mode > 4){
				fprintf(stderr, "ERROR! MODE NOT VALID solSet must be set to matrix[1-4]\n");
				exit(1);
			}
			
			int ngh=std::min(globalNeighborhoodSize, repository.getActualSize());//the global neighborhood size can be larger than the size of the repository
			
			if(mode==1){//only closest neighbor
				neighbor=&repository.getSolution(1);//the closest neighbor
			}
			if(mode==2){//closest neighbor in 90% of the time
				if( (rand()/(double)RAND_MAX) < delta)//uses only the neighborhood
					neighbor=&repository.getSolution(1);//the closest neighbor
				else
					neighbor=&repository.getSolution((rand()%(repository.getActualSize()-1))+1);//a random solution from the population, except for 0 (itself)
			}
			if(mode==3){//only in the neighborhood
				neighbor=&repository.getSolution((rand()%(ngh-1))+1);//a random neighbor, except for 0 (itself)
			}
			if(mode==4){//in the neighborhood in 90% of the time
				if( (rand()/(double)RAND_MAX) < delta)//uses only the neighborhood
					neighbor=&repository.getSolution((rand()%(ngh-1))+1);//a random neighbor, except for 0 (itself)
				else
					neighbor=&repository.getSolution((rand()%(repository.getActualSize()-1))+1);//a random solution from the population, except for 0 (itself)
			}
				
			if(!neighbor->model->init && !sol->model->init){//if both matrices are valid
				
				sol->model->combineMatrices(*neighbor->model);
			
// 				double outputMatrix[decisionNumber][decisionNumber];
// 				for(int i=0;i<decisionNumber;i++){//Calculates: alpha(1-cCov)C+(1-alpha)(1-cCov)C_n
// 					for(int j=0;j<decisionNumber;j++){
// 		// 				outputMatrix[i][j]=alpha*(1-cCov)*sw.C[i][j] + (1-alpha)*(1-cCov)*SwN->C[i][j];
// 						outputMatrix[i][j]=alpha*sol->model->C[i][j] + (1-alpha)*neighbor->model->C[i][j];
// 						sol->model->C[i][j]=outputMatrix[i][j];
// 					}
// 				}
			}
		}
	}
	
	setSize(0);
	
// 	for(int p=0;p<swarmSize;p++){
// 		if(vectorNan(particles[p].solution.objectiveVector, objectiveNumber)){
// 			fprintf(stderr,"\n ERROR! NaN on objective vector after EDA application\n");
// 			exit(1);
// 		}
// 		if(vectorNan(particles[p].solution.decisionVector, decisionNumber)){
// 			fprintf(stderr,"\n ERROR! NaN on decision vector after EDA application\n");
// 			exit(1);
// 		}
// 	}
	if(showTimes){
		gettimeofday(&endTime, NULL); fprintf(stderr, "Pareto CMA-ES processed in ms =\t%03.2f\n", ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
	}
}

void Swarm::DifferentialEvolution(){
	Repository repTemp;
	
	if(decomposition){
		//gathering the solutions to be used to learn and putting them on repTemp
		mergeNeighboringRepositories(repTemp, *this);
	}else{
		fprintf(stderr,"\n ERROR! Differential evolution only available for MOEA/D for now\n");
		exit(1);
	}
	
// 	repTemp.organize();
	int rnd1=-1,rnd2=-1;
	
	for(int p=0;p<swarmSize;p++){
		if(repTemp.getActualSize() > 1){
			while(true){
				rnd1=rand()%repTemp.getActualSize();
				rnd2=rand()%repTemp.getActualSize();
				if(rnd1 != rnd2)
					break;
			}
		}else{
			rnd1=rnd2=0;
		}
		
		const Solution *parent2=&repTemp.getSolution(rnd1);
		const Solution *parent3=&repTemp.getSolution(rnd2);
			
		double f = 0.5;
		//CR=1	
		for(int n=0;n<decisionNumber;n++){
			/*Selected Two Parents*/
			particles[p].solution.decisionVector[n] = repository.getSolution(0).decisionVector[n] + f*(parent2->decisionVector[n] - parent3->decisionVector[n]);
		}
	}
}

void Swarm::SBX(){
	Repository repTemp;
	
	if(decomposition){
		//gathering the solutions to be used to learn and putting them on repTemp
		mergeNeighboringRepositories(repTemp, *this);
	}else{
		fprintf(stderr,"\n ERROR! SBX crossover only available for MOEA/D for now\n");
		exit(1);
	}
	
	int rnd1=-1,rnd2=-1;
	
	for(int p=0;p<swarmSize;p+=2){
		if(repTemp.getActualSize() > 1){
			while(true){
				rnd1=rand()%repTemp.getActualSize();
				rnd2=rand()%repTemp.getActualSize();
				if(rnd1 != rnd2)
					break;
			}
		}else{
			rnd1=rnd2=0;
			printf("Warning! Only one solution in repository for SBX\n");
		}
		
		const Solution *parent1=&repTemp.getSolution(rnd1);
		const Solution *parent2=&repTemp.getSolution(rnd2);
		
		double y1,y2, betaq;
		double crossoverProbability=1.0;//just like MOEA/DD
		double distributionIndex=30.0;//just like MOEA/DD
		
		if(rand()/(double)RAND_MAX <= crossoverProbability){
			for(int i=0;i<decisionNumber;i++){
				double x1=parent1->decisionVector[i];
				double x2=parent2->decisionVector[i];
				if(rand()/(double)RAND_MAX <= 0.5){
					if(abs(x1 - x2) > 1.0e-14){//defines the minimum difference allowed between real values
						if(x1<x2){
							y1=x1;
							y2=x2;
						}else{
							y1=x2;
							y2=x1;
						}
						double rnd=rand()/(double)RAND_MAX;
						double beta = 1.0 + (2.0 * (y1 - inferiorPositionLimit[i]) / (y2 - y1));
						double alpha = 2.0 - pow(beta, -(distributionIndex + 1.0));
						
						if (rnd <= (1.0 / alpha)) {
							betaq = pow(rnd * alpha, (1.0 / (distributionIndex + 1.0)));
						} else {
							betaq = pow(1.0 / (2.0 - rnd * alpha), 1.0 / (distributionIndex + 1.0));
						}
						double c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));

						beta = 1.0 + (2.0 * (superiorPositionLimit[i] - y2) / (y2 - y1));
						alpha = 2.0 - pow(beta, -(distributionIndex + 1.0));

						if (rnd <= (1.0 / alpha)) {
							betaq = pow((rnd * alpha), (1.0 / (distributionIndex + 1.0)));
						} else {
							betaq = pow(1.0 / (2.0 - rnd * alpha), 1.0 / (distributionIndex + 1.0));
						}
						double c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
						
						if (c1 < inferiorPositionLimit[i])
							c1 = inferiorPositionLimit[i];

						if (c2 < inferiorPositionLimit[i])
							c2 = inferiorPositionLimit[i];

						if (c1 > superiorPositionLimit[i])
							c1 = superiorPositionLimit[i];

						if (c2 > superiorPositionLimit[i])
							c2 = superiorPositionLimit[i];
						
						if(rand()/(double)RAND_MAX <= 0.5){
							particles[p].solution.decisionVector[i]=c2;
							if(p+1<swarmSize)
								particles[p+1].solution.decisionVector[i]=c1;
						}else{
							particles[p].solution.decisionVector[i]=c1;
							if(p+1<swarmSize)
								particles[p+1].solution.decisionVector[i]=c2;
						}
					}else{
						particles[p].solution.decisionVector[i]=x1;
						if(p+1<swarmSize)
							particles[p+1].solution.decisionVector[i]=x2;
					}
				}else{
					particles[p].solution.decisionVector[i]=x2;
					if(p+1<swarmSize)
	 					particles[p+1].solution.decisionVector[i]=x1;
				}
			}
		}else{
			printf("NO CROSSOVER\n");
			for(int i=0;i<decisionNumber;i++){
				particles[p].solution.decisionVector[i]=parent1->decisionVector[i];
				if(p+1<swarmSize)
					particles[p+1].solution.decisionVector[i]=parent2->decisionVector[i];
			}
		}
	}
}