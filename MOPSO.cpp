#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <algorithm>
#include <unistd.h>
#include <iostream>

#include <cstdlib>
#include <sstream>
#include <vector>
#include <cassert>

// #include "additionalCode/asa_tandem/asa_cg.c"
#include "additionalCode/cma-es/eda.h"
#include "additionalCode/cobyla/cobyla.c"
#include "variables.cpp"
#include "Problem.cpp"
#include "Solution.cpp"
#include "Particle.cpp"
#include "Repository.cpp"
#include "additionalCode/cma-es/eda.cpp"
#include "Swarm.cpp"
#include "util.cpp"
#include "additionalCode/rBOA/eda.cpp"

//*********************************************************************************************************************************************************************//
bool run(){
	gettimeofday(&startTime, NULL);
	long iteration=0;
	clearBetweenRuns();
	fprintf(stderr, " %s - %d OBJ - %d DEC \n", problemName, objectiveNumber, decisionNumber);
	
	//Initialization of the runs
	for(int s=0;s<swarmNumber;s++){
		swarms[s].initializeParticles();//initialize the particles
		swarms[s].initializeRepository(); //initialize the repository with the best particles found in the initialization
	}

	//Main loop
	while(true){
		if(maxEvals < 0){
			if(iteration >= maxIterations)//stop by max iterations
				break;
		}else
			if(maxEvals <= funcEvals)//stop by max evaluations
				break;
		
// 		gettimeofday(&tmpTime, NULL);
		hmopsoBefore();//if H-MOPSO
		IMultiOperations(iteration); //if imulti or cmulti
		
		
// // 		if(!strcmp(algorithm, "cmaes-mopso")){
// // // 			int pp=4+floor(3*log(decisionNumber));
// // 			int pp=2;
// // 			if(swarms[0].getSize() < pp){
// // 				fprintf(stderr,"set to %d\n",pp);
// // 				for(int sw=0;sw<swarmNumber;sw++){
// // 					swarms[sw].setSize(pp);
// // 				}
// // 			}
// // 		}
		if(!strcmp(algorithm, "moead") || !strcmp(algorithm, "cmaes-mopso")){
// 			int pp=2;
			int pp=4+floor(3*log(decisionNumber));
			if(swarms[0].getSize() < pp){
				fprintf(stderr,"set to %d\n",pp);
				for(int sw=0;sw<swarmNumber;sw++){
					swarms[sw].setSize(pp);
				}
			}
		}
		
		//shuffles the neighborhood so there is no bias
		int neighborhoodShuffled[swarmNumber];
		for(int i=0;i<swarmNumber;i++)
			neighborhoodShuffled[i]=i;//shuffles the neighborhood
		std::random_shuffle( neighborhoodShuffled+0, neighborhoodShuffled+swarmNumber );
		for(int s=0;s<swarmNumber;s++){ //for all the swarms
			int sw=neighborhoodShuffled[s];
			if(diversityPhase || !strcmp(algorithm, "imulti") || !strcmp(algorithm, "smpso") || !strcmp(algorithm, "hmopso")) //if diversity phase or imulti or smpso
				swarms[sw].PSO();
			
			if(!strcmp(algorithm, "cmulti") && !diversityPhase)
				swarms[sw].rBOA();
			
			if(!strcmp(algorithm, "cmaes-mopso") && !diversityPhase)
				swarms[sw].CMAES();
			
			if(!strcmp(algorithm, "pareto-cmaes") && !diversityPhase){
				swarms[sw].pareto_CMAES();
				printf("sz: %d, curr_eval/max_eval: %d/%d\n",swarms->repository.getActualSize(), (int)funcEvals, (int)maxEvals);
			}
			
			if(!strcmp(algorithm, "moead") && !diversityPhase)
// 				swarms[sw].DifferentialEvolution();
				swarms[sw].SBX();
			
			swarms[sw].turbulence();
			
			swarms[sw].evaluation();
			
			//MOEA/DD UPDATE
			if(iteration==0)
				calculateDominationRanks();
			moead_ddUpdate(swarms[sw]);
			
// 			swarms[sw].updateRepository(); //repository update for other algorithms
		}
		
		if(printMOEAD_CMAES && !strcmp(algorithm, "cmaes-mopso") && !diversityPhase){
			int sz=8;
			double tmp[sz];
			tmp[0]=averageWeightedLogLikelihood();
			tmp[1]=averageScalarObjFromPopulation();
			tmp[2]=averageScalarObjFromRepository();
			tmp[3]=likelihoodQuality(1);
			tmp[4]=likelihoodQuality(2);
			tmp[5]=likelihoodQuality(3);
			tmp[6]=likelihoodQuality(4);
			
			repTemp->initialize(archiver, repositorySize*swarmNumber);
			for(int s=0;s<swarmNumber;s++)
				repTemp->add(swarms[s].repository);
			tmp[7]=calculateIGDp(*repTemp);
			
			
// 			if(bestTandem == MAXDOUBLE)
// 				tmp[3]=0;
// 			else
// 				tmp[3]=bestTandem;
			printVectorToFile(tmp, sz, _MOEAD_CMAES);
			// 			printf("%f %f %f\n", averageWeightedLogLikelihood(), averagePbiFromPopulation(), averagePbiFromRepository());
		}
		
		// 			for(int sw=0;sw<swarmNumber;sw++)//number of solutions kept from previous iteration by any subproblem (in our case 1 or 0)
		// 				printf("%d ", swarms[sw].solsKept);
		// 			printf("\n");

// // // // 		if(strcmp(algorithm, "moead-dd")){//if it is not moea/d-dd
// // // // 			for(int sw=0;sw<swarmNumber;sw++){
// // // // 				swarms[sw].updateRepository();
// // // // 			}// END of the swarms loop
// // // // 		}//else{//if it is moea/d-dd

// 			moead_ddUpdate();
			
// 			for(int sw=0;sw<swarmNumber;sw++){
// 				printf("%d ",swarms[sw].repository.getActualSize());
// 			}
// 			printf("\n");

// // // 		}
// 		gettimeofday(&tmpTime, NULL);
		updateNeighborhood(); //update the neighbors of each sub-problem based on different criterions (set as parameter)
// 		updateModelsQuality(); //update the quality of a model, only good models will be used as surrogates
// 		gettimeofday(&endTime, NULL); fprintf(stderr, "Neighborhood updated in ms =\t%03.2f\n",((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);

		
		hmopsoAfter();
		
		mostrador(iteration);
		iteration++;
		
// 		printAllForRoberto();//Print all the solutions at all the iterations. Roberto need this data
		
// 		gettimeofday(&endTime, NULL); fprintf(stderr, "Iteration %d completed in ms =\t%03.2f\n",iteration, ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000);
// 		gettimeofday(&endTime, NULL); fprintf(stderr, "Iteration %d completed in ms =\t%03.2f with %d sols\n",iteration, ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)tmpTime.tv_sec * 1e6 + tmpTime.tv_usec))/1000, repGlobal->getActualSize());
	}// END OF THE MAIN LOOP
	finalizeRun(startTime);

	return true;
}//END of run

//*********************************************************************************************************************************************************************//

void clearBetweenRuns(){
	if(decisionNumber < 0){
		fprintf(stderr,"\nERROR! Number of decision variables is: %d. Are you sure you are using everything correctly? \n", decisionNumber);
		exit(1);
	}
	
	if( ( !strcmp(algorithm, "imulti") || !strcmp(algorithm, "cmulti") || !strcmp(algorithm, "cmaes-mopso") || !strcmp(algorithm, "moead") ) ){//reset the parameters of multi-swarm algorithms
		range=rangeInicial;
		increment=fabs(rangeInicial-rangeFinal)/(numParticoes-1)*-1;//(I-Multi)
		diversityPhase=true;
		decomposition=false;
		swarmNumber=1;
		//sprintf(archiver, "mga");
		//sprintf(archiver, "cd");
		//sprintf(leader, "cd");
	}
	if(!strcmp(algorithm, "hmopso")){//reset the roulette of H-MOPSO
		for(int i=0;i<combinations;i++)
			probabilities[i]=1.0/combinations;
		sprintf(archiver, "cd"); //initialize with smpso, only change to other when the repository is full
		sprintf(leader, "cd");
	}
	if(repTemp != NULL)//a temporary repository to be used in different situations
		delete repTemp;
	repTemp = new Repository;

	if(repGlobal != NULL)
		repGlobal->clear();
	
	if(deltaMin != NULL)
		delete[] deltaMin;
	if(deltaMax != NULL)
		delete[] deltaMax;
	deltaMin = new double[decisionNumber];
	deltaMax = new double[decisionNumber];
	for (int i = 0; i < decisionNumber; i++) {
		deltaMax[i] = (superiorPositionLimit[i] - inferiorPositionLimit[i]) / 2.0;
		deltaMin[i] = -deltaMax[i];
	}
	if(minObjectives != NULL)
		delete[] minObjectives;
	if(maxObjectives != NULL)
		delete[] maxObjectives;
	maxObjectives = new double[objectiveNumber];
	minObjectives = new double[objectiveNumber];
	for(int o=0;o<objectiveNumber;o++){
		maxObjectives[o]=MAXDOUBLE*-1;
		minObjectives[o]=MAXDOUBLE;
	}
	
	if(swarms != NULL) //if there is already memory allocated, re-alocate
		delete[] swarms;
	swarms = new Swarm[swarmNumber];
	if( (!strcmp(archSubSwarms, "pbi")  || !strcmp(archSubSwarms, "tch") || !strcmp(archSubSwarms, "wcp") || !strcmp(archSubSwarms, "wsum") || !strcmp(archSubSwarms, "r-ideal")) && diversityIterations<1)//avoid unecessary particle initializations
		swarms[0].setSize(0);
	else
		swarms[0].setSize(originalPopSize);
}

void IMultiOperations(long iteration){
	if( !strcmp(algorithm, "imulti") || !strcmp(algorithm, "cmulti" ) || !strcmp(algorithm, "cmaes-mopso") || !strcmp(algorithm, "moead") ){
		//end of the diversity phase
		if(iteration == diversityIterations && diversityPhase){
			if(!strcmp(archSubSwarms, "pbi")  || !strcmp(archSubSwarms, "tch") || !strcmp(archSubSwarms, "wcp") || !strcmp(archSubSwarms, "wsum") || !strcmp(archSubSwarms, "r-ideal")){
				decomposition=true;
				if(repGlobal != NULL)
					delete repGlobal;
				repGlobal = new Repository;
				if(useGlobalRepOnDecomposition)
					repGlobal->initialize("cd", 10000);//one thousand solutions
			}
			
			diversityPhase=false; //diversity phase of the IMulti algorithm
			swarmNumber=originalSwarmNumber;
			
			repTemp->initialize(archiver, repositorySize*swarmNumber);
			repTemp->add(swarms[0].repository);
			
// // 			if( !strcmp(algorithm, "imulti"))
// // 				//sprintf(archiver, "ideal");
// // 				strcpy(archiver, archSubSwarms);
// // 			if( !strcmp(algorithm, "cmulti") || !strcmp(algorithm, "cmaes-mopso"))
// // 				strcpy(archiver, archSubSwarms);
// // 			//sprintf(archiver, archSubSwarms);
			
			if(swarms != NULL) //if there is already memory allocated, re-alocate
				delete[] swarms;
			swarms = new Swarm[swarmNumber];
			for(int s=0;s<swarmNumber;s++){ //initialize all the sub-swarms accordingly
				swarms[s].init=true;
				swarms[s].repository.initialize(archSubSwarms, repositorySize);
				// 				swarms[s].repository.initialize(archiver, 1);
				if(decomposition){
					memcpy(swarms[s].repository.weight, reference[s], sizeof(double)*objectiveNumber);
				}
			}
			if(decomposition)
				initializeNeighborhood();
			
			if(repTemp->getActualSize() < swarmNumber){ //if there is not enough solutions for all swarms, reduce the number of swarms
				if(!decomposition)//if not decomposition, limit the number of swarms
					swarmNumber=repTemp->getActualSize();
			}
			
			if(swarmNumber>1){
				subSwarmPopulation=(subSwarmTotalPopulation/swarmNumber);
				if(subSwarmPopulation < 1) subSwarmPopulation=1;
			}else
				subSwarmPopulation=originalPopSize;
			
			for(int s=0;s<swarmNumber;s++){
				if(s < subSwarmTotalPopulation-(swarmNumber*subSwarmPopulation))//deal with non-integer division of subSwarmPopulation by swarmNumber
					swarms[s].setSize(subSwarmPopulation+1);
				else
					swarms[s].setSize(subSwarmPopulation);
			}
			
			if(decomposition)//if decomposition, cluster the solutions according to the weights
				weightClustering(repTemp->getSolutions(), repTemp->getActualSize());
			else
				KMeans(repTemp->getSolutions(), repTemp->getActualSize()); // otherwise cluster using kmeans
				
			for(int s=0;s<swarmNumber;s++){
				//swarms[s].repository.organize();
				swarms[s].initializeParticles(swarms[s].repository.getSolutions(), swarms[s].repository.getActualSize());//initialize the particles
			}
			
			// 			for(int i=0;i<swarmNumber;i++){
			// 				for(int p=0;p<swarms[i].repository.getActualSize();p++){
			// 					printVector(swarms[i].repository.getSolution(p).objectiveVector, objectiveNumber);
			// 					printf("\n");
			// 				}
			// 				printf("\n\n\n");
			// 			}
		}
		//partitioning iteration
		if( !diversityPhase && maxIterations < 0 && iteration>0 && iteration != diversityIterations && ( ((iteration-diversityIterations) % ( (maxIterations-diversityIterations)/numParticoes) ) == 0) && swarmNumber > 1 ){
			range+=increment; //calculate new range
			
			//printf("\n swNumber: %d, swSize: %d \n", swarmNumber, swarmSize);
			repTemp->initialize(archiver, repositorySize*swarmNumber);
			for(int s=0;s<swarmNumber;s++){ //merge all solutions into a single repository
				repTemp->add(swarms[s].repository);
				//printf("\n %d particles", repIntermediate.getActualSize());
			}
			
			swarmNumber=originalSwarmNumber; //restore the swarmNumber, if it was changed previously
			
			for(int s=0;s<swarmNumber;s++){
				//cma-es init
				swarms[s].init=true; //reset the cma-es of each sub-swarm
				swarms[s].repository.clear(); //clear the repositories
			}
			
			//printf("\npartitioning iteration - %d\n",repTemp.getActualSize());
			
			if(repTemp->getActualSize() < swarmNumber){
				if(!decomposition)//if not decomposition, limit the number of swarms
					swarmNumber=repTemp->getActualSize();
			}
			
			// 			printf("%d\n", swarmNumber);
			
			if(swarmNumber>1){
				subSwarmPopulation=(subSwarmTotalPopulation/swarmNumber);
				if(subSwarmPopulation < 1) subSwarmPopulation=1;
			}else
				subSwarmPopulation=originalPopSize;
			
			for(int s=0;s<swarmNumber;s++){
				if(s < subSwarmTotalPopulation-(swarmNumber*subSwarmPopulation))//deal with non-integer division of subSwarmPopulation by swarmNumber
					swarms[s].setSize(subSwarmPopulation+1);
				else
					swarms[s].setSize(subSwarmPopulation);
			}
			
			if(decomposition)//if the weights are in use, cluster the solutions according to them
				weightClustering(repTemp->getSolutions(), repTemp->getActualSize());
			else
				KMeans(repTemp->getSolutions(), repTemp->getActualSize()); // otherwise cluster using kmeans
				//KMeans(repTemp.getSolutions(), repTemp.getActualSize(), swarm, swarmNumber);
				
				for(int s=0;s<swarmNumber;s++){
// 					swarms[s].repository.organize();
					swarms[s].initializeParticles(swarms[s].repository.getSolutions(), swarms[s].repository.getActualSize());//initialize the particles
				}
		}
	}
}

void hmopsoBefore(){
	if(!strcmp(algorithm, "hmopso") && swarms[0].repository.archiverUsed){
		//selects a combination to use
		double rnd=(rand()/(double)RAND_MAX);
		double sum=0;
		used=-1;
		for(int i=0;i<combinations;i++){ //roulette
			sum+=probabilities[i];
			used=i;
			if(rnd <= sum)
				break;
		}
		switch(used){
			case 0:
				sprintf(swarms[0].repository.archiverType, "ideal");//ideal
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "cd");//crowding distance
					break;
			case 1:
				sprintf(swarms[0].repository.archiverType, "ideal");//ideal
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "nwsum");//nwsum
					break;
			case 2:
				sprintf(swarms[0].repository.archiverType, "ideal");//ideal
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "sigma");//sigma
					break;
			case 3:
				sprintf(swarms[0].repository.archiverType, "mga");//mga 
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "cd");//crowding distance
					break;
			case 4:
				sprintf(swarms[0].repository.archiverType, "mga");//mga
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "nwsum");//nwsum
					break;
			case 5:
				sprintf(swarms[0].repository.archiverType, "mga");//mga
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "sigma");//sigma
					break;
			case 6:
				sprintf(swarms[0].repository.archiverType, "cd");//crowding distance
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "cd");//crowding distance
					break;
			case 7:
				sprintf(swarms[0].repository.archiverType, "cd");//crowding distance
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "nwsum");//nwsum
					break;
			case 8:
				sprintf(swarms[0].repository.archiverType, "cd");//crowding distance
				for(int i=0;i<swarms[0].getSize();i++)
					sprintf(swarms[0].particles[i].leaderType, "sigma");//sigma
					break;
			default:
				fprintf(stderr,"Error on switch - H-MOPSO");
				exit(1);
				break;
		}
		
		repTemp->initialize(archiver, swarms[0].repository.getActualSize());
		repTemp->add(swarms[0].repository);
	}	
}

void hmopsoAfter(){
	if(!strcmp(algorithm, "hmopso") && used != -1){
		double vectmp[combinations+1];
		for(int i=0;i<combinations;i++)
			vectmp[i]=probabilities[i];
		vectmp[combinations]=used;
		printVectorToFile(vectmp, combinations+1, _p);
		//double pMinima=0.05;//minimum probability
		//double pMinima=0.01;//minimum probability
		double pMinima=0.005;//minimum probability
		
		//update the probabilities
		double R2_ant=calculateR2(swarms[0].repository, reference, refSize);
		double R2=calculateR2(swarms[0].repository, reference, refSize);
		
		double inc=(1.0/combinations)/10;
		////double inc=(1.0/combinations)/5;
		//double inc=fabs(R2-R2_ant)*10;
		
		if(R2_ant >= R2){ //better or equal
			probabilities[used]+=inc;
			for(int i=0;i<combinations;i++){
				if((probabilities[i]-(inc/combinations)) > pMinima)
					//if((probabilities[i]-(inc/combinations)) > 0)
					probabilities[i]-=(inc/combinations);
				else
					probabilities[used]-=(inc/combinations);
			}
		}else{ //worsen
			//if(R2_ant < R2){ //worsen
			if(probabilities[used]-inc < pMinima){
				inc=probabilities[used]-pMinima;
				probabilities[used]=pMinima;
			}else
				probabilities[used]-=inc;
			for(int i=0;i<combinations;i++){
				if((probabilities[i]+(inc/combinations)) < 1-(pMinima*(combinations-1)) ){
					//if((probabilities[i]+(inc/combinations)) < pMinima*(combinations-1)){ //never enters
					//if((probabilities[i]+(inc/combinations)) < 1) //more than 1
					probabilities[i]+=(inc/combinations);
				}else{
					probabilities[used]+=(inc/combinations);
				}
			}
			//if worsened, restores the repository copy (only improving)
			//memcpy(&swarms[0].repository, &repPrevious, sizeof(Repository));
			
			swarms[0].repository.clear();
			swarms[0].repository.add(repTemp[0]);
			
			// 			//show
			// 			printf("\n%f\t", R2-R2_ant);
			// 			double s=0;
			// 			for(int i=0;i<combinations;i++){
			// 				//printf("%.2f%% ",probabilities[i]*100.0);
			// 				printf("%f ",probabilities[i]);
			// 				s+=probabilities[i];
			// 			}
			// 			printf("(%f) - %d ",s, used);
		}
	}	
}

//show the progress of the optimization
//param - number of the actual iteration of the algorithm
void mostrador(long iteration){
	if(maxEvals < 0){//stop by max iterations
		if( (int)fmod( (iteration+1), (maxIterations/20.0)) == 0){
			fprintf(stderr, "%d%% ", (int)floor( ( (iteration+1.0)/maxIterations) *100) );
		}
	}else{//stop by max evaluations
// // 		printf("%ld / %ld -- current: %d max: %d\n", funcEvals, maxEvals, swarms->repository.getActualSize(), swarms->repository.getMaxSize());
		
		int te=0;//total evaluations
		for(int s=0;s<swarmNumber;s++)
			te+=swarms[s].getSize();
			
		int mi=maxEvals/(double)te;//max iterations
		if( (int)fmod( (iteration+1), (mi/20.0)) == 0){
			fprintf(stderr, "%d%% ", (int)floor( ( (iteration+1.0)/mi) *100) );
		}
	}
	
	//fprintf(stderr,"\n%f ",silhouette(swarm)); //show the average silhouette
	
	
	/* ----------------END OF THE MOSTRADOR FUNCTION---------------------*/
	
	// 	//print the informations of the particles
	// 	for(int p=0;p<swarms[0].getSize();p++){
	// 		printf("\n\nparticle %d: ", p);
	// 		printVector(swarms[0].particles[p].solution.objectiveVector, objectiveNumber);
	// 		printf(" -> ");
	// 		printVector(swarms[0].particles[p].solution.decisionVector, decisionNumber);
	// 		printf(" (vel: ");
	// 		printVector(swarms[0].particles[p].velocity, decisionNumber);
	// 		printf(")\n p_leader:   ");
	// 		printVector(swarms[0].particles[p].localBest.objectiveVector, objectiveNumber);
	// 		printf(" -> ");
	// 		printVector(swarms[0].particles[p].localBest.decisionVector, decisionNumber);
	// 		printf("\n g_leader:   ");
	// 		printVector(swarms[0].particles[p].globalBest.objectiveVector, objectiveNumber);
	// 		printf(" -> ");
	// 		printVector(swarms[0].particles[p].globalBest.decisionVector, decisionNumber);
	// 	}
	// 		printf("\n--------------------\n");
	// 	//END of information printing block
	
	// 	//print the solutions in all the repositories
	// 	if(iteration==199)
	// 	for(int s=0;s<swarmNumber;s++){
	// 		for(int p=0;p<swarms[s].repository.getActualSize();p++){
	// // 			printf("\n rep_%d: ",p);
	// 			printVector(swarms[s].repository.getSolution(p).objectiveVector, objectiveNumber);
	// 			printf("\n");
	// // 			printf(" -> ");
	// // 			printVector(swarms[s].repository.getSolution(p).decisionVector, decisionNumber);
	// 		}
	// 		printf("\n\n\n");
	// 	}
	// 	//END of repository showing block
	/*
	 *	//print the PBI of the best solution in the repository
	 *	printf("\n");
	 *	for(int s=0;s<swarmNumber;s++){
	 *		double minPbi=MAXDOUBLE;
	 *		int bestSol=-1;
	 *		double solNorm[objectiveNumber];
	 *		for(int p=0;p<swarms[s].repository.getActualSize();p++){
	 *			for(int o=0;o<objectiveNumber;o++)
	 * // 				solNorm[o]=normalize(swarms[s].repository.getSolution(p).objectiveVector[o],globalSmallerObjs[o],globalLargerObjs[o]);
	 *				solNorm[o]=swarms[s].repository.getSolution(p).objectiveVector[o];
	 *			
	 *			double pbi=PBI(solNorm, swarms[s].repository.weight);
	 *			if(minPbi > pbi){
	 *				minPbi=pbi;
	 *				bestSol=p;
}
}
printf("%f(%f,%f) sigma: %f",minPbi, solNorm[0], solNorm[1], swarms[s].sigma);

}
//END of best solution PBI showing block*/
	
	// 	//shows the number of solutions in the repository until it is full, once is full show the iteration number
	// 	if(swarms[0].repository.archiverUsed){
	// 		printf("\nit: %d\n", iteration);
	// 		exit(0);
	// 	}else
	// 		printf("\n%d",swarms[0].repository.getActualSize());
	//	//END of showing the number of solutions in the repository until it is full
	
	//	//Shows all the particles of the first six iterations
	// 	if(iteration < 6){
	// 		for(int p=0;p<swarms[0].getSize();p++){
	// 			printf("\nparticle %d: ", p);
	// 			printVector(swarms[0].particles[p].solution.objectiveVector, objectiveNumber);
	// 			printf(" -> ");
	// 			printVector(swarms[0].particles[p].solution.decisionVector, decisionNumber);
	// 		}
	// 		printf("\n\n");
	// 	}
	//	//END of showing the all the particles for the first six iterations
	
	// 	//shows how many particles entered in each front at each iteration
	// 	double tmp_partEnter[swarmNumber];
	// 	for(int i=0;i<swarmNumber;i++){
	// 		tmp_partEnter[i]=swarms[i].repository.particlesEntered;
	// 		printf("%d ", swarms[i].repository.particlesEntered);
	// 	}
	// 	printf("\n");
	// // 	printVectorToFile(tmp_partEnter, swarmNumber, _partEnter);
	// 	//END of the "particles entered" showing block
	
	// 	//calculates the largest positions found considering all swarms
	// 	double max[objectiveNumber];
	// 	for(int i=0;i<swarmNumber;i++)
	// 		for(int o=0;o<objectiveNumber;o++)
	// 			max[o]=std::max(max[o],swarms[i].repository.largerObjs[o]);
	// 	printf("\n max positions found: ");
	// 	printVector(max, objectiveNumber);
	// 	printf("\n");
	// 	//END of largest positions calculation
	
	if(printHV){
		bool hvPerFront=false; //true= shows one hv per swarm. false= combines all the nondominated solutions from all the swarms and shows only one hv
		//shows the hypervolume per front
		double point=2100;//temporary reference point for now
		point=log10(1+ (point*1000) );//transformation of point
		
		double tmp_hv[swarmNumber];
		char refPoints[100]={0};
		for(int o=0;o<objectiveNumber;o++)
			sprintf(refPoints, "%s%f ", refPoints, point);
		char tmpFile[1000]={0};
		sprintf(tmpFile, "/tmp/hv_tmp(%d).txt",getpid());
		char command[1000]={0};
		sprintf(command, "./assessment/metrics/hv/wfg '%s' %s | head -1 | cut -d '=' -f 2",tmpFile, refPoints);
		for(int i=0;i<swarmNumber;i++){
// 			Solution* sols;
// 			int solNumber=0;
// 			if(hvPerFront){
// 				sols=swarms[i].repository.getSolutions();
// 				solNumber=swarms[i].repository.getActualSize();
// 			}else{
// 				repTemp->initialize(archiver, repositorySize*swarmNumber);
// 				for(int s=0;s<swarmNumber;s++)
// 					repTemp->add(swarms[s].repository);
// // 				repTemp->organize();
// 				sols=repTemp->getSolutions();
// 				solNumber=repTemp->getActualSize();
// 				// 				printf("\ncombined sols: %d\t", solNumber);
// 			}
			
			//printf("%d ", swarms[i].repository.particlesEntered);
			char resCommand[100]={0};
			clearFile(tmpFile); //create an empty temp file
			printToFile(tmpFile, "#");
			double hvSign=1;
// 			for(int p=0;p<solNumber;p++){
// 				for(int o=0;o<objectiveNumber;o++){//security check
// 					// 					if(sols[p].objectiveVector[o] > 3){ //always correct
// 					// 					printf("Corrected from:\t%f to \t%f, point:%f\n",sols[p].objectiveVector[o],log10(1+ (sols[p].objectiveVector[o]*1000)), point );
// 					sols[p].objectiveVector[o]=log10(1+ (sols[p].objectiveVector[o]*1000) );
// 					// 						3.0+log10(sols[p].objectiveVector[o]);
// 					// 					}					
// 					if(point < sols[p].objectiveVector[o]){
// 						fprintf(stderr,"\nError! The reference point of the hypervolume is too small (%f, %f)\n", point, sols[p].objectiveVector[o]);
// 						// 						sols[p].objectiveVector[o]=point;
// 						// 						hvSign=0;
// 						exit(1);
// 					}
// 				}
// 				printVectorToFile(sols[p].objectiveVector, objectiveNumber, tmpFile);
// 			}
// 			repGlobal->organize();
			for(int p=0;p<repGlobal->getActualSize();p++)
				printVectorToFile(repGlobal->getSolution(p).objectiveVector, objectiveNumber, tmpFile);
			
			printToFile(tmpFile, "#");
			exec(command, resCommand);
			
			tmp_hv[i]=atof(resCommand)*hvSign;
			
			if(tmp_hv[i] == 0)
				fprintf(stderr,"\n%f -- %s\n", tmp_hv[i], resCommand);
			
			//printf("%f ", atof(resCommand));
			
			if(diversityPhase || !hvPerFront)
				break;
		}
		if(diversityPhase || !hvPerFront)
			printVectorToFile(tmp_hv, 1, _hv);
		else
			printVectorToFile(tmp_hv, swarmNumber, _hv);
		//printf("\n");
		//END of the hypervolume showing block
	}
	
// 	repGlobal->organize();
// 	double value=calculateHypervolume(*repGlobal, 0);
// 	char charValue[100]={0};
// 	sprintf(charValue, "%f", value);
// 
// 	printToFile(_hv, charValue);
	
	// 	//shows the IGD_p per front
	// 	char tmp[5];
	// 	strcpy(tmp, problem);
	// 	toUpperCase(tmp);
	// 	sprintf(arquivo, "assessment/metrics/pareto/%s_%d",tmp,objectiveNumber);
	// 	tamRef=lerArquivos(reference, arquivo);
	// 	
	// 	for(int s=0;s<swarmNumber;s++){
	// 		double gdp = 0;
	// 		double soma = 0;
	// 		//Percorre todos os pontos do conjunto de aproximacao
	// 		for(int i = 0; i<tamRef; i++){
	// 			double menor_distancia = smallestEuclideanDistance(reference[i], swarms[s].repository.getSolutions(), swarms[s].repository.getActualSize()); 
	// 			soma+=menor_distancia*menor_distancia;
	// 		}
	// 		gdp = sqrt( (1.0/tamRef) * soma ); //gdp
	// 		printf("%f ",gdp);
	// 	}
	// 	printf("\n");
	// 	//END of the IGD_p per front showing block
	
	//shows the IGD_p general, considering all fronts
	if(printIGD){
		repTemp->initialize(archiver, repositorySize*swarmNumber);
		for(int s=0;s<swarmNumber;s++){
			repTemp->add(swarms[s].repository);
		}
// 		repTemp->organize();
		fprintf(stderr,"\n%f ",calculateIGDp(*repTemp));
		
		
		
		// 		char tmp[5];
		// 		strcpy(tmp, problem);
		// 		toUpperCase(tmp);
		// 		sprintf(arquivo, "assessment/metrics/pareto/%s_%d",tmp,objectiveNumber);
		// 		tamRef=lerArquivos(reference, arquivo);
		// 
		// 		repTemp.initialize(archiver, repositorySize*swarmNumber);
		// 		for(int s=0;s<swarmNumber;s++){
		// 			repTemp.add(swarms[s].repository);
		// 		}
		// 		repTemp.organize();
		// 		double gdp[1]; gdp[0] = 0;
		// 		double soma = 0;
		// 		//Percorre todos os pontos do conjunto de aproximacao
		// 		for(int i = 0; i<tamRef; i++){
		// 			double menor_distancia = smallestEuclideanDistance(reference[i], repTemp.getSolutions(), repTemp.getActualSize()); 
		// 			soma+=menor_distancia*menor_distancia;
		// 		}
		// 		gdp[0] = sqrt( (1.0/tamRef) * soma ); //gdp
		// 		//printf("%f ",gdp);
		// 		//printf("\n");
		// 		printVectorToFile(gdp, 1, _igd);
		
	}
	// 	//END of the IGD_p general, considering all fronts showing block
	
	// 	//shows the GD_p general, considering all fronts
	// 	char tmp[5];
	// 	strcpy(tmp, problem);
	// 	toUpperCase(tmp);
	// 	sprintf(arquivo, "assessment/metrics/pareto/%s_%d",tmp,objectiveNumber);
	// 	tamRef=lerArquivos(reference, arquivo);
	// 
	// 	repTemp.initialize();
	// 	for(int s=0;s<swarmNumber;s++){
	// 		repTemp.add(swarms[s].repository);
	// 	}
	// 	repTemp.organize();
	// 	double igdp = 0;
	// 	double soma = 0;
	// 	//Percorre todas as solucoes
	// 	for(int i = 0; i<repTemp.getActualSize(); i++){
	// 		double menor_distancia = smallestEuclideanDistance(repTemp.getSolution(i).objectiveVector, reference, tamRef);
	// 		soma+=menor_distancia*menor_distancia;
	// 	}
	// 	igdp  = sqrt( (1.0/tamRef) * soma ); //igdp 
	// 	printf("%f ",igdp );
	// 	printf("\n");
	// 	//END of the GD_p general, considering all fronts showing block

	
// 	//Show the total number of solutions updated considering all swarms and the size of the global repository
// 	printf("size: %d totalUpdates: %d \n",repGlobal->getActualSize(), totalUpdates);
// // 	printf("%d: replacements: %d bestFit: %f\n", iteration, totalUpdates, ( (swarms[0].repository.getSolution(0).objectiveVector[0]*0.5) + swarms[0].repository.getSolution(0).objectiveVector[1]*0.5) );
// 	totalUpdates=0;
// 	//END of the solutions updated and global repository size showing block
// 	
// 	//shows the avg inner distance of the entire population (considering all swarms)
// 	double tmp_obj=0;
// 	double tmp_dec=0;
// 	int cont=0;
// 
// 	for(int s=0;s<swarmNumber;s++){
// 		for(int p=0;p<swarms[s].getSize();p++){
// 			for(int q=p+1;q<swarms[s].getSize();q++){
// 				cont++;
// 				tmp_obj+=calculateEuclideanDistance(swarms[s].particles[p].solution.objectiveVector, swarms[s].particles[q].solution.objectiveVector, objectiveNumber);
// 				tmp_dec+=calculateEuclideanDistance(swarms[s].particles[p].solution.decisionVector, swarms[s].particles[q].solution.decisionVector, decisionNumber);
// 			}
// 		}
// 	}
// 	printf("obj: %f ",tmp_obj/cont);
// 	printf("dec: %f \n",tmp_dec/cont);
// 	fflush(stdout);
// 	//END of the global avg inner distance showing block
	
// 	double sum=0;
// 	for(int o=0;o<objectiveNumber;o++)
// 		sum+=(swarms[0].repository.weight[o]*swarms[0].repository.getSolution(0).objectiveVector[o]);
// 	
// 	printf("WSum: %f\n", sum);
	
	
	
	
	
	
	// 	//shows the avg inner distance per swarm
	// 	double tmp_obj[swarmNumber]={0};
	// 	double tmp_dec[swarmNumber]={0};
	// 
	// 	for(int s=0;s<swarmNumber;s++){
	// 		int cont=0;
	// 		for(int p=0;p<swarms[s].repository.getActualSize();p++){
	// 			for(int q=p+1;q<swarms[s].repository.getActualSize();q++){
	// 				cont++;
	// 				tmp_obj[s]+=distanciaEuclidiana(swarms[s].repository.getSolution(p).objectiveVector, swarms[s].repository.getSolution(q).objectiveVector, objectiveNumber);
	// 				tmp_dec[s]+=distanciaEuclidiana(swarms[s].repository.getSolution(p).decisionVector, swarms[s].repository.getSolution(q).decisionVector, decisionNumber);
	// 			}
	// 		}
	// 		//printf("%f ",tmp_obj[s]/cont);
	// 		//printf("%f ",tmp_dec[s]/cont);
	// 		if(cont>0){
	// 			tmp_obj[s]/=cont;
	// 			tmp_dec[s]/=cont;
	// 		}
	// 	}
	// 	//printf("\n");
	// // 	printVectorToFile(tmp_obj, swarmNumber, _obj);
	// // 	printVectorToFile(tmp_dec, swarmNumber, _dec);
	// 	//END of the avg inner distance showing block
	
	// 	//Shows the clustering quality at each iteration
	// 	repTemp.initialize(archiver, repositorySize*swarmNumber);
	// 	for(int s=0;s<swarmNumber;s++){
	// 		repTemp.add(swarms[s].repository);
	// 	}
	// 	repTemp.organize();
	// 	fprintf(stderr,"\n%f ",calculateIGDp(repTemp));
	// 
	// 	if(swarmNumber > 1)
	// 		fprintf(stderr,"%f %f %f ", silhouette(swarm), dunnIndex(swarm), daviesBouldinIndex(swarm)); //higher, higher, lower
	// 	//END of showing the clustering quality block
	
	
	// 	//shows the number of solutions in each repository
	// 	printf("(%d) ",iteration);
	// 	for(int s=0;s<swarmNumber;s++)
	// 		//printf("swarm %d - %d\n",s,swarms[s].repository.getActualSize());
	// // 		printf("%d ", swarms[s].repository.getActualSize());
	// // 		printf("%d ", swarms[s].neighborhoodSize);
	// 	printf("\n");
	// 	//END of showing the number of solutions in each repository block
	
	// 	//shows the weight vector of each repository
	// 	for(int s=0;s<swarmNumber;s++){
	// 		//printf("swarm %d - %d\n",s,swarms[s].repository.getActualSize());
	// 		printVector(swarms[s].repository.weight, objectiveNumber);
	// 		printf("\n");
	// 	}
	// 	printf("\n\n\n");
	// 	//END of showing weight vector of each repository block
	
	
	// 	//shows the number of particles in each swarm, and the sum (considering all swarms)
	// 	int sum=0;
	// 	for(int s=0;s<swarmNumber;s++){
	// 		printf("%d ",swarms[s].getSize());
	// 		sum+=swarms[s].getSize();
	// 	}
	// 	printf("--- %d\n",sum);
	// 	//END of showing the number of particles in each swarm
	
	
	// 	//print the solutions when reach x iterations
	// 	if(iteration == 180){
	// 		for(int i=0;i<swarmNumber;i++){
	// 			for(int p=0;p<swarms[i].repository.getActualSize();p++){
	// 				printVector(swarms[i].repository.getSolution(p).objectiveVector, objectiveNumber);
	// 				printf("\n");
	// 			}
	// 			printf("\n\n\n");
	// 		}
	// 	}
	// 	//END of printing the solutions block
	
	// 	//printing the iteration, decision and objective vectors to files as Roberto requested
	// 	for(int i=0;i<swarmNumber;i++){
	// 		for(int p=0;p<swarms[i].repository.getActualSize();p++){
	// 			sprintf(line, "%d", iteration);
	// 			
	// 			for(int d=0;d<decisionNumber;d++)
	// 				sprintf(line, "%s %f", line, swarms[i].repository.getSolution(p).decisionVector[d]);
	// 			for(int o=0;o<objectiveNumber;o++)
	// 				sprintf(line, "%s %f", line, swarms[i].repository.getSolution(p).objectiveVector[o]);
	// 			
	// 			sprintf(line, "%s", line);
	// 			//printToFile(_solObj, line);
	// 			printf("%s\n", line);
	// 		}
	// 		printf("\n\n\n");
	// 	}
	// 	//END of printing the iteration, decision and objective vectors
	
	
	// 	//printing the covariance matrix to a separate file as Roberto requested
	// 	if(!diversityPhase){
	// 		for(int s=0;s<swarmNumber;s++){
	// 			printf("\n------------------Final covariance matrix-------------------");
	// 			for(int i=0;i<decisionNumber;i++){
	// 				printf("\n");
	// 				printVector(swarms[s].C[i], decisionNumber);
	// 				//printVectorToFile(swarms[s].C[i], decisionNumber, _matCov);
	// 			}
	// 			//insertBlankLine(_matCov);
	// 			printf("\n-------------------------------------\n");
	// 		}
	// 	}
	// 	//END of printing the covariance matrix
	
	// 	//printing on the screen the objective values vs. the log-likelihood of each solution vs. the euclidean distance from the solution to the mean in the decision space
	// 	printf("\n");
	// 	for(int i=0;i<swarmNumber;i++){
	// 		for(int p=0;p<swarms[i].repository.getActualSize();p++){
	// 			for(int o=0;o<objectiveNumber;o++)
	// 				printf("%f ", swarms[i].repository.getSolution(p).objectiveVector[o]);
	// 			double likelihood_log=logLikelihood(swarms[i].repository.getSolution(p).decisionVector, swarms[i].mean, swarms[i].C, swarms[i].B, swarms[i].D, swarms[i].det);
	// 			printf(" -- %f -> %e -- %f \n", likelihood_log, exp(likelihood_log), calculateEuclideanDistance(swarms[i].repository.getSolution(p).decisionVector, swarms[i].mean, decisionNumber));
	// 		}
	// 	}
	// 	//END of printing the log-likelihood block
}

void algorithmInitializationOperations(const int argc, const char* argv[]){
	unsigned long seed=clock()+time(NULL);
	char run[5]={0};
	if(argc > 2){
		seed=( atof(argv[2]) );
		sprintf(run, ".%d",(int)seed);
	}
	srand(seed);
	
	sprintf(_s, "%s_solutions%s.txt",outputFileName,run);
	sprintf(_f, "%s_fronts%s.txt",outputFileName,run);
	clearFile(_s); //clear the solutions file, so it previous content is deleted
	clearFile(_f); //clear the front file, so it previous content is deleted
	if(printHV){
		sprintf(_hv, "%s_hv%s.txt",outputFileName,run);
		clearFile(_hv); //clear the front file, so it previous content is deleted
	}
	if(printIGD){
		sprintf(_igd, "%s_igd%s.txt",outputFileName,run);
		clearFile(_igd); //clear the front file, so it previous content is deleted
	}
	if(printMOEAD_CMAES){
		sprintf(_MOEAD_CMAES, "%s_moead-cmaes%s.txt",outputFileName,run);
		clearFile(_MOEAD_CMAES); //clear the file, so it previous content is deleted
	}
	if(!strcmp(algorithm, "hmopso")){
		sprintf(_p, "%s_probabilities%s.txt",outputFileName,run);
		clearFile(_p); //clear the front file, so it previous content is deleted
	}
		
	int fileSize=1e4; //3100 should be enough
	reference = new double*[fileSize];
	for(int i=0; i<fileSize;i++)
		reference[i] = new double[objectiveNumber];
	
	if(!strcmp(algorithm, "hmopso") || !strcmp(archSubSwarms, "r2") || !strcmp(archiver, "r2") || !strcmp(archSubSwarms, "pbi") || !strcmp(archSubSwarms, "tch")  || !strcmp(archSubSwarms, "wcp") || !strcmp(archSubSwarms, "wsum") || !strcmp(archSubSwarms, "r-ideal") || !strcmp(cmaRanking, "r2")){
		sprintf(file, "weights/W_%dD",objectiveNumber);
		refSize=readFile(reference, file);
	}
	
	if(!strcmp(archSubSwarms, "pbi")  || !strcmp(archSubSwarms, "tch") || !strcmp(archSubSwarms, "wcp") || !strcmp(archSubSwarms, "wsum") || !strcmp(archSubSwarms, "r-ideal")){
		if(!strcmp(algorithm, "imulti") || !strcmp(algorithm, "cmulti") || !strcmp(algorithm, "cmaes-mopso") || !strcmp(algorithm, "moead") ){
			// 			decomposition=true;
			originalSwarmNumber=swarmNumber=refSize;
			subSwarmPopulation=(subSwarmTotalPopulation/swarmNumber);
			globalTmpWeight = new double[objectiveNumber];//temporary weight vector used in the cobyla local optimizer
		}
	}
}

void finalizeRun(timeval startTime){
	// 	Merge all the nondominated solutions into a single front
	if(decomposition && useGlobalRepOnDecomposition){
		// 		repTemp->initialize(archiver, repGlobal->getActualSize());
		// 		repTemp->add(*repGlobal);
		delete repTemp;
		repTemp=repGlobal;
		for(int s=0;s<swarmNumber;s++){
			for(int sl=0;sl<swarms[s].repository.getActualSize();sl++)
				repTemp->add(swarms[s].repository.getSolutions()[sl]);
		}
	}else{
		if(swarmNumber == 1){
			delete repTemp;
			repTemp=&swarms[0].repository;
		}else{
			repTemp->initialize(archiver, swarmNumber*repositorySize);
			for(int s=0;s<swarmNumber;s++){
				repTemp->add(swarms[s].repository);
			}
		}
	}
// 	repTemp->organize();
	std::sort(repTemp->getSolutions(), repTemp->getSolutions()+repTemp->getActualSize(), compareSolutions(0)); //sort the solutions according to the first objectives, just a cosmetic step
	// END of merging the solutions
	
	//print the solutions found
	for(int p=0;p<repTemp->getActualSize();p++){
		//printVector(front[p].decisionVector, decisionNumber);
		//printf(" -> ");
		//printVector(repUnico.getSolution(p).objectiveVector, objectiveNumber);
		//printf("\n");
		printVectorToFile(repTemp->getSolution(p).decisionVector, decisionNumber, _s);
		printVectorToFile(repTemp->getSolution(p).objectiveVector, objectiveNumber, _f);
	}// END of printing solutions
	
	insertBlankLine(_s);
	insertBlankLine(_f);
	if(printHV){
		insertBlankLine(_hv);
	}
	if(!strcmp(algorithm, "hmopso"))
		insertBlankLine(_p);
	gettimeofday(&endTime, NULL); 
	double diff= ((double)endTime.tv_sec * 1e6 + endTime.tv_usec - ((double)startTime.tv_sec * 1e6 + startTime.tv_usec))/1e6;
	fprintf(stderr, "\nRun in %02.0f:%02.2f s, found %d solutions using %ld func. evals.\n", (diff/60.0), fmod(diff, 60), repTemp->getActualSize(), funcEvals);
	
	//Clear structures	
	if(swarmNumber == 1 || (decomposition && useGlobalRepOnDecomposition) )//in these cases, repTemp only points to other structures that will be freed after
		repTemp=NULL;
	
	if(swarms != NULL){ //if there is already memory allocated, free it
		delete[] swarms;
		swarms=NULL;
	}
	delete repGlobal;
	repGlobal=NULL;
// 	if(!decomposition)//if decomposition it is already freed
	delete repTemp;
	repTemp=NULL;
}

void readParameters(const int argc, const char *inputFile){
	if(argc < 2){
		fprintf(stderr,"uso: mopso <parameters file>\n\n");
		exit(1);
	}
	
	char s[500];
	FILE *stream = fopen(inputFile, "r");
	if( stream == NULL ){
		fprintf(stderr,"\nFail to Open Parameters File!!\n");
		exit(1);
	}
	while(fgets(s,500,stream)){
		if ((s[0]) && (s[0] != '%')) { //if the line is not a comment
			
			if(strlen(s)-1>0)
				s[strlen(s)-1]='\0';
			
			char* var = strtok(s,"=");
			char* value = strtok(NULL,"=");
			
			if(!strncmp(var, "outName", 7))
				strcpy(outputFileName, value);
			if(!strncmp(var, "problem", 7))
				strcpy(problemName, value);
			if(!strncmp(var, "leader", 6))
				strcpy(leader, value);
			if(!strncmp(var, "archiver", 8))
				strcpy(archiver, value);
			if(!strncmp(var, "objectiveNumber", 15))
				objectiveNumber=atoi(value);
			if(!strncmp(var, "population", 10))
				originalPopSize=(int)atof(value);//atoi does not recognize scientific notation, so we use atof and cast to int in some cases
			if(!strncmp(var, "repository", 10))
				repositorySize=(int)atof(value);
			if(!strncmp(var, "swarms", 6))
				originalSwarmNumber=swarmNumber=atoi(value);
			if(!strncmp(var, "iterations", 10))
				maxIterations=(long)atof(value);
			if(!strncmp(var, "runs", 4))
				numExec=atoi(value);
			if(!strncmp(var, "algorithm", 9))
				strcpy(algorithm, value);
			if(!strncmp(var, "partIterations", 14))
				numParticoes=atoi(value);
			if(!strncmp(var, "archSubSwarms", 13))
				strcpy(archSubSwarms, value);
			if(!strncmp(var, "subSwarmTotalPopulation", 23))
				subSwarmTotalPopulation=(int)atof(value);
			if(!strncmp(var, "diversityIterations", 19))
				diversityIterations=atoi(value);
			if(!strncmp(var, "truncType", 9))
				strcpy(truncType, value);
			if(!strncmp(var, "clusteringType", 14))
				strcpy(clusteringType, value);
			if(!strncmp(var, "clusteringDistance", 18))
				strcpy(clusteringDistance, value);
			if(!strncmp(var, "weightDistribution", 18))
				strcpy(weightDistribution, value);
			if(!strncmp(var, "solSet", 6))
				strcpy(solSet, value);
			if(!strncmp(var, "cmaRanking", 10))
				strcpy(cmaRanking, value);
			if(!strncmp(var, "neighborhoodSize", 16))
				globalNeighborhoodSize=atoi(value);
			if(!strncmp(var, "updateNeighborhoodMetric", 24))
				updateNeighborhoodMetric=atoi(value);
		}
	}
	fclose(stream);
	
	if(strcmp(algorithm, "smpso") && strcmp(algorithm, "hmopso") && strcmp(algorithm, "imulti") && strcmp(algorithm, "cmulti") && strcmp(algorithm, "cmaes-mopso") && strcmp(algorithm, "moead") && strcmp(algorithm, "pareto-cmaes")){
		fprintf(stderr,"\nError! Invalid algorithm\n");
		exit(1);
	}
	
	if(strcmp(problemName, "dtlz1") && strcmp(problemName, "dtlz2") && strcmp(problemName, "dtlz3") && strcmp(problemName, "dtlz4") && 
	strcmp(problemName, "dtlz5") && strcmp(problemName, "dtlz6") && strcmp(problemName, "dtlz7") && strcmp(problemName, "wfg1") && 
	strcmp(problemName, "wfg2") && strcmp(problemName, "wfg3") && strcmp(problemName, "wfg4") && strcmp(problemName, "wfg5") &&
	strcmp(problemName, "wfg6") && strcmp(problemName, "wfg7") && strcmp(problemName, "wfg8") && strcmp(problemName, "wfg9") && 
	strcmp(problemName, "uf1") && strcmp(problemName, "uf2") && strcmp(problemName, "uf3") && strcmp(problemName, "uf4") && 
	strcmp(problemName, "uf5") && strcmp(problemName, "uf6") && strcmp(problemName, "uf7") && strcmp(problemName, "uf8") && 
	strcmp(problemName, "uf9") && strcmp(problemName, "uf10") && strncmp(problemName, "tandem",6) && strncmp(problemName, "coco",4)){
		printf("INVALID OPTIMIZATION PROBLEM! (%s)\n", problemName);
		exit(1);
	}
	
	if(!strncmp(problemName, "dtlz", 4)){
		if(!strcmp(problemName, "dtlz1"))
			decisionNumber=objectiveNumber+5-1; //if dtlz1, k=5
			else
				if(!strcmp(problemName, "dtlz7"))
					decisionNumber=objectiveNumber+20-1; //if dtlz7, k=20
				else
					decisionNumber=objectiveNumber+10-1; // otherwise, k=10
						
					inferiorPositionLimit = new double[decisionNumber];
					superiorPositionLimit = new double[decisionNumber];
				
				for (int var = 0; var < decisionNumber; var++) {
					inferiorPositionLimit[var] = 0;
					superiorPositionLimit[var] = 1;
				}
	}
	// 	decisionNumber=3;
	
	if(!strncmp(problemName, "tandem", 6)){//if the problem is from the TANDEM trajectory instances
		int problemId=-1;
		if(strlen(problemName) > 6 && strlen(problemName) < 9){
			char tmp[3]="";
			strncpy(tmp, problemName+6, strlen(problemName)-6);
			problemId=atoi(tmp);
		}else{
			fprintf(stderr,"\nERROR! on problem ID of tandem (%s) (%d)\n", problemName, problemId);
			exit(1);
		}
		switch(problemId){	//best values unconstrained	constrained -- the second objective for constrained problems should be smaller than 3652.5 (3652.5 days == 10 years)
			case 1: {int ar[]={3,2,2,2,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVVVS	1233.49	1087.66
			case 2: {int ar[]={3,2,2,3,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVVES	1412.28	1338.72
			case 3: {int ar[]={3,2,2,4,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVVMS	772.53	664.76
			case 4: {int ar[]={3,2,2,5,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVVJS	976.18	651.55
			case 5: {int ar[]={3,2,3,2,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVEVS	1142.6	938.06
			case 6: {int ar[]={3,2,3,3,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVEES	1673.88	1500.46
			case 7: {int ar[]={3,2,3,4,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVEMS	1163.6	814.51
			case 8: {int ar[]={3,2,3,5,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVEJS	1603.4	740.6
			case 9: {int ar[]={3,2,4,2,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVMVS	657.84	589.03
			case 10: {int ar[]={3,2,4,3,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVMES	1104.51	859.2
			case 11: {int ar[]={3,2,4,4,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVMMS	471.99	336.62
			case 12: {int ar[]={3,2,4,5,6};memcpy(seq,ar,sizeof(int)*5);}break;//EVMJS	603.76	272.13
			case 13: {int ar[]={3,3,2,2,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEVVS	764.74	601.12
			case 14: {int ar[]={3,3,2,3,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEVES	945.4	806.48
			case 15: {int ar[]={3,3,2,4,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEVMS	726.81	460
			case 16: {int ar[]={3,3,2,5,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEVJS	836.46	428.06
			case 17: {int ar[]={3,3,3,2,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEEVS	896.96	558.51
			case 18: {int ar[]={3,3,3,3,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEEES	1242.61	969.96
			case 19: {int ar[]={3,3,3,4,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEEMS	1192.62	506.27
			case 20: {int ar[]={3,3,3,5,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEEJS	1351.53	700
			case 21: {int ar[]={3,3,4,2,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEMVS	812.216	502.03
			case 22: {int ar[]={3,3,4,3,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEMES	1265.44	928.57
			case 23: {int ar[]={3,3,4,4,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEMMS	1077.95	459.33
			case 24: {int ar[]={3,3,4,5,6};memcpy(seq,ar,sizeof(int)*5);}break;//EEMJS	1209.26	281.36
			
			default:
				fprintf(stderr,"\nERROR! on problem ID of tandem (%s) (%d)\n", problemName, problemId);
				exit(1);
		}
		objectiveNumber=2;
		decisionNumber=18;
		inferiorPositionLimit = new double[decisionNumber];
		superiorPositionLimit = new double[decisionNumber];
		//  HERE WE SHOULD SET THE POSITION LIMITS ACCORDINGLY
		for(int i=0;i<decisionNumber;i++){
			inferiorPositionLimit[i] = 0;
			superiorPositionLimit[i] = 1;
		}
	}
	
	if(!strncmp(problemName, "wfg", 3)){//if the problem is from the WFG family
		//l=20 (distance related), k=4 (position related) if M=2, otherwise k=2*(M-1) //from the WFG readme file
		if(objectiveNumber == 2)
			decisionNumber=20+4;
		else
			decisionNumber=20+(2*(objectiveNumber-1));
		
		inferiorPositionLimit = new double[decisionNumber];
		superiorPositionLimit = new double[decisionNumber];
		
		for (int var = 0; var < decisionNumber; var++) {
			inferiorPositionLimit[var] = 0;
			superiorPositionLimit[var] = 2 * (var + 1);
		}
	}
	
	if(!strncmp(problemName, "uf", 2)){//if the problem is from the uf (CEC09) family
		decisionNumber=30;
		inferiorPositionLimit = new double[decisionNumber];
		superiorPositionLimit = new double[decisionNumber];
		inferiorPositionLimit[0]=0;
		superiorPositionLimit[0]=1;
		
		if(!strncmp(problemName, "uf1", 3) || !strncmp(problemName, "uf2", 3) || !strncmp(problemName, "uf5", 3) || !strncmp(problemName, "uf6", 3) || !strncmp(problemName, "uf7", 3)){
			objectiveNumber=2;
			for (int i = 1; i < decisionNumber; i++) {
				inferiorPositionLimit[i]=-1;
				superiorPositionLimit[i]=1;
			}
		}
		
		if(!strncmp(problemName, "uf3", 3)){
			objectiveNumber=2;
			for (int i = 0; i < decisionNumber; i++) {
				inferiorPositionLimit[i]=0;
				superiorPositionLimit[i]=1;
			}
		}
		
		if(!strncmp(problemName, "uf4", 3)){
			objectiveNumber=2;
			for (int i = 1; i < decisionNumber; i++) {
				inferiorPositionLimit[i]=-2;
				superiorPositionLimit[i]=2;
			}
		}
		
		if(!strncmp(problemName, "uf8", 3) || !strncmp(problemName, "uf9", 3) || !strncmp(problemName, "uf10", 4)){
			objectiveNumber=3;
			inferiorPositionLimit[1]=0;
			superiorPositionLimit[1]=1;
			for (int i = 2; i < decisionNumber; i++) {
				inferiorPositionLimit[i]=-2;
				superiorPositionLimit[i]=2;
			}
		}
	}
	
	if(!strncmp(problemName, "coco", 4)){//if the problem is from the coco benchmark
		//initialize with anything because it will be initialized again inside their framework
		decisionNumber=-1;
		objectiveNumber=2;//only bi-objective problems for now
	}
}