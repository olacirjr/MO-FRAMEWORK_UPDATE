#include "MOPSO.cpp"
void Problem::calculateCoco(double* decisionVector, double* objectiveVector){}

int main(const int argc, const char* argv[]){
	
	readParameters(argc, argv[1]);
	algorithmInitializationOperations(argc, argv);
	
// 	sprintf(problem, "dtlz2");
	
	for(int r=0;r<numExec;r++){ //independent runs loop
		fprintf(stderr, "----------------Independent run # %d-------------- ", r);
		
// 		decisionNumber+=r+2;
// 		
// 		delete[] inferiorPositionLimit;
// 		delete[] superiorPositionLimit;
// 		inferiorPositionLimit = new double[decisionNumber];
// 		superiorPositionLimit = new double[decisionNumber];
// 		
// 		for (int var = 0; var < decisionNumber; var++) {
// 			inferiorPositionLimit[var] = 0;
// 			superiorPositionLimit[var] = 1;
// 		}
		
// 		maxEvals=decisionNumber*1e5;//setting the same number of evaluations for all algorithms
// 		maxEvals=50000;//setting the same number of evaluations for all algorithms
		
		run();
		
	} //END of the independent runs loop
	
	exec("notify-send 'End run.'", NULL);
	return 0;
}//END of main