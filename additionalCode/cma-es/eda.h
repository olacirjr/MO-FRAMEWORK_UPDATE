class Repository;
class Solution;
double cmaes_random_Gauss();
static void QLalgo2 (int n, double *d, double *e, double **V);
static void Householder2(int n, double **V, double *d, double *e);
static void Eigen( int N,  double **C, double *diag, double **Q, double *rgtmp);
static double rgdouMax( const double *rgd, int len);
static double rgdouMin( const double *rgd, int len);

//class cmaesModel
class CMAESModel{
public:
	
	CMAESModel();//constructor
	~CMAESModel();//destructor
	CMAESModel(const CMAESModel &source);//copy
	CMAESModel& operator= (const CMAESModel &source);//assignment
	
	void initialize();
	
	void rankMuUpdate(Repository &repository, const bool* injected);
	void rankOneUpdate(const Solution &sol, const bool suc);
	void rankOneOnlySigma(const bool suc);
	bool sample(const double* gaussianMean, double* out);
	
	double logLikelihood(double* sol) const;
	
	double combineMatrices(const CMAESModel &outerModel);
	
// 	const double* getMean(){return mean;}
// 	double** getMatrix(){return C;}
	
	double *prevSol;
	double modelQuality;//stores the quality of the model
	
// private:
	double sigma;
	double cSigma;
	double dSigma;
	double det;
	double pSuc;
	double prevSigma;
	double *pC;
	double *pSigma;
	double *mean;
	double *D;
	double **C;
	double **B;
	bool init;
	int gen;
	
	//gaussian generator
	int flgstored;
	double hold;
}; 
