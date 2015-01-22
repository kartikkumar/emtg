#include "problem.h"

#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"

#ifndef _TESTALG1
#define _TESTALG1

namespace EMTG { namespace Solvers {

	class TestAlg1
	{
	public:

		//default constructor
		TestAlg1(void);

		//constructor with a problem object
		TestAlg1(EMTG::problem* Problem_input);

		//destructor
		virtual ~TestAlg1(void);

		//methods
		void Initialize(EMTG::problem* Problem_input);
		//void GenerateNewPopulation();
		//void EvaluateIndividual(vector<double>& X, double* ObjectiveFunctionValue, double* ConstraintNormValue);
		void Evolve();
		void RandInit(int gens);
		void EvalInd(vector<double>& X, double* ObjectiveFunctionValue, double* ConstraintNormValue);
		void LocOptInd(vector<double>& X);
		void SortBin(vector<vector<double>>& bin,vector<vector<double>>& fit,vector<double>& X,vector<double>& score,int cap,int ind);
		std::vector<vector<vector<double>>> Bins;
		std::vector<vector<vector<double>>> Fits;
		std::vector<double> vars;
		int iter;


		//fields
		int NP;
		int nX;
		int nF;
		//int AgeLimit;
		//double Tau1; //adaptation parameter from jDE
		//double Tau2; //adaptation parameter from jDE
		//double ScaleFactorLowerBound; //adaptation parameter from jDE
		//double ScaleFactorUpperBound; //adaptation parameter from jDE
		//double ScaleFactor;
		//double CrossoverRatio;

		//vector< vector<double> > CurrentGenerationPopulation;
		//vector< vector<double> > NextGenerationPopulation;
		//vector<double> CurrentGenerationObjectiveValue;
		//vector<double> NextGenerationObjectiveValue;
		//vector<double> CurrentGenerationConstraintNorm;
		//vector<double> NextGenerationConstraintNorm;
		std::vector<double> FVector;
		//vector<double> ConstraintViolationVector;
		//vector<double> TrialVector;
		//vector<int> CurrentGenerationAges;
		//vector<int> NextGenerationAges;
		//double TrialObjectiveValue;
		//double TrialConstrainNorm;

		std::vector<double> BestX;
		std::vector<double> xTemp;
		std::vector<double> ConstraintViolationVector;
		double objValTemp;
		double conValTemp;
		double BestObjectiveValue;
		double BestConstraintNorm;
		//int BestIndex;
		

		//pointer to problem object
		EMTG::problem* Problem;

		//random number generator
		boost::mt19937 RNG;
		//boost::uniform_int<> IntegerDistributionNP;
		boost::uniform_int<> capDist;
		boost::uniform_int<> binDist;
		boost::uniform_real<> DoubleDistribution;
		boost::uniform_real<> VarDist;
	};

}}//close namespaces

#endif //_TESTALG1