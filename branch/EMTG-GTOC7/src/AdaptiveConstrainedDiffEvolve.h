//header file for adaptive constrained differential evolution
//a hybrid of jDE (self-adaptive differential evolution) by Brest et al. "Dynamic Optimization using Self-Adaptive Differential Evolution"
//and CDE (constrained differential evolution) by Mezura-Montes, Coello, and Tun-Morales, "Simple Feasibility Rules and Differential Evolution for Constrained Optimization"
//inspired by ACT solution to GTOC6
//created for EMTG by Jacob Englander, 5/6/2013

#include "problem.h"

#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"

#ifndef _ADAPTIVECONSTRAINEDDIFFEVOLVE
#define _ADAPTIVECONSTRAINEDDIFFEVOLVE

namespace EMTG { namespace Solvers {

	class AdaptiveConstrainedDiffEvolve
	{
	public:

		//default constructor
		AdaptiveConstrainedDiffEvolve(void);

		//constructor with a problem object
		AdaptiveConstrainedDiffEvolve(EMTG::problem* Problem_input);

		//destructor
		virtual ~AdaptiveConstrainedDiffEvolve(void);

		//methods
		void Initialize(EMTG::problem* Problem_input);
		void GenerateNewPopulation();
		void EvaluateIndividual(vector<double>& X, double* ObjectiveFunctionValue, double* ConstraintNormValue);
		void Evolve();

		//fields
		int NP;
		int nX;
		int nF;
		int AgeLimit;
		double Tau1; //adaptation parameter from jDE
		double Tau2; //adaptation parameter from jDE
		double ScaleFactorLowerBound; //adaptation parameter from jDE
		double ScaleFactorUpperBound; //adaptation parameter from jDE
		double ScaleFactor;
		double CrossoverRatio;

		vector< vector<double> > CurrentGenerationPopulation;
		vector< vector<double> > NextGenerationPopulation;
		vector<double> CurrentGenerationObjectiveValue;
		vector<double> NextGenerationObjectiveValue;
		vector<double> CurrentGenerationConstraintNorm;
		vector<double> NextGenerationConstraintNorm;
		vector<double> FVector;
		vector<double> ConstraintViolationVector;
		vector<double> TrialVector;
		vector<int> CurrentGenerationAges;
		vector<int> NextGenerationAges;
		double TrialObjectiveValue;
		double TrialConstrainNorm;

		vector<double> BestX;
		double BestObjectiveValue;
		double BestConstraintNorm;
		int BestIndex;
		

		//pointer to problem object
		EMTG::problem* Problem;

		//random number generator
		boost::mt19937 RNG;
		boost::uniform_int<> IntegerDistributionNP;
		boost::uniform_int<> IntegerDistributionnX;
		boost::uniform_real<> DoubleDistribution;
	};

}}//close namespaces

#endif //_ADAPTIVECONSTRAINEDDIFFEVOLVE