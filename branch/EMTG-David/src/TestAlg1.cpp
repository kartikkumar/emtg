#include <ctime>

#include "TestAlg1.h"
#include "EMTG_math.h"

#include <iostream>

namespace EMTG { namespace Solvers {

	//default constructor (never called)
	TestAlg1::TestAlg1(void)
	{
	}

	//constructor with a problem object
	TestAlg1::TestAlg1(EMTG::problem* Problem_input)
	{
		Initialize(Problem_input);
	}

	//destructor
	TestAlg1::~TestAlg1(void)
	{
	}

	//method to initialize the solver
	void TestAlg1::Initialize(EMTG::problem* Problem_input)
	{
		//store a pointer to the problem object
		Problem = Problem_input;

		////extract the dimensions of the problem
		nX = Problem->total_number_of_NLP_parameters;
		nF = Problem->total_number_of_constraints;
		NP = 10 * nX;

	}

	

	//main ACDE loop
	void TestAlg1::Evolve()
	{
		for (int i = 0; i < nX; ++i)
		{
			BestX.push_back(0.0);
		}
		for (int i = 0; i < nF; ++i)
		{
			FVector.push_back(0.0);
		}
		cout << "Do Darwin's Bidding!\n";
		Problem->unscale(BestX.data());
		Problem->evaluate(Problem->X.data(), FVector.data(), Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);
		//BestObjectiveValue = 0;
		//Problem->unscale(&BestX[0]);
		Problem->Xopt = Problem->X; //we store the unscaled Xbest
		//Problem->Xopt = BestX;
		//Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
		//Problem->output();
	}

}}//close namespaces