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
		NP = 100; // Bin Depth -- Make Alg Param
		BestObjectiveValue = EMTG::math::LARGE;
		BestConstraintNorm = EMTG::math::LARGE;
		ConstraintViolationVector.resize(nF - 1);
		int gens = 200; // Make Alg Param
		FVector.resize(nF);
		xTemp.resize(nX);
		DoubleDistribution = boost::uniform_real<>(0.0, 1.0);
		VarDist = boost::uniform_real<>(-1.0, 1.0);
		capDist = boost::uniform_int<>(0, NP-1);
		binDist = boost::uniform_int<>(0, nF-1);
		RNG.seed(time(NULL));
		Bins.resize(nF);
		Fits.resize(nF);

		for (int i = 0; i < nF; i++)
		{
			Bins[i].resize(NP);
			Fits[i].resize(NP);
			for (int j = 0; j < NP; j++)
			{
				Bins[i][j].resize(nX);
				Fits[i][j].resize(nF);
				for (int k = 0; k < nF; k++)
					Fits[i][j][k] = EMTG::math::LARGE;
			}
		}

		iter = 500; // Make Alg Param
		double varPer = 0.1; // Make Alg Param
		vars.resize(nX);
		for (int i = 0; i < nX; i++)
			vars[i] = (varPer/iter)*(Problem->Xupperbounds[i]-Problem->Xlowerbounds[i]);

		RandInit(gens);
	}

	void TestAlg1::RandInit(int gens)
	{
		for (int i = 0; i < gens; i++)
		{
			for (int j = 0; j < nX; j++)
				xTemp[j] = DoubleDistribution(RNG);

			// HERE LOCALLY OPTIMIZE xTemp---------------------------------------
			LocOptInd(xTemp);
			EvalInd(xTemp, &objValTemp, &conValTemp);
			vector<double> score (nF);
			score[0] = objValTemp + conValTemp;
			for (int j = 0; j < nF - 1; ++j)
			{
				score[j+1] = objValTemp + ConstraintViolationVector[j];
			}
			for (int j = 0; j < nF; j++)
			{
				if(score[j] < Fits[j][NP-1][j])
					SortBin(Bins[j],Fits[j],xTemp,score,NP-1,j);
			}
		}
	}

	void TestAlg1::LocOptInd(vector<double>& X)
	{
		vector<double> Xb = X;
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conValTemp;
		for (int j = 0; j < iter; j++)
		{
			vector<double> Xt (nX);
			for (int i = 0; i < nX; i++)
			{
				double temp = Xb[i] + vars[i] * VarDist(RNG);
				if(temp > Problem->Xupperbounds[i])
					temp = 2*Problem->Xupperbounds[i]-temp;
				if(temp < Problem->Xlowerbounds[i])
					temp = 2*Problem->Xlowerbounds[i]-temp;
				Xt[i] = temp;
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conValTemp;
			if(St < Sb)
			{
				Xb = Xt;
				Sb = St;
			}
		}
		X = Xb;
	}

	void TestAlg1::SortBin(vector<vector<double>>& bin,vector<vector<double>>& fit,vector<double>& X,vector<double>& score,int cap,int ind)
	{
		if(cap>0 && score[ind] < fit[cap-1][ind])
			SortBin(bin,fit,X,score,cap-1,ind);
		else
		{
			vector<vector<double>>::iterator it;
			it = bin.begin();
			it+=cap;
			it = bin.insert (it,X);
			vector<vector<double>>::iterator it2;
			it2 = fit.begin();
			it2+=cap;
			it2 = fit.insert (it2,score);
			//bin.insert(cap,X);
			bin.pop_back();
			fit.pop_back();
			if(cap==0)
				cout << "Bin " << ind << " has a new best at Fit = " << score[ind] << "\n";
		}
	}

	void TestAlg1::EvalInd(vector<double>& X, double* ObjectiveFunctionValue, double* ConstraintNormValue)
	{
		//unscale the decision vector
		Problem->unscale(X.data());

		//set the failure flag to false
		bool failflag = false;

		//evaluate the decision vector
		try
		{
			Problem->evaluate(Problem->X.data(), FVector.data(), &Problem->G[0], 0, Problem->iGfun, Problem->jGvar);
		}
		catch (int errorcode) //integration step error
		{
			if (errorcode == 13)
				failflag = true;
		}

		//return the objective function value
		*ObjectiveFunctionValue = FVector[0];

		//compute the constraint violations
		if (failflag)
		{
			*ConstraintNormValue = 1.0e+6;
		}
		else
		{
			for (int j = 0; j < nF - 1; ++j)
			{
				if (!(FVector[j+1] == FVector[j+1]))
					ConstraintViolationVector[j] = EMTG::math::LARGE; //NaN trap
				else if (FVector[j+1] > Problem->Fupperbounds[j+1])
					ConstraintViolationVector[j] = FVector[j+1] - Problem->Fupperbounds[j+1];
				else if (FVector[j+1] < Problem->Flowerbounds[j+1])
					ConstraintViolationVector[j] = Problem->Flowerbounds[j+1] - FVector[j+1];
				else
					ConstraintViolationVector[j] = 0.0;
			}

			//return the norm of the constraints
			*ConstraintNormValue = EMTG::math::norm(ConstraintViolationVector.data(), nF - 1);
		}
	}
	

	//main ACDE loop
	void TestAlg1::Evolve()
	{
		cout << "Evolve Called\n";

		int evoIter = 1000; // Make Alg Param
		double varType;
		double mutCap = 1; // Make Alg Param
		double F_w = 0.85; // Make Alg Param
		vector<double> bn (3);
		vector<double> num (3);
		for (int i = 0; i < evoIter; i++)
		{
			varType = DoubleDistribution(RNG);
			if(varType < mutCap)
			{
				for (int j = 0; j < 3; j++)
				{
					bn[j] = binDist(RNG);
					num[j] = capDist(RNG);
				}
				for (int j = 0; j < nX; j++)
					xTemp[j] = Bins[bn[0]][num[0]][j] + F_w*(Bins[bn[1]][num[1]][j] - Bins[bn[2]][num[2]][j]);
			}
			else
			{
				// GA recomb
			}
			LocOptInd(xTemp);
			EvalInd(xTemp, &objValTemp, &conValTemp);
			vector<double> score (nF);
			score[0] = objValTemp + conValTemp;
			for (int j = 0; j < nF - 1; ++j)
			{
				score[j+1] = objValTemp + ConstraintViolationVector[j];
			}
			for (int j = 0; j < nF; j++)
			{
				if(score[j] < Fits[j][NP-1][j])
					SortBin(Bins[j],Fits[j],xTemp,score,NP-1,j);
			}
		}

		BestX = Bins[1][1];
		Problem->unscale(BestX.data());
		Problem->evaluate(Problem->X.data(), FVector.data(), Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);
		//BestObjectiveValue = 0;
		//Problem->unscale(&BestX[0]);
		Problem->Xopt = Problem->X; //we store the unscaled Xbest
		//Problem->Xopt = BestX;
		Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
		Problem->output();
	}

}}//close namespaces