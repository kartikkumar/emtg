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
		Problem->output_problem_bounds_and_descriptions(Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "XFfile.csv");
		Mission = (mission*)(Problem);
		pM = 0;
		phi = (1.0 + sqrt(5.0)) / 2.0;
		resphi = 2.0 - phi;
		for (int j = 0; j < Problem->options.number_of_journeys; j++)
		{
			for (int p = 0; p < Problem->options.number_of_phases[j]; p++)
			{
				phasePts.push_back(Mission->journeys[j].phases[p].first_X_entry_in_phase);
				pM++;
			}
		}

		////extract the dimensions of the problem
		nX = Problem->total_number_of_NLP_parameters;
		nF = Problem->total_number_of_constraints;
		NP = 100; // Bin Depth -- Make Alg Param
		genC = 0.1; // Fraction of total iterations allowed to serve as stagntion tolerance -- Make Alg Param
		gC = 0;
		gCb = 0;
		BestObjectiveValue = EMTG::math::LARGE;
		BestConstraintNorm = EMTG::math::LARGE;
		ConstraintViolationVector.resize(nF - 1);
		int gens = 500; // Random Init Generations -- Make Alg Param
		FVector.resize(nF);
		xTemp.resize(nX);
		DoubleDistribution = boost::uniform_real<>(0.0, 1.0);
		capDist = boost::uniform_int<>(0, NP-1);
		binDist = boost::uniform_int<>(0, nF-1);
		RNG.seed(time(NULL));
		Bins.resize(nF);
		Fits.resize(nF);
		varPer = 0.000001;
		zeros.resize(nX);
		for (int i = 0; i < nX; i++)
			zeros[i] = 0.0;

		killDist = Problem->options.FD_stepsize; // Make Alg Param
		proxMode = true; // Make Alg Param
		proxKill = false;

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

		iter = 500; // Local optimization iterations -- Make Alg Param

		RandInit(gens);
	}

	void TestAlg1::RandInit(int gens)
	{
		mode = 0;// Rand flag
		for (int i = 0; i < gens; i++)
		{
			for (int j = 0; j < nX; j++)
				xTemp[j] = DoubleDistribution(RNG);

			LocOptInd_SteepDecent(xTemp);
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

	// Locally Optimize Individual -- here by random walk
	void TestAlg1::LocOptInd(vector<double>& X)
	{
		vector<double> Xb = X;
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conValTemp;
		//for (int j = 0; j < iter; j++)
		//{
		//	vector<double> Xt (nX);
		//	for (int i = 0; i < nX; i++)
		//	{
		//		double temp = Xb[i] + vars[i] * VarDist(RNG);
		//		if(temp > Problem->Xupperbounds[i])
		//			temp = 2*Problem->Xupperbounds[i]-temp;
		//		if(temp < Problem->Xlowerbounds[i])
		//			temp = 2*Problem->Xlowerbounds[i]-temp;
		//		Xt[i] = temp;
		//	}
		//	EvalInd(Xt, &objValTemp, &conValTemp);
		//	double St = objValTemp + conValTemp;
		//	if(St < Sb)
		//	{
		//		Xb = Xt;
		//		Sb = St;
		//	}
		//}
		gC = 0;
		int j = 0;
		while (j < iter && gC < genC*iter)
		{
			vector<double> Xt (nX);
			for (int i = 0; i < nX; i++)
			{
				double temp = Xb[i] + varPer * (DoubleDistribution(RNG)-0.5);
				if(temp > 1)
					temp = 2-temp;
				if(temp < 0)
					temp = -1*temp;
				Xt[i] = temp;
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conValTemp;
			if(St < Sb)
			{
				Xb = Xt;
				Sb = St;
				gC = 0;
			}
			j++;
			gC++;
		}
		X = Xb;
	}

	// Locally Optimize Individual -- using steepest decent
	void TestAlg1::LocOptInd_SteepDecent(vector<double>& X)
	{
		vector<double> Xb = X;
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conValTemp;
		vector<double> Xt (nX);
		vector<double> ders (nX);
		double step = 0.0001; // Make Alg Param
		double maxDer = 0;
		double closestDist = 2;
		vector<double> bndDist (nX);
		//for (int i = 0; i < nX; i++)
		//{
		//	Xt[i] = Xb[i] + Problem->options.FD_stepsize;
		//	if(Xt[i] > 1)
		//		Xt[i] = 1-Xt[i];
		//	ders[i] = (Xt[i]-Xb[i])/Problem->options.FD_stepsize;
		//	if(abs(ders[i]) > maxDer)
		//		maxDer = abs(ders[i]);
		//	if(1-Xb[i] < closestDist)
		//		closestDist = 1-Xb[i];
		//	else if(Xb[i] < closestDist)
		//		closestDist = Xb[i];
		//}
		for (int i = 0; i < nX; i++)
		{
			if(1-Xb[i] < Xb[i])
			{
				bndDist[i] = 1-Xb[i];
				if(1-Xb[i] < closestDist)
					closestDist = 1-Xb[i];
			}
			else
			{
				bndDist[i] = Xb[i];
				if(Xb[i] < closestDist)
					closestDist = Xb[i];
			}
			vector<double> Xt (nX);
			for (int j = 0; j < nX; j++)
			{
				if(j==i)
				{
					Xt[j] = Xb[j] + Problem->options.FD_stepsize;
					if(Xt[i] > 1)
						Xt[i] = 2-Xt[i];
					if(Xt[i] < 0 && Xt[i] >= -1)
						Xt[i] = -1*Xt[i];
					if(Xt[i] < -1)
						Xt[i] = DoubleDistribution(RNG);
				}
				else
					Xt[j] = Xb[j];
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conValTemp; 
			ders[i] = (St-Sb)/Problem->options.FD_stepsize;
			if(abs(ders[i]) > maxDer)
				maxDer = abs(ders[i]);
		}
		
		if(maxDer < Problem->options.FD_stepsize)
			maxDer = Problem->options.FD_stepsize;
		//if(maxDer!=maxDer)
		//	maxDer = EMTG::math::LARGE;
		//for(int i = 0; i < nX; i++)
		//{
		//	if(!(ders[i] == ders[i]))
		//	{
		//		ders[i] = EMTG::math::LARGE;
		//		cout << "NaN BS\n";
		//	}
		//}

		// This is to find out which dimension will first hit the boundary along the line made by the derivatives
		double numCl = EMTG::math::LARGE;
		int indCl = -1;
		for (int i = 0; i < nX; i++)
		{
			double numTemp = bndDist[i]/abs((ders[i]/maxDer)*closestDist);
			if(numTemp < numCl)
			{
				numCl = numTemp;
				indCl = i;
			}
		}

		// Find boundary point hit
		vector<double> bndPt (nX);
		for (int i = 0; i < nX; i++)
		{
			bndPt[i] = Xb[i] - numCl * ((ders[i]/maxDer)*closestDist);
			if(bndPt[i] > 1)
				bndPt[i] = 1;
			if(bndPt[i] < 0)
				bndPt[i] = 0;
		}

		// Make directional vector from Xb to boundary (inefficient at the moment since this can be constructed from known data)
		vector<double> dir (nX);
		double tot = 0;
		for (int i = 0; i < nX; i++)
		{
			dir[i] = bndPt[i] - Xb[i];
			tot += pow(dir[i],2);
		}
		tot = sqrt(tot);
		if(tot < Problem->options.FD_stepsize) // Divide by 0 prevention
			tot = Problem->options.FD_stepsize;
		for (int i = 0; i < nX; i++)
			dir[i] = dir[i]/tot;
		//EvalInd(bndPt, &objValTemp, &conValTemp);
		//double Sbnd = objValTemp + conValTemp;
		vector<double> ptB (nX);
		for (int i = 0; i < nX; i++)
			ptB[i] = Xb[i] + (bndPt[i]-Xb[i])/3.0;
		X = GoldenSearch(Xb,ptB,bndPt,Problem->options.FD_stepsize,dir);
	 //   bool rep = true;
		//while (rep)
		//{
		//	for (int i = 0; i < nX; i++)
		//	{
		//		Xt[i] = -1*(ders[i]/maxDer)*step + Xb[i];
		//		if(Xt[i] > 1)
		//			Xt[i] = 2-Xt[i];
		//		if(Xt[i] < 0 && Xt[i] >= -1)
		//			Xt[i] = -1*Xt[i];
		//		if(Xt[i] < -1)
		//			Xt[i] = DoubleDistribution(RNG);
		//	}
		//	EvalInd(Xt, &objValTemp, &conValTemp);
		//	double St = objValTemp + conValTemp;
		//	if(St < Sb)
		//	{
		//		Xb = Xt;
		//		Sb = St;
		//	}
		//	else
		//		rep = false;
		//}
	}

	// Recursive Golden Search Implementation
	vector<double> TestAlg1::GoldenSearch(vector<double>& A, vector<double>& B, vector<double>& C, double tau, vector<double>& dir)
	{
		vector<double> x (nX);
		if(NormOfDif(C,B) > NormOfDif(B,A))
		{
			for (int i = 0; i < nX; i++)
				x[i] = B[i] + resphi * NormOfDif(C,B)*dir[i];
		}
		else
			for (int i = 0; i < nX; i++)
				x[i] = B[i] - resphi * NormOfDif(B,A)*dir[i];
		if(NormOfDif(C,A) < tau*(NormOfDif(B,zeros) + NormOfDif(x,zeros)))
		{
			vector<double> ret (nX);
			for (int i = 0; i < nX; i++)
				ret[i] = A[i] + (NormOfDif(C,A)/2.0)*dir[i];
			return ret;
		}
		EvalInd(x, &objValTemp, &conValTemp);
		double Sx = objValTemp + conValTemp;
		EvalInd(B, &objValTemp, &conValTemp);
		double Sb = objValTemp + conValTemp;
		if(Sx < Sb)
		{
			if(NormOfDif(C,B) > NormOfDif(B,A))
				return GoldenSearch(B,x,C,tau,dir);
			else
				return GoldenSearch(A,x,B,tau,dir);
		}
		else
		{
			if(NormOfDif(C,B) > NormOfDif(B,A))
				return GoldenSearch(A,B,x,tau,dir);
			else
				return GoldenSearch(x,B,C,tau,dir);
		}
	}

	// Locally Optimize Individual -- linear walk with reset
	void TestAlg1::LocOptInd_LinWalk(vector<double>& X)
	{
		vector<double> Xb = X;
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conValTemp;
		gC = 0;
		int j = 0;
		while (j < iter && gC < genC*iter)
		{
			vector<double> xVar (nX);
			for (int i = 0; i < nX; i++)
				xVar[i] = varPer * (DoubleDistribution(RNG)-0.5);
			bool fail = false;
			while (!fail)
			{
				vector<double> Xt (nX);
				for (int i = 0; i < nX; i++)
				{
					Xt[i] = Xb[i] + xVar[i];
					if(Xt[i] > 1)
						Xt[i] = 2-Xt[i];
					if(Xt[i] < 0)
						Xt[i] = -1*Xt[i];
				}
				EvalInd(Xt, &objValTemp, &conValTemp);
				double St = objValTemp + conValTemp;
				if(St < Sb)
				{
					Xb = Xt;
					Sb = St;
					gC = 0;
				}
				else
					fail = true;
			}
			j++;
			gC++;
		}
		X = Xb;
	}

	// Recursive bin sorting method
	void TestAlg1::SortBin(vector<vector<double>>& bin,vector<vector<double>>& fit,vector<double>& X,vector<double>& score,int cap,int ind)
	{
		proxKill = false;
		if(cap>0 && score[ind] < fit[cap-1][ind])
			SortBin(bin,fit,X,score,cap-1,ind);
		else
		{
			if(proxMode && NormOfDif(bin[cap],X) < killDist);
				proxKill = true;
			vector<vector<double>>::iterator it;
			it = bin.begin();
			it+=cap;
			it = bin.insert (it,X);
			vector<vector<double>>::iterator it2;
			it2 = fit.begin();
			it2+=cap;
			it2 = fit.insert (it2,score);
			//bin.insert(cap,X);
			if(proxKill && cap+1 > NP)
			{
				bin.erase(bin.begin()+cap+1);
				fit.erase(fit.begin()+cap+1);
			}
			string tStr;
			switch (mode)
			{
				case 0:
				{
					tStr = "Random";
					break;
				}
				case 1:
				{
					tStr = "DE";
					break;
				}
				case 2:
				{
					tStr = "GA";
					break;
				}
			}
			if(cap==0)
			{
				cout << "Bin " << ind << " has a new best at Fit = " << score[ind] << " by " << tStr << "\n";
				gCb = 0;
			}
		}
		if(cap==NP-1 && !proxKill)
		{
			bin.pop_back();
			fit.pop_back();
		}
	}

	// Evaluate Individual
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
					ConstraintViolationVector[j] = 0.0; // formerly 0.0
			}

			//return the norm of the constraints
			*ConstraintNormValue = EMTG::math::norm(ConstraintViolationVector.data(), nF - 1);
		}
	}
	
	// Norm of the difference of two vectors
	double TestAlg1::NormOfDif(vector<double>& A,vector<double>& B)
	{
		double tot = 0;
		for (int i = 0; i < nX; i++)
			tot += pow((B[i]-A[i]),2.0);
		return sqrt(tot);
	}

	//main algorithm loop
	void TestAlg1::Evolve()
	{
		cout << "Evolve Called\n";

		int evoIter = 30000; // Make Alg Param
		bool deMode = true; // Make Alg Param
		double varType;
		double mutCap = 0.5; // Make Alg Param
		double F_w = 0.85; // Make Alg Param
		vector<double> bn (3);
		vector<double> num (3);
		int i = 0;
		gCb = 0;
		while (i < evoIter && gCb < genC*evoIter)
		{
			varType = DoubleDistribution(RNG);
			if(varType < mutCap) // DE mutation
			{
				mode = 1;// DE flag
				for (int j = 0; j < 3; j++)
				{
					int safety = 0;
					do
					{
						bn[j] = binDist(RNG);
						num[j] = capDist(RNG);
						safety++;
						if(safety > 1000)
							cout << "Infinite Loop -- KILL IT!\n";
					} while(deMode && ((bn[j]==bn[(j+1)%3] && num[j]==num[(j+1)%3]) || (bn[j]==bn[(j+2)%3] && num[j]==num[(j+2)%3])));
				}
				for (int j = 0; j < nX; j++)
				{
					double temp = Bins[bn[0]][num[0]][j] + F_w*(Bins[bn[1]][num[1]][j] - Bins[bn[2]][num[2]][j]);
					if(temp > 1)
						temp = 2-temp;
					if(temp < 0 && temp >= -1)
						temp = -1*temp;
					if(temp < -1)
						temp = DoubleDistribution(RNG);
					xTemp[j] = temp;
				}
			}
			else // GA recombination
			{
				mode = 2;// GA flag
				for (int j = 0; j < pM; j++)
				{
					int f = phasePts[j];
					int e = nX;
					if(j != pM-1)
					{
						e = phasePts[j+1];
					}
					int b = binDist(RNG);
					int c = capDist(RNG);
					for (int k = f; k < e; k++)
						xTemp[k] = Bins[b][c][k];
				}
			}
			// Add bit-wise mutation probability
			LocOptInd_SteepDecent(xTemp);
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
			i++;
			gCb++;
		}
		if (!(gCb < genC*evoIter))
		{
			cout << "Evolution terminated due to stagnation.\n";
		}
		/*for (int i = 0; i < pM; i++) // check phasePts for elements
			cout << "Entry " << i << " is " << phasePts[i] << "\n";*/

		/*for (int i = 0; i < nF; i++) // sanity check to ensure bins didn't change size
			cout << "Bin " << i << " cap = " << Bins[i].size() << "\n";*/

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