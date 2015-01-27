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
		Mission = (mission*)(Problem); // Needed for phase parcing of genome
		pM = 0; // Total phase counter
		phi = (1.0 + sqrt(5.0)) / 2.0; // Golden Ratio
		resphi = 2.0 - phi; // Needed for GoldenSearch
		for (int j = 0; j < Problem->options.number_of_journeys; j++)
		{
			for (int p = 0; p < Problem->options.number_of_phases[j]; p++)
			{
				phasePts.push_back(Mission->journeys[j].phases[p].first_X_entry_in_phase); // Collection of phase indices in genome
				pM++;
			}
		}

		////extract the dimensions of the problem
		nX = Problem->total_number_of_NLP_parameters;
		nF = Problem->total_number_of_constraints;
		NP = 100; // Bin Depth -- Make Alg Param
		genC = 0.1; // Fraction of total iterations allowed to serve as stagntion tolerance -- Make Alg Param
		gC = 0; // stagnation generation counter for local optimization
		gCb = 0; // stagnation generation counter for evolve
		BestObjectiveValue = EMTG::math::LARGE;
		BestConstraintNorm = EMTG::math::LARGE;
		ConstraintViolationVector.resize(nF - 1);
		int gens = 500; // Random Init Generations -- Make Alg Param
		FVector.resize(nF);
		xTemp.resize(nX);
		DoubleDistribution = boost::uniform_real<>(0.0, 1.0);
		capDist = boost::uniform_int<>(0, NP-1); // Used for random member of bin selection
		binDist = boost::uniform_int<>(0, nF-1); // Used for random bin selection
		RNG.seed(time(NULL));
		Bins.resize(nF); // Collection of best trials in bins by constraint
		Fits.resize(nF); // Fitness values corresponding to maintained population
		varPer = 0.000001; // Used as small genome walk value in local optimization
		zeros.resize(nX); // vector of all zeros to allow the norm of vector difference method to also calculate norms
		for (int i = 0; i < nX; i++)
			zeros[i] = 0.0;

		killDist = 1000*Problem->options.FD_stepsize; // Distance within which insertion becomes replacement -- Make Alg Param
		proxMode = true; // Whether or not to kill near duplicates -- Make Alg Param
		proxKill = false; // Used to determine whether or not a deletion has taken place
		elitePer = 0.1*NP; // Index above which ranked placement resets the stagnation counter -- Make Alg Param 
		mutStr = ""; // Used in cmd to notify when the best had undergone mutation

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
		recur = 0; // Local optimization recursion-depth counter 
		recurLim = 7000; // Local optimization recursion-depth limit -- Make Alg Param
		RandInit(gens);
	}

	// Make a sufficient number of random candidates to fill population and allow for differentiation
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
		gC = 0;
		int j = 0;
		while (j < iter && gC < genC*iter)
		{
			vector<double> Xt (nX);
			for (int i = 0; i < nX; i++)
			{
				double temp = Xb[i] + varPer * (DoubleDistribution(RNG)-0.5);
				if(temp > 1.0)
					temp = 2.0-temp;
				if(temp < 0.0)
					temp = -1.0*temp;
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
		vector<double> Xb = X; // Best so far
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conValTemp; // Score of best
		vector<double> ders (nX); // Directional derivatives
		double maxDer = 0; // Max derivative in any direction
		double closestDist = 2; // Distance closest dimension to boundary
		vector<double> bndDist (nX); // Distance to boundary of every dimension
		for (int i = 0; i < nX; i++)
		{
			if(1-Xb[i] < Xb[i]) // Learn boundary
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
			vector<double> Xt (nX); // An single dimension offset of Xb to calculate derivatives
			for (int j = 0; j < nX; j++) // Make Xt
			{
				if(j==i)
				{
					Xt[j] = Xb[j] + Problem->options.FD_stepsize;
					if(Xt[i] > 1.0)
						Xt[i] = 2.0-Xt[i];
					if(Xt[i] < 0.0 && Xt[i] >= -1.0)
						Xt[i] = -1.0*Xt[i];
					if(Xt[i] < -1.0)
						Xt[i] = DoubleDistribution(RNG); // Should never happen;
				}
				else
					Xt[j] = Xb[j];
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conValTemp; // Score of Xt
			ders[i] = (St-Sb)/Problem->options.FD_stepsize; // Simple derivative calculation
			if(abs(ders[i]) > maxDer)
				maxDer = abs(ders[i]);
		}
		
		if(maxDer < Problem->options.FD_stepsize)
			maxDer = Problem->options.FD_stepsize; // Stop divide by 0 error
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
			if(bndPt[i] > 1.0)
				bndPt[i] = 1.0;
			if(bndPt[i] < 0.0)
				bndPt[i] = 0.0;
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
		vector<double> ptB (nX);
		for (int i = 0; i < nX; i++)
			ptB[i] = Xb[i] + (bndPt[i]-Xb[i])/3.0;
		recur = 0;
		X = GoldenSearch(Xb,ptB,bndPt,Problem->options.FD_stepsize,dir);
	}

	// Recursive Golden Search Implementation
	vector<double> TestAlg1::GoldenSearch(vector<double>& A, vector<double>& B, vector<double>& C, double tau, vector<double>& dir)
	{
		recur++;
		vector<double> x (nX);
		if(NormOfDif(C,B) > NormOfDif(B,A))
		{
			for (int i = 0; i < nX; i++)
				x[i] = B[i] + resphi * NormOfDif(C,B)*dir[i];
		}
		else
			for (int i = 0; i < nX; i++)
				x[i] = B[i] - resphi * NormOfDif(B,A)*dir[i];
		if(NormOfDif(C,A) < tau*(NormOfDif(B,zeros) + NormOfDif(x,zeros)) || recur > recurLim)
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
		int offset = 0; // Used for both index of deletion and eraser flag
		if(cap>0 && score[ind] < fit[cap-1][ind])
			SortBin(bin,fit,X,score,cap-1,ind);
		else
		{
			if(proxMode) // True = do not allow near duplicates
			{
				if(NormOfDif(bin[cap],X) < killDist) // Kill below
				{
					proxKill = true;
					offset++;
				}
				if(cap > 0 && NormOfDif(bin[cap-1],X) < killDist) // Kill self
				{
					proxKill = true;
					offset--;
				}
				if(proxKill && offset==0) // Find closer one to kill
				{
					if(NormOfDif(bin[cap],X) < NormOfDif(bin[cap-1],X))
						offset++;
					else
						offset--;
				}
			}
			if(offset!=-1) // Prevents insertion when self is dupicate
			{
				vector<vector<double>>::iterator it;
				it = bin.begin();
				it+=cap;
				it = bin.insert (it,X);
				vector<vector<double>>::iterator it2;
				it2 = fit.begin();
				it2+=cap;
				it2 = fit.insert (it2,score);
			}
			if(offset > 0 && cap+1 > NP) // Kill below
			{
				bin.erase(bin.begin()+cap+offset);
				fit.erase(fit.begin()+cap+offset);
			}
			string tStr;
			switch (mode) // For tracking success
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
				cout << "\nBin " << ind << " has a new best at Fit = " << score[ind] << " by " << tStr << mutStr;
			if(cap < elitePer)
				gCb = 0;
		}
		if(cap==NP-1 && !proxKill) // Only delete last if no one was killed
		{
			bin.pop_back();
			fit.pop_back();
		}
	}

	// Evaluate Individual -- borrowed from ACDE
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
		cout << "\nEvolve Called\n";

		int evoIter = 50000; // Number of candidates to make -- Make Alg Param
		bool deMode = true; // True = unique DE vector choices -- Make Alg Param
		double varType;
		double mutCap = 0.5; // DE/GA trade-off point -- Make Alg Param
		double F_w = 0.85; // DE parameter -- Make Alg Param
		double mutProb = 0.2; // Bit-wise Pareto mutation probability -- Make Alg Param
		double alpha = Problem->options.MBH_Pareto_alpha; // Needed for pareto
		vector<double> bn (3); // Random bin index holder
		vector<double> num (3); // Random bin entry index holder
		for (int i = 0; i < 3; i++)
		{
			bn[i] = -1;
			num[i] = -1;
		}

		int i = 0;
		gCb = 0;
		while (i < evoIter && gCb < genC*evoIter)
		{
			mutStr = "";
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
					} while(deMode && ((bn[j]==bn[(j+1)%3] && num[j]==num[(j+1)%3]) || (bn[j]==bn[(j+2)%3] && num[j]==num[(j+2)%3]))); // Forces reselection if choice is duplicate
				}
				for (int j = 0; j < nX; j++)
				{
					double temp = Bins[bn[0]][num[0]][j] + F_w*(Bins[bn[1]][num[1]][j] - Bins[bn[2]][num[2]][j]);
					if(temp > 1.0)
						temp = 2.0-temp;
					if(temp < 0.0 && temp >= -1.0)
						temp = -1.0*temp;
					if(temp < -1.0)
						temp = DoubleDistribution(RNG);
					xTemp[j] = temp;
				}
			}
			else // GA recombination
			{
				mode = 2;// GA flag
				for (int j = 0; j < pM; j++) // Construct trial by phase combination
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
			if(DoubleDistribution(RNG) < mutProb) // Pareto mutation
			{
				mutStr = " with Mut";
				for (int j = 0; j < nX; j++)
				{
					double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
					int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;
					double perturb = s * Problem->options.MBH_max_step_size * r;
					double temp = xTemp[j] + perturb;
					if(temp > 1.0)
					{
						if(temp < 2.0 && abs(xTemp[j]-1) < Problem->options.FD_stepsize)
						{
							temp = 1.0;
							//cout << "Bound Upper";
						}
						else
						{
							temp = DoubleDistribution(RNG);
							//cout << "\nReset Upper";
						}
					}
					else if(temp < 0.0)
					{
						if(temp > -1.0 && abs(xTemp[j]) < Problem->options.FD_stepsize)
						{
							temp = 0.0;
							//cout << "Bound Lower";
						}
						else
						{
							temp = DoubleDistribution(RNG);
							//cout << "\nReset Lower";
						}
					}
					xTemp[j] = temp;
				}
			}
			if(gCb == 1)
				cout << "\n";
			cout <<".";
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
			cout << "\nEvolution terminated due to stagnation.\n";
		else
			cout << "\n";
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