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
		genC = 2000; // Percentage of total iterations allowed to serve as stagntion tolerance -- Make Alg Param
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
		varPer = Problem->options.MBH_max_step_size;
		zeros.resize(nX); // vector of all zeros to allow the norm of vector difference method to also calculate norms
		
		for (int i = 0; i < nX; i++)
			zeros[i] = 0.0;
		trashFits.resize(nF);
		for (int i = 0; i < nF; i++)
			trashFits[i] = EMTG::math::LARGE;
		killDist = 40000000*Problem->options.FD_stepsize; // Distance within which insertion becomes replacement -- Make Alg Param
		proxMode = true; // Whether or not to kill near duplicates -- Make Alg Param
		proxKill = false; // Used to determine whether or not a deletion has taken place
		elitePer = 0.05*NP; // Index above which ranked placement resets the stagnation counter -- Make Alg Param 
		mutStr = ""; // Used in cmd to notify when the best had undergone mutation

		for (int i = 0; i < nF; i++)
		{
			Bins[i].resize(NP);
			Fits[i].resize(NP);
			for (int j = 0; j < NP; j++)
			{
				Bins[i][j].resize(nX);
				Fits[i][j].resize(nF);
				Fits[i][j] = trashFits;
			}
		}

		iter = 1000; // Local optimization iterations -- Make Alg Param
		recur = 0; // Local optimization recursion-depth counter 
		recurLim = 20; // Local optimization recursion-depth limit -- Make Alg Param
		optLim = 1;
		conWeight = 50; // Weighting of constraints -- Make Alg Param
		conLW = 5;
		steepStep = 10^-5; // Max steep step -- Make Alg Param
		ders.resize(nX);
		RandInit(gens,false);
	}

	// Make a sufficient number of random candidates to fill population and allow for differentiation
	void TestAlg1::RandInit(int gens,bool sed)
	{
		mode = 0;// Rand flag
		double alpha = Problem->options.MBH_Pareto_alpha; // Needed for pareto
		vector<double> xEph (nX);
		for (int i = 0; i < gens; i++)
		{
			if(!sed)
			{
				for (int j = 0; j < nX; j++)
					xEph[j] = DoubleDistribution(RNG);
			}
			else
			{
				for (int j = 0; j < nX; j++)
				{
					for (int k = nF-1; k >= 0; k--)
					{
						double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
						int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;
						double perturb = s * varPer * r;
						double temp = Bins[k][0][j] + perturb;
						if(temp > 1.0)
						{
							//if(temp < 2.0 && abs(xEph[j]-1) < Problem->options.FD_stepsize)
							//{
							//	temp = 1.0;
							//	//cout << "Bound Upper";
							//}
							//else
							//{
							//	temp = DoubleDistribution(RNG);
							//	//cout << "\nReset Upper";
							//}
							temp = 1.0;
						}
						else if(temp < 0.0)
						{
							//if(temp > -1.0 && abs(xEph[j]) < Problem->options.FD_stepsize)
							//{
							//	temp = 0.0;
							//	//cout << "Bound Lower";
							//}
							//else
							//{
							//	temp = DoubleDistribution(RNG);
							//	//cout << "\nReset Lower";
							//}
							temp = 0.0;
						}
						xEph[j] = temp;
					}
				}
			}
			SubmitInd(xEph, &objValTemp, &conValTemp);
		}
	}

	// Local Optimization method holder to reduce difficulty in scheme change
	void TestAlg1::LocalOpt(vector<double>& X,int b)
	{
		for (int i = 0; i < optLim; i++)
		{
			LocOptInd_SteepDescentv2(X,b);
		}
	}

	// Locally Optimize Individual -- here by random walk
	void TestAlg1::LocOptInd_RandWalk(vector<double>& X,int b)
	{
		vector<double> Xb = X;
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conWeight * conValTemp; // Score of best
		if(b!=0)
			Sb = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
		double alpha = Problem->options.MBH_Pareto_alpha; // Needed for pareto
		
		gC = 0;
		int j = 0;
		while (j < iter && gC < recurLim)
		{
			vector<double> Xt (nX);
			for (int i = 0; i < nX; i++)
			{
				//double temp = Xb[i] + varPer * (DoubleDistribution(RNG)-0.5);
				double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
				int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;
				double perturb = s * varPer * r;
				double temp = Xb[i] + perturb;
				if(temp > 1.0)
					temp = 1.0;
				if(temp < 0.0)
					temp = 0.0;
				Xt[i] = temp;
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conWeight * conValTemp;
			if(b!=0)
				St = objValTemp + conWeight * ConstraintViolationVector[j] + conLW*conValTemp;
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

	void TestAlg1::LocOptInd_SteepDescentv2(vector<double>& X,int b)
	{
		int safe = 0;
		vector<double> Xb (nX);
		for (int i = 0; i < nX; i++)
			Xb[i] = X[i];
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conWeight * conValTemp; // Score of best
		if(b!=0)
			Sb = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
		vector<double> X2 (nX);
		vector<double> X3 (nX);
		vector<double> dir (nX);
		bool reset = false;
		for (int i = 0; i < nX; i++)
		{
			vector<double> Xt (nX);
			for (int j = 0; j < nX; j++)
			{
				if(j==i)
				{
					Xt[j] = Xb[j] + Problem->options.FD_stepsize;
					if(Xt[i] > 1.0)
						Xt[i] = 1.0;
					if(Xt[i] < 0.0)
						Xt[i] = 0.0;
				}
				else
					Xt[j] = Xb[j];
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conWeight * conValTemp; // Score of Xt
			if(b!=0)
				St = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
			ders[i] = (St-Sb)/Problem->options.FD_stepsize; // Simple derivative calculation
		}
		double derNorm = NormOfDif(ders,zeros);
		while (derNorm > 10*Problem->options.FD_stepsize && safe < iter)
		{
			for (int i = 0; i < nX; i++)
			{
				dir[i] = ders[i]/derNorm;
				X2[i] = Xb[i] - steepStep * dir[i];
				if(X2[i] > 1.0)
				{
					X2[i] = 1.0;
					//cout << "]";
				}
				if(X2[i] < 0.0)
				{
					X2[i] = 0.0;
					//cout << "[";
				}
				double nTemp = NormOfDif(Xb,X2);
				for (int j = 0; j < nX; j++)
					dir[i] = (X2[i] - Xb[i])/nTemp;
				X3[i] = Xb[i] + resphi*(X2[i]-Xb[i]);
				if(X3[i] > 1.0)
				{
					X3[i] = 1.0;
					cout << "\nSteep Bound Impossibility!\n";
				}
				if(X3[i] < 0.0)
				{
					X3[i] = 0.0;
					cout << "\nSteep Bound Impossibility!\n";
				}
			}
			Xb = GoldenSearch(Xb,X3,X2,Problem->options.FD_stepsize,dir,b);
			EvalInd(Xb, &objValTemp, &conValTemp);
			double Sb = objValTemp + conWeight * conValTemp; // Score of best
			if(b!=0)
				Sb = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
			reset = !(NormOfDif(Xb,X2) < Problem->options.FD_stepsize);
			if(reset)
			{
				for (int i = 0; i < nX; i++)
				{
					vector<double> Xt (nX);
					for (int j = 0; j < nX; j++)
					{
						if(j==i)
						{
							Xt[j] = Xb[j] + Problem->options.FD_stepsize;
							if(Xt[i] > 1.0)
								Xt[i] = 1.0;
							if(Xt[i] < 0.0)
								Xt[i] = 0.0;
						}
						else
							Xt[j] = Xb[j];
					}
					EvalInd(Xt, &objValTemp, &conValTemp);
					double St = objValTemp + conWeight * conValTemp; // Score of Xt
					if(b!=0)
						St = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
					ders[i] = (St-Sb)/Problem->options.FD_stepsize; // Simple derivative calculation
				}
				derNorm = NormOfDif(ders,zeros);
			}
			safe++;
		} 
		if(!(safe < iter))
			cout << "\nSteepDescentv2 LOOP STOP\n";
		X = Xb;
	}

	// Locally Optimize Individual -- using steepest decent
	void TestAlg1::LocOptInd_SteepDescent(vector<double>& X,int b)
	{
		vector<double> Xb = X; // Best so far
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conWeight * conValTemp; // Score of best
		if(b!=0)
			Sb = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
		vector<double> ders (nX); // Directional derivatives
		double maxDer = 0; // Max derivative in any direction
		double closestDist = 5; // Distance closest dimension to boundary
		vector<double> bndDist (nX); // Distance to boundary of every dimension
		for (int i = 0; i < nX; i++)
		{
			//if(1-Xb[i] < Xb[i]) // Learn boundary
			//{
			//	bndDist[i] = 1-Xb[i];
			//	if(1-Xb[i] < closestDist)
			//		closestDist = 1-Xb[i];
			//}
			//else
			//{
			//	bndDist[i] = Xb[i];
			//	if(Xb[i] < closestDist)
			//		closestDist = Xb[i];
			//}
			vector<double> Xt (nX); // An single dimension offset of Xb to calculate derivatives
			for (int j = 0; j < nX; j++) // Make Xt
			{
				if(j==i)
				{
					Xt[j] = Xb[j] + Problem->options.FD_stepsize;
					if(Xt[i] > 1.0)
						Xt[i] = 1.0;
					if(Xt[i] < 0.0)
						Xt[i] = 0.0;
				}
				else
					Xt[j] = Xb[j];
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conWeight * conValTemp; // Score of Xt
			if(b!=0)
				St = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
			ders[i] = (St-Sb)/Problem->options.FD_stepsize; // Simple derivative calculation
			if(abs(ders[i]) > maxDer)
				maxDer = abs(ders[i]);
			if(ders[i] < 0) // Suspect direction
				bndDist[i] = 1-Xb[i];
			else if(ders[i] > 0)
				bndDist[i] = Xb[i];
			else
				bndDist[i] = 5;
			if(bndDist[i] < closestDist)
				closestDist = bndDist[i];
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
			bndPt[i] = Xb[i] - numCl * ((ders[i]/maxDer)*closestDist); // Notice sign suspicion
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
			ptB[i] = Xb[i] + (bndPt[i]-Xb[i])*resphi;
		recur = 0;
		X =  GoldenSearch(Xb,ptB,bndPt,Problem->options.FD_stepsize,dir,b);
	}

	void TestAlg1::LocOptInd_ParetoDescent(vector<double>& X,int b)
	{
		vector<double> Xb = X; // Best so far
		double alpha = Problem->options.MBH_Pareto_alpha; // Needed for pareto
		EvalInd(Xb, &objValTemp, &conValTemp);
		double Sb = objValTemp + conWeight * conValTemp; // Score of best
		if(b!=0)
			Sb = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
		vector<double> ders (nX); // Directional derivatives
		double maxDer = 0; // Max derivative in any direction
		double closestDist = 5; // Distance closest dimension to boundary
		vector<double> bndDist (nX); // Distance to boundary of every dimension
		vector<double> Xt (nX); // An single dimension offset of Xb to calculate derivatives
		for (int i = 0; i < nX; i++)
		{
			vector<double> Xt (nX); // An single dimension offset of Xb to calculate derivatives
			for (int j = 0; j < nX; j++) // Make Xt
			{
				if(j==i)
				{
					Xt[j] = Xb[j] + Problem->options.FD_stepsize;
					if(Xt[i] > 1.0)
						Xt[i] = 1.0;
					if(Xt[i] < 0.0)
						Xt[i] = 0.0;
				}
				else
					Xt[j] = Xb[j];
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conWeight * conValTemp; // Score of Xt
			if(b!=0)
				St = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
			ders[i] = (St-Sb)/Problem->options.FD_stepsize; // Simple derivative calculation
		}
		int j = 0;
		gC = 0;
		double nDers = NormOfDif(ders,zeros);
		for (int i = 0; i < nX; i++)
			ders[i] = ders[i]/nDers;
		while (j < iter && gC < recurLim)
		{
			double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
			int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;
			double perturb = s * varPer * r;
			for (int k = 0; k < nX; k++)
			{
				Xt[k] = Xb[k] + perturb*ders[k];
				if(Xt[k] > 1.0)
				{
					Xt[k] = 1.0;
					//cout << "]";
				}
				if(Xt[k] < 0.0)
				{
					Xt[k] = 0.0;
					//cout << "[";
				}
			}
			EvalInd(Xt, &objValTemp, &conValTemp);
			double St = objValTemp + conWeight * conValTemp; // Score of Xt
			if(b!=0)
				St = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
			if(St < Sb)
			{
				Xb = Xt;
				Sb = St;
				gC = -1;
			}
			j++;
			gC++;
		}
		X = Xb;
	}

	void TestAlg1::LocOptInd_HybridOpt(vector<double>& X,int b)
	{
		LocOptInd_ParetoDescent(X,b);
		LocOptInd_RandWalk(X,b);
	}

	// Recursive Golden Search Implementation
	std::vector<double> TestAlg1::GoldenSearch(vector<double>& A, vector<double>& B, vector<double>& C, double tau, vector<double>& dir, int b)
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
		double Sx = objValTemp + conWeight * conValTemp;
		if(b!=0)
			Sx = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
		EvalInd(B, &objValTemp, &conValTemp);
		double Sb = objValTemp + conWeight * conValTemp;
		if(b!=0)
			Sb = objValTemp + conWeight * ConstraintViolationVector[b-1] + conLW*conValTemp;
		if(Sx < Sb)
		{
			if(NormOfDif(C,B) > NormOfDif(B,A))
				return GoldenSearch(B,x,C,tau,dir,b);
			else
				return GoldenSearch(A,x,B,tau,dir,b);
		}
		else
		{
			if(NormOfDif(C,B) > NormOfDif(B,A))
				return GoldenSearch(A,B,x,tau,dir,b);
			else
				return GoldenSearch(x,B,C,tau,dir,b);
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
		while (j < iter && gC < genC)
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
			if(ind == 0 && cap < elitePer && offset!=-1)
				gCb = -1;				
			if(cap==0)
				cout << "\nBin " << ind << " new best at Fit = " << score[ind] << " => " << objValTemp << " by " << tStr << mutStr << " with Dist = " << NormOfDif(Bins[ind][0],Bins[ind][1]);
			if(offset > 0 && cap!=NP-1) // Kill below
			{
				bin.erase(bin.begin()+cap+offset);
				fit.erase(fit.begin()+cap+offset);
			}
		}
		if(cap==NP-1 && bin.size() > NP) // Only delete last if no one was killed
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
					ConstraintViolationVector[j] = abs(FVector[j+1] - Problem->Fupperbounds[j+1]);
				else if (FVector[j+1] < Problem->Flowerbounds[j+1])
					ConstraintViolationVector[j] = abs(Problem->Flowerbounds[j+1] - FVector[j+1]);
				else
					ConstraintViolationVector[j] = 0.0;
			}

			//return the norm of the constraints
			*ConstraintNormValue = EMTG::math::norm(ConstraintViolationVector.data(), nF - 1);
		}
	}
	
	// Does bin-wise optimization and sorting
	void TestAlg1::SubmitInd(vector<double> X, double* ObjectiveFunctionValue, double* ConstraintNormValue)
	{
		for (int i = 0; i < nF; i++)
		{
			//xTemp = X;
			for (int j = 0; j < nX; j++)
				xTemp[j] = X[j];
			LocalOpt(xTemp,i);
			EvalInd(xTemp,ObjectiveFunctionValue,ConstraintNormValue);
			vector<double> score (nF);
			score[0] = objValTemp + conWeight * conValTemp;
			if(score[0] < 0.1)
			{
				cout << "\nANOMALY DETECTED!\n";
				ofstream outputfile;
				outputfile.open(Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "_Anomaly.emtg", ios::out | ios::app);
				outputfile.precision(20);
				for (int j = 0; j < xTemp.size(); j++)
				{
					outputfile << xTemp[j];
					outputfile << ", ";
				}
				outputfile.close();
			}
			for (int j = 0; j < nF - 1; ++j)
				score[j+1] = objValTemp + conWeight * ConstraintViolationVector[j] + conLW*conValTemp;
			if(score[i] < Fits[i][NP-1][i])
			{
					SortBin(Bins[i],Fits[i],xTemp,score,NP-1,i);
					if(proxKill)
						cout << "x";
			}
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
		proxKill = false;
		double alpha = Problem->options.MBH_Pareto_alpha; // Needed for pareto
		int evoIter = 500000; // Number of candidates to make -- Make Alg Param
		bool deMode = true; // True = unique DE vector choices -- Make Alg Param
		double varType;
		double mutCap = 0.5; // DE/GA trade-off point -- Make Alg Param
		double F_w = 0.85; // DE parameter -- Make Alg Param
		bool flush = true; // Population diversity action -- Make Alg Param
		int flushStrat = 1; // Population diversity strategy -- Make Alg Param
		int flushInd = 5; // Strat 1 = kill mod index ||| Strat 2 = Index past which fitnesses are nullified -- Make Alg Param
		int flushKW = 3; // Strat 1 kill-band width -- Make Alg Param
		int flushLim = 500; // Stagnation history -- Make Alg Param
		int flushGen = 1000; // New randoms post-flush -- Make Alg Param
		int flushC = 0; // Count of flushes occured
		int flushCLim = 3; // Limit of flushes -- Make Alg Param
		//fitThresh = 0.001; // Fitness minimum change -- Make Alg Param
		//hist.resize(flushLim); // Fit history
		//for (int i = 0; i< hist.size(); i++)
		//	hist[i] = EMTG::math::LARGE;
		double flushWait = 0.1;
		bool seeded = true; // Seed rand reset with best -- Make Alg Param 
		int mutType = 2; // Mutation type employed -- Make Alg Param
		double mutProb = 0.02; // Pareto mutation probability -- Make Alg Param
		vector<int> bn (3); // Random bin index holder
		vector<int> num (3); // Random bin entry index holder
		int writeFreq = 25000; // Writing frequency in case of crash -- Make Alg Param
		vector<double> xTri (nX);
		for (int i = 0; i < 3; i++)
		{
			bn[i] = -1;
			num[i] = -1;
		}

		int i = 0;
		gCb = 0;
		mutStr = "";
		while (i < evoIter && gCb < genC)
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
					} while(deMode && (num[j]==num[(j+1)%3] || num[j]==num[(j+2)%3])); // Forces reselection if choice is duplicate -- now using duplicate bin knowledge
					//} while(deMode && ((bn[j]==bn[(j+1)%3] && num[j]==num[(j+1)%3]) || (bn[j]==bn[(j+2)%3] && num[j]==num[(j+2)%3]))); // Forces reselection if choice is duplicate
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
					xTri[j] = temp;
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
						xTri[k] = Bins[b][c][k];
				}
			}
			bool muted = false;
			switch (mutType)
			{
			case 1:
				{
					double mutUse = DoubleDistribution(RNG);
					if(mutUse < mutProb) // Pareto mutation
					{
						muted = true;
						mutStr = " with Mut";
						for (int j = 0; j < nX; j++)
						{
							double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
							int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;
							double perturb = s * varPer * r;
							double temp = xTri[j] + perturb;
							if(temp > 1.0)
							{
								if(temp < 2.0 && abs(xTri[j]-1) < Problem->options.FD_stepsize)
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
								if(temp > -1.0 && abs(xTri[j]) < Problem->options.FD_stepsize)
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
							xTri[j] = temp;
						}
					}
					break;
				}
			case 2:
				{
					for (int j = 0; j < nX; j++)
					{
						double mutUse = DoubleDistribution(RNG);
						if(mutUse < mutProb)
						{
							muted = true;
							mutStr = " with Mut";
							double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
							int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;
							double perturb = s * varPer * r;
							double temp = xTri[j] + perturb;
							if(temp > 1.0)
							{
								if(temp < 2.0 && abs(xTri[j]-1) < Problem->options.FD_stepsize)
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
								if(temp > -1.0 && abs(xTri[j]) < Problem->options.FD_stepsize)
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
							xTri[j] = temp;
						}
					}
					break;
				}
			}
			string mdStr = "_"; // Operation display string (mode string)
			if(varType < mutCap)
			{
				if(muted)
					mdStr = ":";
				else
					mdStr = ".";
			}
			else
			{
				if(muted)
					mdStr = "+";
				else
					mdStr = "-";
			}
			if(gCb == 0)
			{
				double tI = i;
				cout << "\n" << (tI/evoIter)*100 << "%\n";
			}
			cout << mdStr;
			SubmitInd(xTri, &objValTemp, &conValTemp);
			//if(proxKill)
			//	cout << "x";
			mutStr = "";
			//hist[(i%flushLim)] = Fits[0][0][0];
			gCb++;
			// !(hist[(i%flushLim)] - hist[((i+1)%flushLim)] < -1*fitThresh) // old
			if(flush && flushC < flushCLim && i > flushWait*evoIter && gCb > flushLim) // Flush procedure
			{
				if(flushC == 0)
				{
					flushLim = flushLim*2;
					optLim = optLim*2;
				}
				int delC = 0; // Deletion counter
				killDist = killDist/2;
				cout << "\nFLUSH\n";
				switch (flushStrat)
				{
				case 1: // Banded deletion
					{
						for (int j = 0; j < nF; j++)
						{
							for (int m = Bins[j].size()-1; m >= flushInd; m--)
							{
								if(m % flushInd < flushKW)
								{
									Bins[j].erase(Bins[j].begin()+m);
									Fits[j].erase(Fits[j].begin()+m);
								}
							}
							while (Bins[j].size() < NP)
							{
								Bins[j].push_back(zeros);
								Fits[j].push_back(trashFits);
							}
						}
						break;
					}
				case 2: // Back-end voiding
					{
						for (int j = 0; j < nF; j++)
						{
							for (int k = flushInd; k < NP; k++)
							{
								Fits[j][k] = trashFits;
							}
						}
					}
				}
				gCb = 0;
				flushC++;
				RandInit(flushGen,seeded);
			}
			if(i % writeFreq == 0 && i!=0)
			{
				cout << "\nsaving\n";
				for (int i = 0; i < nF; i++)
				{
					BestX = Bins[i][0];
					Problem->unscale(BestX.data());
					Problem->evaluate(Problem->X.data(), FVector.data(), Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);
					Problem->Xopt = Problem->X; //we store the unscaled Xbest
					Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "_Bin_" + to_string((_ULonglong)i) + ".emtg";
					Problem->output();
				}
			}
			i++;
		}
		RandInit(flushGen,true);
		if (!(gCb < genC))
			cout << "\nEvolution terminated due to stagnation.\n";
		else
			cout << "\n";
		/*for (int i = 0; i < pM; i++) // check phasePts for elements
			cout << "Entry " << i << " is " << phasePts[i] << "\n";*/

		for (int i = 0; i < nF; i++) // sanity check to ensure bins didn't change size
		{
			cout << "Bin " << i << " cap = " << Bins[i].size() << "\n";
			BestX = Bins[i][0];
			Problem->unscale(BestX.data());
			Problem->evaluate(Problem->X.data(), FVector.data(), Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);
			//BestObjectiveValue = 0;
			//Problem->unscale(&BestX[0]);
			Problem->Xopt = Problem->X; //we store the unscaled Xbest
			//Problem->Xopt = BestX;
			Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "_Bin_" + to_string((_ULonglong)i) + ".emtg";
			Problem->output();
		}
	}

}}//close namespaces