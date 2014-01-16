//adaptive constrained differential evolution
//a hybrid of jDE (self-adaptive differential evolution) by Brest et al. "Dynamic Optimization using Self-Adaptive Differential Evolution"
//and CDE (constrained differential evolution) by Mezura-Montes, Coello, and Tun-Morales, "Simple Feasibility Rules and Differential Evolution for Constrained Optimization"
//inspired by ACT solution to GTOC6
//created for EMTG by Jacob Englander, 5/6/2013

#include "AdaptiveConstrainedDiffEvolve.h"
#include "EMTG_math.h"

namespace EMTG { namespace Solvers {

	//default constructor (never called)
	AdaptiveConstrainedDiffEvolve::AdaptiveConstrainedDiffEvolve(void)
	{
	}

	//constructor with a problem object
	AdaptiveConstrainedDiffEvolve::AdaptiveConstrainedDiffEvolve(EMTG::problem* Problem_input)
	{
		Initialize(Problem_input);
	}

	//destructor
	AdaptiveConstrainedDiffEvolve::~AdaptiveConstrainedDiffEvolve(void)
	{
	}

	//method to initialize the solver
	void AdaptiveConstrainedDiffEvolve::Initialize(EMTG::problem* Problem_input)
	{
		//store a pointer to the problem object
		Problem = Problem_input;

		//extract the dimensions of the problem
		nX = Problem->total_number_of_NLP_parameters;
		nF = Problem->total_number_of_constraints;
		NP = 10 * nX;

		//size the storage arrays
		CurrentGenerationPopulation.clear();
		NextGenerationPopulation.clear();
		vector<double> Xdummy(nX);
		
		for (int i = 0; i < NP; ++i)
		{
			CurrentGenerationPopulation.push_back(Xdummy);
			NextGenerationPopulation.push_back(Xdummy);
		}

		CurrentGenerationObjectiveValue.resize(NP);
		NextGenerationObjectiveValue.resize(NP);
		CurrentGenerationConstraintNorm.resize(NP);
		NextGenerationConstraintNorm.resize(NP);
		CurrentGenerationAges.resize(NP);
		NextGenerationAges.resize(NP);

		TrialVector.resize(nX);
		FVector.resize(nF);
		ConstraintViolationVector.resize(nF - 1);

		//initialize the "best" values
		BestObjectiveValue = EMTG::math::LARGE;
		BestConstraintNorm = EMTG::math::LARGE;

		//initialize the adaptive coefficients
		//values from Brest et al. "Dynamic Optimization using Self-Adaptive Differential Evolution"
		Tau1 = 0.1;
		Tau2 = 0.1;
		ScaleFactorLowerBound = 1.0 / NP;
		ScaleFactorUpperBound = 0.9;

		ScaleFactor = 0.5; //RNG.rand() * (ScaleFactorUpperBound - ScaleFactorLowerBound) + ScaleFactorLowerBound;
		CrossoverRatio = 0.9; //RNG.rand();

		//age at which an individual starts being checked for "death"
		AgeLimit = 20;

		//initialize the random number generator
		DoubleDistribution = boost::uniform_real<> (0.0, 1.0);
		IntegerDistributionNP = boost::uniform_int<> (0, NP - 1);
		IntegerDistributionnX = boost::uniform_int<> (0, nX - 1);

		//seed the random number generator based on the node's clock
		RNG.seed(std::time(0));

		//finally, generate the population
		GenerateNewPopulation();
	}

	//method to generate a new population
	void AdaptiveConstrainedDiffEvolve::GenerateNewPopulation()
	{
		for (int i = 0; i < NP; ++i)
		{
			for (int j = 0; j < nX; ++j)
				CurrentGenerationPopulation[i][j] = DoubleDistribution(RNG);

			EvaluateIndividual(CurrentGenerationPopulation[i], &CurrentGenerationObjectiveValue[i], &CurrentGenerationConstraintNorm[i]);
			CurrentGenerationAges[i] = 0;

			//if the newest point is better than the known best, adopt it as the new best
			if (CurrentGenerationConstraintNorm[i] < Problem->options.snopt_feasibility_tolerance && CurrentGenerationObjectiveValue[i] < BestObjectiveValue)
			{
				BestX = CurrentGenerationPopulation[i];
				BestConstraintNorm = CurrentGenerationConstraintNorm[i];
				BestObjectiveValue = CurrentGenerationObjectiveValue[i];
				BestIndex = i;

				cout << "New point found with C = " << BestConstraintNorm << ", J = " << BestObjectiveValue << endl;

				//Problem->unscale(&BestX[0]);
				//Problem->Xopt = Problem->X; //we store the unscaled Xbest
				//Problem->Xopt = BestX;
				//Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
				//Problem->output();
			}
			else if (CurrentGenerationConstraintNorm[i] < BestConstraintNorm)
			{
				BestX = CurrentGenerationPopulation[i];
				BestConstraintNorm = CurrentGenerationConstraintNorm[i];
				BestObjectiveValue = CurrentGenerationObjectiveValue[i];
				BestIndex = i;

				cout << "New point found with C = " << BestConstraintNorm << ", J = " << BestObjectiveValue << endl;
				
				//Problem->unscale(&BestX[0]);
				//Problem->Xopt = Problem->X; //we store the unscaled Xbest
				//Problem->Xopt = BestX;
				//Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
				//Problem->output();
			}
		}
	}

	//method to evaluate one individual in the population
	void AdaptiveConstrainedDiffEvolve::EvaluateIndividual(vector<double>& X, double* ObjectiveFunctionValue, double* ConstraintNormValue)
	{
		//unscale the decision vector
		Problem->unscale(X.data());

		//set the failure flag to false
		bool failflag = false;

		//evaluate the decision vector
		try
		{
			Problem->evaluate(&(Problem->X[0]), FVector.data(), &Problem->G[0], 0, Problem->iGfun, Problem->jGvar);
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
	void AdaptiveConstrainedDiffEvolve::Evolve()
	{
		//Loop until converged
		int halt = 0;
		int Generation = 0;
		int nEvals = NP;
		int StallGeneration = 0;
		while (halt == 0)
		{
			//step 1: increment the generation counter
			++Generation;
			++StallGeneration;

			//step 2: update ScaleFactor and CrossOverRatio
			if (DoubleDistribution(RNG) < Tau1)
				ScaleFactor = ScaleFactorLowerBound + DoubleDistribution(RNG) * ScaleFactorUpperBound;

			if (DoubleDistribution(RNG) < Tau2)
				CrossoverRatio = DoubleDistribution(RNG);

			//step 3: iterate through the population
			for (int i = 0; i < NP; ++i)
			{
				//step 3.1: pick three random indices from the population
				int r1, r2, r3;
				
				r1 = IntegerDistributionNP(RNG);

				while(true)
				{
					r2 = IntegerDistributionNP(RNG);

					if (!(r2 == r1))
						break;
				}
				
				while(true)
				{
					r3 = IntegerDistributionNP(RNG);
					if (!(r3 == r2) && !(r3 == r1))
						break;
				}

				//step 3.2 construct the trial vector
				int randj = IntegerDistributionnX(RNG);
				for (int j = 0; j < nX; ++j)
				{
					if (DoubleDistribution(RNG) < CrossoverRatio || j == randj)
						TrialVector[j] = CurrentGenerationPopulation[r3][j] + ScaleFactor * (CurrentGenerationPopulation[r1][j] - CurrentGenerationPopulation[r2][j]);
					else
						TrialVector[j] = CurrentGenerationPopulation[r3][j];
				}

				//step 3.3 evaluate the trial vector
				EvaluateIndividual(TrialVector, &TrialObjectiveValue, &TrialConstrainNorm);
				++nEvals;

				//step 3.4 compare the trial vector to the current element of the population
				
				//check feasibility of the new solution
				if (TrialConstrainNorm > Problem->options.snopt_feasibility_tolerance)
				{
					//if the trial point is infeasible, check the current point
					if (CurrentGenerationConstraintNorm[i] > Problem->options.snopt_feasibility_tolerance)
					{
						//if both solutions are infeasible then discard the one which is least feasible
						if (TrialConstrainNorm < CurrentGenerationConstraintNorm[i])
						{
							NextGenerationPopulation[i] = TrialVector;
							NextGenerationConstraintNorm[i] = TrialConstrainNorm;
							NextGenerationObjectiveValue[i] = TrialObjectiveValue;
							

							if (TrialConstrainNorm < BestConstraintNorm)
							{
								if ((BestConstraintNorm - TrialConstrainNorm) / BestConstraintNorm > 0.00001)
									StallGeneration = 0;

								BestX = TrialVector;
								BestConstraintNorm = TrialConstrainNorm;
								BestObjectiveValue = TrialObjectiveValue;
								BestIndex = i;

								StallGeneration = 0;

								cout << "New point found with C = " << BestConstraintNorm << ", J = " << BestObjectiveValue << endl;

								//Problem->unscale(&BestX[0]);
								//Problem->Xopt = Problem->X; //we store the unscaled Xbest
								//Problem->Xopt = BestX;
								//Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
								//Problem->output();
							}
						}
						//otherwise you don't have to do anything and the current member of the population is retained
						//however we need to check to see if the member has "aged out"
						else
						{
							if (!(BestIndex == i) && CurrentGenerationAges[i] > AgeLimit)
							{
								if (DoubleDistribution(RNG) < 0.1 || (CurrentGenerationConstraintNorm[i] - BestConstraintNorm < 0.0001))
								{
									for (int j = 0; j < nX; ++j)
										NextGenerationPopulation[i][j] = DoubleDistribution(RNG);

									EvaluateIndividual(NextGenerationPopulation[i], &NextGenerationObjectiveValue[i], &NextGenerationConstraintNorm[i]);
									NextGenerationAges[i] = 0;
								}
							}
							else
								NextGenerationAges[i] = CurrentGenerationAges[i] + 1;
						}
					}
				}
				else
				{
					//check to see if the current point is feasible
					if (CurrentGenerationConstraintNorm[i] > Problem->options.snopt_feasibility_tolerance)
					{
						//if the current point is infeasible and the trial point is feasible, retain the trial point and discard the current point
						NextGenerationPopulation[i] = TrialVector;
						NextGenerationConstraintNorm[i] = TrialConstrainNorm;
						NextGenerationObjectiveValue[i] = TrialObjectiveValue;
					}
					else
					{
						//both points are feasible, so it is necessary to compare their objective function values
						if (TrialObjectiveValue < CurrentGenerationObjectiveValue[i])
						{
							NextGenerationPopulation[i] = TrialVector;
							NextGenerationConstraintNorm[i] = TrialConstrainNorm;
							NextGenerationObjectiveValue[i] = TrialObjectiveValue;


							if (TrialObjectiveValue < BestObjectiveValue)
							{
								if ((BestObjectiveValue - TrialObjectiveValue) / BestObjectiveValue > 1.0e-3)
									StallGeneration = 0;

								BestX = TrialVector;
								BestConstraintNorm = TrialConstrainNorm;
								BestObjectiveValue = TrialObjectiveValue;
								BestIndex = i;

								cout << "New point found with C = " << BestConstraintNorm << ", J = " << BestObjectiveValue << endl;

								Problem->unscale(&BestX[0]);
								Problem->Xopt = Problem->X; //we store the unscaled Xbest
								Problem->Xopt = BestX;
								Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
								Problem->output();
							}
						}
						//otherwise you don't have to do anything and the current member of the population is retained
						//however we need to check to see if the member has "aged out"
						else
						{
							if (!(BestIndex == i) && CurrentGenerationAges[i] > AgeLimit)
							{
								if (DoubleDistribution(RNG) < 0.1  || (CurrentGenerationObjectiveValue[i] - BestObjectiveValue < 0.0001))
								{
									for (int j = 0; j < nX; ++j)
										NextGenerationPopulation[i][j] = DoubleDistribution(RNG);

									EvaluateIndividual(NextGenerationPopulation[i], &NextGenerationObjectiveValue[i], &NextGenerationConstraintNorm[i]);
									NextGenerationAges[i] = 0;
								}
							}
							else
								NextGenerationAges[i] = CurrentGenerationAges[i] + 1;
						}
					}
				}
			}

			//step 4: the current population is replaced by the new population
			CurrentGenerationPopulation = NextGenerationPopulation;
			CurrentGenerationObjectiveValue = NextGenerationObjectiveValue;
			CurrentGenerationConstraintNorm = NextGenerationConstraintNorm;
			CurrentGenerationAges = NextGenerationAges;
			
			//step 5: check convergence
			if (Generation == 1e+4)
			{
				cout << "ACDE halted due to number of generations reaching " << Generation << endl;
				halt = 1;
			}
			else if (nEvals == 1e+6)
			{
				cout << "ACDE halted due to number of function evaluations reaching " << nEvals << endl;
				halt = 2;
			}
			else if (StallGeneration == 1000)
			{
				cout << "ACDE halted due to number of stall generations reaching " << StallGeneration << endl;
				halt = 3;
			}
		}

		Problem->unscale(&BestX[0]);
		Problem->Xopt = Problem->X; //we store the unscaled Xbest
		Problem->Xopt = BestX;
		Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
		Problem->output();
	}

}}//close namespaces