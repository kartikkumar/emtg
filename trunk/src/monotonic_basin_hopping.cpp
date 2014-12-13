//Monotonic Basin Hopping
//for EMTG version 8
//Jacob Englander 7-27-2012
#include <iostream>
#include <fstream>
#include <ctime>

#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"

#include "problem.h"
#include "monotonic_basin_hopping.h"
#include "EMTG_math.h"
#include "SNOPT_user_function.h"

#include "snoptProblemExtension.h"

#ifdef _use_WORHP
#include "EMTG_WORHP_interface.h"
#endif

using namespace std;

namespace EMTG { namespace Solvers {
	//constructors
	MBH::MBH() {}

	MBH::MBH(EMTG::problem* Problem_input)
	{
		//initialize the MBH variables
		initialize(Problem_input);

		//create the dynamic SNOPT variables
		neF = Problem->total_number_of_constraints;

		lenA  = Problem->total_number_of_NLP_parameters*Problem->total_number_of_constraints;

		iAfun = new SNOPT_INT_TYPE[lenA];
		jAvar = new SNOPT_INT_TYPE[lenA];
		A  = new SNOPT_DOUBLE_TYPE[lenA];

		lenG   = Problem->total_number_of_NLP_parameters*Problem->total_number_of_constraints;
		iGfun = new SNOPT_INT_TYPE[lenG];
		jGvar = new SNOPT_INT_TYPE[lenG];

		x      = new SNOPT_DOUBLE_TYPE[Problem->total_number_of_NLP_parameters];
		xlow   = new SNOPT_DOUBLE_TYPE[Problem->total_number_of_NLP_parameters];
		xupp   = new SNOPT_DOUBLE_TYPE[Problem->total_number_of_NLP_parameters];
		xmul   = new SNOPT_DOUBLE_TYPE[Problem->total_number_of_NLP_parameters];
		xstate = new    SNOPT_INT_TYPE[Problem->total_number_of_NLP_parameters];

		F      = new SNOPT_DOUBLE_TYPE[neF];
		Flow   = new SNOPT_DOUBLE_TYPE[neF];
		Fupp   = new SNOPT_DOUBLE_TYPE[neF];
		Fmul   = new SNOPT_DOUBLE_TYPE[neF];
		Fstate = new SNOPT_INT_TYPE[neF];

#ifdef Heritage_SNOPT7
		nxnames = 1;
		nFnames = 1;
		xnames = new char[nxnames*8];
		Fnames = new char[nFnames*8];
#endif

		ObjRow = 0;
		ObjAdd = 0;

		//initialize a few more SNOPT variables whose quantities do not change with an MBH run
		for (int k=0; k < Problem->total_number_of_NLP_parameters; ++k)
		{
			xlow[k] = 0.0;
			xupp[k] = 1.0;
		}
		for (int k=0; k < Problem->total_number_of_constraints; ++k)
		{
			Flow[k] = Problem->Flowerbounds[k];
			Fupp[k] = Problem->Fupperbounds[k];
		}

		//we don't have any real intuition for what the intial "states" of the F functions will be
		//the same applies for the constraint Lagrange multipliers, so set all of the Fstates = 0 and Fmuls to 0.0
		for (int k=0; k < Problem->total_number_of_constraints; ++k)
		{
			Fstate[k] = 0;
			Fmul[k] = 0.0;
		}

		this->SNOPTproblem = new snoptProblemExtension(true);

		this->SNOPTproblem->setProblemSize(Problem->total_number_of_NLP_parameters, neF);
		this->SNOPTproblem->setObjective(ObjRow, ObjAdd);
		this->SNOPTproblem->setUserspace((SNOPT_INT_TYPE*)&SNOPT_start_time, 500, (SNOPT_DOUBLE_TYPE*)Problem, 500);
		this->SNOPTproblem->setA(lenA, iAfun, jAvar, A);
		this->SNOPTproblem->setG(lenG, iGfun, jGvar);
		//this->SNOPTproblem->setIntParameter("Total real workspace", lenA + lenG);
		//this->SNOPTproblem->setIntParameter("Total integer workspace", lenA + lenG);
		//this->SNOPTproblem->setIntParameter("Total character workspace", lenA + lenG);
#ifdef Heritage_SNOPT7
		this->SNOPTproblem->setXNames(xnames, nxnames);
		this->SNOPTproblem->setFNames(Fnames, nFnames);
#endif
		this->SNOPTproblem->setProbName("EMTG");
		this->SNOPTproblem->setUserFun(SNOPT_user_function);
		this->SNOPTproblem->setIntParameter("Iterations limit", 100 * Problem->options.snopt_major_iterations);
		this->SNOPTproblem->setIntParameter("Major iterations limit", Problem->options.snopt_major_iterations);
        //if (Problem->options.mission_type == 5)
            //this->SNOPTproblem->setIntParameter("Scale option", 2);
        if (Problem->options.derivative_type == 4)
            this->SNOPTproblem->setIntParameter("Derivative option", 1);
        else
		    this->SNOPTproblem->setIntParameter("Derivative option", 0);
		this->SNOPTproblem->setIntParameter("Minor print level", 0);
		this->SNOPTproblem->setRealParameter("Major feasibility tolerance", Problem->options.snopt_feasibility_tolerance);
		
		this->SNOPTproblem->setIntParameter("Major Print Level", 1);
		this->SNOPTproblem->setRealParameter("Optimality tolerance", 1.0e-6);
		if (Problem->options.check_derivatives)
		{
			this->SNOPTproblem->setIntParameter("Print file", 1);
			this->SNOPTproblem->setIntParameter("Summary file", 1);
            this->SNOPTproblem->setIntParameter("System information", 1);
			this->SNOPTproblem->setIntParameter("Verify level", 3); //0 = cheap test 1 = individual gradients checked (OK or BAD) 2 = Individual columns of the Jacobian are checked 3 = 1 and 2 happen -1 = Derivative checking is disabled
		}
		if (Problem->options.quiet_NLP)
		{
			this->SNOPTproblem->setIntParameter("Major Print Level", 0);
			this->SNOPTproblem->setIntParameter("Minor print level", 0);
			this->SNOPTproblem->setIntParameter("Print No", 0);
			this->SNOPTproblem->setIntParameter("Summary file", 0);
			this->SNOPTproblem->setParameter("Suppress parameters");
		}

		if (Problem->options.NLP_solver_mode)
			this->SNOPTproblem->setParameter("Minimize");
		else
			this->SNOPTproblem->setParameter("Feasible point");

		//set up the random number generator
		DoubleDistribution = boost::uniform_real<>(0.0, 1.0);

		//search through the problem object and identify which decision variables are flight time variables
		if (Problem->options.MBH_time_hop_probability > 0.0)
		{
			for (int entry = 0; entry < Problem->total_number_of_NLP_parameters; ++entry)
				if ( Problem->Xdescriptions[entry].find("flight time") < 1024)
				{
					this->time_variable_indices.push_back(entry);
				}
		}

		//search through the problem object and identify which decision variables are significant (i.e. non-control)
		if (Problem->options.MBH_time_hop_probability > 0.0)
		{
			for (int entry = 0; entry < Problem->total_number_of_NLP_parameters; ++entry)
				if ( !(Problem->Xdescriptions[entry].find("u_x") < 1024 || Problem->Xdescriptions[entry].find("u_y") < 1024 || Problem->Xdescriptions[entry].find("u_z") < 1024 || Problem->Xdescriptions[entry].find("Isp") < 1024) )
				{
					this->significant_variable_indices.push_back(entry);
				}
		}
	}

	//destructor
	MBH::~MBH()
	{
		//delete the SNOPT object
		delete SNOPTproblem;

		//delete all temporary SNOPT arrays
		delete [] iAfun;
		delete [] jAvar;
		delete [] A;

		delete [] iGfun;
		delete [] jGvar;
		
		delete [] x;
		delete [] xlow;
		delete [] xupp;
		delete [] xmul;
		delete [] xstate;

		delete [] F;
		delete [] Flow;
		delete [] Fupp;
		delete [] Fmul;
		delete [] Fstate;

#ifdef Heritage_SNOPT7
		delete [] xnames;
		delete [] Fnames;
#endif
	}

	
	//method to initialize the MBH solver
	//resets all storage fields
	int MBH::initialize(EMTG::problem* Problem_input)
	{
		Problem = Problem_input;


		//size the storage vectors
		Xtrial_scaled.resize(Problem->total_number_of_NLP_parameters, 0);
		Xlocalbest_scaled.resize(Problem->total_number_of_NLP_parameters, 0);
		Xcurrent_scaled.resize(Problem->total_number_of_NLP_parameters, 0);
		Xbest_scaled.resize(Problem->total_number_of_NLP_parameters, 0);
		Xbest.resize(Problem->total_number_of_NLP_parameters, 0);

		archive.clear();
		archive_scores.clear();

		//clear the scores
		ftrial = EMTG::math::LARGE;
		fcurrent = EMTG::math::LARGE;
		fbest = EMTG::math::LARGE;
		most_feasible = EMTG::math::LARGE;
		
		//clear the counters
		number_of_solutions = 0;
		number_of_improvements = 0;
		number_of_failures_since_last_improvement = 0;
		number_of_resets = 0;

		//reset the Jacobian computed flag
		computed_Jacobian = false;
		jacfullrankflag = true;

		//reset the Jacobian printed flag
		this->printed_sparsity = false;

		//seed the random number generator based on the node's clock
		RNG.seed(time(0));

		return 0;
	}

	//function to reset to a new point, either at the beginning of the problem or when a certain number of failures is reached
	int MBH::reset_point()
	{
		//reset the number of failures
		this->number_of_failures_since_last_improvement = 0;

		//increment the number of global resets
		++this->number_of_resets;

		//generate a new random trial point
		if (Problem->options.MBH_zero_control_initial_guess > 0)
		{
			for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
				this->Xtrial_scaled[k] = 0.5;
			for (size_t sigk = 0; sigk < this->significant_variable_indices.size(); ++sigk)
				this->Xtrial_scaled[this->significant_variable_indices[sigk]] = DoubleDistribution(RNG);
		}
		else
		{
			for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
				this->Xtrial_scaled[k] = DoubleDistribution(RNG);
		}

		fcurrent = EMTG::math::LARGE;

		return 0;
	}

	//function to seed a new point based on an input trial point
	int MBH::seed(vector<double> seed_vector)
	{
		//reset the number of failures
		number_of_failures_since_last_improvement = 0;

		//read a trial point
		if (seed_vector.size() > Xtrial_scaled.size())
		{
			cout << "Seed vector is longer than problem decision vector. Truncating." << endl;
		}
		for (size_t k = 0; k < Xtrial_scaled.size(); ++k)
			Xtrial_scaled[k] = (seed_vector[k] - Problem->Xlowerbounds[k]) / (Problem->Xupperbounds[k] - Problem->Xlowerbounds[k]);
		if (seed_vector.size() < Problem->total_number_of_NLP_parameters)
		{
			cout << "Seed vector is shorter than problem decision vector. Extending with random values." << endl;
			for (size_t k = 0; k < Problem->total_number_of_NLP_parameters - seed_vector.size(); ++k)
				Xtrial_scaled[k] = DoubleDistribution(RNG);
		}

        this->Problem->evaluate(seed_vector.data(), this->F, Problem->G.data(), 0, this->Problem->iGfun, this->Problem->jGvar);
        double Ccurrent = this->check_feasibility();

		fcurrent = EMTG::math::LARGE;

		return 0;
	}

	//function to perform a "hop" operation
	int MBH::hop()
	{
		if (Problem->options.MBH_hop_distribution == 0)
		{
			//perform a uniform "hop"
			if (Problem->options.MBH_zero_control_initial_guess > 1)
			{
				for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
					this->Xtrial_scaled[k] = 0.5;
				for (int sigk = 0; sigk < this->significant_variable_indices.size(); ++sigk)
				{
					double r = DoubleDistribution(RNG);

					step_size = 2 * (r - 0.5) * Problem->options.MBH_max_step_size;

					this->Xtrial_scaled[this->significant_variable_indices[sigk]] = Xcurrent_scaled[this->significant_variable_indices[sigk]] + step_size;
				}
			}
			else
			{
				for (int k=0; k < Problem->total_number_of_NLP_parameters; ++k)
				{
					double r = DoubleDistribution(RNG);

					step_size = 2 * (r - 0.5) * Problem->options.MBH_max_step_size;

					Xtrial_scaled[k] = Xcurrent_scaled[k] + step_size;
				}
			}
		}
		else if (Problem->options.MBH_hop_distribution == 1)
		{
			//perform a Cauchy "hop"
			if (Problem->options.MBH_zero_control_initial_guess > 1)
			{
				for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
					this->Xtrial_scaled[k] = 0.5;
				for (int sigk = 0; sigk < this->significant_variable_indices.size(); ++sigk)
				{
					//Cauchy generator v2 as per instructions from Arnold Englander 1-12-2014 using Cauchy CDF

					double r = 0.0;

					while (r < math::SMALL)
						r = DoubleDistribution(RNG);

					step_size = Problem->options.MBH_max_step_size * tan(math::PI*(r - 0.5));

					this->Xtrial_scaled[this->significant_variable_indices[sigk]] = Xcurrent_scaled[this->significant_variable_indices[sigk]] + step_size;
				}
			}
			else
			{
				for (int k=0; k < Problem->total_number_of_NLP_parameters; ++k)
				{
					//Cauchy generator v2 as per instructions from Arnold Englander 1-12-2014 using Cauchy CDF

					double r = 0.0;

					while (r < math::SMALL)
						r = DoubleDistribution(RNG);

					step_size = Problem->options.MBH_max_step_size * tan(math::PI*(r - 0.5));

					Xtrial_scaled[k] = Xcurrent_scaled[k] + step_size;
				}
			}
		}
		else if (Problem->options.MBH_hop_distribution == 2)
		{
			double alpha = Problem->options.MBH_Pareto_alpha;
			//perform a Pareto hop
			if (Problem->options.MBH_zero_control_initial_guess > 1)
			{
				for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
					this->Xtrial_scaled[k] = 0.5;
				for (int sigk = 0; sigk < this->significant_variable_indices.size(); ++sigk)
				{
					//Pareto distribution from x = 1.0 with alpha, modified to start from x = 0.0

					//s*(((alpha-1)/xU_min)/power((xU_min/u),-alpha))
					double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
					int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;

					step_size = s * Problem->options.MBH_max_step_size * r;

					this->Xtrial_scaled[this->significant_variable_indices[sigk]] = Xcurrent_scaled[this->significant_variable_indices[sigk]] + step_size;
				}
			}
			else
			{
				for (int k=0; k < Problem->total_number_of_NLP_parameters; ++k)
				{
					//Pareto distribution from x = 1.0 with alpha, modified to start from x = 0.0

					//s*(((alpha-1)/xU_min)/power((xU_min/u),-alpha))
					double r = ( (alpha - 1.0) / math::SMALL ) / pow((math::SMALL / (math::SMALL + DoubleDistribution(RNG))), -alpha);
					int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;

					step_size = s * Problem->options.MBH_max_step_size * r;

					Xtrial_scaled[k] = Xcurrent_scaled[k] + step_size;
				}
			}
		}
		else if (Problem->options.MBH_hop_distribution == 3)
		{
			//perform a Gaussian hop
			double sigma = Problem->options.MBH_max_step_size;
			double sigma2 = sigma * sigma;
			if (Problem->options.MBH_zero_control_initial_guess > 1)
			{
				for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
					this->Xtrial_scaled[k] = 0.5;
				for (int sigk = 0; sigk < this->significant_variable_indices.size(); ++sigk)
				{
					double r = DoubleDistribution(RNG);
					int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;

					step_size = s / (sigma * sqrt(math::TwoPI)) * exp(-r*r / (2*sigma2));

					this->Xtrial_scaled[this->significant_variable_indices[sigk]] = Xcurrent_scaled[this->significant_variable_indices[sigk]] + step_size;
				}
			}
			else
				{
				for (int k=0; k < Problem->total_number_of_NLP_parameters; ++k)
				{
					double r = DoubleDistribution(RNG);
					int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;

					step_size = s / (sigma * sqrt(math::TwoPI)) * exp(-r*r / (2*sigma2));

					Xtrial_scaled[k] = Xcurrent_scaled[k] + step_size;
				}
			}
		}

		//MBH clipping
        for (size_t k = 0; k < Xtrial_scaled.size(); ++k)
        {
            if (Xtrial_scaled[k] > 1.0)
                Xtrial_scaled[k] = 1.0;
            else if (Xtrial_scaled[k] < 0.0)
                Xtrial_scaled[k] = 0.0;
        }


		return 0;
	}

	//function to perform a "time hop" operation
	int MBH::time_hop()
	{
		//loop through any time variables and if (uniform random < threshold) then add/subtract a synodic period
		for (int timeindex = 0; timeindex < time_variable_indices.size(); ++timeindex)
		{
            if (DoubleDistribution(RNG) < Problem->options.MBH_time_hop_probability)
			{
				int k = time_variable_indices[timeindex];
				int s = DoubleDistribution(RNG) > 0.5 ? 1 : -1;
				Xtrial_scaled[k] = Xcurrent_scaled[k] + s * Problem->synodic_periods[timeindex] / (Problem->Xupperbounds[k] - Problem->Xlowerbounds[k]);
			}
		}

		return 0;
	}

	//function to perform a "slide" operation, i.e. run SNOPT
	int MBH::slide()
	{
		//Step 1: set the current state equal to the initial guess
		for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
		{
			xstate[k] = 0;
			xmul[k] = 0.0;
			x[k] = Xtrial_scaled[k];
		}

		//print the sparsity file and XF files if this is the first pass, otherwise don't to save time and hard drive cycles
		if (!this->printed_sparsity)
		{
			this->printed_sparsity = true;
			if (!Problem->options.quiet_basinhopping)
			{
				Problem->output_Jacobian_sparsity_information(Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "_SparsityDescriptions.csv");
				Problem->output_problem_bounds_and_descriptions(Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "XFfile.csv");
			}
		}

		if (this->Problem->options.NLP_solver_type == 1)
		{
			//run WORHP
#ifdef _use_WORHP
			this->WORHP_interface = new EMTG_WORHP_interface(Problem);
			this->WORHP_interface->SetInitialGuess(Xtrial_scaled.data());
			this->WORHP_interface->Solve();
			delete this->WORHP_interface;
#else
			cout << "WORHP not installed" << endl;
#endif
		}
		else
		{
			//run SNOPT
			for (int k=0; k < Problem->total_number_of_constraints; ++k)
			{
				Fstate[k] = 0;
				Fmul[k] = 0.0;
				F[k] = 0.0;
			}
			//Step 2: attempt to calculate the Jacobian
			SNOPTproblem->setX          ( x, xlow, xupp, xmul, xstate );
			SNOPTproblem->setF          ( F, Flow, Fupp, Fmul, Fstate );

			if (!computed_Jacobian || !jacfullrankflag)
			{
				vector<bool> cjacflag(neF);
				if (true /*Problem->options.derivative_type > 0*/)
				{
					if (!Problem->options.quiet_basinhopping)
						Problem->output_Jacobian_sparsity_information(Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "_SparsityDescriptions.csv");

					for (size_t entry = 0; entry < Problem->Adescriptions.size(); ++entry)
					{
						iAfun[entry] = Problem->iAfun[entry];
						jAvar[entry] = Problem->jAvar[entry];
						A[entry] = Problem->A[entry];
					}
					for (size_t entry = 0; entry < Problem->Gdescriptions.size(); ++entry)
					{
						iGfun[entry] = Problem->iGfun[entry];
						jGvar[entry] = Problem->jGvar[entry];
					}
				
					SNOPTproblem->setNeA(Problem->Adescriptions.size());
					SNOPTproblem->setNeG(Problem->Gdescriptions.size());

					computed_Jacobian = true;
				}
				else
				{
					
					inform = SNOPTproblem->computeJac    ();
					if (!Problem->options.quiet_basinhopping)
					{
						write_sparsity("sparsitySNJac.txt");
						Problem->output_Jacobian_sparsity_information(Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "_SparsityDescriptions_SNJAC.csv");
					}

					//Step 3: If the Jacobian calculation succeeded, then check if it is full rank
					if (inform == 102)
					{
						computed_Jacobian = true;
			
						for (int k=0; k<neF; ++k)
							cjacflag[k] = false;

						for (int i=0; i<SNOPTproblem->getNeG(); ++i)
							cjacflag[iGfun[i]-1] = true;

						for (int i=0; i<SNOPTproblem->getNeA(); ++i)
							cjacflag[iAfun[i]-1] = true;
			
						jacfullrankflag = true;
						for (int i=0; i<neF-1; ++i)
						{
							if (!cjacflag[i])
								jacfullrankflag = false;
						}
					}
				}	
			}

			//If the Jacobian is full rank, then run SNOPT
			if (jacfullrankflag)
			{
				SNOPT_start_time = time(NULL);
				inform = SNOPTproblem->solve( 0 );
			}
			else if (!Problem->options.quiet_basinhopping)
			{
				cout << "Jacobian is not full rank. SNOPT cannot be run." << endl;
			}
		}

		//Step 4: store the results
		for (int k = 0; k < Problem->total_number_of_NLP_parameters; ++k)
		{
			Xtrial_scaled[k] = x[k];
		}

		Problem->unscale(x);

		if (computed_Jacobian || this->Problem->options.NLP_solver_type == 1)
		{
			try
			{
				Problem->evaluate(Problem->X.data(), this->F, Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);
			}
			catch (int errorcode) //integration step error
			{
				if (!Problem->options.quiet_basinhopping)
					cout << "EMTG::Invalid initial point or failure in objective function." << endl;
				F[0] = EMTG::math::LARGE;
			}
		}
		else
		{
			if (!Problem->options.quiet_basinhopping)
				cout << "EMTG::Jacobian was not computed successfully." << endl;
			F[0] = EMTG::math::LARGE;
		}

		for (int k = 0; k < Problem->total_number_of_constraints; ++k)
		{
			Problem->F[k] = F[k];
		}

		return (int) jacfullrankflag;
	}

	
	//function to run MBH
	int MBH::run()
	{
		bool new_point;
		bool seeded_step = false;
		double best_feasibility = math::LARGE;

		//If we have seeded MBH, start from that seed. Otherwise we will need to generate a new random point.
		if ( !(Problem->options.seed_MBH) )
			new_point = true;
		else
		{
			seeded_step = true;
			new_point = false;
		}

		int number_of_attempts = 1;
		bool continue_flag = true;

		//print the archive header
		string archive_file = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + "archive.emtg_archive";
		if (!Problem->options.quiet_basinhopping)
			print_archive_header(archive_file);

		time_t tstart = time(NULL);

		do
		{
			++number_of_attempts;
			if (new_point)
			{
				//Step 1: generate a random new point
				this->reset_point();

				this->number_of_failures_since_last_improvement = 0;
			}
			else if (!seeded_step)
			{
				//Step 1 (alternate): perturb the existing point
				this->hop();

                if (Problem->options.MBH_time_hop_probability > 0.0  && best_feasibility >= this->Problem->options.snopt_feasibility_tolerance)
					this->time_hop();
			}
			
			//if seeding MBH, only the first step runs from the seed. After that hopping occurs.
			seeded_step = false;

			//Step 2: apply the slide operator
			//Step 2.1: If we are in two-step mode, make the finite differencing step coarse
			int temporary_derivative_code;
			if (Problem->options.MBH_two_step)
			{
				if (!Problem->options.quiet_basinhopping)
					cout << "Performing coarse optimization step" << endl;
				SNOPTproblem->setRealParameter("Difference interval", Problem->options.FD_stepsize_coarse);
				temporary_derivative_code = Problem->options.derivative_type;
				Problem->options.derivative_type = 1;
			}

			int ransnopt = slide();

			//Step 3: determine if the new trial point is feasible and if so, operate on it
			double feasibility = this->check_feasibility();

			//Step 3.01 if this is an MGA or MGA-DSM problem, make the finite differencing step tight and run again
			if (Problem->options.MBH_two_step && (inform <= 3 || feasibility < Problem->options.snopt_feasibility_tolerance))
			{
				if (!Problem->options.quiet_basinhopping)
				{
					cout << "Coarse optimization step succeeded with J = " << F[0] << endl;
					cout << "Performing fine optimization step" << endl;
				}
				SNOPTproblem->setRealParameter("Difference interval", Problem->options.FD_stepsize);
				Problem->options.derivative_type = temporary_derivative_code;
				vector<double> coarseX = this->Xtrial_scaled;
				double coarseF = this->F[0];
				slide();

				feasibility = check_feasibility();

				if (inform <= 3 || feasibility < Problem->options.snopt_feasibility_tolerance)
				{
					if (!Problem->options.quiet_basinhopping)
						if (!Problem->options.quiet_basinhopping)
							cout << "Fine optimization step succeeded with J = " << F[0] << endl;
					if (this->F[0] > coarseF)
					{
						if (!Problem->options.quiet_basinhopping)
							cout << "Coarse solution was better, retaining coarse solution." << endl;
						this->Xtrial_scaled = coarseX;
						try
						{
							Problem->evaluate(this->Xtrial_scaled.data(), this->F, Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);
						}
						catch (int e)
						{
							if (!Problem->options.quiet_basinhopping)
								std::cout << "Failure to evaluate " << Problem->options.description << std::endl;
						}
					}
				}
			}

			//note: I do not trust SNOPT's "requested accuracy could not be achieved" return state - I prefer my own feasibility check
			if ((inform <= 2 && inform >= 0)|| feasibility < Problem->options.snopt_feasibility_tolerance)
			{
				//Step 3.1: if the trial point is feasible, add it to the archive
				Problem->unscale(Xtrial_scaled.data());
				archive.push_back(Problem->X);
				archive_scores.push_back(F[0]);
				archive_timestamps.push_back(time(NULL) - tstart);
				archive_step_count.push_back(number_of_attempts);
				archive_reset_count.push_back(number_of_resets);

				if (!Problem->options.quiet_basinhopping)
				{
					print_archive_line(archive_file, number_of_solutions);
					cout << "Hop evaluated mission " << Problem->options.description << " with fitness " << F[0] << endl;
				}

				++number_of_solutions;

				//Step 3.2: if the point came from a hop operation and is superior to the current point, adopt it as the new current point
				if (!(new_point))
				{
					if (F[0] < fcurrent)
					{
						fcurrent = F[0];
						Xcurrent_scaled = Xtrial_scaled;

						if (!Problem->options.quiet_basinhopping)
							cout << "New local best" << endl;

						++number_of_improvements;

						number_of_failures_since_last_improvement = 0;
					}
					else
						++number_of_failures_since_last_improvement;
				}
				else //if this is from a reset, adopt it as the current point and turn off the "new point" flag
				{
					fcurrent = F[0];
					Xcurrent_scaled = Xtrial_scaled;
					new_point = false;
				}

				//Step 3.3 update the global best if applicable
				if (F[0] < fbest)
				{
					fbest = F[0];
					Xbest_scaled = this->Xtrial_scaled; 
					Problem->unscale(this->Xbest_scaled.data());
					Problem->Xopt = Problem->X; //we store the unscaled Xbest
					Problem->best_cost = fbest;
					

					if (!Problem->options.quiet_basinhopping)
						cout << "New global best" << endl;

					//Write out a results file for the current global best
					try
					{
						Problem->evaluate(Problem->X.data(), this->F, Problem->G.data(), 0, Problem->iGfun, Problem->jGvar);
					}
					catch (int errorcode)
					{
						std::cout << "Failure to evaluate " << Problem->options.description << std::endl;
						F[0] = EMTG::math::LARGE;
					}
					Problem->options.outputfile = Problem->options.working_directory + "//" + Problem->options.mission_name + "_" + Problem->options.description + ".emtg";
					Problem->output();
				}

			}
			else if (this->Problem->options.ACE_feasible_point_finder && best_feasibility >= this->Problem->options.snopt_feasibility_tolerance)
			{
				//if we have not yet found our first feasible point and the ACE feasible point finder is enabled
				//then we should see if this point is "more feasible" than best one we have so far
				if (feasibility < best_feasibility)
				{
					if (!Problem->options.quiet_basinhopping)
						std::cout << "Acquired slightly less infeasible point with feasibility " << feasibility << std::endl;
					fcurrent = F[0];
					Xcurrent_scaled = Xtrial_scaled;
					best_feasibility = feasibility;
					new_point = false;
				}
				++number_of_failures_since_last_improvement;
			}
			else
			{
				++number_of_failures_since_last_improvement;

				if (feasibility < most_feasible)
				{
					most_feasible = feasibility;
					Problem->unscale(x);
					closest_to_feasible_solution = Problem->X;
				}
			}

			if (number_of_failures_since_last_improvement >= Problem->options.MBH_max_not_improve)
				new_point = true;

		} while (number_of_attempts < Problem->options.MBH_max_trials && (time(NULL) - tstart) < Problem->options.MBH_max_run_time);

		if (!Problem->options.quiet_basinhopping)
		{
			cout << endl;
			if (number_of_solutions > 0)
				cout << "Best value found was " << fbest << endl;
			else
				cout << "No feasible solutions found." << endl;
		}

		return number_of_solutions;
	}

	//function to print the archive header to a text file
	int MBH::print_archive_header(string filename)
	{
		//print the archive in CSV format
		//the final entry in each line is the fitness value

		ofstream outputfile (filename.c_str(), ios::trunc);

		//archive column headers
		for (int entry = 0; entry < Problem->total_number_of_NLP_parameters; ++entry)
			outputfile << Problem->Xdescriptions[entry] << ",";
		outputfile << "reset count,";
		outputfile << "step count,";
		outputfile << "solution timestamp,";
		outputfile << "Objective function" << endl;



		outputfile.close();

		return 0;
	}

	//function to print the archive line to a text file
	int MBH::print_archive_line(string filename, int linenumber)
	{
		//print the archive in CSV format
		//the final entry in each line is the fitness value

		ofstream outputfile (filename.c_str(), ios::app);

		//archive lines
		for (int entry = 0; entry < Problem->total_number_of_NLP_parameters; ++entry)
		{
			if (this->Problem->Xdescriptions[entry].find("epoch") < 1024 || this->Problem->Xdescriptions[entry].find("time") < 1024)
			{
				outputfile << this->archive[linenumber][entry] / 86400.0 << ",";
			}
			else
				outputfile << this->archive[linenumber][entry] << ",";
		}
		outputfile << archive_reset_count[linenumber] << ",";
		outputfile << archive_step_count[linenumber] << ",";
		outputfile << archive_timestamps[linenumber] << ",";
		outputfile << archive_scores[linenumber] << endl;

		outputfile.close();

		return number_of_solutions;
	}

	//function to write out the Jacobian sparsity
	int MBH::write_sparsity(string filename)
	{
		ofstream sparsityfile(filename.c_str(), ios::trunc);

		for (int k = 0; k < SNOPTproblem->getNeA(); ++k)
			sparsityfile << "Linear, " << iAfun[k] << "," << jAvar[k] << endl;

		for (int k = 0; k < SNOPTproblem->getNeG(); ++k)
			sparsityfile << "Nonlinear, " << iGfun[k] << "," << jGvar[k] << endl;

		sparsityfile.close();

		return 0;
	}

	//function to check feasibility of a solution
	double MBH::check_feasibility()
	{
		double max_constraint_violation = 0.0;
		int worst_constraint;
		for (int k = 1; k < Problem->total_number_of_constraints; ++k)
		{
			if (F[k] > Fupp[k])
			{
				if (F[k] - Fupp[k] > max_constraint_violation)
				{
					max_constraint_violation = F[k] - Fupp[k];
					worst_constraint = k;
				}
			}
			else if (F[k] < Flow[k])
			{
				if (Flow[k] - F[k] > max_constraint_violation)
				{
					max_constraint_violation = Flow[k] - F[k];
					worst_constraint = k;
				}
			}
		}

		return max_constraint_violation / EMTG::math::norm(this->Xtrial_scaled.data(), Problem->total_number_of_NLP_parameters);
	}

}} //close namespace