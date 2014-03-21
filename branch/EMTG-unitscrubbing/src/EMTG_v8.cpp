//============================================================================
// Name        : EMTG_v8.cpp
// Author      : Jacob Englander
// Version     :
// Copyright   : 
// Description : Main launch function for EMTG_v8
// Description : EMTG_v8 is a generic optimizer that hands MGA, MGADSM, MGALT, and FBLT mission types
//============================================================================

#include "missionoptions.h"
#include "mission.h"
#include "integerGA.h"
#include "outerloop_NSGAII.h"


#include "universe.h"
#include "body.h"

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"

#include "SpiceUsr.h"

#include <iostream>
#include <fstream>
#include <sstream>

#ifdef EMTG_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#endif

using namespace std;
using namespace boost::filesystem;
using namespace boost::gregorian;
using namespace boost::posix_time;

int main(int argc, char* argv[]) 
{
	//delete the fort if present
#ifndef _STONEAGECplusplus
	fs::path fort(L"fort.1"); 
	fs::remove(fort);
#endif


	cout << "program starting" << endl;

	//parse the options file
	string options_file_name;
	if (argc == 1)
		options_file_name = "default.emtgopt";
	else if (argc == 2) 
		options_file_name.assign(argv[1]);

	cout << options_file_name << endl;

	EMTG::missionoptions options(options_file_name);
	cout << options.error_message << endl;
	if (!(options.file_status == 0))
	{
		cout << "Aborting program run." << endl;
		cin.ignore();
		return 0;
	}

	//create a working directory for the problem
	//TODO this may change when we develop a parallel version
	//*****************************************************************
	ptime now = second_clock::local_time();
	std::stringstream timestream;
	timestream << static_cast<int>(now.date().month()) << now.date().day() << now.date().year() << "_" << now.time_of_day().hours() << now.time_of_day().minutes() << now.time_of_day().seconds();

		
	//define a new working directory
	options.working_directory = "..//EMTG_v8_results//" + options.mission_name + "_" + timestream.str();
	//create the working directory
	try 
	{ 
		path p(options.working_directory);
		path puniverse(options.working_directory + "/Universe");
        boost::filesystem::create_directories (p); 
		boost::filesystem::create_directories (puniverse);
    } 
	catch (std::exception &e) 
	{ 
		std::cerr << "Error " << e.what() << ": Directory creation failed" << std::endl; 
    } 


	//print the options file to the new directory
	options.print_options_file(options.working_directory + "//" + options.mission_name + ".emtgopt");
	

	//load all ephemeris data
	vector<fs::path> SPICE_files;
	EMTG::filesystem::get_all_files_with_extension(fs::path(options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);
	EMTG::filesystem::get_all_files_with_extension(fs::path(options.universe_folder + "/ephemeris_files/"), ".cmt", SPICE_files);

	string filestring;
	for (size_t k = 0; k < SPICE_files.size(); ++k)
	{
		filestring = options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();
		furnsh_c(filestring.c_str());
	}

	//SPICE reference frame kernel
	string leapsecondstring = options.universe_folder + "/ephemeris_files/" + options.SPICE_leap_seconds_kernel;
	string referenceframestring = options.universe_folder + "/ephemeris_files/" + options.SPICE_reference_frame_kernel;
	furnsh_c(leapsecondstring.c_str());
	furnsh_c(referenceframestring.c_str());

	//disable SPICE errors. This is because we can, and will often, go off the edge of an ephemeris file.
	errprt_c ("SET", 100, "NONE");
	erract_c ("SET", 100, "RETURN");

	//create a vector of universes for each journey
	boost::ptr_vector<EMTG::Astrodynamics::universe> TheUniverse;
	options.TU = 0;
	for (int j = 0; j < options.number_of_journeys; ++j)
	{
		TheUniverse.push_back(new EMTG::Astrodynamics::universe(j, options.universe_folder + "//" + options.journey_central_body[j] + ".emtg_universe", &options));
		stringstream universenamestream;

		universenamestream << options.journey_central_body[j] + "_Journey_" << j << ".universe_output";

		TheUniverse[j].print_universe(options.working_directory + "//" +"universe//" + universenamestream.str(), &options);

		cout << options.working_directory + "//" +"universe//" + universenamestream.str() << " written" << endl;

		if (TheUniverse[j].TU > options.TU)
			options.TU = TheUniverse[j].TU;
	}

	//*****************************************************************


	//next, it is time to start the outer-loop
	//if we are running in parallel, start MPI
#ifdef EMTG_MPI
	boost::mpi::environment MPIEnvironment;
	boost::mpi::communicator MPIWorld;
#endif

	if (options.run_outerloop && options.outerloop_objective_function_choices.size() == 0)
	{
		// TODO call to outer-loop GA

		//Step 1: instantiate a random number generator
		boost::mt19937 RNG;
		RNG.seed(std::time(0));

		//Step 2: instantiate a GA object
		integerGA outerloop(&options, TheUniverse, RNG);

		//Step 3: make sure that output files are defined
		if (options.solutions_file == "")
			options.solutions_file = options.working_directory + "//solutions.txt";
		if (options.population_file == "")
			options.population_file = options.working_directory + "//population.txt";
		if (options.convergence_file == "")
			options.convergence_file = options.working_directory + "//convergence.txt";

		//Step 4: run the GA
		outerloop.evolve(&options, TheUniverse);
	}
	else if (options.run_outerloop && options.outerloop_objective_function_choices.size() > 0)
	{
		//Step 1: instantiate an NSGA-II object
#ifdef EMTG_MPI
		GeneticAlgorithm::outerloop_NSGAII NSGAII(options, &MPIEnvironment, &MPIWorld);
#else
		GeneticAlgorithm::outerloop_NSGAII NSGAII(options);
#endif

		//Step 2: generate a random population
		NSGAII.set_populationsize(options.outerloop_popsize);
		NSGAII.set_mutationrate(options.outerloop_mu);
		NSGAII.set_max_generations(options.outerloop_genmax);
		NSGAII.generatepop();
		NSGAII.startclock();
		NSGAII.evolve(options, TheUniverse);

		//Step 3: evaluate the population, in this case meaning just get the names of the files
		//NSGAII.evaluatepop(options, TheUniverse);

		//Step 4: write out the population
		NSGAII.writepop(options.working_directory + "//NSGAII_final_population.csv");

		//Step 5: write out the archive
		NSGAII.write_archive(options.working_directory + "//NSGAII_archive.csv");
	}
	else
	{
		//set up a batch output file
		string outputfilestring = options.working_directory + "//" + options.mission_name + "_batch_summary.emtgbatch";
		std::ofstream outputfile(outputfilestring.c_str(), ios::trunc);
		outputfile.width(30); outputfile << left << "Sequence";
		outputfile.width(3); outputfile << " | ";

		switch (options.objective_type)
		{
		case 0: //deltaV
			outputfile.width(17); outputfile << left << "deltaV (km/s)";
			break;
		case 2:
			outputfile.width(17); outputfile << left << "delivered mass (kg)";
			break;
		}
		
		outputfile << endl;
		for (int k = 0; k < 50; ++k)
		outputfile << "-";
		outputfile << endl;
		outputfile.close();

		//run the trial outer-loop vectors

		for (int trial = 0; trial < options.number_of_trial_sequences; ++trial)
		{
			//first we need to create an outer-loop decision vector based on the parameters from the options file
			vector<int> Xouterloop_trial;

			//if we have specified the mission type (all MGA, MGA-DSM, MGA-LT, FBLT, FBLT-S, MGA-NDSM) then it only takes one decision variable to encode a phase
			//if we have NOT specified the mission type, it takes two decision variables to encode a phase
			int phase_encode_length = (options.mission_type > 6 ? 2 : 1);
			
			for (int j = 0; j < options.number_of_journeys; ++j)
			{
				for (int p = 0; p < options.max_phases_per_journey; ++p)
				{
					Xouterloop_trial.push_back(options.sequence_input[trial][j][p]);

					//encode the phase type in the decision vector only if we are allowing the outer-loop to choose phase type
					if (options.sequence_input[trial][j][p] > 0 && phase_encode_length == 2)
							Xouterloop_trial.push_back(options.phase_type_input[j][p]);
				}

				if (phase_encode_length == 2)
					Xouterloop_trial.push_back(options.phase_type_input[j][options.phase_type_input[j].size() - 1]);

				options.number_of_phases = options.number_of_phases_input[trial];
			}

			//next, instantiate and optimize a problem object
			try
			{
				if (options.problem_type == 0) //regular EMTG missions
				{
					EMTG::mission TrialMission(&Xouterloop_trial[0], &options, TheUniverse, 0, 0);

					//and now, as a demo, print the mission tree
					TrialMission.output_mission_tree(options.working_directory + "//" + TrialMission.options.mission_name + "_" + TrialMission.options.description + "_missiontree.emtgtree");

					//copy the appropriate trialX, if necessary
					if (options.run_inner_loop == 0 || options.run_inner_loop == 4)
					{
						TrialMission.options.current_trialX = options.trialX[trial];
						
						//if we are optimizing an FBLT-S mission with an MGA-LT or FBLT initial guess
						if (options.mission_type == 4 && (TrialMission.options.current_trialX.size() == TrialMission.Xdescriptions.size() - TrialMission.options.total_number_of_phases))
						{
							EMTG::missionoptions temp_options = options;
							temp_options.mission_type = 3;

							EMTG::mission temp_mission(Xouterloop_trial.data(), &temp_options, TheUniverse, 0, 0);

							temp_mission.evaluate(TrialMission.options.current_trialX.data(), temp_mission.F.data(), temp_mission.G.data(), 0, temp_mission.iGfun, temp_mission.jGvar);

							TrialMission.options.current_trialX = temp_mission.create_initial_guess(TrialMission.options.current_trialX, TrialMission.Xdescriptions);
						}

						//if we are interpolating an initial guess to change the resolution
						if (options.interpolate_initial_guess && options.run_inner_loop > 0)
						{
							TrialMission.interpolate(Xouterloop_trial.data(), options.current_trialX);
						}
						
					}
					else if (options.run_inner_loop == 2)
					{
						TrialMission.options.current_trialX = options.trialX[trial];

						if (options.interpolate_initial_guess && options.seed_MBH)
							TrialMission.interpolate(Xouterloop_trial.data(), options.current_trialX);

					}

					

					//evaluate the mission
					TrialMission.optimize();

					//output the mission
					TrialMission.output();

					//output GMAT files
					//temporarily we can't output MGA or MGA-DSM missions
					if (options.mission_type > 1 && options.create_GMAT_script)
					{
						TrialMission.output_GMAT_preamble();
						TrialMission.output_GMAT_mission();
					}

					//store the results in a database file
					string outputfilestring = options.working_directory + "//" + options.mission_name + "_batch_summary.emtgbatch";
					std::ofstream outputfile(outputfilestring.c_str(), ios::app);
					outputfile.width(30); outputfile << left << TrialMission.options.description;
					outputfile.width(3); outputfile << " | ";

					switch (options.objective_type)
					{
					case 0: //deltaV
						outputfile.width(17); outputfile << left << TrialMission.F[0];
						break;
					case 2:
						outputfile.width(17); outputfile << left << TrialMission.F[0] * -options.maximum_mass;
						break;
					}
					outputfile << endl;
					outputfile.close();
				}
			}
			catch (std::exception &e)
			{
				cout << "Error " << e.what() << ": Failure to run inner-loop solver" << endl;
			}
		}//end loop over trial sequences


	}

	//unload SPICE


	for (size_t k = 0; k < SPICE_files.size(); ++k)
	{
		filestring = options.universe_folder + "ephemeris_files/" + SPICE_files[k].string();
		unload_c(filestring.c_str());
	}

	unload_c((options.universe_folder + "ephemeris_files/" + options.SPICE_leap_seconds_kernel).c_str());
	unload_c((options.universe_folder + "ephemeris_files/" + options.SPICE_reference_frame_kernel).c_str());	

#ifndef BACKGROUND_MODE
	std::cout << "EMTG run complete. Press enter to close window." << endl;
	std::cin.ignore();
#endif

	return 0;
}
