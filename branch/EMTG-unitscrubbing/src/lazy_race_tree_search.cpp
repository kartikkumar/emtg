#include "lazy_race_tree_search.h"

namespace EMTG{
	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in)
	{
		//call lazy race-tree search

		//I suggest putting this in a separate source file but here is the general outline:

		//1. pick a starting body (maybe this is an input? maybe you arrange it so you always start from the first body in the list?)

		//2. farm out and evaluate earliest arrival date trajectories with constrained dry mass to every body in the list except the one that you started from

		//3. pick the cheapest

		//4. repeat from (2), checking to make sure there are no duplicates. Only check for duplicates if the options->lazy_race_tree_allow_duplicates is set to 0.
		//(this is so that we can use lazy race-tree for a GTOC6 type problem that allows duplicates at some later date)

		//if you need to add stuff in the options file search for lazy race-tree in missionoptions->cpp

		//set up a batch output file

		
		//SPECIFY STARTING ASTEROID ID

		double time_left = 6.0 * 365.25 * 86400.0;
		double current_cost;
		int sequence_length = 3;
		int num_bodies_left_to_evaluate = 3;
		int starting_body_ID;

		std::vector <int> target_sequence(sequence_length, 0);
		std::vector <double> time_to_get_to_each_body_in_branch(sequence_length, 0.0);

		//Populate the sequence vector ("number" column from the Universe file
		starting_body_ID = 54;
		target_sequence[0] = 62;
		target_sequence[1] = 13;
		target_sequence[2] = 65;
		//target_sequence[3] = 3;
		//target_sequence[4] = 101;

		//Loop over levels in the tree
		do
		{
			std::fill(time_to_get_to_each_body_in_branch.begin(), time_to_get_to_each_body_in_branch.end(), 0.0);

			//Loop over branches in the level
			for (size_t body = 0; body < num_bodies_left_to_evaluate; ++body)
			{
				//CHANGE DESTINATION LIST (WE ONLY EVER HAVE ONE JOURNEY, ONE PHASE)
				missionoptions branch_options = *options;

				branch_options.destination_list[0][0] = target_sequence[0];
				branch_options.destination_list[0][1] = target_sequence[1];
				
				vector<int> journey_sequence;
				journey_sequence.push_back(branch_options.destination_list[0][0]);
				journey_sequence.push_back(branch_options.destination_list[0][1]);
				branch_options.sequence.push_back(journey_sequence);
				branch_options.number_of_phases[0] = journey_sequence.size() - 1;
				vector<int> journey_phase_type(journey_sequence.size() - 1, branch_options.mission_type);
				branch_options.phase_type.push_back(journey_phase_type);

				ostringstream mission_name_stream;
				mission_name_stream << options->destination_list[0][0] << "_" << options->destination_list[0][1];
				branch_options.mission_name = mission_name_stream.str();

				//current time
				ptime now = second_clock::local_time();
				std::stringstream timestream;
				timestream << static_cast<int>(now.date().month()) << now.date().day() << now.date().year() << "_" << now.time_of_day().hours() << now.time_of_day().minutes() << now.time_of_day().seconds();


				//define a new working directory
				branch_options.working_directory = "..//EMTG_v8_results//" + options->mission_name + "_" + timestream.str();
				//create the working directory
				try
				{
					path p(branch_options.working_directory);
					path puniverse(branch_options.working_directory + "/Universe");
					boost::filesystem::create_directories(p);
					boost::filesystem::create_directories(puniverse);
				}
				catch (std::exception &e)
				{
					std::cerr << "Error " << e.what() << ": Directory creation failed" << std::endl;
				}


				//print the options file to the new directory
				branch_options.print_options_file(branch_options.working_directory + "//" + options->mission_name + ".emtgopt");
				string outputfilestring = branch_options.working_directory + "//" + options->mission_name + "_batch_summary.emtgbatch";
				std::ofstream outputfile(outputfilestring.c_str(), ios::trunc);
				outputfile.width(30); outputfile << left << "Sequence";
				outputfile.width(3); outputfile << " | ";

				outputfile << endl;
				for (int k = 0; k < 50; ++k)
					outputfile << "-";
				outputfile << endl;
				outputfile.close();

				//TRY CATCH - IF MISSION IS INFEASIBLE JUST MOVE TO THE NEXT MISSION
				
				EMTG::mission branch_mission(&branch_options, TheUniverse_in);

				try{
					branch_mission.optimize();
				}
				catch (int err)
				{
					getchar();
				}
				current_cost = branch_mission.best_cost;
				
				//PLACE CURRENT FLIGHT TIME INTO THE TIME TO GET TO EACH BODY VECTOR
				

			}//end branch loop

			//DETERMINE WHICH BODY WAS "CHEAPEST" TO GET TO (MINIMUM OF TIME VECTOR)
			//REMOVE THAT BODY FROM THE SEQUENCE LIST

			//REDUCE TIME_LEFT
			//REMOVE STARTING BODY FROM BODY LIST
			//MAKE NEW BODY LIST THE CURRENT BODY LIST
			--num_bodies_left_to_evaluate;

		}while (time_left > 0.0 && num_bodies_left_to_evaluate > 0);//end level loop


	}

}//end EMTG namespace