#include "lazy_race_tree_search.h"

namespace EMTG{
	void LRTS(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in)
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

		//Instantiate a mission object
		EMTG::mission branch_mission(options, TheUniverse_in);

		double time_left = 6.0 * 31557600.0;
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
			for (size_t j = 0; j < num_bodies_left_to_evaluate; ++j)
			{
				//CHANGE DESTINATION LIST (WE ONLY EVER HAVE ONE JOURNEY, ONE PHASE)
				options->destination_list[0][0] = target_sequence[0];
				options->destination_list[0][1] = target_sequence[1];

				ostringstream mission_name_stream;
				mission_name_stream << options->destination_list[0][0] << "_" << options->destination_list[0][1];
				options->mission_name = mission_name_stream.str();

				//current time
				ptime now = second_clock::local_time();
				std::stringstream timestream;
				timestream << static_cast<int>(now.date().month()) << now.date().day() << now.date().year() << "_" << now.time_of_day().hours() << now.time_of_day().minutes() << now.time_of_day().seconds();


				//define a new working directory
				options->working_directory = "..//EMTG_v8_results//" + options->mission_name + "_" + timestream.str();
				//create the working directory
				try
				{
					path p(options->working_directory);
					path puniverse(options->working_directory + "/Universe");
					boost::filesystem::create_directories(p);
					boost::filesystem::create_directories(puniverse);
				}
				catch (std::exception &e)
				{
					std::cerr << "Error " << e.what() << ": Directory creation failed" << std::endl;
				}


				//print the options file to the new directory
				options->print_options_file(options->working_directory + "//" + options->mission_name + ".emtgopt");
				string outputfilestring = options->working_directory + "//" + options->mission_name + "_batch_summary.emtgbatch";
				std::ofstream outputfile(outputfilestring.c_str(), ios::trunc);
				outputfile.width(30); outputfile << left << "Sequence";
				outputfile.width(3); outputfile << " | ";

				outputfile << endl;
				for (int k = 0; k < 50; ++k)
					outputfile << "-";
				outputfile << endl;
				outputfile.close();

				//TRY CATCH - IF MISSION IS INFEASIBLE JUST MOVE TO THE NEXT MISSION
				double cost;
				try{
					branch_mission.optimize();
				}
				catch (int err)
				{
					getchar();
				}
				cost = branch_mission.F[0];

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