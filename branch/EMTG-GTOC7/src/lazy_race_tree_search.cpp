//Donald Ellison and Ryne Beeson
//May 22nd 2014


#include "lazy_race_tree_search.h"

namespace EMTG{
	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in, std::vector <int> & asteroid_list, std::vector <int> & best_sequence, std::string & branch_directory, std::string & tree_summary_file_location)
	{
		
		//make a copy of the options structure
		EMTG::missionoptions branch_options = *options;

		

		double time_left = 5.5 * 365.25 * 86400.0; //POSSIBLE OPTIONS STRUCTURE INCLUSION
		double current_cost;
		double current_flight_time;
		double current_mass;
		double current_wait_time;
		int num_asteroids = asteroid_list.size();
		int starting_body_ID;
		int failures_in_current_level = 0;
		int tree_level = 0;
		
		//This flag controls a hack objective function since the true "arrive as early as possible" objective function does not work well
		//This DOES NOT include the wait time in the SNOPT objective function, it just adds it at the end so that the 
		//tree search makes the decision to go with the asteroid that it gets to on the earliest calendar date
		bool include_wait_time_in_cost = true;

		std::vector <double> cost_to_get_to_each_body_in_level(num_asteroids, 0.0);
		std::vector <double> time_to_get_to_each_body_in_level(num_asteroids, 0.0);
		std::vector <double> final_mass_for_each_body_in_level(num_asteroids, 0.0);
		std::vector <double> wait_time_for_each_body_in_level(num_asteroids, 0.0);

		vector<int> journey_sequence;

		//POSSIBLE OPTIONS STRUCTURE INCLUSION
		//The first asteroid in the list is assumed to be the one you are starting from
		starting_body_ID = asteroid_list[0];
		asteroid_list.erase(asteroid_list.begin());
		best_sequence.push_back(starting_body_ID);

		int number_of_branches_in_current_level = num_asteroids - 1; //after we remove the starting body from consideration, this is how many banches are in the first level

		//Loop over levels in the tree
		do
		{
			std::fill(cost_to_get_to_each_body_in_level.begin(), cost_to_get_to_each_body_in_level.end(), 1.0e+20);
			std::fill(time_to_get_to_each_body_in_level.begin(), time_to_get_to_each_body_in_level.end(), 0.0);
			std::fill(cost_to_get_to_each_body_in_level.begin(), cost_to_get_to_each_body_in_level.end(), 1.0e+20);
			std::fill(wait_time_for_each_body_in_level.begin(), wait_time_for_each_body_in_level.end(), 1.0e+20);

			//Loop over branches in the level
			for (size_t branch = 0; branch < number_of_branches_in_current_level; ++branch)
			{
				//CHANGE DESTINATION LIST (WE ONLY EVER HAVE ONE JOURNEY, ONE PHASE)
				journey_sequence.clear();
				branch_options.sequence.clear();
				branch_options.phase_type.clear();

				//specify which two asteroids you're flying between
				branch_options.destination_list[0][0] = starting_body_ID;
				branch_options.destination_list[0][1] = asteroid_list[branch];
				
				
				journey_sequence.push_back(branch_options.destination_list[0][0]);
				journey_sequence.push_back(branch_options.destination_list[0][1]);
				branch_options.sequence.push_back(journey_sequence);
				branch_options.number_of_phases[0] = journey_sequence.size() - 1;

				vector<int> journey_phase_type(journey_sequence.size() - 1, branch_options.mission_type);
				branch_options.phase_type.push_back(journey_phase_type);

				ostringstream mission_name_stream;
				mission_name_stream << branch_options.destination_list[0][0] << "_" << branch_options.destination_list[0][1];
				branch_options.mission_name = mission_name_stream.str();

				

				//define a new working directory
				branch_options.working_directory = branch_directory + branch_options.mission_name;
				//create the working directory
				try
				{
					path p(branch_options.working_directory);
					boost::filesystem::create_directories(p);
				}
				catch (std::exception &e)
				{
					std::cerr << "Error " << e.what() << ": Directory creation failed" << std::endl;
				}


				//print the options file to the new directory
				branch_options.print_options_file(branch_options.working_directory + "//" + branch_options.mission_name + ".emtgopt");
				
				

				EMTG::mission branch_mission(&branch_options, TheUniverse_in);

			
				
				branch_mission.optimize();
				
				//if this branch was not feasible, then make its cost really big
				if (branch_mission.number_of_solutions == 0)
				{
					//std::cout << "No feasible LRTS solution" << std::endl;
					//getchar();
					current_cost = 1.0e+20;
					++failures_in_current_level;

					//if we fail to find any feasible solutions for every branch in the level
					//then this is as far as this tree can go
					if (failures_in_current_level == number_of_branches_in_current_level)
						return;
				}
				
				current_flight_time = branch_mission.Xopt[1];
				current_mass = branch_mission.Xopt[2];
				current_wait_time = branch_mission.Xopt[0] - branch_options.launch_window_open_date;

				if (branch_options.objective_type == 1 && include_wait_time_in_cost)
					current_cost = branch_mission.best_cost + current_wait_time/TheUniverse_in[0].TU; //only one journey, therefore, always first entry in universe vector
				else
					current_cost = branch_mission.best_cost;

				cost_to_get_to_each_body_in_level[branch] = current_cost;
				time_to_get_to_each_body_in_level[branch] = current_flight_time;
				final_mass_for_each_body_in_level[branch] = current_mass;
				wait_time_for_each_body_in_level[branch] = current_wait_time;
				
				write_branch_summary(branch_mission, branch_options, tree_summary_file_location, tree_level, branch);

			}//end branch loop

			//DETERMINE WHICH BRANCH WAS "CHEAPEST" 
			std::vector<double>::iterator best_cost_of_level = std::min_element(cost_to_get_to_each_body_in_level.begin(), cost_to_get_to_each_body_in_level.end());
			
			//std::cout << *best_cost_of_level << std::endl;

			//find index of body that had the best cost
			int next_starting_body_index = best_cost_of_level - cost_to_get_to_each_body_in_level.begin();

			//assign new starting body for the next level in the tree
			//this is the body with the best cost function in the current level
			starting_body_ID = asteroid_list[next_starting_body_index];
			best_sequence.push_back(starting_body_ID);

			//the starting body that was just assigned should now be removed from the list of available targets
			asteroid_list.erase(asteroid_list.begin() + next_starting_body_index);

			//REDUCE TIME_LEFT
			//we are letting the probe go for 5.5 years after which it should probably stop
			time_left -= wait_time_for_each_body_in_level[next_starting_body_index] + time_to_get_to_each_body_in_level[next_starting_body_index];
			
			//reduce number of bodies for the next branch
			--number_of_branches_in_current_level;

			//ALTER WET MASS
			branch_options.maximum_mass = final_mass_for_each_body_in_level[next_starting_body_index];

			//ALTER STARTING EPOCH
			//time you waited after required 30 days + time it took to get to next body + 30 days required at the destination asteroid
			branch_options.launch_window_open_date += wait_time_for_each_body_in_level[next_starting_body_index] + time_to_get_to_each_body_in_level[next_starting_body_index] + 30.0*86400.0;

			++tree_level;

		}while (time_left > 0.0 && number_of_branches_in_current_level > 0);//end level loop


		
	}

}//end EMTG namespace