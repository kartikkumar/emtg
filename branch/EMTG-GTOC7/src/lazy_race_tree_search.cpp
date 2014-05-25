//Donald Ellison and Ryne Beeson
//May 22nd 2014


#include "lazy_race_tree_search.h"

namespace EMTG{
	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in, std::vector <int> & asteroid_list, std::vector <int> & best_sequence, std::string & branch_directory, std::string & tree_summary_file_location)
	{
		
		//make a copy of the options structure
		EMTG::missionoptions branch_options = *options;
		std::vector <int> & asteroid_sublist = asteroid_list;
		

		double time_left = options->lazy_race_tree_maximum_duration * 86400.0; //POSSIBLE OPTIONS STRUCTURE INCLUSION
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

		int number_of_branches_in_current_level = num_asteroids - 1; //after we remove the starting body from consideration, this is how many banches are in the first level

		vector<int> journey_sequence;

		//POSSIBLE OPTIONS STRUCTURE INCLUSION
		//The first asteroid in the list is assumed to be the one you are starting from
		starting_body_ID = options->lazy_race_tree_start_location_ID;
		asteroid_list.erase(std::find(asteroid_list.begin(), asteroid_list.end(), starting_body_ID));
		best_sequence.push_back(starting_body_ID);


		//Loop over levels in the tree
		do
		{
			//filter from the full list to the sublist
			
			asteroid_sublist = filter_asteroid_list(starting_body_ID, branch_options.launch_window_open_date / 86400.0, asteroid_list, TheUniverse_in, options);
			
			number_of_branches_in_current_level = asteroid_sublist.size(); 
			
			//we've run out of asteroids/hit a dead-end based on the current search 'ball'
			if (number_of_branches_in_current_level <= 0) {
				break;
			}

			//resize them to the subfilter
			cost_to_get_to_each_body_in_level.resize(number_of_branches_in_current_level);
			time_to_get_to_each_body_in_level.resize(number_of_branches_in_current_level);
			final_mass_for_each_body_in_level.resize(number_of_branches_in_current_level);
			wait_time_for_each_body_in_level.resize(number_of_branches_in_current_level);

			//refill them, since sadly, resize even with a val argument doesn't guarantee old values are changed
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
				branch_options.destination_list[0][1] = asteroid_sublist[branch];
				
				
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
			starting_body_ID = asteroid_sublist[next_starting_body_index];
			best_sequence.push_back(starting_body_ID);

			//the starting body that was just assigned should now be removed from the list of available targets
			asteroid_list.erase(std::find(asteroid_list.begin(), asteroid_list.end(), starting_body_ID));

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


	//epoch needs to be in MJD
	std::vector<int> filter_asteroid_list(int const & current_asteroid, double const & epoch, std::vector<int> & asteroid_list, boost::ptr_vector<Astrodynamics::universe> & myUniverse, missionoptions * options) {

		std::vector<int> filtered_asteroid_list;
		double current_asteroid_state[6], target_asteroid_state[6];
		double euclidean_distance = 0, velocity_difference = 0;


		//reserve as much room as the other list, but don't actually allocate
		filtered_asteroid_list.reserve(asteroid_list.size());


		//determine where we are
		myUniverse[0].bodies[current_asteroid].locate_body(epoch, current_asteroid_state, false, options);

		for (std::vector<int>::iterator asteroid = asteroid_list.begin(); asteroid != asteroid_list.end(); ++asteroid) {

			myUniverse[0].bodies[*asteroid].locate_body(epoch, target_asteroid_state, false, options);

			euclidean_distance = std::sqrt((current_asteroid_state[0] - target_asteroid_state[0])*(current_asteroid_state[0] - target_asteroid_state[0]) + (current_asteroid_state[1] - target_asteroid_state[1])*(current_asteroid_state[1] - target_asteroid_state[1]) + (current_asteroid_state[2] - target_asteroid_state[2])*(current_asteroid_state[2] - target_asteroid_state[2]));
			velocity_difference = std::abs(std::sqrt((current_asteroid_state[3])*(current_asteroid_state[3]) + (current_asteroid_state[4])*(current_asteroid_state[4]) + (current_asteroid_state[5])*(current_asteroid_state[5])) - std::sqrt((target_asteroid_state[3])*(target_asteroid_state[3]) + (target_asteroid_state[4])*(target_asteroid_state[4]) + (target_asteroid_state[5])*(target_asteroid_state[5])));

			//should we include this point?
			if (euclidean_distance < 100000000 && velocity_difference < 2.0) //update this to be pulling from the options file
				filtered_asteroid_list.push_back(*asteroid);

		}

		return filtered_asteroid_list;

	}
}//end EMTG namespace