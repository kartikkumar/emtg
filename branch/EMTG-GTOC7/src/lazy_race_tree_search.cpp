//Donald Ellison and Ryne Beeson
//May 22nd 2014



#include "lazy_race_tree_search.h"
#ifdef EMTG_MPI
namespace boost {
	namespace mpi {
		template <> struct is_commutative<EMTG::pair_min, std::pair<int, double> > : mpl::true_{};
	}
}
#endif

namespace EMTG{
#ifdef EMTG_MPI
	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in, std::vector <int> & asteroid_list, std::vector <int> & best_sequence, std::vector <double> & epoch_sequence, std::vector <double> & mass_sequence, std::string & branch_directory, std::string & tree_summary_file_location, boost::mpi::environment & MPIenvironment, boost::mpi::communicator & world)
#else
	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in, std::vector <int> & asteroid_list, std::vector <int> & best_sequence, std::vector <double> & epoch_sequence, std::vector <double> & mass_sequence, std::string & branch_directory, std::string & tree_summary_file_location)
#endif
	{
		
		//make a copy of the options structure
		EMTG::missionoptions branch_options = *options;
		std::vector <int> & asteroid_sublist = asteroid_list;
		

		double time_left = options->lazy_race_tree_maximum_duration * 86400.0; //POSSIBLE OPTIONS STRUCTURE INCLUSION
		double time_to_remove = 0.0;
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


		//THESE SHOULD GO IN THE OPTIONS STRUCTURE/GUI
		//WITH A NOTE THAT THEY ONLY APPLY FOR MIN. PROP WITH THE TIME EXTENSION FEATURE
		
		double max_flight_time = branch_options.lazy_race_tree_final_flight_time_bound;
		double flight_time_increment = branch_options.lazy_race_tree_flight_time_increment;

		//THIS SHOULD BE FLIPPED TO TRUE IMMEDIATELY IF WE ARE DOING MIN. PROP WITHOUT THE TIME EXTENSION
		bool reached_max_upper_flighttime_bound = false;
		bool we_are_repeating_the_level = false;

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
		epoch_sequence.push_back(branch_options.launch_window_open_date);
		mass_sequence.push_back(branch_options.maximum_mass);

#ifdef EMTG_MPI //this block of code lets each core pick its own subset of asteroids to play with
		std::vector<int> my_asteroidlist;
		std::pair<int, double> my_best, absolute_best; //need these later, declared here because it is the first #ifdef available

		for (int index = 0; index < asteroid_list.size(); ++index) {
			if (index%(world.size()) == world.rank()) //this one belongs to me
				my_asteroidlist.push_back(asteroid_list[index]);
		}

#endif

		//We will most likely be receiving initial mothership arrival epochs, so we want to advance
		//the launch_window_open_date ahead by 30 days to force the requisit stay time
		branch_options.launch_window_open_date += 2592000.0;


		//Loop over levels in the tree
		do
		{

			if(!we_are_repeating_the_level)
				branch_options.total_flight_time_bounds[1] = options->total_flight_time_bounds[1];
			//EMTG::missionoptions branch_options = *options;

			//filter from the full list to the sublist
#ifdef EMTG_MPI			
			asteroid_sublist = filter_asteroid_list(starting_body_ID, branch_options.launch_window_open_date / 86400.0, my_asteroidlist, TheUniverse_in, options);
#else
			asteroid_sublist = filter_asteroid_list(starting_body_ID, branch_options.launch_window_open_date / 86400.0, asteroid_list, TheUniverse_in, options);
#endif
			number_of_branches_in_current_level = asteroid_sublist.size(); 


#ifndef EMTG_MPI //we only want to break if we're running solo; if we are MPI, we just want to recognize that we're "out"
			//we've run out of asteroids/hit a dead-end based on the current search 'ball'
			if (number_of_branches_in_current_level <= 0) {
				break;
			}
#endif
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

			
			failures_in_current_level = 0; //reset this for the new level

			//Loop over branches in the level
			for (size_t branch = 0; branch < number_of_branches_in_current_level && number_of_branches_in_current_level != 0; ++branch)
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
				
				if (options->enable_emtg_output_files)
				{
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
				}		
				

					

				//Instantiate a new mission, then optimize it
				EMTG::mission branch_mission(&branch_options, TheUniverse_in);
				branch_mission.optimize();


				current_flight_time = branch_mission.Xopt[1];
				current_mass = branch_mission.Xopt[2];
				current_wait_time = branch_mission.Xopt[0] - branch_options.launch_window_open_date;

				if (branch_options.objective_type == 1 && include_wait_time_in_cost)
					current_cost = branch_mission.best_cost + current_wait_time / TheUniverse_in[0].TU; //only one journey, therefore, always first entry in universe vector
				else
					current_cost = branch_mission.best_cost;

				

				//if this branch was not feasible, then make its cost really big
				if (branch_mission.number_of_solutions == 0)
				{
					//std::cout << "No feasible LRTS solution" << std::endl;
					//getchar();
					current_cost = 1.0e+20;
					++failures_in_current_level;
				}

				cost_to_get_to_each_body_in_level[branch] = current_cost;
				time_to_get_to_each_body_in_level[branch] = current_flight_time;
				final_mass_for_each_body_in_level[branch] = current_mass;
				wait_time_for_each_body_in_level[branch] = current_wait_time;
				
				write_branch_summary(branch_mission, branch_options, tree_summary_file_location, tree_level, branch);

			}//end branch loop
		
//#ifndef EMTG_MPI			
//			//if we fail to find any feasible solutions for every branch in the level
//			//then this is as far as this tree can go
//			if (failures_in_current_level == number_of_branches_in_current_level) {
//				return; //if we are not in MPI mode, then we've bottomed out this tree
//			}
//#endif

			//DETERMINE WHICH BRANCH WAS "CHEAPEST" 
			std::vector<double>::iterator best_cost_of_level = std::min_element(cost_to_get_to_each_body_in_level.begin(), cost_to_get_to_each_body_in_level.end());
			
#ifdef EMTG_MPI //now share our minimum across everyone and find the true minimum
			if (number_of_branches_in_current_level == 0) { //we didn't have any asteroids, so we skipped this round
				my_best = std::make_pair(world.rank(), 1.0e+20);
			} else {
				my_best = std::make_pair(world.rank(), *best_cost_of_level);
			}

#endif
			/*if this node did not produce ANY feasible results, we are going to pretend that it was a zero-sized node.  We do this AFTER we've registered its 
			correct objective function, etc.  This lines up with the later test that if *ALL* branches are zero-sized, we ran out of targets and we're done.  By
			pretending we are a zero-sized branch, if ALL branches PRETEND too, that really means they all failed out, and we've equally hit a dead-end. */
			int number_of_successes_in_current_level = number_of_branches_in_current_level - failures_in_current_level;

#ifdef EMTG_MPI
			//we make the number_of_branches_in_current_level be the TOTAL number across ALL satellites so the WHILE loop can exit if we sum to zero.  Also doubles as a global check of failure due to a previous if-statement
			number_of_successes_in_current_level = boost::mpi::all_reduce<int>(world, number_of_successes_in_current_level, std::plus<int>());
#endif

			
			if (number_of_successes_in_current_level == 0 )
			{
				//if we are running min. propellant and had no successes in this level
				//try the level again with a slightly extended max flight time
				//otherwise, if we are running another objective function we are actually done with this level
				//so we will just abandon the lazy race tree search
				if (branch_options.objective_type == 2)
				{
					if (reached_max_upper_flighttime_bound) {
						return;//need to exit out at this stage, for fear of otherwise having ambiguous reduce next and multi-point broadcast
					}
					//extend the flight time
					branch_options.total_flight_time_bounds[1] += flight_time_increment*86400.0;

					//we have reached the allowed upper limit of our flight time extension tactic
					//the next time throug the do-while will be the last
					if (branch_options.total_flight_time_bounds[1] >= max_flight_time*86400.0)
						reached_max_upper_flighttime_bound = true;

					we_are_repeating_the_level = true;
					//try the level again
					continue;
				}
				else
					return; 
			}
			
			//if we've gone this far, then we have found one feasible solution
			//we need to reset we_are_repeating_the_level if it has been turned on
			we_are_repeating_the_level = false;

#ifdef EMTG_MPI
			absolute_best = boost::mpi::all_reduce<std::pair<int, double> >(world, my_best, pair_min());
			
			//now we know who has the absolute best!
			if (absolute_best.first == world.rank()) { //I have the best.  Note this is a HANGING if statement across code
#endif

				//find index of body that had the best cost
				int next_starting_body_index = best_cost_of_level - cost_to_get_to_each_body_in_level.begin();
				time_to_remove = wait_time_for_each_body_in_level[next_starting_body_index] + time_to_get_to_each_body_in_level[next_starting_body_index];
				
				//ALTER WET MASS
				branch_options.maximum_mass = final_mass_for_each_body_in_level[next_starting_body_index];

				//ALTER STARTING EPOCH
				//time you waited after required 30 days + time it took to get to next body + 30 days required at the destination asteroid
				branch_options.launch_window_open_date += wait_time_for_each_body_in_level[next_starting_body_index] + time_to_get_to_each_body_in_level[next_starting_body_index] + 30.0*86400.0;

				//assign new starting body for the next level in the tree
				//this is the body with the best cost function in the current level
				starting_body_ID = asteroid_sublist[next_starting_body_index];
#ifdef EMTG_MPI
			}  //close the if statement, only when in MPI
				
			//share the time to remove
			boost::mpi::broadcast(world, time_to_remove, absolute_best.first);
			
			//broadcast updated final wetmass
			boost::mpi::broadcast(world, branch_options.maximum_mass, absolute_best.first);

			//broadcast updated launch window opening date
			boost::mpi::broadcast(world, branch_options.launch_window_open_date, absolute_best.first);

			//now need to broadcast out the best answer asteroid
			boost::mpi::broadcast(world, starting_body_ID, absolute_best.first);

			

			
			//the starting body that was just assigned should now be removed from the list of available targets
			std::vector<int>::iterator tobeerased = std::find(my_asteroidlist.begin(), my_asteroidlist.end(), starting_body_ID);
			if (tobeerased != my_asteroidlist.end())
				my_asteroidlist.erase(tobeerased);
#else
			asteroid_list.erase(std::find(asteroid_list.begin(), asteroid_list.end(), starting_body_ID));
#endif
			//append the best sequence
			best_sequence.push_back(starting_body_ID);
			//append the starting mass at each target
			mass_sequence.push_back(branch_options.maximum_mass);
			//append the epoch at each target
			epoch_sequence.push_back(branch_options.launch_window_open_date - 30.0*86400.0);
			
			//REDUCE TIME_LEFT
			time_left -= time_to_remove;

			
			++tree_level;

		}while (time_left > 0.0 && branch_options.maximum_mass > branch_options.minimum_dry_mass);//end level loop


		
	}//end lazy race tree search function


	//epoch needs to be in MJD

	std::vector<int> filter_asteroid_list(int const & current_asteroid, double const & epoch, std::vector<int> & asteroid_list, boost::ptr_vector<Astrodynamics::universe> & myUniverse, missionoptions * options) {

		std::vector<int> filtered_asteroid_list;
		double current_asteroid_state[6], target_asteroid_state[6];
		double euclidean_distance = 0, velocity_difference = 0;


		//reserve as much room as the other list, but don't actually allocate
		filtered_asteroid_list.reserve(asteroid_list.size());


		//determine where we are
		myUniverse[0].bodies[current_asteroid-1].locate_body(epoch, current_asteroid_state, false, options);

		for (std::vector<int>::iterator asteroid = asteroid_list.begin(); asteroid != asteroid_list.end(); ++asteroid) {

			myUniverse[0].bodies[(*asteroid)-1].locate_body(epoch, target_asteroid_state, false, options);

			euclidean_distance = std::sqrt((current_asteroid_state[0] - target_asteroid_state[0])*(current_asteroid_state[0] - target_asteroid_state[0]) + (current_asteroid_state[1] - target_asteroid_state[1])*(current_asteroid_state[1] - target_asteroid_state[1]) + (current_asteroid_state[2] - target_asteroid_state[2])*(current_asteroid_state[2] - target_asteroid_state[2]));
			velocity_difference = std::abs(std::sqrt((current_asteroid_state[3])*(current_asteroid_state[3]) + (current_asteroid_state[4])*(current_asteroid_state[4]) + (current_asteroid_state[5])*(current_asteroid_state[5])) - std::sqrt((target_asteroid_state[3])*(target_asteroid_state[3]) + (target_asteroid_state[4])*(target_asteroid_state[4]) + (target_asteroid_state[5])*(target_asteroid_state[5])));

			//should we include this point?
			if (euclidean_distance < options->lazy_race_tree_radius && velocity_difference < options->lazy_race_tree_velocity_difference) //update this to be pulling from the options file
				filtered_asteroid_list.push_back(*asteroid);

		}

		return filtered_asteroid_list;

	}


	void write_branch_summary(EMTG::mission & branch_mission, EMTG::missionoptions & branch_options, std::string & tree_summary_file_location, int & tree_level, const int & branch)
	{
		//open the summary file for appending
		std::ofstream outputfile(tree_summary_file_location.c_str(), std::ios::app);

		//Write things of interest for this branch
		outputfile.width(15); outputfile << left << tree_level;
		outputfile.width(15); outputfile << branch;
		outputfile.width(35); outputfile << branch_options.destination_list[0][0];
		outputfile.width(35); outputfile << branch_options.destination_list[0][1];
		outputfile.width(30); outputfile.precision(15); outputfile << branch_options.launch_window_open_date;
		outputfile.width(25); outputfile.precision(15); outputfile << branch_mission.Xopt[0];
		outputfile.width(25); outputfile.precision(15); outputfile << branch_mission.Xopt[1] / 86400.0; //DEBUGGING
		outputfile.width(25); outputfile.precision(15); outputfile << (branch_mission.Xopt[0] + branch_mission.Xopt[1] - branch_options.launch_window_open_date) / 86400.0; // DEBUGGING
		outputfile.width(25); outputfile.precision(15); outputfile << branch_mission.Xopt[0] + branch_mission.Xopt[1];
		outputfile.width(20); outputfile << branch_options.maximum_mass;
		outputfile.width(20); outputfile << branch_mission.Xopt[2];


		outputfile << std::endl;
		//close the tree summary file
		outputfile.close();
	}


	void load_asteroid_list(std::string & asteroid_filename, std::vector <int> & asteroid_list)
	{

		std::ifstream inputfile(asteroid_filename.c_str());

		int number;

		if (!inputfile.is_open())
		{
			std::cout << "Cannot find asteroid file for lazy race tree search: " + asteroid_filename << std::endl;
		}

		while (!inputfile.eof())
		{
			inputfile >> number;
			asteroid_list.push_back(number);
		}

	}


}//end EMTG namespace