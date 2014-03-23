//outer-loop NSGA-II
//written to solve systems optimization as formulated by J. Englander, M. Vavrina, and D. Ellison
//collaborative effort by J. Englander and M. Vavrina based on A. Ghosh's abstract GA spec and M. Vavrina's NSGA-II spec
#include "outerloop_NSGAII.h"
#include "mission.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "boost/algorithm/string.hpp"

namespace GeneticAlgorithm
{
		EMTG_outerloop_solution::EMTG_outerloop_solution(const std::vector<int>& X_in, const int& number_of_objectives)
	{
		this->X = X_in;
		this->fitness_values.resize(number_of_objectives, 1.0e+100);
	}

	bool EMTG_outerloop_solution::operator>(EMTG_outerloop_solution& OtherSolution)
	{
		return this->compare_objective_greaterthan(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::operator<(EMTG_outerloop_solution& OtherSolution)
	{
		return this->compare_objective_lessthan(OtherSolution, 0, 1.0e-6);
	}
	
	bool EMTG_outerloop_solution::operator>=(EMTG_outerloop_solution& OtherSolution)
	{
		return this->compare_objective_greaterthanorequalto(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::operator<=(EMTG_outerloop_solution& OtherSolution)
	{
		return this->compare_objective_lessthanorequalto(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::operator==(EMTG_outerloop_solution& OtherSolution)
	{
		return this->compare_objective_equalto(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::compare_objective_greaterthan(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] > tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_greaterthanorequalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] >= tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_lessthan(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] < tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_lessthanorequalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] <= tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_equalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (fabs(this->fitness_values[index] - OtherSolution.fitness_values[index]) < tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_description(EMTG_outerloop_solution& OtherSolution)
	{
		return (this->description == OtherSolution.description ? true : false);
	}

	//function to evaluate a member of the population
	void EMTG_outerloop_solution::set_data_pointers(void* optionspointer, void* Universepointer)
	{
		options_base = (EMTG::missionoptions*) optionspointer;
		Universe = (boost::ptr_vector< EMTG::Astrodynamics::universe>*) Universepointer;
	}

	EMTG::missionoptions EMTG_outerloop_solution::parse_outer_loop_decision_vector()
	{
		//Step 1: initialize a counter and a description string
		int Xindex = 0;
		stringstream descriptionstream;
		descriptionstream << options_base->mission_name;

		//Step 2: create a new options structure as a copy of the original
		EMTG::missionoptions options = *options_base;
				
		//Step 3: decode decision vector and create an options structure
		//we must decode the decision vector in the order that we encoded it
		
		//global mission choices
		if (options.outerloop_vary_power && options.outerloop_power_choices.size() > 1)
		{
			options.power_at_1_AU = options.outerloop_power_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_" << options.power_at_1_AU << "kW";
		}
		if (options.outerloop_vary_launch_epoch && options.outerloop_launch_epoch_choices.size() > 1)
		{
			options.launch_window_open_date = options.outerloop_launch_epoch_choices[X[Xindex]] * 86400.0;
			++Xindex;
			descriptionstream << "_LD" << (int) (options.launch_window_open_date / 86400.0);
		}
		if (options.outerloop_vary_flight_time_upper_bound && options.outerloop_flight_time_upper_bound_choices.size() > 1)
		{
			options.total_flight_time_bounds[1] = options.outerloop_flight_time_upper_bound_choices[X[Xindex]] * 86400;

			if (X[Xindex] > 0 && options.outerloop_restrict_flight_time_lower_bound)
				options.total_flight_time_bounds[0] = options.outerloop_flight_time_upper_bound_choices[X[Xindex] - 1] * 86400.0;
			else
				options.total_flight_time_bounds[0] = 0.0;

			++Xindex;
			descriptionstream << "_" << (int)(options.total_flight_time_bounds[1] / 86400.0) << "d";
		}
		if (options.outerloop_vary_thruster_type && options.outerloop_thruster_type_choices.size() > 1)
		{
			options.engine_type = options.outerloop_thruster_type_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_" << options.thruster_names[options.engine_type];
		}
		if (options.outerloop_vary_number_of_thrusters && options.outerloop_number_of_thrusters_choices.size() > 1)
		{
			options.number_of_engines = options.outerloop_number_of_thrusters_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_nTh" << options.number_of_engines;
		}
		if (options.outerloop_vary_launch_vehicle && options.outerloop_launch_vehicle_choices.size() > 1)
		{
			options.LV_type = options.outerloop_launch_vehicle_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_" << options.LV_names[options.LV_type];
		}
		if (options.outerloop_vary_departure_C3 && options.outerloop_departure_C3_choices.size() > 1)
		{
			if (X[Xindex] > 0)
				options.journey_initial_impulse_bounds[0][0] = sqrt(options.outerloop_departure_C3_choices[X[Xindex] - 1]);
			else
				options.journey_initial_impulse_bounds[0][0] = 0.0;
			options.journey_initial_impulse_bounds[0][1] = sqrt(options.outerloop_departure_C3_choices[X[Xindex]]);
			
			descriptionstream << "_C3d" << options.outerloop_departure_C3_choices[X[Xindex]];
			++Xindex;
		}
		if (options.outerloop_vary_arrival_C3 && options.outerloop_arrival_C3_choices.size() > 1)
		{
			if (X[Xindex] > 0)
				options.journey_final_velocity.back()[0] = sqrt(options.outerloop_arrival_C3_choices[X[Xindex] - 1]);
			else
				options.journey_final_velocity.back()[0] = 0.0;
			options.journey_final_velocity.back()[1] = sqrt(options.outerloop_arrival_C3_choices[X[Xindex]]);
			
			descriptionstream << "_C3a" << options.outerloop_arrival_C3_choices[X[Xindex]];
			++Xindex;
		}
		//choices for each journey
		for (size_t j = 0; j < options.number_of_journeys; ++j)
		{
			vector<int> journey_sequence;
			descriptionstream << "_" << (*Universe)[j].central_body_name << "(";

			//journey destination
			if (options.outerloop_vary_journey_destination[j] && options.outerloop_journey_destination_choices[j].size() > 1)
			{
				//for successive journeys, starting location is the previous journey's destination
				if (j > 0)
				{
					options.destination_list[j][0] = options.destination_list[j-1][1];
				}

				options.destination_list[j][1] = options.outerloop_journey_destination_choices[j][X[Xindex]];
				++Xindex;
			}

			//encode journey destination 1
			journey_sequence.push_back(options.destination_list[j][0]);
			switch (options.destination_list[j][0])
			{
				case -3: //begin at point on orbit
					{
						descriptionstream << "o";
						break;
					}
				case -2: //begin at fixed point
					{
						descriptionstream << "p";
						break;
					}
				case -1: //begin at SOI
					{
						descriptionstream << "s";
						break;
					}
				default:
					descriptionstream << (*Universe)[j].bodies[options.destination_list[j][0]-1].short_name;
			}

			//flyby sequence using Englander and Conway "null gene" method
			if (options.outerloop_vary_journey_flyby_sequence[j] && options.outerloop_journey_flyby_sequence_choices[j].size() > 1)
			{
				for (size_t p = 0; p < options.outerloop_journey_maximum_number_of_flybys[j]; ++p)
				{
					//is this a valid flyby or has the outer-loop chosen a null code?
					if (X[Xindex] < options.outerloop_journey_flyby_sequence_choices[j].size())
					{
						journey_sequence.push_back(options.outerloop_journey_flyby_sequence_choices[j][X[Xindex]]);
						descriptionstream << (*Universe)[j].bodies[journey_sequence.back()-1].short_name;
					}

					++Xindex;
					
				}
			}

			//encode journey destination 2
			journey_sequence.push_back(options.destination_list[j][1]);
			switch (options.destination_list[j][1])
			{
				case -3: //begin at point on orbit
					{
						descriptionstream << "o";
						break;
					}
				case -2: //begin at fixed point
					{
						descriptionstream << "p";
						break;
					}
				case -1: //begin at SOI
					{
						descriptionstream << "s";
						break;
					}
				default:
					descriptionstream << (*Universe)[j].bodies[options.destination_list[j][1]-1].short_name;
			}
			
			descriptionstream << ")";

			options.sequence.push_back(journey_sequence);
			options.number_of_phases[j] = journey_sequence.size() - 1;
			vector<int> journey_phase_type(journey_sequence.size() - 1, options.mission_type);
			options.phase_type.push_back(journey_phase_type);
		}

		vector< vector<int> > sequence_input;
		vector< vector<int> > phase_type_input;

		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			vector<int> journey_sequence_input;
			vector<int> journey_phase_type_input;
			for (int p = 1; p < options.max_phases_per_journey + 1; ++p)
			{
				if (p < options.number_of_phases[j])
				{
					journey_sequence_input.push_back(options.sequence[j][p]);
					journey_phase_type_input.push_back(options.phase_type[j][p]);
				}
				else
				{
					journey_sequence_input.push_back(0);
					journey_phase_type_input.push_back(options.mission_type);
				}
			}
			sequence_input.push_back(journey_sequence_input);
			phase_type_input.push_back(journey_phase_type_input);
		}
		options.sequence_input.clear();
		options.phase_type_input.clear();
		options.sequence_input.push_back(sequence_input);
		options.number_of_trial_sequences = 1;
		options.phase_type_input = phase_type_input;

		options.description = descriptionstream.str();
		options.mission_name = options.description;
		this->description = options.description;

		//if this is a previously evaluated solution that we want to take another cut at, seed the optimizer from its decision vector
		if (this->Xinner.size() > 0)
		{
			options.seed_MBH = 1;
			options.current_trialX = this->Xinner;
		}
		else
			options.seed_MBH = 0;

		return options;
	}

	void EMTG_outerloop_solution::evaluate()
	{
		//Step 1: create the options structure
		EMTG::missionoptions options = this->parse_outer_loop_decision_vector();

		//Step 2: instantiate a mission
		EMTG::mission TrialMission(&options, *Universe);

		//Step 3: run the inner-loop
		options.print_options_file(options.working_directory + "//" + this->description + ".emtgopt");
		if (!options.quiet_outerloop)
		{
			TrialMission.output_mission_tree(options.working_directory + "//" + this->description + ".emtgtree");
			cout << "Optimizing mission: " << TrialMission.options.description << endl;
		}

		//Step 3.1: Do a quick check - if any journey destination is visited twice (later make this a switch) then assign poor fitness values and return
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			for (int jj = 0; jj < options.number_of_journeys; ++jj)
			{
				if ( !(jj == j) && (options.destination_list[j][1] == options.destination_list[jj][1]) )
				{
					std::cout << "Mission " <<TrialMission.options.description << " contains duplicate journeys. Discarding." << std::endl;
					this->innerloop_fitness = 1.0e+100;
					for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
						this->fitness_values[objective] = 1.0e+100;
					return;
				}
			}
		}
		//Step 3.2 if we made it this far then optimize
		TrialMission.optimize();

		//Step 4: extract the fitness values
		//if the trial mission has no feasible solution then set the fitness values very high
		if (TrialMission.number_of_solutions == 0)
		{
			this->innerloop_fitness = 1.0e+100;
			for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
				this->fitness_values[objective] = 1.0e+100;
			return;
		}
		else
		{
			this->Xinner = TrialMission.Xopt;
			TrialMission.extract_objective_function_values(this->fitness_values);
			this->innerloop_fitness = TrialMission.F[0];
		}

		return;
	}

	void EMTG_outerloop_solution::write_csv_line(std::string filename)
	{
		std::ofstream outputfile(filename.c_str(), ios::out | ios::app);

		for (size_t gene = 0; gene < this->X.size(); ++gene)
			outputfile << this->X[gene] << ",";
		outputfile << this->description << "," << this->generation_found << "," << this->timestamp;
		for (size_t objective = 0; objective < this->fitness_values.size(); ++objective)
			outputfile << "," << this->fitness_values[objective];
		for (size_t Xinner_entry = 0; Xinner_entry < this->Xinner.size(); ++Xinner_entry)
			outputfile << "," << Xinner[Xinner_entry];
		outputfile << std::endl;

		outputfile.close();
	}

	//*****************************************methods for NSGA-II
	//constructors
	outerloop_NSGAII::outerloop_NSGAII(const EMTG::missionoptions& options)
	{
		//seed the random number generator
		this->RNG.seed(time(NULL));

		//calculate upper and lower bounds on the decision vector
		this->calcbounds(options);
	}
	outerloop_NSGAII::outerloop_NSGAII(int boost_random_seed, const EMTG::missionoptions& options)
	{
		//seed the random number generator
		this->RNG.seed(boost_random_seed);
		
		//calculate upper and lower bounds on the decision vector
		this->calcbounds(options);
	}
	//MPI-enabled constructors
#ifdef EMTG_MPI
	outerloop_NSGAII::outerloop_NSGAII(const EMTG::missionoptions& options, boost::mpi::environment* MPIEnvironmentin, boost::mpi::communicator* MPIWorldin)
	{
		//seed the random number generator
		this->RNG.seed(time(NULL));

		//calculate upper and lower bounds on the decision vector
		this->calcbounds(options);

		//set pointers to MPI communication objects
		this->MPIEnvironment = MPIEnvironmentin;
		this->MPIWorld = MPIWorldin;
	}
	outerloop_NSGAII::outerloop_NSGAII(int boost_random_seed, const EMTG::missionoptions& options, boost::mpi::environment* MPIEnvironmentin, boost::mpi::communicator* MPIWorldin)
	{
		//seed the random number generator
		this->RNG.seed(boost_random_seed);
		
		//calculate upper and lower bounds on the decision vector
		this->calcbounds(options);

		//set pointers to MPI communication objects
		this->MPIEnvironment = MPIEnvironmentin;
		this->MPIWorld = MPIWorldin;
	}
#endif

		//method to seed the initial population
	void outerloop_NSGAII::supply_pop(const std::vector< EMTG_outerloop_solution >& input_population)
	{
		this->this_generation = input_population;
	}

	//method to seed the archive
	void outerloop_NSGAII::supply_archive(std::vector< EMTG_outerloop_solution > const & input_archive)
	{
		this->archive_of_solutions = input_archive;
	}

	//method to retrieve the current generation
	const std::vector< EMTG_outerloop_solution >& outerloop_NSGAII::get_population()
	{
		return this->this_generation;
	}

	//method to retrieve the archive
	const std::vector< EMTG_outerloop_solution >& outerloop_NSGAII::get_archive()
	{
		return this->archive_of_solutions;
	}

	//methods to set GA parameters
	void outerloop_NSGAII::set_CR(double const & CR_in)
	{
		this->CR = CR_in;
	}

	void outerloop_NSGAII::set_mutationrate(double const & mu_in)
	{
		this->mu = mu_in;
	}

	void outerloop_NSGAII::set_tournament_size(unsigned int const & tournamentsize_in)
	{
		this->tournamentsize = tournamentsize_in;
	}

	void outerloop_NSGAII::set_stall_generations(unsigned int const & stall_generations_in)
	{
		this->stallmax = stall_generations_in;
	}

	void outerloop_NSGAII::set_populationsize(unsigned int const & populationsize_in)
	{
		this->popsize = populationsize_in;
	}

	void outerloop_NSGAII::set_max_generations(unsigned int const & max_generations_in)
	{
		this->genmax = max_generations_in;
	}

	void outerloop_NSGAII::set_tolfit(std::vector<double> const & tolfit_in)
	{
		this->tolfit = tolfit_in;
	}

	void outerloop_NSGAII::set_current_generation(const int& gen_in)
	{
		this->current_generation = gen_in;
	}

	void outerloop_NSGAII::set_objectives(const int& number_of_objectives_in, const std::vector<std::string>& objective_descriptions_in)
	{
		this->number_of_objectives = number_of_objectives_in;
		this->objective_descriptions = objective_descriptions_in;
	}

	//*****************private methods

	//NSGAII crowded tournament selection function (Matt)
	std::vector< EMTG_outerloop_solution >  outerloop_NSGAII::select(const std::vector<EMTG_outerloop_solution> &parent_population_in)
	{
		std::vector< EMTG_outerloop_solution > parent_pool;

		for (int i = 0; i < 2; i++)   // two rounds in tournament to create parent_pool of same size as parent_population
		{
			// shuffle population
			this->parent_population = parent_population_in;
#ifndef _STONEAGECplusplus
			std::shuffle(this->parent_population.begin(), this->parent_population.end(), this->RNG); // TODO: make sure shuffle works with RNG input
#else
			std::random_shuffle(parent_population.begin(), parent_population.end()); // TODO: make sure random_shuffle works with rand() function; modify to use mersenne twister
#endif

			// compare solutions to generate half of parent pool
			for (int j = 0; j < parent_population.size()/2; ++j)
			{
				if (parent_population[j].pareto_rank < parent_population[parent_population.size() - j - 1].pareto_rank)
					parent_pool.push_back(parent_population[j]);
				else if (parent_population[j].pareto_rank > parent_population[parent_population.size() - j - 1].pareto_rank)
					parent_pool.push_back(parent_population[parent_population.size() - j - 1]);
				else
				{
					// individuals belong to the same local front so use crowding distance measure to determine winner (choose individual w/ larger crowding_distnance)
					if (parent_population[j].crowding_distance > parent_population[parent_population.size() - j - 1].crowding_distance)
						parent_pool.push_back(parent_population[j]);
					else
						parent_pool.push_back(parent_population[parent_population.size() - j - 1]);
				}
			}
		}
		return parent_pool; // parent_pool (vector of NSGAII individuals) is returned
	}

	//GA crossover function (Matt)
	std::vector< EMTG_outerloop_solution >  outerloop_NSGAII::crossover_uniformInt(std::vector<EMTG_outerloop_solution> &parent_pool)
	{
		// uniform crossover for an integer GA crossover (crossover at each gene position)

		// initialize
		std::vector< EMTG_outerloop_solution > children_population;
		std::vector< int> X_child1;
		std::vector< int> X_child2;
		IntegerDistribution = boost::uniform_int<>(0, 1);

		// shuffle population
#ifndef _STONEAGECplusplus
		std::shuffle(this->parent_pool.begin(), this->parent_pool.end(), this->RNG);
#else
		std::random_shuffle(parent_pool.begin(), parent_pool.end()); // TODO: make sure random_shuffle works with rand() function; modify to use mersenne twister
#endif

		for (int i = 0; i < this->popsize/2; ++i)
		{
			X_child1.clear();
			X_child2.clear();
			for (int gene = 0 ; gene < parent_pool[i].X.size(); ++gene)
			{
				int coinflip = IntegerDistribution(RNG);
				if (coinflip == 0)
				{
					//if coinflip==0 child 1 chromosome receives gene from parent 1, child 2 receives gene from parent 2 
					X_child1.push_back(parent_pool[i].X[gene]);
					X_child2.push_back(parent_pool[this->popsize/2 + i].X[gene]);
				}
				else
				{
					//if coinflip==1 child 2 chromosome receives gene from parent 1, child 1 receives gene from parent 2
					X_child1.push_back(parent_pool[i].X[gene]);
					X_child2.push_back(parent_pool[this->popsize/2 + i].X[gene]);
				}
			}

			// add fitness value place holders
			
			// add children to offspring population
			children_population.push_back(EMTG_outerloop_solution(X_child1, this->number_of_objectives));
			children_population.push_back(EMTG_outerloop_solution(X_child2, this->number_of_objectives));
		}
		return children_population;  //next generation
	}

	//NSGA-II non-dominated sort (Matt)
	void outerloop_NSGAII::non_dominated_sort(std::vector<EMTG_outerloop_solution> &population)
	{
		//based on fast-non-dominated sort by Deb

		//declare sorting variables
		int obj_compare_lt_i;
		int obj_compare_lteq2_i;
		int obj_compare_lt_j;
		int obj_compare_lteq2_j;
		int front;
		std::vector< EMTG_outerloop_solution > local_front;

		this->nondominated_fronts.clear();

		//loop for calculating individuals comprising 1st nondominated front
		for (size_t i = 0; i < population.size(); ++i)
		{
			//initialize individual solution attributes
			population[i].dominated_by.clear(); 
			population[i].dominates.clear(); 
			population[i].pareto_rank = 0;
			population[i].ndom = population.size();

			//compare solution i with all other solutions to find if it is dominated
			for (size_t j = 0; j < population.size(); ++j)
			{
				if (i != j)
				{	 
					//check if i is dominated by j
					obj_compare_lt_i = 0;
					obj_compare_lteq2_i = 0;
					obj_compare_lt_j = 0;
					obj_compare_lteq2_j = 0;
					for (int obj_ind = 0; obj_ind < this->number_of_objectives; ++obj_ind)
					{
						obj_compare_lt_i = obj_compare_lt_i + population[j].compare_objective_lessthan(population[i], obj_ind, 0);
						obj_compare_lteq2_i = obj_compare_lteq2_i + population[j].compare_objective_lessthanorequalto(population[i], obj_ind, 1.0e-06);
						obj_compare_lt_j = obj_compare_lt_j + population[i].compare_objective_lessthan(population[j], obj_ind, 0);
						obj_compare_lteq2_j = obj_compare_lteq2_j + population[i].compare_objective_lessthanorequalto(population[j], obj_ind, 1.0e-06);
					}
					
					if ((obj_compare_lteq2_i == this->number_of_objectives) && (obj_compare_lt_i > 1))
						population[i].dominated_by.push_back(j); // individual i is dominated by j
					else if ((obj_compare_lteq2_j == this->number_of_objectives) && (obj_compare_lt_j > 1))
						population[i].dominates.push_back(j); // individual j is dominated by i
				}			
			}
			population[i].ndom = population[i].dominated_by.size(); // number of times individual i is dominated
			if (population[i].ndom == 0)
			{
				population[i].pareto_rank = 1;
				local_front.push_back(population[i]); // first non-dominated front
			}
		}

		// create top ranking local front
		this->nondominated_fronts.push_back(local_front);

		// loop for determining remaining local non-dominated fronts
		int frontcount = 0;
		while (nondominated_fronts[frontcount].size() != 0)
		{
			local_front.clear();
			// for each individual in the current front vector (@frontcount)
			for (int individual = 0; individual < this->nondominated_fronts[frontcount].size(); ++individual)
			{
				// for each individual in the 'dominates' vector
				for (int dom_by_ind = 0; dom_by_ind < this->nondominated_fronts[frontcount][individual].dominates.size(); ++dom_by_ind)
				{
					population[this->nondominated_fronts[frontcount][individual].dominates[dom_by_ind]].ndom = population[this->nondominated_fronts[frontcount][individual].dominates[dom_by_ind]].ndom  - 1; // reduce domination counter by 1
					if (population[this->nondominated_fronts[frontcount][individual].dominates[dom_by_ind]].ndom == 0)
					{
						population[this->nondominated_fronts[frontcount][individual].dominates[dom_by_ind]].pareto_rank = frontcount+2; // assign pareto rank to the individual
						local_front.push_back(population[this->nondominated_fronts[frontcount][individual].dominates[dom_by_ind]]);  // add solution object to current local front vector
					}
				}
			}
			frontcount = frontcount + 1;
			if (local_front.size() != 0)
				this->nondominated_fronts.push_back(local_front);
			else
				break;
		}

		// remove extra vector from this->nondominated_fronts
		if (this->nondominated_fronts.back().size() == 0)
			this->nondominated_fronts.pop_back();


		// assign crowding distance to member of each front
		for (int front = 0; front < this->nondominated_fronts.size(); ++front)
		{
			this->assign_crowding_distance(this->nondominated_fronts[front]);
		}

		// rebuild population 
		//TODO: ensure members of parent population were updated with pareto_rank value
		population.clear();
		for (int front = 0; front < this->nondominated_fronts.size(); ++front)
		{
			for (int individual = 0; individual < this->nondominated_fronts[front].size(); ++individual)
			{
				population.push_back(this->nondominated_fronts[front][individual]);
			}
		}
	}


	//NSGA-II assign crowding distance (Matt)
	void outerloop_NSGAII::assign_crowding_distance(std::vector<EMTG_outerloop_solution> &local_front)
	{
		// assign crowding distance measure described by Deb
		// updates .crowding_distance in each individual of local_front input

		// declare variables
		std::vector< std::vector< std::pair< double, int > > > obj_value;
		std::vector< std::pair< double, int > > obj_value_ind;
		std::vector< double > min_obj_value;
		std::vector< double > max_obj_value;

		// create vector of objective function value and index pairs (order from local_front order)
		for (int obj_ind = 0; obj_ind < this->number_of_objectives; ++obj_ind)
		{
			obj_value_ind.clear();
			for (int i = 0; i < local_front.size(); ++i)  // for each individual, i, in local_front
			{
				obj_value_ind.push_back(make_pair(local_front[i].fitness_values[obj_ind], i)); // create vector of objective function value and index pairs (order from local_front order)
			}
			obj_value.push_back(obj_value_ind);
		}

		// sort obj_value vector of vector pairs (.second of pair is associated with order of local_front)
		for (int obj_ind = 0; obj_ind < this->number_of_objectives; ++obj_ind)
		{
			std::sort(obj_value[obj_ind].begin(), obj_value[obj_ind].end());
		}

		// assign large crowding distance measure to solutions at boundaries for each objective function
		for (int obj_ind = 0; obj_ind < this->number_of_objectives; ++obj_ind)
		{
			local_front[obj_value[obj_ind][0].second].crowding_distance = 9.9e99; // individual with smallest objective value
			local_front[obj_value[obj_ind][local_front.size()-1].second].crowding_distance = 9.9e98; // individual with largest objective value (set smaller than individual w/ smallest obj value)

			// set min and max of objective function values for all objectives of the vector of local_front
			min_obj_value.push_back(local_front[obj_value[obj_ind][0].second].fitness_values[obj_ind]);
			max_obj_value.push_back(local_front[obj_value[obj_ind][local_front.size()-1].second].fitness_values[obj_ind]);
		}

		// assign crowding distance measure to remaining individuals according to "distance" of neighbors in objective space
		for (int i = 1; i < local_front.size()-1; ++i)  // for each individual, i, in local_front
		{
			for (int obj_ind = 0; obj_ind < this->number_of_objectives; ++obj_ind)
			{
				if (obj_ind == 0)
					local_front[obj_value[obj_ind][i].second].crowding_distance  = 0;  // initialize crowding_distance to 0
				double prev_crowd_dist = local_front[obj_value[obj_ind][i].second].crowding_distance;
				double obj_range = max_obj_value[obj_ind] - min_obj_value[obj_ind]; // range of objective function values for current objective
				double current_crowd_dist = (local_front[obj_value[obj_ind][i+1].second].fitness_values[obj_ind] - local_front[obj_value[obj_ind][i-1].second].fitness_values[obj_ind]) / obj_range;
				local_front[obj_value[obj_ind][i].second].crowding_distance = prev_crowd_dist + current_crowd_dist; // add crowding distance measure of current objective
			}
		}
	}

	//GA mutation function (Matt)

	void outerloop_NSGAII::mutate(std::vector<EMTG_outerloop_solution> &population)
	{
		// mutate genes of members of population based on mu value

		//declare/initialize variables
		this->DoubleDistribution = boost::uniform_real<>(0.0, 1.0);
		int genes_in_X = population[0].X.size();

		for (int individual = 0; individual < population.size(); ++individual)
		{
			for (int gene = 0; gene < genes_in_X; ++gene)
			{
				if (this->DoubleDistribution(RNG) < this->mu)
				{
					this->IntegerDistribution = boost::uniform_int<>(0, this->Xupperbounds[gene]); 
					population[individual].X[gene] = IntegerDistribution(RNG);
				}
			}
		}
	}


	//GA population_filter function (Matt)
	void outerloop_NSGAII::filter_population()
	{
		// reduce parent population to size this->popsize (size n), taking the top performing individuals
		// this is the elitism operation of the NSGAII

		// declare/initialize
		std::vector< EMTG_outerloop_solution > population;
		int slots_remaining = this->popsize; // number of remaining open slots of population 
		std::vector< std::pair < double, int > > crowd_dist_vec; // vector crowding_distance,index pairs

		this->parent_population.clear(); // clear parent population (will update at end of function)

		for (int front = 0; front < this->nondominated_fronts.size(); ++front)
		{
			if (slots_remaining > this->nondominated_fronts[front].size())
			{
				// add all members of current front
				for (int individual = 0; individual < this->nondominated_fronts[front].size(); ++individual)
				{
					population.push_back(this->nondominated_fronts[front][individual]);
					slots_remaining = this->popsize - population.size();
				}
			}
			else
			{
				// add least crowded individuals of current front until population is size n
				// create vector of crowding_distance, index pairs
				for (int i = 0; i < this->nondominated_fronts[front].size(); ++i)
				{
					crowd_dist_vec.push_back(make_pair(this->nondominated_fronts[front][i].crowding_distance, i));
				}

				//rearrange crowd_dist_vec in ascending order according to crowding distance
				std::sort(crowd_dist_vec.begin(), crowd_dist_vec.end()); 

				for (int i = 0; i < slots_remaining; ++i)
				{
					population.push_back(this->nondominated_fronts[front][crowd_dist_vec[crowd_dist_vec.size() - i -1].second]);
				}
				break;
			}
		}

		// update parent population (now size n)
		this->parent_population = population;
	}

	//*****************public methods

	//evaluate the population (Jacob)
	void outerloop_NSGAII::evaluatepop(const EMTG::missionoptions& options, const boost::ptr_vector<EMTG::Astrodynamics::universe>& Universe)
	{
		//the first step is to create subpopulation which does not appear in the archive, i.e. has not been evaluated
		//we really only need to keep track of their indices because later we will just place them back into the main population anyway
		std::vector< int > unevaluated_individuals_indices;
		bool NewSolution;

		if (options.outerloop_reevaluate_full_population) //this adds EVERYTHING to the list of things to be evaluated
		{
			for (size_t individual = 0; individual < this->popsize; ++individual)
				unevaluated_individuals_indices.push_back(individual);
		}
		else
		{
			for (size_t individual = 0; individual < this->popsize; ++individual)
			{
				//check if this solution has been evaluated before
				if ( !(this->this_generation[individual].description == "") )
					NewSolution = false;
				else
				{
					NewSolution = true;
					this->this_generation[individual].set_data_pointers((void*) &options, (void*) &(Universe));
					this->this_generation[individual].parse_outer_loop_decision_vector();

					for (size_t archived_solution = 0; archived_solution < this->archive_of_solutions.size(); ++archived_solution)
					{
						if ( this->archive_of_solutions[archived_solution].description == this->this_generation[individual].description )
						{
							if (!options.quiet_outerloop)
							{
								cout << "Solution " << this->this_generation[individual].description << " has already been evaluated with fitnesses [";

								for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
								{
									this->this_generation[individual].fitness_values[objective] = this->archive_of_solutions[archived_solution].fitness_values[objective];

									if (objective > 0)
										cout << ", ";

									cout << this->this_generation[individual].fitness_values[objective];
								}

								cout << "]" << endl;
							}
							NewSolution = false;

							this->this_generation[individual].timestamp = this->archive_of_solutions[archived_solution].timestamp;
							this->this_generation[individual].generation_found = this->archive_of_solutions[archived_solution].generation_found;
							this->this_generation[individual].Xinner = this->archive_of_solutions[archived_solution].Xinner;
							break;
						}
					}
				}

				if (NewSolution)
					unevaluated_individuals_indices.push_back(individual);
			}
		}

#ifdef EMTG_MPI
		//if MPI is enabled, conduct a parallel evaluation of the population
		//first announce that this processor is ready to go
		if (this->MPIWorld->rank() == 0)
		{
			if (!options.quiet_outerloop)
				std::cout << "Processor " << this->MPIWorld->rank() << " ready to broadcast population of " << unevaluated_individuals_indices.size() << " individuals." << std::endl;
		}
		
		//distribute the population into subvectors

		//vector of subpopulations
		std::vector< std::vector < EMTG_outerloop_solution > > subpopulations;

		//wrap the population into the vector of subpopulations
		if (this->MPIWorld->rank() == 0)
		{
			int populationcounter = 0;
			for (size_t Rank = 0; Rank < this->MPIWorld->size(); ++Rank)
			{
				std::vector< EMTG_outerloop_solution > temppop;
				for (size_t entry = 0; entry < unevaluated_individuals_indices.size() / this->MPIWorld->size(); ++entry)
				{
					temppop.push_back(this->this_generation[unevaluated_individuals_indices[populationcounter]]);
					++populationcounter;
				}
				if (Rank < unevaluated_individuals_indices.size() % this->MPIWorld->size())
				{
					temppop.push_back(this->this_generation[unevaluated_individuals_indices[populationcounter]]);
					++populationcounter;
				}
				subpopulations.push_back(temppop);
			}

			if (!options.quiet_outerloop)
				std::cout << "Processor " << this->MPIWorld->rank() << " is ready to scatter" << std::endl;
		}
		
		//then scatter the subvectors
		std::vector< EMTG_outerloop_solution > my_subpopulation;

		boost::mpi::scatter <vector <EMTG_outerloop_solution> > (*(this->MPIWorld), subpopulations, my_subpopulation, 0);

		//announce how many problems are to be evaluated by each processor
		if (!options.quiet_outerloop)
		{
			try
			{
				if (this->MPIWorld->rank() == 0)
					std::cout << "Attempting to scatter the population for generation " << this->current_generation << std::endl;
				std::cout << "Processor " << this->MPIWorld->rank() << " evaluating " << my_subpopulation.size() << " cases." << std::endl;
			}
			catch (int e)
			{
				
				std::cout << "Failure to scatter population! Exiting..." << std::endl;

				throw e;
			}
		}
		
		
		//evaluate the local subvector
		for (size_t individual = 0; individual < my_subpopulation.size(); ++individual)
		{
			if (!options.quiet_outerloop)
				std::cout << "Processor #" << this->MPIWorld->rank() << " " << my_subpopulation[individual].description << std::endl;
			my_subpopulation[individual].set_data_pointers((void*) &options, (void*) &(Universe));
			try
			{
				my_subpopulation[individual].evaluate();
			}
			catch (int e)
			{
				std::cout << "Processor #" << this->MPIWorld->rank() << " CRASHED evaluating " << my_subpopulation[individual].description << std::endl;
				std::cout << "Exception " << e << std::endl;
				for (int obj = 0; obj < my_subpopulation[individual].fitness_values.size(); ++ obj)
					my_subpopulation[individual].fitness_values[obj] = 1.0e+100;
			}
			my_subpopulation[individual].timestamp = std::time(NULL) - this->tstart;
			my_subpopulation[individual].generation_found = this->current_generation;
		}

		//gather the results
		try
		{
			if (!options.quiet_outerloop)
			{
				if (this->MPIWorld->rank() == 0)
					std::cout << "Attempting to gather the population for generation " << this->current_generation << std::endl;
			}
			boost::mpi::gather <vector <EMTG_outerloop_solution> > (*(this->MPIWorld), my_subpopulation, subpopulations, 0);
		}
		catch (int e)
		{
			std::cout << "Failure to gather population! Exiting..." << std::endl;

			throw e;
		}
#else
		//serial evaluation of the population
		for (size_t individual = 0; individual < unevaluated_individuals_indices.size(); ++individual)
		{
			this->this_generation[unevaluated_individuals_indices[individual]].set_data_pointers((void*) &options, (void*) &(Universe));
			this->this_generation[unevaluated_individuals_indices[individual]].evaluate();
			this->this_generation[unevaluated_individuals_indices[individual]].timestamp = std::time(NULL) - this->tstart;
			this->this_generation[unevaluated_individuals_indices[individual]].generation_found = this->current_generation;
		}
#endif
		
		//unwrap the results into the main population vector and add each solution to the archive IFF it is unique
#ifdef EMTG_MPI
		if (this->MPIWorld->rank() == 0)
		{
			int population_counter = 0;
			vector< EMTG_outerloop_solution > temppop;
			for (size_t Rank = 0; Rank < this->MPIWorld->size(); ++Rank)
			{
				for (size_t entry = 0; entry < subpopulations[Rank].size(); ++entry)
				{
					//store this solution
					this->this_generation[unevaluated_individuals_indices[population_counter]] = subpopulations[Rank][entry];
					++population_counter;

					//determine if this solution is a member of the archive. If not, add it
					//note this is NOT redundant. If the solution is in the archive but the new one is better, replace it.
					//If the old one is better, replace the new one with the old one.
					bool NewSolution = true;
					for (size_t ArchiveEntry = 0; ArchiveEntry < this->archive_of_solutions.size(); ++ArchiveEntry)
					{
						if (subpopulations[Rank][entry].description == this->archive_of_solutions[ArchiveEntry].description)
						{
							//if this solution is better than the one in the archive with the same name, over-write the archive
							if (subpopulations[Rank][entry].innerloop_fitness < this->archive_of_solutions[ArchiveEntry].innerloop_fitness)
							{
								std::cout << "Solution to " << subpopulations[Rank][entry].description << " is superior to its value in the archive. Overwriting." << std::endl;
								this->archive_of_solutions[ArchiveEntry] = subpopulations[Rank][entry];
							}
							//otherwise overwrite the current solution with what is in the archive
							else
								this->this_generation[unevaluated_individuals_indices[population_counter-1]] = this->archive_of_solutions[ArchiveEntry];

							NewSolution = false;
							break;
						}
					}
					if (NewSolution)
						this->archive_of_solutions.push_back(subpopulations[Rank][entry]);
				}
			}
		}
#else
		
		for (size_t entry = 0; entry < unevaluated_individuals_indices.size(); ++entry)
		{
			//determine if this solution is a member of the archive. If not, add it
			//note this is NOT redundant. While it is true that we do not evaluate any solution that already existed
			//in the archive, it is possible for us to evaluate several copies of a new solution at the same time
			//perhaps in the future there might be another (better?) way to do this by checking to see if any of the
			//identical solutions dominate the others. For now we'll just keep the first one.
			bool NewSolution = true;
			for (size_t ArchiveEntry = 0; ArchiveEntry < this->archive_of_solutions.size(); ++ArchiveEntry)
			{
				if (this->this_generation[unevaluated_individuals_indices[entry]].description == this->archive_of_solutions[ArchiveEntry].description)
				{
					NewSolution = false;
					break;
				}
			}
			if (NewSolution || options.outerloop_reevaluate_full_population)
				this->archive_of_solutions.push_back(this->this_generation[unevaluated_individuals_indices[entry]]);
		}
#endif
	}

	//generate a random population (Jacob)
	void outerloop_NSGAII::generatepop()
	{
		//first clear the starting population
		this->this_generation.clear();

		//clear the generation count
		this->current_generation = 0;
		this->stall_generation = 0;

		//now create a new population
		for (size_t individual = 0; individual < this->popsize; ++ individual)
		{
			vector<int> Xindividual;
			for (size_t gene = 0; gene < this->Xdescriptions.size(); ++gene)
			{
				Xindividual.push_back(this->random_integer[gene](this->RNG));
			}

			this->this_generation.push_back(EMTG_outerloop_solution(Xindividual, this->number_of_objectives));

			this->this_generation[individual].generation_found = 0;
			this->this_generation[individual].timestamp = 0;
		}
	}

	//read a population file (Jacob)
	void outerloop_NSGAII::readpop(std::string filename)
	{
		std::ifstream inputfile(filename.c_str());
		int linenumber = 0;
		char line_buffer[2048];
		std::vector<std::string> linecell;
		int number_of_genes = 0;

		if (!inputfile.is_open())
		{
			std::cout << "Cannot find population file " << filename << std::endl;
			throw 0;
		}
		
		while (!inputfile.eof()) 
		{
			++linenumber;
			char peek = inputfile.peek();
			if (peek == '#' || peek == '\r' || peek == '\n') 
			{
				//comment or blank line, do not parse
				inputfile.getline(line_buffer, 2048);	
			}
			else 
			{
				inputfile.getline(line_buffer, 2048);
				boost::split(linecell, line_buffer, boost::is_any_of(","));

				//the number of fields with "Gene" as the first four letters is the number of decision variables
				for (size_t entry = 0; entry < linecell.size(); ++entry)
				{
					if (linecell[entry].find("Gene") < 1024)
					{
						++number_of_genes;
					}
					else
						break;
				}

				if ( !(number_of_genes == this->Xdescriptions.size()) )
				{
					cout << "Number of genes in population file " << filename << " does not match number of genes for this mission script. Perhaps you have loaded the wrong file?" << std::endl;
					throw 0;
				}

				//skip the second header line
				inputfile.getline(line_buffer, 2048);
				
				break;
			}
		}
				
			

		//if we made it this far, it is time to read the rest of the input file and create the new population vector
		this->this_generation.clear();
		while (!inputfile.eof())
		{
			++linenumber;
			char peek = inputfile.peek();
			if (peek == '#' || peek == '\r' || peek == '\n') 
			{
				//comment or blank line, do not parse
				inputfile.getline(line_buffer, 2048);
				break;
			}
			else
			{
				//grab a line
				inputfile.getline(line_buffer, 2048);
				
				//split the line by commas
				boost::split(linecell, line_buffer, boost::is_any_of(","));

				if (linecell.size() > 1)
				{
					//the first N entries of the line are genes
					std::vector<int> Xtemp;
					for (size_t gene = 0; gene < number_of_genes; ++gene)
						Xtemp.push_back( atoi(linecell[gene].c_str()) );
				
					//add this vector to the population
					this->this_generation.push_back(EMTG_outerloop_solution(Xtemp, this->number_of_objectives));

					//fill in the description, generation found, and timestamp
					this->this_generation.back().description = linecell[number_of_genes];
					this->this_generation.back().generation_found = atoi(linecell[number_of_genes + 1].c_str());
					this->this_generation.back().timestamp = atoi(linecell[number_of_genes + 2].c_str());
				
					//fill in the objective values
					for (size_t objective = 0; objective < this->number_of_objectives; ++objective)
						this->this_generation.back().fitness_values[objective] = atof(linecell[number_of_genes + 3 + objective].c_str());

					//fill in the inner-loop vector
					for (size_t index = number_of_genes + 3 + this->number_of_objectives; index < linecell.size(); ++index)
						this->this_generation.back().Xinner.push_back(atof(linecell[index].c_str()));
				}
			}
		}
	}

	//write a population file (Jacob)
	void outerloop_NSGAII::writepop(std::string filename)
	{
#ifdef EMTG_MPI
		//if we are in parallel mode and NOT the head node, don't write anything
		if (!(this->MPIWorld->rank() == 0))
			return;
#endif

		ofstream outputfile(filename.c_str(), ios::trunc);

		outputfile << "#GA population file" << std::endl;
		outputfile << "##Written by EMTG_v8 core program compiled " << __DATE__ << " " << __TIME__ << std::endl;
		outputfile << std::endl;
		
		//output header rows
		for (size_t gene = 0; gene < this->Xdescriptions.size(); ++gene)
			outputfile << "Gene " << gene << ",";
		outputfile << "Description, Generation found, timestamp";
		for (size_t objective = 0; objective < this->number_of_objectives; ++objective)
			outputfile << "," << this->objective_descriptions[objective];
		outputfile << ",Inner-loop decision vector";
		outputfile << std::endl;
		for (size_t gene = 0; gene < this->Xdescriptions.size(); ++gene)
			outputfile << Xdescriptions[gene] << ",";
		outputfile << endl;

		outputfile.close();

		//output every member of the population
		for (size_t member = 0; member < this->this_generation.size(); ++member)
			this->this_generation[member].write_csv_line(filename);

		
	}

	//read an archive file (Jacob)
	void outerloop_NSGAII::read_archive(std::string filename)
	{
		std::ifstream inputfile(filename.c_str());
		int linenumber = 0;
		char line_buffer[2048];
		std::vector<std::string> linecell;
		int number_of_genes = 0;

		if (!inputfile.is_open())
		{
			std::cout << "Cannot find archive file " << filename << std::endl;
			throw 0;
		}
		
		while (!inputfile.eof()) 
		{
			++linenumber;
			char peek = inputfile.peek();
			if (peek == '#' || peek == '\r' || peek == '\n') 
			{
				//comment or blank line, do not parse
				inputfile.getline(line_buffer, 2048);	
			}
			else 
			{
				inputfile.getline(line_buffer, 2048);
				boost::split(linecell, line_buffer, boost::is_any_of(","));

				//the number of fields with "Gene" as the first four letters is the number of decision variables
				for (size_t entry = 0; entry < linecell.size(); ++entry)
				{
					if (linecell[entry].find("Gene") < 1024)
					{
						++number_of_genes;
					}
					else
						break;
				}

				if ( !(number_of_genes == this->Xdescriptions.size()) )
				{
					cout << "Number of genes in archive file " << filename << " does not match number of genes for this mission script. Perhaps you have loaded the wrong file?" << std::endl;
					throw 0;
				}

				//skip the second header line
				inputfile.getline(line_buffer, 2048);

				break;
			}
		}
				
			

		//if we made it this far, it is time to read the rest of the input file and create the new population vector
		this->archive_of_solutions.clear();
		while (!inputfile.eof())
		{
			++linenumber;
			char peek = inputfile.peek();
			if (peek == '#' || peek == '\r' || peek == '\n') 
			{
				//comment or blank line, do not parse
				inputfile.getline(line_buffer, 2048);
				break;
			}
			else
			{
				//grab a line
				inputfile.getline(line_buffer, 2048);
				
				//split the line by commas
				boost::split(linecell, line_buffer, boost::is_any_of(","));

				if (linecell.size() > 1)
				{
					//the first N entries of the line are genes
					std::vector<int> Xtemp;
					for (size_t gene = 0; gene < number_of_genes; ++gene)
						Xtemp.push_back( atoi(linecell[gene].c_str()) );
				
					//add this vector to the population
					this->archive_of_solutions.push_back(EMTG_outerloop_solution(Xtemp, this->number_of_objectives));

					//fill in the description, generation found, and timestamp
					this->archive_of_solutions.back().description = atoi(linecell[number_of_genes].c_str());
					this->archive_of_solutions.back().generation_found = atoi(linecell[number_of_genes + 1].c_str());
					this->archive_of_solutions.back().timestamp = atoi(linecell[number_of_genes + 2].c_str());
				
					//fill in the objective values
					for (size_t objective = 0; objective < this->number_of_objectives; ++objective)
						this->archive_of_solutions.back().fitness_values[objective] = atof(linecell[number_of_genes + 3 + objective].c_str());

					//fill in the inner-loop vector
					for (size_t index = number_of_genes + 3 + this->number_of_objectives; index < linecell.size(); ++index)
						this->archive_of_solutions.back().Xinner.push_back(atof(linecell[index].c_str()));
				}
			}
		}
	}

	//write an archive file (Jacob)
	void outerloop_NSGAII::write_archive(std::string filename)
	{
#ifdef EMTG_MPI
		//if we are in parallel mode and NOT the head node, don't write anything
		if (!(this->MPIWorld->rank() == 0))
			return;
#endif

		ofstream outputfile(filename.c_str(), ios::trunc);

		outputfile << "#GA archive file" << std::endl;
		outputfile << "##Written by EMTG_v8 core program compiled " << __DATE__ << " " << __TIME__ << std::endl;
		outputfile << std::endl;
		
		//output header rows
		for (size_t gene = 0; gene < this->Xdescriptions.size(); ++gene)
			outputfile << "Gene " << gene << ",";
		outputfile << "Description, Generation found, timestamp";
		for (size_t objective = 0; objective < this->number_of_objectives; ++objective)
			outputfile << "," << this->objective_descriptions[objective];
		outputfile << ",Inner-loop decision vector";
		outputfile << std::endl;
		for (size_t gene = 0; gene < this->Xdescriptions.size(); ++gene)
			outputfile << Xdescriptions[gene] << ",";
		
		outputfile << endl;

		//output every member of the population
		for (size_t member = 0; member < this->archive_of_solutions.size(); ++member)
			this->archive_of_solutions[member].write_csv_line(filename);

		outputfile.close();
	}

	//run the main GA evolution loop (Matt and Jacob)
	void outerloop_NSGAII::evolve(const EMTG::missionoptions& options, const boost::ptr_vector<EMTG::Astrodynamics::universe>& Universe)
	{
		std::srand(time(NULL));
		// evolves this_generation population towards Pareto front (main NSGA-II function)

		//Initialize
		//************
		//1. generate random initial parent population (size n)
		//1a. evaluate objectives for initial parent population
		this->evaluatepop(options, Universe); //evaluates this->this_generation vector of solution indiviudals

		//1b. assign pareto rank and crowding distance value to initial parent population using non-dominated sorting
		this->parent_population = this->this_generation; //random initial parent population of size n based on copy of this->this_generation
		std::stringstream popfilestream;
		popfilestream << options.working_directory << "//NSGAII_initial_population.csv";
		this->writepop(popfilestream.str());
		popfilestream.clear();
		this->write_archive(options.working_directory + "//NSGAII_archive.csv");

		this->non_dominated_sort(this->parent_population); // assigns local non-dominated front value to parent_population members

		//2. generate initial child population from initial parent population (size n)
		//2a. crowded tournment selection to generate parent pool
		this->parent_pool = this->select(this->parent_population);
		
		//2b. apply crossover to generate offspring from parent pool
		this->children_population = this->crossover_uniformInt(this->parent_pool);

		//2c. mutate offspring to finalize initial child population
		this->mutate(this->children_population);

		//2d. evaluate objectives for initial child population
		this->this_generation = this->children_population;
		this->evaluatepop(options, Universe);
		this->write_archive(options.working_directory + "//NSGAII_archive.csv");

		this->children_population = this->this_generation;
		this->current_generation = 1;
		

		// Main loop
		//**************
		while (this->current_generation <= this->genmax)
		{
			//3. combine parent and child populations to create combined population
			for (int individual = 0; individual < this->children_population.size(); ++individual)
			{
				this->parent_population.push_back(children_population[individual]);
			}

			//4. assign fitness to parent population (size 2n), using non-dominated sorting
			this->non_dominated_sort(this->parent_population); // assigns local non-dominated front value to parent_population members

			//5. generate new parent population from combined population by filling N slots with the best N designs from combine population
			this->filter_population(); //updates parent population
			this->this_generation = this->parent_population;

			// print the current population to a file
			std::stringstream popfilestream;
			popfilestream << options.working_directory << "//NSGAII_population_gen_" << this->current_generation << ".csv";
			this->writepop(popfilestream.str());

			//6. generate new offspring via genetic operators
			//6a. crowded tournment selection to generate parent pool from parent population
			this->parent_pool = this->select(this->parent_population);

			//6b. apply crossover to generate offspring from parent pool
			this->children_population = this->crossover_uniformInt(this->parent_pool);

			//6c. mutate offspring to finalize new child population
			this->mutate(this->children_population); // update child population

			//6d. evaluate objectives for new child population
			this->last_generation = this->parent_population;
			this->this_generation = this->children_population;
			this->evaluatepop(options, Universe);
			this->children_population = this->this_generation;  

			//7. write out the current archive in case we need to warm-start another GA later
			this->write_archive(options.working_directory + "//NSGAII_archive.csv");

			++this->current_generation;  // increase generation counter
		}
	}



	//reset the GA (Jacob)
	void outerloop_NSGAII::reset()
	{
		//clear the generation count
		this->current_generation = 0;
		this->stall_generation = 0;
		
		//clear the populations
		this->this_generation.clear();
		this->last_generation.clear();
		
		//clear the archive
		this->archive_of_solutions.clear();
	}

	//GA upper and lower bounds calculation (Jacob)
	void outerloop_NSGAII::calcbounds(const EMTG::missionoptions& options)
	{
		//set up the objectives
		std::vector< std::string > objective_menu_descriptions;
		objective_menu_descriptions.push_back("BOL power at 1 AU (kW)");
		objective_menu_descriptions.push_back("Launch epoch (MJD)");
		objective_menu_descriptions.push_back("Flight time (days)");
		objective_menu_descriptions.push_back("Thruster preference");
		objective_menu_descriptions.push_back("Number of thrusters");
		objective_menu_descriptions.push_back("Launch vehicle preference");
		objective_menu_descriptions.push_back("Delivered mass to final target");
		objective_menu_descriptions.push_back("Final journey mass increment (for maximizing sample return)");
		objective_menu_descriptions.push_back("First journey departure C3 (km^2/s^2)");
		objective_menu_descriptions.push_back("Final journey arrival C3 (km^2/s^2)");

		for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
		{
			this->objective_descriptions.push_back(objective_menu_descriptions[options.outerloop_objective_function_choices[objective]]);
		}
		this->number_of_objectives = objective_descriptions.size();

		//construct the outer-loop decision vector
		//for each decision vector entry we need upper bounds, lower bounds, a description, and a random integer generator
		//order of entries is:
		//power
		//launch window open date
		//flight time upper bound
		//thruster type
		//number of thrusters
		//launch vehicle
		//for each journey:
		//	destination
		//  flyby 1
		//  flyby 2
		//  ...
		//  flyby n
		//
		//note that any of these can be turned on or off

		//vary global mission settings
		if (options.outerloop_vary_power && options.outerloop_power_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_power_choices.size() - 1);
			this->Xdescriptions.push_back("BOL Power at 1 AU (kW)");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying power" << std::endl;
		}
		if (options.outerloop_vary_launch_epoch && options.outerloop_launch_epoch_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_launch_epoch_choices.size() - 1);
			this->Xdescriptions.push_back("Launch window open date (MJD)");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying launch window open date" << std::endl;
		}
		if (options.outerloop_vary_flight_time_upper_bound && options.outerloop_flight_time_upper_bound_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_flight_time_upper_bound_choices.size() - 1);
			this->Xdescriptions.push_back("Flight time upper bound (d)");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying flight time upper bound" << std::endl;
		}
		if (options.outerloop_vary_thruster_type && options.outerloop_thruster_type_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_thruster_type_choices.size() - 1);
			this->Xdescriptions.push_back("Thruster type");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying thruster type" << std::endl;
		}
		if (options.outerloop_vary_number_of_thrusters && options.outerloop_number_of_thrusters_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_number_of_thrusters_choices.size() - 1);
			this->Xdescriptions.push_back("Number of thrusters");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying number of thrusters" << std::endl;
		}
		if (options.outerloop_vary_launch_vehicle && options.outerloop_launch_vehicle_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_launch_vehicle_choices.size() - 1);
			this->Xdescriptions.push_back("Launch vehicle");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying launch vehicle" << std::endl;
		}
		if (options.outerloop_vary_departure_C3 && options.outerloop_departure_C3_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_departure_C3_choices.size() - 1);
			this->Xdescriptions.push_back("Departure C3 upper bound");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying first journey departure C3" << std::endl;
		}
		if (options.outerloop_vary_arrival_C3 && options.outerloop_arrival_C3_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_arrival_C3_choices.size() - 1);
			this->Xdescriptions.push_back("Arrival C3 upper bound");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			if (!options.quiet_outerloop)
				std::cout << "varying last journey arrival C3" << std::endl;
		}

		//options for each journey
		for (size_t j = 0; j < options.number_of_journeys; ++j)
		{
			//journey destination
			if (options.outerloop_vary_journey_destination[j] && options.outerloop_journey_destination_choices[j].size() > 1)
			{
				this->Xlowerbounds.push_back(0);
				this->Xupperbounds.push_back(options.outerloop_journey_destination_choices[j].size() - 1);
				std::stringstream descriptionstream;
				descriptionstream << "Journey " << j << " destination";
				this->Xdescriptions.push_back(descriptionstream.str());
				this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
				if (!options.quiet_outerloop)
					std::cout << "varying journey " << j << " destination" << std::endl;
			}

			//flyby sequence using Englander and Conway "null gene" method
			if (options.outerloop_vary_journey_flyby_sequence[j] && options.outerloop_journey_flyby_sequence_choices[j].size() > 1)
			{
				for (size_t p = 0; p < options.outerloop_journey_maximum_number_of_flybys[j]; ++p)
				{
					this->Xlowerbounds.push_back(0);
					this->Xupperbounds.push_back( (options.outerloop_journey_flyby_sequence_choices[j].size() - 1) * 2);
					std::stringstream descriptionstream;
					descriptionstream << "Journey " << j << " potential flyby target " << p;
					this->Xdescriptions.push_back(descriptionstream.str());
					this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
					if (!options.quiet_outerloop)
						std::cout << "varying journey " << j << " potential flyby target " << p << std::endl;
				}
			}
		}
	}
}

