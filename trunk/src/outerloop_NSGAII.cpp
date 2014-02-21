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
			options.launch_window_open_date = options.outerloop_launch_epoch_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_LD" << (int) options.launch_window_open_date;
		}
		if (options.outerloop_vary_flight_time_upper_bound && options.outerloop_flight_time_upper_bound_choices.size() > 1)
		{
			options.total_flight_time_bounds[1] = options.outerloop_flight_time_upper_bound_choices[X[Xindex]];

			if (X[Xindex] > 0)
				options.total_flight_time_bounds[0] = options.outerloop_flight_time_upper_bound_choices[X[Xindex] - 1];
			else
				options.total_flight_time_bounds[0] = 0.0;

			++Xindex;
			descriptionstream << "_" << (int)(options.total_flight_time_bounds[1]) << "d";
		}
		if (options.outerloop_vary_thruster_type && options.outerloop_thruster_type_choices.size() > 1)
		{
			options.engine_type = options.outerloop_thruster_type_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_" << options.thruster_names[options.outerloop_thruster_type_choices[X[Xindex]]];
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
			descriptionstream << "_" << options.LV_names[options.outerloop_launch_vehicle_choices[X[Xindex]]];
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
					journey_sequence.push_back(options.destination_list[j][0]);
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
			vector<int> phase_type(journey_sequence.size() - 1, options.mission_type);
			options.phase_type.push_back(phase_type);
		}
		options.description = descriptionstream.str();
		this->description = options.description;

		return options;
	}

	void EMTG_outerloop_solution::evaluate()
	{
		//Step 1: create the options structure
		EMTG::missionoptions options = this->parse_outer_loop_decision_vector();

		//Step 2: instantiate a mission
		EMTG::mission TrialMission(&options, *Universe);

		//Step 3: run the inner-loop
		TrialMission.output_mission_tree(options.working_directory + "//" + this->description + ".emtgtree");
		cout << "Optimizing mission: " << TrialMission.options.description << endl;
		TrialMission.optimize();

		//Step 4: extract the fitness values
		//if the trial mission has no feasible solution then set the fitness values very high
		if (TrialMission.number_of_solutions == 0)
		{
			for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
				this->fitness_values[objective] = 1.0e+100;
			return;
		}
		else
		{
			this->Xinner = TrialMission.Xopt;
			TrialMission.extract_objective_function_values(this->fitness_values);
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

	//GA selection function (Matt)
	void outerloop_NSGAII::select()
	{
	}

	//GA mutation function (Matt)
	void outerloop_NSGAII::mutate()
	{
	}

	//GA crossover function (Matt)
	void outerloop_NSGAII::crossover()
	{
	}

	//NSGA-II non-dominated sort (Matt)
	void outerloop_NSGAII::non_dominated_sort()
	{
		
	}

	//*****************public methods

	//evaluate the population (Jacob)
	void outerloop_NSGAII::evaluatepop(const EMTG::missionoptions& options, const boost::ptr_vector<EMTG::Astrodynamics::universe>& Universe)
	{
		//the first step is to create subpopulation which does not appear in the archive, i.e. has not been evaluated
		//we really only need to keep track of their indices because later we will just place them back into the main population anyway
		std::vector< int > unevaluated_individuals_indices;
		bool NewSolution;
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
						cout << "Solution " << this->this_generation[individual].description << " has already been evaluated with fitnesses [";

						for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
						{
							this->this_generation[individual].fitness_values[objective] = this->archive_of_solutions[archived_solution].fitness_values[objective];

							if (objective > 0)
								cout << ", ";

							cout << this->this_generation[individual].fitness_values[objective];
						}

						cout << "]" << endl;
						NewSolution = false;
						break;
					}
				}
			}

			if (NewSolution)
				unevaluated_individuals_indices.push_back(individual);
		}

#ifdef EMTG_MPI
		//if MPI is enabled, conduct a parallel evaluation of the populatio
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
				if (Rank < this->popsize % this->MPIWorld->size())
				{
					temppop.push_back(this->this_generation[unevaluated_individuals_indices[populationcounter]]);
					++populationcounter;
				}
				subpopulations.push_back(temppop);
			}
		}
		
		//then scatter the subvectors
		std::vector< EMTG_outerloop_solution > my_subpopulation;

		boost::mpi::scatter(*(this->MPIWorld), subpopulations, &my_subpopulation, this->MPIWorld->size(), 0);
		
		//evaluate the local subvector
		for (size_t individual = 0; individual < my_subpopulation.size(); ++individual)
		{
			std::cout << "Proc #" << this->MPIWorld->rank() << " " << my_subpopulation[individual].description << std::endl;
			my_subpopulation[individual].set_data_pointers((void*) &options, (void*) &(Universe));
			my_subpopulation[individual].evaluate();
			my_subpopulation[individual].timestamp = std::time(NULL) - this->tstart;
			my_subpopulation[individual].generation_found = this->current_generation;
		}

		//gather the results
		boost::mpi::gather(*(this->MPIWorld), my_subpopulation, subpopulations, 0);
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
					//note this is NOT redundant. While it is true that we do not evaluate any solution that already existed
					//in the archive, it is possible for us to evaluate several copies of a new solution at the same time
					//perhaps in the future there might be another (better?) way to do this by checking to see if any of the
					//identical solutions dominate the others. For now we'll just keep the first one.
					bool NewSolution = true;
					for (size_t ArchiveEntry = 0; ArchiveEntry < this->archive_of_solutions.size(); ++ArchiveEntry)
					{
						if (subpopulations[Rank][entry].description == this->archive_of_solutions[ArchiveEntry].description)
						{
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
			if (NewSolution)
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
	void outerloop_NSGAII::evolve(const EMTG::missionoptions& options, const EMTG::Astrodynamics::universe& TheUniverse)
	{
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
			std::cout << "varying power" << std::endl;
		}
		if (options.outerloop_vary_launch_epoch && options.outerloop_launch_epoch_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_launch_epoch_choices.size() - 1);
			this->Xdescriptions.push_back("Launch window open date (MJD)");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			std::cout << "varying launch window open date" << std::endl;
		}
		if (options.outerloop_vary_flight_time_upper_bound && options.outerloop_flight_time_upper_bound_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_flight_time_upper_bound_choices.size() - 1);
			this->Xdescriptions.push_back("Flight time upper bound (d)");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			std::cout << "varying flight time upper bound" << std::endl;
		}
		if (options.outerloop_vary_thruster_type && options.outerloop_thruster_type_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_thruster_type_choices.size() - 1);
			this->Xdescriptions.push_back("Thruster type");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			std::cout << "varying thruster type" << std::endl;
		}
		if (options.outerloop_vary_number_of_thrusters && options.outerloop_number_of_thrusters_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_number_of_thrusters_choices.size() - 1);
			this->Xdescriptions.push_back("Number of thrusters");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			std::cout << "varying number of thrusters" << std::endl;
		}
		if (options.outerloop_vary_launch_vehicle && options.outerloop_launch_vehicle_choices.size() > 1)
		{
			this->Xlowerbounds.push_back(0);
			this->Xupperbounds.push_back(options.outerloop_launch_vehicle_choices.size() - 1);
			this->Xdescriptions.push_back("Launch vehicle");
			this->random_integer.push_back(boost::uniform_int<>(this->Xlowerbounds[this->Xlowerbounds.size() - 1], this->Xupperbounds[this->Xupperbounds.size() - 1]));
			std::cout << "varying launch vehicle" << std::endl;
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
					std::cout << "varying journey " << j << " potential flyby target " << p << std::endl;
				}
			}
		}
	}
}

