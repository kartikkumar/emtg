//outer-loop single objective genetic algorithm (SGA) class
//can solve full systems optimization problem for single objective
//integer encoding
//J. Englander 5-10-2014
//based heavily on outer-loop NSGAII by M. Vavrina and J. Englander

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "outerloop_SGA.h"
#include "mission.h"

#include "boost/algorithm/string.hpp"

namespace GeneticAlgorithm
{
	//constructors
	outerloop_SGA::outerloop_SGA(const EMTG::missionoptions& options)
	{
		//seed the random number generator
		this->RNG.seed(time(NULL));

		//calculate upper and lower bounds on the decision vector
		this->calcbounds(options);
	}
	outerloop_SGA::outerloop_SGA(int boost_random_seed, const EMTG::missionoptions& options)
	{
		//seed the random number generator
		this->RNG.seed(boost_random_seed);

		//calculate upper and lower bounds on the decision vector
		this->calcbounds(options);
	}
	//MPI-enabled constructors
#ifdef EMTG_MPI
	outerloop_SGA::outerloop_SGA(const EMTG::missionoptions& options, boost::mpi::environment* MPIEnvironmentin, boost::mpi::communicator* MPIWorldin)
	{
		//seed the random number generator
		this->RNG.seed(time(NULL));

		//calculate upper and lower bounds on the decision vector
		this->calcbounds(options);

		//set pointers to MPI communication objects
		this->MPIEnvironment = MPIEnvironmentin;
		this->MPIWorld = MPIWorldin;
	}
	outerloop_SGA::outerloop_SGA(int boost_random_seed, const EMTG::missionoptions& options, boost::mpi::environment* MPIEnvironmentin, boost::mpi::communicator* MPIWorldin)
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
	void outerloop_SGA::supply_pop(const std::vector< EMTG_outerloop_solution >& input_population)
	{
		this->this_generation = input_population;
	}

	//method to seed the archive
	void outerloop_SGA::supply_archive(std::vector< EMTG_outerloop_solution > const & input_archive)
	{
		this->archive_of_solutions = input_archive;
	}

	//method to retrieve the current generation
	const std::vector< EMTG_outerloop_solution >& outerloop_SGA::get_population()
	{
		return this->this_generation;
	}

	//method to retrieve the archive
	const std::vector< EMTG_outerloop_solution >& outerloop_SGA::get_archive()
	{
		return this->archive_of_solutions;
	}

	//methods to set GA parameters
	void outerloop_SGA::set_CR(double const & CR_in)
	{
		this->CR = CR_in;
	}

	void outerloop_SGA::set_mutationrate(double const & mu_in)
	{
		this->mu = mu_in;
	}

	void outerloop_SGA::set_elitecount(unsigned int const & elitecount_in)
	{
		this->elitecount = elitecount_in;
	}

	void outerloop_SGA::set_tournament_size(unsigned int const & tournamentsize_in)
	{
		this->tournamentsize = tournamentsize_in;
	}

	void outerloop_SGA::set_stall_generations(unsigned int const & stall_generations_in)
	{
		this->stallmax = stall_generations_in;
	}

	void outerloop_SGA::set_populationsize(unsigned int const & populationsize_in)
	{
		this->popsize = populationsize_in;
	}

	void outerloop_SGA::set_max_generations(unsigned int const & max_generations_in)
	{
		this->genmax = max_generations_in;
	}

	void outerloop_SGA::set_tolfit(std::vector<double> const & tolfit_in)
	{
		this->tolfit = tolfit_in;
	}

	void outerloop_SGA::set_current_generation(const int& gen_in)
	{
		this->current_generation = gen_in;
	}

	void outerloop_SGA::set_objectives(const int& number_of_objectives_in, const std::vector<std::string>& objective_descriptions_in)
	{
		this->number_of_objectives = number_of_objectives_in;
		this->objective_descriptions = objective_descriptions_in;
	}

	//*****************private methods
	//GA crossover function
	void outerloop_SGA::crossover_uniformInt(EMTG_outerloop_solution* Parent1, EMTG_outerloop_solution* Parent2, int index_in_child_population)
	{
		// uniform crossover for an integer GA crossover (crossover at each gene position)
		// based on original integer GA concept from EMTGv6
		int ChromosomeLength = Parent1->X.size();
		vector<int> Xbaby(ChromosomeLength);

		//choose two crossover points
		int coinflip = IntegerDistribution(RNG);
		for (int gene = 0; gene < ChromosomeLength; ++gene)
		if (coinflip == 0)
			Xbaby[gene] = Parent1->X[gene];
		else
			Xbaby[gene] = Parent2->X[gene];

		//put the baby into the child population
		this->children_population[index_in_child_population] = EMTG_outerloop_solution(Xbaby, 1);
	}


	//GA mutation function
	void outerloop_SGA::mutate(EMTG_outerloop_solution* Individual)
	{
		// mutate genes of members of population based on mu value

		//declare/initialize variables
		this->DoubleDistribution = boost::uniform_real<>(0.0, 1.0);
		int ChromosomeLength = Individual->X.size();
		bool mutated_flag = false;

		for (int gene = 0; gene < ChromosomeLength; ++gene)
		{
			if (this->DoubleDistribution(RNG) < this->mu)
			{
				mutated_flag = true;
				this->IntegerDistribution = boost::uniform_int<>(0, this->Xupperbounds[gene]);
				Individual->X[gene] = this->random_integer[gene](this->RNG);
			}
		}

		if (mutated_flag)
			Individual->clear_innerloop_solution();
	}

	//GA selection function
	EMTG_outerloop_solution outerloop_SGA::select()
	{
		//first fill the parent pool
		
		IntegerDistribution = boost::uniform_int<>(0, this->popsize - 1);
		
		std::vector<EMTG_outerloop_solution> parent_pool;
		std::vector<int> parent_pool_indices;
		int candidate_parent_index;

		parent_pool_indices.push_back(IntegerDistribution(this->RNG));
		parent_pool.push_back(this->this_generation[parent_pool_indices[0]]);

		for (int k = 1; k < this->tournamentsize; ++k)
		{
			bool flag; 
			do
			{
				flag = false;
				candidate_parent_index = IntegerDistribution(RNG);

				for (int q = 0; q < parent_pool_indices.size(); ++q) //check to see if we're in the parent pool
				{
					if (candidate_parent_index == parent_pool_indices[q]) //if this index is already in the parent pool, turn on the flag
						flag = 1;
				}
			} while (flag);
			parent_pool.push_back(this->this_generation[candidate_parent_index]);
		}

		//sort the parent pool
		std::sort(parent_pool.begin(), parent_pool.end());

		//return the best in the sort
		return parent_pool[0];
	}

	//*****************public methods

	//evaluate the population (Jacob)
	void outerloop_SGA::evaluatepop(const EMTG::missionoptions& options, const boost::ptr_vector<EMTG::Astrodynamics::universe>& Universe)
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
				if (!(this->this_generation[individual].description == ""))
					NewSolution = false;
				else
				{
					NewSolution = true;
					this->this_generation[individual].set_data_pointers((void*)&options, (void*)&(Universe));
					this->this_generation[individual].parse_outer_loop_decision_vector();

					for (size_t archived_solution = 0; archived_solution < this->archive_of_solutions.size(); ++archived_solution)
					{
						if (this->archive_of_solutions[archived_solution].description == this->this_generation[individual].description)
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

		boost::mpi::scatter <vector <EMTG_outerloop_solution> >(*(this->MPIWorld), subpopulations, my_subpopulation, 0);

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
			my_subpopulation[individual].set_data_pointers((void*)&options, (void*)&(Universe));
			try
			{
				my_subpopulation[individual].evaluate();
			}
			catch (int e)
			{
				std::cout << "Processor #" << this->MPIWorld->rank() << " CRASHED evaluating " << my_subpopulation[individual].description << std::endl;
				std::cout << "Exception " << e << std::endl;
				for (int obj = 0; obj < my_subpopulation[individual].fitness_values.size(); ++obj)
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
			boost::mpi::gather <vector <EMTG_outerloop_solution> >(*(this->MPIWorld), my_subpopulation, subpopulations, 0);
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
			this->this_generation[unevaluated_individuals_indices[individual]].set_data_pointers((void*)&options, (void*)&(Universe));
			try
			{
				this->this_generation[unevaluated_individuals_indices[individual]].evaluate();
			}
			catch (int e)
			{
				std::cout << "CRASHED evaluating " << this->this_generation[unevaluated_individuals_indices[individual]].description << std::endl;
				std::cout << "Exception " << e << std::endl;
				for (int obj = 0; obj < this->this_generation[unevaluated_individuals_indices[individual]].fitness_values.size(); ++obj)
					this->this_generation[unevaluated_individuals_indices[individual]].fitness_values[obj] = 1.0e+100;
			}

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
					//If the solution is in the archive but the new one is better, replace it
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
								this->this_generation[unevaluated_individuals_indices[population_counter - 1]] = this->archive_of_solutions[ArchiveEntry];

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
			//If the solution is in the archive but the new one is better, replace it
			//If the old one is better, replace the new one with the old one.
			bool NewSolution = true;
			for (size_t ArchiveEntry = 0; ArchiveEntry < this->archive_of_solutions.size(); ++ArchiveEntry)
			{
				if (this->this_generation[unevaluated_individuals_indices[entry]].description == this->archive_of_solutions[ArchiveEntry].description)
				{

					//if this solution is better than the one in the archive with the same name, over-write the archive
					if (this->this_generation[unevaluated_individuals_indices[entry]].innerloop_fitness < this->archive_of_solutions[ArchiveEntry].innerloop_fitness)
					{
						std::cout << "Solution to " << this->this_generation[unevaluated_individuals_indices[entry]].description << " is superior to its value in the archive. Overwriting." << std::endl;
						this->archive_of_solutions[ArchiveEntry] = this->this_generation[unevaluated_individuals_indices[entry]];
					}
					//otherwise overwrite the current solution with what is in the archive
					else
						this->this_generation[unevaluated_individuals_indices[entry]] = this->archive_of_solutions[ArchiveEntry];

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
	void outerloop_SGA::generatepop()
	{
		//first clear the starting population
		this->this_generation.clear();

		//clear the generation count
		this->current_generation = 0;
		this->stall_generation = 0;

		//now create a new population
		for (size_t individual = 0; individual < this->popsize; ++individual)
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
	void outerloop_SGA::readpop(std::string filename)
	{
		std::ifstream inputfile(filename.c_str());
		int linenumber = 0;
		char line_buffer[65536];
		std::vector<std::string> linecell;
		int number_of_genes = 0;

		this->current_generation = 0;
		this->stall_generation = 0;

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
				inputfile.getline(line_buffer, 65536);
			}
			else
			{
				inputfile.getline(line_buffer, 65536);
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

				if (!(number_of_genes == this->Xdescriptions.size()))
				{
					cout << "Number of genes in population file " << filename << " does not match number of genes for this mission script. Perhaps you have loaded the wrong file?" << std::endl;
					throw 0;
				}

				//skip the second header line
				inputfile.getline(line_buffer, 65536);

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
				inputfile.getline(line_buffer, 65536);
				break;
			}
			else
			{
				//grab a line
				inputfile.getline(line_buffer, 65536);

				//split the line by commas
				boost::split(linecell, line_buffer, boost::is_any_of(","));

				if (linecell.size() > 1)
				{
					//the first N entries of the line are genes
					std::vector<int> Xtemp;
					for (size_t gene = 0; gene < number_of_genes; ++gene)
						Xtemp.push_back(atoi(linecell[gene].c_str()));

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
	void outerloop_SGA::writepop(std::string filename)
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
	void outerloop_SGA::read_archive(std::string filename)
	{
		std::ifstream inputfile(filename.c_str());
		int linenumber = 0;
		char line_buffer[65536];
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
			if (peek == '#' || peek == '\r' || peek == '\n' || peek == ',')
			{
				//comment or blank line, do not parse
				inputfile.getline(line_buffer, 65536);
			}
			else
			{
				inputfile.getline(line_buffer, 65536);
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

				if (!(number_of_genes == this->Xdescriptions.size()))
				{
					cout << "Number of genes in archive file " << filename << " does not match number of genes for this mission script. Perhaps you have loaded the wrong file?" << std::endl;
					throw 0;
				}

				//skip the second header line
				inputfile.getline(line_buffer, 65536);

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
				inputfile.getline(line_buffer, 65536);
				break;
			}
			else
			{
				//grab a line
				inputfile.getline(line_buffer, 65536);

				//split the line by commas
				boost::split(linecell, line_buffer, boost::is_any_of(","));

				if (linecell.size() > 1)
				{
					//the first N entries of the line are genes
					std::vector<int> Xtemp;
					for (size_t gene = 0; gene < number_of_genes; ++gene)
						Xtemp.push_back(atoi(linecell[gene].c_str()));

					//add this vector to the population
					this->archive_of_solutions.push_back(EMTG_outerloop_solution(Xtemp, this->number_of_objectives));

					//fill in the description, generation found, and timestamp
					this->archive_of_solutions.back().description = linecell[number_of_genes].c_str();
					this->archive_of_solutions.back().generation_found = atoi(linecell[number_of_genes + 1].c_str());
					this->archive_of_solutions.back().timestamp = atoi(linecell[number_of_genes + 2].c_str());

					//fill in the objective values
					for (size_t objective = 0; objective < this->number_of_objectives; ++objective)
						this->archive_of_solutions.back().fitness_values[objective] = atof(linecell[number_of_genes + 3 + objective].c_str());

					//fill in the inner-loop vector
					for (size_t index = number_of_genes + 3 + this->number_of_objectives; index < linecell.size(); ++index)
					if (!(linecell[index] == ""))
						this->archive_of_solutions.back().Xinner.push_back(atof(linecell[index].c_str()));
				}
			}
		}
	}

	//write an archive file (Jacob)
	void outerloop_SGA::write_archive(std::string filename)
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
	void outerloop_SGA::evolve(const EMTG::missionoptions& options, const boost::ptr_vector<EMTG::Astrodynamics::universe>& Universe)
	{
		std::srand(time(NULL));
		// evolves this_generation population towards Pareto front (main NSGA-II function)

		//Initialize
		//************
		//1. generate random initial parent population (size n)
		//1a. evaluate objectives for initial parent population
		this->evaluatepop(options, Universe); //evaluates this->this_generation vector of solution indiviudals

		//1b. sort the initial population
		std::stringstream popfilestream;
		popfilestream << options.working_directory << "//SGA_initial_population.SGA";
		this->writepop(popfilestream.str());
		popfilestream.clear();
		this->write_archive(options.working_directory + "//SGA_archive.SGA");
		std::sort(this->this_generation.begin(), this->this_generation.end());
		this->children_population.resize(this->popsize);

		this->current_generation = options.outerloop_warmstart;

		// Main loop
		//**************
		while (this->current_generation <= this->genmax)
		{
			//2. write out the current population
			std::stringstream popfilestream;
			popfilestream << options.working_directory << "//SGA_population_gen_" << this->current_generation << ".SGA";
			this->writepop(popfilestream.str());

			//3. increase generation counter
			++this->current_generation;

			//4. generate initial child population from initial parent population (size n)
			for (int individual = 0; individual < this->popsize; ++individual)
			{
				//produce (CR * popsize) individuals via crossover and tournament selection
				if (individual < this->CR * this->popsize)
				{
					EMTG_outerloop_solution Parent1;
					EMTG_outerloop_solution Parent2;
					Parent1 = this->select();
					Parent2 = this->select();
					this->crossover_uniformInt(&Parent1, &Parent2, individual);
				}

				//populate the remaining ((1-CR)*popsize - elitecount) individuals by grabbing from the previous generation with tournament selection
				else if (individual < this->popsize - this->elitecount)
				{
					this->children_population[individual] = this->select();
					this->mutate(&(this->children_population[individual]));
				}

				//add the elite individuals
				else
					this->children_population[individual] = this->archive_of_solutions[this->popsize - 1 - individual];
			}

#ifdef EMTG_MPI
			std::cout << "Processor " << this->MPIWorld->rank() << " has completed crossover and mutation and is ready to evaluate the population in generation " << this->current_generation << std::endl;
#endif
			//5. evaluate objectives for child population
			this->this_generation = this->children_population;
			this->evaluatepop(options, Universe);

			//6. sort the population and the archive
			std::sort(this->this_generation.begin(), this->this_generation.end());
			std::sort(this->archive_of_solutions.begin(), this->archive_of_solutions.end());

			//7. write out the archive
			this->write_archive(options.working_directory + "//SGA_archive.SGA");
		}
	}



	//reset the GA (Jacob)
	void outerloop_SGA::reset()
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
	void outerloop_SGA::calcbounds(const EMTG::missionoptions& options)
	{
		//set up the objectives
		std::vector< std::string > objective_menu_descriptions;
		objective_menu_descriptions.push_back("BOL power at 1 AU (kW)");
		objective_menu_descriptions.push_back("Launch epoch (MJD)");
		objective_menu_descriptions.push_back("Flight time (days)");
		objective_menu_descriptions.push_back("Thruster preference");
		objective_menu_descriptions.push_back("Number of thrusters");
		objective_menu_descriptions.push_back("Launch vehicle preference");
		objective_menu_descriptions.push_back("Delivered mass to final target (kg)");
		objective_menu_descriptions.push_back("Final journey mass increment (for maximizing sample return)");
		objective_menu_descriptions.push_back("First journey departure C3 (km^2/s^2)");
		objective_menu_descriptions.push_back("Final journey arrival C3 (km^2/s^2)");
		objective_menu_descriptions.push_back("Total delta-v (km/s)");
		objective_menu_descriptions.push_back("Inner-loop objective (whatever it was)");

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
			if (options.outerloop_vary_journey_flyby_sequence[j] && options.outerloop_journey_flyby_sequence_choices[j].size() > 0)
			{
				for (size_t p = 0; p < options.outerloop_journey_maximum_number_of_flybys[j]; ++p)
				{
					this->Xlowerbounds.push_back(0);
					this->Xupperbounds.push_back((options.outerloop_journey_flyby_sequence_choices[j].size() - 1) * 2);
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

