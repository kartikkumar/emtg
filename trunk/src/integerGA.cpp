//integer GA
//for use with EMTG v0.2 and later
//Jacob Englander 11/20/2010

#include "missionoptions.h"
#include "integerGA.h"
#include "quicksort.h"
#include "mission.h"

#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"

#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;


//constructor
integerGA::integerGA(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse, boost::mt19937 RNG_in)
{
	mu = options->outerloop_mu;
	CR = options->outerloop_CR;
	genmax = options->outerloop_genmax;
	popsize = options->outerloop_popsize;
	tolfit = options->outerloop_tolfit;
	tournamentsize = options->outerloop_tournamentsize;
	
	elitecount = options->outerloop_elitecount;
	stallmax = options->outerloop_stallmax;
		
	chromosomelength = calcbounds(options, TheUniverse);

	phase_encode_length = (options->mission_type > 2 ? 2 : 1);

	//size the population
	vector<int> individual(chromosomelength,0);

	for (int k = 0; k < popsize; ++k)
	{
		population.push_back(individual);
		indiceslist.push_back(k);
	}

	scores.resize(popsize);
	times_generation.resize(popsize);
	descriptions.resize(popsize);
	parentpool.resize(tournamentsize);
	parentscores.resize(tournamentsize);

	RNG = RNG_in;
	DoubleDistribution = boost::uniform_real<>(0.0, 1.0);


}

//reset function
int integerGA::reset()
{
	population.clear();
	indiceslist.clear();
	scores.clear();
	descriptions.clear();
	parentpool.clear();
	parentscores.clear();

	//size the population
	vector<int> individual(chromosomelength,0);

	for (int k=0;k<popsize;++k)
	{
		population.push_back(individual);
		indiceslist.push_back(k);
	}

	scores.resize(popsize);
	descriptions.resize(popsize);
	parentpool.resize(tournamentsize);
	parentscores.resize(tournamentsize);

	return 0;
}

//bounds generation function
int integerGA::calcbounds(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse)
{
	upper_bounds.clear();
	lower_bounds.clear();

	for (int j = 0; j < options->number_of_journeys; ++j)
	{
		for (int p = 0; p < options->max_phases_per_journey; ++p)
		{
			lower_bounds.push_back(0);
			upper_bounds.push_back(TheUniverse[j].size_of_flyby_menu);

			if (options->mission_type == 4 || options->mission_type == 5 || options->mission_type == 6) //optimizer chooses between two choices offered
			{
				lower_bounds.push_back(0);
				upper_bounds.push_back(1);
			}
			else if (options->mission_type == 7) //optimizer chooses MGA, MGA-DSM, or MGA-LT
			{
				lower_bounds.push_back(0);
				upper_bounds.push_back(2);
			}
		}
	}

	return upper_bounds.size();
}

//population generation function
int integerGA::generatepop()
{
	optfit = 1.0e10;

	//place a random integer in [0,encodelength) in each bin
	for (int i=0;i<popsize;++i)
	{
		for (int j=0;j<chromosomelength;++j)
		{
			IntegerDistribution = boost::uniform_int<>(0, upper_bounds[j] - 1);
			population[i][j] = IntegerDistribution(RNG);
		}
	}

	return 0;
}

//function to read in a population and solutions database for a "warm start"
int integerGA::readpop()
{
	//first read in the current population
	string tempstring;
	int tempint;
	double tempdouble;
	ifstream inputfile("population.txt", ios::in);
	for (int i=0;i<popsize;++i)
	{
		inputfile >> tempstring;

		for (int j=0;j<chromosomelength;++j)
			inputfile >> population[i][j];

		inputfile >> descriptions[i];

		inputfile >> scores[i];
	}
	inputfile.close();

	//then read in and parse the solutions database
	inputfile.open("solutions.txt", ios::in);
	int count = 0;
	while (!inputfile.eof())
	{
		vector<int> tempintvec(chromosomelength);
		inputfile >> tempstring;
		solutions_database_indices.push_back(count);
		for (int i=0; i<chromosomelength; ++i)
			inputfile >> tempintvec[i];

		solutions_database.push_back(tempintvec);

		inputfile >> tempstring;
		descriptions_database.push_back(tempstring);
		
		inputfile >> tempint;
		solutions_database_generations.push_back(tempint);

		inputfile >> tempdouble;
		scores_database.push_back(tempdouble);

		inputfile >> tempint;
		times_database.push_back(tempint);

		++count;
	}

	return 0;
}

//selection function
int integerGA::select(int* parent)
{
	//first fill the parent pool
	int flag;

	IntegerDistribution = boost::uniform_int<>(0, popsize-1);

	parentpool[0] = IntegerDistribution(RNG);

	for (int k = 1; k < tournamentsize; ++k)
	{
		do
		{
			flag = 0;
			parentpool[k] = IntegerDistribution(RNG);

			for (int q = 0; q < k-1; ++q) //check to see if we're in the parent pool
			{
				if (parentpool[k] == parentpool[q]) //if this index is already in the parent pool, turn on the flag
					flag = 1;
			}
		} while (flag);
	}

	//record the scores of the parents
	for (int k = 0; k < tournamentsize; ++k)
		parentscores[k] = scores[parentpool[k]];

	//then sort the parent pool
	EMTG::quicksort(parentpool, parentscores, 0, tournamentsize-1);
	*parent = parentpool[0];

	return 0;
}

//mutation function
int integerGA::mutate(int index)
{
	for (int k = 0; k < chromosomelength; ++k)
	{
		if (DoubleDistribution(RNG) < mu)
		{
			IntegerDistribution = boost::uniform_int<>(0, upper_bounds[k] - 1);
			population[index][k] = IntegerDistribution(RNG);
		}
	}

	return 0;
}

//crossover function
int integerGA::crossover(int parent1, int parent2, int index)
{
	//choose two crossover points
	IntegerDistribution = boost::uniform_int<>(0, ((chromosomelength/phase_encode_length - 1) * phase_encode_length));

	int crosspoint1 = IntegerDistribution(RNG); 
	int crosspoint2 = IntegerDistribution(RNG); 

	if (crosspoint1 > crosspoint2)
	{
		int temp = crosspoint2;
		crosspoint2 = crosspoint1;
		crosspoint1 = temp;
	}

	for (int k = 0; k < crosspoint1; ++k)
		population[index][k] = old_population[parent1][k];
	for (int k = crosspoint1; k<crosspoint2; ++k)
		population[index][k] = old_population[parent2][k];
	for (int k = crosspoint2; k<chromosomelength; ++k)
		population[index][k] = old_population[parent1][k];

	return 0;
}

//evaluate function
int integerGA::evaluate(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse)
{
	if (options->outerloop_useparallel)
	{

	}
	else //run the outer-loop in serial
	{
		vector<double> augmented_scores_database = scores_database;
		vector<string> augmented_descriptions_database = descriptions_database;
		

		for (int k=0; k < popsize; ++k)
		{
			time_t tstart = time(NULL);

			vector<int> Xouter = population[k];

			EMTG::mission TrialMission(&Xouter[0], options, TheUniverse, 0, 0);

			bool optimize_flag = true;
			for (size_t i = 0; i < augmented_scores_database.size(); ++i)
			{
				if (augmented_descriptions_database[i] == TrialMission.options.description) //this sequence has already been evaluated
				{
					cout << "Mission " << TrialMission.options.description << " has already been evaluated with fitness " << augmented_scores_database[i] << endl;
					scores[k] = augmented_scores_database[i];
					descriptions[k] = augmented_descriptions_database[i];
					times_generation[k] = (time(NULL) - tstart);
					optimize_flag = false;
					break;
				}
			}
		
			if (optimize_flag)
			{
				cout << "Optimizing mission: " << TrialMission.options.description << endl;
				TrialMission.output_mission_tree(options->working_directory + "//" + TrialMission.options.mission_name + "_" + TrialMission.options.description + "_missiontree.emtgtree");
				TrialMission.optimize();
				scores[k] = TrialMission.F[0];
				descriptions[k] = TrialMission.options.description;
				times_generation[k] = (time(NULL) - tstart);
				augmented_scores_database.push_back(scores[k]);
				augmented_descriptions_database.push_back(descriptions[k]);
			}
		}
	}

	return 0;
}

//evolve function
int integerGA::evolve(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse)
{
	tstart = time(NULL);

	//first, generate an initial population
	if (options->outerloop_warmstart)
	{
		readpop();
		
		//set generations to zero
		generation = options->outerloop_warmstart;
	}
	else
	{
		generatepop();
		
		//set generations to zero
		generation = 0;
	}

	//set stall generations to zero
	stallgen = 0;

	//then evaluate the initial population
	evaluate(options, TheUniverse);

	//sort the initial population
	old_population = population;
	old_descriptions = descriptions;
	vector<int> tempindices = indiceslist;
	EMTG::quicksort(tempindices, scores, 0, popsize-1);
	for (int k=0;k<popsize;++k)
	{
		population[k] = old_population[tempindices[k]];
		descriptions[k] = old_descriptions[tempindices[k]];
	}

	//print the population, sort and print the database of solutions found so far (remove duplicates)
	writepop(options);
	//print to the convergence history
	ofstream outputfile;
	outputfile.open (options->convergence_file.c_str(), ios::out | ios::trunc);
	outputfile << "Convergence history" << endl;
	outputfile << endl;

	//now iterate
	int parent1, parent2;

	while (stallgen < stallmax && generation < genmax)
	{	
		
		//increment the generation count
		++generation;

		//save the old population
		old_population = population;

		//temporary code, zero out the new population
		for (int i=0;i<popsize;++i)
		{
			for (int j=0;j<chromosomelength;++j)
				population[i][j] = 0;
		}


		for (int k=0; k<popsize; ++k)
		{
			//produce (CR*popsize) individuals via crossover and tournament selection
			if (k < CR*popsize)
			{
				select(&parent1);
				select(&parent2);
				crossover(parent1, parent2, k);
			}

			//populate the remaining ((1-CR)*popsize - elitecount) individuals by grabbing from the previous generation with tournament selection
			else if (k<popsize - elitecount)
			{
				select(&parent1);
				population[k] = old_population[parent1];
				mutate(parent1);
			}

			//add the elite individuals
			else
				population[k] = old_population[popsize - 1 - k];
		}

		//evaluate everybody
		evaluate(options, TheUniverse);

		//sort the population
		old_population = population;
		old_descriptions = descriptions;
		old_times_generation = times_generation;

		tempindices = indiceslist;
		EMTG::quicksort(tempindices, scores, 0, popsize-1);
		for (int k=0;k<popsize;++k)
		{
			population[k] = old_population[tempindices[k]];
			descriptions[k] = old_descriptions[tempindices[k]];
			times_generation[k] = old_times_generation[tempindices[k]];
		}
		//check convergence
		if (optfit - scores[0] <= tolfit) //this is a stall generation
			++stallgen;
		else
		{
			optfit = scores[0];
			stallgen = 0;
		}

		//print the population, sort and print the database of solutions found so far (remove duplicates)
		writepop(options);

		//write generation count and best solution to the screen
		outputfile << "Generation " << generation <<  ", t="<<time(NULL)-tstart<<" - best sequence is " <<descriptions_database[0]<< " with fitness " << scores_database[0] << endl;
	}

	tfinish = time(NULL);

	return 0;
}

//sort and write population to the disk
int integerGA::writepop(EMTG::missionoptions* options)
{
	//first write out the current population
	ofstream outputfile;

	//generate the population file
	outputfile.open (options->population_file.c_str(), ios::out | ios::app);
	outputfile << "Population at generation "<< generation << endl;
	outputfile << endl;

	for (int i = 0; i < popsize; ++i)
	{
		outputfile << "[" << i+1 << "]";
		for (int j = 0; j < chromosomelength; ++j)
		{
			outputfile.width(3); outputfile << population[i][j];
		}
		outputfile.width(chromosomelength + 2); outputfile << descriptions[i];
		outputfile.width(3); outputfile <<" "<< scores[i] << endl;
	}

	outputfile.close();

	//now update the solutions database
	int flag;
	for (int i = 0; i < popsize; ++i)
	{
		flag = 1;
		for (size_t j = 0; j<(solutions_database_indices.size()); ++j)
		{
			if (descriptions[i] == descriptions_database[j])
			{
				if (scores[i] < scores_database[j])
					scores_database[j] = scores[i];
				flag = 0;
				break;
			}
		}

		if (flag)
		{
			solutions_database_indices.push_back(solutions_database_indices.size());
			solutions_database.push_back(population[i]);
			scores_database.push_back(scores[i]);
			times_database.push_back(times_generation[i]);
			descriptions_database.push_back(descriptions[i]);
			solutions_database_generations.push_back(generation);
		}
	}
	
	//sort the solutions database
	vector< vector<int> > old_database = solutions_database;
	old_descriptions = descriptions_database;
	vector<time_t> old_times_database = times_database;
	vector<int> tempindices = solutions_database_indices;
	vector<int> old_solutions_database_generations = solutions_database_generations;
	EMTG::quicksort(tempindices, scores_database, 0, solutions_database_indices.size()-1);
	
	for (size_t k = 0; k < solutions_database.size() ; ++k)
	{
		solutions_database[k] = old_database[tempindices[k]];
		solutions_database_generations[k] = old_solutions_database_generations[tempindices[k]];
		descriptions_database[k] = old_descriptions[tempindices[k]];
		times_database[k] = old_times_database[tempindices[k]];
	}
	
	//write out the solutions database
	outputfile.open (options->solutions_file.c_str(), ios::out | ios::trunc);
	outputfile << "Solutions database as of generation "<< generation << endl;
	outputfile << endl;

	for (size_t i=0; i < solutions_database_indices.size(); ++i)
	{
		outputfile.width(6);
		stringstream tempstream;
		tempstream << "[" << i+1 << "]";
		outputfile << tempstream.str();
		for (int j=0;j<chromosomelength;++j)
		{
			outputfile.width(3); outputfile << solutions_database[i][j];
		}
		outputfile.width(80); outputfile << descriptions_database[i];
		outputfile.width(5); outputfile << solutions_database_generations[i];
		outputfile.width(14); outputfile << scores_database[i];
		outputfile.width(8); outputfile << times_database[i] << endl;
	}

	outputfile.close();
		
	return 0;
}
