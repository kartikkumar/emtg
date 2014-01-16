//header file for integer GA
//for use with EMTG v0.2 and later
//Jacob Englander 11/20/2010
//Modified to understand Universe flyby menus 11/16/2012

#ifndef _GA
#define _GA

#include "universe.h"

#include <vector>
#include <string>
#include "missionoptions.h"
#include <ctime>

#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"

using namespace std;

class integerGA
{
private:
	//variables
	vector< vector<int> > population, old_population;
	vector<int> upper_bounds;
	vector<int> lower_bounds;
	vector<double> scores;
	vector<string> descriptions, old_descriptions;
	vector< vector<int> > solutions_database;
	vector<int> solutions_database_indices;
	vector<int> solutions_database_generations;
	vector<double> scores_database;
	vector<time_t> times_generation, old_times_generation;
	vector<time_t> times_database;
	vector<string> descriptions_database;
	vector<int> indiceslist; //always goes 0:(popsize-1)
	double CR; //crossover rate
	double mu; //mutation rate
	int tournamentsize; //tournament size for selection
	int stallmax; //maximum number of stall generations
	int popsize; //population size
	int genmax; //maximum number of generations
	double tolfit; //fitness tolerance
	int encodelength; //maximum value for any element
	int chromosomelength; //length of a chromosome
	int elitecount; //how many elite individuals to retain
	int phase_encode_length; //a phase is transcribed by two values if the phase type is free to be chosen by the optimizer, one otherwise

	time_t tstart, tfinish;

	int generation, stallgen; //counters to track during the optimization

	vector<int> parentpool;
	vector<double> parentscores;
	


	boost::mt19937 RNG;
	boost::uniform_int<> IntegerDistribution;
	boost::uniform_real<> DoubleDistribution;

	//methods
	int select(int* parent);
	int mutate(int index);
	int crossover(int parent1, int parent2, int index);
	int evaluate(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse);
	int calcbounds(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse);

public:
	//variables
	double optfit;

	//constructor
	integerGA(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse, boost::mt19937 RNG_in);

	//methods
	int generatepop();
	int readpop();
	int evolve(EMTG::missionoptions* options, boost::ptr_vector<EMTG::Astrodynamics::universe>& TheUniverse);
	int writepop(EMTG::missionoptions* options);
	int reset();
};

#endif //_GA
