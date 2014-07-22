//header file for outer-loop single objective genetic algorithm (SGA) class
//can solve full systems optimization problem for single objective
//integer encoding
//J. Englander 5-10-2014
//based heavily on outer-loop NSGAII by M. Vavrina and J. Englander

#include "missionoptions.h"
#include "universe.h"
#include "EMTG_outerloop_solution.h"

#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/ptr_container/ptr_vector.hpp"

#include <vector>
#include <string>
#include <ctime>

#ifdef EMTG_MPI
#include "boost/serialization/string.hpp"
#include "boost/serialization/vector.hpp"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#endif

#ifndef OUTERLOOP_SGA
#define OUTERLOOP_SGA

namespace GeneticAlgorithm
{
	class outerloop_SGA
	{
	private:
		//methods
		EMTG_outerloop_solution select();
		void crossover_uniformInt(EMTG_outerloop_solution* Parent1, EMTG_outerloop_solution* Parent2, int index_in_child_population);
		void mutate(EMTG_outerloop_solution* Individual);

		//distributions
		std::vector< boost::uniform_int<> > random_integer;

		//fields
		boost::mt19937 RNG;
		boost::uniform_int<> IntegerDistribution;
		boost::uniform_real<> DoubleDistribution;
		double CR; //crossover rate
		double mu; //mutation rate
		int elitecount;
		int tournamentsize; //tournament size for selection
		int stallmax; //maximum number of stall generations
		int genmax; //maximum number of generations
		int popsize; //number of individuals in the population
		int current_generation;
		int stall_generation;
		int number_of_objectives;
		std::vector<std::string> objective_descriptions;
		std::vector<double> tolfit; //tolerance for each fitness function
		std::vector< EMTG_outerloop_solution > this_generation, last_generation;
		std::vector< EMTG_outerloop_solution > parent_population, children_population;
		std::vector< EMTG_outerloop_solution > archive_of_solutions;
		std::vector<int> Xupperbounds, Xlowerbounds;
		std::vector<std::string> Xdescriptions;

		time_t tstart, tfinish;

		//communication code - pointers to World and Environment
#ifdef EMTG_MPI
		boost::mpi::environment* MPIEnvironment;
		boost::mpi::communicator* MPIWorld;
#endif

	public:
		//constructor
		outerloop_SGA(const EMTG::missionoptions& options); //seedless
		outerloop_SGA(int boost_random_seed, const EMTG::missionoptions& options);
#ifdef EMTG_MPI
		outerloop_SGA(const EMTG::missionoptions& options, boost::mpi::environment* MPIEnvironment, boost::mpi::communicator* MPIWorld);
		outerloop_SGA(int boost_random_seed, const EMTG::missionoptions& options, boost::mpi::environment* MPIEnvironment, boost::mpi::communicator* MPIWorld);
#endif

		//methods
		void evaluatepop(const EMTG::missionoptions& options, const boost::ptr_vector<EMTG::Astrodynamics::universe>& Universe);
		void generatepop();
		void readpop(std::string filename);
		void writepop(std::string filename);
		void read_archive(std::string filename);
		void write_archive(std::string filename);
		void evolve(const EMTG::missionoptions& options, const boost::ptr_vector<EMTG::Astrodynamics::universe>& Universe);
		void reset();
		void calcbounds(const EMTG::missionoptions& options);
		void startclock() { this->tstart = time(NULL); }

		void supply_pop(std::vector< EMTG_outerloop_solution > const & input_population);
		void supply_archive(std::vector< EMTG_outerloop_solution > const & input_archive);
		const std::vector< EMTG_outerloop_solution >& get_population();
		const std::vector< EMTG_outerloop_solution >& get_archive();
		void set_CR(double const & CR_in);
		void set_mutationrate(double const & mu_in);
		void set_elitecount(unsigned int const & elitecount_in);
		void set_tournament_size(unsigned int const & tournamentsize_in);
		void set_stall_generations(unsigned int const & stall_generations_in);
		void set_populationsize(unsigned int const & populationsize_in);
		void set_max_generations(unsigned int const & max_generations_in);
		void set_tolfit(std::vector<double> const & tolfit_in);
		void set_current_generation(const int& gen_in);
		void set_objectives(const int& number_of_objectives_in, const std::vector<std::string>& objective_descriptions_in);
	};
} //close namespace

#endif //OUTERLOOP_SGA