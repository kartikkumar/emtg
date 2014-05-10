//header file for EMTG outer-loop solution class
//for use with NSGAII and SGA
//M. Vavrina and J. Englander

#include "missionoptions.h"
#include "universe.h"

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

#ifndef EMTG_OUTERLOOP_SOLUTION
#define EMTG_OUTERLOOP_SOLUTION

namespace GeneticAlgorithm
{
	class EMTG_outerloop_solution
	{
	public:
		//constructor
		EMTG_outerloop_solution() {};
		EMTG_outerloop_solution(const std::vector<int>& X_in, const int& number_of_objectives);
		//destructor
		~EMTG_outerloop_solution() {};

		//methods
		void evaluate();
		EMTG::missionoptions parse_outer_loop_decision_vector();
		void set_data_pointers(void* optionspointer, void* Universepointer);
		void write_csv_line(std::string outputfile);

		//fields
		std::vector<int> X;
		std::vector<double> Xinner;
		time_t timestamp;
		int generation_found;
		std::vector<double> fitness_values;
		std::string description;
		std::vector< int > dominated_by;	// vector of other EMTG_outerloop_solution objects that the solution is dominated by
		std::vector< int > dominates;		// vector of other EMTG_outerloop_solution objects that the solution dominates;
		double crowding_distance;								// relative measure of distance to others along non-dominated front
		int pareto_rank;										// level of non-dominated front that the solution belongs to (the lower the closer to the pareto front)
		int ndom;												// dominated-by counter for non-dpminated sorting
		double innerloop_fitness;

		//members

		//comparators
		bool operator>(EMTG_outerloop_solution OtherSolution);
		bool operator<(EMTG_outerloop_solution OtherSolution);
		bool operator>=(EMTG_outerloop_solution OtherSolution);
		bool operator<=(EMTG_outerloop_solution OtherSolution);
		bool operator==(EMTG_outerloop_solution OtherSolution);
		bool compare_objective_greaterthan(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit);
		bool compare_objective_greaterthanorequalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit);
		bool compare_objective_lessthan(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit);
		bool compare_objective_lessthanorequalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit);
		bool compare_objective_equalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit);
		bool compare_description(EMTG_outerloop_solution& OtherSolution);

	private:
		//serialization code, if applicable
#ifdef EMTG_MPI
		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & this->X;
			ar & this->Xinner;
			ar & this->timestamp;
			ar & this->generation_found;
			ar & this->fitness_values;
			ar & this->description;
			ar & this->crowding_distance;
			ar & this->pareto_rank;
			ar & this->ndom;
			ar & this->dominated_by;
			ar & this->dominates;
			ar & this->innerloop_fitness;
		}
#endif

		EMTG::missionoptions* options_base;
		boost::ptr_vector<EMTG::Astrodynamics::universe>* Universe;

	};
}// close namespace

#endif //EMTG_OUTERLOOP_SOLUTION