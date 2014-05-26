
#include "missionoptions.h"
#include "mission.h"
#include "outerloop_NSGAII.h"
#include "outerloop_SGA.h"
#include <cstdlib>

#include "universe.h"
#include "body.h"

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"

#include "SpiceUsr.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>



#ifdef EMTG_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#endif

using namespace std;
using namespace boost::filesystem;
using namespace boost::gregorian;
using namespace boost::posix_time;

#ifndef _LAZY_RACE_TREE_SEARCH
#define _LAZY_RACE_TREE_SEARCH

namespace EMTG{
#ifdef EMTG_MPI
	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in, std::vector <int> & asteroid_list, std::vector <int> & best_sequence, std::string & branch_directory, std::string & tree_summary_file_location, boost::mpi::environment MPIenvironment, boost::mpi::communicator world);
#else
	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in, std::vector <int> & asteroid_list, std::vector <int> & best_sequence, std::string & branch_directory, std::string & tree_summary_file_location);
#endif
	
	std::vector<int> filter_asteroid_list(int const & current_asteroid, double const & epoch, std::vector<int> & asteroid_list, boost::ptr_vector<Astrodynamics::universe> & myUniverse, missionoptions * options);
	
	
	struct pair_min {

		std::pair<int, double> operator() (const std::pair<int, double> & left, std::pair<int, double> right) { return ((left.second < right.second) ? (left) : (right)); }

	};

	
	
	//ugh; why are we header declaring functions that aren't in the matching .cpp??? - Alex
	void load_asteroid_list(std::string & asteroid_filename, std::vector <int> & asteroid_list);
	void write_branch_summary(EMTG::mission & branch_mission, EMTG::missionoptions & branch_options, std::string & tree_summary_file_location, int & tree_level, const int & branch);
	
	
}

#endif