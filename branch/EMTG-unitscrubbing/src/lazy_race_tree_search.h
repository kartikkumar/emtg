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


#ifdef EMTG_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#endif

using namespace std;
using namespace boost::filesystem;
using namespace boost::gregorian;
using namespace boost::posix_time;

namespace EMTG{

	void lazy_race_tree_search(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in);
}