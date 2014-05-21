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

#ifdef _STONEAGECplusplus
#include <execinfo.h>
#include <signal.h>
void handler(int sig) {
	void *array[10];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, 10);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(1);
}
#endif

namespace EMTG{

	void LRTS(missionoptions * options, boost::ptr_vector<Astrodynamics::universe> & TheUniverse_in);
}