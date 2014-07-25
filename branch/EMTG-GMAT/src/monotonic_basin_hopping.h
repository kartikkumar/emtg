//header file for Monotonic Basin Hopping
//for EMTG version 8
//Jacob Englander 7-27-2012

#include "problem.h"

#include "boost/random/uniform_int.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/mersenne_twister.hpp"

#include "snopt.hh"
#include "snoptProblemExtension.h"
#include "snfilewrapper.hh"



#ifdef _use_WORHP
#include "EMTG_WORHP_interface.h"
#endif


#ifndef _EMTG_SOLVERS
#define _EMTG_SOLVERS

namespace EMTG { namespace Solvers {

class MBH
{
public:
	//constructor
	MBH();
	MBH(EMTG::problem* Problem_input);

	//destructor
	virtual ~MBH();

	//methods
	int initialize(EMTG::problem* Problem_input);
	int reset_point();
	int seed(vector<double> seed_vector);
	int hop();
	int time_hop();
	int slide();
	int run();
	double check_feasibility();
	int print_archive_header(string filename);
	int print_archive_line(string filename, int linenumber);
	int write_sparsity(string filename);
	int read_sparsity(string filename);
	
	
	//fields

	//pointer to problem object
	EMTG::problem* Problem;

	//variable to track the amount of time for a given SNOPT run
	time_t SNOPT_start_time;

	//dummy real array
	doublereal* DummyReal;

	//archive of feasible decision vectors
	vector<double> Xtrial_scaled;
	double ftrial;
	vector<double> Xlocalbest_scaled;
	double flocalbest;
	vector<double> Xcurrent_scaled;
	double fcurrent;
	vector<double> Xbest_scaled;
	double fbest;
	vector<double> Xbest;
	vector< vector<double> > archive;
	vector<double> archive_scores;
	vector<time_t> archive_timestamps;
	vector<int> archive_step_count;
	vector<int> archive_reset_count;
	double most_feasible;
	vector<double> closest_to_feasible_solution;
	double fcrit;

	//helper arrays
	vector<int> time_variable_indices;
	vector<int> significant_variable_indices;

	//counters
	int number_of_solutions; //how many feasible solutions found so far
	int number_of_improvements; //how many times have we improved within a basin
	int number_of_resets;
	int number_of_failures_since_last_improvement;
	int Jacobian_offset;

	//track the step size
	double step_size;

	//track whether or not a Jacobian sparsity pattern has yet been identified
	bool computed_Jacobian;
	bool jacfullrankflag;

	//track whether or not the sparsity file and XFfile have been printed
	bool printed_sparsity;

	//random number generator
	boost::mt19937 RNG;
	boost::uniform_real<> DoubleDistribution;

	//pointer to SNOPT object
	snoptProblemExtension* SNOPTproblem;
	
	//other fields for SNOPT
	integer neF;

	integer lenA;

	integer *iAfun;
	integer *jAvar;
	doublereal *A;

	integer lenG;
	integer *iGfun;
	integer *jGvar;

	doublereal *x;
	doublereal *xlow;
	doublereal *xupp;
	doublereal *xmul;
	integer    *xstate;

	doublereal *F;
	doublereal *Flow;
	doublereal *Fupp;
	doublereal *Fmul;
	integer    *Fstate;

	integer nxnames;
	integer nFnames;
	char *xnames;
	char *Fnames;

	integer    ObjRow;
	doublereal ObjAdd;
	
	//pointer to WORHP object

#ifdef _use_WORHP
	EMTG_WORHP_interface* WORHP_interface;
#endif
};


}}//close namespace


#endif //_EMTG_SOLVERS