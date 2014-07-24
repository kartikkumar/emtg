/*
 * problem.h
 *
 *  Created on: Jul 17, 2012
 *      Author: Jacob
 */

#include <string>
#include <vector>

#include "missionoptions.h"


#ifndef PROBLEM_H_
#define PROBLEM_H_

namespace EMTG {

class problem {
public:
	//constructor
	problem();

	//destructor
	virtual ~problem();

	//methods

	//optimize function
	//return 0 for success, 1 for failure
	int optimize();

	//bounds calculation function
	//return 0 for success, 1 for failure
	virtual int calcbounds() = 0;

	//functions to scale and unscale decision vectors
	int unscale(double* Xscaled);
	int scale(double* Xscaled);

	//outer-loop parse function
	//return 0 for success, 1 for failure
	virtual int parse_outer_loop(int* Xouter, int n_outer_loop);

	//function to output X and F bounds, descriptions
	int output_problem_bounds_and_descriptions();
	int output_problem_bounds_and_descriptions(string filestring);
	

	//function to output the Jacobian sparsity information
	virtual int output_Jacobian_sparsity_information(string filestring);

	//function to check the derivatives via central differencing
	int check_and_print_derivatives();
	int check_and_print_derivatives(string filestring);
	

	//virtual function templates
	virtual int evaluate(double* X, double* F, double* G, int needG, const vector<int>& iGfun, const vector<int>& jGvar) = 0;
	virtual int output() = 0;
	virtual vector<double> create_initial_guess(vector<double> XFBLT, vector<string>& NewXDescriptions) = 0;
	virtual void interpolate(int* Xouter, const vector<double>& initialguess) = 0;
	virtual void convert_cartesian_solution_to_polar(const vector<double>& initialguess) = 0;
	virtual void convert_polar_solution_to_cartesian(const vector<double>& initialguess) = 0;
	//performance characteristics function
	//used to extract various pieces of mission data for a multi-objective GA
	virtual void extract_objective_function_values(std::vector<double>& objective_functions) {};

	//fields

	//identifying information
	int thread_ID; //what thread are we assigned to
	int problem_ID; //what problem are we

	//encapsulated options structure
	missionoptions options;

	//solver parameters
	vector<double> X0; //initial guess, may be supplied or generated at random
	vector<double> X; //current decision vector
	vector<double> Xopt; //optimal decision vector
	vector<double> F; //constraint vector
	vector<double> G; //nonlinear Jacobian vector
	vector<double> A; //linear Jacobian vector
	double best_cost;
	int total_number_of_constraints; //total number of nonlinear constraints
	int total_number_of_NLP_parameters; //total number of NLP parameters
	vector<double> Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds;
	vector<string> Xdescriptions, Fdescriptions;
	vector<int> iAfun, jAvar, iGfun, jGvar;
	vector<string> Adescriptions, Gdescriptions;
	int number_of_solutions;

	//vector of synodic periods for any periodic variables
	vector<double> synodic_periods;
};

} /* namespace EMTG */
#endif /* PROBLEM_H_ */
