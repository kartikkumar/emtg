/*
 * mission.h
 *
 *  Created on: Jul 17, 2012
 *      Author: Jacob
 */

#ifndef MISSION_H_
#define MISSION_H_

#include <vector>

#include "problem.h"
#include "journey.h"
#include "file_utilities.h"

#include "boost/ptr_container/ptr_vector.hpp"

namespace EMTG {

class mission: public EMTG::problem 
{
public:
	//constructor
	mission();
    mission(int* Xouter, const missionoptions& options_in, const boost::ptr_vector<Astrodynamics::universe>& TheUniverse_in);
    mission(const missionoptions& options_in, const boost::ptr_vector<Astrodynamics::universe>& TheUniverse_in);

	//destructor
	virtual ~mission();

	//methods
	//evaluate function
	//return 0 if successful, 1 if failure
	virtual int evaluate(double* X, 
                        double* F,
                        double* G, 
                        int needG, 
                        const vector<int>& iGfun,
                        const vector<int>& jGvar);

	//output function
	//return 0 if successful, 1 if failure
	virtual void output();

    //method to output a forward-integrated ephemeris
    virtual void write_ephemeris_file();

	//performance characteristics function
	//used to extract various pieces of mission data for a multi-objective GA
	void extract_objective_function_values(std::vector<double>& objective_functions);

	//output mission structure
	virtual void output_mission_tree(string filename);

	//bounds calculation function
	//return 0 for success, 1 for failure
	virtual void calcbounds();

	//outer-loop parse function
	//return 0 for success, 1 for failure
	virtual int parse_outer_loop(int* Xouter);

	//function to create an initial guess for another mission type
	virtual void create_initial_guess(const int& desired_mission_type, const bool& VSI);

	//function to interpolate an initial guess
	virtual void interpolate(int* Xouter, const vector<double>& initialguess);

	//functions to convert initial guess between cartesian and polar coordinates
	virtual void convert_cartesian_solution_to_polar(const vector<double>& initialguess);
	virtual void convert_polar_solution_to_cartesian(const vector<double>& initialguess);


	//function to find constraint/objective function dependencies due to an spiral anywhere in the mission
	void find_dependencies_due_to_escape_spiral(vector<double>* Xupperbounds,
												vector<double>* Xlowerbounds,
												vector<double>* Fupperbounds,
												vector<double>* Flowerbounds,
												vector<string>* Xdescriptions,
												vector<string>* Fdescriptions,
												vector<int>* iAfun,
												vector<int>* jAvar,
												vector<int>* iGfun,
												vector<int>* jGvar,
												vector<string>* Adescriptions,
												vector<string>* Gdescriptions,
												missionoptions* options,
												const int& Findex);

	void find_dependencies_due_to_capture_spiral(vector<double>* Xupperbounds,
												vector<double>* Xlowerbounds,
												vector<double>* Fupperbounds,
												vector<double>* Flowerbounds,
												vector<string>* Xdescriptions,
												vector<string>* Fdescriptions,
												vector<int>* iAfun,
												vector<int>* jAvar,
												vector<int>* iGfun,
												vector<int>* jGvar,
												vector<string>* Adescriptions,
												vector<string>* Gdescriptions,
												missionoptions* options,
												const int& Findex);

	void find_dependencies_due_to_capture_spiral_in_final_journey(vector<double>* Xupperbounds,
																vector<double>* Xlowerbounds,
																vector<double>* Fupperbounds,
																vector<double>* Flowerbounds,
																vector<string>* Xdescriptions,
																vector<string>* Fdescriptions,
																vector<int>* iAfun,
																vector<int>* jAvar,
																vector<int>* iGfun,
																vector<int>* jGvar,
																vector<string>* Adescriptions,
																vector<string>* Gdescriptions,
																missionoptions* options,
																const int& Findex);

	//vector of journeys
	boost::ptr_vector<journey> journeys;

	//vector of universes
	boost::ptr_vector<Astrodynamics::universe> TheUniverse;

	//fields
	int number_of_journeys;
	double current_deltaV; //value to hold the current deltaV as we progress through the mission
	double current_epoch; //value to hold the current epoch as we progress through the mission
	double current_state[7]; //array to hold the current spacecraft state as we progress through the mission
	double dry_mass; //in kg
	double total_propellant_mass; //in kg

	//derivative information
	vector<int> objective_function_G_indices;
	vector<int> timeconstraints_G_indices;
	vector<double> timeconstraints_X_scale_ranges;
	int derivative_of_flight_time_with_respect_to_launch_date_G_index;
	vector<int> derivative_of_flight_time_with_respect_to_journey_initial_mass_increment_ratios_for_spirals;
	vector<int> dry_mass_constraint_G_indices;
	vector<int> dry_mass_constraint_X_indices;
	vector<double> dry_mass_constraint_X_ranges;
    vector<int> final_mass_constraint_G_indices;
    vector<int> final_mass_constraint_X_indices;
    vector<double> final_mass_constraint_X_ranges;
	vector<int> propellant_mass_constraint_G_indices;
	vector<int> propellant_mass_constraint_X_indices;
	vector<double> propellant_mass_constraint_X_ranges;
	vector<int> objectivefunction_X_indices;
	vector<int> objectivefunction_G_indices;
	vector<double> objectivefunction_X_scale_ranges;

	//time information
	double max_TU; //largest time unit, for constraint scaling
};

} /* namespace EMTG */
#endif /* MISSION_H_ */
