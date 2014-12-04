/*
 * journey.h
 *
 *  Created on: Jul 17, 2012
 *      Author: Jacob
 */
#include <vector>
#include <memory>

#include "missionoptions.h"
#include "phase.h"
#include "universe.h"

#include "boost/ptr_container/ptr_vector.hpp"

#ifndef JOURNEY_H_
#define JOURNEY_H_



namespace EMTG 
{

    class journey 
    {
    public:
	    //constructor
	    journey();
	    journey(const int& j, const missionoptions& options);

	    //destructor
	    virtual ~journey();

	    //methods
	    //evaluate function
	    //return 0 if successful, 1 if failure
	    int evaluate(   const double* X,
                        int* Xindex, 
                        double* F, 
                        int* Findex, 
                        double* G, 
                        int* Gindex, 
                        const int& needG, 
                        const int& j, 
                        double* current_epoch,
                        double* current_state,
                        double* current_deltaV, 
                        EMTG::Astrodynamics::universe& Universe,
                        missionoptions* options);

	    //output functions
	    //main output function
	    void output(missionoptions* options,
                   const double& launchdate,
                   const int& j, 
                   int& jprint,
                   EMTG::Astrodynamics::universe& Universe, 
                   int* eventcount);

        //method to output "journey and a half" information that occurs while the spacecraft is "hanging out" at a body prior to departure
        void output_journey_prologue(missionoptions* options,
                                    const double& launchdate,
                                    const int& j,
                                    int& jprint,
                                    EMTG::Astrodynamics::universe& Universe,
                                    int* eventcount);

        //method to output additional stay time at the target body
        void output_journey_postlogue(  missionoptions* options,
                                        const double& launchdate,
                                        const int& j,
                                        int& jprint,
                                        EMTG::Astrodynamics::universe& Universe,
                                        int* eventcount);

        //method to output journey header
        void output_journey_header(missionoptions* options,
                                    EMTG::Astrodynamics::universe& Universe,
                                    const int& j,
                                    int& jprint,
                                    const int& waiting);

	    //bounds calculation function
	    //return 0 if successful, 1 if failure
	    void calcbounds(vector<double>* Xupperbounds,
                       vector<double>* Xlowerbounds,
                       vector<double>* Fupperbounds, 
                       vector<double>* Flowerbounds, 
                       vector<string>* Xdescriptions, 
                       vector<string>* Fdescriptions, 
                       vector<int>* iAfun, 
                       vector<int>* jAvar, 
                       vector<int>* iGfun, 
                       vector<int>* jGvar, 
                       vector<double>* A, 
                       vector<string>* Adescriptions,
                       vector<string>* Gdescriptions, 
                       vector<double>* synodic_periods,
                       int j, 
                       EMTG::Astrodynamics::universe& Universe,
                       missionoptions* options);

	    //function to find constraint dependecies due to an escape or capture spiral in this or a previous journey
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
												    int j,
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
												    int j,
												    missionoptions* options,
												    const int& Findex);

	    //method to create an initial guess of another mission type
	    void create_initial_guess(const int& desired_mission_type,
							    const bool& VSI,
							    double& current_epoch,
							    const int& j, 
							    vector<double>& NewX,
							    int& NewXIndex,
							    const vector<string>& NewXDescriptions,
							    const missionoptions& options,
                                const Astrodynamics::universe& Universe);

        //method to output a .e ephemeris file
        void write_ephemeris_file(  const missionoptions& options,
                                    const EMTG::Astrodynamics::universe& Universe,
                                    const double& launch_epoch,
                                    int j);

	    //vector of phases
	    boost::ptr_vector<phase> phases;

	    //vector of boundary states
	    vector< vector <double> > boundary_states;

	    //fields
	    int number_of_phases;
	    int journey_index;
	    double mu; //gravitational parameter of central body in km^3/s^2
	    string central_body_name; //name of central body
	    double LU; //characteristic length unit of central body (1 AU for sun, or 1 body radius for others)
	    double TU; //characteristic time unit of central body
	    double journey_initial_mass_increment_scale_factor;

	    //derivative information
	    int first_X_entry_in_journey;
	    vector<int> timeconstraints_G_indices;
	    vector<double> timeconstraints_X_scale_ranges;
    };

} /* namespace EMTG */
#endif /* JOURNEY_H_ */
