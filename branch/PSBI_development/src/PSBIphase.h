//header file for Parallel-Shooting Bounded Impulse (PSBI) phase
//a transcription which draws heritage from both parallel shooting and Sims-Flanagan
//Jacob Englander 10-3-2014

#include <vector>

#include "phase.h"
#include "journey.h"
#include "missionoptions.h"
#include "universe.h"
#include "STM.h"

#ifndef _PSBI_PHASE
#define _PSBI_PHASE

namespace EMTG
{
    class PSBIphase : public phase
    {
    public:
        //constructor
        PSBIphase();
        PSBIphase(int j, int p, missionoptions* options);

        //destructor
        virtual ~PSBIphase();

        //evaluate function
        //return 0 if successful, 1 if failure
        int evaluate(double* X, 
                    int* Xindex, 
                    double* F,
                    int* Findex,
                    double* G, 
                    int* Gindex,
                    int needG, 
                    double* current_epoch, 
                    double* current_state,
                    double* current_deltaV,
                    double* boundary1_state, 
                    double* boundary2_state, 
                    int j,
                    int p,
                    EMTG::Astrodynamics::universe* Universe,
                    missionoptions* options);

        //output function
        //return 0 if successful, 1 if failure
        int output(missionoptions* options,
                    const double& launchdate, 
                    int j,
                    int p,
                    EMTG::Astrodynamics::universe* Universe, 
                    int* eventcount);

        //bounds calculation function
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
                        vector<string>* Adescriptions, 
                        vector<string>* Gdescriptions, 
                        vector<double>* synodic_periods,
                        int j,
                        int p,
                        EMTG::Astrodynamics::universe* Universe,
                        missionoptions* options);

        //time information
        vector <double> event_epochs;

        //state information and constraint information
        vector< vector<double> > defect_constraint;
        vector<double> throttle;
        vector<double> dVmax;
        vector< vector<double> > dV;
        vector< vector<double> > ForceVector;

        //****************************************************************
        //derivative information
        
        //Kepler derivative information
            vector< Kepler::STM > STM; //vector of state transition matrices
            vector<double> Kepler_F, Kepler_Fdot, Kepler_G, Kepler_Gdot, Kepler_Fdotdot, Kepler_Gdotdot;
            vector<double> Propagation_Step_Time_Fraction;
            vector<double> Propagation_Step_Time_Fraction_Derivative;

        //derivatives for the radial distance constraint
            vector< vector<int> > radus_constraint_G_indicies;
            vector< vector <double> > radius_constraint_X_scale_ranges;

        //derivatives for all left-handed defect constraints
            vector< vector< vector<int> > > G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
            vector< vector< vector<double> > > X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;

            vector< vector<int> > G_index_of_derivative_of_defect_with_respect_to_BOL_power;

        //derivatives for the FIRST left-handed defect constraint only
            vector<int> G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_mission_initial_mass_multiplier;

            vector<int> G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_journey_initial_mass_increment_multiplier;

            vector< vector<int> > G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition;
            vector< vector<double> > X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition;

            vector< vector<int> > G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity;
            vector< vector<double> > X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity;

            vector<int> G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_previous_phase_arrival_mass;
            double X_scale_range_previous_phase_arrival_mass;

        //derivatives for SUCCESSIVE left-handed defect constraints only
            vector< vector< vector<int> > > G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;
            vector< vector< vector<double> > > X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;

        //derivatives of the right-handed defect constraint
            vector< vector<int> > G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables;
            vector< vector<double> > X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables;

            vector<int> G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_current_phase_arrival_mass;
            double X_scale_range_current_phase_arrival_mass;

            vector< vector<int> > G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity;
            vector< vector<double> > X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity;

            vector< vector<int> > G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition;
            vector< vector<double> > X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition;

        

    };
} // end namespace EMTG

#endif //_PSBI_PHASE