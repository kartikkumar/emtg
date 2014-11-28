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
        int evaluate(const double* X,
                    int* Xindex,
                    double* F,
                    int* Findex,
                    double* G,
                    int* Gindex,
                    const int& needG,
                    double* current_epoch,
                    double* current_state,
                    double* current_deltaV,
                    double* boundary1_state,
                    double* boundary2_state,
                    const int& j,
                    const int& p,
                    EMTG::Astrodynamics::universe* Universe,
                    missionoptions* options);

        //output function
        //return 0 if successful, 1 if failure
        int output(missionoptions* options,
                    const double& launchdate,
                    const int& j,
                    const int& p,
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
                        const int& j,
                        const int& p,
                        EMTG::Astrodynamics::universe* Universe,
                        missionoptions* options);

        //method to calculate the derivatives of defect constraints
    protected:
        void calculate_defect_derivatives(double* G,
                                        const int& j,
                                        const int& p,
                                        missionoptions* options,
                                        EMTG::Astrodynamics::universe* Universe);

        //function to calculate the derivative of a match point constraint with respect to a decision variable in the forward propagation
        void calculate_derivative_of_right_hand_state(missionoptions* options,
                                                    EMTG::Astrodynamics::universe* Universe,
                                                    const int& step,
                                                    double& dxdu,
                                                    double& dydu,
                                                    double& dzdu,
                                                    double& dxdotdu,
                                                    double& dydotdu,
                                                    double& dzdotdu,
                                                    double& dmdu,
                                                    double& dtdu,
                                                    double& dtotal_available_thrust_time_du,
                                                    double& dPdu);
        //function to create an initial guess for another phase type
        virtual void create_initial_guess(const int& desired_mission_type,
            const bool& VSI,
            double& current_epoch,
            const int& j,
            const int& p,
            vector<double>& NewX,
            int& NewXIndex,
            const vector<string>& NewXDescriptions,
            const missionoptions& options,
            const Astrodynamics::universe& Universe);

        //time information
        vector <double> event_epochs;

        //state information and constraint information
        double minimum_radial_distance;
        vector< vector<double> > left_hand_state;
        vector< vector<double> > right_hand_state;
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
            Kepler::STM* initial_coast_STM;
            double initial_coast_Kepler_F, initial_coast_Kepler_Fdot, initial_coast_Kepler_G, initial_coast_Kepler_Gdot, initial_coast_Kepler_Fdotdot, initial_coast_Kepler_Gdotdot;
            Kepler::STM* terminal_coast_STM;
            double terminal_coast_Kepler_F, terminal_coast_Kepler_Fdot, terminal_coast_Kepler_G, terminal_coast_Kepler_Gdot, terminal_coast_Kepler_Fdotdot, terminal_coast_Kepler_Gdotdot;
            bool initial_coast;
            bool terminal_coast;
            vector<double> Propagation_Step_Time_Fraction;
            vector<double> Propagation_Step_Time_Fraction_Derivative;

        //derivatives for the radial distance constraint
            vector< vector<int> > radius_constraint_G_indices;
            vector< vector <double> > radius_constraint_X_scale_ranges;

        //derivatives for all left-handed defect constraints
            vector< vector<int> > G_index_of_derivative_of_defect_constraints_with_respect_to_current_state;
            vector< vector<double> > X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state;
            vector< vector< vector<int> > > G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
            vector< vector< vector<double> > > X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;

            vector< vector<int> > G_index_of_derivative_of_defect_with_respect_to_BOL_power;

        //derivatives for the FIRST left-handed defect constraint only
            int G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_mission_initial_mass_multiplier;
            double X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_mission_initial_mass_multiplier;

            int G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_journey_initial_mass_increment_multiplier;

            vector< vector<int> > G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition;
            vector< vector<double> > X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition;

            vector< vector<int> > G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity;
            vector< vector<double> > X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity;

            int G_index_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_phase_arrival_mass;
            double X_scale_range_previous_phase_arrival_mass;

            vector < vector<int> > G_index_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition;
            vector < vector<double> > X_scale_range_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition;

        //derivatives for SUCCESSIVE left-handed defect constraints only
            vector< vector< vector<int> > > G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state;
            vector< vector< vector<double> > > X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state;
            vector< vector< vector<int> > > G_index_of_derivative_of_defect_constraints_with_respect_to_previous_control;
            vector< vector< vector<double> > > X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_control;

        //derivatives of the right-handed defect constraint
            vector< vector<int> > G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control;
            vector< vector<double> > X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control;

            vector< vector<int> > G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables;
            vector< vector<double> > X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables;

            int G_index_of_derivative_of_rightmost_mass_defect_constraint_with_respect_to_current_phase_arrival_mass;
            double X_scale_range_of_derivative_of_rightmost_mass_defect_constraint_with_respect_to_current_phase_arrival_mass;

            vector< vector<int> > G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity;
            vector< vector<double> > X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity;

            vector< vector<int> > G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition;
            vector< vector<double> > X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition;

        

    };
} // end namespace EMTG

#endif //_PSBI_PHASE