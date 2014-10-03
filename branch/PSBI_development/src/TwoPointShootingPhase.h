//header file for TwoPointShootingPhase
//this is an intermediate phase type that one should never instantiate - it contains features common to MGALT, FBLT, MGANDSM (for now)

#include "phase.h"
#include "journey.h"
#include "missionoptions.h"
#include "universe.h"
#include "STM.h"


#ifndef TWOPOINTSHOOTINGPHASE
#define TWOPOINTSHOOTINGPHASE

namespace EMTG {

    class TwoPointShootingPhase : public EMTG::phase {
    public:
        //constructor
        TwoPointShootingPhase();

        //destructor
        virtual ~TwoPointShootingPhase();

        //evaluate function
        //return 0 if successful, 1 if failure
        virtual int evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options) = 0;

        //output function
        //return 0 if successful, 1 if failure
        virtual int output(missionoptions* options, const double& launchdate, int j, int p, EMTG::Astrodynamics::universe* Universe, int* eventcount) = 0;

        //bounds calculation function
        //return 0 if successful, 1 if failure
        virtual void calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options) = 0;

        //top-level function to calculate the match point derivatives
        virtual void calculate_match_point_derivatives(double* G,
                                                        int* Gindex,
                                                        const int& j,
                                                        const int& p,
                                                        missionoptions* options,
                                                        EMTG::Astrodynamics::universe* Universe) {};

        //function to calculate the derivative of a match point constraint with respect to a decision variable in the forward propagation
        virtual void calculate_match_point_forward_propagation_derivatives(double* G,
                                                                        int* Gindex,
                                                                        const int& j,
                                                                        const int& p,
                                                                        missionoptions* options,
                                                                        EMTG::Astrodynamics::universe* Universe,
                                                                        const int& step,
                                                                        const int& stepnext,
                                                                        double& dxdu,
                                                                        double& dydu,
                                                                        double& dzdu,
                                                                        double& dxdotdu,
                                                                        double& dydotdu,
                                                                        double& dzdotdu,
                                                                        double& dmdu,
                                                                        double& dtdu,
                                                                        double& dtotal_available_thrust_time_du,
                                                                        double& dPdu) {};

        //function to calculate the derivative of a match point constraint with respect to a decision variable in the backward propagation
        virtual void calculate_match_point_backward_propagation_derivatives(double* G,
                                                                            int* Gindex,
                                                                            const int& j,
                                                                            const int& p,
                                                                            missionoptions* options,
                                                                            EMTG::Astrodynamics::universe* Universe,
                                                                            const int& backstep,
                                                                            const int& stepnext,
                                                                            double& dxdu,
                                                                            double& dydu,
                                                                            double& dzdu,
                                                                            double& dxdotdu,
                                                                            double& dydotdu,
                                                                            double& dzdotdu,
                                                                            double& dmdu,
                                                                            double& dtdu,
                                                                            double& dtotal_available_thrust_time_du,
                                                                            double& dPdu) {};

        //function to create an initial guess for another phase type
        virtual void create_initial_guess(const int& desired_mission_type,
                                        const bool& VSI,
                                        double& current_epoch,
                                        const int& j,
                                        const int& p,
                                        vector<double>& NewX,
                                        int& NewXIndex,
                                        const vector<string>& NewXDescriptions,
                                        const missionoptions& options)		{};

        //method to output a .e ephemeris file
        virtual void write_ephemeris_file(const missionoptions& options,
                                        const EMTG::Astrodynamics::universe& Universe,
                                        const double& launch_epoch,
                                        const double& journey_starting_epoch,
                                        vector< vector<string> >& output_line_array,
                                        const int& j,
                                        const int& p)	{};


        //functions that only exist for phases of this type
    protected:
        void calcbounds_phase_thruster_parameters(const string& prefix,
                                                int first_X_entry_in_phase,
                                                vector<double>* Xupperbounds,
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
                                                int p,
                                                EMTG::Astrodynamics::universe* Universe,
                                                missionoptions* options);

        void calcbounds_LT_match_points(const string& prefix,
                                        int first_X_entry_in_phase,
                                        vector<double>* Xupperbounds,
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
                                        int p,
                                        EMTG::Astrodynamics::universe* Universe,
                                        missionoptions* options);

        //time information
        vector <double> event_epochs;

        //state information
        vector<double> match_point_state;
        vector<double> throttle;
        vector<double> dVmax;
        vector< vector<double> > dV;
        vector< vector<double> > ForceVector;

        //match point constraints
        vector< vector< vector <int> > > match_point_constraint_G_indices;
        vector<double> match_point_constraint_X_scale_ranges;
        vector< Kepler::STM > Forward_STM; //vector of state transition matrices
        vector< Kepler::STM > Backward_STM; //vector of state transition matrices
        Kepler::STM Current_STM;
        vector<double> Kepler_F_Forward, Kepler_Fdot_Forward, Kepler_G_Forward, Kepler_Gdot_Forward, Kepler_Fdotdot_Forward, Kepler_Gdotdot_Forward;
        vector<double> Kepler_F_Backward, Kepler_Fdot_Backward, Kepler_G_Backward, Kepler_Gdot_Backward, Kepler_Fdotdot_Backward, Kepler_Gdotdot_Backward;
        double Kepler_F_Current, Kepler_Fdot_Current, Kepler_G_Current, Kepler_Gdot_Current, Kepler_Fdotdot_Current, Kepler_Gdotdot_Current;
        vector<double> Propagation_Step_Time_Fraction_Forward, Propagation_Step_Time_Fraction_Backward;
        vector<double> Propagation_Step_Time_Fraction_Derivative_Forward, Propagation_Step_Time_Fraction_Derivative_Backward;
        vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass;
        vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass;
        vector< vector<int> > G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables;
        vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier;
        vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier;
        vector<int> G_index_of_derivative_of_match_point_with_respect_to_BOL_power;
    };
}

#endif TWOPOINTSHOOTINGPHASE