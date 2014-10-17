//Parallel - Shooting Bounded Impulse(PSBI) phase
//a transcription which draws heritage from both parallel shooting and Sims-Flanagan
//Jacob Englander 10-3-2014


#include <sstream>
#include <fstream>

#include "PSBIphase.h"
#include "Astrodynamics.h"
#include "missionoptions.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "universe.h"
#include "EMTG_Matrix.h"
#include "Kepler_Lagrange_Laguerre_Conway_Der.h"

#include "SpiceUsr.h"

namespace EMTG
{
    //default constructor - never used
    PSBIphase::PSBIphase() {}

    //main constructor
    PSBIphase::PSBIphase(int j, int p, missionoptions* options)
    {
        //must resize all data vectors to the correct length
        vector<double> state_dummy(7);
        vector<double> dV_or_control_dummy(3);

        for (size_t step = 0; step < options->num_timesteps; ++step)
        {
            this->left_hand_state.push_back(state_dummy);
            this->spacecraft_state.push_back(state_dummy);
            this->right_hand_state.push_back(state_dummy);
            this->dV.push_back(dV_or_control_dummy);
            this->control.push_back(dV_or_control_dummy);
            this->ForceVector.push_back(dV_or_control_dummy);
            this->dagravdRvec.push_back(dV_or_control_dummy);
            this->dagravdtvec.push_back(dV_or_control_dummy);
        }

        this->event_epochs.resize(options->num_timesteps);
        this->dVmax.resize(options->num_timesteps);
        this->available_power.resize(options->num_timesteps);
        this->available_mass_flow_rate.resize(options->num_timesteps);
        this->available_thrust.resize(options->num_timesteps);
        this->available_Isp.resize(options->num_timesteps);
        this->active_power.resize(options->num_timesteps);
        this->number_of_active_engines.resize(options->num_timesteps);
        this->throttle.resize(options->num_timesteps);

        //size the vectors that will be used to calculate the b-plane
        this->V_infinity_in.resize(3, 1);
        this->V_infinity_out.resize(3, 1);
        this->BoundaryR.resize(3, 1);
        this->BoundaryV.resize(3, 1);

        //set the bodies
        this->boundary1_location_code = options->sequence[j][p];
        this->boundary2_location_code = options->sequence[j][p + 1];

        //size the vectors of state transition matrices
        this->STM.resize(options->num_timesteps * 2);
        this->Kepler_F.resize(options->num_timesteps * 2);
        this->Kepler_Fdot.resize(options->num_timesteps * 2);
        this->Kepler_G.resize(options->num_timesteps * 2);
        this->Kepler_Gdot.resize(options->num_timesteps * 2);
        this->Kepler_Fdotdot.resize(options->num_timesteps * 2);
        this->Kepler_Gdotdot.resize(options->num_timesteps * 2);
        this->Propagation_Step_Time_Fraction.resize(options->num_timesteps * 2);
        this->Propagation_Step_Time_Fraction.resize(options->num_timesteps * 2);

        //if there are initial and/or terminal coasts, instantiate the relevant STMs
        this->initial_coast = false;
        this->terminal_coast = false;
        if (j == 0 && p == 0 && options->forced_post_launch_coast > math::SMALL)
        {
            this->initial_coast_STM = new Kepler::STM;
            this->initial_coast = true;
        }
        else if ((p > 0 || options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4 || options->journey_departure_type[j] == 6)
            && options->forced_flyby_coast > math::SMALL)
        {
            this->initial_coast_STM = new Kepler::STM;
            this->initial_coast = true;
        }
        if ((p < options->number_of_phases[j] - 1 || (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5))
            && options->forced_flyby_coast > math::SMALL)
        {
            this->terminal_coast_STM = new Kepler::STM;
            this->terminal_coast = true;
        }

        this->current_mass_increment = 0.0;
        this->journey_initial_mass_increment_scale_factor = 1.0;

        //size the time step vector
        time_step_sizes.resize(options->num_timesteps * 2);

        this->dTdP.resize(options->num_timesteps);
        this->dmdotdP.resize(options->num_timesteps);
        this->dTdIsp.resize(options->num_timesteps);
        this->dmdotdIsp.resize(options->num_timesteps);
        this->dPdr.resize(options->num_timesteps);
        this->dPdt.resize(options->num_timesteps);
        this->dFSRPdr.resize(options->num_timesteps);

        //set derivatives for spirals
        this->spiral_escape_dm_after_dm_before = 1.0;
    }

    //destructor
    PSBIphase::~PSBIphase()
    {
        //if there are initial and/or terminal coasts, delete the relevant STMs
        if (this->initial_coast)
            delete this->initial_coast_STM;
        if (this->terminal_coast)
            delete this->terminal_coast_STM;
    }

    //bounds calculation method
    void PSBIphase::calcbounds(vector<double>* Xupperbounds,
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
                                missionoptions* options)
    {
        //this function calculates the upper and lower bounds for the decision and constraint vectors for PSBI
        //create a prefix string with journey and phase information
        stringstream prefixstream;
        prefixstream << "j" << j << "p" << p << ": ";
        string prefix = prefixstream.str();
        int first_X_entry_in_phase = Xupperbounds->size();

        //initialize an array of strings describing the defect constraints
        vector<string> statename;
        statename.push_back("x");
        statename.push_back("y");
        statename.push_back("z");
        statename.push_back("xdot");
        statename.push_back("ydot");
        statename.push_back("zdot");
        statename.push_back("mass");

        //**************************************************************************
        //calculate bounds on variables and constraints governing the left boundary
        this->calcbounds_left_boundary(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

        //**************************************************************************
        //if EMTG is choosing an input power or Isp for the phase (for REP/NEP models), then this information must be encoded
        this->calcbounds_phase_thruster_parameters(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

        //**************************************************************************
        //next, we need to encode the phase flight time
        this->calcbounds_flight_time(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, synodic_periods, j, p, Universe, options);

        //**************************************************************************
        //calculate bounds on variables and constraints governing the right boundary
        this->calcbounds_right_boundary(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

        //**************************************************************************
        //if we are using a variable scale/standard deviation distribution to pick the control points
        if (options->step_size_distribution == 3 || options->step_size_distribution == 4)
        {
            Xlowerbounds->push_back(1.0);
            Xupperbounds->push_back(options->step_size_stdv_or_scale);
            Xdescriptions->push_back(prefix + "step size distribution standard deviation or scale factor");
        }

        //**************************************************************************
        //now encode the bounds on the state and control variables, plus the defect constraints, for each control point
        //points are given in spherical coordinates
        for (size_t step = 0; step < options->num_timesteps; ++step)
        {
            stringstream stepstream;
            stepstream << step;

            //define the maximum distance from the central body
            //first check to see if either boundary orbit is a hyperbola
            //if so, use the pseudoa for that boundary orbit
            //if not, use twice the pseudoa (i.e. twice the apoapse distance)
            double maxd1 = this->e1 >= 1.0 ? this->pseudoa1 : 2.0 * this->pseudoa1;
            double maxd2 = this->e2 >= 1.0 ? this->pseudoa2 : 2.0 * this->pseudoa2;
            //then use the larger of these values betweeen boundary 1 and boundary 2
            double max_distance = max(maxd1, maxd2);

            //define the maximum velocity
            
            //first check to see if either boundary orbit is a hyperbola
            //if so the maximum velocity for that boundary is the hyperbola periapse velocity
            //if not use twice the periapse velocity (i.e. for bounded orbits)
            double maxv1 = this->e1 >= 1.0 ? sqrt(Universe->mu * (2.0 / (this->a1 * (1 - this->e1)) - 1.0 / this->a1)) : 2.0 * sqrt(Universe->mu * (2.0 / (this->a1 * (1 - this->e1)) - 1.0 / this->a1));
            double maxv2 = this->e2 >= 1.0 ? sqrt(Universe->mu * (2.0 / (this->a2 * (1 - this->e2)) - 1.0 / this->a2)) : 2.0 * sqrt(Universe->mu * (2.0 / (this->a2 * (1 - this->e2)) - 1.0 / this->a2));
            //then use the larger of these values betweeen boundary 1 and boundary 2
            double max_velocity = max(maxv1, maxv2);

            //x
            Xlowerbounds->push_back(-max_distance);
            Xupperbounds->push_back(max_distance);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " x");

            //y
            Xlowerbounds->push_back(-max_distance);
            Xupperbounds->push_back(max_distance);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " y");

            //z
            Xlowerbounds->push_back(-max_distance);
            Xupperbounds->push_back(max_distance);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " z");

            //xdot
            Xlowerbounds->push_back(-max_velocity);
            Xupperbounds->push_back(max_velocity);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " xdot");

            //ydot
            Xlowerbounds->push_back(-max_velocity);
            Xupperbounds->push_back(max_velocity);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " ydot");

            //zdot
            Xlowerbounds->push_back(-max_velocity);
            Xupperbounds->push_back(max_velocity);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " zdot");

            //mass
            //count up the mass increments preceding this phase
            current_mass_increment = 0.0;
            for (int jj = 0; jj <= j; ++jj)
                current_mass_increment += options->journey_starting_mass_increment[jj];

            if (options->journey_arrival_type[j] == 0) //chemical orbit insertion, so we need to drop the EP dry mass (and therefore have it available to drop)
                Xlowerbounds->push_back(options->EP_dry_mass);
            else
                Xlowerbounds->push_back(EMTG::math::SMALL);
            Xupperbounds->push_back(options->maximum_mass + current_mass_increment);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " mass");

            //u_x
            Xlowerbounds->push_back(-1.0);
            Xupperbounds->push_back(1.0);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_x");

            //u_y
            Xlowerbounds->push_back(-1.0);
            Xupperbounds->push_back(1.0);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_y");

            //u_z
            Xlowerbounds->push_back(-1.0);
            Xupperbounds->push_back(1.0);
            Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_z");

            //for variable specific impulse propulsion systems, we must also encode the Isp for this time step
            int state_and_control_elements_this_step = 10;
            if (options->engine_type == 4 || options->engine_type == 13)
            {
                Xlowerbounds->push_back(options->IspLT_minimum);
                Xupperbounds->push_back(options->IspLT);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
                ++state_and_control_elements_this_step;
            }
            else if (options->engine_type == 12)
            {
                Xlowerbounds->push_back(3000.0);
                Xupperbounds->push_back(5000.0);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
                ++state_and_control_elements_this_step;
            }

            //position magnitude constraint
            //radial distance must always be greater than central body minimum
            //or, later, whatever user-specified minimum distance we have
            this->minimum_radial_distance = Universe->minimum_safe_distance;
            Flowerbounds->push_back(-math::LARGE);
            Fupperbounds->push_back(0.0);
            Fdescriptions->push_back(prefix + "step " + stepstream.str() + " radial distance constraint");
            //radial distance constraint is dependent ONLY on x, y, z for THIS step
            vector<int> step_radius_constraint_G_indices;
            vector<double> step_radius_constraint_X_scale_ranges;
            for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
            {
                if (((*Xdescriptions)[entry].find(prefix + "step " + stepstream.str() + " x") < 1024 && !((*Xdescriptions)[entry].find("dot") < 1024))
                    || ((*Xdescriptions)[entry].find(prefix + "step " + stepstream.str() + " y") < 1024 && !((*Xdescriptions)[entry].find("dot") < 1024))
                    || ((*Xdescriptions)[entry].find(prefix + "step " + stepstream.str() + " z") < 1024 && !((*Xdescriptions)[entry].find("dot") < 1024)))
                {
                    iGfun->push_back(Fdescriptions->size() - 1);
                    jGvar->push_back(entry);
                    stringstream EntryNameStream;
                    EntryNameStream << "Derivative of " << prefix + "step " << step << " radial distance constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                    Gdescriptions->push_back(EntryNameStream.str());
                    step_radius_constraint_G_indices.push_back(Gdescriptions->size() - 1);
                    step_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                }
            }
            this->radius_constraint_G_indices.push_back(step_radius_constraint_G_indices);
            this->radius_constraint_X_scale_ranges.push_back(step_radius_constraint_X_scale_ranges);


            //left-hand defect constraints: all steps have this!
            //temporary arrays for G index tracking
            
            vector<int> step_G_index_of_derivative_of_defect_constraints_with_respect_to_current_state;
            vector<double> step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state;
            vector< vector<int> > step_G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
            vector< vector<double> > step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
            vector<int> step_G_index_of_derivative_of_defect_with_respect_to_BOL_power(7);
            vector< vector<int> > step_G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;
            vector< vector<double> >step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;

            for (size_t state = 0; state < 7; ++state)
            {
                Flowerbounds->push_back(-math::SMALL);
                Fupperbounds->push_back(math::SMALL);
                Fdescriptions->push_back(prefix + "step " + stepstream.str() + " defect " + statename[state]);

                //the defect constraint is ALWAYS dependent on:
                //for all defects
                //  CURRENT left-hand state variable matching this defect
                vector<int> state_G_index_of_derivative_of_defect_constraints_with_respect_to_current_state;
                vector<double> state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state;
                for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                {
                    if ((*Xdescriptions)[entry] == prefix + "step " + stepstream.str() + " " + statename[state])
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << "step " << step << " " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        step_G_index_of_derivative_of_defect_constraints_with_respect_to_current_state.push_back(Gdescriptions->size() - 1);
                        step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }
                }
                
                //  ALL previous time variables including the current phase flight time
                //  note that the FIRST entry in the list of time derivatives is with respect to the CURRENT phase flight time
                vector<int> state_G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
                vector<double> state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
                for (int entry = Xdescriptions->size() - 1; entry >= 0; --entry)
                {
                    if (((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024) && (step > 0 || state < 6))
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << "step " << step << " " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        state_G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables.push_back(Gdescriptions->size() - 1);
                        state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }
                }
                step_G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables.push_back(state_G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables);
                step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables.push_back(state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables);

                //  mission/journey global Isp and/or power variables
                if (options->objective_type == 13)
                {
                    for (size_t Xentry = 0; Xentry < Xdescriptions->size(); ++Xentry)
                    {
                        if ((*Xdescriptions)[Xentry].find("engine input power (kW)") < 1024)
                        {
                            bool duplicateflag = false;
                            stringstream entry_tag_stream;
                            entry_tag_stream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
                            for (int Gentry = Gdescriptions->size() - 1; Gentry >= 0; --Gentry)
                            {
                                if ((*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
                                {
                                    duplicateflag = true;
                                    step_G_index_of_derivative_of_defect_with_respect_to_BOL_power[state] = Gentry;
                                    this->power_range = (*Xupperbounds)[Xentry] - math::SMALL;
                                    break;
                                }
                            }
                            if (!duplicateflag)
                            {
                                iGfun->push_back(Fdescriptions->size() - 1);
                                jGvar->push_back(Xentry);
                                stringstream EntryNameStream;
                                EntryNameStream << "Derivative of " << prefix << "step " << step << " " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
                                Gdescriptions->push_back(EntryNameStream.str());
                                step_G_index_of_derivative_of_defect_with_respect_to_BOL_power[state] = Gdescriptions->size() - 1;
                                this->power_range = (*Xupperbounds)[Xentry] - math::SMALL;
                                break;
                            }
                        }
                    }
                }
                //  if applicable, spiral things affect every defect constraint
                this->find_dependencies_due_to_escape_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
                this->find_dependencies_due_to_capture_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);

                //for the first defect in the phase only
                
                vector<int> state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition;
                vector<double> state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition;
                vector<int> state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity;
                vector<double> state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity;
                vector<int> state_G_index_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition;
                vector<double> state_X_scale_range_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition;
                
                if (step == 0)
                {
                    //if applicable, variable left hand boundary condition
                    if (p == 0 && state < 6 &&
                                    (options->journey_departure_elements_vary_flag[j][0]
                                    || options->journey_departure_elements_vary_flag[j][1]
                                    || options->journey_departure_elements_vary_flag[j][2]
                                    || options->journey_departure_elements_vary_flag[j][3]
                                    || options->journey_departure_elements_vary_flag[j][4]
                                    || options->journey_departure_elements_vary_flag[j][5]))
                    {
                        for (int entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                        {
                            if ((*Xdescriptions)[entry].find("left boundary point") < 1024)
                            {
                                iGfun->push_back(Fdescriptions->size() - 1);
                                jGvar->push_back(entry);
                                stringstream EntryNameStream;
                                EntryNameStream << "Derivative of leftmost defect constraint on state " << statename[state] << "F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                Gdescriptions->push_back(EntryNameStream.str());
                                state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition.push_back(Gdescriptions->size() - 1);
                                state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                            }
                        }
                        this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition.push_back(state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition);
                        this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition.push_back(state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_variable_left_boundary_condition);
                    }//end  block for variable left-hand boundary

                    //if applicable, previous journey variable right-hand boundary condition
                    //this occurs for first phase of successive journeys where previous journey ended in a variable boundary condition
                    if (p == 0 && j > 0)
                    {
                        if (options->destination_list[j][0] == -1 && options->destination_list[j - 1][1] == -1 && state < 6
                            && (options->journey_arrival_elements_vary_flag[j - 1][0]
                                || options->journey_arrival_elements_vary_flag[j - 1][1]
                                || options->journey_arrival_elements_vary_flag[j - 1][2]
                                || options->journey_arrival_elements_vary_flag[j - 1][3]
                                || options->journey_arrival_elements_vary_flag[j - 1][4]
                                || options->journey_arrival_elements_vary_flag[j - 1][5]))
                        {
                            for (int entry = first_X_entry_in_phase; entry >=0; --entry)
                            {
                                stringstream previous_journey_prefix_stream;
                                previous_journey_prefix_stream << "j" << j - 1 << "p" << options->number_of_phases[j - 1] - 1 << ": ";
                                if ((*Xdescriptions)[entry].find(previous_journey_prefix_stream.str()) < 1024 && (*Xdescriptions)[entry].find("right boundary point") < 1024)
                                {
                                    iGfun->push_back(Fdescriptions->size() - 1);
                                    jGvar->push_back(entry);
                                    stringstream EntryNameStream;
                                    EntryNameStream << "Derivative of leftmost defect constraint on state " << statename[state] << "F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                    Gdescriptions->push_back(EntryNameStream.str());
                                    state_G_index_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition.push_back(Gdescriptions->size() - 1);
                                    state_X_scale_range_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                                }
                            }
                            this->G_index_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition.push_back(state_G_index_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition);
                            this->X_scale_range_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition.push_back(state_X_scale_range_of_derivative_of_leftmost_defect_constraint_with_respect_to_previous_journey_variable_right_boundary_condition);
                        }
                    }//end block for previous journey variable right-hand boundary condition

                    //if applicable, current phase initial v-infinity (which for the first phase is in polar, successive phases in cartesian)
                    //note this is applicable whenever we are not in the first phase of a journey beginning with a free direct departure or departure spiral
                    if (!(p == 0 && (options->journey_departure_type[j] == 2 || options->journey_departure_type[j] == 5)))
                    {
                        //the first phase of the journey has velocity given in polar coordinates
                        if (p == 0)
                        {
                            for (int entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                            {
                                if ((*Xdescriptions)[entry].find("magnitude of outgoing velocity asymptote") < 1024)
                                {
                                    iGfun->push_back(Fdescriptions->size() - 1);
                                    jGvar->push_back(entry);
                                    stringstream EntryNameStream;
                                    EntryNameStream << "Derivative of "<< prefix << statename[state] << " leftmost defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                    Gdescriptions->push_back(EntryNameStream.str());
                                    state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back(Gdescriptions->size() - 1);
                                    state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                                }
                                else if ((*Xdescriptions)[entry].find("RA of departure asymptote") < 1024 && state < 6)
                                {
                                    iGfun->push_back(Fdescriptions->size() - 1);
                                    jGvar->push_back(entry);
                                    stringstream EntryNameStream;
                                    EntryNameStream << "Derivative of "<< prefix << statename[state] << " leftmost defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                    Gdescriptions->push_back(EntryNameStream.str());
                                    state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back(Gdescriptions->size() - 1);
                                    state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                                }
                                else if ((*Xdescriptions)[entry].find("DEC of departure asymptote") < 1024 && state < 6)
                                {
                                    iGfun->push_back(Fdescriptions->size() - 1);
                                    jGvar->push_back(entry);
                                    stringstream EntryNameStream;
                                    EntryNameStream << "Derivative of " << prefix << statename[state] << " leftmost defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                    Gdescriptions->push_back(EntryNameStream.str());
                                    state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back(Gdescriptions->size() - 1);
                                    state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                                }
                            }
                        }
                        //successive phases have velocity given in cartesian coordinates
                        else
                        {
                            for (int entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                            {
                                if ((*Xdescriptions)[entry].find("initial velocity increment") < 1024 && state < 6)
                                {
                                    iGfun->push_back(Fdescriptions->size() - 1);
                                    jGvar->push_back(entry);
                                    stringstream EntryNameStream;
                                    EntryNameStream << "Derivative of " << prefix << statename[state] << " leftmost defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                    Gdescriptions->push_back(EntryNameStream.str());
                                    state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back(Gdescriptions->size() - 1);
                                    state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                                }
                            }
                        }
                        this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back(state_G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity);
                        this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity.push_back(state_X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity);
                    }//end block for current phase initial v-infinity

                    //previous phase arrival mass, if this is NOT the first phase in the mission
                    if (j > 0 || p > 0)
                    {
                        for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
                        {
                            if ((*Xdescriptions)[entry].find("arrival_mass") < 1024 && state == 6)
                            {
                                iGfun->push_back(Fdescriptions->size() - 1);
                                jGvar->push_back(entry);
                                stringstream EntryNameStream;
                                EntryNameStream << "Derivative of " << prefix << statename[state] << " leftmost defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                Gdescriptions->push_back(EntryNameStream.str());
                                this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_previous_phase_arrival_mass.push_back(Gdescriptions->size() - 1);
                                this->X_scale_range_previous_phase_arrival_mass = (*Xupperbounds)[entry] - (*Xlowerbounds)[entry];
                                break;
                            }
                        }
                    }//end block for previous phase arrival mass

                    //if applicable, variable mission initial mass
                    //note that this can only happen in the first journey and phase because phases are separable
                    if (p == 0 && j == 0 && options->allow_initial_mass_to_vary && state == 6)
                    {
                        //step forward through the decision vector until you hit the initial mass multiplier
                        for (size_t entry = 0; entry < Xdescriptions->size() - 1; ++entry)
                        {
                            if ((*Xdescriptions)[entry].find("initial mass multiplier") < 1024)
                            {
                                iGfun->push_back(Fdescriptions->size() - 1);
                                jGvar->push_back(entry);
                                stringstream EntryNameStream;
                                EntryNameStream << "Derivative of " << prefix << statename[state] << " leftmost defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                Gdescriptions->push_back(EntryNameStream.str());
                                this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_mission_initial_mass_multiplier.push_back(Gdescriptions->size() - 1);
                                break;
                            }
                        }
                    }//close if block for variable initial mass

                    //if applicable, variable journey initial mass increment scale factor
                    if (options->journey_variable_mass_increment[j] && state == 6)
                    {
                        //note that the current phases's initial mass scales linearly with ALL previous journey initial mass multipliers
                        for (int pj = 0; pj <= j; ++pj)
                        {
                            for (int pp = 0; pp < (pj == j ? p : options->number_of_phases[pj]); ++pp)
                            {
                                stringstream pprefix_stream;
                                pprefix_stream << "j" << pj << "p" << pp;
                                string pprefix = pprefix_stream.str();

                                for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
                                {
                                    if ((*Xdescriptions)[entry].find(pprefix) < 1024 && (*Xdescriptions)[entry].find("journey initial mass scale factor") < 1024)
                                    {
                                        iGfun->push_back(Fdescriptions->size() - 1);
                                        jGvar->push_back(entry);
                                        stringstream EntryNameStream;
                                        EntryNameStream << "Derivative of " << prefix << statename[state] << " leftmost defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                        Gdescriptions->push_back(EntryNameStream.str());
                                        this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_journey_initial_mass_increment_multiplier.push_back(Gdescriptions->size() - 1);
                                    }
                                }
                            }
                        }
                    }//end leftmost defect constraint dependence on journey initial mass multiplier
                }
                //for successive defects only
                //  the previous control point's state and control
                else
                {
                    vector<int> state_G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;
                    vector<double> state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;
                    for (size_t entry = Xdescriptions->size() - 2 * state_and_control_elements_this_step; entry < Xdescriptions->size() - state_and_control_elements_this_step; ++entry)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix + "step " << step << " " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        state_G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control.push_back(Gdescriptions->size() - 1);
                        state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }

                    //put the sparsity entries into the step holder
                    step_G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control.push_back(state_G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control);
                    step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control.push_back(state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control);
                }
            }

            //now append the sparsity entries for this step's defect constraints to the phase's archive of sparsity entries
            this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state.push_back(step_G_index_of_derivative_of_defect_constraints_with_respect_to_current_state);
            this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state.push_back(step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state);
            this->G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables.push_back(step_G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables);
            this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables.push_back(step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables);
            this->G_index_of_derivative_of_defect_with_respect_to_BOL_power.push_back(step_G_index_of_derivative_of_defect_with_respect_to_BOL_power);
            this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control.push_back(step_G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control);
            this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control.push_back(step_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control);

            //throttle magnitude constraint
            //throttle = 0
            Flowerbounds->push_back(0.0);
            Fupperbounds->push_back(1.0);
            Fdescriptions->push_back(prefix + "step " + stepstream.str() + " throttle magnitude constraint");

            //Jacobian entries for the throttle magnitude constraint
            //The throttle magnitude constraint is dependent only on the preceding three throttle components
            vector<int> step_G_indices;
            int vary_Isp_flag = (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13) ? 1 : 0;
            for (size_t entry = Xdescriptions->size() - 1 - vary_Isp_flag; entry > Xdescriptions->size() - 4 - vary_Isp_flag; --entry)
            {
                iGfun->push_back(Fdescriptions->size() - 1);
                jGvar->push_back(entry);
                stringstream EntryNameStream;
                EntryNameStream << "Derivative of " << prefix + "step " << step << " throttle magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                Gdescriptions->push_back(EntryNameStream.str());

                //store the position in the G vector
                step_G_indices.push_back(Gdescriptions->size() - 1);
            }
            control_vector_G_indices.push_back(step_G_indices);
        }

        //**************************************************************************
        //right hand side defect constraint - the last segment only must connect to the right-hand boundary condition for the phase
        
        for (size_t state = 0; state < 7; ++state)
        {
            Flowerbounds->push_back(-math::SMALL);
            Fupperbounds->push_back(math::SMALL);
            Fdescriptions->push_back(prefix + "rightmost " + statename[state] + " defect_constraint");

            //derivatives
            //rightmost step state and control variables
            //this will just be the last N decision vectors
            vector<int> state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control;
            vector<double> state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control;
            int state_and_control_elements = options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13 ? 11 : 10;
            for (size_t entry = Xdescriptions->size() - state_and_control_elements; entry < Xdescriptions->size(); ++entry)
            {
                iGfun->push_back(Fdescriptions->size() - 1);
                jGvar->push_back(entry);
                stringstream EntryNameStream;
                EntryNameStream << "Derivative of " << prefix + "rightmost " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                Gdescriptions->push_back(EntryNameStream.str());
                state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control.push_back(Gdescriptions->size() - 1);
                state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
            }
            this->G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control.push_back(state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control);
            this->X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control.push_back(state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_rightmost_state_and_control);

            //  flight time (all previous, first entry is CURRENT phase flight time)
            vector<int> state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables;
            vector<double> state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables;
            for (int entry = Xdescriptions->size() - 1; entry >= 0; --entry)
            {
                if ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024)
                {
                    iGfun->push_back(Fdescriptions->size() - 1);
                    jGvar->push_back(entry);
                    stringstream EntryNameStream;
                    EntryNameStream << "Derivative of " << prefix << " rightmost " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                    Gdescriptions->push_back(EntryNameStream.str());
                    state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables.push_back(Gdescriptions->size() - 1);
                    state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                }
            }
            this->G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables.push_back(state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables);
            this->X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables.push_back(state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_flight_time_variables);

            //CURRENT phase arrival mass - only relevant for mass defect
            if (state == 6)
            {
                for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                {
                    if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << " rightmost " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        this->G_index_of_derivative_of_rightmost_mass_defect_constraint_with_respect_to_current_phase_arrival_mass == Gdescriptions->size() - 1;
                        this->X_scale_range_of_derivative_of_rightmost_mass_defect_constraint_with_respect_to_current_phase_arrival_mass = (*Xupperbounds)[entry] - (*Xlowerbounds)[entry];
                        break;
                    }
                }
            }

            //final v-infinity vector (no dependency exists for the mass defect
            //there is always a final v-infinity vector in the decision vector unless:
            //a) this phase ends in a low-thrust rendezvous (option 3)
            //b) this phase ends in a low-thrust fixed v-infinity intercept (option 5)
            //c) this phase ends in a capture spiral (option 7)
            if (!(p == options->number_of_phases[j] - 1 && state < 6 && (options->journey_arrival_type[j] == 3
                                                                        || options->journey_arrival_type[j] == 5
                                                                        || options->journey_arrival_type[j] == 7)))
            {
                vector<int> state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity;
                vector<double> state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity;
                for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                {
                    if ((*Xdescriptions)[entry].find("terminal velocity increment") < 1024)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << " rightmost " << statename[state] << " defect constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity.push_back(Gdescriptions->size() - 1);
                        state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }
                }
                this->G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity.push_back(state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity);
                this->X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity.push_back(state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_phase_terminal_velocity);
            }

            //if applicable, variable right hand boundary condition
            //(no dependency exists for the mass defect)
            if (p == options->number_of_phases[j] - 1 && state < 6 && (options->journey_arrival_elements_vary_flag[j][0]
                                                                        || options->journey_arrival_elements_vary_flag[j][1]
                                                                        || options->journey_arrival_elements_vary_flag[j][2]
                                                                        || options->journey_arrival_elements_vary_flag[j][3]
                                                                        || options->journey_arrival_elements_vary_flag[j][4]
                                                                        || options->journey_arrival_elements_vary_flag[j][5]))
            {
                vector<int> state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition;
                vector<double> state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition;

                for (int entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                {
                    if ((*Xdescriptions)[entry].find("right boundary point") < 1024)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of rightmost defect constraint on state " << statename[state] << "F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition.push_back(Gdescriptions->size() - 1);
                        state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }
                }
                this->G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition.push_back(state_G_index_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition);
                this->X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition.push_back(state_X_scale_range_of_derivative_of_rightmost_defect_constraints_with_respect_to_variable_right_boundary_condition);
            }//end  block for variable right-hand boundary

            //if applicable, spiral things affect the right-hand defect constraint (time is mass and mass is time)
            this->find_dependencies_due_to_escape_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
            this->find_dependencies_due_to_capture_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
        }//end loop over state for right-hand defect constraint

        
    }

    //evaluate function
    //return 0 if successful, 1 if failure
    int PSBIphase::evaluate(double* X,
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
                            missionoptions* options)
    {
        //declare some local variables
        int errcode = 0;
        double temp_available_thrust;
        double temp_available_mass_flow_rate;
        double temp_available_Isp;
        double temp_available_power;
        double temp_active_power;
        int temp_number_of_active_engines;
        double local_throttle, impulse_magnitude;

        //******************************************************************
        //Steps 1-4: Process the left boundary condition
        process_left_boundary_condition(X, Xindex, F, Findex, G, Gindex, needG, current_epoch, current_state, current_deltaV, boundary1_state, boundary2_state, j, p, Universe, options);

        //******************************************************************
        //Step 5: Process the right boundary condition
        process_right_boundary_condition(X, Xindex, F, Findex, G, Gindex, needG, current_epoch, current_state, current_deltaV, boundary1_state, boundary2_state, j, p, Universe, options);

        //******************************************************************
        //Step 6: propagate and thrust, filling in defect constraints as we go

        //Step 6.1: determine the length of a timestep
        double step_size_normalization_coefficient = 0.0;
        if (options->step_size_distribution == 3 || options->step_size_distribution == 4)
        {
            time_step_distribution_scale_or_stdv = X[*Xindex];
            ++(*Xindex);
        }
        else
            time_step_distribution_scale_or_stdv = options->step_size_stdv_or_scale;

        for (int step = 0; step < options->num_timesteps; ++step)
        {
            time_step_sizes[step] = compute_timestep_width_from_distribution(step + 0.5, options, time_step_distribution_scale_or_stdv);
            step_size_normalization_coefficient += time_step_sizes[step];
        }

        double total_available_thrust_time = TOF;
        if (j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6)
            total_available_thrust_time -= options->forced_post_launch_coast;
        else if ((p > 0 || p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4 || options->journey_departure_type[j] == 6)) && options->forced_flyby_coast > 1.0e-6)
            total_available_thrust_time -= options->forced_flyby_coast;
        if ((p < options->number_of_phases[j] - 1 || (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5)) && options->forced_flyby_coast > 1.0e-6)
            total_available_thrust_time -= options->forced_flyby_coast;

        for (int step = 0; step < options->num_timesteps; ++step)
        {
            time_step_sizes[step] *= total_available_thrust_time / step_size_normalization_coefficient;
        }

        //Step 6.2: Propagate across an initial coast if applicable
        double state_after_initial_coast[7];
        phase_time_elapsed_forward = 0.0;
        if (j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6)
        {
            //if this is a launch AND we are doing a forced post-launch initial coast
            Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_beginning_of_phase,
                                                        state_after_initial_coast,
                                                        Universe->mu,
                                                        Universe->LU,
                                                        options->forced_post_launch_coast,
                                                        this->initial_coast_Kepler_F,
                                                        this->initial_coast_Kepler_G,
                                                        this->initial_coast_Kepler_Fdot,
                                                        this->initial_coast_Kepler_Gdot,
                                                        this->initial_coast_Kepler_Fdotdot,
                                                        this->initial_coast_Kepler_Gdotdot,
                                                        this->initial_coast_STM,
                                                        (options->derivative_type > 1 && needG ? true : false));

            this->phase_time_elapsed_forward += options->forced_post_launch_coast;
            state_after_initial_coast[6] = this->state_at_beginning_of_phase[6];
        }
        else if ((p > 0 || p == 0 && (options->journey_departure_type[j] == 3
                                        || options->journey_departure_type[j] == 4
                                        || options->journey_departure_type[j] == 6))
                                    && options->forced_flyby_coast > 1.0e-6)
        {
            //if we are coming out of a flyby and we are doing a forced post-flyby coast
            Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_beginning_of_phase,
                                                        state_after_initial_coast,
                                                        Universe->mu,
                                                        Universe->LU,
                                                        options->forced_flyby_coast,
                                                        this->initial_coast_Kepler_F,
                                                        this->initial_coast_Kepler_G,
                                                        this->initial_coast_Kepler_Fdot,
                                                        this->initial_coast_Kepler_Gdot,
                                                        this->initial_coast_Kepler_Fdotdot,
                                                        this->initial_coast_Kepler_Gdotdot,
                                                        this->initial_coast_STM,
                                                        (options->derivative_type > 1 && needG ? true : false));

            this->phase_time_elapsed_forward += options->forced_flyby_coast;
            state_after_initial_coast[6] = this->state_at_beginning_of_phase[6];
        }

        //Step 6.3 propagate forward over the steps
        for (size_t step = 0; step < options->num_timesteps; ++step)
        {
            //Step 6.3.1 extract the left-hand state variables for this phase
            for (size_t state = 0; state < 7; ++state)
                this->left_hand_state[step][state] = X[*Xindex + state];
            (*Xindex) += 7;

            //Step 6.3.2 apply the radial distance constraint and its derivatives
            double r = math::norm(this->left_hand_state[step].data(), 3);
            F[*Findex] = (this->minimum_radial_distance - r) / this->minimum_radial_distance;
            ++(*Findex);

            if (needG && options->derivative_type > 0)
            {
                G[this->radius_constraint_G_indices[step][0]] = this->radius_constraint_X_scale_ranges[step][0] * (-this->left_hand_state[step][0] / r) / this->minimum_radial_distance;
                G[this->radius_constraint_G_indices[step][1]] = this->radius_constraint_X_scale_ranges[step][1] * (-this->left_hand_state[step][1] / r) / this->minimum_radial_distance;
                G[this->radius_constraint_G_indices[step][2]] = this->radius_constraint_X_scale_ranges[step][2] * (-this->left_hand_state[step][2] / r) / this->minimum_radial_distance;
            }

            //Step 6.3.3 apply the left-hand defect constraint and its derivatives
            //if this is the first step and there was an initial coast then we apply relative to the initial coast
            if (step == 0 && this->initial_coast)
            {
                for (size_t k = 0; k<3; ++k)
                {
                    //position
                    F[*Findex + k] = (this->left_hand_state[step][k] - state_after_initial_coast[k]) / Universe->LU;

                    //velocity
                    F[*Findex + k + 3] = (this->left_hand_state[step][k + 3] - state_after_initial_coast[k + 3]) / Universe->LU * Universe->TU;
                }
                //mass
                F[*Findex + 6] = (this->left_hand_state[step][6] - state_after_initial_coast[6]) / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
                (*Findex) += 7;
            }
            //if this is the first step and there was not an initial coast, apply relative to state_at_beginning_of_phase
            else if (step == 0)
            {
                for (size_t k = 0; k<3; ++k)
                {
                    //position
                    F[*Findex + k] = (this->left_hand_state[step][k] - this->state_at_beginning_of_phase[k]) / Universe->LU;

                    //velocity
                    F[*Findex + k + 3] = (this->left_hand_state[step][k + 3] - this->state_at_beginning_of_phase[k + 3]) / Universe->LU * Universe->TU;
                }
                //mass
                F[*Findex + 6] = (this->left_hand_state[step][6] - this->state_at_beginning_of_phase[6]) / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
                (*Findex) += 7;
            }
            //if this is NOT the first step then apply relative to the right-hand boundary of the previous step
            else
            {
                for (size_t k = 0; k<3; ++k)
                {
                    //position
                    F[*Findex + k] = (this->left_hand_state[step][k] - this->right_hand_state[step - 1][k]) / Universe->LU;

                    //velocity
                    F[*Findex + k + 3] = (this->left_hand_state[step][k + 3] - this->right_hand_state[step - 1][k + 3]) / Universe->LU * Universe->TU;
                }
                //mass
                F[*Findex + 6] = (this->left_hand_state[step][6] - this->right_hand_state[step - 1][6]) / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
                (*Findex) += 7;
            }

            //Step 6.3.4 propagate to the midpoint of the step
            Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->left_hand_state[step].data(),
                                                        this->spacecraft_state[step].data(),
                                                        Universe->mu,
                                                        Universe->LU,
                                                        this->time_step_sizes[step] / 2.0,
                                                        this->Kepler_F[2 * step],
                                                        this->Kepler_G[2 * step],
                                                        this->Kepler_Fdot[2 * step],
                                                        this->Kepler_Gdot[2 * step],
                                                        this->Kepler_Fdotdot[2 * step],
                                                        this->Kepler_Gdotdot[2 * step],
                                                        &(this->STM[2 * step]),
                                                        (options->derivative_type > 1 && needG ? true : false));
            this->spacecraft_state[step][6] = this->left_hand_state[step][6];
            this->phase_time_elapsed_forward += this->time_step_sizes[step] / 2.0;


            //Step 6.3.5 extract the control unit vector and, if applicable, the variable Isp for the step
            control[step][0] = X[*Xindex];
            control[step][1] = X[*Xindex + 1];
            control[step][2] = X[*Xindex + 2];
            (*Xindex) += 3;

            if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
            {
                available_Isp[step] = X[*Xindex];
                ++(*Xindex);
            }

            //Step 6.3.6 apply the control unit vector magnitude constraint and its derivatives
            //"mass leak" throttle of 1.0e-10 to ensure that we always have at least some thrusting (prevents derivative from being zero)
            local_throttle = math::norm(control[step].data(), 3) + 1.0e-10;
            throttle[step] = local_throttle;
            if (options->derivative_type > 1 && needG)
            {
                G[control_vector_G_indices[step][0]] = 2.0 * control[step][2] / local_throttle;
                G[control_vector_G_indices[step][1]] = 2.0 * control[step][1] / local_throttle;
                G[control_vector_G_indices[step][2]] = 2.0 * control[step][0] / local_throttle;
            }

            F[*Findex] = local_throttle;
            ++(*Findex);

            //Step 6.3.7 compute thruster performance for this step
            this->event_epochs[step] = phase_start_epoch + phase_time_elapsed_forward;
            EMTG::Astrodynamics::force_model(options,
                                            Universe,
                                            spacecraft_state[step].data(),
                                            &event_epochs[step],
                                            X,
                                            control[step].data(),
                                            &available_thrust[step],
                                            &available_mass_flow_rate[step],
                                            &available_Isp[step],
                                            &available_power[step],
                                            &active_power[step],
                                            &number_of_active_engines[step],
                                            ForceVector[step].data(),
                                            (options->derivative_type > 1 && needG ? 1 : 0),
                                            &dTdP[step],
                                            &dmdotdP[step],
                                            &dTdIsp[step],
                                            &dmdotdIsp[step],
                                            &dPdr[step],
                                            &dPdt[step],
                                            &dFSRPdr[step],
                                            dagravdRvec[step],
                                            dagravdtvec[step]);

            double effective_mass = spacecraft_state[step][6] > 1.0e-3 ? spacecraft_state[step][6] : 1.0e-3;
            dVmax[step] = options->engine_duty_cycle * available_thrust[step] / effective_mass * (time_step_sizes[step]);

            //Step 6.3.8 apply the control for this step
            for (size_t k = 0; k < 3; ++k)
            {
                dV[step][k] = ForceVector[step][k] / spacecraft_state[step][6] * (time_step_sizes[step]);
                spacecraft_state[step][k + 3] += dV[step][k];
            }
            impulse_magnitude = local_throttle * dVmax[step];
            *current_deltaV += impulse_magnitude;
            spacecraft_state[step][6] -= local_throttle * options->engine_duty_cycle * available_mass_flow_rate[step] * (time_step_sizes[step]);

            //Step 6.3.9 propagate to the right-hand side of the step
            Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->spacecraft_state[step].data(),
                                                        this->right_hand_state[step].data(),
                                                        Universe->mu,
                                                        Universe->LU,
                                                        this->time_step_sizes[step] / 2.0,
                                                        this->Kepler_F[2 * step + 1],
                                                        this->Kepler_G[2 * step + 1],
                                                        this->Kepler_Fdot[2 * step + 1],
                                                        this->Kepler_Gdot[2 * step + 1],
                                                        this->Kepler_Fdotdot[2 * step + 1],
                                                        this->Kepler_Gdotdot[2 * step + 1],
                                                        &(this->STM[2 * step + 1]),
                                                        (options->derivative_type > 1 && needG ? true : false));
            this->right_hand_state[step][6] = this->spacecraft_state[step][6];
            this->phase_time_elapsed_forward += this->time_step_sizes[step] / 2.0;
        }

        //Step 6.4 apply the terminal coast if applicable
        double state_before_terminal_coast[7];
        if ((p < options->number_of_phases[j] - 1
            || (options->journey_arrival_type[j] == 2
            || options->journey_arrival_type[j] == 5)) && options->forced_flyby_coast > 1.0e-6)
        {
            //if we are going into a flyby or intercept
            Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_end_of_phase,
                                                        state_before_terminal_coast,
                                                        Universe->mu,
                                                        Universe->LU,
                                                        -options->forced_flyby_coast,
                                                        this->terminal_coast_Kepler_F,
                                                        this->terminal_coast_Kepler_G,
                                                        this->terminal_coast_Kepler_Fdot,
                                                        this->terminal_coast_Kepler_Gdot,
                                                        this->terminal_coast_Kepler_Fdotdot,
                                                        this->terminal_coast_Kepler_Gdotdot,
                                                        this->terminal_coast_STM,
                                                        (options->derivative_type > 1 && needG ? true : false));

            this->phase_time_elapsed_backward += options->forced_flyby_coast;
            state_before_terminal_coast[6] = this->state_at_end_of_phase[6];
        }

        //Step 6.5 apply the right-hand defect constraint
        //this is applied with respect to the pre-terminal coast state if a coast occured
        
        if (terminal_coast)
        {
            for (size_t k = 0; k<3; ++k)
            {
                //position
                F[*Findex + k] = (this->right_hand_state[options->num_timesteps - 1][k] - state_before_terminal_coast[k]) / Universe->LU;

                //velocity
                F[*Findex + k + 3] = (this->right_hand_state[options->num_timesteps - 1][k + 3] - state_before_terminal_coast[k + 3]) / Universe->LU * Universe->TU;
            }
            //mass
            F[*Findex + 6] = (this->right_hand_state[options->num_timesteps - 1][6] - state_before_terminal_coast[6]) / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
            (*Findex) += 7;
        }
        //otherwise it is applied to the state at the end of the phase
        else 
        {
            for (size_t k = 0; k<3; ++k)
            {
                //position
                F[*Findex + k] = (this->right_hand_state[options->num_timesteps - 1][k] - this->state_at_end_of_phase[k]) / Universe->LU;

                //velocity
                F[*Findex + k + 3] = (this->right_hand_state[options->num_timesteps - 1][k + 3] - this->state_at_end_of_phase[k + 3]) / Universe->LU * Universe->TU;
            }
            //mass
            F[*Findex + 6] = (this->right_hand_state[options->num_timesteps - 1][6] - this->state_at_end_of_phase[6]) / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
            (*Findex) += 7;
        }

        //Step 6.6 calculate defect constraint derivatives
        if (needG && options->derivative_type > 1)
            this->calculate_defect_derivatives(G, j, p, options, Universe);

        //******************************************************************
        //Step 7: process the arrival, if applicable
        if (p == options->number_of_phases[j] - 1)
        {
            if (options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 5)
                dV_arrival_magnitude = 0.0;

            //note that "3" for journey_arrival_type indicates a "low-thrust rendezvous," which means we are there already and we don't need to do anything
            else
            {
                //compute the arrival deltaV
                if (boundary2_location_code > 0) //ending at body
                    dV_arrival_magnitude = process_arrival(state_at_end_of_phase + 3,
                    boundary2_state,
                    current_state + 3,
                    current_epoch,
                    Body2->mu,
                    Body2->r_SOI,
                    F,
                    Findex,
                    j,
                    options,
                    Universe);
                else //ending at point on central body SOI, fixed point, or fixed orbit
                    dV_arrival_magnitude = process_arrival(state_at_end_of_phase + 3,
                    boundary2_state,
                    current_state + 3,
                    current_epoch,
                    Universe->mu,
                    Universe->r_SOI,
                    F,
                    Findex,
                    j,
                    options,
                    Universe);

                //drop the electric propulsion stage
                state_at_end_of_phase[6] -= options->EP_dry_mass;

                //apply the arrival burn
                state_at_end_of_phase[6] *= exp(-dV_arrival_magnitude * 1000 / ((options->IspChem > 0 ? options->IspChem : 1e-6)* options->g0));
                *current_deltaV += dV_arrival_magnitude;

                //if this is a terminal intercept, we need to compute derivatives
                if (needG && options->journey_arrival_type[j] == 2)
                {
                    double C3_desired = options->journey_final_velocity[j][1] * options->journey_final_velocity[j][1];
                    //derivative with respect to x component of terminal velocity
                    G[terminal_velocity_constraint_G_indices[0]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[0] * X[terminal_velocity_constraint_X_indices[0]] / C3_desired;

                    //derivative with respect to y component of terminal velocity
                    G[terminal_velocity_constraint_G_indices[1]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[1] * X[terminal_velocity_constraint_X_indices[1]] / C3_desired;

                    //derivative with respect to z component of terminal velocity
                    G[terminal_velocity_constraint_G_indices[2]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[2] * X[terminal_velocity_constraint_X_indices[2]] / C3_desired;
                }
                else if (needG && options->journey_arrival_type[j] == 6)
                {
                    double r = math::norm(boundary2_state, 3) / Universe->LU;
                    double denominator = r * r * r * Universe->LU * Universe->LU;
                    double LUTU2 = (Universe->TU / Universe->LU) * (Universe->TU / Universe->LU);
                    double sqLU = sqrt(Universe->LU);

                    //derivative with respect to x component of position
                    G[escape_constraint_G_indices[0]] = escape_constraint_X_scale_ranges[0] * X[escape_constraint_X_indices[0]] / denominator;

                    //derivative with respect to x component of position
                    G[escape_constraint_G_indices[1]] = escape_constraint_X_scale_ranges[1] * X[escape_constraint_X_indices[1]] / denominator;

                    //derivative with respect to x component of position
                    G[escape_constraint_G_indices[2]] = escape_constraint_X_scale_ranges[2] * X[escape_constraint_X_indices[2]] / denominator;

                    //derivative with respect to x component of terminal velocity
                    G[escape_constraint_G_indices[3]] = escape_constraint_X_scale_ranges[3] * X[escape_constraint_X_indices[3]] * LUTU2;

                    //derivative with respect to y component of terminal velocity
                    G[escape_constraint_G_indices[4]] = escape_constraint_X_scale_ranges[4] * X[escape_constraint_X_indices[4]] * LUTU2;

                    //derivative with respect to z component of terminal velocity
                    G[escape_constraint_G_indices[5]] = escape_constraint_X_scale_ranges[5] * X[escape_constraint_X_indices[5]] * LUTU2;
                }
            }

        }

        //******************************************************************
        //Step 8: update the current epoch
        *current_epoch += this->TOF;

        //******************************************************************
        //Step 9: update the current state
        if (p == (options->number_of_phases[j] - 1) && options->journey_arrival_type[j] == 7)
        {//for terminal phases ending in a capture spiral
            for (int k = 0; k < 7; ++k)
                current_state[k] = this->spiral_capture_state_after_spiral[k];
            *current_epoch += this->spiral_capture_time;
        }
        else
        {
            for (int k = 0; k < 7; ++k)
                current_state[k] = this->state_at_end_of_phase[k];
        }

        return errcode;
    }

    //output function
    //return 0 if successful, 1 if failure
    int PSBIphase::output(missionoptions* options,
        const double& launchdate,
        int j,
        int p,
        EMTG::Astrodynamics::universe* Universe,
        int* eventcount)
    {
        //Step 1: store data that will be used for the printing
        double empty_vector[] = { 0, 0, 0 };
        double phase_time_elapsed = 0.0;
        string event_type;
        double temp_power, temp_thrust, temp_mdot, temp_Isp,
            temp_active_power, temp_dTdP, temp_dmdotdP,
            temp_dTdIsp, temp_dmdotdIsp, temp_dPdr, temp_dPdt;
        int temp_active_thrusters;

        if (p > 0 || (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4))
        {
            event_type = "upwr_flyby";
            math::Matrix<double> periapse_state = calculate_flyby_periapse_state(V_infinity_in, V_infinity_out, flyby_altitude, *Body1);
            math::Matrix<double> periapse_R(3, 1);
            periapse_R(0) = periapse_state(0);
            periapse_R(1) = periapse_state(1);
            periapse_R(2) = periapse_state(2);
            Bplane.define_bplane(V_infinity_in, BoundaryR, BoundaryV);
            Bplane.compute_BdotR_BdotT_from_periapse_position(Body1->mu, V_infinity_in, periapse_R, &BdotR, &BdotT);

            //compute RA and DEC in the frame of the target body
            this->Body1->J2000_body_equatorial_frame.construct_rotation_matrices(this->phase_start_epoch / 86400.0 + 2400000.5);
            math::Matrix<double> rot_out_vec = this->Body1->J2000_body_equatorial_frame.R_from_ICRF_to_local * V_infinity_in;

            this->RA_departure = atan2(rot_out_vec(1), rot_out_vec(0));

            this->DEC_departure = asin(rot_out_vec(2) / V_infinity_in.norm());
        }
        else if (j == 0 && boundary1_location_code > 0 && options->LV_type >= 0)
            event_type = "launch";
        else if (options->journey_departure_type[j] == 6)
        {
            event_type = "zeroflyby";
            //compute RA and DEC in the frame of the target body
            this->Body1->J2000_body_equatorial_frame.construct_rotation_matrices(this->phase_start_epoch / 86400.0 + 2400000.5);
            math::Matrix<double> rot_out_vec = this->Body1->J2000_body_equatorial_frame.R_from_ICRF_to_local * V_infinity_in;

            this->RA_departure = atan2(rot_out_vec(1), rot_out_vec(0));

            this->DEC_departure = asin(rot_out_vec(2) / V_infinity_in.norm());
        }
        else
            event_type = "departure";

        string boundary1_name;
        string boundary2_name;

        switch (boundary1_location_code)
        {
        case -1:
        {
            boundary2_name = "free point";
            break;
        }
        case -2: //begin at SOI
        {
            boundary1_name = "Hyp-arrival";
            break;
        }
        default:
            boundary1_name = (Universe->bodies[boundary1_location_code - 1].name);
        }

        switch (boundary2_location_code)
        {
        case -1:
        {
            boundary2_name = "free point";
            break;
        }
        default:
            boundary2_name = (Universe->bodies[boundary2_location_code - 1].name);
        }
        double initial_Isp;
        if (j == 0 && p == 0)
        {
            if (options->journey_departure_type[j] == 1 || options->LV_type == -1)
            {
                initial_Isp = options->IspDS;
            }
            else
            {
                initial_Isp = -1;
            }
        }
        else
        {
            initial_Isp = options->IspChem;
        }

        //Step 2: all phases have at least 2+n events: departure, arrival, and n thrust/coast
        //*****************************************************************************
        //first let's print the departure/flyby

        //Step 2.01 if this phase starts with an escape spiral, we have to encode the spiral start date
        if (options->journey_departure_type[j] == 5 && p == 0)
        {
            event_type = "departure";

            this->write_summary_line(options,
                Universe,
                eventcount,
                (this->phase_start_epoch - this->spiral_escape_time) / 86400.0,
                "begin_spiral",
                this->Body1->name,
                this->spiral_escape_time / 86400.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                this->spiral_escape_state_before_spiral,
                empty_vector,
                empty_vector,
                this->spiral_escape_dv,
                this->spiral_escape_thrust * 1000.0, //kN to N conversion
                this->spiral_escape_Isp,
                this->spiral_escape_power,
                this->spiral_escape_mdot,
                this->spiral_escape_number_of_engines,
                this->spiral_escape_active_power);

        }

        //then print the departure itself
        for (int k = 0; k < 3; ++k)
            dVdeparture[k] = V_infinity_out(k);

        //we have to calculate the available power at departure
        Astrodynamics::find_engine_parameters(options,
            math::norm(this->state_at_beginning_of_phase, 3) / Universe->LU,
            (this->phase_start_epoch),
            &temp_thrust,
            &temp_mdot,
            &temp_Isp,
            &temp_power,
            &temp_active_power,
            &temp_active_thrusters,
            false,
            &temp_dTdP,
            &temp_dmdotdP,
            &temp_dTdIsp,
            &temp_dmdotdIsp,
            &temp_dPdr,
            &temp_dPdt);

        write_summary_line(options,
            Universe,
            eventcount,
            this->phase_start_epoch / 86400.0,
            event_type,
            boundary1_name,
            0,
            (p > 0 ? flyby_altitude : Bradius),
            this->Btheta,
            (p > 0 ? this->flyby_turn_angle : -1),
            this->RA_departure,
            this->DEC_departure,
            this->C3_departure,
            this->state_at_beginning_of_phase,
            this->dVdeparture,
            empty_vector,
            (p == 0 ? (options->journey_departure_type[j] == 5 ? 0.0 : this->dV_departure_magnitude) : this->flyby_outgoing_v_infinity),
            -1,
            initial_Isp,
            temp_power,
            0,
            0,
            0);

        //*****************************************************************************
        //next, if there was an initial coast, we must print it
        //we'll have to start by propagating to the halfway point of the initial coast step
        if ((j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6) || ((p > 0 || p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4 || options->journey_departure_type[j] == 6)) && options->forced_flyby_coast > 1.0e-6))
        {
            double state_at_initial_coast_midpoint[7];
            double initial_coast_duration;
            state_at_initial_coast_midpoint[6] = state_at_beginning_of_phase[6];

            if (j == 0 && p == 0 && options->forced_post_launch_coast > 0.0)
            {
                //initial coast after launch
                initial_coast_duration = options->forced_post_launch_coast;
            }
            else
            {
                //initial coast after flyby
                initial_coast_duration = options->forced_flyby_coast;
            }

            Kepler::Kepler_Lagrange_Laguerre_Conway_Der(state_at_beginning_of_phase,
                state_at_initial_coast_midpoint,
                Universe->mu,
                Universe->LU,
                initial_coast_duration / 2.0);

            //we have to calculate the available power at the midpoint of the forced coast
            Astrodynamics::find_engine_parameters(options,
                math::norm(state_at_initial_coast_midpoint, 3) / Universe->LU,
                phase_start_epoch + initial_coast_duration / 2.0,
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                false,
                &temp_dTdP,
                &temp_dmdotdP,
                &temp_dTdIsp,
                &temp_dmdotdIsp,
                &temp_dPdr,
                &temp_dPdt);

            write_summary_line(options,
                Universe,
                eventcount,
                (phase_start_epoch + 0.5 * initial_coast_duration) / 86400.0,
                "force-coast",
                "deep-space",
                initial_coast_duration / 86400.0,
                -1,
                -1,
                -1,
                0.0,
                0.0,
                0.0,
                state_at_initial_coast_midpoint,
                empty_vector,
                empty_vector,
                0.0,
                0.0,
                0.0,
                temp_power,
                0.0,
                0,
                0.0);

            phase_time_elapsed += initial_coast_duration;
        }

        //*****************************************************************************
        //next, we must print each thrust arc


        for (int step = 0; step < options->num_timesteps; ++step)
        {
            double thrust_vector[3];
            double impulse_magnitude = sqrt(dV[step][0] * dV[step][0] + dV[step][1] * dV[step][1] + dV[step][2] * dV[step][2]);
            double throttle_magnitude = sqrt(control[step][0] * control[step][0] + control[step][1] * control[step][1] + control[step][2] * control[step][2]);
            double current_mass_flow_rate;

            double current_Isp, current_thrust, current_power;
            if (options->engine_type == 0) //fixed thrust/Isp
            {
                current_thrust = options->Thrust * 1000.0;
                current_Isp = options->IspLT;
                current_power = -1;
                current_mass_flow_rate = current_thrust / current_Isp / options->g0;
            }
            else //thrust, Isp are functions of power
            {
                current_thrust = available_thrust[step] * 1000.0; //kN to N conversion
                current_Isp = available_Isp[step];
                current_power = available_power[step];
                current_mass_flow_rate = available_mass_flow_rate[step];// current_thrust / current_Isp / options->g0;
            }

            if (EMTG::math::norm(control[step].data(), 3) > 1.0e-2 && fabs(current_thrust) > 1.0e-6)
            {
                event_type = "PSBIthrust";
            }
            else
            {
                event_type = "coast";
            }

            for (int k = 0; k < 3; ++k)
                thrust_vector[k] = control[step][k] * options->engine_duty_cycle * current_thrust;

            double state[7];
            state[6] = spacecraft_state[step][6];
            for (int k = 0; k < 6; ++k)
                state[k] = spacecraft_state[step][k];
            if (step < options->num_timesteps / 2)
            {
                for (int k = 0; k < 3; ++k)
                    state[k + 3] -= dV[step][k];
            }

            write_summary_line(options,
                Universe,
                eventcount,
                (this->phase_start_epoch + phase_time_elapsed + 0.5 * this->time_step_sizes[step]) / 86400.0,
                event_type,
                "deep-space",
                (this->time_step_sizes[step]) / 86400.0,
                -1,
                -1,
                -1,
                atan2(this->dV[step][1] + math::SMALL, this->dV[step][0]) * EMTG::math::PI / 180.0,
                asin(this->dV[step][2] / impulse_magnitude) * EMTG::math::PI / 180.0,
                0,
                state,
                this->dV[step].data(),
                thrust_vector,
                impulse_magnitude,
                current_thrust,
                current_Isp,
                current_power,
                math::norm(this->control[step].data(), 3) * current_mass_flow_rate,
                this->number_of_active_engines[step],
                this->active_power[step]);

            phase_time_elapsed += this->time_step_sizes[step];
        }

        //*****************************************************************************
        //next, if there was an terminal coast, we must print it
        //we'll have to start by propagating to the halfway point of the initial coast step
        if ((p < options->number_of_phases[j] - 1 || (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5)) && options->forced_flyby_coast > 1.0e-6)
        {
            double state_at_terminal_coast_midpoint[7];
            double terminal_coast_duration = options->forced_flyby_coast;
            state_at_terminal_coast_midpoint[6] = state_at_end_of_phase[6];

            Kepler::Kepler_Lagrange_Laguerre_Conway_Der(state_at_end_of_phase,
                state_at_terminal_coast_midpoint,
                Universe->mu,
                Universe->LU,
                -terminal_coast_duration / 2.0);

            //we have to calculate the available power at the midpoint of the forced coast
            Astrodynamics::find_engine_parameters(options,
                math::norm(state_at_terminal_coast_midpoint, 3) / Universe->LU,
                this->phase_start_epoch + phase_time_elapsed + 0.5 * terminal_coast_duration,
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                false,
                &temp_dTdP,
                &temp_dmdotdP,
                &temp_dTdIsp,
                &temp_dmdotdIsp,
                &temp_dPdr,
                &temp_dPdt);

            write_summary_line(options,
                Universe,
                eventcount,
                (this->phase_start_epoch + phase_time_elapsed + 0.5 * terminal_coast_duration) / 86400.0,
                "force-coast",
                "deep-space",
                terminal_coast_duration / 86400.0,
                -1,
                -1,
                -1,
                0.0,
                0.0,
                0.0,
                state_at_terminal_coast_midpoint,
                empty_vector,
                empty_vector,
                0.0,
                0.0,
                0.0,
                temp_power,
                0.0,
                0,
                0.0);

            phase_time_elapsed += terminal_coast_duration;
        }

        //*****************************************************************************
        //finally, terminal phases have an arrival maneuver

        if (p == options->number_of_phases[j] - 1)
        {
            if (options->journey_arrival_type[j] == 0)
                event_type = "insertion";
            else if (options->journey_arrival_type[j] == 1)
                event_type = "rendezvous";
            else if (options->journey_arrival_type[j] == 2)
                event_type = "intercept";
            else if (options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 7)
                event_type = "LT_rndzvs";
            else if (options->journey_arrival_type[j] == 5 || options->journey_arrival_type[j] == 4)
                event_type = "match-vinf";

            //compute RA and DEC in the frame of the target body
            if (options->destination_list[j][1] > 0)
            {
                this->Body2->J2000_body_equatorial_frame.construct_rotation_matrices((this->phase_start_epoch + this->TOF) / 86400.0 + 2400000.5);
                math::Matrix<double> rot_in_vec(3, 1, this->dVarrival);
                math::Matrix<double> rot_out_vec = this->Body2->J2000_body_equatorial_frame.R_from_ICRF_to_local * rot_in_vec;

                this->RA_arrival = atan2(rot_out_vec(1), rot_out_vec(0));

                this->DEC_arrival = asin(rot_out_vec(2) / sqrt(this->C3_arrival));
            }
            else
            {
                this->RA_arrival = 0.0;
                this->DEC_arrival = 0.0;
            }


            double dV_arrival_mag;
            if (options->journey_arrival_type[j] == 2)
            {
                dV_arrival_mag = sqrt(this->C3_arrival);
            }
            else if (options->journey_arrival_type[j] == 4 || options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 7)
            {
                dV_arrival_mag = 0;
                this->dVarrival[0] = 0;
                this->dVarrival[1] = 0;
                this->dVarrival[2] = 0;

                this->RA_arrival = 0.0;
                this->DEC_arrival = 0.0;
            }
            else
            {
                dV_arrival_mag = this->dV_arrival_magnitude;
            }

            //we have to calculate the available power at boundary
            Astrodynamics::find_engine_parameters(options,
                math::norm(this->state_at_end_of_phase, 3) / Universe->LU,
                (this->phase_start_epoch + this->TOF),
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                false,
                &temp_dTdP,
                &temp_dmdotdP,
                &temp_dTdIsp,
                &temp_dmdotdIsp,
                &temp_dPdr,
                &temp_dPdt);

            write_summary_line(options,
                Universe,
                eventcount,
                (this->phase_start_epoch + this->TOF) / 86400.0,
                event_type,
                boundary2_name,
                0,
                -1,
                -1,
                -1,
                this->RA_arrival,
                this->DEC_arrival,
                this->C3_arrival,
                this->state_at_end_of_phase,
                this->dVarrival,
                empty_vector,
                dV_arrival_mag,
                -1,
                options->IspChem,
                temp_power,
                0,
                0,
                0);

            //if the phase ends in a capture spiral, print it
            if (options->journey_arrival_type[j] == 7)
            {
                this->write_summary_line(options,
                    Universe,
                    eventcount,
                    (this->phase_start_epoch + this->TOF + this->spiral_capture_time) / 86400.0,
                    "end_spiral",
                    this->Body2->name,
                    this->spiral_capture_time / 84600.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    this->spiral_capture_state_after_spiral,
                    empty_vector,
                    empty_vector,
                    this->spiral_capture_dv,
                    this->spiral_capture_thrust * 1000.0, //kN to N conversion
                    this->spiral_capture_Isp,
                    this->spiral_capture_power,
                    this->spiral_capture_mdot,
                    this->spiral_capture_number_of_engines,
                    this->spiral_capture_active_power);
            }
        }

        return 0;
    }

    //method to calculate defect constraint derivatives
    void PSBIphase::calculate_defect_derivatives(double* G,
                                                const int& j,
                                                const int& p,
                                                missionoptions* options,
                                                EMTG::Astrodynamics::universe* Universe)
    {
        //declare variables which will be generally useful
        double dxdu, dydu, dzdu, dxdotdu, dydotdu, dzdotdu, dmdu, dtdu;

        //derivatives for the left-hand defect constraints
        //each constraint has a left side and a right side, and derivatives are posed as:
        //dC/du = dRHS/du - dLHS/du
        //RHS is ALWAYS a decision variable and if dRHS/du is nonzero then dLHS/du is zero and vice versa
        //therefore it is practical to start with the RHS derivatives
        for (size_t step = 0; step < options->num_timesteps; ++step)
        {
            //derivative with respect to current state variable
            //RHS is defined, LHS is zero
            G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state[step][0]] = this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state[step][0] / Universe->LU;
            G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state[step][1]] = this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state[step][1] / Universe->LU;
            G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state[step][2]] = this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state[step][2] / Universe->LU;
            G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state[step][3]] = this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state[step][3] / Universe->LU * Universe->TU;
            G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state[step][4]] = this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state[step][4] / Universe->LU * Universe->TU;
            G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state[step][5]] = this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state[step][5] / Universe->LU * Universe->TU;
            G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_current_state[step][6]] = this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_current_state[step][6] / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);

            //all of the other derivatives are of the LHS, which requires propagating an STM
            //there are different cases:
            //1. first step, initial coast present
            //2. first step, no initial coast
            //3. successive step (have to propagate STMs from left hand side of previous step to right hand side of current step)
            if (step == 0)
            {
                //if this journey starts with a launch then the initial v-infinity is given in polar coordinates
                if (p == 0 && options->journey_departure_type[j] == 0)
                {
                    double vinf = sqrt(this->C3_departure);
                    double cosRA = cos(this->RA_departure);
                    double sinRA = sin(this->RA_departure);
                    double cosDEC = cos(this->DEC_departure);
                    double sinDEC = sin(this->DEC_departure);
                    double dVxdVinf = cosRA * cosDEC;
                    double dVydVinf = sinRA * cosDEC;
                    double dVzdVinf = sinDEC;
                    double dVxdRA = vinf * (-sinRA * cosDEC);
                    double dVydRA = vinf * (cosRA * cosDEC);
                    double dVzdRA = 0.0;
                    double dVxdDEC = vinf * (-cosRA * sinDEC);
                    double dVydDEC = vinf * (-sinRA * sinDEC);
                    double dVzdDEC = vinf * cosDEC;

                    if (initial_coast)
                    {
                        //derivatives with respect to vinf
                        for (size_t state = 0; state < 3; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][0]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][0]
                            * ((*this->initial_coast_STM)(state, 3) * dVxdVinf
                            + (*this->initial_coast_STM)(state, 4) * dVydVinf
                            + (*this->initial_coast_STM)(state, 5) * dVzdVinf) / Universe->LU;
                        for (size_t state = 3; state < 6; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][0]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][0]
                            * ((*this->initial_coast_STM)(state, 3) * dVxdVinf
                            + (*this->initial_coast_STM)(state, 4) * dVydVinf
                            + (*this->initial_coast_STM)(state, 5) * dVzdVinf) / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[6][0]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[6][0] * this->dmdvinf / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
                        //derivatives with respect to RA
                        for (size_t state = 0; state < 3; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][1]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][1]
                            * ((*this->initial_coast_STM)(state, 3) * dVxdRA
                            + (*this->initial_coast_STM)(state, 4) * dVydRA
                            + (*this->initial_coast_STM)(state, 5) * dVzdRA) / Universe->LU;
                        for (size_t state = 3; state < 6; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][1]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][1]
                            * ((*this->initial_coast_STM)(state, 3) * dVxdRA
                            + (*this->initial_coast_STM)(state, 4) * dVydRA
                            + (*this->initial_coast_STM)(state, 5) * dVzdRA) / Universe->LU * Universe->TU;
                        //derivatives with respect to DEC
                        for (size_t state = 0; state < 3; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][2]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][2]
                            * ((*this->initial_coast_STM)(state, 3) * dVxdDEC
                            + (*this->initial_coast_STM)(state, 4) * dVydDEC
                            + (*this->initial_coast_STM)(state, 5) * dVzdDEC) / Universe->LU;
                        for (size_t state = 3; state < 6; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][2]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][2]
                            * ((*this->initial_coast_STM)(state, 3) * dVxdDEC
                            + (*this->initial_coast_STM)(state, 4) * dVydDEC
                            + (*this->initial_coast_STM)(state, 5) * dVzdDEC) / Universe->LU * Universe->TU;
                    }
                    else //if no initial coast
                    {
                        //no dependence of position on initial velocity increment
                        //derivatives with respect to vinf
                        for (size_t state = 0; state < 3; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][0]] = 0.0;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[3][0]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[3][0] * dVxdVinf / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[4][0]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[4][0] * dVydVinf / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[5][0]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[5][0] * dVzdVinf / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[6][0]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[6][0] * this->dmdvinf / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
                        //derivatives with respect to RA
                        for (size_t state = 0; state < 3; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][1]] = 0.0;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[3][1]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[3][1] * dVxdRA / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[4][1]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[4][1] * dVydRA / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[5][1]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[5][1] * dVzdRA / Universe->LU * Universe->TU;
                        //derivatives with respect to DEC
                        for (size_t state = 0; state < 3; ++state)
                            G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[state][2]] = 0.0;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[3][2]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[3][2] * dVxdDEC / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[4][2]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[4][2] * dVydDEC / Universe->LU * Universe->TU;
                        G[this->G_index_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[5][2]] = -this->X_scale_range_of_derivative_of_leftmost_defect_constraints_with_respect_to_phase_initial_velocity[5][2] * dVzdDEC / Universe->LU * Universe->TU;

                    }
                    
                }
                //otherwise if this phase does not start with a launch then the initial v-infinity is given in cartesian coordinates
                //(i.e. flyby v-infinity out)
                //and the calculation is much simpler
                else if (p > 0 || (p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 6)))
                {

                }
            }
            else //for successive steps
            {
                //TODO derivative with respect to previous step state variables
               
                //derivative with respect to previous step control variables
                double umag = math::norm(control[step - 1].data(), 3);
                double deltat = time_step_sizes[step];

                for (size_t c = 0; c < 3; ++c)
                {
                    //first we need the derivative of the right hand side of the previous step with respect to the control in that step
                    dxdu = this->STM[2 * step - 1](0, c + 3) * dVmax[step - 1];
                    dydu = this->STM[2 * step - 1](1, c + 3) * dVmax[step - 1];
                    dzdu = this->STM[2 * step - 1](2, c + 3) * dVmax[step - 1];
                    dxdotdu = this->STM[2 * step - 1](3, c + 3) * dVmax[step - 1];
                    dydotdu = this->STM[2 * step - 1](4, c + 3) * dVmax[step - 1];
                    dzdotdu = this->STM[2 * step - 1](5, c + 3) * dVmax[step - 1];
                    dmdu = -(available_mass_flow_rate[step - 1] * deltat * options->engine_duty_cycle) * ((control[step - 1][c] / (umag + 1.0e-10)));

                    //now we can fill that into the constraint equation
                    G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][0][7 + c]] = -this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][0][7 + c] * dxdu / Universe->LU;
                    G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][1][7 + c]] = -this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][1][7 + c] * dydu / Universe->LU;
                    G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][2][7 + c]] = -this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][2][7 + c] * dzdu / Universe->LU;
                    G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][3][7 + c]] = -this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][3][7 + c] * dxdotdu / Universe->LU * Universe->TU;
                    G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][4][7 + c]] = -this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][4][7 + c] * dydotdu / Universe->LU * Universe->TU;
                    G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][5][7 + c]] = -this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][5][7 + c] * dzdotdu / Universe->LU * Universe->TU;
                    G[this->G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][6][7 + c]] = -this->X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control[step][6][7 + c] * dmdu / (options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
                }
            }
        }

        //derivatives for the right-hand defect constraint
    }

}//close namespace EMTG