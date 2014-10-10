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
            this->spacecraft_state.push_back(state_dummy);
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
        this->STM.resize(options->num_timesteps + 1);
        this->Kepler_F.resize(options->num_timesteps + 1);
        this->Kepler_Fdot.resize(options->num_timesteps + 1);
        this->Kepler_G.resize(options->num_timesteps + 1);
        this->Kepler_Gdot.resize(options->num_timesteps + 1);
        this->Kepler_Fdotdot.resize(options->num_timesteps + 1);
        this->Kepler_Gdotdot.resize(options->num_timesteps + 1);
        this->Propagation_Step_Time_Fraction.resize(options->num_timesteps + 1);
        this->Propagation_Step_Time_Fraction.resize(options->num_timesteps + 1);

        this->current_mass_increment = 0.0;
        this->journey_initial_mass_increment_scale_factor = 1.0;

        //size the time step vector
        time_step_sizes.resize(options->num_timesteps + 1);

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
        //destructor does not have to do anything
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
            Xdescriptions->push_back(prefix + "arrival mass");

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
            Flowerbounds->push_back(-math::LARGE);
            Fupperbounds->push_back(0.0);
            Fdescriptions->push_back(prefix + " step " + stepstream.str() + " radial defect constraint");
            //radial distance constraint is dependent ONLY on x, y, z for THIS step
            vector<int> step_radius_constraint_G_indices;
            vector<double> step_radius_constraint_X_scale_ranges;
            for (size_t entry = Xdescriptions->size() - 1 - state_and_control_elements_this_step; entry < Xdescriptions->size() - 1 - state_and_control_elements_this_step + 3; ++entry)
            {
                iGfun->push_back(Fdescriptions->size() - 1);
                jGvar->push_back(entry);
                stringstream EntryNameStream;
                EntryNameStream << "Derivative of " << prefix + "step " << step << " radial distance constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                Gdescriptions->push_back(EntryNameStream.str());
                step_radius_constraint_G_indices.push_back(Gdescriptions->size() - 1);
                step_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
            }
            this->radus_constraint_G_indicies.push_back(step_radius_constraint_G_indices);
            this->radius_constraint_X_scale_ranges.push_back(step_radius_constraint_X_scale_ranges);


            //left-hand defect constraints: all steps have this!
            vector<string> statename;
            statename.push_back("x");
            statename.push_back("y");
            statename.push_back("z");
            statename.push_back("xdot");
            statename.push_back("ydot");
            statename.push_back("zdot");
            statename.push_back("m");

            //temporary arrays for G index tracking
            
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
                //  ALL previous time variables including the current phase flight time
                //  note that the FIRST entry in the list of time derivatives is with respect to the CURRENT phase flight time
                vector<int> state_G_index_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
                vector<double> state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_flight_time_variables;
                for (size_t entry = Xdescriptions->size() - 1; entry >= 0; --entry)
                {
                    if ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
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
                                EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
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

                //for the first defect only
                //  if applicable, variable mission initial mass
                //  if applicable, variable initial mass increment scale factor
                //  if applicable, variable left hand boundary condition
                //  if applicable, phase v-infinity (which for the first phase is in polar, successive phases in cartesian)
                //  previous phase arrival mass, if this is NOT the first phase in the mission
                // 
                if (step == 0)
                {
                    if (p == 0)
                    {

                    }
                }
                //for successive defects only
                //  the previous control point's state and control
                else
                {
                    vector<int> state_G_index_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;
                    vector<double> state_X_scale_range_of_derivative_of_defect_constraints_with_respect_to_previous_state_and_control;
                    for (size_t entry = Xdescriptions->size() - 1 - 2 * state_and_control_elements_this_step; entry < Xdescriptions->size() - 1 - state_and_control_elements_this_step; ++entry)
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
                step_G_indices.push_back(iGfun->size() - 1);
            }
            control_vector_G_indices.push_back(step_G_indices);
        }

        //**************************************************************************
        //right hand side defect constraint - the last segment only must connect to the right-hand boundary condition for the phase
        //has derivatives with respect to:
        //  rightmost step state and control variables
        //  flight time (all previous, first entry is CURRENT phase flight time)
        //  arrival mass
        //  final v-infinity vector
        //  variable right hand boundary condition
        //  if applicable, spiral things affect every defect constraint
        this->find_dependencies_due_to_escape_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
        this->find_dependencies_due_to_capture_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);

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

        return errcode;
    }

    //output function
    //return 0 if successful, 1 if failure
    int PSBIphase::output(  missionoptions* options,
                            const double& launchdate,
                            int j,
                            int p,
                            EMTG::Astrodynamics::universe* Universe,
                            int* eventcount)
    {
        return 0;
    }


}//close namespace EMTG