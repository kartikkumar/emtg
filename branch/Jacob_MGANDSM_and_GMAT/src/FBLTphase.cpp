/*
 * FBLTphase.cpp
 *
 *  Created on: September 17, 2012
 *      Author: Jacob
 */
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "FBLTphase.h"
#include "Astrodynamics.h"
#include "missionoptions.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "equations_of_motion.h"
#include "universe.h"
#include "EMTG_string_utilities.h"

#include "SpiceUsr.h"

//#define ENABLE_FBLT_INSTRUMENTATION



namespace EMTG 
{
    FBLT_phase::FBLT_phase() 
    {
    //default constructor does nothing

    }

    FBLT_phase::FBLT_phase(const int& j, const int& p, const missionoptions& options)
    {
        //call phase initialize method
        this->initialize(j, p, options);
        
	    //must resize all data vectors to the correct length
	    std::vector<double> vector7_dummy(7, 0.0);
        std::vector<double> vector3_dummy(3, 0.0);
        this->match_point_state.resize(7);

		//state vector containers
		this->spacecraft_state_forward = std::vector <double> (11 + 11 * 11, 0.0);
		this->spacecraft_state_forward_prop = std::vector <double>(11 + 11 * 11, 0.0);
		this->spacecraft_state_end_coast = std::vector <double>(11 + 11 * 11, 0.0);
		this->spacecraft_state_backward = std::vector <double>(11 + 11 * 11, 0.0);
		this->spacecraft_state_backward_prop = std::vector <double>(11 + 11 * 11, 0.0);
		this->spacecraft_state_propagate = std::vector <double>(11 + 11 * 11, 0.0);
		this->spacecraft_state_propagate_next = std::vector <double>(11 + 11 * 11, 0.0);
		this->match_point_state = std::vector <double>(7, 0.0);

		//phase TOF derivative containers
		this->dspacecraft_state_forwarddTOF = EMTG::math::Matrix <double> (7, 2, 0.0);
		this->dspacecraft_state_forward_propdTOF = EMTG::math::Matrix <double>(7, 2, 0.0);
		this->dspacecraft_state_backwarddTOF = EMTG::math::Matrix <double>(7, 2, 0.0);
		this->dspacecraft_state_backward_propdTOF = EMTG::math::Matrix <double>(7, 2, 0.0);
		this->dspacecraft_state_end_coastdTOF = EMTG::math::Matrix <double>(7, 2, 0.0);

		//How does the current FBLT segment's starting epoch change w.r.t. changes in the phase flight times
		this->dcurrent_epochdTOF = std::vector <double>(2, 0.0);


		//containers for the output method
		this->augmented_state_at_initial_coast_midpoint = std::vector <double>(11 + 11 * 11, 0.0);
		this->augmented_state_at_terminal_coast_midpoint = std::vector <double>(11 + 11 * 11, 0.0);
		this->dummy_state_TOF_derivatives = EMTG::math::Matrix <double>(7, 2, 0.0);



		//empty control vector for when we are coasting
	    for (int step = 0; step < options.num_timesteps; ++step) 
        {
            this->spacecraft_state.push_back(vector7_dummy);
            this->control.push_back(vector3_dummy);
	    }


        this->event_epochs.resize(options.num_timesteps);
        this->available_power.resize(options.num_timesteps);
        this->available_thrust.resize(options.num_timesteps);
        this->available_mass_flow_rate.resize(options.num_timesteps);
        this->available_Isp.resize(options.num_timesteps);
        this->active_power.resize(options.num_timesteps);
        this->number_of_active_engines.resize(options.num_timesteps);

        //vector to track the state and derivatives of the central body
        std::vector<double> central_body_state_dummy(options.derivative_type > 2 ? 12 : 6);
        for (size_t step = 0; step < options.num_timesteps; ++step)
            this->central_body_state_mks.push_back(central_body_state_dummy);

	    //set up the integrator
		//for analytical FBLT derivatives, we need to integrate STM
		//entries as states as well, so we must instantiate the integrator 
		//with the normal s/c states but also slots for the STM entries
		this->STMrows = 11;
		this->STMcolumns = 11;
		this->num_states = 11 + 11 * 11;
	    integrator = new EMTG::integration::rk8713M(num_states, options.number_of_phases[j]);

	    //size the time step vector
	    this->time_step_sizes.resize(options.num_timesteps);

	    //set derivatives for spirals
	    this->spiral_escape_dm_after_dm_before = 1.0;

        //set up STMS
        for (size_t step = 0; step < options.num_timesteps / 2; ++step)
        {
            this->STM_archive_forward.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
            this->STM_archive_backward.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
            this->forward_cumulative_STM_archive.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
            this->backward_cumulative_STM_archive.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
        }

        if ((j == 0 && p == 0 && options.forced_post_launch_coast > 1.0e-6) || ((p > 0 || p == 0 && (options.journey_departure_type[j] == 3 || options.journey_departure_type[j] == 4 || options.journey_departure_type[j] == 6)) && options.forced_flyby_coast > 1.0e-6))
        {
            this->detect_initial_coast = true;
            this->initial_coast_STM = EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0);
        }
        else
            this->detect_initial_coast = false;

        if ((p < options.number_of_phases[j] - 1 || (options.journey_arrival_type[j] == 2 || options.journey_arrival_type[j] == 5)) && options.forced_flyby_coast > 1.0e-6)
        {
            this->detect_terminal_coast = true;
            this->terminal_coast_STM = EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0);
        }
        else
            this->detect_terminal_coast = false;
		
        //vector to track the distance from the central body for each step
        std::vector<double> dummy_distance_vector(options.journey_distance_constraint_number_of_bodies[j]);
        std::vector< std::vector <double> > body_position_dummy(options.journey_distance_constraint_number_of_bodies[j], vector3_dummy);
        for (int step = 0; step < options.num_timesteps; ++step)
        {
            this->distance_from_body.push_back(dummy_distance_vector);
            this->distance_constraint_relative_position.push_back(body_position_dummy);
        }



        //size the empty control vector
        this->empty_vector.resize((options.engine_type == 4 || options.engine_type == 12 || options.engine_type == 13 ? 4 : 3), 0.0);
    }

    FBLT_phase::~FBLT_phase() {
	    //delete the integrator object
		    delete integrator;
    }

    //evaluate function
    //return 0 if successful, 1 if failure
    int FBLT_phase::evaluate(const double* X,
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
	    //Step 5: For FBLT, we need to know the state of the spacecraft at the right hand side (end) of the phase in order to propagate backward
	    process_right_boundary_condition(X, Xindex, F, Findex, G, Gindex, needG, current_epoch, current_state, current_deltaV, boundary1_state, boundary2_state, j, p, Universe, options);

	    //******************************************************************
	    //Step 6: thrust and propagate forward and back

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
	    if ( (p < options->number_of_phases[j] - 1 ||  (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5) ) && options->forced_flyby_coast > 1.0e-6)
		    total_available_thrust_time -= options->forced_flyby_coast;

	    for (int step = 0; step < options->num_timesteps; ++step)
	    {
		    time_step_sizes[step] *= total_available_thrust_time / step_size_normalization_coefficient;
	    }

	    //Step 6.2: propagate forward
	    phase_time_elapsed_forward = 0.0;
	    //store the initial prefered integration step size
	    double resumeH = time_step_sizes[0]/2.0/Universe->TU;
	    double resumeError = 1.0e-13;

	    //initialize the forward integration
	    for (int k = 0; k < 7; ++k)
		    spacecraft_state_forward[k] = state_at_beginning_of_phase[k];

        //in the FBLT STM test bed I passed ux,uy,uz and TOF around as part of the state
        //to keep things consistent, let's just set these to zero for all time
        //we will know something is wrong if, in fact, they are ever not zero
        for (int k = 7; k < this->STMrows; ++k)
            spacecraft_state_forward[k] = 0.0;


		//WE NEED TO FILL THIS WITH PHASE LEFT HAND BOUNDARY STATE DERIVATIVES (SPICE FINITE DIFFERENCED)
        if (options->derivative_type > 2 && needG)
        {
            for (size_t state = 0; state < 3; ++state)
            {
                dspacecraft_state_forwarddTOF(state, 0) = 0.0;
				dspacecraft_state_forwarddTOF(state + 3, 0) = 0.0;
				dspacecraft_state_forwarddTOF(state, 1) = this->left_boundary_state_derivative[state] / Universe->LU * Universe->TU;
				dspacecraft_state_forwarddTOF(state + 3, 1) = this->left_boundary_state_derivative[state + 3] / Universe->LU * Universe->TU * Universe->TU;
            }
            dspacecraft_state_forwarddTOF(6, 0) = 0.0;
            dspacecraft_state_forwarddTOF(6, 1) = 0.0;
        }
		else
            dspacecraft_state_forwarddTOF.assign_zeros();

		//How does the temporal length of the current FBLT segment change w.r.t. changes in the phase flight times
		double dsegment_timedTOF;

	    //scale the integration state array to LU and TU
	    for (int k = 0; k < 6; ++k)
	    {
		    spacecraft_state_forward[k] /= Universe->LU;
	    }

	    for (int k = 3; k < 6; ++k)
	    {
		    spacecraft_state_forward[k] *= Universe->TU;
	    }
	    //scale the mass
	    spacecraft_state_forward[6] /= options->maximum_mass;

		//the partial derivative of the current epoch w.r.t. phase flight times is 
		//seeded with either a 0.0 if we are considering a left-hand boundary epoch w.r.t. the current
		//phase flight time or a 1.0 in all other situations
		//since we are about to seed a forward propagation, the seeding is as follows
		dcurrent_epochdTOF[0] = 0.0;
		dcurrent_epochdTOF[1] = 1.0;

	    //Step 6.2.0.1 if there is an initial coast, propagate through it
	    if (detect_initial_coast)
	    {
			//A coast is a user-specified FIXED time and is not affected by the phase TOF decision variables
			dsegment_timedTOF = 0.0;

		    //if this is a launch AND we are doing a forced post-launch initial coast
		    if (j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6)
		    {
			    //initial coast after launch
			    initial_coast_duration = options->forced_post_launch_coast;
		    }
		    else
		    {	
			    //initial coast after flyby
			    initial_coast_duration = options->forced_flyby_coast;
		    }

            double resumeH = initial_coast_duration * 86400 / Universe->TU;

            //The initial coast STM entries of the state vector must be initialized to the identity
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_forward[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_forward[i] = 1.0;

		    integrator->adaptive_step_int(	spacecraft_state_forward,
											dspacecraft_state_forwarddTOF,
                                            spacecraft_state_end_coast,
											dspacecraft_state_end_coastdTOF,
										    empty_vector, 
										    (phase_start_epoch) / Universe->TU,
											dcurrent_epochdTOF,
										    X[0],
                                            initial_coast_duration / Universe->TU,
											dsegment_timedTOF,
										    &resumeH,
										    &resumeError,
											options->integrator_tolerance,
										    EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &temp_available_thrust,
                                            &temp_available_mass_flow_rate,
                                            &temp_available_Isp,
                                            &temp_available_power,
                                            &temp_active_power,
                                            &temp_number_of_active_engines,
                                            this->STMrows,
                                            this->STMcolumns,
										    (void*)options,
										    (void*)Universe,
										    DummyControllerPointer      );	

            //Store the initial coast's first half STM
            int statecount = this->STMrows;
            for (size_t i = 0; i < this->STMrows; ++i)
            {
                for (size_t j = 0; j < this->STMcolumns; ++j)
                {
                    initial_coast_STM(i, j) = spacecraft_state_end_coast[statecount];
                    ++statecount;
                }
            }


            phase_time_elapsed_forward += initial_coast_duration;

			for (size_t i = 0; i < 2; ++i)
				dcurrent_epochdTOF[i] += dsegment_timedTOF;

			//pass on the states and their phase TOF derivatives
            for (int k = 0; k < 7; ++k)
				spacecraft_state_forward = spacecraft_state_end_coast;

			dspacecraft_state_forwarddTOF = dspacecraft_state_end_coastdTOF;
        }


		//////////////////////////////////////////////////////////////////
		//
		// PROPAGATE THROUGH FORWARD FBLT SEGMENTS
		//
		//////////////////////////////////////////////////////////////////

		dsegment_timedTOF = 1.0 / options->num_timesteps;

        for (int step = 0; step < options->num_timesteps / 2; ++step)
        {
            //step 6.2.1 extract the control unit vector from the decision vector
            control[step][0] = X[*Xindex];
            control[step][1] = X[*Xindex + 1];
            control[step][2] = X[*Xindex + 2];


            //extract the specific impulse for this step (VSI only)
            if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
            {
                available_Isp[step] = X[*Xindex + 3];
                ++(*Xindex);
            }

            double throttle = math::norm(control[step].data(), 3) + 1.0e-10;
            if (options->derivative_type > 0 && needG)
            {
                G[control_vector_G_indices[step][0]] = 2.0 * control[step][2] / throttle;
                G[control_vector_G_indices[step][1]] = 2.0 * control[step][1] / throttle;
                G[control_vector_G_indices[step][2]] = 2.0 * control[step][0] / throttle;
                (*Gindex) += 3;
            }
            (*Xindex) += 3;

            //step 6.2.2 apply the control unit vector magnitude constraint
            F[*Findex] = throttle;
            ++(*Findex);

            //step 6.2.3 propagate the spacecraft to the end of the FBLT step using the control unit vector
 
            //The STM entries of the state vector must be initialized to the identity before every step
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_forward[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_forward[i] = 1.0;


#ifdef ENABLE_FBLT_INSTRUMENTATION

            //before we propagate with derivatives, let's do a "mini finite differencing" check
			//NOTE: NOT set up to work with coasts presently
            double control_perturbation = 1.0e-6;
            double dstate_du[7][3]; //indexed as [state][control]
            vector<double> reference_control = this->control[step];
            std::vector <double> spacecraft_state_forward_prop_plus(11 + 11 * 11, 0.0);
            std::vector <double> spacecraft_state_forward_prop_minus(11 + 11 * 11, 0.0);

			EMTG::math::Matrix <double> dspacecraft_state_forward_prop_plusdTOF(7, 2, 0.0);
			EMTG::math::Matrix <double> dspacecraft_state_forward_prop_minusdTOF(7, 2, 0.0);

            for (size_t cindex = 0; cindex < 3; ++cindex)
            {
                vector<double> control_minus = reference_control;
                vector<double> control_plus = reference_control;
                control_minus[cindex] -= control_perturbation;
                control_plus[cindex] += control_perturbation;

                //get the forward step evaluation
                for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_forward[i] = 0.0;
                for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_forward[i] = 1.0;

                integrator->adaptive_step_int(spacecraft_state_forward,
											dspacecraft_state_forwarddTOF,
                                            spacecraft_state_forward_prop_plus,
											dspacecraft_state_forward_prop_plusdTOF,
                                            control_plus,
                                            (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
											dcurrent_epochdTOF,
                                            X[0],
                                            time_step_sizes[step] / Universe->TU,
											dsegment_timedTOF,
                                            &resumeH,
                                            &resumeError,
											options->integrator_tolerance,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &available_thrust[step],
                                            &available_mass_flow_rate[step],
                                            &available_Isp[step],
                                            &available_power[step],
                                            &active_power[step],
                                            &number_of_active_engines[step],
                                            this->STMrows,
                                            this->STMcolumns,
                                            (void*)options,
                                            (void*)Universe,
                                            DummyControllerPointer);

                //get the backward step evaluation
                for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_forward[i] = 0.0;
                for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_forward[i] = 1.0;

                integrator->adaptive_step_int(spacecraft_state_forward,
											dspacecraft_state_forwarddTOF,
                                            spacecraft_state_forward_prop_minus,
											dspacecraft_state_forward_prop_minusdTOF,
                                            control_minus,
                                            (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
											dcurrent_epochdTOF,
                                            X[0],
                                            time_step_sizes[step] / Universe->TU,
											dsegment_timedTOF,
                                            &resumeH,
                                            &resumeError,
											options->integrator_tolerance,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &available_thrust[step],
                                            &available_mass_flow_rate[step],
                                            &available_Isp[step],
                                            &available_power[step],
                                            &active_power[step],
                                            &number_of_active_engines[step],
                                            this->STMrows,
                                            this->STMcolumns,
                                            (void*)options,
                                            (void*)Universe,
                                            DummyControllerPointer);

                for (size_t stateindex = 0; stateindex < 7; ++stateindex)
					dstate_du[stateindex][cindex] = (spacecraft_state_forward_prop_plus[stateindex] - spacecraft_state_forward_prop_minus[stateindex]) / (2.0 * control_perturbation);
            }

            //now generate the actual STM
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_forward[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_forward[i] = 1.0;
#endif
			 
			

            integrator->adaptive_step_int( spacecraft_state_forward,
				                           dspacecraft_state_forwarddTOF,
                                           spacecraft_state_forward_prop,
										   dspacecraft_state_forward_propdTOF,
                                           control[step],
                                           (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
										   dcurrent_epochdTOF,
                                           X[0],
                                           time_step_sizes[step] / Universe->TU,
										   dsegment_timedTOF,
                                           &resumeH,
                                           &resumeError,
										   options->integrator_tolerance,
                                           EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                           &available_thrust[step],
                                           &available_mass_flow_rate[step],
                                           &available_Isp[step],
                                           &available_power[step],
                                           &active_power[step],
                                           &number_of_active_engines[step],
                                           this->STMrows,
                                           this->STMcolumns,
                                           (void*)options,
                                           (void*)Universe,
                                           DummyControllerPointer);


            //Store this step's STM in the forward archive
            int statecount = this->STMrows;
            for (size_t i = 0; i < this->STMrows; ++i)
            {
                for (size_t j = 0; j < this->STMcolumns; ++j)
                {
                    STM_archive_forward[step](i, j) = spacecraft_state_forward_prop[statecount];
                    ++statecount;
                }
            }
            
#ifdef ENABLE_FBLT_INSTRUMENTATION
            //finally, compare the propagated STM to the central differenced STM
            std::cout << "Difference between propagated STM and central differenced STM, forward step " << step << std::endl;
            std::cout << "Written as: index, propagated, differenced, absolute error, relative error" << std::endl;
            vector<string> statenames;
            vector<string> cnames;
            statenames.push_back("x");
            statenames.push_back("y");
            statenames.push_back("z");
            statenames.push_back("vx");
            statenames.push_back("vy");
            statenames.push_back("vz");
            statenames.push_back("m");
            cnames.push_back("ux");
            cnames.push_back("uy");
            cnames.push_back("uz");
            for (size_t state = 0; state < 7; ++state)
            {
                for (size_t cindex = 0; cindex < 3; ++cindex)
                {
                    std::cout << "[" << statenames[state] << "," << cnames[cindex] << "] ";
                    std::cout << string_utilities::convert_number_to_formatted_string(STM_archive_forward[step](state, 7 + cindex), 2) << "   "
                    << string_utilities::convert_number_to_formatted_string(dstate_du[state][cindex], 2) << "    "
                    << string_utilities::convert_number_to_formatted_string(STM_archive_forward[step](state, 7 + cindex) - dstate_du[state][cindex], 2) << "    "
                    << string_utilities::convert_number_to_formatted_string((STM_archive_forward[step](state, 7 + cindex) - dstate_du[state][cindex]) / dstate_du[state][cindex], 2)
                    << std::endl;
                }
            }
            getchar();
            //end local STM check
#endif

            //copy the state over
            for (size_t state = 0; state < 11; ++state)
                spacecraft_state_forward[state] = spacecraft_state_forward_prop[state];

			//copy the state phase TOF derivatives over
			dspacecraft_state_forwarddTOF = dspacecraft_state_forward_propdTOF;

            //step 6.2.4 encode the epoch of the step midpoint
            event_epochs[step] = phase_start_epoch + phase_time_elapsed_forward + 0.5 * time_step_sizes[step];

			//advance the clock for the next FBLT segment
            phase_time_elapsed_forward += time_step_sizes[step];
			
			for (size_t i = 0; i < 2; ++i)
				dcurrent_epochdTOF[i] += dsegment_timedTOF;

            //step 6.2.5 apply the forward body distance constraints
            for (int body = 0; body < options->journey_distance_constraint_number_of_bodies[j]; ++body)
            {
                //step 6.2.5.1 compute the distance from the spacecraft to the body
                //(this is a different procedure if the body is the central body)
                if (options->journey_distance_constraint_bodies[j][body] == -2)
                {
                    this->distance_constraint_relative_position[step][body][0] = this->spacecraft_state[step][0];
                    this->distance_constraint_relative_position[step][body][1] = this->spacecraft_state[step][1];
                    this->distance_constraint_relative_position[step][body][2] = this->spacecraft_state[step][2];
                }
                else
                {
                    double body_state[6]; //TODO THIS MUST CHANGE IF WE WANT TIME DERIVATIVES FOR DISTANCE CONSTRAINT
                    Universe->bodies[options->journey_distance_constraint_bodies[j][body]].locate_body(this->event_epochs[step], body_state, false, options);
                    this->distance_constraint_relative_position[step][body][0] = this->spacecraft_state[step][0] - body_state[0];
                    this->distance_constraint_relative_position[step][body][1] = this->spacecraft_state[step][1] - body_state[1];
                    this->distance_constraint_relative_position[step][body][2] = this->spacecraft_state[step][2] - body_state[2];
                }
                this->distance_from_body[step][body] = math::norm(this->distance_constraint_relative_position[step][body].data(), 3);

                //step 6.2.5.2 apply the constraint
                F[*Findex] = this->distance_from_body[step][body] / Universe->LU;
                ++(*Findex);
            }
        }
		
#ifdef ENABLE_FBLT_INSTRUMENTATION
        //now let's repeat the integration inside a finite difference loop to check the multiplication of the STMs
        double control_perturbation = 1.0e-6;
        std::vector<double> temp_control;
        double dstate_du[7][3]; //indexed as [state][control]
		std::vector<double> spacecraft_state_forward_prop_plus(11 + 11 * 11, 0.0);
		std::vector<double> spacecraft_state_forward_prop_plus_2(11 + 11 * 11, 0.0);
		std::vector<double> spacecraft_state_forward_prop_minus(11 + 11 * 11, 0.0);
		std::vector<double> spacecraft_state_forward_prop_minus_2(11 + 11 * 11, 0.0);

		EMTG::math::Matrix <double> dspacecraft_state_forward_prop_plusdTOF(7, 2, 0.0);
		EMTG::math::Matrix <double> dspacecraft_state_forward_prop_plus_2dTOF(7, 2, 0.0);
		EMTG::math::Matrix <double> dspacecraft_state_forward_prop_minusdTOF(7, 2, 0.0);
		EMTG::math::Matrix <double> dspacecraft_state_forward_prop_minus_2dTOF(7, 2, 0.0);

        for (size_t cindex = 0; cindex < 3; ++cindex)
        {
            //the following only works when there is NO initial coast
            for (size_t state = 0; state < 3; ++state)
            {
                spacecraft_state_forward_prop_plus[state] = this->state_at_beginning_of_phase[state] / Universe->LU;
                spacecraft_state_forward_prop_minus[state] = this->state_at_beginning_of_phase[state] / Universe->LU;
                spacecraft_state_forward_prop_plus[state + 3] = this->state_at_beginning_of_phase[state + 3] / Universe->LU * Universe->TU;
                spacecraft_state_forward_prop_minus[state + 3] = this->state_at_beginning_of_phase[state + 3] / Universe->LU * Universe->TU;
            }
            spacecraft_state_forward_prop_plus[6] = this->state_at_beginning_of_phase[6] / options->maximum_mass;
            spacecraft_state_forward_prop_minus[6] = this->state_at_beginning_of_phase[6] / options->maximum_mass;
            for (size_t placeholder = 7; placeholder < 11; ++placeholder)
            {
                spacecraft_state_forward_prop_plus[placeholder] = 0.0;
                spacecraft_state_forward_prop_minus[placeholder] = 0.0;
            }

            //first the forward step of the central difference
            for (int step = 0; step < options->num_timesteps / 2; ++step)
            {
                //clear the local STM

                for (size_t i = this->STMrows; i < this->num_states; ++i)
                    spacecraft_state_forward_prop_plus[i] = 0.0;
                for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                    spacecraft_state_forward_prop_plus[i] = 1.0;

                //set the perturbed control
                temp_control = control[step];
                if (step == 0)
                    temp_control[cindex] += control_perturbation;

                for (size_t uindex = 0; uindex < 3; ++uindex)
                    spacecraft_state_forward_prop_plus[7 + uindex] = temp_control[uindex];

                //propagate the spacecraft to the end of the FBLT step using the control unit vector

                integrator->adaptive_step_int(spacecraft_state_forward_prop_plus,
												dspacecraft_state_forward_prop_plusdTOF,
                                                spacecraft_state_forward_prop_plus_2,
												dspacecraft_state_forward_prop_plus_2dTOF,
                                                temp_control,
                                                (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
												dcurrent_epochdTOF,
                                                X[0],
                                                time_step_sizes[step] / Universe->TU,
												dsegment_timedTOF,
                                                &resumeH,
                                                &resumeError,
												options->integrator_tolerance,
                                                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                                &available_thrust[step],
                                                &available_mass_flow_rate[step],
                                                &available_Isp[step],
                                                &available_power[step],
                                                &active_power[step],
                                                &number_of_active_engines[step],
                                                this->STMrows,
                                                this->STMcolumns,
                                                (void*)options,
                                                (void*)Universe,
                                                DummyControllerPointer);

                //copy the state over
                for (size_t state = 0; state < 11; ++state)
                    spacecraft_state_forward_prop_plus[state] = spacecraft_state_forward_prop_plus_2[state];
            }

            //then the backward step of the central difference
            for (int step = 0; step < options->num_timesteps / 2; ++step)
            {
                //clear the local STM

                for (size_t i = this->STMrows; i < this->num_states; ++i)
                    spacecraft_state_forward_prop_minus[i] = 0.0;
                for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                    spacecraft_state_forward_prop_minus[i] = 1.0;

                //set the perturbed control
                temp_control = control[step];
                if (step == 0)
                    temp_control[cindex] -= control_perturbation;

                for (size_t uindex = 0; uindex < 3; ++uindex)
                    spacecraft_state_forward_prop_minus[7 + uindex] = temp_control[uindex];

                //propagate the spacecraft to the end of the FBLT step using the control unit vector

                integrator->adaptive_step_int(spacecraft_state_forward_prop_minus,
												dspacecraft_state_forward_prop_minusdTOF,
                                                spacecraft_state_forward_prop_minus_2,
												dspacecraft_state_forward_prop_minus_2dTOF,
                                                temp_control,
                                                (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
												dcurrent_epochdTOF,
                                                X[0],
                                                time_step_sizes[step] / Universe->TU,
												dsegment_timedTOF,
                                                &resumeH,
                                                &resumeError,
                                                options->integrator_tolerance,
                                                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                                &available_thrust[step],
                                                &available_mass_flow_rate[step],
                                                &available_Isp[step],
                                                &available_power[step],
                                                &active_power[step],
                                                &number_of_active_engines[step],
                                                this->STMrows,
                                                this->STMcolumns,
                                                (void*)options,
                                                (void*)Universe,
                                                DummyControllerPointer);

                //copy the state over
                for (size_t state = 0; state < 11; ++state)
                    spacecraft_state_forward_prop_minus[state] = spacecraft_state_forward_prop_minus_2[state];
            }

            for (size_t state = 0; state < 7; ++state)
                dstate_du[state][cindex] = (spacecraft_state_forward_prop_plus[state] - spacecraft_state_forward_prop_minus[state]) / (2.0 * control_perturbation);
        }

        //finally, compare the propagated STM to the central differenced STM
        std::cout << "Difference between propagated STM and central differenced STM after " << options->num_timesteps / 2 << " steps" << std::endl;
        std::cout << "Written as: index, propagated, differenced, absolute error, relative error" << std::endl;
        std::vector<string> statenames;
        std::vector<string> cnames;
        statenames.push_back("x");
        statenames.push_back("y");
        statenames.push_back("z");
        statenames.push_back("vx");
        statenames.push_back("vy");
        statenames.push_back("vz");
        statenames.push_back("m");
        cnames.push_back("ux");
        cnames.push_back("uy");
        cnames.push_back("uz");
        EMTG::math::Matrix <double> forward_cumulative_STM(this->STMrows, this->STMcolumns, 0.0);
        for (size_t i = 0; i < this->STMrows; ++i)
            forward_cumulative_STM(i, i) = 1.0;

        //build the derivatives matrix through successive STM multiplication
        //every time step we add one more STM to the chain
        for (int step = options->num_timesteps / 2 - 1; step >= 0; --step)
        {
            EMTG::math::Matrix <double> stripped_step_STM = STM_archive_forward[step];
            if (step > 0)
            {
                for (size_t row = 0; row < 7; ++row)
                {
                    for (size_t column = 7; column < 10; ++column)
                        stripped_step_STM(row, column) = 0.0;
                }
            }
            forward_cumulative_STM *= stripped_step_STM;
        }

        for (size_t state = 0; state < 7; ++state)
        {
            for (size_t cindex = 0; cindex < 3; ++cindex)
            {
                std::cout << "[" << statenames[state] << "," << cnames[cindex] << "] ";
                std::cout << string_utilities::convert_number_to_formatted_string(forward_cumulative_STM(state, 7 + cindex), 2) << "   "
                    << string_utilities::convert_number_to_formatted_string(dstate_du[state][cindex], 2) << "    "
                    << string_utilities::convert_number_to_formatted_string(forward_cumulative_STM(state, 7 + cindex) - dstate_du[state][cindex], 2) << "    "
                    << string_utilities::convert_number_to_formatted_string((forward_cumulative_STM(state, 7 + cindex) - dstate_du[state][cindex]) / dstate_du[state][cindex], 2)
                    << std::endl;
            }
        }

        cout << endl;
        cout << "Check of propagated state vs propagated state+ and propagated state-" << endl;
        for (size_t state = 0; state < 7; ++state)
        {
            std::cout << "[" << statenames[state] << "] ";
            cout << string_utilities::convert_number_to_formatted_string(spacecraft_state_forward[state], 2) << "   "
                << string_utilities::convert_number_to_formatted_string(spacecraft_state_forward_prop_plus[state], 2) << "   "
                << string_utilities::convert_number_to_formatted_string(spacecraft_state_forward_prop_minus[state], 2) << endl;
        }

        getchar();
        
        //end STM checker

#endif
        
	    //Step 6.3: propagate backward
	    phase_time_elapsed_backward = 0.0;
	    //store the initial prefered integration step size
	    resumeH = -time_step_sizes[options->num_timesteps - 1] / Universe->TU;
	    resumeError = 1.0e-13;
	
	   
		//WE NEED TO FILL THIS WITH PHASE RIGHT HAND BOUNDARY STATE DERIVATIVES (SPICE FINITE DIFFERENCED)
        if (options->derivative_type > 2 && needG)
        {
            for (size_t state = 0; state < 3; ++state)
            {
                dspacecraft_state_backwarddTOF(state, 0) = this->right_boundary_state_derivative[state] / Universe->LU * Universe->TU;
				dspacecraft_state_backwarddTOF(state + 3, 0) = this->right_boundary_state_derivative[state + 3] / Universe->LU * Universe->TU * Universe->TU;
				dspacecraft_state_backwarddTOF(state, 1) = this->right_boundary_state_derivative[state] / Universe->LU * Universe->TU;
				dspacecraft_state_backwarddTOF(state + 3, 1) = this->right_boundary_state_derivative[state + 3] / Universe->LU * Universe->TU * Universe->TU;
            }
            dspacecraft_state_backwarddTOF(6, 0) = 0.0;
            dspacecraft_state_backwarddTOF(6, 1) = 0.0;
        }
        else
            dspacecraft_state_backwarddTOF.assign_zeros();

	    //first initialize the backward integration
	    for (int k = 0; k < 7; ++k)
		    spacecraft_state_backward[k] = state_at_end_of_phase[k];

		for (int k = 7; k < this->STMrows; ++k)
			spacecraft_state_backward[k] = 0.0;

	    //scale the integration state array to LU and TU
	    for (int k = 0; k < 6; ++k)
	    {
		    spacecraft_state_backward[k] /= Universe->LU;
	    }

	    for (int k = 3; k < 6; ++k)
	    {
		    spacecraft_state_backward[k] *= Universe->TU;
	    }

	    //scale the mass
	    spacecraft_state_backward[6] /= options->maximum_mass;

		//the partial derivative of the current epoch w.r.t. phase flight times is as follows...
		//same as the forward, except now the current phase flight time DOES impact the current epoch
		dcurrent_epochdTOF[0] = 1.0;
		dcurrent_epochdTOF[1] = 1.0;
		



	    //Step 6.3.0.1 if there is an terminal coast, propagate through it
	    if (this->detect_terminal_coast)
	    {
			//A coast is a user-specified FIXED time and is not affected by the decision variables
			dsegment_timedTOF = 0.0;

			//zero out the coast RHS state container
		    std::fill(spacecraft_state_end_coast.begin(), spacecraft_state_end_coast.end(), 0.0);

		    //initial coast after flyby
		    terminal_coast_duration = options->forced_flyby_coast;

		    double resumeH = -terminal_coast_duration / 2.0 / Universe->TU;

			//The terminal coast STM entries of the state vector must be initialized to the identity
			for (size_t i = this->STMrows; i < this->num_states; ++i)
				spacecraft_state_backward[i] = 0.0;
			for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
				spacecraft_state_backward[i] = 1.0;


		    integrator->adaptive_step_int(	spacecraft_state_backward,
											dspacecraft_state_backwarddTOF,
                                            spacecraft_state_end_coast,
											dspacecraft_state_end_coastdTOF,
										    empty_vector, 
										    (phase_end_epoch) / Universe->TU,
											dcurrent_epochdTOF,
										    X[0],
										    -terminal_coast_duration / Universe->TU, 
											dsegment_timedTOF,
										    &resumeH,
										    &resumeError,
                                            options->integrator_tolerance,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &temp_available_thrust,
                                            &temp_available_mass_flow_rate,
                                            &temp_available_Isp,
                                            &temp_available_power,
                                            &temp_active_power,
                                            &temp_number_of_active_engines,
											this->STMrows,
											this->STMcolumns,
										    (void*)options,
										    (void*)Universe,
										    DummyControllerPointer      );	

			//Store the terminal coast's first-half STM 
			int statecount = this->STMrows;
			for (size_t i = 0; i < this->STMrows; ++i)
			{
				for (size_t j = 0; j < this->STMcolumns; ++j)
				{
                    terminal_coast_STM(i, j) = spacecraft_state_end_coast[statecount];
					++statecount;
				}
			}

		    phase_time_elapsed_backward += terminal_coast_duration;

			for (size_t i = 0; i < 2; ++i)
				dcurrent_epochdTOF[i] += dsegment_timedTOF;

			//pass on the states and their phase TOF derivatives
		    for (int k = 0; k < 7; ++k)
			    spacecraft_state_backward[k] = spacecraft_state_end_coast[k];

			dspacecraft_state_backwarddTOF = dspacecraft_state_end_coastdTOF;
	    }


		//////////////////////////////////////////////////////////////////
		//
		// PROPAGATE THROUGH BACKWARD FBLT SEGMENTS
		//
		//////////////////////////////////////////////////////////////////
		dsegment_timedTOF = 1.0 / options->num_timesteps;

	    for (int step = 0; step < options->num_timesteps/2; ++step)
	    {
		    //translate into backward steps
		    size_t backstep = options->num_timesteps - 1 - step;

		    //step 6.3.1 extract the burn parameters from the decision vector
		    //extract the specific impulse for this step (VSI only)
		    if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
		    {
			    control[backstep][0] = X[*Xindex+4*(backstep - options->num_timesteps/2)];
			    control[backstep][1] = X[*Xindex+1+4*(backstep - options->num_timesteps/2)];
			    control[backstep][2] = X[*Xindex+2+4*(backstep - options->num_timesteps/2)];
			    available_Isp[backstep] = X[*Xindex+3+4*(backstep - options->num_timesteps/2)];
		    }
		    else
		    {
			    control[backstep][0] = X[*Xindex+3*(backstep - options->num_timesteps/2)];
			    control[backstep][1] = X[*Xindex+1+3*(backstep - options->num_timesteps/2)];
			    control[backstep][2] = X[*Xindex+2+3*(backstep - options->num_timesteps/2)];
		    }

			double throttle = math::norm(control[backstep].data(), 3) + 1.0e-10;
		    if (options->derivative_type > 0 && needG)
		    {
			    G[control_vector_G_indices[backstep][0]] = 2.0 * control[backstep][2] / throttle;
			    G[control_vector_G_indices[backstep][1]] = 2.0 * control[backstep][1] / throttle;
			    G[control_vector_G_indices[backstep][2]] = 2.0 * control[backstep][0] / throttle;
			    (*Gindex) += 3;
		    }

		    //step 6.3.2 apply the control unit vector magnitude constraint
            F[*Findex + (backstep - options->num_timesteps / 2) * (1 + options->journey_distance_constraint_number_of_bodies[j])] = throttle;
		

			//The STM entries of the state vector must be initialized to the identity before every step
			for (size_t i = this->STMrows; i < this->num_states; ++i)
				spacecraft_state_backward[i] = 0.0;
			for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
				spacecraft_state_backward[i] = 1.0;

            
#ifdef ENABLE_FBLT_INSTRUMENTATION
            //before we propagate with derivatives, let's do a "mini finite differencing" check
            double control_perturbation = 1.0e-6;
            double dstate_du[7][3]; //indexed as [state][control]
            std::vector<double> reference_control = this->control[backstep];
            std::vector <double> spacecraft_state_backward_prop_plus(11 + 11 * 11, 0.0);
            std::vector <double> spacecraft_state_backward_prop_minus(11 + 11 * 11, 0.0);

			EMTG::math::Matrix <double> dspacecraft_state_backward_prop_plusdTOF(7, 2, 0.0);
			EMTG::math::Matrix <double> dspacecraft_state_backward_prop_minusdTOF(7, 2, 0.0);

            for (size_t cindex = 0; cindex < 3; ++cindex)
            {
                vector<double> control_minus = reference_control;
                vector<double> control_plus = reference_control;
                control_minus[cindex] -= control_perturbation;
                control_plus[cindex] += control_perturbation;

                //get the forward step evaluation
                for (size_t i = this->STMrows; i < this->num_states; ++i)
                    spacecraft_state_backward[i] = 0.0;
                for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                    spacecraft_state_backward[i] = 1.0;

                integrator->adaptive_step_int(spacecraft_state_backward,
											dspacecraft_state_backwarddTOF,
                                            spacecraft_state_backward_prop_plus,
											dspacecraft_state_backward_prop_plusdTOF,
                                            control_plus,
                                            (phase_end_epoch - phase_time_elapsed_backward) / Universe->TU,
											dcurrent_epochdTOF,
                                            X[0],
                                            -time_step_sizes[backstep] / Universe->TU,
											dsegment_timedTOF,
                                            &resumeH,
                                            &resumeError,
                                            options->integrator_tolerance,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &available_thrust[backstep],
                                            &available_mass_flow_rate[backstep],
                                            &available_Isp[backstep],
                                            &available_power[backstep],
                                            &active_power[backstep],
                                            &number_of_active_engines[backstep],
                                            this->STMrows,
                                            this->STMcolumns,
                                            (void*)options,
                                            (void*)Universe,
                                            DummyControllerPointer);

                //get the backward step evaluation
                for (size_t i = this->STMrows; i < this->num_states; ++i)
                    spacecraft_state_backward[i] = 0.0;
                for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                    spacecraft_state_backward[i] = 1.0;

                integrator->adaptive_step_int(spacecraft_state_backward,
											dspacecraft_state_backwarddTOF,
                                            spacecraft_state_backward_prop_minus,
											dspacecraft_state_backward_prop_minusdTOF,
                                            control_minus,
                                            (phase_end_epoch - phase_time_elapsed_backward) / Universe->TU,
											dcurrent_epochdTOF,
                                            X[0],
                                            -time_step_sizes[backstep] / Universe->TU,
											dsegment_timedTOF,
                                            &resumeH,
                                            &resumeError,
                                            options->integrator_tolerance,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &available_thrust[backstep],
                                            &available_mass_flow_rate[backstep],
                                            &available_Isp[backstep],
                                            &available_power[backstep],
                                            &active_power[backstep],
                                            &number_of_active_engines[backstep],
                                            this->STMrows,
                                            this->STMcolumns,
                                            (void*)options,
                                            (void*)Universe,
                                            DummyControllerPointer);

                for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    dstate_du[stateindex][cindex] = (spacecraft_state_backward_prop_plus[stateindex] - spacecraft_state_backward_prop_minus[stateindex]) / (2.0 * control_perturbation);
            }
#endif
		
		    //step 6.3.3 propagate the spacecraft to the midpoint of the step using the control unit vector
		    integrator->adaptive_step_int(	spacecraft_state_backward,
											dspacecraft_state_backwarddTOF,
										    spacecraft_state_backward_prop,
											dspacecraft_state_backward_propdTOF,
											control[backstep],  
										    (phase_end_epoch - phase_time_elapsed_backward) / Universe->TU,
											dcurrent_epochdTOF,
										    X[0],
										    -time_step_sizes[backstep] / Universe->TU, 
											dsegment_timedTOF,
										    &resumeH,
										    &resumeError,
											options->integrator_tolerance,
										    EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
										    &available_thrust[backstep],
										    &available_mass_flow_rate[backstep],
										    &available_Isp[backstep],
										    &available_power[backstep],
										    &active_power[backstep],
										    &number_of_active_engines[backstep],
											this->STMrows,
											this->STMcolumns,
										    (void*)options,
										    (void*)Universe,
										    DummyControllerPointer                                                  );
		

			//Store this half-step's STM in the backwards archive
			int statecount = this->STMrows;
			for (size_t i = 0; i < this->STMrows; ++i)
			{
				for (size_t j = 0; j < this->STMcolumns; ++j)
				{
                    STM_archive_backward[step](i, j) = spacecraft_state_backward_prop[statecount];
                    spacecraft_state_backward[statecount] = spacecraft_state_backward_prop[statecount];
					++statecount;
				}
			}
            
#ifdef ENABLE_FBLT_INSTRUMENTATION
            //finally, compare the propagated STM to the central differenced STM
            std::cout << "Difference between propagated STM and central differenced STM, backward step " << step << std::endl;
            std::cout << "Written as: index, propagated, differenced, absolute error, relative error" << std::endl;
            vector<string> statenames;
            vector<string> cnames;
            statenames.push_back("x");
            statenames.push_back("y");
            statenames.push_back("z");
            statenames.push_back("vx");
            statenames.push_back("vy");
            statenames.push_back("vz");
            statenames.push_back("m");
            cnames.push_back("ux");
            cnames.push_back("uy");
            cnames.push_back("uz");
            for (size_t state = 0; state < 7; ++state)
            {

                for (size_t cindex = 0; cindex < 3; ++cindex)
                {
                    std::cout << "[" << statenames[state] << "," << cnames[cindex] << "] ";
                    std::cout << string_utilities::convert_number_to_formatted_string(STM_archive_backward[step](state, 7 + cindex), 2) << "   "
                        << string_utilities::convert_number_to_formatted_string(dstate_du[state][cindex], 2) << "    "
                        << string_utilities::convert_number_to_formatted_string(STM_archive_backward[step](state, 7 + cindex) - dstate_du[state][cindex], 2) << "    "
                        << string_utilities::convert_number_to_formatted_string((STM_archive_backward[step](state, 7 + cindex) - dstate_du[state][cindex]) / dstate_du[state][cindex], 2)
                        << std::endl;
                }
            }
            getchar();
            //end local STM check
#endif

            //copy the state over
            for (size_t state = 0; state < 11; ++state)
                spacecraft_state_backward[state] = spacecraft_state_backward_prop[state];

			dspacecraft_state_backwarddTOF = dspacecraft_state_backward_propdTOF;

		    //step 6.3.4 encode the epoch of the step midpoint
		    event_epochs[backstep] = phase_end_epoch - phase_time_elapsed_backward - 0.5 * time_step_sizes[backstep];


			//advance the clock and the current epoch derivatives
		    phase_time_elapsed_backward += time_step_sizes[backstep];

			for (size_t i = 0; i < 2; ++i)
				dcurrent_epochdTOF[i] += dsegment_timedTOF;

            //step 6.3.5 apply the backward body distance constraints
            for (int body = 0; body < options->journey_distance_constraint_number_of_bodies[j]; ++body)
            {
                //step 6.3.5.1 compute the distance from the spacecraft to the body
                //(this is a different procedure if the body is the central body)
                if (options->journey_distance_constraint_bodies[j][body] == -2)
                {
                    this->distance_constraint_relative_position[backstep][body][0] = this->spacecraft_state[backstep][0];
                    this->distance_constraint_relative_position[backstep][body][1] = this->spacecraft_state[backstep][1];
                    this->distance_constraint_relative_position[backstep][body][2] = this->spacecraft_state[backstep][2];
                }
                else
                {
                    double body_state[6]; //TODO THIS MUST CHANGE IF WE WANT TIME DERIVATIVES FOR DISTANCE CONSTRAINT
                    Universe->bodies[options->journey_distance_constraint_bodies[j][body]].locate_body(this->event_epochs[backstep], body_state, false, options);
                    this->distance_constraint_relative_position[backstep][body][0] = this->spacecraft_state[backstep][0] - body_state[0];
                    this->distance_constraint_relative_position[backstep][body][1] = this->spacecraft_state[backstep][1] - body_state[1];
                    this->distance_constraint_relative_position[backstep][body][2] = this->spacecraft_state[backstep][2] - body_state[2];
                }
                this->distance_from_body[backstep][body] = math::norm(this->distance_constraint_relative_position[backstep][body].data(), 3);

                //step 6.3.5.2 apply the constraint
                F[*Findex + (backstep - options->num_timesteps / 2) * (1 + options->journey_distance_constraint_number_of_bodies[j]) + body + 1] = this->distance_from_body[backstep][body] / Universe->LU;
            }
	    }
		
	    //step Xindex back to the end of the arc
	    if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
		    (*Xindex) += 4 * options->num_timesteps/2;
	    else
		    (*Xindex) += 3 * options->num_timesteps/2;

	    //step Findex back to the end of the arc
        (*Findex) += options->num_timesteps / 2 * (1 + options->journey_distance_constraint_number_of_bodies[j]);

	    //Step 6.4: enforce match point constraint
	    for (size_t k = 0; k < 3; ++k)
	    {
		    //position
            F[*Findex + k] = (spacecraft_state_backward_prop[k] - spacecraft_state_forward_prop[k]);
		
		    //velocity
            F[*Findex + k + 3] = (spacecraft_state_backward_prop[k + 3] - spacecraft_state_forward_prop[k + 3]);

		    //unscale the match point state (normalized to mks)
            this->match_point_state[k] = spacecraft_state_forward_prop[k] * Universe->LU;
            this->match_point_state[k + 3] = spacecraft_state_forward_prop[k + 3] * Universe->LU / Universe->TU;
	    }
	    //mass
        F[*Findex + 6] = (spacecraft_state_backward_prop[6] - spacecraft_state_forward_prop[6]);
	    (*Findex) += 7;

		//unscale the match point mass
        this->match_point_state[6] = spacecraft_state_forward_prop[6] * options->maximum_mass;


		//CALCULATE MATCH POINT DERIVATIVES HERE

		if (options->derivative_type > 1 && needG)
			this->calculate_match_point_derivatives(G, Gindex, j, p, options, Universe);

	    //******************************************************************
	    //Step 7: process the arrival, if applicable
        if (p == options->number_of_phases[j] - 1)
            this->process_arrival(  current_state,
                                    current_deltaV,
                                    boundary2_state,
                                    current_epoch,
                                    X,
                                    Xindex,
                                    F,
                                    Findex,
                                    G,
                                    j,
                                    p,
                                    needG,
                                    options,
                                    Universe);
	
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

 	    return 0;
    }


    //bounds calculation function
    //return 0 if successful, 1 if failure
    void FBLT_phase::calcbounds(vector<double>* Xupperbounds,
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
        missionoptions* options)
    {
	    //this function calculates the upper and lower bounds for the decision and constraint vectors for FBLT
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
        this->calcbounds_step_distribution_scale_factor(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	    //**************************************************************************
	    //next, we need to include the decision variables and constraints for each burn
        this->calcbounds_LT_controls(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);
	    
	    //**************************************************************************
	    //finally, we encode the match point continuity constraints and their Jacobian entries,
	    //noting that every patch point constraint in the phase has a derivative with respect to every variable in the phase
	    //in addition, the patch point constraints have a derivative with respect to the previous phase's arrival mass
	    //and the patch point constraints have a derivative with respect to all previous time variables, including the launch date
        this->calcbounds_LT_match_points(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	    //***************************************************************************
	    //if this is the last phase, encode any constraints for the arrival processing
        this->calcbounds_arrival_constraints(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);
    }

    //output function
    //return 0 if successful, 1 if failure
    void FBLT_phase::output(missionoptions* options,
        const double& launchdate,
        const int& j,
        const int& p,
        EMTG::Astrodynamics::universe* Universe,
        int* eventcount)
    {
        //Step 1: store data that will be used for the printing
        double phase_time_elapsed = 0.0;
        string event_type;
        double temp_power, temp_thrust, temp_mdot, temp_Isp,
            temp_active_power, temp_dTdP, temp_dmdotdP,
            temp_dTdIsp, temp_dmdotdIsp, temp_dPdr, temp_dPdt;
        int temp_active_thrusters;
        std::vector <double> spacecraft_state_propagate (11 + 11 * 11, 0.0), spacecraft_state_propagate_next (11 + 11 * 11, 0.0);
        double resumeH;
        double resumeError = 1.0e-13;
        double dummy_parameter = 0.0;
		double dsegment_timedTOF = 0.0; //dummy variable for output

		std::fill(spacecraft_state_propagate.begin(), spacecraft_state_propagate.end(), 0.0);
		std::fill(spacecraft_state_propagate_next.begin(), spacecraft_state_propagate_next.end(), 0.0);

        //scale the integration state array to LU and TU
        for (int k = 0; k < 6; ++k)
        {
            spacecraft_state_propagate[k] = this->state_at_beginning_of_phase[k] / Universe->LU;
        }

        for (int k = 3; k < 6; ++k)
        {
            spacecraft_state_propagate[k] *= Universe->TU;
        }
        //scale the mass
        spacecraft_state_propagate[6] = state_at_beginning_of_phase[6] / options->maximum_mass;

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
        else if (j == 0 && boundary1_location_code > 0 && (options->LV_type >= 0 || options->LV_type == -2))
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
                &empty_vector[0],
                &empty_vector[0],
                this->spiral_escape_dv,
                this->spiral_escape_thrust * 1000.0, //kN to N conversion
                this->spiral_escape_Isp,
                this->spiral_escape_power,
                this->spiral_escape_mdot,
                this->spiral_escape_number_of_engines,
                this->spiral_escape_active_power);
        }

        //print the departure/flyby
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
            phase_start_epoch / 86400.0,
            event_type,
            boundary1_name,
            0,
            ( (p > 0 || options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4) ? flyby_altitude : Bradius),
            (Btheta),
            ( (p > 0 || options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4) ? flyby_turn_angle : -1),
            RA_departure,
            DEC_departure,
            C3_departure,
            state_at_beginning_of_phase,
            dVdeparture,
            &empty_vector[0],
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
        if (this->detect_initial_coast)
        {
            //if this is a launch AND we are doing a forced post-launch initial coast
            
            if (j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6)
            {
                //initial coast after launch
                initial_coast_duration = options->forced_post_launch_coast;
            }
            else
            {
                //initial coast after flyby
                initial_coast_duration = options->forced_flyby_coast;
            }



            //initialize the propagator state
            for (size_t i = 7; i < 11; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_propagate[i] = 1.0;
            
			std::fill(augmented_state_at_initial_coast_midpoint.begin(), augmented_state_at_initial_coast_midpoint.end(), 0.0);

			//these are just dummy variables for the output
			dcurrent_epochdTOF[0] = 0.0;
			dcurrent_epochdTOF[1] = 1.0;
			dsegment_timedTOF = 0.0;
            
            resumeH = initial_coast_duration * 86400 / 2.0 / Universe->TU;

            integrator->adaptive_step_int(spacecraft_state_propagate,
											dummy_state_TOF_derivatives,
                augmented_state_at_initial_coast_midpoint,
											dummy_state_TOF_derivatives,
											empty_vector,
                (phase_start_epoch) / Universe->TU,
											dcurrent_epochdTOF,
                launchdate,
                initial_coast_duration / 2.0 / Universe->TU,
											dsegment_timedTOF,
                &resumeH,
                &resumeError,
                options->integrator_tolerance,
                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                this->STMrows,
                this->STMcolumns,
                (void*)options,
                (void*)Universe,
                DummyControllerPointer);



            //initialize the propagator state
            for (size_t i = 7; i < 11; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_propagate[i] = 1.0;

			dummy_state_TOF_derivatives.assign_zeros();

            //propagate forward to the end of the initial coast
            integrator->adaptive_step_int(augmented_state_at_initial_coast_midpoint,
											dummy_state_TOF_derivatives,
                spacecraft_state_propagate,
											dummy_state_TOF_derivatives,
                empty_vector,
                (phase_start_epoch) / Universe->TU,
											dcurrent_epochdTOF,
                launchdate,
                initial_coast_duration / 2.0 / Universe->TU,
											dsegment_timedTOF,
                &resumeH,
                &resumeError,
                options->integrator_tolerance,
                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                this->STMrows,
                this->STMcolumns,
                (void*)options,
                (void*)Universe,
                DummyControllerPointer);

            for (int k = 0; k < 3; ++k)
            {
                state_at_initial_coast_midpoint[k] = augmented_state_at_initial_coast_midpoint[k] * Universe->LU;
                state_at_initial_coast_midpoint[k + 3] = augmented_state_at_initial_coast_midpoint[k + 3] * Universe->LU / Universe->TU;
            }
            state_at_initial_coast_midpoint[6] = augmented_state_at_initial_coast_midpoint[6] * options->maximum_mass;

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
                &empty_vector[0],
                &empty_vector[0],
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
            //initialize the propagator state
            for (size_t i = 7; i < 11; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_propagate[i] = 1.0;

			dummy_state_TOF_derivatives.assign_zeros();

            resumeH = time_step_sizes[0] / 2.0 / Universe->TU;

            //first propagate to get the mid-point state of each step
            integrator->adaptive_step_int(spacecraft_state_propagate,
				                            dummy_state_TOF_derivatives,
                spacecraft_state_propagate_next,
				                            dummy_state_TOF_derivatives,
                                            this->control[step],
                (phase_start_epoch + phase_time_elapsed) / Universe->TU,
				                            dcurrent_epochdTOF,
                launchdate,
                this->time_step_sizes[step] / 2.0 / Universe->TU,
				                            dsegment_timedTOF,
                &resumeH,
                &resumeError,
                options->integrator_tolerance,
                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                this->STMrows,
                this->STMcolumns,
                (void*)options,
                (void*)Universe,
                DummyControllerPointer);

            //save the midpoint state
            for (size_t state = 0; state < 7; ++state)
                this->spacecraft_state[step][state] = spacecraft_state_propagate_next[state];

            //initialize the propagator state
            for (size_t i = 7; i < 11; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_propagate[i] = 1.0;

			dummy_state_TOF_derivatives.assign_zeros();

            //then propagate to the end of the step
            integrator->adaptive_step_int(spacecraft_state_propagate_next,
											dummy_state_TOF_derivatives,
                spacecraft_state_propagate,
											dummy_state_TOF_derivatives,
                                            this->control[step],
                (phase_start_epoch + phase_time_elapsed + 0.5 * this->time_step_sizes[step]) / Universe->TU,
											dcurrent_epochdTOF,
                launchdate,
                this->time_step_sizes[step] / 2.0 / Universe->TU,
											dsegment_timedTOF,
                &resumeH,
                &resumeError,
                options->integrator_tolerance,
                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                this->STMrows,
                this->STMcolumns,
                (void*)options,
                (void*)Universe,
                DummyControllerPointer);

            double angle1, angle2;
            if (step >= (options->num_timesteps / 2))
            {
                angle1 = atan2(control[step][1] + math::SMALL, control[step][0]) * EMTG::math::PI / 180.0;
                angle2 = asin(control[step][2] / EMTG::math::norm(control[step].data(), 3)) * EMTG::math::PI / 180.0;
            }
            else
            {
                angle1 = atan2(-control[step][1] + math::SMALL, -control[step][0]) * EMTG::math::PI / 180.0;
                angle2 = asin(-control[step][2] / EMTG::math::norm(control[step].data(), 3)) * EMTG::math::PI / 180.0;
            }

            //get the power and propulsion parameters at this half-step
            double dTdP, dmdotdP, dTdIsp, dmdotdIsp, dPdr, dPdt;
            EMTG::Astrodynamics::find_engine_parameters(options, EMTG::math::norm(spacecraft_state[step].data(), 3),
                event_epochs[step] - launchdate,
                &this->available_thrust[step],
                &this->available_mass_flow_rate[step],
                &this->available_Isp[step],
                &this->available_power[step],
                &this->active_power[step],
                &this->number_of_active_engines[step],
                false,
                &dTdP,
                &dmdotdP,
                &dTdIsp,
                &dmdotdIsp,
                &dPdr,
                &dPdt);

            double thrust_vector[3];
            for (int k = 0; k < 3; ++k)
                thrust_vector[k] = control[step][k] * available_thrust[step] * 1000.0; //kN to N conversion

            if (EMTG::math::norm(control[step].data(), 3) > 1.0e-2 && fabs(available_thrust[step]) > 1.0e-6)
            {
                event_type = "FBLTthrust";
            }
            else
            {
                event_type = "coast";
            }

            //scale from LU, LU/TU, MU to km, km/s, kg
            double scaled_state[7];
            for (int k = 0; k < 3; ++k)
            {
                scaled_state[k] = spacecraft_state[step][k] * Universe->LU;
                scaled_state[k + 3] = spacecraft_state[step][k + 3] * Universe->LU / Universe->TU;
            }
            scaled_state[6] = spacecraft_state[step][6] * options->maximum_mass;

            write_summary_line(options,
                Universe,
                eventcount,
                (phase_start_epoch + phase_time_elapsed + 0.5 * time_step_sizes[step]) / 86400.0,
                event_type,
                "deep-space",
                time_step_sizes[step] / 86400.0,
                -1,
                -1,
                -1,
                angle1,
                angle2,
                0,
                scaled_state,
                this->control[step].data(),
                thrust_vector,
                EMTG::math::norm(this->control[step].data(), 3),
                this->available_thrust[step] * 1000.0, //kN to N conversion
                this->available_Isp[step],
                this->available_power[step],
                math::norm(this->control[step].data(), 3) * this->available_mass_flow_rate[step],
                this->number_of_active_engines[step],
                this->active_power[step]);

            phase_time_elapsed += time_step_sizes[step];

            //if we have stepped halfway through, insert the match point line
            if (step == options->num_timesteps / 2 - 1)
                write_summary_line(options,
                Universe,
                eventcount,
                (phase_start_epoch + phase_time_elapsed) / 86400.0,
                "match_point",
                "deep-space",
                0.0,
                -1,
                -1,
                -1,
                0,
                0,
                0,
                match_point_state.data(),
                &empty_vector[0],
                &empty_vector[0],
                0,
                0,
                0,
                0,
                0,
                0,
                0);
        }

        //*****************************************************************************
        //next, if there was an terminal coast, we must print it
        if (this->detect_terminal_coast)
        {
			std::fill(augmented_state_at_terminal_coast_midpoint.begin(), augmented_state_at_terminal_coast_midpoint.end(), 0.0);

            //initial coast after flyby
            terminal_coast_duration = options->forced_flyby_coast;

            double resumeH = -terminal_coast_duration / 2.0 / Universe->TU;
            double resumeError = 1.0e-13;

			dummy_state_TOF_derivatives.assign_zeros();

            //The terminal coast STM entries of the state vector must be initialized to the identity
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;

            integrator->adaptive_step_int(spacecraft_state_propagate,
											dummy_state_TOF_derivatives,
                augmented_state_at_terminal_coast_midpoint,
											dummy_state_TOF_derivatives,
                                            empty_vector,
                (phase_end_epoch + phase_time_elapsed) / Universe->TU,
											dcurrent_epochdTOF,
                launchdate,
                terminal_coast_duration / 2.0 / Universe->TU,
											dsegment_timedTOF,
                &resumeH,
                &resumeError,
                options->integrator_tolerance,
                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                &temp_thrust,
                &temp_mdot,
                &temp_Isp,
                &temp_power,
                &temp_active_power,
                &temp_active_thrusters,
                this->STMrows,
                this->STMcolumns,
                (void*)options,
                (void*)Universe,
                DummyControllerPointer);

            for (int k = 0; k < 3; ++k)
            {
                state_at_terminal_coast_midpoint[k] = augmented_state_at_terminal_coast_midpoint[k] * Universe->LU;
                state_at_terminal_coast_midpoint[k + 3] = augmented_state_at_terminal_coast_midpoint[k + 3] * Universe->LU / Universe->TU;
            }
            state_at_terminal_coast_midpoint[6] = augmented_state_at_terminal_coast_midpoint[6] * options->maximum_mass;

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
                (phase_start_epoch + phase_time_elapsed + 0.5 * terminal_coast_duration) / 86400.0,
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
                &empty_vector[0],
                &empty_vector[0],
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
            else if (options->journey_arrival_type[j] == 3)
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
                dV_arrival_mag = sqrt(C3_arrival);
            }
            else if (options->journey_arrival_type[j] == 4 || options->journey_arrival_type[j] == 3)
            {
                dV_arrival_mag = 0;
                dVarrival[0] = 0;
                dVarrival[1] = 0;
                dVarrival[2] = 0;
                this->RA_arrival = 0.0;
                this->DEC_arrival = 0.0;
            }
            else
            {
                dV_arrival_mag = dV_arrival_magnitude;
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
                (phase_start_epoch + TOF) / 86400.0,
                event_type,
                boundary2_name,
                0,
                -1,
                -1,
                -1,
                RA_arrival,
                DEC_arrival,
                C3_arrival,
                state_at_end_of_phase,
                dVarrival,
                &empty_vector[0],
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
                    this->spiral_capture_time / 86400.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    this->spiral_capture_state_after_spiral,
                    &empty_vector[0],
                    &empty_vector[0],
                    this->spiral_capture_dv,
                    this->spiral_capture_thrust * 1000.0, //kN to N conversion
                    this->spiral_capture_Isp,
                    this->spiral_capture_power,
                    this->spiral_capture_mdot,
                    this->spiral_capture_number_of_engines,
                    this->spiral_capture_active_power);
            }
        }
    }


    //function to write a forward-integrated ephemeris file
    void FBLT_phase::write_ephemeris_file(  const missionoptions& options,
                                            const EMTG::Astrodynamics::universe& Universe,
                                            const double& launch_epoch,
                                            const double& journey_starting_epoch,
                                            vector< vector<string> >& output_line_array,
                                            const int& j,
                                            const int& p)
    {
        //Step 1: get the initial state of the spacecraft and write it out
        //recall that STK wants m and m/s, not km and km/s
        vector<string> output_line(7);
        vector<double> current_state(11*11+11);
        double current_epoch;

        //The STM entries of the current_state vector must be initialized to the identity
        for (size_t i = this->STMrows; i < this->num_states; ++i)
            current_state[i] = 0.0;
        for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
            current_state[i] = 1.0;

        //Step 1a: if there is a spiral
        if (p == 0 && options.journey_departure_type[j] == 5)
        {
            current_epoch = this->phase_start_epoch - this->spiral_escape_time;
            output_line[0] = string_utilities::convert_number_to_formatted_string(current_epoch - journey_starting_epoch, 2);
            for (size_t k = 0; k < 6; ++k)
            {
                output_line[k + 1] = string_utilities::convert_number_to_formatted_string(this->spiral_escape_state_before_spiral[k] * 1000.0, 2);
                current_state[k] = this->spiral_escape_state_before_spiral[k];
            }               
            current_state[6] = this->spiral_escape_state_before_spiral[6];
        }
        //Step 1b if there is not a spiral
        else
        {
            current_epoch = this->phase_start_epoch;
            output_line[0] = string_utilities::convert_number_to_formatted_string(current_epoch - journey_starting_epoch, 2);
            for (size_t k = 0; k < 6; ++k)
            {
                output_line[k + 1] = string_utilities::convert_number_to_formatted_string(this->state_at_beginning_of_phase[k] * 1000.0, 2);
                current_state[k] = this->state_at_beginning_of_phase[k];
            }    
            current_state[6] = this->state_at_beginning_of_phase[6];
        }
        output_line_array.push_back(output_line);

        //Step 2a: if this phase starts with a spiral, propagate the spacecraft with the planet until it leaves
        if (p == 0 && options.journey_departure_type[j] == 5)
        {
            double temp_body_state[6];
            double time_remaining_in_spiral = this->spiral_escape_time;

            while (time_remaining_in_spiral > 0.0)
            {
                if (time_remaining_in_spiral > 86400.0)
                {
                    time_remaining_in_spiral -= 86400.0;
                    current_epoch += 86400.0;
                }
                else
                {
                    current_epoch += time_remaining_in_spiral;
                    time_remaining_in_spiral = 0.0;
                }
                this->Body1->locate_body(current_epoch,
                                        temp_body_state,
                                        false,
                                        (missionoptions*)&options);
                
                output_line[0] = string_utilities::convert_number_to_formatted_string(current_epoch, 2);
                for (size_t k = 0; k < 6; ++k)
                {
                    output_line[k + 1] = string_utilities::convert_number_to_formatted_string(temp_body_state[k] * 1000.0, 2);
                    current_state[k] = temp_body_state[k];
                }
                current_state[6] = this->state_at_beginning_of_phase[6];

                output_line_array.push_back(output_line);
            }
        }
        //Step 2b: alternatively, if this phase starts with a forced post-launch coast, propagate along that coast
        else if (j == 0 && p == 0 && options.forced_post_launch_coast > 0.0)
        {
            for (size_t k = 0; k < 7; ++k)
                current_state[k] = state_at_beginning_of_phase[k];

            this->propagate_forward_ephemeris(  options,
                                                Universe,
                                                launch_epoch,
                                                current_epoch,
                                                current_state,
                                                -1,
                                                output_line_array,
                                                j,
                                                p,
                                                options.forced_post_launch_coast,
                                                journey_starting_epoch);
                                                
        }
        //Step 2c: if this phase starts with a post-flyby coast, propagate along that coast
        else if ((p > 0 || (p == 0 && (options.journey_departure_type[j] == 3 || options.journey_departure_type[j] == 4 || options.journey_departure_type[j] == 6)))
                 && options.forced_flyby_coast > 0.0)
        {
            for (size_t k = 0; k < 7; ++k)
                current_state[k] = state_at_beginning_of_phase[k];

            this->propagate_forward_ephemeris(  options,
                                                Universe,
                                                launch_epoch,
                                                current_epoch,
                                                current_state,
                                                -1,
                                                output_line_array,
                                                j,
                                                p,
                                                options.forced_flyby_coast,
                                                journey_starting_epoch);
        }

        //Step 3: propagate along the thrust/coast arcs
        for (size_t step = 0; step < options.num_timesteps; ++step)
        {
            this->propagate_forward_ephemeris(  options,
                                                Universe,
                                                launch_epoch,
                                                current_epoch,
                                                current_state,
                                                step,
                                                output_line_array,
                                                j,
                                                p,
                                                this->time_step_sizes[step],
                                                journey_starting_epoch);
        }

        //Step 4a: if this phase ends with a spiral, propagate the spacecraft with the planet until the end of the spiral
        if (p == options.number_of_phases[j] - 1 && options.journey_arrival_type[j] == 7)
        {
            double temp_body_state[6];
            double time_remaining_in_spiral = this->spiral_capture_time;

            while (time_remaining_in_spiral > 0.0)
            {
                if (time_remaining_in_spiral > 86400.0)
                {
                    time_remaining_in_spiral -= 86400.0;
                    current_epoch += 86400.0;
                }
                else
                {
                    current_epoch += time_remaining_in_spiral;
                    time_remaining_in_spiral = 0.0;
                }
                this->Body2->locate_body(current_epoch,
                                        temp_body_state,
                                        false,
                                        (missionoptions*)&options);

                output_line[0] = string_utilities::convert_number_to_formatted_string(current_epoch, 2);
                for (size_t k = 0; k < 6; ++k)
                {
                    output_line[k + 1] = string_utilities::convert_number_to_formatted_string(temp_body_state[k] * 1000.0, 2);
                    current_state[k] = temp_body_state[k];
                }

                output_line_array.push_back(output_line);
            }
        }
        //Step 4b: if this phase ends with a pre-flyby coast, propagate along that coast
        if (p < options.number_of_phases[j] - 1
            || (j < options.number_of_journeys - 1 && p == options.number_of_phases[j] - 1)
            && (options.journey_arrival_type[j] == 2 || options.journey_arrival_type[j] == 5))
        {
            this->propagate_forward_ephemeris(  options,
                                                Universe,
                                                launch_epoch,
                                                current_epoch,
                                                current_state,
                                                -1,
                                                output_line_array,
                                                j,
                                                p,
                                                options.forced_flyby_coast,
                                                journey_starting_epoch);
        }

        std::cout << "Forward propagation error for journey " << j << ", phase " << p << std::endl;
        std::cout << "t: " << this->phase_end_epoch - current_epoch << std::endl;
        std::cout << "x: " << this->state_at_end_of_phase[0] - current_state[0] << std::endl;
        std::cout << "y: " << this->state_at_end_of_phase[1] - current_state[1] << std::endl;
        std::cout << "z: " << this->state_at_end_of_phase[2] - current_state[2] << std::endl;
        std::cout << "vx: " << this->state_at_end_of_phase[3] - current_state[3] << std::endl;
        std::cout << "vy: " << this->state_at_end_of_phase[4] - current_state[4] << std::endl;
        std::cout << "vz: " << this->state_at_end_of_phase[5] - current_state[5] << std::endl;
        std::cout << "m: " << this->state_at_end_of_phase[6] - current_state[6] << std::endl;
            
    }

    //method to forward propagate a state for making a daily ephemeris
    void FBLT_phase::propagate_forward_ephemeris(const missionoptions& options,
                                                const EMTG::Astrodynamics::universe& Universe,
                                                const double& launch_epoch,
                                                double& current_epoch,
                                                vector<double>& current_state,
                                                const int& control_step,
                                                vector< vector<string> >& output_line_array,
                                                const int& j,
                                                const int& p,
                                                const double& propagation_time,
                                                const double& journey_starting_epoch)
    {
        vector<string> output_line(7);
        double time_remaining = propagation_time;
        std::vector <double> empty_control (3, 0.0);
        std::vector <double> temp_state(11 * 11 + 11);
		double dsegment_timedTOF = 0.0;

        double available_thrust, available_mass_flow_rate, available_Isp, available_power, active_power;
        int number_of_active_engines;

        for (size_t k = 0; k < 3; ++k)
        {
            current_state[k] /= Universe.LU;
            current_state[k+3] *= Universe.TU / Universe.LU;
        }
        current_state[6] /= options.maximum_mass;

		dummy_state_TOF_derivatives.assign_zeros();

        while (time_remaining > 0.0)
        {
            double step_time = time_remaining > 86400.0 ? 86400.0 : time_remaining;
            double resumeError = options.integrator_tolerance;
            double resumeH = step_time / Universe.TU;

			//I DON'T THINK THIS WILL WORK IN IT'S CURRENT FORM....CURRENT STATE MUST BE LENGTH 132
            this->integrator->adaptive_step_int(current_state,
												  dummy_state_TOF_derivatives,
                                                temp_state,
												  dummy_state_TOF_derivatives,
                                                  control_step >= 0 ? this->control[control_step] : empty_control,
                                                  current_epoch / Universe.TU,
												  dcurrent_epochdTOF,
                                                (current_epoch - launch_epoch) / Universe.TU,
                                                step_time / Universe.TU,
												  dsegment_timedTOF,
                                                &resumeH,
                                                &resumeError,
                                                options.integrator_tolerance,
                                                EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                                &available_thrust,
                                                &available_mass_flow_rate,
                                                &available_Isp,
                                                &available_power,
                                                &active_power,
                                                &number_of_active_engines,
												this->STMrows,
												this->STMcolumns,
                                                (void*)&options,
                                                (void*)&Universe,
                                                this->DummyControllerPointer);

            current_epoch += step_time;
            time_remaining -= step_time;

            //copy temp_state to current_state
            for (size_t k = 0; k < 11*11+11; ++k)
                current_state[k] = temp_state[k];

            //create a text output line
            output_line[0] = string_utilities::convert_number_to_formatted_string(current_epoch - journey_starting_epoch, 2);
            for (size_t k = 0; k < 3; ++k)
            {
                output_line[k + 1] = string_utilities::convert_number_to_formatted_string(current_state[k] * Universe.LU * 1000.0, 2);
                output_line[k + 4] = string_utilities::convert_number_to_formatted_string(current_state[k + 3] * Universe.LU / Universe.TU * 1000.0, 2);
            }

            output_line_array.push_back(output_line);
        }

        //re-scale current_state so it can be used on the other end
        for (size_t k = 0; k < 3; ++k)
        {
            current_state[k] *= Universe.LU;
            current_state[k + 3] *= Universe.LU / Universe.TU;
        }
        current_state[6] *= options.maximum_mass;
    }


	//function to calculate the match point derivatives
	void FBLT_phase::calculate_match_point_derivatives(double* G,
		int* Gindex,
		const int& j,
		const int& p,
		missionoptions* options,
		EMTG::Astrodynamics::universe* Universe)
	{
        //build the derivatives matrix through successive STM multiplication
        //every time step we add one more STM to the chain
        //this is symmetric for forward and backward propagation
        //and before we can do that, we need to create a vector of STMs with the control entries zeroed out
        std::vector< EMTG::math::Matrix <double> > forward_stripped_STMs = this->STM_archive_forward;
        std::vector< EMTG::math::Matrix <double> > backward_stripped_STMs = this->STM_archive_backward;
        static EMTG::math::Matrix<double> forward_cumulative_stripped_STM(this->STMrows, this->STMcolumns);
        static EMTG::math::Matrix<double> backward_cumulative_stripped_STM(this->STMrows, this->STMcolumns);
        forward_cumulative_stripped_STM.assign_zeros();
        backward_cumulative_stripped_STM.assign_zeros();

        //initialize
        for (size_t i = 0; i < this->STMrows; ++i)
        {
            forward_cumulative_stripped_STM(i, i) = 1.0;
            backward_cumulative_stripped_STM(i, i) = 1.0;
        }

        //step through and create STM chains
        for (int step = options->num_timesteps / 2 - 1; step >= 0; --step)
        {
            //create a stripped version of this step's STM
            for (int row = 0; row < 7; ++row)
            {
                for (int column = 7; column < 10; ++column)
                {
                    forward_stripped_STMs[step](row, column) = 0.0;
                    backward_stripped_STMs[step](row, column) = 0.0;
                }
            }

            //compute the cumulative STM for this step
            if (step < options->num_timesteps / 2 - 1)
            {
                forward_cumulative_stripped_STM *= forward_stripped_STMs[step + 1];
                backward_cumulative_stripped_STM *= backward_stripped_STMs[step + 1];
            }

            //multiply in this step's STM
            this->forward_cumulative_STM_archive[step] = forward_cumulative_stripped_STM * STM_archive_forward[step];
            this->backward_cumulative_STM_archive[step] = backward_cumulative_stripped_STM * STM_archive_backward[step];
        }
        //advanced the cumulative stripped STMs all the way to the endpoints
        forward_cumulative_stripped_STM *= forward_stripped_STMs[0];
        backward_cumulative_stripped_STM *= backward_stripped_STMs[0];

        //create initial and terminal coast STMs if appropriate
        static EMTG::math::Matrix <double> cumulative_initial_coast_STM;
        static EMTG::math::Matrix <double> cumulative_terminal_coast_STM;
        if (detect_initial_coast)
        {
            //strip the initial coast STM because there is no control applied
            for (int row = 0; row < 7; ++row)
            {
                for (int column = 7; column < 10; ++column)
                {
                    this->initial_coast_STM(row, column) = 0.0;
                }
            }
            
            //construct the cumulative STM for the initial coast
            cumulative_initial_coast_STM = forward_cumulative_stripped_STM * this->initial_coast_STM;
        }

        if (detect_terminal_coast)
        {
            //strip the terminal coast STM because there is no control applied
            for (int row = 0; row < 7; ++row)
            {
                for (int column = 7; column < 10; ++column)
                {
                    terminal_coast_STM(row, column) = 0.0;
                }
            }

            //construct the cumulative STM for the terminal coast
            cumulative_terminal_coast_STM = backward_cumulative_stripped_STM * this->terminal_coast_STM;
        }

        
		//Match point constraint derivatives with respect to forward controls
        for (int step = options->num_timesteps / 2 - 1; step >= 0; --step)
		{
			//place the derivatives in the Jacobian
			//since the match point constraint is defined backward - forward, 
			//the forward derivatives will have a negative sign

			//match point with respect to control of the current step
			for (size_t cindex = 0; cindex < 3; ++cindex)
			{
				for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    G[match_point_constraint_G_indices[step + 1][stateindex][cindex]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[step + 1][stateindex][cindex]]] * forward_cumulative_STM_archive[step](stateindex, 7 + cindex);
			}
		}


			

		//derivative with respect to Isp for forward propagation with VSI thruster
		if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
		{
			 
		}

		if (options->objective_type == 13) //derivatives with respect to power for minimum power problems
		{
			
		}


        //derivatives with respect to phase departure from the left boundary
		if (p == 0)
		{
			if (j == 0 && options->allow_initial_mass_to_vary && !(options->journey_departure_type[j] == 5))
			{
                if (detect_initial_coast)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier[stateindex]]] * this->unscaled_phase_initial_mass / options->maximum_mass * cumulative_initial_coast_STM(stateindex, 6);
                    }
                }
                else
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier[stateindex]]] * this->unscaled_phase_initial_mass / options->maximum_mass * forward_cumulative_stripped_STM(stateindex, 6);
                    }
                }
			}
			if (options->journey_variable_mass_increment[j] && !(options->journey_departure_type[j] == 5))
			{
                if (detect_initial_coast)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier[stateindex]]] * options->journey_starting_mass_increment[j] / options->maximum_mass * cumulative_initial_coast_STM(stateindex, 6);
                    }
                }
                else
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier[stateindex]]] * options->journey_starting_mass_increment[j] / options->maximum_mass * forward_cumulative_stripped_STM(stateindex, 6);
                    }
                }
			}
			if (options->journey_departure_type[j] == 0)// && options->LV_type > 0)
			{
                
				//initial asymptote derivatives
                double cRA = cos(RA_departure);
                double sRA = sin(RA_departure);
                double cDEC = cos(DEC_departure);
                double sDEC = sin(DEC_departure);
                double v_infinity = sqrt(C3_departure);
                
                //derivatives with respect to v-infinity
                double dvx0_dvinf = cRA * cDEC;
                double dvy0_dvinf = sRA * cDEC;
                double dvz0_dvinf = sDEC;

                //derivatives with respect to RA
                double dvx0_dRA = v_infinity * (-sRA*cDEC);
                double dvy0_dRA = v_infinity * (cRA*cDEC);
                double dvz0_dRA = 0.0;

                //derivatives with respect to DEC
                double dvx0_dDEC = v_infinity * (-cRA*sDEC);
                double dvy0_dDEC = v_infinity * (-sRA*sDEC);
                double dvz0_dDEC = v_infinity * (cDEC);
 
                if (detect_initial_coast)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[0][stateindex][0]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][0]]] * Universe->TU / Universe->LU
                            * (cumulative_initial_coast_STM(stateindex, 3) * dvx0_dvinf
                            + cumulative_initial_coast_STM(stateindex, 4) * dvy0_dvinf
                            + cumulative_initial_coast_STM(stateindex, 5) * dvz0_dvinf
                            + cumulative_initial_coast_STM(stateindex, 6) * this->dmdvinf * this->mission_initial_mass_multiplier * this->unscaled_phase_initial_mass / options->maximum_mass * Universe->TU / Universe->LU);

                        G[match_point_constraint_G_indices[0][stateindex][1]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][1]]] * Universe->TU / Universe->LU
                            * (cumulative_initial_coast_STM(stateindex, 3) * dvx0_dRA
                            + cumulative_initial_coast_STM(stateindex, 4) * dvy0_dRA
                            + cumulative_initial_coast_STM(stateindex, 5) * dvz0_dRA);

                        G[match_point_constraint_G_indices[0][stateindex][2]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][2]]] * Universe->TU / Universe->LU
                            * (cumulative_initial_coast_STM(stateindex, 3) * dvx0_dDEC
                            + cumulative_initial_coast_STM(stateindex, 4) * dvy0_dDEC
                            + cumulative_initial_coast_STM(stateindex, 5) * dvz0_dDEC);
                    }
                }   
                else
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[0][stateindex][0]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][0]]] * Universe->TU / Universe->LU
                            * (forward_cumulative_stripped_STM(stateindex, 3) * dvx0_dvinf
                            + forward_cumulative_stripped_STM(stateindex, 4) * dvy0_dvinf
                            + forward_cumulative_stripped_STM(stateindex, 5) * dvz0_dvinf
                            + forward_cumulative_stripped_STM(stateindex, 6) * this->dmdvinf * this->mission_initial_mass_multiplier * this->unscaled_phase_initial_mass / options->maximum_mass * Universe->TU / Universe->LU);

                        G[match_point_constraint_G_indices[0][stateindex][1]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][1]]] * Universe->TU / Universe->LU
                            * (forward_cumulative_stripped_STM(stateindex, 3) * dvx0_dRA
                            + forward_cumulative_stripped_STM(stateindex, 4) * dvy0_dRA
                            + forward_cumulative_stripped_STM(stateindex, 5) * dvz0_dRA);

                        G[match_point_constraint_G_indices[0][stateindex][2]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][2]]] * Universe->TU / Universe->LU
                            * (forward_cumulative_stripped_STM(stateindex, 3) * dvx0_dDEC
                            + forward_cumulative_stripped_STM(stateindex, 4) * dvy0_dDEC
                            + forward_cumulative_stripped_STM(stateindex, 5) * dvz0_dDEC);
                    }
                }

                
			}

			if (j > 0) //for successive journeys the mass at the beginning of the phase affects the following patch point
			{
				//this must be disabled for phases that start with spirals
				if (!(options->journey_arrival_type[j - 1] == 7))
				{
                    if (detect_initial_coast)
                    {
                        for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                        {
                            G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]]] * cumulative_initial_coast_STM(stateindex, 6) / options->maximum_mass;
                        }
                    }
                    else
                    {
                        for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                        {
                            G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]]] * forward_cumulative_stripped_STM(stateindex, 6) / options->maximum_mass;
                        }
                    }
				}
			}
		}
		else//for other phases other than the first, i.e. coming out of a flyby. Velocity is cartesian and do not forget initial mass!
		{
            if (detect_initial_coast)
            {
                //initial velocity
                for (size_t vindex = 0; vindex < 3; ++vindex)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[0][stateindex][vindex]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][vindex]]] * cumulative_initial_coast_STM(stateindex, vindex + 3) * Universe->TU / Universe->LU;
                    }
                }
                //initial mass
                for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                {
                    G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]]] * cumulative_initial_coast_STM(stateindex, 6) / options->maximum_mass;
                }
            }
            else
            {
                //initial velocity
                for (size_t vindex = 0; vindex < 3; ++vindex)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[0][stateindex][vindex]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][stateindex][vindex]]] * forward_cumulative_stripped_STM(stateindex, vindex + 3) * Universe->TU / Universe->LU;
                    }
                }
                //initial mass
                for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                {
                    G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]] = -options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[stateindex]]] * forward_cumulative_stripped_STM(stateindex, 6) / options->maximum_mass;
                }
            }
		}

		//derivatives of the constraints with respect to the forward flight time variables
		//this in turn is dependent on the velocity and acceleration of the left-most body
		if (options->derivative_type > 2)
		{
            for (int stateindex = 0; stateindex < 7; ++stateindex)
            {
				for (int timevar = 0; timevar < this->G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables[0].size() - 1; ++timevar)
					G[this->G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex][timevar]] = -X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex][timevar] * this->dspacecraft_state_forwarddTOF(stateindex, 1) / Universe->TU;				
            }


			//the following derivatives are for the current phase flight time ONLY
			if (options->derivative_type > 3)
			{
				for (int stateindex = 0; stateindex < 7; ++stateindex)
				{
					G[this->G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex].back()] = -X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex].back() * this->dspacecraft_state_forwarddTOF(stateindex, 0) / Universe->TU;
				}
			}
		}
		
		//compute and store the backward derivatives of the match point constraints with respect to the control unit vector
		//and (if applicable) variable Isp
		//build the derivatives matrix through successive STM multiplication
		//every time step we add one more STM to the chain
        for (int step = 0; step < options->num_timesteps / 2; ++step)
		{
			//place the derivatives in the Jacobian
			//since the match point constraint is defined backward - forward, 
			//the backward derivatives will have a positive sign

			//match point with respect to u_x of the current step
            //match point with respect to control of the current step
            for (size_t cindex = 0; cindex < 3; ++cindex)
            {
                for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                {
                    G[match_point_constraint_G_indices[options->num_timesteps - step][stateindex][cindex]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[options->num_timesteps - step][stateindex][cindex]]] * backward_cumulative_STM_archive[step](stateindex, 7 + cindex);
                }
            }
		}

		//derivative with respect to Isp for backward propagation with VSI thruster
		if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
		{
			for (int step = 0; step < options->num_timesteps / 2; ++step)
			{
				//translate into backward steps
				int backstep = options->num_timesteps - 1 - step;

				
			}
		}
		 //end backward stuff

		//derivative of match point constraints with respect to arrival mass
		//disable this for a phase ending in a capture spiral
		if (!(p == options->number_of_phases[j] && options->journey_arrival_type[j] == 7))
		{
            if (detect_terminal_coast)
            {
                for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                {
                    G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass[stateindex]] = options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass[stateindex]]] * cumulative_terminal_coast_STM(stateindex, 6) / options->maximum_mass;
                }
            }
            else
            {
                for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                {
                    G[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass[stateindex]] = options->X_scale_ranges[options->jGvar[this->G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass[stateindex]]] * backward_cumulative_stripped_STM(stateindex, 6) / options->maximum_mass;
                }
            }
		}

		//derivatives of the constraints with respect to the backward flight time variables
		//this in turn is dependent on the velocity and acceleration of the right-most body
		if (options->derivative_type > 2)
		{
            for (int stateindex = 0; stateindex < 7; ++stateindex)
            {
                for (int timevar = 0; timevar < this->G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables[0].size() - 1; ++timevar)
                    G[this->G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex][timevar]] = -X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex][timevar] * this->dspacecraft_state_backwarddTOF(stateindex, 1) / Universe->TU;
            }

			//the following derivatives are for the current phase flight time
			if (options->derivative_type > 3)
			{
                for (int stateindex = 0; stateindex < 7; ++stateindex)
                    G[this->G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex].back()] = -X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables[stateindex].back() * this->dspacecraft_state_backwarddTOF(stateindex, 0) / Universe->TU;
			}
		}

        
        // derivative with respect to arrival velocity
        //only evaluated for phases that are not terminal intercepts
        if (!(p == options->number_of_phases[j] - 1 && ((options->journey_arrival_type[j] == 1) || options->journey_arrival_type[j] == 3) || options->journey_arrival_type[j] == 5 || options->journey_arrival_type[j] == 7))
        {
            if (detect_terminal_coast)
            {
                for (size_t vindex = 0; vindex < 3; ++vindex)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]]] * cumulative_terminal_coast_STM(stateindex, vindex + 3) * Universe->TU / Universe->LU;
                    }
                }
            }
            else
            {
                for (size_t vindex = 0; vindex < 3; ++vindex)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]]] * backward_cumulative_stripped_STM(stateindex, vindex + 3) * Universe->TU / Universe->LU;
                    }
                }
            }
        }//end code with respect to arrival velocity
        
	}

} /* namespace EMTG */
