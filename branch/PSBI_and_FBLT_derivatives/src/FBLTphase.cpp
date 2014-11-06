/*
 * FBLTphase.cpp
 *
 *  Created on: September 17, 2012
 *      Author: Jacob
 */
#include <sstream>
#include <fstream>
#include <iostream>

#include "FBLTphase.h"
#include "Astrodynamics.h"
#include "missionoptions.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "equations_of_motion.h"
#include "universe.h"
#include "EMTG_string_utilities.h"

#include "SpiceUsr.h"



namespace EMTG {

FBLT_phase::FBLT_phase() {
//default constructor does nothing

}

    FBLT_phase::FBLT_phase(int j, int p, missionoptions* options)
    {
	    //must resize all data vectors to the correct length
		vector<double> state_dummy(11 + 11*11);
	    vector<double> dV_or_control_dummy(3);

	    for (int step = 0; step < options->num_timesteps; ++step) {
		    spacecraft_state.push_back(state_dummy);
		    control.push_back(dV_or_control_dummy);
	    }

	    match_point_state.resize(7);


	    event_epochs.resize(options->num_timesteps);
	    available_power.resize(options->num_timesteps);
	    available_thrust.resize(options->num_timesteps);
	    available_mass_flow_rate.resize(options->num_timesteps);
	    available_Isp.resize(options->num_timesteps);
	    active_power.resize(options->num_timesteps);
	    number_of_active_engines.resize(options->num_timesteps);

        //vector to track the state and derivatives of the central body
        vector<double> central_body_state_dummy(options->derivative_type > 2 ? 12 : 6);
        for (size_t step = 0; step < options->num_timesteps; ++step)
            this->central_body_state_mks.push_back(central_body_state_dummy);

	    //set the bodies
	    boundary1_location_code = options->sequence[j][p];
	    boundary2_location_code = options->sequence[j][p+1];

	    //size the vectors that will be used to calculate the b-plane
	    V_infinity_in.resize(3, 1);
	    V_infinity_out.resize(3, 1);
	    BoundaryR.resize(3, 1);
	    BoundaryV.resize(3, 1);

	    //set up the integrator
		//for analytical FBLT derivatives, we need to integrate STM
		//entries as states as well, so we must instantiate the integrator 
		//with the normal s/c states but also slots for the STM entries
		this->STMrows = 11;
		this->STMcolumns = 11;
		this->num_states = 11 + 11 * 11;
	    integrator = new EMTG::integration::rk8713M(num_states);

	    current_mass_increment = 0.0;
	    journey_initial_mass_increment_scale_factor = 1.0;

	    //size the time step vector
	    time_step_sizes.resize(options->num_timesteps);

	    //set derivatives for spirals
	    this->spiral_escape_dm_after_dm_before = 1.0;

		//set up STMS
		for (size_t step = 0; step < options->num_timesteps / 2; ++step)
		{
			this->STM_archive_forward.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
			this->STM_archive_backward.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
            this->forward_cumulative_STM_archive.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
            this->backward_cumulative_STM_archive.push_back(EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0));
		}
        if ((j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6) || ((p > 0 || p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4 || options->journey_departure_type[j] == 6)) && options->forced_flyby_coast > 1.0e-6))
        {

            this->detect_initial_coast = true;
            this->initial_coast_STM = EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0);
        }
        else
            this->detect_initial_coast = false;
        
        if ((p < options->number_of_phases[j] - 1 || (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5)) && options->forced_flyby_coast > 1.0e-6)
        {
            this->detect_terminal_coast = true;
            this->terminal_coast_STM = EMTG::math::Matrix< double >(this->STMrows, this->STMcolumns, 0.0);
        }
        else
            this->detect_terminal_coast = false;

    }

    FBLT_phase::~FBLT_phase() {
	    //delete the integrator object
		    delete integrator;
    }

    //evaluate function
    //return 0 if successful, 1 if failure
    int FBLT_phase::evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
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
        if ((p < options->number_of_phases[j] - 1 || (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5)) && options->forced_flyby_coast > 1.0e-6)
            total_available_thrust_time -= options->forced_flyby_coast;

        for (int step = 0; step < options->num_timesteps; ++step)
        {
            time_step_sizes[step] *= total_available_thrust_time / step_size_normalization_coefficient;
        }

        //Step 6.2: propagate forward
        phase_time_elapsed_forward = 0.0;
        //store the initial prefered integration step size
        double resumeH = time_step_sizes[0] / 2.0 / Universe->TU;
        double resumeError = 1.0e-13;

        //first initialize the forward integration
        //the following array holds the spacecraft state at "half steps," i.e. halfway through each integration segment

        //WE MUST CHANGE THIS INTO A C++ VECTOR!!
        double spacecraft_state_forward[11 + 11 * 11];
        double spacecraft_state_forward_prop[11 + 11 * 11];
        for (int k = 0; k < 7; ++k)
            spacecraft_state_forward[k] = state_at_beginning_of_phase[k];

        //in the FBLT STM test bed I passed ux,uy,uz and TOF around as part of the state
        //to keep things consistent, let's just set these to zero for all time
        //we will know something is wrong if, in fact, they are ever not zero
        for (int k = 7; k < this->STMrows; ++k)
            spacecraft_state_forward[k] = 0.0;

        //we want to calculate the STM for every time step, but since we are propagating in half-steps let's do that for now
        //first make an archive to store all of the STMs
        //we also need to made a place to store the initial coast STMs if there is an initial coast


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

        //Step 6.2.0.1 if there is an initial coast, propagate through it
        if (detect_initial_coast)
        {
            //if this is a launch AND we are doing a forced post-launch initial coast
            double spacecraft_state_end_coast[11 + 11*11];
            double empty_vector[] = { 0.0, 0.0, 0.0 };
            double dummy_parameter = 0.0;

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


            integrator->adaptive_step_int(spacecraft_state_forward,
                                            spacecraft_state_end_coast,
                                            empty_vector,
                                            (phase_start_epoch) / Universe->TU,
                                            X[0],
                                            initial_coast_duration / Universe->TU,
                                            &resumeH,
                                            &resumeError,
                                            1.0e-8,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            (int*)&dummy_parameter,
                                            this->STMrows,
                                            this->STMcolumns,
                                            (void*)options,
                                            (void*)Universe,
                                            DummyControllerPointer);

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

            for (int k = 0; k < 7; ++k)
                spacecraft_state_forward[k] = spacecraft_state_end_coast[k];
        }

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


            /*
            //before we propagate with derivatives, let's do a "mini finite differencing" check
            double control_perturbation = 1.0e-6;
            double dstate_du[7][3]; //indexed as [state][control]
            vector<double> reference_control = this->control[step];
            double spacecraft_state_forward_prop_plus[11 + 11 * 11];
            double spacecraft_state_forward_prop_minus[11 + 11 * 11];
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
            spacecraft_state_forward_prop_plus,
            control_plus.data(),
            (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
            X[0],
            time_step_sizes[step] / Universe->TU,
            &resumeH,
            &resumeError,
            1.0e-8,
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
            spacecraft_state_forward_prop_minus,
            control_minus.data(),
            (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
            X[0],
            time_step_sizes[step] / Universe->TU,
            &resumeH,
            &resumeError,
            1.0e-8,
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
            */
            //now generate the actual STM
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_forward[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_forward[i] = 1.0;

            integrator->adaptive_step_int(spacecraft_state_forward,
                spacecraft_state_forward_prop,
                control[step].data(),
                (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
                X[0],
                time_step_sizes[step] / Universe->TU,
                &resumeH,
                &resumeError,
                1.0e-8,
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
            /*
            //finally, compare the propagated STM to the central differenced STM
            std::cout << "Difference between propagated STM and central differenced STM" << std::endl;
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
            */
            //getchar();

            //copy the state over
            for (size_t state = 0; state < 11; ++state)
                spacecraft_state_forward[state] = spacecraft_state_forward_prop[state];

            //step 6.2.4 encode the epoch of the step midpoint
            event_epochs[step] = phase_start_epoch + phase_time_elapsed_forward + 0.5 * time_step_sizes[step];
            phase_time_elapsed_forward += time_step_sizes[step];
        }

        /*
        //now let's repeat the integration inside a finite difference loop to check the multiplication of the STMs
        double control_perturbation = 1.0e-6;
        vector<double> temp_control;
        double dstate_du[7][3]; //indexed as [state][control]
        double spacecraft_state_forward_prop_plus[11 + 11 * 11];
        double spacecraft_state_forward_prop_plus_2[11 + 11 * 11];
        double spacecraft_state_forward_prop_minus[11 + 11 * 11];
        double spacecraft_state_forward_prop_minus_2[11 + 11 * 11];

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
                                                spacecraft_state_forward_prop_plus_2,
                                                temp_control.data(),
                                                (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
                                                X[0],
                                                time_step_sizes[step] / Universe->TU,
                                                &resumeH,
                                                &resumeError,
                                                1.0e-8,
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
                                                spacecraft_state_forward_prop_minus_2,
                                                temp_control.data(),
                                                (phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
                                                X[0],
                                                time_step_sizes[step] / Universe->TU,
                                                &resumeH,
                                                &resumeError,
                                                1.0e-8,
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
        */
        //end STM checker
        
	    //Step 6.3: propagate backward
	    phase_time_elapsed_backward = 0.0;
	    //store the initial prefered integration step size
	    resumeH = -time_step_sizes[options->num_timesteps - 1] / Universe->TU;
	    resumeError = 1.0e-6;
	
	    //first initialize the backward integration
	    double spacecraft_state_backward[11 + 11*11];
        double spacecraft_state_backward_prop[11 + 11 * 11];

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

	    //Step 6.3.0.1 if there is an terminal coast, propagate through it
	    if (this->detect_terminal_coast)
	    {
		    double spacecraft_state_end_coast[11 + 11*11];
		    double empty_vector[] = {0.0,0.0,0.0};
		    double dummy_parameter = 0.0;

		    //initial coast after flyby
		    terminal_coast_duration = options->forced_flyby_coast;

		    double resumeH = -terminal_coast_duration / 2.0 / Universe->TU;

			
			
			//The terminal coast STM entries of the state vector must be initialized to the identity
			for (size_t i = this->STMrows; i < this->num_states; ++i)
				spacecraft_state_backward[i] = 0.0;
			for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
				spacecraft_state_backward[i] = 1.0;


		    integrator->adaptive_step_int(	spacecraft_state_backward,
                                            spacecraft_state_end_coast,
										    empty_vector, 
										    (phase_end_epoch) / Universe->TU,
										    X[0],
										    -terminal_coast_duration / Universe->TU, 
										    &resumeH,
										    &resumeError,
										    1.0e-8,
										    EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
										    &dummy_parameter,
										    &dummy_parameter,
										    &dummy_parameter,
										    &dummy_parameter,
										    &dummy_parameter,
										    (int*) &dummy_parameter,
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
					terminal_coast_STM(i, j) = state_at_terminal_coast_midpoint[statecount];
					++statecount;
				}
			}

		    phase_time_elapsed_backward += terminal_coast_duration;

		    for (int k = 0; k < 7; ++k)
			    spacecraft_state_backward[k] = spacecraft_state_end_coast[k];
	    }

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
		    F[*Findex + (backstep - options->num_timesteps/2)] = throttle;


			//The STM entries of the state vector must be initialized to the identity before every step
			for (size_t i = this->STMrows; i < this->num_states; ++i)
				spacecraft_state_backward[i] = 0.0;
			for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
				spacecraft_state_backward[i] = 1.0;
		
		    //step 6.3.3 propagate the spacecraft to the midpoint of the step using the control unit vector
		    integrator->adaptive_step_int(	spacecraft_state_backward,
										    spacecraft_state_backward_prop,
										    control[backstep].data(),  
										    (phase_end_epoch - phase_time_elapsed_backward) / Universe->TU,
										    X[0],
										    -time_step_sizes[backstep] / Universe->TU, 
										    &resumeH,
										    &resumeError,
										    1.0e-8,
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

            //copy the state over
            for (size_t state = 0; state < 11; ++state)
                spacecraft_state_backward[state] = spacecraft_state_backward_prop[state];

		    //step 6.3.4 encode the epoch of the step midpoint
		    event_epochs[backstep] = phase_end_epoch - phase_time_elapsed_backward - 0.5 * time_step_sizes[backstep];
		    phase_time_elapsed_backward += time_step_sizes[backstep];
	    }

	    //step Xindex back to the end of the arc
	    if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
		    (*Xindex) += 4 * options->num_timesteps/2;
	    else
		    (*Xindex) += 3 * options->num_timesteps/2;

	    //step Findex back to the end of the arc
	    (*Findex) += options->num_timesteps/2;

	    //Step 6.4: enforce match point constraint
	    for (size_t k=0; k<3; ++k)
	    {
		    //position
            F[*Findex + k] = (spacecraft_state_backward_prop[k] - spacecraft_state_forward_prop[k]);
		
		    //velocity
            F[*Findex + k + 3] = (spacecraft_state_backward_prop[k + 3] - spacecraft_state_forward_prop[k + 3]);

		    //unscale the match point state (normalized to mks)
            match_point_state[k] = spacecraft_state_forward_prop[k] * Universe->LU;
            match_point_state[k + 3] = spacecraft_state_forward_prop[k + 3] * Universe->LU / Universe->TU;
	    }
	    //mass
        F[*Findex + 6] = (spacecraft_state_backward_prop[6] - spacecraft_state_forward_prop[6]);
	    (*Findex) += 7;

		//unscale the match point mass
        match_point_state[6] = spacecraft_state_forward_prop[6] * options->maximum_mass;


		//CALCULATE MATCH POINT DERIVATIVES HERE
		if (options->derivative_type > 1 && needG)
			this->calculate_match_point_derivatives(G, Gindex, j, p, STM_archive_forward, STM_archive_backward, options, Universe);

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
				    dV_arrival_magnitude = process_arrival(	state_at_end_of_phase+3,
														    boundary2_state,
														    current_state+3,
														    current_epoch,
														    Body2->mu,
														    Body2->r_SOI,
														    F,
														    Findex, 
														    j, 
														    options,
														    Universe);
			    else //arriving at a boundary point in free space
				    dV_arrival_magnitude = process_arrival(	state_at_end_of_phase + 3, 
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
			    state_at_end_of_phase[6] *= exp(-dV_arrival_magnitude * 1000/ ((options->IspChem > 0 ? options->IspChem : 1e-6)* options->g0));
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
			    /*if ( needG && options->journey_arrival_declination_constraint_flag[j] && (options->journey_arrival_type[j] == 0 || options->journey_arrival_type[j] == 2) ) //intercept with bounded v-infinity or orbit insertion
			    {
				    double Vx = X[arrival_declination_constraint_X_indices[0]];
				    double Vy = X[arrival_declination_constraint_X_indices[1]];
				    double Vz = X[arrival_declination_constraint_X_indices[2]];
				    double A = sqrt(Vx*Vx + Vy*Vy);
				    double B = Vx*Vx + Vy*Vy + Vz*Vz;

				    //derivative with respect to x component of terminal velocity
				    G[arrival_declination_constraint_G_indices[0]] = arrival_declination_constraint_X_scale_ranges[0] * Vx*Vz / (A*B) / options->journey_arrival_declination_bounds[j][1];

				    //derivative with respect to y component of terminal velocity
				    G[arrival_declination_constraint_G_indices[1]] = arrival_declination_constraint_X_scale_ranges[0] * Vy*Vz / (A*B) / options->journey_arrival_declination_bounds[j][1];

				    //derivative with respect to z component of terminal velocity
				    G[arrival_declination_constraint_G_indices[2]] = arrival_declination_constraint_X_scale_ranges[0] * -A/B / options->journey_arrival_declination_bounds[j][1];
			    }*/
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

	    return 0;
    }


    //bounds calculation function
    //return 0 if successful, 1 if failure
    void FBLT_phase::calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
    {
	    //this function calculates the upper and lower bounds for the decision and constraint vectors for FBLT
	    //create a prefix string with journey and phase information
	    stringstream prefixstream;
	    prefixstream << "j" << j << "p" << p << ": ";
	    string prefix = prefixstream.str();
	    int first_X_entry_in_phase = Xupperbounds->size();

	    //**************************************************************************
	    //calculate bounds on variables and constraints governing the left boundary
	    calcbounds_left_boundary(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	    //**************************************************************************
	    //if EMTG is choosing an input power or Isp for the phase (for REP/NEP models), then this information must be encoded
	    calcbounds_phase_thruster_parameters(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	    //**************************************************************************
	    //next, we need to encode the phase flight time
	    calcbounds_flight_time(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, synodic_periods, j, p, Universe, options);

	    //**************************************************************************
	    //calculate bounds on variables and constraints governing the right boundary
	    calcbounds_right_boundary(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	    //**************************************************************************
	    //if we are using a variable scale/standard deviation distribution to pick the control points
	    if (options->step_size_distribution == 3 || options->step_size_distribution == 4)
	    {
		    Xlowerbounds->push_back(1.0);
		    Xupperbounds->push_back(options->step_size_stdv_or_scale);
		    Xdescriptions->push_back(prefix + "step size distribution standard deviation or scale factor");
	    }

	    //**************************************************************************
	    //next, we need to include the decision variables and constraints for each burn
	    //now, for each timestep
	    if (options->control_coordinate_system == 0) //Cartesian control
	    {
		    for (int w = 0; w < options->num_timesteps; ++w)
		    {
			    stringstream stepstream;
			    stepstream << w;
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
			    if (options->engine_type == 4 || options->engine_type == 13)
			    {
				    Xlowerbounds->push_back(options->IspLT_minimum);
				    Xupperbounds->push_back(options->IspLT);
				    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
			    }
			    else if (options->engine_type == 12)
			    {
				    Xlowerbounds->push_back(3000.0);
				    Xupperbounds->push_back(5000.0);
				    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
			    }

			    //and the throttle magnitude constraint
			    //throttle = 0
			    Flowerbounds->push_back(0.0);
			    Fupperbounds->push_back(1.0);
			    Fdescriptions->push_back(prefix + "step " + stepstream.str() + " throttle magnitude constraint");

			    //Jacobian entries for the throttle magnitude constraint
			    //The throttle magnitude constraint is dependent only on the preceding three throttle components
			    vector<int> step_G_indices;
			    vector<double> dummyvalues(3);
			    int vary_Isp_flag = (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13) ? 1 : 0;
			    for (size_t entry = Xdescriptions->size() - 1 - vary_Isp_flag; entry > Xdescriptions->size() - 4 - vary_Isp_flag; --entry)
			    {
				    iGfun->push_back(Fdescriptions->size() - 1);
				    jGvar->push_back(entry);
				    stringstream EntryNameStream;
				    EntryNameStream << "Derivative of " << prefix + "step " << w << " throttle magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
				    Gdescriptions->push_back(EntryNameStream.str());

				    //store the position in the G vector
				    step_G_indices.push_back(iGfun->size() - 1);
			    }
			    control_vector_G_indices.push_back(step_G_indices);
		    }
	    }
	    else if (options->control_coordinate_system == 1) //Polar control
	    {
		    for (int w = 0; w < options->num_timesteps; ++w)
		    {
			    stringstream stepstream;
			    stepstream << w;
			    //u_x
			    Xlowerbounds->push_back(1.0e-10);
			    Xupperbounds->push_back(1.0);
			    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_Throttle");

			    //u_y
			    Xlowerbounds->push_back(0.0);
			    Xupperbounds->push_back(1.0);
			    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_u");

			    //u_z
			    Xlowerbounds->push_back(0.0);
			    Xupperbounds->push_back(1.0);
			    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_v");

			    //for variable specific impulse propulsion systems, we must also encode the Isp for this time step
			    if (options->engine_type == 4 || options->engine_type == 13)
			    {
				    Xlowerbounds->push_back(options->IspLT_minimum);
				    Xupperbounds->push_back(options->IspLT);
				    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
			    }
			    else if (options->engine_type == 12)
			    {
				    Xlowerbounds->push_back(3000.0);
				    Xupperbounds->push_back(5000.0);
				    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
			    }
		    }
	    }



	    //**************************************************************************
	    //finally, we encode the match point continuity constraints and their Jacobian entries,
	    //noting that every patch point constraint in the phase has a derivative with respect to every variable in the phase
	    //in addition, the patch point constraints have a derivative with respect to the previous phase's arrival mass
	    //and the patch point constraints have a derivative with respect to all previous time variables, including the launch date
	    calcbounds_LT_match_points(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	    //***************************************************************************
	    //Construct a helper array "rectangular prism" of G indices for match point derivatives with respect to state variables within the phase
	    vector < vector <int> > timeslice;
	    vector <int> scanline;
	    vector<string> constraint_type;
	    constraint_type.push_back(" patch point x constraint");
	    constraint_type.push_back(" patch point y constraint");
	    constraint_type.push_back(" patch point z constraint");
	    constraint_type.push_back(" patch point xdot constraint");
	    constraint_type.push_back(" patch point ydot constraint");
	    constraint_type.push_back(" patch point zdot constraint");
	    constraint_type.push_back(" patch point m constraint");

	    for (int step = 0; step < options->num_timesteps + 2; ++step)
	    {
		    timeslice.clear();
		
		    //loop over constraint type dimension
		    for (int c = 0; c < 7; ++c)
		    {
			    scanline.clear();

			    string constraintname = prefix + constraint_type[c];

			    for (size_t entry = 0; entry < Gdescriptions->size(); ++entry)
			    {
				    if ( (*Gdescriptions)[entry].find(constraintname) < 1024)
				    {
					    if (step == 0) //derivatives with respect to initial velocity
					    {
						    //the first step of the first phase of a journey is abnormal because instead of encoding an XYZ vector, we have encoded a magnitude and two angles
						    if (p == 0)
						    {
							    if ( (*Gdescriptions)[entry].find("magnitude of outgoing velocity asymptote") < 1024)
							    {
								    scanline.push_back(entry);     //derivative with respect to magnitude
								    scanline.push_back(entry + 1); //derivative with respect to RA
								    scanline.push_back(entry + 2); //derivative with respect to DEC
							    }
						    }
						    //otherwise look for an XYZ initial velocity increment
						    else
						    {
							    if ( (*Gdescriptions)[entry].find("initial velocity increment x") < 1024)
							    {
								    scanline.push_back(entry);     //derivative with respect to initial velocity increment x
								    scanline.push_back(entry + 1); //derivative with respect to initial velocity increment y
								    scanline.push_back(entry + 2); //derivative with respect to initial velocity increment z
							    }
						    }
					    }
					    //derivatives with respect to final velocity - these do not exist for a terminal rendezvous phase
					    else if (step == options->num_timesteps + 1 && (!(p == options->number_of_phases[j] - 1 && ((options->journey_arrival_type[j] == 1) || options->journey_arrival_type[j] == 3) || options->journey_arrival_type[j] == 5)))
					    {
						    if ( (*Gdescriptions)[entry].find("terminal velocity increment x") < 1024)
						    {
							    scanline.push_back(entry);	   //derivative with respect to terminal velocity increment x
							    scanline.push_back(entry + 1); //derivative with respect to terminal velocity increment y
							    scanline.push_back(entry + 2); //derivative with respect to terminal velocity increment z
						    }
					    }
					    //derivatives with respect to the control parameters
					    else
					    {
						    stringstream stepstream;
						    stepstream << step - 1;

						    string controlname = prefix + "step " + stepstream.str() + " u_x";

						    if ( (*Gdescriptions)[entry].find(controlname) < 1024)
						    {
							    scanline.push_back(entry);	   //derivative with respect to u_x
							    scanline.push_back(entry + 1); //derivative with respect to u_y
							    scanline.push_back(entry + 2); //derivative with respect to u_z
						    }
					    }

					
				    } //end if ( (*Gdescriptions)[entry].find(constraintname) )
			    } //end loop over Gdescriptions
			
			    timeslice.push_back(scanline);
			
		    } //end loop over constraint entries
		    match_point_constraint_G_indices.push_back(timeslice);
	
	    } //end loop over time steps

	    //***************************************************************************
	    //if this is the last phase, encode any constraints for the arrival processing
	    if (p == options->number_of_phases[j] - 1)
	    {
		    if (options->journey_arrival_type[j] == 2) //intercept with bounded v-infinity
		    {
			    Flowerbounds->push_back((options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1])*(options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1]) - 1);
			    Fupperbounds->push_back(0.0);
			    Fdescriptions->push_back(prefix + "arrival C3 constraint");

			    //Jacobian entry for a bounded v-infinity intercept
			    //this is a nonlinear constraint dependent only on the terminal velocity vector for this phase
			    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
			    {
				    if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
				    {
					    iGfun->push_back(Fdescriptions->size() - 1);
					    jGvar->push_back(entry);
					    stringstream EntryNameStream;
					    EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					    Gdescriptions->push_back(EntryNameStream.str());
					    terminal_velocity_constraint_G_indices.push_back(iGfun->size() - 1);
					    terminal_velocity_constraint_X_indices.push_back(entry);
					    terminal_velocity_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
				    }
			    }
		    }

		    if ( options->journey_arrival_declination_constraint_flag[j] && (options->journey_arrival_type[j] == 0 || options->journey_arrival_type[j] == 2) ) //intercept with bounded v-infinity or orbit insertion
		    {
			    Flowerbounds->push_back(options->journey_arrival_declination_bounds[j][0] / options->journey_arrival_declination_bounds[j][1] - 1.0);
			    Fupperbounds->push_back(0.0);
			    Fdescriptions->push_back(prefix + "arrival declination constraint");

			    //Jacobian entry for arrival declination constraint
			    //this is a nonlinear constraint dependent only on the terminal velocity vector for this phase
			    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
			    {
				    if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
				    {
					    iGfun->push_back(Fdescriptions->size() - 1);
					    jGvar->push_back(entry);
					    stringstream EntryNameStream;
					    EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					    Gdescriptions->push_back(EntryNameStream.str());
					    arrival_declination_constraint_G_indices.push_back(iGfun->size() - 1);
					    arrival_declination_constraint_X_indices.push_back(entry);
					    arrival_declination_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
				    }
			    }
		    }

            else if (options->journey_arrival_type[j] == 6) //enforce escape constraint, E = 0
		    {
			    Flowerbounds->push_back(-math::SMALL);
			    Fupperbounds->push_back(math::LARGE);
			    Fdescriptions->push_back(prefix + "arrival escape condition (zero energy constraint)");

			    //Jacobian entry for zero energy condition
			    //this is a nonlinear constraint dependent on the terminal velocity vector
			    //note that it is NOT dependent on position because these phases always end at the SOI and SOI distance is constant
			    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
			    {
				    if ((*Xdescriptions)[entry].find("terminal velocity") < 1024 || (*Xdescriptions)[entry].find("right boundary point") < 1024)
				    {
					    iGfun->push_back(Fdescriptions->size() - 1);
					    jGvar->push_back(entry);
					    stringstream EntryNameStream;
					    EntryNameStream << "Derivative of " << prefix << " escape constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					    Gdescriptions->push_back(EntryNameStream.str());
					    escape_constraint_G_indices.push_back(iGfun->size() - 1);
					    escape_constraint_X_indices.push_back(entry);
					    escape_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
				    }
			    }
		    }	
	    }
    }

    //output function
    //return 0 if successful, 1 if failure
    int FBLT_phase::output(missionoptions* options, const double& launchdate, int j, int p, EMTG::Astrodynamics::universe* Universe, int* eventcount)
    {
	    //Step 1: store data that will be used for the printing
	    double empty_vector[] = {0,0,0};
	    double phase_time_elapsed = 0.0;
	    string event_type;
	    double temp_power, temp_thrust, temp_mdot, temp_Isp,
		    temp_active_power, temp_dTdP, temp_dmdotdP,
		    temp_dTdIsp, temp_dmdotdIsp, temp_dPdr, temp_dPdt;
	    int temp_active_thrusters;
        double spacecraft_state_propagate[11 + 11 * 11];
        double resumeH = initial_coast_duration * 86400 / Universe->TU;
        double resumeError = 1.0e-8;
        double dummy_parameter = 0.0;

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
			    initial_Isp =  options->IspDS;
		    }
		    else
		    {
			    initial_Isp =  -1;
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

	    //print the departure/flyby
	    for (int k = 0; k < 3; ++k)
		    dVdeparture[k] = V_infinity_out(k);

	    //we have to calculate the available power at departure
	    Astrodynamics::find_engine_parameters(	options,
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
						    (p > 0 ? flyby_altitude : Bradius),
						    (Btheta),
						    (p > 0 ? flyby_turn_angle : -1),
						    RA_departure,
						    DEC_departure,
						    C3_departure,
						    state_at_beginning_of_phase,
						    dVdeparture,
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
	    if (this->detect_initial_coast)
	    {
            //if this is a launch AND we are doing a forced post-launch initial coast
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

            double augmented_state_at_initial_coast_midpoint[11 + 11 * 11];

            //initialize the propagator state
            for (size_t i = 7; i < 11; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_propagate[i] = 1.0;

            integrator->adaptive_step_int(spacecraft_state_propagate,
                                        augmented_state_at_initial_coast_midpoint,
                                        empty_vector,
                                        (phase_start_epoch) / Universe->TU,
                                        launchdate,
                                        initial_coast_duration / 2.0 / Universe->TU,
                                        &resumeH,
                                        &resumeError,
                                        1.0e-8,
                                        EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        (int*)&dummy_parameter,
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

            //propagate forward to the end of the initial coast
            integrator->adaptive_step_int(augmented_state_at_initial_coast_midpoint,
                                            spacecraft_state_propagate,
                                            empty_vector,
                                            (phase_start_epoch) / Universe->TU,
                                            launchdate,
                                            initial_coast_duration / 2.0 / Universe->TU,
                                            &resumeH,
                                            &resumeError,
                                            1.0e-8,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            (int*)&dummy_parameter,
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
		    Astrodynamics::find_engine_parameters(	options,
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

		    write_summary_line(	options,
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
            //initialize the propagator state
            for (size_t i = 7; i < 11; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;
            for (size_t i = this->STMrows; i < this->num_states; i = i + STMrows + 1)
                spacecraft_state_propagate[i] = 1.0;

            //first propagate to get the mid-point state of each step
            integrator->adaptive_step_int(spacecraft_state_propagate,
                                        this->spacecraft_state[step].data(),
                                        this->control[step].data(),
                                        (phase_start_epoch + phase_time_elapsed) / Universe->TU,
                                        launchdate,
                                        this->time_step_sizes[step] / 2.0 / Universe->TU,
                                        &resumeH,
                                        &resumeError,
                                        1.0e-8,
                                        EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        (int*)&dummy_parameter,
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

            //then propagate to the end of the step
            integrator->adaptive_step_int(this->spacecraft_state[step].data(),
                                        spacecraft_state_propagate,
                                        this->control[step].data(),
                                        (phase_start_epoch + phase_time_elapsed + 0.5 * this->time_step_sizes[step]) / Universe->TU,
                                        launchdate,
                                        this->time_step_sizes[step] / 2.0 / Universe->TU,
                                        &resumeH,
                                        &resumeError,
                                        1.0e-8,
                                        EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        &dummy_parameter,
                                        (int*)&dummy_parameter,
                                        this->STMrows,
                                        this->STMcolumns,
                                        (void*)options,
                                        (void*)Universe,
                                        DummyControllerPointer);

		    double angle1, angle2;
		    if (step >= (options->num_timesteps / 2))
		    {
			    angle1 =  atan2(control[step][1] + math::SMALL, control[step][0]) * EMTG::math::PI / 180.0;
			    angle2 =  asin(control[step][2] / EMTG::math::norm(control[step].data(), 3)) * EMTG::math::PI / 180.0;
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
			    scaled_state[k+3] = spacecraft_state[step][k+3] * Universe->LU / Universe->TU;
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
						    math::norm(this->control[step].data(),3) * this->available_mass_flow_rate[step],
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
						        empty_vector,
						        empty_vector,
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
            double augmented_state_at_terminal_coast_midpoint[11 + 11 * 11];
            
            double empty_vector[] = { 0.0, 0.0, 0.0 };
            double dummy_parameter = 0.0;

            //initial coast after flyby
            terminal_coast_duration = options->forced_flyby_coast;

            double resumeH = -terminal_coast_duration / 2.0 / Universe->TU;
            double resumeError = 1.0e-8;


            //The terminal coast STM entries of the state vector must be initialized to the identity
            for (size_t i = this->STMrows; i < this->num_states; ++i)
                spacecraft_state_propagate[i] = 0.0;

            integrator->adaptive_step_int(spacecraft_state_propagate,
                                            augmented_state_at_terminal_coast_midpoint,
                                            empty_vector,
                                            (phase_end_epoch + phase_time_elapsed) / Universe->TU,
                                            launchdate,
                                            terminal_coast_duration / 2.0 / Universe->TU,
                                            &resumeH,
                                            &resumeError,
                                            1.0e-8,
                                            EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            &dummy_parameter,
                                            (int*)&dummy_parameter,
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
		    Astrodynamics::find_engine_parameters(	options,
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


		    write_summary_line(	options,
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
		    Astrodynamics::find_engine_parameters(	options,
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
								    this->spiral_capture_time / 86400.0,
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
        double current_state[7];
        double current_epoch;

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
        else if ((p > 0 || (p == 0 && (options.journey_departure_elements_type[j] == 3 || options.journey_departure_elements_type[j] == 4 || options.journey_departure_elements_type[j] == 6)))
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
                                                double* current_state,
                                                const int& control_step,
                                                vector< vector<string> >& output_line_array,
                                                const int& j,
                                                const int& p,
                                                const double& propagation_time,
                                                const double& journey_starting_epoch)
    {
        vector<string> output_line(7);
        double time_remaining = propagation_time;
        double empty_control[] = {0.0, 0.0, 0.0};
        double temp_state[7];

        double available_thrust, available_mass_flow_rate, available_Isp, available_power, active_power;
        int number_of_active_engines;

        for (size_t k = 0; k < 3; ++k)
        {
            current_state[k] /= Universe.LU;
            current_state[k+3] *= Universe.TU / Universe.LU;
        }
        current_state[6] /= options.maximum_mass;

        while (time_remaining > 0.0)
        {
            double step_time = time_remaining > 86400.0 ? 86400.0 : time_remaining;
            double resumeError = 1.0e-8;
            double resumeH = step_time / Universe.TU;

            this->integrator->adaptive_step_int(current_state,
                                                temp_state,
                                                control_step >= 0 ? this->control[control_step].data() : empty_control,
                                                current_epoch / Universe.TU,
                                                (current_epoch - launch_epoch) / Universe.TU,
                                                step_time / Universe.TU,
                                                &resumeH,
                                                &resumeError,
                                                1.0e-8,
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
            for (size_t k = 0; k < 7; ++k)
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
		std::vector < EMTG::math::Matrix< double > > & STM_archive_forward,
		std::vector < EMTG::math::Matrix< double > > & STM_archive_backward,
		missionoptions* options,
		EMTG::Astrodynamics::universe* Universe)
	{
        //build the derivatives matrix through successive STM multiplication
        //every time step we add one more STM to the chain
        //this is symmetric for forward and backward propagation
        //and before we can do that, we need to create a vector of STMs with the control entries zeroed out
        vector<EMTG::math::Matrix <double>> forward_stripped_STMs = this->STM_archive_forward;
        vector<EMTG::math::Matrix <double>> backward_stripped_STMs = this->STM_archive_backward;
        EMTG::math::Matrix<double> forward_cumulative_stripped_STM(this->STMrows, this->STMcolumns, 0.0);
        EMTG::math::Matrix<double> backward_cumulative_stripped_STM(this->STMrows, this->STMcolumns, 0.0);

        //initialize
        for (size_t i = 0; i < this->STMrows; ++i)
        {
            forward_cumulative_stripped_STM(i, i) = 1.0;
            backward_cumulative_stripped_STM(i, i) = 1.0;
        }

        //step through and create STM chains
        for (int step = options->num_timesteps / 2 - 1; step >=0; --step)
        {
            //create a stripped version of this step's STM
            for (size_t row = 0; row < 7; ++row)
            {
                for (size_t column = 7; column < 10; ++column)
                {
                    forward_stripped_STMs[step](row, column) = 0.0;
                    backward_stripped_STMs[step](row, column) = 0.0;
                }
            }

            //compute the cumulative STM for this step
            for (int stepnext = options->num_timesteps / 2 - 1; stepnext > step; --stepnext)
            {
                forward_cumulative_stripped_STM *= forward_stripped_STMs[stepnext];
                backward_cumulative_stripped_STM *= backward_stripped_STMs[stepnext];
            }

            //multiply in this step's STM
            this->forward_cumulative_STM_archive[step] = forward_cumulative_stripped_STM * STM_archive_forward[step];
            this->backward_cumulative_STM_archive[step] = backward_cumulative_stripped_STM * STM_archive_backward[step];
        }
        //advanced the cumulative stripped STMs all the way to the endpoints
        forward_cumulative_stripped_STM *= forward_stripped_STMs[0];
        backward_cumulative_stripped_STM *= backward_stripped_STMs[0];

		//Match point constraint derivatives with respect to forward controls
		//
        for (int step = 0; step < options->num_timesteps / 2; ++step)
		{
			//place the derivatives in the Jacobian
			//since the match point constraint is defined backward - forward, 
			//the forward derivatives will have a negative sign

			//match point with respect to control of the current step
			for (size_t cindex = 0; cindex < 3; ++cindex)
			{
				for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    G[match_point_constraint_G_indices[step + 1][stateindex][cindex]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[step + 1][stateindex][cindex]]] * this->forward_cumulative_STM_archive[step](stateindex, 7 + cindex);
			}
		}


			

		//derivative with respect to Isp for forward propagation with VSI thruster
		if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
		{
			 
		}

		if (options->objective_type == 13) //derivatives with respect to power for minimum power problems
		{
			
		}

		if (p == 0)
		{
			if (j == 0 && options->allow_initial_mass_to_vary && !(options->journey_departure_type[j] == 5))
			{
				
			}
			if (options->journey_variable_mass_increment[j] && !(options->journey_departure_type[j] == 5))
			{
				
			}
			if (options->journey_departure_type[j] == 0)// && options->LV_type > 0)
			{
				
			}

			if (j > 0) //for successive journeys the mass at the beginning of the phase affects the following patch point
			{
				//this must be disabled for phases that start with spirals
				if (!(options->journey_arrival_type[j - 1] == 7))
				{
					
				}
			}
		}
		else//for other phases other than the first
		{
			
		}

		//derivatives of the constraints with respect to the forward flight time variables
		//this in turn is dependent on the velocity and acceleration of the left-most body
		if (options->derivative_type > 2)
		{
			

			//the following derivatives are for the current phase flight time ONLY
			if (options->derivative_type > 3)
			{
				
			}
		}
		
		//compute and store the backward derivatives of the match point constraints with respect to the control unit vector
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
			
		}

		//derivatives of the constraints with respect to the backward flight time variables
		//this in turn is dependent on the velocity and acceleration of the right-most body
		if (options->derivative_type > 2)
		{
			

			//the following derivatives are for the current phase flight time
			if (options->derivative_type > 3)
			{
				
			}
		}


		//derivative with respect to arrival velocity
		//only evaluated for phases that are not terminal intercepts
		if (!(p == options->number_of_phases[j] - 1 && ((options->journey_arrival_type[j] == 1) || options->journey_arrival_type[j] == 3) || options->journey_arrival_type[j] == 5 || options->journey_arrival_type[j] == 7))
		{
            /*
            if (detect_terminal_coast)
            {
                //first we need to construct the cumulative STM for the terminal coast
                EMTG::math::Matrix <double> cumulative_terminal_coast_STM = this->terminal_coast_STM * backward_cumulative_stripped_STM;

                //then fill out the derivatives for the terminal coast
                for (size_t vindex = 0; vindex < 3; ++vindex)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]]] * cumulative_terminal_coast_STM(stateindex, vindex + 3);
                    }
                }
            }
            else
            { 
                for (size_t vindex = 0; vindex < 3; ++vindex)
                {
                    for (size_t stateindex = 0; stateindex < 7; ++stateindex)
                    {
                        G[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[options->num_timesteps + 1][stateindex][vindex]]] * backward_cumulative_stripped_STM(stateindex, vindex + 3);
                    }
                }
            }
            */
		}

	}

} /* namespace EMTG */
