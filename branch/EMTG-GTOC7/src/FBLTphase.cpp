/*
 * FBLTphase.cpp
 *
 *  Created on: September 17, 2012
 *      Author: Jacob
 */

#include "FBLTphase.h"
#include "Astrodynamics.h"
#include "missionoptions.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "equations_of_motion.h"
#include "universe.h"

#include "SpiceUsr.h"

#include <sstream>
#include <fstream>

#include "EMTG_string_utilities.h"

namespace EMTG {

FBLT_phase::FBLT_phase() {
//default constructor does nothing

}

FBLT_phase::FBLT_phase(int j, int p, missionoptions* options)
{
	//must resize all data vectors to the correct length
	vector<double> state_dummy(7);
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

	//set the bodies
	boundary1_location_code = options->sequence[j][p];
	boundary2_location_code = options->sequence[j][p+1];

	//size the vectors that will be used to calculate the b-plane
	V_infinity_in.resize(3, 1);
	V_infinity_out.resize(3, 1);
	BoundaryR.resize(3, 1);
	BoundaryV.resize(3, 1);

	//set up the integrator
	integrator = new EMTG::integration::rk8713M(7);

	current_mass_increment = 0.0;
	journey_initial_mass_increment_scale_factor = 1.0;

	//size the time step vector
	time_step_sizes.resize(options->num_timesteps);

	//set derivatives for spirals
	this->spiral_escape_dm_after_dm_before = 1.0;
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
	else if ( (p > 0 || p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4) ) && options->forced_flyby_coast > 1.0e-6)
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

	//first initialize the forward integration
	//the following array holds the spacecraft state at "half steps," i.e. halfway through each integration segment
	double spacecraft_state_forward[7];
	for (int k = 0; k < 7; ++k)
		spacecraft_state_forward[k] = state_at_beginning_of_phase[k];

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
	if ( (j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6) || ( (p > 0 || p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4) ) && options->forced_flyby_coast > 1.0e-6) )
	{
		//if this is a launch AND we are doing a forced post-launch initial coast
		double spacecraft_state_end_coast[7];
		double empty_vector[] = {0.0,0.0,0.0};
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

		double resumeH = initial_coast_duration/2.0 * 86400/Universe->TU;

		integrator->adaptive_step_int(	spacecraft_state_forward,
										state_at_initial_coast_midpoint,
										empty_vector, 
										(phase_start_epoch) / Universe->TU,
										X[0],
										initial_coast_duration / 2.0 /Universe->TU, 
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
										(void*)options,
										(void*)Universe,
										DummyControllerPointer      );	

		integrator->adaptive_step_int(	state_at_initial_coast_midpoint,
										spacecraft_state_end_coast,
										empty_vector, 
										(phase_start_epoch + initial_coast_duration / 2.0) / Universe->TU,
										X[0],
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
										(int*) &dummy_parameter,
										(void*)options,
										(void*)Universe,
										DummyControllerPointer      );	


		phase_time_elapsed_forward += initial_coast_duration;

		for (int k = 0; k < 7; ++k)
			spacecraft_state_forward[k] = spacecraft_state_end_coast[k];
	}

	for (int step = 0; step < options->num_timesteps/2; ++step)
	{
		//step 6.2.1 extract the control unit vector from the decision vector
		control[step][0] = X[*Xindex];
		control[step][1] = X[*Xindex+1];
		control[step][2] = X[*Xindex+2];
		

		//extract the specific impulse for this step (VSI only)
		if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
		{
			available_Isp[step] = X[*Xindex+3];
			++(*Xindex);
		}

		double throttle = math::norm(control[step].data(), 3);
		if (options->derivative_type > 0 && needG)
		{
			G[control_vector_G_indices[step][0]] = 2.0 * control[step][2] / throttle;
			G[control_vector_G_indices[step][1]] = 2.0 * control[step][1] / throttle;
			G[control_vector_G_indices[step][2]] = 2.0 * control[step][0] / throttle;
			(*Gindex) += 3;
		}
		(*Xindex) +=3;

		//step 6.2.2 apply the control unit vector magnitude constraint
		F[*Findex] = throttle;
		++(*Findex);
					
		//step 6.2.3 propagate the spacecraft to the midpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state_forward,
										spacecraft_state[step].data(),
										control[step].data(), 
										(phase_start_epoch + phase_time_elapsed_forward) / Universe->TU,
										X[0],
										time_step_sizes[step] / 2.0 / Universe->TU, 
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
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );

		
		//step 6.2.4 encode the epoch of the step midpoint
		event_epochs[step] = phase_start_epoch + phase_time_elapsed_forward + 0.5 * time_step_sizes[step];
		phase_time_elapsed_forward += time_step_sizes[step];
		
		//step 6.2.5 propagate the spacecraft to the endpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state[step].data(),
										spacecraft_state_forward,
										control[step].data(),  
										(event_epochs[step]) / options->TU,
										X[0],
										time_step_sizes[step]/2.0 / Universe->TU, 
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
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );

	}

	//Step 6.3: propagate backward
	phase_time_elapsed_backward = 0.0;
	//store the initial prefered integration step size
	resumeH = -time_step_sizes[options->num_timesteps - 1] / Universe->TU;
	resumeError = 1.0e-6;
	
	//first initialize the backward integration
	//the following array holds the spacecraft state at "half steps," i.e. integration steps where a burn is not applied
	double spacecraft_state_backward[7];
	for (int k = 0; k < 7; ++k)
		spacecraft_state_backward[k] = state_at_end_of_phase[k];

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
	if ( (p < options->number_of_phases[j] - 1 ||  (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5) ) && options->forced_flyby_coast > 1.0e-6)
	{
		double spacecraft_state_end_coast[7];
		double empty_vector[] = {0.0,0.0,0.0};
		double dummy_parameter = 0.0;

		//initial coast after flyby
		terminal_coast_duration = options->forced_flyby_coast;

		double resumeH = -terminal_coast_duration / 2.0 / Universe->TU;

		integrator->adaptive_step_int(	spacecraft_state_backward,
										state_at_terminal_coast_midpoint,
										empty_vector, 
										(phase_end_epoch) / Universe->TU,
										X[0],
										-terminal_coast_duration / 2.0 / Universe->TU, 
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
										(void*)options,
										(void*)Universe,
										DummyControllerPointer      );	

		integrator->adaptive_step_int(	state_at_terminal_coast_midpoint,
										spacecraft_state_end_coast,
										empty_vector, 
										(phase_end_epoch - terminal_coast_duration / 2.0) / Universe->TU,
										X[0],
										-terminal_coast_duration / 2.0 / Universe->TU, 
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
										(void*)options,
										(void*)Universe,
										DummyControllerPointer      );	


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

		double throttle = math::norm(control[backstep].data(), 3);
		if (options->derivative_type > 0 && needG)
		{
			G[control_vector_G_indices[backstep][0]] = 2.0 * control[backstep][2] / throttle;
			G[control_vector_G_indices[backstep][1]] = 2.0 * control[backstep][1] / throttle;
			G[control_vector_G_indices[backstep][2]] = 2.0 * control[backstep][0] / throttle;
			(*Gindex) += 3;
		}

		//step 6.3.2 apply the control unit vector magnitude constraint
		F[*Findex + (backstep - options->num_timesteps/2)] = throttle;
		
		//step 6.3.3 propagate the spacecraft to the midpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state_backward,
										spacecraft_state[backstep].data(),
										control[backstep].data(),  
										(phase_end_epoch - phase_time_elapsed_backward) / Universe->TU,
										X[0],
										-time_step_sizes[backstep] / 2.0 / Universe->TU, 
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
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );
		
		//step 6.3.4 encode the epoch of the step midpoint
		event_epochs[backstep] = phase_end_epoch - phase_time_elapsed_backward - 0.5 * time_step_sizes[backstep];
		phase_time_elapsed_backward += time_step_sizes[backstep];
		
		//step 6.3.5 propagate the spacecraft to the endpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state[backstep].data(),
										spacecraft_state_backward,
										control[backstep].data(),  
										(event_epochs[backstep]) / Universe->TU,
										X[0],
										-time_step_sizes[backstep] / 2.0 / Universe->TU, 
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
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );
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
		F[*Findex+k] = (spacecraft_state_backward[k] - spacecraft_state_forward[k]);
		
		//velocity
		F[*Findex+k+3] = (spacecraft_state_backward[k+3] - spacecraft_state_forward[k+3]);

		//match point state
		match_point_state[k] = spacecraft_state_forward[k] * Universe->LU;
		match_point_state[k+3] = spacecraft_state_forward[k+3] * Universe->LU / Universe->TU;
	}
	//mass
	F[*Findex+6] = (spacecraft_state_backward[6] - spacecraft_state_forward[6]);
	(*Findex) += 7;

	match_point_state[6] = spacecraft_state_forward[6] * options->maximum_mass;


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
int FBLT_phase::calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
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
	//and the patch point constraints have a derivatives with respect to all previous time variables, including the launch date
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
		
	return 0;
}

//output function
//return 0 if successful, 1 if failure
int FBLT_phase::output(missionoptions* options, const double& launchdate, int j, int p, EMTG::Astrodynamics::universe* Universe, int* eventcount)
{
	//Step 1: store data that will be used for the printing
	double empty_vector[] = {0,0,0};
	double phase_time_elapsed = 0.0;
	string event_type;
	if (p == 0)
	{
		if (j == 0 && boundary1_location_code > 0 && options->LV_type >= 0)
			event_type = "launch";
		else
			event_type = "departure";
	}
	else
	{
		event_type = "upwr_flyby";
		math::Matrix<double> periapse_state = calculate_flyby_periapse_state(V_infinity_in, V_infinity_out, flyby_altitude, *Body1);
		math::Matrix<double> periapse_R(3,1);
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
						-1,
						0,
						0,
						0);

	//*****************************************************************************
	//next, if there was an initial coast, we must print it
	//we'll have to start by propagating to the halfway point of the initial coast step
	if ( (j == 0 && p == 0 && options->forced_post_launch_coast > 1.0e-6) || ( (p > 0 || p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4) ) && options->forced_flyby_coast > 1.0e-6) )
	{
		for (int k = 0; k < 3; ++k)
		{
			state_at_initial_coast_midpoint[k] *= Universe->LU;
			state_at_initial_coast_midpoint[k+3] *= Universe->LU / Universe->TU;
		}
		state_at_initial_coast_midpoint[6] *= options->maximum_mass;

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
							0.0,
							0.0,
							0,
							0.0);

		phase_time_elapsed += initial_coast_duration;
	}


	//*****************************************************************************
	//next, we must print each thrust arc

	for (int step = 0; step < options->num_timesteps; ++step)
	{
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

		/*if (options->engine_type == 0) //fixed thrust/Isp
		{
			current_thrust = options->Thrust;
			current_Isp = options->IspLT;
			current_power = -1;
			current_mass_flow_rate = current_thrust / current_Isp / options->g0;
		}
		else //thrust, Isp are functions of power
		{
			double mdot;
			double dTdP, dmdotdP, dTdIsp, dmdotdIsp, dPdr, dPdt;
			EMTG::Astrodynamics::find_engine_parameters(options, EMTG::math::norm(spacecraft_state[step].data(), 3),
														event_epochs[step] - launchdate,
														&available_thrust[step],
														&mdot,
														&current_Isp,
														&available_power[step],
														&active_power[step],
														&number_of_active_engines[step],
														false,
														&dTdP,
														&dmdotdP,
														&dTdIsp,
														&dmdotdIsp,
														&dPdr,
														&dPdt);

			current_thrust = available_thrust[step];
			current_power = available_power[step];
			current_mass_flow_rate = current_thrust / current_Isp / options->g0;
		}*/
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
						math::norm(control[step].data(),3) * available_mass_flow_rate[step],
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
	if ( (p < options->number_of_phases[j] - 1 ||  (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 5) ) && options->forced_flyby_coast > 1.0e-6)
	{
		for (int k = 0; k < 3; ++k)
		{
			state_at_terminal_coast_midpoint[k] *= Universe->LU;
			state_at_terminal_coast_midpoint[k+3] *= Universe->LU / Universe->TU;
		}
		state_at_terminal_coast_midpoint[6] *= options->maximum_mass;

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
							0.0,
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
						-1,
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

//function to output in GTOC7 format
void FBLT_phase::output_GTOC7_format(missionoptions* options, EMTG::Astrodynamics::universe* Universe, const std::string& GTOC_output_file, int j, int p)
{
	//this is a daily format. We will not interpolate, rather we will integrate with piecewise constant control
	//in this format we assume that each journey has only one phase, which is reasonable for GTOC7 where all journeys end in rendezvous and no flybys are allowed

	ofstream GTOC7file;

	if (j == 0)
	{
		GTOC7file.open(GTOC_output_file.c_str(), ios::trunc);
		GTOC7file << "# PROBE NUMBER:         X" << endl;
	}
	else
		GTOC7file.open(GTOC_output_file.c_str(), ios::app);

	

	if (j == 0)
	{
		GTOC7file << "# PHASE NUMBER:         " << j + 1 << endl;
		GTOC7file << "# DESCRIPTION: From Mother ship to Asteroid " << this->Body1->spice_ID << endl;
		GTOC7file << "#  Time (MJD)             x (km)                 y (km)                 z (km)                 vx (km/s)              vy (km/s)              vz (km/s)              mass (kg)              Thrust_x (N)           Thrust_y (N)           Thrust_z (N)" << endl;

		//find the position and velocity of the first body at the launch date minus thirty days
		double probe_release_state[7];
		this->Body1->locate_body(this->phase_start_epoch - 30.0 * 86400.0, probe_release_state, false, options);
		probe_release_state[6] = 2000.0;

		//write out the departure state
		GTOC7file << " ";
		GTOC7file.precision(14);
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->phase_start_epoch / 86400.0 - 30.0, 2) << " ";
		for (int k = 0; k < 7; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(probe_release_state[k], 2) << " ";
		for (int k = 0; k < 3; ++k)
		{
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(0.0, 2);
			if (k < 2)
				GTOC7file << " ";
		}

		GTOC7file << endl;
	}
		
	GTOC7file << "# PHASE NUMBER:         " << j + 2 << endl;
	GTOC7file << "# DESCRIPTION: From Asteroid " << this->Body1->spice_ID<< " to Asteroid " << this->Body2->spice_ID << endl;

	GTOC7file << "#  Time (MJD)             x (km)                 y (km)                 z (km)                 vx (km/s)              vy (km/s)              vz (km/s)              mass (kg)              Thrust_x (N)           Thrust_y (N)           Thrust_z (N)" << endl;

	//write out the departure state
	GTOC7file << " ";
	GTOC7file.precision(14);
	GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->phase_start_epoch / 86400.0, 2) << " ";
	for (int k = 0; k < 7; ++k)
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_beginning_of_phase[k], 2) << " ";
	for (int k = 0; k < 3; ++k)
	{
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->available_thrust[0] * 1000.0 * this->control[0][k], 2);
		if (k < 2)
			GTOC7file << " ";
	}
	
	GTOC7file << endl;	

	//now loop over days
	//assume all time-steps are the same width, i.e. no funky distributions
	double initial_state[7], propagated_state[7];
	for (int k = 0; k < 3; ++k)
		propagated_state[k] = this->state_at_beginning_of_phase[k] / Universe->LU;
	for (int k = 3; k < 6; ++k)
		propagated_state[k] = this->state_at_beginning_of_phase[k] / Universe->LU * Universe->TU;
	propagated_state[6] = 1.0;

	double time_step_width = this->time_step_sizes[0];
	for (int day = 1; day <= (int) floor(this->TOF / 86400.0); ++day)
	{
		for (int k = 0; k < 7; ++k)
			initial_state[k] = propagated_state[k];

		int control_index = day / ((int) ceil(time_step_width / 86400.0) );
		double day_control[3];
		day_control[0] = this->control[control_index][0];
		day_control[1] = this->control[control_index][1];
		day_control[2] = this->control[control_index][2];

		double resumeH = 86400.0 / Universe->TU;
		double resumeError = 1.0e-13;
		double available_thrust, available_mass_flow_rate, available_Isp, available_power, active_power;
		int number_of_active_engines;

		this->integrator->adaptive_step_int(	initial_state,
												propagated_state,
												day_control,
												0.0,
												86400.0 / Universe->TU,
												86400.0 / Universe->TU,
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
												(void*)options,
												(void*)Universe,
												this->DummyControllerPointer);

		GTOC7file << " ";
		GTOC7file.precision(14);
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->phase_start_epoch / 86400.0 + day, 2) << " ";
		for (int k = 0; k < 3; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(propagated_state[k] * Universe->LU, 2) << " ";
		for (int k = 3; k < 6; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(propagated_state[k] * Universe->LU / Universe->TU, 2) << " ";
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(propagated_state[6] * this->state_at_beginning_of_phase[6], 2) << " ";
		for (int k = 0; k < 3; ++k)
		{
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->available_thrust[0] * 1000.0 * day_control[k], 2);
			if (k < 2)
				GTOC7file << " ";
		}

		GTOC7file << endl;
	}

	//write out the arrival state
	GTOC7file << " ";
	GTOC7file.precision(14);
	GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->phase_end_epoch / 86400.0, 2) << " ";
	for (int k = 0; k < 7; ++k)
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_end_of_phase[k], 2) << " ";
	for (int k = 0; k < 3; ++k)
	{
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(0.0, 2);// this->available_thrust.back() * 1000.0 * this->control.back()[k], 2);
		if (k < 2)
			GTOC7file << " ";
	}

	GTOC7file << endl;
}

} /* namespace EMTG */
