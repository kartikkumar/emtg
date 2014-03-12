/*
 * FBLTSphase.cpp
 *
 *  Created on: March 22nd, 2013
 *      Author: Jacob
 */

#include "FBLTSphase.h"
#include "Astrodynamics.h"
#include "missionoptions.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "equations_of_motion.h"
#include "universe.h"

#include "SpiceUsr.h"

#include <sstream>
#include <fstream>

namespace EMTG {

FBLTS_phase::FBLTS_phase() {
//default constructor does nothing

}

FBLTS_phase::FBLTS_phase(int j, int p, missionoptions* options)
{
	//must resize all data vectors to the corrrect length
	vector<double> state_dummy(8);
	vector<double> dV_or_control_dummy(3);

	for (int step = 0; step < options->num_timesteps; ++step) {
		spacecraft_state.push_back(state_dummy);
		control.push_back(dV_or_control_dummy);
	}

	match_point_state.resize(8);


	event_epochs.resize(options->num_timesteps);
	time_step_width.resize(options->num_timesteps);
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
	integrator = new EMTG::integration::rk8713M(8);

	current_mass_increment = 0.0;
	journey_initial_mass_increment_scale_factor = 1.0;
}

FBLTS_phase::~FBLTS_phase() {
	//delete the integrator object
		delete integrator;
}

//evaluate function
//return 0 if successful, 1 if failure
int FBLTS_phase::evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options) 
{
	//declare some local variables
	int errcode = 0;

	//******************************************************************
	//Steps 1-4: Process the left boundary condition
	process_left_boundary_condition(X, Xindex, F, Findex, G, Gindex, needG, current_epoch, current_state, current_deltaV, boundary1_state, boundary2_state, j, p, Universe, options);
		
	//******************************************************************
	//Step 5: For FBLTS, we need to know the state of the spacecraft at the right hand side (end) of the phase in order to propagate backward
	process_right_boundary_condition(X, Xindex, F, Findex, G, Gindex, needG, current_epoch, current_state, current_deltaV, boundary1_state, boundary2_state, j, p, Universe, options);

	//******************************************************************
	//Step 6: thrust and propagate forward and back

	//Step 6.1: determine the length of an s-step
	double s_step = s_final / options->num_timesteps;

	//Step 6.2: propagate forward

	//store the initial prefered integration step size
	double resumeH = s_step / 2.0;
	double resumeError = 1.0e-6;

	//first initialize the forward integration
	//the following array holds the spacecraft state at "half steps," i.e. halfway through each integration segment
	double spacecraft_state_forward[8];
	for (int k = 0; k < 7; ++k)
		spacecraft_state_forward[k] = state_at_beginning_of_phase[k];

	spacecraft_state_forward[7] = *current_epoch;

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

	//scale the time
	spacecraft_state_forward[7] *= 86400 / Universe->TU;

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

		//step 6.2.3 record the time at the beginning of the step
		double time_at_beginning_of_step = spacecraft_state_forward[7];
					
		//step 6.2.4 propagate the spacecraft to the midpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state_forward,
										spacecraft_state[step].data(),
										control[step].data(), 
										step * s_step,
										X[0],
										s_step / 2.0, 
										&resumeH,
										&resumeError,
										1.0e-10,
										EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust_sundman,
										&available_thrust[step],
										&available_mass_flow_rate[step],
										&available_Isp[step],
										&available_power[step],
										&active_power[step],
										&number_of_active_engines[step],
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );

		
		//step 6.2.5 encode the epoch of the step midpoint
		event_epochs[step] = spacecraft_state[step][7] * Universe->TU / 86400;
		
		//step 6.2.6 propagate the spacecraft to the endpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state[step].data(),
										spacecraft_state_forward,
										control[step].data(),  
										(step + 0.5) * s_step,
										X[0],
										s_step / 2.0, 
										&resumeH,
										&resumeError,
										1.0e-10,
										EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust_sundman,
										&available_thrust[step],
										&available_mass_flow_rate[step],
										&available_Isp[step],
										&available_power[step],
										&active_power[step],
										&number_of_active_engines[step],
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );

		//step 6.2.7 compute the width of the time step
		time_step_width[step] = (spacecraft_state_forward[7] - time_at_beginning_of_step) * Universe->TU / 86400;

	}

	//Step 6.3: propagate backward

	//store the initial prefered integration step size
	resumeH = -s_step / 2.0;
	resumeError = 1.0e-6;
	
	//first initialize the backward integration
	//the following array holds the spacecraft state at "half steps," i.e. integration steps where a burn is not applied
	double spacecraft_state_backward[8];
	for (int k = 0; k < 7; ++k)
		spacecraft_state_backward[k] = state_at_end_of_phase[k];
	spacecraft_state_backward[7] = (*current_epoch + TOF);

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

	//scale the time
	spacecraft_state_backward[7] *= 86400 / Universe->TU;

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

		//step 6.3.3 record the time at the end of the step
		double time_at_end_of_step = spacecraft_state_backward[7];
		
		//step 6.3.4 propagate the spacecraft to the midpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state_backward,
										spacecraft_state[backstep].data(),
										control[backstep].data(),  
										step * s_step,
										X[0],
										-s_step / 2.0, 
										&resumeH,
										&resumeError,
										1.0e-10,
										EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust_sundman,
										&available_thrust[backstep],
										&available_mass_flow_rate[backstep],
										&available_Isp[backstep],
										&available_power[backstep],
										&active_power[backstep],
										&number_of_active_engines[backstep],
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );
		
		//step 6.3.5 encode the epoch of the step midpoint
		event_epochs[backstep] = spacecraft_state[backstep][7] * Universe->TU / 86400;
		
		//step 6.3.6 propagate the spacecraft to the endpoint of the phase using the control unit vector
		integrator->adaptive_step_int(	spacecraft_state[backstep].data(),
										spacecraft_state_backward,
										control[backstep].data(),  
										(step + 0.5) * s_step,
										X[0],
										-s_step / 2.0, 
										&resumeH,
										&resumeError,
										1.0e-10,
										EMTG::Astrodynamics::EOM::EOM_inertial_continuous_thrust_sundman,
										&available_thrust[backstep],
										&available_mass_flow_rate[backstep],
										&available_Isp[backstep],
										&available_power[backstep],
										&active_power[backstep],
										&number_of_active_engines[backstep],
										(void*)options,
										(void*)Universe,
										DummyControllerPointer                                                  );

		//step 6.3.7 compute the width of the time step
		time_step_width[backstep] = (time_at_end_of_step - spacecraft_state_backward[7]) * Universe->TU / 86400;
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
		match_point_state[k] = spacecraft_state_forward[k];
		match_point_state[k+3] = spacecraft_state_forward[k+3];
	}
	//mass
	F[*Findex+6] = (spacecraft_state_backward[6] - spacecraft_state_forward[6]);
	
	//time
	F[*Findex+7] = (spacecraft_state_backward[7] - spacecraft_state_forward[7]);
	(*Findex) += 8;

	match_point_state[6] = spacecraft_state_forward[6] * options->maximum_mass;
	match_point_state[7] = spacecraft_state_forward[7] * Universe->TU / 86400;


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
				dV_arrival_magnitude = process_arrival(state_at_end_of_phase+3, boundary2_state, current_state+3, Body2->mu, Body2->r_SOI, F, Findex, j, options, Universe);
			else //arriving at a boundary point in free space
				dV_arrival_magnitude = process_arrival(state_at_end_of_phase+3, boundary2_state, current_state+3, Universe->mu, Universe->r_SOI, F, Findex, j, options, Universe);

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
			if ( needG && options->journey_arrival_declination_constraint_flag[j] && (options->journey_arrival_type[j] == 0 || options->journey_arrival_type[j] == 2) ) //intercept with bounded v-infinity or orbit insertion
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
	*current_epoch += TOF;

	//******************************************************************
	//Step 9: update the current state
	for (int k = 0; k < 7; ++k)
		current_state[k] = state_at_end_of_phase[k];

	return 0;
}


//bounds calculation function
//return 0 if successful, 1 if failure
int FBLTS_phase::calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
{
	//this function calculates the upper and lower bounds for the decision and constraint vectors for FBLTS
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
	//next, we need to include the decision variables and constraints for each burn
	//now, for each timestep
	for (int w=0; w < options->num_timesteps; ++w)
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
			EntryNameStream << "Derivative of " << prefix + "step" << " throttle magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
			Gdescriptions->push_back(EntryNameStream.str());

			//store the position in the G vector
			step_G_indices.push_back(iGfun->size() - 1);
		}
		control_vector_G_indices.push_back(step_G_indices);
	}




	//**************************************************************************
	//finally, we encode the match point continuity constraints and their Jacobian entries,
	//noting that every patch point constraint in the phase has a derivative with respect to every variable in the phase
	//in addition, the patch point constraints have a derivative with respect to the previous phase's arrival mass
	//and the patch point constraints have a derivatives with respect to all previous time variables, including the launch date
	calcbounds_LT_match_points(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

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

        if (options->journey_arrival_type[j] == 6) //enforce escape constraint, E = 0
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
int FBLTS_phase::output(missionoptions* options, const double& launchdate, int j, int p, EMTG::Astrodynamics::universe* Universe, int* eventcount)
{
	//Step 1: store data that will be used for the printing
	double empty_vector[] = {0,0,0};
	string event_type;
	if (p == 0)
	{
		if (j == 0 && boundary1_location_code > 0 && options->LV_type >= 0)
			event_type = "launch";
		else
		{
			event_type = "departure";
			math::Matrix<double> periapse_state = calculate_flyby_periapse_state(V_infinity_in, V_infinity_out, flyby_altitude, *Body1);
			Bplane.define_bplane(V_infinity_in, BoundaryR, BoundaryV);
			Bplane.compute_Bradius_Btheta_from_periapse_position(Body1->mu, V_infinity_in, periapse_state, &Bradius, &Btheta);
		}
	}
	else
		event_type = "upwr_flyby";

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

	//Step 2: all phases have at least two events: departure/flyby and burn
	//*****************************************************************************
	//first let's print the departure/flyby
	
	for (int k = 0; k < 3; ++k)
		dVdeparture[k] = V_infinity_out(k);

	write_summary_line(options,
						Universe,
						eventcount,
						phase_start_epoch,
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
						(p == 0 ? dV_departure_magnitude : flyby_outgoing_v_infinity),
						-1,
						initial_Isp,
						-1,
						0,
						0,
						0);


	//*****************************************************************************
	//next, we must print each thrust arc

	for (int step = 0; step < options->num_timesteps; ++step)
	{
		if (EMTG::math::norm(control[step].data(), 3) > 1.0e-2)
		{
			event_type = "FBLTSthrust";
		}
		else
		{
			event_type = "coast";
		}

		double angle1, angle2;
		if (step >= (options->num_timesteps / 2))
		{
			angle1 =  atan2(control[step][1], control[step][0]) * EMTG::math::PI / 180.0;
			angle2 =  asin(control[step][2] / EMTG::math::norm(control[step].data(), 3)) * EMTG::math::PI / 180.0;
		}
		else
		{
			angle1 = atan2(-control[step][1], -control[step][0]) * EMTG::math::PI / 180.0;
			angle2 = asin(-control[step][2] / EMTG::math::norm(control[step].data(), 3)) * EMTG::math::PI / 180.0;
		}

		double current_Isp, current_thrust, current_power, current_mass_flow_rate;
		if (options->engine_type == 0) //fixed thrust/Isp
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
			EMTG::Astrodynamics::find_engine_parameters(options,
														EMTG::math::norm(spacecraft_state[step].data(), 3),
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
			//current_Isp = available_thrust[step] / mdot / options->g0;
			current_power = available_power[step];
			current_mass_flow_rate = current_thrust / current_Isp / options->g0;
		}
		double thrust_vector[3];
		for (int k = 0; k < 3; ++k)
			thrust_vector[k] = control[step][k] * current_thrust;

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
						event_epochs[step],
						event_type,
						"deep-space",
						time_step_width[step],
						-1,
						-1,
						-1,
						angle1,
						angle2,
						0,
						scaled_state,
						control[step].data(),
						thrust_vector,
						EMTG::math::norm(control[step].data(), 3),
						current_thrust,
						current_Isp,
						current_power,
						math::norm(control[step].data(),3) * current_mass_flow_rate,
						number_of_active_engines[step],
						active_power[step]);

		//if we have stepped halfway through, insert the match point line
		if (step == options->num_timesteps / 2 - 1)
			write_summary_line(options,
						Universe,
						eventcount,
						match_point_state[7],
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
		}
		else
		{
			dV_arrival_mag = dV_arrival_magnitude;
		}

		write_summary_line(options,
						Universe,
						eventcount,
						phase_start_epoch + TOF,
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
	}

	return 0;
}

} /* namespace EMTG */
