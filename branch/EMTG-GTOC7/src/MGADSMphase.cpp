/*
 * MGADSMphase.cpp
 *
 *  Created on: Jul 15, 2012
 *      Author: Jacob
 */

#include "MGADSMphase.h"
#include "missionoptions.h"
#include "Astrodynamics.h"

#ifdef _EMTG_Lambert
#include "Lambert.h"
#endif

#include "Kepler_Lagrange_Laguerre_Conway_Der.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "universe.h"

#include "SpiceUsr.h"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "EMTG_string_utilities.h"

using namespace std;


namespace EMTG {

MGA_DSM_phase::MGA_DSM_phase() :
	time_before_burn(0),
	time_after_burn(0),
	b_plane_insertion_angle(0)
{
	//default constructor doesn't do anything
}

MGA_DSM_phase::MGA_DSM_phase(int j, int p, missionoptions* options) :
	time_before_burn(0),
	time_after_burn(0),
	b_plane_insertion_angle(0)
{
	//every phase has at least one burn
	int size = 1;
	//if this phase begins with a departure, then there is an additional burn
	if (p == 0 && !(options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4 || options->journey_departure_type[j] == 6))
		++size;

	//if this phase ends with an arrival, then there is an additional burn
	if (p == options->number_of_phases[j] - 1)
		++size;

	dVmag.resize(size);

	//set the bodies
	boundary1_location_code = options->sequence[j][p];
	boundary2_location_code = options->sequence[j][p+1];

	current_mass_increment = 0.0;
	journey_initial_mass_increment_scale_factor = 1.0;

	V_infinity_in.resize(3,1);
	V_infinity_out.resize(3,1);
	BoundaryR.resize(3,1);
	BoundaryV.resize(3,1);
}

MGA_DSM_phase::~MGA_DSM_phase() 
{
	//destructor does not need to do anything special
}

//evaluate function
//return 0 if successful, 1 if failure
int MGA_DSM_phase::evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
{
	//declare some local variables
	int errcode = 0;
	double lambert_v1[3];
	double lambert_v2[3];
	double vinf;
	double xaxis[] = {1,0,0};
	double hz1, hz2; //z component of angular momentum
	int dVindex = 0;

	//set pointers to the bodies
	if (boundary1_location_code > 0)
		Body1 = &Universe->bodies[boundary1_location_code - 1];
	if (boundary2_location_code > 0)
		Body2 = &Universe->bodies[boundary2_location_code - 1];

	//we need to know if we are the first phase of the journey and the journey does not start with a flyby
	if (p == 0 && !(options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 6))
	{
		//Step 1: extract the journey start epoch
		//if this is the first journey, the first decision variable is the starting epoch
		//otherwise, the first decision variable is the journey wait time
		if (j == 0)
		{
			*current_epoch = X[*Xindex];
			++(*Xindex);
		}
		else
		{
			//if we are not the first journey, are we starting from a hyperbolic arrival? or the boundary of the sphere of influence?
			//if so, then there is no wait time. If not, then the first decision variable is the stay time at the first body in the journey (i.e. at the asteroid for sample return)
			if (!(options->sequence[j-1][p+1] == -1 || boundary1_location_code == -1))
			{
				*current_epoch += X[*Xindex];
				++(*Xindex);
			}
		}

		//Step 2: locate the first boundary point
		locate_boundary_point(boundary1_location_code, options->journey_departure_type[j], true, Universe, boundary1_state, current_state+3, *current_epoch, X, Xindex, F, Findex, G, Gindex, false, j, p, options);
		

		//Step 3: compute the departure asymptote
		//Step 3.1 extract the departure parameters
		if (j == 0 || !(options->journey_departure_type[j] == 3))
		{
			if (boundary1_location_code == -2) //if we are starting at periapse of a hyperbola
			{
				//find the velocity periapse of the incoming hyperbola
				//note this is a scalar quantity - we already chose, in phase::locate_boundary_point(), b-plane parameters that guarantee that the incoming hyperbola is coplanar with the capture orbit
				double rp_inbound_hyperbola = math::norm(boundary1_state, 3);

				//current_state+3 holds the V-infinity after phase::locate_boundary_point() is called
				double v_infinity_inbound_hyperbola = math::norm(current_state+3, 3);
				double vp_inbound_hyperbola = sqrt(2 * Universe->mu / rp_inbound_hyperbola - v_infinity_inbound_hyperbola * v_infinity_inbound_hyperbola);

				//compute vinf - in this context it is not actually v-infinity, rather it is the magnitude of the velocity difference between the hyperbola and the capture ellipse
				double vp_capture_orbit = math::norm(boundary1_state+3, 3);
				double dv_capture = vp_capture_orbit - vp_inbound_hyperbola;
				vinf = fabs(dv_capture);

				//compute RA and DEC of the delta-V vector - this will be opposite from RA and DEC of the capture orbit state since the burn is tangential
				double tempvec[3];
				tempvec[0] = -boundary1_state[3];
				tempvec[1] = -boundary1_state[4];
				tempvec[2] = -boundary1_state[5];
				DEC_departure = asin(tempvec[2] / vp_capture_orbit);
				RA_departure = atan2(tempvec[1], tempvec[0]);

				//Step 3.2 compute the outgoing velocity vector
				dVdeparture[0] = vinf * cos(RA_departure)*cos(DEC_departure);
				dVdeparture[1] = vinf * sin(RA_departure)*cos(DEC_departure);
				dVdeparture[2] = vinf * sin(DEC_departure);
			}
			else
			{
				vinf = X[*Xindex];
				C3_departure = vinf*vinf;
				RA_departure = X[*Xindex + 1];
				DEC_departure = X[*Xindex + 2];
				*Xindex += 3;

				//Step 3.2 compute the outgoing velocity vector
				dVdeparture[0] = vinf * cos(RA_departure)*cos(DEC_departure);
				dVdeparture[1] = vinf * sin(RA_departure)*cos(DEC_departure);
				dVdeparture[2] = vinf * sin(DEC_departure);
			}
		}
		else //if this is a successive phase starting from an SOI, then we do NOT want to have a burn at the SOI
		{
			vinf = 0;
			C3_departure = 0;
			RA_departure = 0;
			DEC_departure = 0;

			//Step 3.2 compute the outgoing velocity vector
			dVdeparture[0] = 0;
			dVdeparture[1] = 0;
			dVdeparture[2] = 0;
		}
		
		//Step 4: compute the state post-departure
		//Step 4.1: compute position and velocity
		if (boundary1_location_code == -2)
		{
			for (int k = 0; k < 3; ++k)
			{
				state_at_beginning_of_phase[k] = boundary1_state[k];
				state_at_beginning_of_phase[k+3] = boundary1_state[k+3];
				boundary1_state[k+3] -= dVdeparture[k];
			}
		}
		else
		{
			for (int k = 0; k < 3; ++k)
			{
				state_at_beginning_of_phase[k] = boundary1_state[k];
				state_at_beginning_of_phase[k+3] = boundary1_state[k+3] + dVdeparture[k];
			}
		}
		//Step 4.2 compute the mass post-departure
		if (options->journey_departure_type[j] == 0 || options->journey_departure_type[j] == 2) //this is a direct insertion or launch
		{
			dVmag[dVindex] = vinf;
			if (options->include_initial_impulse_in_cost)
				*current_deltaV += vinf;
			++dVindex;
			
			if (j == 0) //this is a launch
			{

				//figure out how much mass we can launch
				if (options->LV_type > 0 || options->LV_type == -2)
				{
					double launch_mass;
					EMTG::Astrodynamics::find_mass_to_orbit(C3_departure, DEC_departure * 180.0 / math::PI, options->LV_type, &launch_mass, &dmdvinf, options);
					state_at_beginning_of_phase[6] = launch_mass > options->maximum_mass ? options->maximum_mass : launch_mass;
					
					//apply the mass margin
					state_at_beginning_of_phase[6] *= 1.0 - options->LV_margin;
				}
				else if (options->LV_type == 0)
				{
					state_at_beginning_of_phase[6] = options->maximum_mass;

					//add the starting mass increment
					state_at_beginning_of_phase[6] += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					dmdvinf = 0.0;
				}
				else //chemical burn using departure stage engine
				{
					double expfun = exp(-vinf * 1000 / (options->IspDS * options->g0));

					double initialmass = options->maximum_mass;

					//add the starting mass increment
					initialmass += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					state_at_beginning_of_phase[6] = initialmass * expfun;
					dmdvinf = -initialmass * 1000 / (options->IspDS * options->g0) * expfun;
				}
			}
			else if (options->journey_departure_type[j] == 0)
			{
				double expfun = exp(-vinf * 1000 / (options->IspChem * options->g0));
				if (j > 0)
				{
					double initialmass = current_state[6] < 1.0e-6 ? 1.0e-6 : current_state[6];

					//add the starting mass increment
					initialmass += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					state_at_beginning_of_phase[6] = initialmass * expfun;
					dmdvinf = -initialmass * 1000 / (options->IspChem * options->g0) * expfun;
				}
				else
				{
					double initialmass = options->maximum_mass;

					//add the starting mass increment
					initialmass += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					state_at_beginning_of_phase[6] = initialmass * expfun;
					dmdvinf = -options->maximum_mass * 1000 / (options->IspChem * options->g0) * expfun;
				}
			}
			else //"free" departure
			{
				state_at_beginning_of_phase[6] = (j == 0 ? options->maximum_mass : current_state[6]);

				//add the starting mass increment
				state_at_beginning_of_phase[6] += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

				dmdvinf = 0.0;
			}
		}
		else if (options->journey_departure_type[j] == 1) //depart from a parking orbit
		{
			double vinf;
			dVmag[dVindex] = EMTG::Astrodynamics::insertion_burn(state_at_beginning_of_phase+3, boundary1_state+3, Universe->bodies[boundary1_location_code-1].mu, Universe->bodies[boundary1_location_code-1].r_SOI, options->journey_departure_elements[j][0], options->journey_departure_elements[j][1], &vinf);
			*current_deltaV += dVmag[dVindex];
			C3_departure = vinf*vinf;
			++dVindex;

			if (boundary1_location_code < 0) //cannot do a simplified parking orbit departure unless the departure point is a body!
			{
				cout << "Cannot do a simplified parking orbit departure unless the departure point is a body. For other departure point types, use a 'direct' departure" << endl;
				throw 20;
			}

			if (j == 0) //if this is the first journey, no mass is assigned yet
				state_at_beginning_of_phase[6] = options->maximum_mass * exp(-dVmag[0] * 1000/ (options->IspDS * options->g0));
			else
				state_at_beginning_of_phase[6] = current_state[6] * exp(-dVmag[0] * 1000/ (options->IspDS * options->g0));
		}
	}
	else if (options->journey_departure_type[j] == 6)
	{
		//for first-phase zero-turn flybys
		//there are no decision variables or constraints

		//step 3.1 compute incoming v_infinity at flyby
		if (p == 0)
			this->locate_boundary_point(this->boundary1_location_code,
			options->journey_departure_type[j],
			true,
			Universe,
			boundary1_state,
			current_state + 3,
			*current_epoch,
			X,
			Xindex,
			F,
			Findex,
			G,
			Gindex,
			needG,
			j,
			p,
			options);

		for (int k = 0; k < 3; ++k)
			this->V_infinity_in(k) = current_state[k + 3] - boundary1_state[k + 3];
		this->V_infinity_out = this->V_infinity_in;
		
		this->flyby_outgoing_v_infinity = this->V_infinity_out.norm();
		this->C3_departure = this->flyby_outgoing_v_infinity*this->flyby_outgoing_v_infinity;

		this->flyby_altitude = 0.0;
		this->flyby_turn_angle = 0.0;
		this->dV_departure_magnitude = 0.0;

		//calculate the b-plane parameters, check the periapse altitude
		this->BoundaryR.assign_all(boundary1_state);
		this->BoundaryV.assign_all(boundary1_state + 3);

		dVdeparture[0] = 0.0;
		dVdeparture[1] = 0.0;
		dVdeparture[2] = 0.0;



		//store the state at the beginning of the phase, post-flyby
		for (int k = 0; k < 6; ++k)
		{
			this->state_at_beginning_of_phase[k] = current_state[k];
		}

		//store the mass
		if (p == 0)
			this->state_at_beginning_of_phase[6] = current_state[6] + options->journey_starting_mass_increment[j];
		else
			this->state_at_beginning_of_phase[6] = current_state[6];
	}
	else
	{
		//there is no alternate step 2

		//Step 3 (alternate): successive states start with a flyby
		//we need to process the flyby

		if (p == 0)
			locate_boundary_point(boundary1_location_code, options->journey_departure_type[j], true, Universe, boundary1_state, current_state+3, *current_epoch, X, Xindex, F, Findex, G, Gindex, false, j, p, options);
		BoundaryR.assign_all(boundary1_state);
		BoundaryV.assign_all(boundary1_state+3);

		//Step 3.1 extract the parameters of the flyby
		double rp_ratio = X[*Xindex];
		flyby_altitude = (rp_ratio - 1) * Body1->radius;
		b_plane_insertion_angle = X[*Xindex + 1];
		*Xindex += 2;

		//Step 3.2 process the flyby
		double flyby_orbit_energy;
		EMTG::Astrodynamics::unpowered_flyby(current_state+3, boundary1_state+3, Body1->mu, Body1->radius, Body1->r_SOI,
											rp_ratio, b_plane_insertion_angle, state_at_beginning_of_phase+3,
											&C3_departure, &flyby_turn_angle, &flyby_orbit_energy);

		//Step 3.2.1 calculate the b-plane parameters, check the periapse altitude
		for (int k = 0; k < 3; ++k)
		{
			V_infinity_in(k) = current_state[k+3] - boundary1_state[k+3];
			V_infinity_out(k) = state_at_beginning_of_phase[k+3] - boundary1_state[k+3];
		}

		//Step 3.3 apply the flyby hyperbola constraint
		F[*Findex] = (flyby_orbit_energy == 0 ? 0 : flyby_orbit_energy > 0 ? log10(flyby_orbit_energy) : -log10(-flyby_orbit_energy));
		++(*Findex);

		//Step 4: compute the position post-departure
		for (int k = 0; k < 3; ++k)
		{
			state_at_beginning_of_phase[k] = boundary1_state[k];
		}
		state_at_beginning_of_phase[6] = current_state[6];

		double dV_flyby[3];
		for (int k = 0; k < 3; ++k)
			dV_flyby[k] = state_at_beginning_of_phase[k] - boundary1_state[k];
		flyby_outgoing_v_infinity = V_infinity_in.norm();
	}

	

	//Step 5 propagate to the burn point
	//Step 5.1 extract parameters of the burn and flight time
	this->TOF = X[*Xindex];
	this->burn_index = X[*Xindex + 1];
	this->phase_start_epoch = *current_epoch;
	this->phase_end_epoch = *current_epoch + TOF;
	*Xindex += 2;
	
	//Step 5.2 compute the time before and after the burn
	this->time_before_burn = this->burn_index * this->TOF;
	this->time_after_burn = this->TOF - this->time_before_burn;

	//Step 5.3 propagate
	Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_beginning_of_phase,
												this->state_before_burn, 
												Universe->mu, 
												Universe->LU,
												this->time_before_burn,
												this->Kepler_F_Current,
												this->Kepler_Fdot_Current,
												this->Kepler_G_Current, 
												this->Kepler_Gdot_Current, 
												this->Kepler_Fdotdot_Current, 
												this->Kepler_Gdotdot_Current,
												this->Current_STM,
												false);

	this->state_before_burn[6] = this->state_at_beginning_of_phase[6];

	//Step 6: locate the second boundary point
	locate_boundary_point(boundary2_location_code,
							options->journey_arrival_type[j],
							false,
							Universe,
							boundary2_state,
							current_state+3, 
							*current_epoch + this->TOF, 
							X, 
							Xindex, 
							F, 
							Findex,
							G, 
							Gindex,
							false, 
							j,
							p,
							options);

	//Step 7: solve Lambert's problem to the right hand boundary point
#ifdef _EMTG_Lambert
	EMTG::Astrodynamics::Lambert (	this->state_before_burn,
									boundary2_state,
									time_after_burn,
									Universe->mu,
									1, 
									0, 
									lambert_v1, 
									lambert_v2);

	//check the angular momentum vector
	//we do not want to go retrograde because of a burn. It's OK to go retrograde because of a flyby, but we don't want to do it propulsively
	//similarly if we are already retrograde, the burn shouldn't drive us prograde
	//in other words, orbit direction switches should not be done by burns, only by flybys (this lets us pick Lambert "long way" or "short way"
	hz1 = this->state_before_burn[0]*this->state_before_burn[3]-this->state_before_burn[1]*this->state_before_burn[4]; //angular momentum before burn
	hz2 = this->state_before_burn[0]*lambert_v1[1]-this->state_before_burn[1]*lambert_v1[0]; //angular momentum after burn

	if (!(hz1 == hz2))
		EMTG::Astrodynamics::Lambert(	state_before_burn,
										boundary2_state,
										time_after_burn,
										Universe->mu,
										0,
										0,
										lambert_v1,
										lambert_v2);
#endif

	//Step 8: compute the state after the burn
	for (int k = 0; k < 3; ++k)
	{
		state_after_burn[k] = state_before_burn[k];
		state_after_burn[k+3] = lambert_v1[k];
		midcourse_dV[k] = state_after_burn[k+3] - state_before_burn[k+3];
	}
	dVmag[dVindex] = EMTG::math::norm(midcourse_dV, 3);
	state_after_burn[6] = state_before_burn[6] * exp(-dVmag[dVindex] * 1000/ (options->IspChem * options->g0));
	*current_deltaV += dVmag[dVindex];
	++dVindex;

	//Step 8.1 apply the minimum distance from the central body constraint
	//F[*Findex] = math::norm(state_before_burn, 3) / Universe->LU;
	//++(*Findex);

	//Step 9: Process the end of the phase and compute the final state vector. Is this an arrival?
	if (p == options->number_of_phases[j] - 1)
	{
		if (boundary2_location_code > 0) //ending at body
			dVmag[dVindex] = process_arrival(	lambert_v2,
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
			dVmag[dVindex] = process_arrival(	lambert_v2, 
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
		
		state_at_end_of_phase[6] = state_after_burn[6] * exp(-dVmag[dVindex] * 1000/ (options->IspChem * options->g0));
		*current_deltaV += dVmag[dVindex];
		++dVindex;	
	}
	else
	{
		state_at_end_of_phase[6] = state_after_burn[6];
	}
	for (int k = 0; k < 3; ++k)
	{
		state_at_end_of_phase[k] = boundary2_state[k];
		state_at_end_of_phase[k+3] = lambert_v2[k];
	}

	//Step 10: advance the current epoch
	*current_epoch += TOF;
	
	//Step 11: update the current state
	for (int k = 0; k < 7; ++k)
		current_state[k] = state_at_end_of_phase[k];

	return 0;
}

//output function
//return 0 if successful, 1 if failure
int MGA_DSM_phase::output(missionoptions* options, const double& launchdate, int j, int p, EMTG::Astrodynamics::universe* Universe, int* eventcount)
{
	
	//Step 1: store data that will be used for the printing
	double empty_vector[] = {0,0,0};
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

	string boundary1_name;
	string boundary2_name;

	switch (boundary1_location_code)
	{
		case -1:
			{
				boundary2_name = "free point";
				break;
			}
		case -2: //begin at hyperbolic periapse
			{
				boundary1_name = "Hyp-arrival";
				break;
			}
		default:
			boundary1_name = (Body1->name);
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
						(p == 0 && !(options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4 || options->journey_departure_type[j] == 6)) ? dVmag[0] : flyby_outgoing_v_infinity,
						-1,
						initial_Isp,
						-1,
						0,
						0,
						0);

	//*****************************************************************************
	//Next, all phases have a deep-space maneuver (DSM). Coast until the DSM
	double timestep = time_before_burn / options->num_timesteps;
	double output_state[7];
	output_state[6] = state_at_beginning_of_phase[6];
	
	for (int step = 0; step < options->num_timesteps; ++step)
	{
		//compute the current epoch
		double epoch = phase_start_epoch + timestep * (step + 0.5);

		//propagate the spacecraft
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_beginning_of_phase, 
													output_state,
													Universe->mu,
													Universe->LU,
													(epoch - this->phase_start_epoch), 
													this->Kepler_F_Current,
													this->Kepler_Fdot_Current, 
													this->Kepler_G_Current,
													this->Kepler_Gdot_Current,
													this->Kepler_Fdotdot_Current,
													this->Kepler_Gdotdot_Current, 
													this->Current_STM,
													false);

		//write the summary line
		write_summary_line(options,
							Universe,
							eventcount,
							epoch / 86400.0,
							"coast",
							"deep-space",
							timestep / 86400.0,
							-1,
							-1,
							-1,
							0,
							0,
							0,
							output_state,
							empty_vector,
							empty_vector,
							0,
							-1,
							-1,
							-1,
							0,
							0,
							0);
	}

	//then print the DSM
	write_summary_line(options,
						Universe,
						eventcount,
						(phase_start_epoch + time_before_burn) / 86400.0,
						"chem_burn",
						"deep-space",
						0,
						-1,
						-1,
						-1,
						atan2(midcourse_dV[1], midcourse_dV[2]),
						asin(midcourse_dV[2] / EMTG::math::norm(midcourse_dV, 3)),
						0,
						state_after_burn,
						midcourse_dV,
						empty_vector,
						EMTG::math::norm(midcourse_dV, 3),
						-1,
						options->IspChem,
						-1,
						0,
						0,
						0);

	//now print the coast after the DSM
	timestep = time_after_burn / options->num_timesteps;
	output_state[6] = state_after_burn[6];

	for (int step = 0; step < options->num_timesteps; ++step)
	{
		//compute the current epoch
		double epoch = phase_start_epoch + time_before_burn + timestep * (step + 0.5);

		//propagate the spacecraft
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_after_burn,
													output_state,
													Universe->mu,
													Universe->LU,
													timestep * (step + 0.5),
													this->Kepler_F_Current,
													this->Kepler_Fdot_Current, 
													this->Kepler_G_Current, 
													this->Kepler_Gdot_Current,
													this->Kepler_Fdotdot_Current,
													this->Kepler_Gdotdot_Current, 
													this->Current_STM,
													false);

		//write the summary line
		write_summary_line(options,
							Universe,
							eventcount,
							epoch / 86400.0,
							"coast",
							"deep-space",
							timestep / 86400.0,
							-1,
							-1,
							-1,
							0,
							0,
							0,
							output_state,
							empty_vector,
							empty_vector,
							0,
							-1,
							-1,
							-1,
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
		else if (options->journey_arrival_type[j] == 1 || options->journey_arrival_type[j] == 3)
			event_type = "rendezvous";
		else if (options->journey_arrival_type[j] == 2)
			event_type = "intercept";
		else if (options->journey_arrival_type[j] == 4)
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
		else if (options->journey_arrival_type[j] == 4)
		{
			dV_arrival_mag = 0;
			dVarrival[0] = 0;
			dVarrival[1] = 0;
			dVarrival[2] = 0;
		}
		else
		{
			if (p == 0)
			{
				dV_arrival_mag = dVmag[2];
			}
			else
			{
				dV_arrival_mag = dVmag[1];
			}
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
	}

	return 0;
}

//bounds calculation function
//return 0 if successful, 1 if failure
int MGA_DSM_phase::calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
{
	//this function calculates the upper and lower bounds for the decision and constraint vectors for MGA-DSM

	//create a prefix string with journey and phase information
	stringstream prefixstream;
	prefixstream << "j" << j << "p" << p << ": ";
	string prefix = prefixstream.str();
	int first_X_entry_in_phase = Xupperbounds->size();

	//first, we need to know if we are the first phase in the journey
	if (p == 0)
	{
		//do not encode time variables for phases starting in a flyby
		if (!(options->journey_departure_type[j] == 6))
		{

			//if we are the first phase, we also need to know if we are the first journey
			if (j == 0)
			{
				//if so, then the first decision variable is the launch epoch in MJD
				Xlowerbounds->push_back(options->launch_window_open_date + options->journey_wait_time_bounds[j][0]);
				Xupperbounds->push_back(options->launch_window_open_date + options->journey_wait_time_bounds[j][1]);
				Xdescriptions->push_back(prefix + "launch epoch (MJD)");
			}
			else
			{
				//if we are not the first journey, are we starting from the boundary of a sphere of influence?
				//if so, then there is no wait time. If not, then the first decision variable is the stay time at the first body in the journey (i.e. at the asteroid for sample return)
				if (!(options->sequence[j - 1][p + 1] == -1 || boundary1_location_code == -1))
				{
					Xlowerbounds->push_back(options->journey_wait_time_bounds[j][0]);
					Xupperbounds->push_back(options->journey_wait_time_bounds[j][1]);
					Xdescriptions->push_back(prefix + "stay time (days)");
				}
			}

			if (boundary1_location_code == -1) //if this boundary point is at a free point in space, with the various elements either fixed or free
			{
				vector<string> CartesianElementNames;
				CartesianElementNames.push_back("x (km)");
				CartesianElementNames.push_back("y (km)");
				CartesianElementNames.push_back("z (km)");
				CartesianElementNames.push_back("xdot (km/s)");
				CartesianElementNames.push_back("ydot (km/s)");
				CartesianElementNames.push_back("zdot (km/s)");

				vector<string> ClassicalOrbitElementNames;
				ClassicalOrbitElementNames.push_back("SMA (km)");
				ClassicalOrbitElementNames.push_back("ECC (km)");
				ClassicalOrbitElementNames.push_back("INC (rad)");
				ClassicalOrbitElementNames.push_back("RAAN (rad)");
				ClassicalOrbitElementNames.push_back("AOP (rad)");
				ClassicalOrbitElementNames.push_back("TA (rad)");

				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_departure_elements_vary_flag[j][k])
					{
						Xlowerbounds->push_back(options->journey_departure_elements_bounds[j][k][0]);
						Xupperbounds->push_back(options->journey_departure_elements_bounds[j][k][1]);
						if (options->journey_departure_elements_type[j])
							Xdescriptions->push_back(prefix + " left boundary point " + ClassicalOrbitElementNames[k]);
						else
							Xdescriptions->push_back(prefix + " left boundary point " + CartesianElementNames[k]);
					}
				}

				//if it is possible for the optimizer to select a point inside the exclusion zone of the central body
				//then there must be a nonlinear constraint to prevent this
				if (options->journey_departure_elements_type[j]) //classical orbit elements
				{
					//this constraint is applied if we are varying SMA or ECC
					if (options->journey_departure_elements_vary_flag[j][0] || options->journey_departure_elements_vary_flag[j][1])
					{
						Flowerbounds->push_back(0.0);
						Fupperbounds->push_back(math::LARGE);
						Fdescriptions->push_back(prefix + " left boundary central body exclusion radius constraint");

						//this constraint has derivatives with respect to SMA, ECC, and TA
						//only create a derivative entry with respect to an orbit element if that element is being varied
						if (options->journey_departure_elements_vary_flag[j][0]) //SMA
						{
							for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							{
								if ((*Xdescriptions)[entry].find("left boundary point SMA") < 1024)
								{
									iGfun->push_back(Fdescriptions->size() - 1);
									jGvar->push_back(entry);
									stringstream EntryNameStream;
									EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									Gdescriptions->push_back(EntryNameStream.str());
									left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								}
							}
						}
						if (options->journey_departure_elements_vary_flag[j][1]) //ECC
						{
							for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							{
								if ((*Xdescriptions)[entry].find("left boundary point ECC") < 1024)
								{
									iGfun->push_back(Fdescriptions->size() - 1);
									jGvar->push_back(entry);
									stringstream EntryNameStream;
									EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									Gdescriptions->push_back(EntryNameStream.str());
									left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								}
							}
						}
						if (options->journey_departure_elements_vary_flag[j][5]) //TA
						{
							for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							{
								if ((*Xdescriptions)[entry].find("left boundary point TA") < 1024)
								{
									iGfun->push_back(Fdescriptions->size() - 1);
									jGvar->push_back(entry);
									stringstream EntryNameStream;
									EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									Gdescriptions->push_back(EntryNameStream.str());
									left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								}
							}
						}
					}
				}
				else //cartesian orbit elements
				{
					//this constraint is applied if we are varying x, y, or z
					if (options->journey_departure_elements_vary_flag[j][0] || options->journey_departure_elements_vary_flag[j][1] || options->journey_departure_elements_vary_flag[j][2])
					{
						Flowerbounds->push_back(0.0);
						Fupperbounds->push_back(math::LARGE);
						Fdescriptions->push_back(prefix + " left boundary central body exclusion radius constraint");
					}

					//this constraint has derivatives with respect to x, y, and z
					//only create a derivative entry with respect to an orbit element if that element is being varied
					if (options->journey_departure_elements_vary_flag[j][0]) //x
					{
						for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						{
							if ((*Xdescriptions)[entry].find("left boundary point x (km)") < 1024)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(entry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								Gdescriptions->push_back(EntryNameStream.str());
								left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							}
						}
					}
					if (options->journey_departure_elements_vary_flag[j][1]) //y
					{
						for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						{
							if ((*Xdescriptions)[entry].find("left boundary point y (km)") < 1024)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(entry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								Gdescriptions->push_back(EntryNameStream.str());
								left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							}
						}
					}
					if (options->journey_departure_elements_vary_flag[j][2]) //z
					{
						for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						{
							if ((*Xdescriptions)[entry].find("left boundary point z (km)") < 1024)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(entry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								Gdescriptions->push_back(EntryNameStream.str());
								left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							}
						}
					}
				}
			}
			//if we are starting at periapse of an arrival hyperbola, we must choose the orbit elements of the initial capture orbit
			//this will be done as [rp, ra, inc, raan, aop, 0]
			//note the true anomaly is always zero because we always capture into the periapse of our desired orbit
			else if (boundary1_location_code == -2)
			{
				Xlowerbounds->push_back(Universe->minimum_safe_distance);
				Xupperbounds->push_back(10 * Universe->minimum_safe_distance);
				Xdescriptions->push_back(prefix + "left boundary periapse distance");

				Xlowerbounds->push_back(Universe->minimum_safe_distance);
				Xupperbounds->push_back(Universe->r_SOI);
				Xdescriptions->push_back(prefix + "left boundary apoapse distance");

				Xlowerbounds->push_back(-2 * EMTG::math::PI);
				Xupperbounds->push_back(2 * EMTG::math::PI);
				Xdescriptions->push_back(prefix + "left boundary inclination");

				Xlowerbounds->push_back(-2 * EMTG::math::PI);
				Xupperbounds->push_back(2 * EMTG::math::PI);
				Xdescriptions->push_back(prefix + "left boundary RAAN");

				Xlowerbounds->push_back(-2 * EMTG::math::PI);
				Xupperbounds->push_back(2 * EMTG::math::PI);
				Xdescriptions->push_back(prefix + "left boundary AOP");

				Flowerbounds->push_back(-EMTG::math::LARGE);
				Fupperbounds->push_back(0.0);
				Fdescriptions->push_back(prefix + "left boundary rp < ra");

				Flowerbounds->push_back(-EMTG::math::LARGE);
				Fupperbounds->push_back(0.0);
				Fdescriptions->push_back(prefix + "incoming orbit must be a hyperbola, i.e. v_periapse > v_escape");
			}

			//then we have three variables to parameterize the departure velocity, only for phases where there is a departure velocity
			//we do NOT encode a departure velocity for successive phases which start at the boundary of an SOI (that would be too many degrees of freedom)
			//we do NOT encode a departure velocity for phases which start from an inbound hyperbola because we are encoding the capture state instead
			if ((j == 0 || !(options->journey_departure_type[j] == 3)) && !(boundary1_location_code == -2))
			{
				//First, we have outgoing velocity
				Xlowerbounds->push_back(options->journey_initial_impulse_bounds[j][0]);
				Xupperbounds->push_back(options->journey_initial_impulse_bounds[j][1]);
				Xdescriptions->push_back(prefix + "magnitude of outgoing velocity asymptote");

				//then we have the angles defining the departure asymptote
				Xlowerbounds->push_back(-math::PI);
				Xupperbounds->push_back(math::PI);
				Xdescriptions->push_back(prefix + "RA of departure asymptote");

				if (j == 0 && boundary1_location_code > 0) //if this is the first journey and we are leaving from a planet, i.e. if this is a launch
				{
					Xlowerbounds->push_back(options->DLA_bounds[0] * math::PI / 180.0);
					Xupperbounds->push_back(options->DLA_bounds[1] * math::PI / 180.0);
				}
				else
				{
					Xlowerbounds->push_back(-math::PI / 2.0);
					Xupperbounds->push_back(math::PI / 2.0);
				}

				Xdescriptions->push_back(prefix + "DEC of departure asymptote");
			}
		}
	}
	else
	{
		//if we are NOT the first phase then we need only two variables, defining the flyby geometry
		//we must also implement the flyby hyperbolic constraint

		//flyby hyperbolic constraint
		Flowerbounds->push_back(0.0);
		Fupperbounds->push_back(math::LARGE);
		Fdescriptions->push_back(prefix + "flyby orbit energy plus 10% fudge factor must be positive");
		//constraint has a derivative with respect to all preceeding variables
		for (size_t entry = 0; entry < Xdescriptions->size() - 2; ++entry)
		{
			iGfun->push_back(Fdescriptions->size() - 1);
			jGvar->push_back(entry);
			stringstream EntryNameStream;
			EntryNameStream << "Derivative of " << prefix << " flyby orbit energy constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
			Gdescriptions->push_back(EntryNameStream.str());
		}

		//altitude ratio (as a fraction of minimum safe altitude)
		Xlowerbounds->push_back(1.0 + Universe->bodies[boundary1_location_code-1].minimum_safe_flyby_altitude / Universe->bodies[boundary1_location_code-1].radius);
		if (Universe->bodies[boundary1_location_code-1].mass < 1.0e+25)
			Xupperbounds->push_back(10.0);
		else
			Xupperbounds->push_back(300.0);
		Xdescriptions->push_back(prefix + "flyby altitude ratio");

		//b-plane insertion angle
		Xlowerbounds->push_back(-math::PI);
		Xupperbounds->push_back(math::PI);
		Xdescriptions->push_back(prefix + "b-plane insertion angle");
	}
	
	//**************************************************************************
	//next, we need to encode the phase flight time
	calcbounds_flight_time(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, synodic_periods, j, p, Universe, options);

	//all MGA-DSM phases encode a burn index
	Xlowerbounds->push_back(0.05);
	Xupperbounds->push_back(0.95);
	Xdescriptions->push_back(prefix + "burn index");

	/*//all MGA-DSM phases encode a no-collision constraint
	Flowerbounds->push_back(Universe->minimum_safe_distance / Universe->LU);
	Fupperbounds->push_back(Universe->r_SOI / Universe->LU);
	Fdescriptions->push_back(prefix + "minimum safe distance constraint");
	//constraint has a derivative with respect to all preceeding variables
	for (size_t entry = 0; entry < Xdescriptions->size(); ++entry)
	{
		iGfun->push_back(Fdescriptions->size() - 1);
		jGvar->push_back(entry);
		stringstream EntryNameStream;
		EntryNameStream << "Derivative of " << prefix << " minimum safe distance constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
		Gdescriptions->push_back(EntryNameStream.str());
	}*/

	//******************
	//if we are the last journey, then encode any variables necessary for the right hand boundary condition
	if (p == (options->number_of_phases[j] - 1))
	{
		if  (boundary2_location_code == -1) //if this boundary point is at a free point in space, with the various elements either fixed or free
			{
				vector<string> CartesianElementNames;
				CartesianElementNames.push_back("x (km)");
				CartesianElementNames.push_back("y (km)");
				CartesianElementNames.push_back("z (km)");
				CartesianElementNames.push_back("xdot (km/s)");
				CartesianElementNames.push_back("ydot (km/s)");
				CartesianElementNames.push_back("zdot (km/s)");

				vector<string> ClassicalOrbitElementNames;
				ClassicalOrbitElementNames.push_back("SMA (km)");
				ClassicalOrbitElementNames.push_back("ECC (km)");
				ClassicalOrbitElementNames.push_back("INC (rad)");
				ClassicalOrbitElementNames.push_back("RAAN (rad)");
				ClassicalOrbitElementNames.push_back("AOP (rad)");
				ClassicalOrbitElementNames.push_back("TA (rad)");

				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_arrival_elements_vary_flag[j][k])
					{
						Xlowerbounds->push_back(options->journey_arrival_elements_bounds[j][k][0]);
						Xupperbounds->push_back(options->journey_arrival_elements_bounds[j][k][1]);
						if (options->journey_arrival_elements_type[j])
							Xdescriptions->push_back(prefix + " right boundary point " + ClassicalOrbitElementNames[k]);
						else
							Xdescriptions->push_back(prefix + " right boundary point " + CartesianElementNames[k]);
					}
				}

				//if it is possible for the optimizer to select a point inside the exclusion zone of the central body
				//then there must be a nonlinear constraint to prevent this
				if (options->journey_arrival_elements_type[j]) //classical orbit elements
				{
					//this constraint is applied if we are varying SMA or ECC
					if (options->journey_arrival_elements_vary_flag[j][0] || options->journey_arrival_elements_vary_flag[j][1])
					{
						Flowerbounds->push_back(0.0);
						Fupperbounds->push_back(math::LARGE);
						Fdescriptions->push_back(prefix + " right boundary central body exclusion radius constraint");

						//this constraint has derivatives with respect to SMA, ECC, and TA
						//only create a derivative entry with respect to an orbit element if that element is being varied
						if (options->journey_arrival_elements_vary_flag[j][0]) //SMA
						{
							for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							{
								if ((*Xdescriptions)[entry].find("right boundary point SMA") < 1024)
								{
									iGfun->push_back(Fdescriptions->size() - 1);
									jGvar->push_back(entry);
									stringstream EntryNameStream;
									EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									Gdescriptions->push_back(EntryNameStream.str());
									right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								}
							}
						}
						if (options->journey_arrival_elements_vary_flag[j][1]) //ECC
						{
							for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							{
								if ((*Xdescriptions)[entry].find("right boundary point ECC") < 1024)
								{
									iGfun->push_back(Fdescriptions->size() - 1);
									jGvar->push_back(entry);
									stringstream EntryNameStream;
									EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									Gdescriptions->push_back(EntryNameStream.str());
									right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								}
							}
						}
						if (options->journey_arrival_elements_vary_flag[j][5]) //TA
						{
							for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							{
								if ((*Xdescriptions)[entry].find("right boundary point TA") < 1024)
								{
									iGfun->push_back(Fdescriptions->size() - 1);
									jGvar->push_back(entry);
									stringstream EntryNameStream;
									EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									Gdescriptions->push_back(EntryNameStream.str());
									right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								}
							}
						}
					}
				}
				else //cartesian orbit elements
				{
					//this constraint is applied if we are varying x, y, or z
					if (options->journey_arrival_elements_vary_flag[j][0] || options->journey_arrival_elements_vary_flag[j][1] || options->journey_arrival_elements_vary_flag[j][2])
					{
						Flowerbounds->push_back(0.0);
						Fupperbounds->push_back(math::LARGE);
						Fdescriptions->push_back(prefix + " right boundary central body exclusion radius constraint");
					}

					//this constraint has derivatives with respect to x, y, and z
					//only create a derivative entry with respect to an orbit element if that element is being varied
					if (options->journey_arrival_elements_vary_flag[j][0]) //x
					{
						for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						{
							if ((*Xdescriptions)[entry].find("right boundary point x (km)") < 1024)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(entry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								Gdescriptions->push_back(EntryNameStream.str());
								right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							}
						}
					}
					if (options->journey_arrival_elements_vary_flag[j][1]) //y
					{
						for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						{
							if ((*Xdescriptions)[entry].find("right boundary point y (km)") < 1024)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(entry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								Gdescriptions->push_back(EntryNameStream.str());
								right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							}
						}
					}
					if (options->journey_arrival_elements_vary_flag[j][2]) //z
					{
						for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						{
							if ((*Xdescriptions)[entry].find("right boundary point z (km)") < 1024)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(entry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								Gdescriptions->push_back(EntryNameStream.str());
								right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							}
						}
					}
				}
			}

		//variables and constraints necessary to process the arrival type
		if (options->journey_arrival_type[j] == 2) //flyby with bounded v-infinity
		{
			Flowerbounds->push_back((options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1])*(options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1]) - 1);
			Fupperbounds->push_back(0.0);
			Fdescriptions->push_back(prefix + "arrival C3 constraint");

			//Jacobian entry for a bounded v-infinity intercept
			//this is a nonlinear constraint dependent on all previous variables
			for (int entry = Xdescriptions->size() - 1; entry >= 0; --entry)
			{
				iGfun->push_back(Fdescriptions->size() - 1);
				jGvar->push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
				Gdescriptions->push_back(EntryNameStream.str());
			}
		}
	}


	return 0;
}

void MGA_DSM_phase::output_GTOC7_format(missionoptions* options, EMTG::Astrodynamics::universe* Universe, const std::string& GTOC_output_file, int j, int p)
{
	//mothership file format
	//state before and after all events
	ofstream GTOC7file;

	if (j == 0)
	{
		GTOC7file.open(GTOC_output_file.c_str(), ios::trunc);
		GTOC7file << "# MOTHER SHIP:         X" << endl;
	}
	else
		GTOC7file.open(GTOC_output_file.c_str(), ios::app);

	

	if (j == 0)
	{
		GTOC7file << "# EVENT NUMBER:         " << 1 << endl;
		GTOC7file << "# DESCRIPTION: DEPARTURE" << endl;
		GTOC7file << "#  Time (MJD)             x (km)                 y (km)                 z (km)                 vx (km/s)              vy (km/s)              vz (km/s)              mass (kg)" << endl;
		GTOC7file << " ";
		GTOC7file.precision(14);
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->phase_start_epoch / 86400.0, 2) << " ";
		for (int k = 0; k < 7; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_beginning_of_phase[k], 2) << " ";
		GTOC7file << endl;
	}

	GTOC7file << "# EVENT NUMBER:         " << j+2 << endl;
	GTOC7file << "# DESCRIPTION: Impulse " << j+1 << endl;
	GTOC7file << "#  Time (MJD)             x (km)                 y (km)                 z (km)                 vx (km/s)              vy (km/s)              vz (km/s)              mass (kg)" << endl;
	GTOC7file << " ";
	GTOC7file.precision(14);
	GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string((this->phase_start_epoch + this->time_before_burn) / 86400.0, 2) << " ";
	for (int k = 0; k < 7; ++k)
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_before_burn[k], 2) << " ";
	GTOC7file << endl;
	GTOC7file << " ";
	GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string((this->phase_start_epoch + this->time_before_burn) / 86400.0, 2) << " ";
	for (int k = 0; k < 7; ++k)
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_after_burn[k], 2) << " ";
	GTOC7file << endl;

	if (options->journey_arrival_type[j] == 1)
	{
		GTOC7file << "# EVENT NUMBER:         " << 1 << endl;
		GTOC7file << "# DESCRIPTION: Impulse " << j+2 << endl;
		GTOC7file << "#  Time (MJD)             x (km)                 y (km)                 z (km)                 vx (km/s)              vy (km/s)              vz (km/s)              mass (kg)" << endl;
		GTOC7file << " ";
		GTOC7file.precision(14);
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string((this->phase_end_epoch) / 86400.0, 2) << " ";
		for (int k = 0; k < 3; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_end_of_phase[k], 2) << " ";
		for (int k = 0; k < 3; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_end_of_phase[k+3], 2) << " ";
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_end_of_phase[6] * exp(this->dVmag.back() * 1000 / (options->IspChem * options->g0)), 2) << " ";
		GTOC7file << endl;
		GTOC7file << " ";
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string((this->phase_end_epoch) / 86400.0, 2) << " ";
		for (int k = 0; k < 3; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_end_of_phase[k], 2) << " ";
		for (int k = 0; k < 3; ++k)
			GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_end_of_phase[k + 3] + this->dVarrival[k], 2) << " ";
		GTOC7file << EMTG::string_utilities::convert_number_to_formatted_string(this->state_at_end_of_phase[6], 2) << " ";
		GTOC7file << endl;
	}
}

} /* namespace EMTG */
