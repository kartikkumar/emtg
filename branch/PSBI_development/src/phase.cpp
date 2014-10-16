/*
 * phase.cpp
 *
 *  Created on: Jul 15, 2012
 *      Author: Jacob
 */

#include <iostream>
#include <fstream>
#include <sstream>

#include "phase.h"
#include "journey.h"
#include "Astrodynamics.h"
#include "EMTG_math.h"
#include "mjd_to_mdyhms.h"
#include "universe.h"
#include "kepler_lagrange_laguerre_conway.h"

#include "SpiceUsr.h"

using namespace std;

namespace EMTG {

	phase::phase() :
		boundary1_location_code(0),
		boundary2_location_code(0),
		TOF(0),
		phase_end_epoch(0),
		phase_start_epoch(0),
		flyby_turn_angle(0),
		flyby_altitude(0),
		RA_departure(0),
		RA_arrival(0),
		DEC_departure(0),
		DEC_arrival(0),
		C3_departure(0),
		C3_arrival(0)
	{
		// default constructor is never used

	}

	phase::~phase()
	{
	// default destructor is never used (I think it is superceded by the daughter class destructors)
	}

	void phase::process_left_boundary_condition(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, const int& needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
	{
		double vinf_out;
		double vinf_in;

		//set pointers to the bodies
		if (this->boundary1_location_code > 0)
			this->Body1 = &Universe->bodies[boundary1_location_code - 1];
		if (this->boundary2_location_code > 0)
			this->Body2 = &Universe->bodies[boundary2_location_code - 1];

		//we need to know if we are the first phase of the journey and the journey does not start with a flyby
		if (p == 0 && !(options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 6))
		{
			//Step 1: extract the journey start epoch
			//if this is the first journey, the first decision variable is the starting epoch
			//otherwise, the first decision variable is the journey wait time
			if (j == 0)
			{
				*current_epoch = X[*Xindex];
				this->phase_wait_time = X[*Xindex] - options->launch_window_open_date;
				++(*Xindex);
			}
			else
			{
				*current_epoch += X[*Xindex];
				this->phase_wait_time = X[*Xindex];
				++(*Xindex);
			}

			//Step 2: locate the first body
			this->locate_boundary_point(this->boundary1_location_code,
										options->journey_departure_type[j],
										true,
										Universe,
										boundary1_state,
										current_state+3,
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

			//Step 3: compute the departure asymptote
			//Step 3.1 extract the departure parameters
			if (!(options->journey_departure_type[j] == 5 || options->journey_departure_type[j] == 2)) //if this journey starts with an impulse
			{
				vinf_out = X[*Xindex];
				this->C3_departure = vinf_out*vinf_out;
				this->RA_departure = X[*Xindex + 1];
				this->DEC_departure = X[*Xindex + 2];
				*Xindex += 3;

				//Step 3.2 compute the outgoing velocity vector
				double cDEC = cos(DEC_departure);
				this->V_infinity_out(0) = vinf_out * cos(this->RA_departure)*cDEC;
				this->V_infinity_out(1) = vinf_out * sin(this->RA_departure)*cDEC;
				this->V_infinity_out(2) = vinf_out * sin(this->DEC_departure);			

				//Step 4: compute the state post-departure
				//Step 4.1: compute position and velocity
				for (int k = 0; k < 3; ++k)
				{
					this->state_at_beginning_of_phase[k] = boundary1_state[k];
					this->state_at_beginning_of_phase[k+3] = boundary1_state[k+3] + this->V_infinity_out(k);
				}
				//Step 4.2 compute the mass post-departure
				if (options->journey_departure_type[j] == 0 || options->journey_departure_type[j] == 2) //this is a direct insertion or launch
				{
					this->dV_departure_magnitude = vinf_out;

					if (options->include_initial_impulse_in_cost)
						*current_deltaV += vinf_out;
			
					if (j == 0) //this is a launch
					{

						//figure out how much mass we can launch
						if (options->LV_type > 0 || options->LV_type == -2)
						{
							double launch_mass;
							EMTG::Astrodynamics::find_mass_to_orbit(this->C3_departure,
																	this->DEC_departure * 180.0 / math::PI,
																	options->LV_type,
																	&launch_mass,
																	&dmdvinf,
																	options);

							this->state_at_beginning_of_phase[6] = launch_mass > options->maximum_mass ? options->maximum_mass : launch_mass;
						
							//apply the mass margin
							this->state_at_beginning_of_phase[6] *= 1.0 - options->LV_margin;

							//convert from dm/dC3 to dm/dvinf
							this->dmdvinf *= (2 * vinf_out)*(1.0 - options->LV_margin);
						}
						else if (options->LV_type == 0)
						{
							this->state_at_beginning_of_phase[6] = options->maximum_mass;

							//add the starting mass increment
							this->state_at_beginning_of_phase[6] += this->journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

							this->dmdvinf = 0.0;
						}
						else //chemical burn using departure stage engine
						{
							double expfun = exp(-vinf_out * 1000 / (options->IspDS * options->g0));

							double initialmass = (options->maximum_mass + this->journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j]) * expfun;

							this->state_at_beginning_of_phase[6] = initialmass;
							this->dmdvinf = -initialmass * 1000 / (options->IspDS * options->g0) * expfun;
						}
					}
					else if (options->journey_departure_type[j] == 0)
					{
						double expfun = exp(-vinf_out * 1000 / (options->IspChem * options->g0));
						if (j > 0)
						{
							double initialmass = (current_state[6] + this->journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j]) * expfun;

							this->state_at_beginning_of_phase[6] = initialmass;
							this->dmdvinf = -initialmass * 1000 / (options->IspChem * options->g0) * expfun;
						}
						else
						{
							double initialmass = (current_state[6] + this->journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j]) * expfun;

							this->state_at_beginning_of_phase[6] = initialmass;
							this->dmdvinf = -options->maximum_mass * 1000 / (options->IspChem * options->g0) * expfun;
						}
					}
					else //"free" departure
					{
						state_at_beginning_of_phase[6] = (j == 0 ? options->maximum_mass : current_state[6]);

						//add the starting mass increment
						this->state_at_beginning_of_phase[6] += this->journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];
						this->C3_departure = 0.0;
						this->RA_departure = 0.0;
						this->DEC_departure = 0.0;

						this->dmdvinf = 0.0;
					}
				}
				else if (options->journey_departure_type[j] == 1) //depart from a parking orbit
				{
					double vinf;
					this->dV_departure_magnitude = EMTG::Astrodynamics::insertion_burn(state_at_beginning_of_phase+3,
																						boundary1_state+3,
																						Body1->mu,
																						Body1->r_SOI,
																						options->journey_departure_elements[j][0],
																						options->journey_departure_elements[j][1],
																						&vinf);
					*current_deltaV += this->dV_departure_magnitude;
					this->C3_departure = vinf*vinf;

					if (this->boundary1_location_code < 0) //cannot do a simplified parking orbit departure unless the departure point is a body!
					{
						cout << "Cannot do a simplified parking orbit departure unless the departure point is a body. For other departure point types, use a 'direct' departure" << endl;
						throw 20;
					}

					if (j == 0) //if this is the first journey, no mass is assigned yet
						this->state_at_beginning_of_phase[6] = options->maximum_mass * exp(-this->dV_departure_magnitude * 1000/ (options->IspDS * options->g0));
					else
						this->state_at_beginning_of_phase[6] = current_state[6] * exp(-this->dV_departure_magnitude * 1000/ (options->IspDS * options->g0));
				}

				if (j == 0 && options->allow_initial_mass_to_vary)
				{
					//if we have enabled varying the initial mass, then pass through a mass multiplier
					this->mission_initial_mass_multiplier = X[*Xindex];
					++(*Xindex);
					this->unscaled_phase_initial_mass = this->state_at_beginning_of_phase[6];
					this->state_at_beginning_of_phase[6] *= this->mission_initial_mass_multiplier;
				}
			}//end code for journeys that start with an impulse
			else if (options->journey_departure_type[j] == 2)//free direct departure
			{
				//no initial impulse
				this->C3_departure = 0;
				this->RA_departure = 0;
				this->DEC_departure = 0;
				this->dV_departure_magnitude = 0.0;

				//Step 3.2 compute the outgoing velocity vector
				this->V_infinity_out.assign_zeros();

				
				//*******************************************************
				//Step 4: compute the state post-departure

				for (int k = 0; k < 6; ++k)
					this->state_at_beginning_of_phase[k] = boundary1_state[k];

				double initialmass = j == 0 ? options->maximum_mass : current_state[6];

				//add the starting mass increment
				initialmass += this->journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

				this->state_at_beginning_of_phase[6] = initialmass;
				this->dmdvinf = 0.0;

				if (j == 0 && options->allow_initial_mass_to_vary)
				{
					//if we have enabled varying the initial mass, then pass through a mass multiplier
					this->mission_initial_mass_multiplier = X[*Xindex];
					++(*Xindex);
					this->unscaled_phase_initial_mass = this->state_at_beginning_of_phase[6];
					this->state_at_beginning_of_phase[6] *= this->mission_initial_mass_multiplier;
				}
			}
			else if (options->journey_departure_type[j] == 5)//for journeys starting with a spiral
			{
				//journeys which start from a spiral have no initial impulse
				this->C3_departure = 0;
				this->RA_departure = 0;
				this->DEC_departure = 0;

				//Step 3.2 compute the outgoing velocity vector
				this->V_infinity_out.assign_zeros();

				//*******************************************************
				//Step 4: compute the state post-departure

				//Step 4.1: compute the necessary inputs to the spiral model

				//Step 4.1.1 if we are using a VSI thruster then we need to extract the spiral Isp from the decision vector 
				if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
				{
					this->spiral_escape_Isp = X[*Xindex];
					++(*Xindex);
				}

				//Step 4.1.2 first we will need to compute the mass before the spiral
				//if this is the first journey, the initial mass is the maximum mass modified by the initial mass multiplier
				if (j == 0)
				{
					if (options->allow_initial_mass_to_vary)
					{
						this->mission_initial_mass_multiplier = X[*Xindex];
						++(*Xindex);
					}
					else
						this->mission_initial_mass_multiplier = 1.0;

					this->unscaled_phase_initial_mass = options->maximum_mass;
					this->spiral_escape_state_before_spiral[6] = this->unscaled_phase_initial_mass * this->mission_initial_mass_multiplier;
				}
				//if this is not the first journey, the initial mass is the pass-through "current mass"
				//and if applicable, the journey initial mass increment is added
				else
				{
					this->spiral_escape_state_before_spiral[6] = current_state[6] + this->journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];
				}
				//if this is the first journey, the mass before the spiral is the user-defined maximum mass modified by the initial mass multiplier
				//if this is NOT the first journey, the mass before the spiral is the mass from the end of the previous journey
				this->spiral_escape_mass_before = (j == 0 ? options->maximum_mass : current_state[6]);

				//Step 4.1.3 the rest of the state before the spiral is equal to the boundary state
				for (int k = 0; k < 6; ++k)
					this->spiral_escape_state_before_spiral[k] = boundary1_state[k];

				

				//Step 4.1.4 compute the starting body's distance from the Sun at the beginning of the spiral
				double distance_from_sun_in_AU;
				if (!(Universe->central_body_SPICE_ID == 10))
				{
					double central_body_state_in_km[6];
					double position_relative_to_sun_in_AU[3];
					Universe->locate_central_body(*current_epoch, central_body_state_in_km, options);

					position_relative_to_sun_in_AU[0] = (central_body_state_in_km[0] + this->spiral_escape_state_before_spiral[0]) / options->AU;
					position_relative_to_sun_in_AU[1] = (central_body_state_in_km[1] + this->spiral_escape_state_before_spiral[1]) / options->AU;
					position_relative_to_sun_in_AU[2] = (central_body_state_in_km[2] + this->spiral_escape_state_before_spiral[2]) / options->AU;

					distance_from_sun_in_AU = math::norm(position_relative_to_sun_in_AU, 3);
				}
				else //if we are orbiting the sun, it's quite silly to look up the position of the sun relative to itself, so don't bother
					distance_from_sun_in_AU = math::norm(boundary1_state, 3)/Universe->LU;

				//Step 4.2: process the spiral
				if (options->spiral_model_type == 0)
				{
					Astrodynamics::Battin_spiral(this->spiral_escape_state_before_spiral[6],
														options->journey_escape_spiral_starting_radius[j],
														distance_from_sun_in_AU,
														*current_epoch - X[0],
														this->spiral_escape_Isp,
														this->spiral_escape_thrust,
														this->spiral_escape_mdot,
														this->spiral_escape_number_of_engines,
														this->spiral_escape_power,
														this->spiral_escape_active_power,
														this->Body1->mu,
														(options->derivative_type > 1 && needG),
														options,
														&this->spiral_escape_time,
														&this->spiral_escape_mass_after,
														&this->spiral_escape_dv,
														&this->spiral_escape_dm_after_dm_before,
														&this->spiral_escape_dt_spiral_dm_before,
														&this->spiral_escape_dm_after_dIsp,
														&this->spiral_escape_dt_spiral_dIsp);
				}
				else if (options->spiral_model_type == 1)
				{
					Astrodynamics::Edelbaum_spiral(this->spiral_escape_state_before_spiral[6],
													options->journey_escape_spiral_starting_radius[j],
													distance_from_sun_in_AU,
													Body1->r_SOI,
													*current_epoch - X[0],
													this->spiral_escape_Isp,
													this->spiral_escape_thrust,
													this->spiral_escape_mdot,
													this->spiral_escape_number_of_engines,
													this->spiral_escape_power,
													this->spiral_escape_active_power,
													this->Body1->mu,
													(options->derivative_type > 1 && needG),
													options,
													&this->spiral_escape_time,
													&this->spiral_escape_mass_after,
													&this->spiral_escape_dv,
													&this->spiral_escape_dm_after_dm_before,
													&this->spiral_escape_dt_spiral_dm_before,
													&this->spiral_escape_dm_after_dIsp,
													&this->spiral_escape_dt_spiral_dIsp);
				}

				//Step 4.3 construct the post-spiral state
				//Step 4.3.1 advance time
				*current_epoch += this->spiral_escape_time;

				//Step 4.3.2 find the position of the body at the new phase start time and store it in state_at_beginning_of_phase
				this->locate_boundary_point(this->boundary1_location_code,
											options->journey_departure_type[j],
											true,
											Universe,
											this->state_at_beginning_of_phase,
											current_state+3,
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

				//Step 4.3.3 fill in the new mass
				this->state_at_beginning_of_phase[6] = this->spiral_escape_mass_after;
			}

		}//end departure code
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

			//Step 3.2 extract V_infinity_out from decision vector and apply equal v-infinity constraint
			for (int k = 0; k < 3; ++k)
			{
				this->V_infinity_out(k) = X[*Xindex];
				++(*Xindex);
				
				double v_infinity_in_k = this->V_infinity_in(k);

				F[*Findex] = (this->V_infinity_out(k) - v_infinity_in_k) / v_infinity_in_k;
				++(*Findex);

				if (options->derivative_type > 0 && needG)
				{
					G[flyby_velocity_magnitude_constraint_G_indices[k * 2]] = 1.0 / v_infinity_in_k * flyby_constraints_X_scale_ranges[k * 2];
					G[flyby_velocity_magnitude_constraint_G_indices[k * 2 + 1]] = -1.0 / v_infinity_in_k * flyby_constraints_X_scale_ranges[k * 2 + 1];
				}
			}

			this->flyby_outgoing_v_infinity = this->V_infinity_out.norm();
			this->C3_departure = this->flyby_outgoing_v_infinity*this->flyby_outgoing_v_infinity;

			this->flyby_altitude = 0.0;
			this->flyby_turn_angle = 0.0;
			this->dV_departure_magnitude = 0.0;

			//calculate the b-plane parameters, check the periapse altitude
			this->BoundaryR.assign_all(boundary1_state);
			this->BoundaryV.assign_all(boundary1_state + 3);



			//store the state at the beginning of the phase, post-flyby
			for (int k = 0; k < 3; ++k)
			{
				this->state_at_beginning_of_phase[k] = boundary1_state[k];
				this->state_at_beginning_of_phase[k + 3] = boundary1_state[k + 3] + V_infinity_out(k);
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
		
			//step 3.1 compute incoming v_infinity at flyby
			if (p == 0)
				this->locate_boundary_point(this->boundary1_location_code,
											options->journey_departure_type[j],
											true,
											Universe, 
											boundary1_state, 
											current_state+3, 
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
				this->V_infinity_in(k) = current_state[k+3] - boundary1_state[k+3];

			//store the left boundary state derivative if we are using time derivatives
			if (options->derivative_type > 2 && needG)
			{
				left_boundary_state_derivative[0] = boundary1_state[3];
				left_boundary_state_derivative[1] = boundary1_state[4];
				left_boundary_state_derivative[2] = boundary1_state[5];
				left_boundary_state_derivative[3] = boundary1_state[6];
				left_boundary_state_derivative[4] = boundary1_state[7];
				left_boundary_state_derivative[5] = boundary1_state[8];
			}

			vinf_in = this->V_infinity_in.norm();

			//step 3.2 compute outgoing v_infinity at flyby
			//this is drawn from the decision vector unless we are in the first phase of a journey that starts with a fixed v-infinity-out flyby
			if (p == 0 && options->journey_departure_type[j] == 4)
			{
				for (int k = 0; k < 3; ++k)
				{
					this->V_infinity_out(k) = options->journey_initial_velocity[j][k];
				}
			}
			else
			{
				for (int k = 0; k < 3; ++k)
				{
					this->V_infinity_out(k) = X[*Xindex];
					++(*Xindex);
				}
			}
			this->flyby_outgoing_v_infinity = this->V_infinity_out.norm();
			this->C3_departure = this->flyby_outgoing_v_infinity*this->flyby_outgoing_v_infinity;

			//step 3.3 enforce flyby velocity magnitude match constraint
			//Step 3.3.1 constraint value
			F[*Findex] = (this->flyby_outgoing_v_infinity - vinf_in) / Universe->LU * Universe->TU;
		
			//Step 3.3.2 constraint derivatives
			if (options->derivative_type > 0 && needG && !(p == 0 && options->journey_departure_type[j] == 4))
			{
				G[flyby_velocity_magnitude_constraint_G_indices[0]] = flyby_constraints_X_scale_ranges[0] * X[flyby_constraints_X_indices[0]] / this->flyby_outgoing_v_infinity / Universe->LU * Universe->TU;
				G[flyby_velocity_magnitude_constraint_G_indices[1]] = flyby_constraints_X_scale_ranges[1] * X[flyby_constraints_X_indices[1]] / this->flyby_outgoing_v_infinity / Universe->LU * Universe->TU;
				G[flyby_velocity_magnitude_constraint_G_indices[2]] = flyby_constraints_X_scale_ranges[2] * X[flyby_constraints_X_indices[2]] / this->flyby_outgoing_v_infinity / Universe->LU * Universe->TU;
				(*Gindex) += 3;

				bool compute_flyby_terminal_velocity_derivative = true;
				if (j > 0)
				{
					if (options->journey_departure_type[j] == 3 && !(options->journey_arrival_type[j-1] == 2))
						compute_flyby_terminal_velocity_derivative = false;
				}

				if (compute_flyby_terminal_velocity_derivative)
				{
					G[flyby_velocity_magnitude_constraint_G_indices[3]] = -flyby_constraints_X_scale_ranges[3] * X[flyby_constraints_X_indices[3]] / vinf_in / Universe->LU * Universe->TU;
					G[flyby_velocity_magnitude_constraint_G_indices[4]] = -flyby_constraints_X_scale_ranges[4] * X[flyby_constraints_X_indices[4]] / vinf_in / Universe->LU * Universe->TU;
					G[flyby_velocity_magnitude_constraint_G_indices[5]] = -flyby_constraints_X_scale_ranges[5] * X[flyby_constraints_X_indices[5]] / vinf_in / Universe->LU * Universe->TU;
					(*Gindex) += 3;
				}
			}
			++(*Findex);

			//step 3.4 enforce minimum flyby altitude constraint
			//step 3.4.1 compute the flyby turn angle
			this->flyby_turn_angle = acos(this->V_infinity_in.dot(this->V_infinity_out) / (vinf_in * this->flyby_outgoing_v_infinity));

			//step 3.4.2 compute the required periapse altitude
			if (this->flyby_turn_angle == 0.0 || this->C3_departure < sqrt(2*this->Body1->mu/(0.9*0.9*this->Body1->r_SOI)))
				this->flyby_altitude = 0.0;
			else
				this->flyby_altitude = this->Body1->mu/(this->C3_departure) * (1/sin(this->flyby_turn_angle/2.0) - 1) - this->Body1->radius;
		
			//step 3.4.3 constraint magnitude
			F[*Findex] = (this->Body1->minimum_safe_flyby_altitude - this->flyby_altitude) / (this->Body1->minimum_safe_flyby_altitude + this->Body1->radius);
		
			//Step 3.4.4 flyby altitude constraint derivatives
			if (options->derivative_type > 0 && needG && !(p == 0 && options->journey_departure_type[j] == 4))
			{
				double Vfx = this->V_infinity_out(0);
				double Vfy = this->V_infinity_out(1);
				double Vfz = this->V_infinity_out(2);
				double V0x = this->V_infinity_in(0);
				double V0y = this->V_infinity_in(1);
				double V0z = this->V_infinity_in(2);

				double VfdotVf = this->V_infinity_out.dot(V_infinity_out);
				double VfdotV0 = this->V_infinity_in.dot(V_infinity_out);
				double V0dotV0 = this->V_infinity_in.dot(V_infinity_in);
	
				double term1 = VfdotV0 / sqrt(VfdotVf * V0dotV0);
				double acos_term1_over2 = 0.5 * acos(term1);

				double mu = Body1->mu;

				double termAcoeff = 2.0 * mu * (1 / sin(acos_term1_over2) - 1) / (VfdotVf*VfdotVf);
				double termBcoeff = mu * cos(acos_term1_over2) / ((term1 - 1) * sqrt(1 - VfdotV0*VfdotV0/(V0dotV0*VfdotVf)) * sqrt(V0dotV0 * VfdotVf*VfdotVf*VfdotVf*VfdotVf*VfdotVf));
				double termA, termB;

				//Vfx
				termA = termAcoeff * Vfx;
				termB = termBcoeff * (V0x*Vfy*Vfy - V0y*Vfx*Vfy + V0x*Vfz*Vfz  - V0z*Vfx*Vfz);
	
				G[flyby_altitude_constraint_G_indices[0]] = flyby_constraints_X_scale_ranges[0] * (termA + termB) / (this->Body1->minimum_safe_flyby_altitude + this->Body1->radius);
			
				//Vfy
				termA = termAcoeff * Vfy;
				termB = termBcoeff * (V0y*Vfx*Vfx - V0x*Vfy*Vfx + V0y*Vfz*Vfz  - V0z*Vfy*Vfz);
			
				G[flyby_altitude_constraint_G_indices[1]] = flyby_constraints_X_scale_ranges[1] * (termA + termB) / (this->Body1->minimum_safe_flyby_altitude + this->Body1->radius);
			
				//Vfz
				termA = termAcoeff * Vfz;
				termB = termBcoeff * (V0z*Vfx*Vfx - V0x*Vfz*Vfx + V0z*Vfy*Vfy  - V0y*Vfz*Vfy);
			
				G[flyby_altitude_constraint_G_indices[2]] = flyby_constraints_X_scale_ranges[2] * (termA + termB) / (this->Body1->minimum_safe_flyby_altitude + this->Body1->radius);
			
				(*Gindex) += 3;

				bool compute_flyby_terminal_velocity_derivative = true;
				if (j > 0)
				{
					if (options->journey_departure_type[j] == 3 && !(options->journey_arrival_type[j-1] == 2))
						compute_flyby_terminal_velocity_derivative = false;
				}

				if (compute_flyby_terminal_velocity_derivative)
				{
					//derivatives with respect to V0
					double term;
					double Coeff = mu * cos(acos_term1_over2) / ((term1 - 1) * sqrt((1 - term1*term1) * V0dotV0*V0dotV0*V0dotV0 * VfdotVf*VfdotVf*VfdotVf));
			
					//V0x
					term = Coeff * (Vfx*V0y*V0y - V0x*Vfy*V0y + Vfx*V0z*V0z - V0x*Vfz*V0z);
					G[flyby_altitude_constraint_G_indices[5]] = flyby_constraints_X_scale_ranges[5] * term / (this->Body1->minimum_safe_flyby_altitude + this->Body1->radius);
			
					//V0y
					term = Coeff * (Vfy*V0x*V0x - V0y*Vfx*V0x + Vfy*V0z*V0z - V0y*Vfz*V0z);
					G[flyby_altitude_constraint_G_indices[4]] = flyby_constraints_X_scale_ranges[4] * term / (this->Body1->minimum_safe_flyby_altitude + this->Body1->radius);

					//V0z
					term = Coeff * (Vfz*V0x*V0x - V0z*Vfx*V0x + Vfz*V0y*V0y - V0z*Vfy*V0y);
					G[flyby_altitude_constraint_G_indices[3]] = flyby_constraints_X_scale_ranges[3] * term / (this->Body1->minimum_safe_flyby_altitude + this->Body1->radius);

					(*Gindex) += 3;
				}
			}

			++(*Findex);

			//calculate the b-plane parameters, check the periapse altitude
			this->BoundaryR.assign_all(boundary1_state);
			this->BoundaryV.assign_all(boundary1_state+3);

		

			//store the state at the beginning of the phase, post-flyby
			for (int k = 0; k < 3; ++k)
			{
				this->state_at_beginning_of_phase[k] = boundary1_state[k];
				this->state_at_beginning_of_phase[k+3] = boundary1_state[k+3] + V_infinity_out(k);
			}

			//store the mass
			if (p == 0)
				this->state_at_beginning_of_phase[6] = current_state[6] + options->journey_starting_mass_increment[j];
			else
				this->state_at_beginning_of_phase[6] = current_state[6];
		}// end flyby code

		//store the current epoch
		this->phase_start_epoch = *current_epoch;
	}

	void phase::process_right_boundary_condition(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, const int& needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
	{		
		//Step 5.1: if EMTG is choosing an input power or Isp for the phase (for REP/NEP models), then this information must be encoded
		if (!(options->mission_type == 0 || options->mission_type == 1 || options->mission_type == 4))
		{
			if (options->engine_type == 1)
			{
				//constant Isp, efficiency, EMTG computes input power
				if (j == 0 && p == 0)
					options->power_at_1_AU = X[*Xindex];

				for (int step = 0; step < options->num_timesteps; ++step)
					available_power[step] = options->power_at_1_AU;
				++(*Xindex);
			}
			else if (options->engine_type == 2)
			{
				//constant power, EMTG chooses Isp
				if (j == 0 && p == 0)
					options->IspLT = X[*Xindex];

				for (int step = 0; step < options->num_timesteps; ++step)
					available_Isp[step] = options->IspLT;
				++(*Xindex);
			}
			else if (options->objective_type == 13 && j == 0 && p == 0)
			{
				//constant Isp, efficiency, EMTG computes input power
				options->power_at_1_AU = X[*Xindex];
				++(*Xindex);
			}
		}

		//Step 5.2: extract the time of flight
		TOF = X[*Xindex];
		++(*Xindex);
		phase_start_epoch = *current_epoch;
		phase_end_epoch = *current_epoch + TOF;

		//Step 5.3: locate the second body
		locate_boundary_point(boundary2_location_code, options->journey_arrival_type[j], false, Universe, boundary2_state, current_state+3, *current_epoch + TOF, X, Xindex, F, Findex, G, Gindex, needG, j, p, options);
	
		//Step 5.4: if this is not a terminal rendezvous, extract the terminal velocity increment
		//otherwise, the terminal state is the body state or the terminal v-infinity
		if (!(p == options->number_of_phases[j] - 1 && (options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 5 || options->journey_arrival_type[j] == 6 || options->journey_arrival_type[j] == 7)))
		{
			dVarrival[0] = X[*Xindex];
			dVarrival[1] = X[*Xindex+1];
			dVarrival[2] = X[*Xindex+2];
			(*Xindex) += 3;

			for (int k = 0; k < 3; ++k)
			{
				state_at_end_of_phase[k] = boundary2_state[k];
				state_at_end_of_phase[k+3] = boundary2_state[k+3] + dVarrival[k];
			}
		}
		else if (options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 6 || options->journey_arrival_type[j] == 7)//for a rendezvous, escape, or capture spiral there is no terminal velocity increment
		{
			for (int k = 0; k < 6; ++k)
				state_at_end_of_phase[k] = boundary2_state[k];
		}
		else if (options->journey_arrival_type[j] == 5)//for matching a terminal velocity vector
		{
			for (int k = 0; k < 3; ++k)
			{
				state_at_end_of_phase[k] = boundary2_state[k];
				state_at_end_of_phase[k+3] = boundary2_state[k+3] + options->journey_final_velocity[j][k];
			}
		}


		//step 5.5 extract the mass at the end of the phase
		state_at_end_of_phase[6] = X[*Xindex];
		++(*Xindex);


		if (options->journey_arrival_type[j] == 7)//for journeys which end in a spiral
		{
			//journeys which end with a spiral have no terminal velocity vector
			this->C3_arrival = 0.0;
			this->RA_arrival = 0.0;
			this->DEC_arrival = 0.0;

			//*******************************************************
			//Step 4: compute the state pre-arrival

			//Step 4.1: compute the necessary inputs to the spiral model

			//Step 4.1.1 if we are using a VSI thruster then we need to extract the spiral Isp from the decision vector 
			if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
			{
				this->spiral_capture_Isp = X[*Xindex];
				++(*Xindex);
			}

			//Step 4.1.2 the state before the spiral is equal to the boundary state
			for (int k = 0; k < 7; ++k)
				this->spiral_capture_state_before_spiral[k] = state_at_end_of_phase[k];

				

			//Step 4.1.4 compute the starting body's distance from the Sun at the beginning of the spiral
			double distance_from_sun_in_AU;
			if (!(Universe->central_body_SPICE_ID == 10))
			{
				double central_body_state_in_km[6];
				double position_relative_to_sun_in_AU[3];
				Universe->locate_central_body(*current_epoch + this->TOF, central_body_state_in_km, options);

				position_relative_to_sun_in_AU[0] = (central_body_state_in_km[0] + this->spiral_capture_state_before_spiral[0]) / options->AU;
				position_relative_to_sun_in_AU[1] = (central_body_state_in_km[1] + this->spiral_capture_state_before_spiral[1]) / options->AU;
				position_relative_to_sun_in_AU[2] = (central_body_state_in_km[2] + this->spiral_capture_state_before_spiral[2]) / options->AU;

				distance_from_sun_in_AU = math::norm(position_relative_to_sun_in_AU, 3);
			}
			else //if we are orbiting the sun, it's quite silly to look up the position of the sun relative to itself, so don't bother
				distance_from_sun_in_AU = math::norm(boundary2_state, 3)/Universe->LU;

			//Step 4.2: process the spiral
			if (options->spiral_model_type == 0)
			{
				Astrodynamics::Battin_spiral(this->spiral_capture_state_before_spiral[6],
											options->journey_capture_spiral_final_radius[j],
											distance_from_sun_in_AU,
											*current_epoch - X[0] + this->TOF,
											this->spiral_capture_Isp,
											this->spiral_capture_thrust,
											this->spiral_capture_mdot,
											this->spiral_capture_number_of_engines,
											this->spiral_capture_power,
											this->spiral_capture_active_power,
											this->Body2->mu,
											(options->derivative_type > 1 && needG),
											options,
											&this->spiral_capture_time,
											&this->spiral_capture_mass_after,
											&this->spiral_capture_dv,
											&this->spiral_capture_dm_after_dm_before,
											&this->spiral_capture_dt_spiral_dm_before,
											&this->spiral_capture_dm_after_dIsp,
											&this->spiral_capture_dt_spiral_dIsp);
			}
			else if (options->spiral_model_type == 1)
			{
				Astrodynamics::Edelbaum_spiral(this->spiral_capture_state_before_spiral[6],
												options->journey_capture_spiral_final_radius[j],
												distance_from_sun_in_AU,
												Body2->r_SOI,
												*current_epoch - X[0] + this->TOF,
												this->spiral_capture_Isp,
												this->spiral_capture_thrust,
												this->spiral_capture_mdot,
												this->spiral_capture_number_of_engines,
												this->spiral_capture_power,
												this->spiral_capture_active_power,
												this->Body2->mu,
												(options->derivative_type > 1 && needG),
												options,
												&this->spiral_capture_time,
												&this->spiral_capture_mass_after,
												&this->spiral_capture_dv,
												&this->spiral_capture_dm_after_dm_before,
												&this->spiral_capture_dt_spiral_dm_before,
												&this->spiral_capture_dm_after_dIsp,
												&this->spiral_capture_dt_spiral_dIsp);
			}

			//Step 4.3 construct the post-spiral state
			//Step 4.3.1 advance time (not required)

			//Step 4.3.2 find the position of the body at the new phase start time and store it in spiral_capture_state_after_spiral
			locate_boundary_point(boundary2_location_code, options->journey_departure_type[j], true, Universe, this->spiral_capture_state_after_spiral, current_state+3, *current_epoch + this->TOF + this->spiral_capture_time, X, Xindex, F, Findex, G, Gindex, needG, j, p, options);

			//Step 4.3.3 fill in the new mass
			this->spiral_capture_state_after_spiral[6] = this->spiral_capture_mass_after;
		}//end capture spirals

	}

	double phase::process_arrival(	double* incoming_velocity,
									double* boundary_state,
									double* X_infinity,
									double* current_epoch,
									double mu,
									double r_SOI, 
									double* F,
									int* Findex, 
									int j, 
									missionoptions* options, 
									Astrodynamics::universe* Universe)
	{
		switch (options->journey_arrival_type[j])
			{
			case 0: //orbit insertion
				{
					//make sure that the user is not trying to use the simple orbit insertion with an actual fixed orbit
					if (boundary2_location_code == -2 || boundary2_location_code == -1)
					{
						cout << "When ending at a fixed point or orbit, must encode a 'rendezvous' in the options file. Cannot use simplified orbit insertion when the actual final orbit is modeled." << endl;
						throw 20;
					}

					double vinf;
					double dV = EMTG::Astrodynamics::insertion_burn(incoming_velocity, boundary_state+3, mu, r_SOI, options->journey_arrival_elements[j][0], options->journey_arrival_elements[j][1], &vinf);

					this->C3_arrival = vinf*vinf;
					this->dVarrival[0] = incoming_velocity[0] - boundary_state[3];
					this->dVarrival[1] = incoming_velocity[1] - boundary_state[4];
					this->dVarrival[2] = incoming_velocity[2] - boundary_state[5];

					if (options->journey_arrival_declination_constraint_flag[j])
					{
						//first we need to rotate the incoming velocity vector into the reference frame of the arrival body
						this->Body2->J2000_body_equatorial_frame.construct_rotation_matrices(*current_epoch / 86400.0 + 2400000.5);
						math::Matrix<double> rot_in_vec(3, 1, this->dVarrival);
						math::Matrix<double> rot_out_vec = this->Body2->J2000_body_equatorial_frame.R_from_ICRF_to_local * rot_in_vec;

						this->RA_arrival = atan2(rot_out_vec(1), rot_out_vec(0));

						this->DEC_arrival = asin(rot_out_vec(2) / vinf);

						F[*Findex] = this->DEC_arrival / options->journey_arrival_declination_bounds[j][1] - 1.0;
						++(*Findex);
					}

					return dV;
				}
			case 1: //rendezvous
				{
					for (int k = 0; k < 3; ++k)
						this->dVarrival[k] = boundary_state[k+3] - incoming_velocity[k];

					double vinf = math::norm(dVarrival, 3);
					this->C3_arrival = vinf*vinf;

					return vinf;
				}
			case 2: //flyby with bounded v-infinity
				{
					for (int k = 0; k < 3; ++k)
						this->dVarrival[k] = incoming_velocity[k] - boundary_state[k + 3];
					this->C3_arrival = dVarrival[0]*dVarrival[0] + dVarrival[1]*dVarrival[1] + dVarrival[2]*dVarrival[2];
					double C3_max = options->journey_final_velocity[j][1] * options->journey_final_velocity[j][1];
					F[*Findex] = C3_arrival / C3_max - 1;
					++(*Findex);

					if (options->journey_arrival_declination_constraint_flag[j])
					{
						//first we need to rotate the incoming velocity vector into the reference frame of the arrival body
						this->Body2->J2000_body_equatorial_frame.construct_rotation_matrices(*current_epoch / 86400.0 + 2400000.5);
						math::Matrix<double> rot_in_vec(3, 1, this->dVarrival);
						math::Matrix<double> rot_out_vec = this->Body2->J2000_body_equatorial_frame.R_from_ICRF_to_local * rot_in_vec;

						this->RA_arrival = atan2(rot_out_vec(1), rot_out_vec(0));

						this->DEC_arrival = asin(rot_out_vec(2) / sqrt(this->C3_arrival));

						F[*Findex] = this->DEC_arrival / options->journey_arrival_declination_bounds[j][1] - 1.0;
						++(*Findex);
					}

					return 0;
				}
			case 3: //low-thrust rendezvous
				{
					//the low-thrust rendezvous constraints are implicity satisfied because no terminal velocity vector is added to the target's velocity to find the spacecraft velocity
					this->C3_arrival = 0;
					this->DEC_arrival = 0.0;
					this->RA_arrival = 0.0;
					return 0;
				}
			case 4: //match v-infinity vector (impulsive)
				{
					for (int k = 0; k < 3; ++k)
						this->dVarrival[k] = (boundary_state[k+3] + options->journey_final_velocity[j][k]) - incoming_velocity[k];

					double vinf = math::norm(dVarrival, 3);
					this->C3_arrival = vinf*vinf;
				
					return vinf;
				}
			case 5: //match v-infinity vector (low-thrust)
				{
					//the low-thrust rendezvous constraints are implicity satisfied because the terminal velocity vector is fixed

					for (int k = 0; k < 3; ++k)
						dVarrival[k] = boundary_state[k+3] - incoming_velocity[k];
					this->C3_arrival = dVarrival[0]*dVarrival[0] + dVarrival[1]*dVarrival[1] + dVarrival[2]*dVarrival[2];

					return 0;
				}
			case 6: //escape (enforce zero energy constraint)
				{
					double v = math::norm(incoming_velocity, 3) * Universe->TU / Universe->LU;
					double r = math::norm(boundary_state, 3) / Universe->LU;

					//zero energy constraint
					F[*Findex] = v*v / 2.0 - 1.0 / r;
					++(*Findex);

					//derivatives of zero energy constraint

					//dVarrival
					this->dVarrival[0] = 0;
					this->dVarrival[1] = 0;
					this->dVarrival[2] = 0;

					return 0;
				}
			case 7: //low-thrust capture spiral (after low-thrust rendezvous) rendezvous
				{
					//the low-thrust rendezvous constraints are implicity satisfied because no terminal velocity vector is added to the target's velocity to find the spacecraft velocity
					return 0;
				}
			}
		std::cout << "Invalid arrival type" << endl;
		return 1e+10;
	}

	//function to locate boundary points
	int phase::locate_boundary_point(int location, int boundary_type, bool left_boundary, EMTG::Astrodynamics::universe* Universe, double* boundary_state, double* V_infinity, double epoch, double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, const int& needG, int j, int p,  missionoptions* options)
	{
		if (left_boundary) //this is the left boundary of the phase
		{
			if (location > 0) //if this boundary point is at a body
			{
				if (needG && options->derivative_type > 2)
				{
					Universe->bodies[location-1].locate_body(epoch,
															boundary_state,
															true,
															options);

					left_boundary_state_derivative[0] = boundary_state[3];
					left_boundary_state_derivative[1] = boundary_state[4];
					left_boundary_state_derivative[2] = boundary_state[5];
					left_boundary_state_derivative[3] = boundary_state[6];
					left_boundary_state_derivative[4] = boundary_state[7];
					left_boundary_state_derivative[5] = boundary_state[8];
				}
				else
					Universe->bodies[location-1].locate_body(epoch,
															boundary_state,
															false,
															options);
			}
			else if (location == -1 && j > 0)
			{
				//if we are starting this journey at a non-body point AND this is not the first journey then the starting location is wherever we ended the last journey
				//here we take advantage of the fact that the "V_infinity" field is actually the 4th, 5th, and 6th entries of a 7-vector containing the spacecraft state at the end of the previous journey
				double boundary_state_copy[6];
				for (int k = 0; k < 6; ++k)
					boundary_state_copy[k] = V_infinity[k - 3];

				double Fdummy, Gdummy, Ftdummy, Gtdummy;
				
				Kepler::KeplerLagrangeLaguerreConway(boundary_state_copy,
													boundary_state, 
													Universe->mu,
													this->phase_wait_time,
													&Fdummy, 
													&Gdummy,
													&Ftdummy,
													&Gtdummy);
			}
			else if (location == -1 && j == 0) //if this boundary point is at a free point in space, with the various elements either fixed or free
			{
				//For each element, either extract from the options structure or from the decision vector, depending on whether the options
				//structure specifies that the element should be varied

				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_departure_elements_vary_flag[j][k])
					{
						//if the user selected to vary this orbit element, extract it from the decision vector
						if (options->journey_departure_elements_type[j])
							this->left_boundary_orbit_elements[k] = X[*Xindex];
						else
							this->left_boundary_local_frame_state[k] = X[*Xindex];
						++(*Xindex);
					}
					else
					{
						//if the user selected to specify this orbit element
						if (options->journey_departure_elements_type[j])
							this->left_boundary_orbit_elements[k] = options->journey_departure_elements[j][k];
						else
							this->left_boundary_local_frame_state[k] = options->journey_departure_elements[j][k];
					}
				}

				//if the orbit was specified in COE, it needs to be converted to inertial coordinates
				if (options->journey_departure_elements_type[j])
					Astrodynamics::COE2inertial(left_boundary_orbit_elements, Universe->mu, this->left_boundary_local_frame_state);

				//if it is possible for the optimizer to select a point inside the exclusion zone of the central body
				//then there must be a nonlinear constraint to prevent this
				if (options->journey_departure_elements_type[j]) //classical orbit elements
				{
					//this constraint is applied if we are varying SMA or ECC
					if (options->journey_departure_elements_vary_flag[j][0] || options->journey_departure_elements_vary_flag[j][1])
					{
						double a = this->left_boundary_orbit_elements[0];
						double e = this->left_boundary_orbit_elements[1];
						double f = this->left_boundary_orbit_elements[5];
						double cosf = cos(f);
						double sinf = sqrt(1.0 - cosf*cosf);
						double r = a * (1 - e*e) / (1 + e * cosf);

						F[*Findex] = (r - Universe->minimum_safe_distance) / Universe->minimum_safe_distance;
						++(*Findex);

						if (needG)
						{
							//this constraint has derivatives with respect to SMA, ECC, and TA
							//only create a derivative entry with respect to an orbit element if that element is being varied
							int WhichDerivative = 0;
							if (options->journey_departure_elements_vary_flag[j][0]) //SMA
							{
								G[left_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance * (1 - e*e) / (1 + e * cosf);
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_departure_elements_vary_flag[j][1]) //ECC
							{
								G[left_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance * -(2*a*e + a*cosf * (1 + e*e)) / ((1 + e*cosf)*(1 + e*cosf));
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_departure_elements_vary_flag[j][5]) //TA
							{
								G[left_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance * a*e*sinf*(1 - e*e) / ((1 + e*cosf)*(1 + e*cosf));
								++(*Gindex);
								++WhichDerivative;
							}
						}
					}
				}
				else //cartesian orbit elements
				{
					//this constraint is applied if we are varying x, y, or z
					if (options->journey_departure_elements_vary_flag[j][0] || options->journey_departure_elements_vary_flag[j][1] || options->journey_departure_elements_vary_flag[j][2])
					{
						double x = this->left_boundary_local_frame_state[0];
						double y = this->left_boundary_local_frame_state[1];
						double z = this->left_boundary_local_frame_state[2];
						double r = sqrt(x*x + y*y + z*z);

						F[*Findex] = (r - Universe->minimum_safe_distance) / Universe->minimum_safe_distance;
						++(*Findex);

						if (needG)
						{
							//this constraint has derivatives with respect to x, y, and z
							//only create a derivative entry with respect to an orbit element if that element is being varied
							int WhichDerivative = 0;
							if (options->journey_departure_elements_vary_flag[j][0]) //x
							{
								G[left_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance / r * X[left_boundary_central_body_exclusion_radius_constraint_X_indices[WhichDerivative]];
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_departure_elements_vary_flag[j][1]) //y
							{
								G[left_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance / r * X[left_boundary_central_body_exclusion_radius_constraint_X_indices[WhichDerivative]];
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_departure_elements_vary_flag[j][2]) //z
							{
								G[left_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance / r * X[left_boundary_central_body_exclusion_radius_constraint_X_indices[WhichDerivative]];
								++(*Gindex);
								++WhichDerivative;
							}
						}
					}
				}

				//the orbit must be rotated into EMTG's internal frame, J2000 Earth Equatorial
				Universe->LocalFrame.construct_rotation_matrices(epoch / 86400.0 + 2400000.5);
				math::Matrix<double> rot_in_vec(3,1,this->left_boundary_local_frame_state);
				math::Matrix<double> rot_out_vec = Universe->LocalFrame.R_from_local_to_ICRF * rot_in_vec;
				for (int k = 0; k < 3; ++k)
				{
					boundary_state[k] = rot_out_vec(k);
					rot_in_vec(k) = this->left_boundary_local_frame_state[k+3];
				}
				rot_out_vec = Universe->LocalFrame.R_from_local_to_ICRF * rot_in_vec;
				for (int k = 0; k < 3; ++k)
				{
					boundary_state[k+3] = rot_out_vec(k);
				}

				//finally, we have the opportunity to propagate the state forward in time if there is an associated journey wait time
				double Fdummy, Gdummy, Ftdummy, Gtdummy;
				double boundary_copy[6];
				for (int k = 0; k < 6; ++k)
					boundary_copy[k] = boundary_state[k];

				Kepler::KeplerLagrangeLaguerreConway(boundary_copy,
													boundary_state,
													Universe->mu,
													this->phase_wait_time,
													&Fdummy,
													&Gdummy,
													&Ftdummy,
													&Gtdummy);
			}
		}
		else //this is the right boundary of the phase, so this is an arrival
		{

			if (location > 0) //if this boundary point is at a body
			{
				if (needG && options->derivative_type > 2)
				{
					Universe->bodies[location-1].locate_body(epoch, boundary_state, true, options);
					right_boundary_state_derivative[0] = boundary_state[3];
					right_boundary_state_derivative[1] = boundary_state[4];
					right_boundary_state_derivative[2] = boundary_state[5];
					right_boundary_state_derivative[3] = boundary_state[6];
					right_boundary_state_derivative[4] = boundary_state[7];
					right_boundary_state_derivative[5] = boundary_state[8];
				}
				else
					Universe->bodies[location-1].locate_body(epoch, boundary_state, false, options);
			}
			else if (location == -1) //if this boundary point is at a free point in space, with the various elements either fixed or free
			{
				//note: boundary point is supplied in the local Universe frame

				//For each element, either extract from the options structure or from the decision vector, depending on whether the options
				//structure specifies that the element should be varied
				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_arrival_elements_vary_flag[j][k])
					{
						//if the user selected to vary this orbit element, extract it from the decision vector
						if (options->journey_arrival_elements_type[j])
							this->right_boundary_orbit_elements[k] = X[*Xindex];
						else
							this->right_boundary_local_frame_state[k] = X[*Xindex];
						++(*Xindex);
					}
					else
					{
						//if the user selected to specify this orbit element
						if (options->journey_arrival_elements_type[j])
							this->right_boundary_orbit_elements[k] = options->journey_arrival_elements[j][k];
						else
							this->right_boundary_local_frame_state[k] = options->journey_arrival_elements[j][k];
					}
				}

				//if the orbit was specified in COE, it needs to be converted to inertial coordinates
				if (options->journey_arrival_elements_type[j])
					Astrodynamics::COE2inertial(this->right_boundary_orbit_elements, Universe->mu, this->right_boundary_local_frame_state);
				//if it is possible for the optimizer to select a point inside the exclusion zone of the central body
				//then there must be a nonlinear constraint to prevent this
				if (options->journey_arrival_elements_type[j]) //classical orbit elements
				{
					//this constraint is applied if we are varying SMA or ECC
					if (options->journey_arrival_elements_vary_flag[j][0] || options->journey_arrival_elements_vary_flag[j][1])
					{
						double a = this->right_boundary_orbit_elements[0];
						double e = this->right_boundary_orbit_elements[1];
						double f = this->right_boundary_orbit_elements[5];
						double cosf = cos(f);
						double sinf = sqrt(1.0 - cosf*cosf);
						double r = a * (1 - e*e) / (1 + e * cosf);

						F[*Findex] = (r - Universe->minimum_safe_distance) / Universe->minimum_safe_distance;
						++(*Findex);

						if (needG)
						{
							//this constraint has derivatives with respect to SMA, ECC, and TA
							//only create a derivative entry with respect to an orbit element if that element is being varied
							int WhichDerivative = 0;
							if (options->journey_arrival_elements_vary_flag[j][0]) //SMA
							{
								G[right_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance * (1 - e*e) / (1 + e * cosf);
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_arrival_elements_vary_flag[j][1]) //ECC
							{
								G[right_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance * -(2*a*e + a*cosf * (1 + e*e)) / ((1 + e*cosf)*(1 + e*cosf));
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_arrival_elements_vary_flag[j][5]) //TA
							{
								G[right_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance * a*e*sinf*(1 - e*e) / ((1 + e*cosf)*(1 + e*cosf));
								++(*Gindex);
								++WhichDerivative;
							}
						}
					}
				}
				else //cartesian orbit elements
				{
					//this constraint is applied if we are varying x, y, or z
					if (options->journey_arrival_elements_vary_flag[j][0] || options->journey_arrival_elements_vary_flag[j][1] || options->journey_arrival_elements_vary_flag[j][2])
					{
						double x = this->right_boundary_local_frame_state[0];
						double y = this->right_boundary_local_frame_state[1];
						double z = this->right_boundary_local_frame_state[2];
						double r = sqrt(x*x + y*y + z*z);

						F[*Findex] = (r - Universe->minimum_safe_distance) / Universe->minimum_safe_distance;
						++(*Findex);

						if (needG)
						{
							//this constraint has derivatives with respect to x, y, and z
							//only create a derivative entry with respect to an orbit element if that element is being varied
							int WhichDerivative = 0;
							if (options->journey_arrival_elements_vary_flag[j][0]) //x
							{
								G[right_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance / r * X[right_boundary_central_body_exclusion_radius_constraint_X_indices[WhichDerivative]];
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_arrival_elements_vary_flag[j][1]) //y
							{
								G[right_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance / r * X[right_boundary_central_body_exclusion_radius_constraint_X_indices[WhichDerivative]];
								++(*Gindex);
								++WhichDerivative;
							}
							if (options->journey_arrival_elements_vary_flag[j][2]) //z
							{
								G[right_boundary_central_body_exclusion_radius_constraint_G_indices[WhichDerivative]] = right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges[WhichDerivative] / Universe->minimum_safe_distance / r * X[right_boundary_central_body_exclusion_radius_constraint_X_indices[WhichDerivative]];
								++(*Gindex);
								++WhichDerivative;
							}
						}
					}
				}

				//finally, the orbit must be rotated into EMTG's internal frame, J2000 Earth Equatorial
				Universe->LocalFrame.construct_rotation_matrices(epoch / 86400.0 + 2400000.5);
				math::Matrix<double> rot_in_vec(3,1,this->right_boundary_local_frame_state);
				math::Matrix<double> rot_out_vec = Universe->LocalFrame.R_from_local_to_ICRF * rot_in_vec;
				for (int k = 0; k < 3; ++k)
				{
					boundary_state[k] = rot_out_vec(k);
					rot_in_vec(k) = this->right_boundary_local_frame_state[k+3];
				}
				rot_out_vec = Universe->LocalFrame.R_from_local_to_ICRF * rot_in_vec;
				for (int k = 0; k < 3; ++k)
				{
					boundary_state[k+3] = rot_out_vec(k);
				}

			}
		}

		return 0;
	}

	//function to calculate flyby periapse state given v_infinity)in and v_infinity_out
	//this is the basic version that assumes a symmetric, unpowered flyby
	//returns state vector in MJ2000Eq / ICRF
	//credit to Max Schadegg summer 2013
	math::Matrix<double> phase::calculate_flyby_periapse_state(math::Matrix<double>& Vinf_in, math::Matrix<double>& Vinf_out, const double& flyby_altitude, EMTG::Astrodynamics::body& TheBody)
	{
        //calculate unit vector pointed towards periapse
        math::Matrix<double> periapse_position_unit_vector = (Vinf_in - Vinf_out).unitize();

        //calculate angular momentum unit vector
        math::Matrix<double> angular_momentum_vector = Vinf_in.cross(Vinf_out);            

        //calculate velocity unit vector at periapse
        math::Matrix<double> periapse_velocity_unit_vector = (angular_momentum_vector.cross(periapse_position_unit_vector)).unitize();
    
        //calculate velocity magnitude at periapse
        double periapse_velocity_magnitude = sqrt(2 * TheBody.mu / (TheBody.radius + flyby_altitude) + Vinf_in.dot(Vinf_in));

		//transform from unit vector space to state space
		math::Matrix<double> periapse_position_vector = periapse_position_unit_vector * (TheBody.radius + flyby_altitude);
		math::Matrix<double> periapse_velocity_vector = periapse_velocity_unit_vector * periapse_velocity_magnitude;

		return periapse_position_vector.vert_cat(periapse_velocity_vector);
	}

	math::Matrix<double> phase::calculate_periapse_state_from_asymptote_and_parking_orbit(math::Matrix<double>& V_infinity, double parking_orbit_incination, double parking_orbit_altitude, double &epoch, EMTG::Astrodynamics::universe* Universe, EMTG::Astrodynamics::body &TheBody)
	{
		//adapted from David Eagle
		//http://www.cdeagle.com/interplanetary/hyper1_matlab.pdf
		//http://www.cdeagle.com/interplanetary/hyper_ftn.pdf


		//construct the unit vector in the direction of the velocity asymptote
		math::Matrix<double> Shat = V_infinity.unitize();

		//construct the unit vector in the direction of the planet's pole, expressed in Earth J2000Eq
		TheBody.J2000_body_equatorial_frame.construct_rotation_matrices(epoch / 86400.0 + 2400000.5);
		double zhatl[] = {0,0,1};
		math::Matrix<double> zhat_local(3,1, zhatl);

		math::Matrix<double> Phat = TheBody.J2000_body_equatorial_frame.R_from_local_to_ICRF * zhat_local;

		//compute the angle between the outgoing asymptote and the spin axis of the planet
		double Beta = acos(Shat.dot(Phat));

		//compute declination of launch asymptote (DLA)
		double DLA = math::PIover2 - Beta;

		//compute right ascension of asymptote
		double RLA = atan2(V_infinity(1), V_infinity(0));

		//if absolute value of the DLA is larger than the inclination of the parking orbit then this function will not work
		if (fabs(DLA) > fabs(parking_orbit_incination))
		{
			cout << "fabs(DLA) > fabs(parking_orbit_incination, cannot generate periapse state using existing code" << endl;
			throw 1711;
		}

		//compute the magnitude of the velocity asymptote
		double v_infinity = V_infinity.norm();

		//compute the distance from the planet at periapse
		double r_periapse = TheBody.radius + parking_orbit_altitude;

		//compute the orbit elements of the parking orbit
		//inclination was an input
		double INC = parking_orbit_incination;

		//compute candidate RAAN
		double RAAN1 = math::PIover2 + RLA + asin(1.0/ (tan(Beta) * tan(INC)));
		double RAAN2 = math::PI + RLA - asin(1.0/ (tan(Beta) * tan(INC)));

		//compute candidate TA
		double Eta = asin( 1.0 / ( 1.0 + r_periapse * v_infinity*v_infinity / TheBody.mu ) );
		double TA1 = acos( cos(Beta) / sin(INC) ) - Eta;
		double TA2 = -acos( cos(Beta) / sin(INC) ) - Eta;

		//there are two equal solutions for RAAN and TA, but we will only use the first one

		//compute the unit vector in the direction of the spacecraft's position at the moment of injection
		math::Matrix<double> Rhat(3,1);
		Rhat(0) = cos(RAAN1) * cos(TA1) - sin(RAAN1) * cos(INC) * sin(TA1);
		Rhat(1) = sin(RAAN1) * cos(TA1) + cos(RAAN1) * cos(INC) * sin(TA1);
		Rhat(2) = sin(INC) * sin(TA1);

		//compute the position vector at the moment of injection
		math::Matrix<double> R = Rhat * r_periapse;

		//compute the cosine of the angle between the spacecraft's position vector and the departure asymptote unit vector
		double cos_psi = Shat.dot(Rhat);

		//compute the velocity vector at periapse required to achieve the desired hyperbolic asymptote
		double d = sqrt(TheBody.mu / ( ( 1 + cos_psi ) * r_periapse ) + v_infinity*v_infinity / 4);

		math::Matrix<double> V = Shat * (d + v_infinity / 2.0) + Rhat * (d - v_infinity / 2.0);

		return R.vert_cat(V);
	}

	double phase::compute_timestep_width_from_distribution(double step, missionoptions* options, double& scale_or_stdv)
	{
		switch (options->step_size_distribution)
		{
		case 0: //uniform distribution
			return 1.0 / options->num_timesteps;
		/*case 1: case 3://Gaussian distribution with or without variable standard deviation
			{
				double t = (double) step / options->num_timesteps;
				return 1.0 / (scale_or_stdv * sqrt(math::TwoPI)) * exp(-(t - 0.5)*(t - 0.5)/(2*scale_or_stdv*scale_or_stdv));
			}
		case 2: case 4://Cauchy distribution with or without variable scale width
			{
				double t = (double) step / options->num_timesteps;
				return scale_or_stdv / ( math::PI * ( scale_or_stdv*scale_or_stdv + (t - 0.5)*(t - 0.5) ) );
			}*/
		}

		return 0.0;
	}

	void phase::write_summary_line(missionoptions* options,
									EMTG::Astrodynamics::universe* Universe,
									int* eventcount,
									double current_epoch_MJD,
									string event_type,
									string boundary_name,
									double timestep_size,
									double flyby_altitude,
									double bdotr,
									double bdott,
									double angle1,
									double angle2,
									double C3,
									double* state,
									double* dV,
									double* ThrustVector,
									double dVmag,
									double Thrust,
									double Isp,
									double AvailPower,
									double mdot,
									int number_of_active_engines,
									double active_power)
	{
		//create a 3-element storage vector that will be used every time something needs to be rotated to the local frame
		math::Matrix<double> display_vector(3,1);
		math::Matrix<double> rot_in_vec(3,1);

		//construct the rotation matrix from ICRF to the local frame
		Universe->LocalFrame.construct_rotation_matrices(current_epoch_MJD + 2400000.5);

		//open the output file
		ofstream outputfile(options->outputfile.c_str(), ios::out | ios::app);

		outputfile.width(5); outputfile << *eventcount;
		++(*eventcount);
		outputfile.width(3); outputfile << " | ";

		//output the event epoch in both MJD and MM/DD/YYYY
		outputfile.width(16); outputfile.setf(ios::fixed, ios::floatfield); outputfile.precision(8); outputfile << current_epoch_MJD + 2400000.5;
		outputfile.width(3); outputfile << " | ";
		int month, day, year, hrs, mins;
		double secs;
		mjd_to_mdyhms(current_epoch_MJD, &month, &day, &year, &hrs, &mins, &secs);
		stringstream datestream;
		datestream << month << "/" << day << "/" << year;
		outputfile.width(11); outputfile << datestream.str();
		outputfile.width(3); outputfile << " | ";

		outputfile.width(12); outputfile << event_type;
		outputfile.width(3); outputfile << " | ";

		outputfile.width(25); outputfile << boundary_name;
		outputfile.width(3); outputfile << " | ";

		outputfile.width(15); outputfile << timestep_size;
		outputfile.width(3); outputfile << " | ";

		//rp, BdotR, and BdotT
		if (event_type == "upwr_flyby" || event_type == "pwr_flyby")
		{
			outputfile.width(19); outputfile << flyby_altitude;
			outputfile.width(3); outputfile << " | ";
			outputfile.width(19); outputfile << BdotR;
			outputfile.width(3); outputfile << " | ";
			outputfile.width(19); outputfile << BdotT;
			outputfile.width(3); outputfile << " | ";
		}
		else if (boundary_name == "Hyp-arrival")
		{
			outputfile.width(19); outputfile << math::norm(state, 3);
			outputfile.width(3); outputfile << " | ";
			outputfile.width(19); outputfile << BdotR;
			outputfile.width(3); outputfile << " | ";
			outputfile.width(19); outputfile << BdotT;
			outputfile.width(3); outputfile << " | ";
		}
		else
		{
			outputfile.width(19); outputfile << "-";
			outputfile.width(3); outputfile << " | ";
			outputfile.width(19); outputfile << "-";
			outputfile.width(3); outputfile << " | ";
			outputfile.width(19); outputfile << "-";
			outputfile.width(3); outputfile << " | ";
		}

		//control angles
		if (!(event_type == "coast"))
		{
			outputfile.precision(3); outputfile.width(8); outputfile << angle1 * 180.0/math::PI;
			outputfile.width(3); outputfile << " | ";
			outputfile.precision(3); outputfile.width(8); outputfile << angle2 * 180.0/math::PI;
			outputfile.width(3); outputfile << " | ";
		}
		else
		{
			outputfile.precision(3); outputfile.width(8); outputfile << 0.0;
			outputfile.width(3); outputfile << " | ";
			outputfile.precision(3); outputfile.width(8); outputfile << 0.0;
			outputfile.width(3); outputfile << " | ";
		}

		//C3
		if (event_type == "upwr_flyby" || event_type == "pwr_flyby" || (event_type == "launch" && options->LV_type >= 0) || event_type == "intercept" || event_type == "insertion" || event_type == "departure" || event_type == "interface" || event_type == "zeroflyby")
		{
			outputfile.precision(5);
			if (boundary_name == "Hyp-arrival")
			{
				double v0 = math::norm(state+3,3);
				outputfile.width(14); outputfile << v0*v0;
			}
			else
			{
				outputfile.width(14); outputfile << C3;
			}
		}
		else
		{
			outputfile.width(14); outputfile << "-";
		}
		outputfile.width(3); outputfile << " | ";


		//state at event +
		if (options->output_units == 0)
		{
			outputfile.precision(8);
			for (size_t k = 0; k < 3; ++k)
				rot_in_vec(k) = state[k];
			display_vector = Universe->LocalFrame.R_from_ICRF_to_local * rot_in_vec;
			for (size_t k=0;k<3;++k)
			{
				outputfile.width(19); outputfile << display_vector(k);
				outputfile.width(3); outputfile << " | ";
			}
			outputfile.precision(8);
			for (size_t k = 0; k < 3; ++k)
				rot_in_vec(k) = state[k+3];
			display_vector = Universe->LocalFrame.R_from_ICRF_to_local * rot_in_vec;
			for (size_t k=0;k<3;++k)
			{
				outputfile.width(19); outputfile << display_vector(k);
				outputfile.width(3); outputfile << " | ";
			}
		}
		else if (options->output_units == 1) //LU
		{
			outputfile.precision(8);
			for (size_t k = 0; k < 3; ++k)
				rot_in_vec(k) = state[k];
			display_vector = Universe->LocalFrame.R_from_ICRF_to_local * rot_in_vec;
			for (size_t k=0;k<3;++k)
			{
				outputfile.width(19); outputfile << display_vector(k) / Universe->LU;
				outputfile.width(3); outputfile << " | ";
			}
			for (size_t k = 0; k < 3; ++k)
				rot_in_vec(k) = state[k+3];
			display_vector = Universe->LocalFrame.R_from_ICRF_to_local * rot_in_vec;
			for (size_t k=0;k<3;++k)
			{
				outputfile.width(19); outputfile << display_vector(k) / Universe->LU * 86400;
				outputfile.width(3); outputfile << " | ";
			}
		}

		//deltaV - we can't track deltaV vectors for flybys
		if (event_type == "upwr_flyby" || event_type == "pwr_flyby" || event_type == "LT_rndzvs" || event_type == "zeroflyby")
		{
			for (int k = 0; k < 3; ++k)
			{
				outputfile.width(19); outputfile <<  "-";
				outputfile.width(3); outputfile << " | ";
			}
		}
		else
		{
			outputfile.precision(12);
			for (size_t k = 0; k < 3; ++k)
				rot_in_vec(k) = dV[k];
			display_vector = Universe->LocalFrame.R_from_ICRF_to_local * rot_in_vec;
			if (options->output_units == 0) //dV in km/s
			{
				for (int k = 0; k < 3; ++k)
				{
					outputfile.width(19); outputfile << display_vector(k);
					outputfile.width(3); outputfile << " | ";
				}
			}
			else
			{
				for (size_t k=0;k<3;++k) //dV in LU/day
				{
					outputfile.width(19); outputfile << display_vector(k) / Universe->LU * 86400;
					outputfile.width(3); outputfile << " | ";
				}
			}
		}

		//Thrust
		outputfile.precision(12);
		if (event_type == "SFthrust" || event_type == "FBLTthrust" || event_type == "PSBIthrust")
		{
			for (size_t k = 0; k < 3; ++k)
				rot_in_vec(k) = ThrustVector[k];
			display_vector = Universe->LocalFrame.R_from_ICRF_to_local * rot_in_vec;
		
			for (int k = 0; k < 3; ++k)
			{
				outputfile.width(19); outputfile << display_vector(k);
				outputfile.width(3); outputfile << " | ";
			}
		}
		else
		{
			for (int k = 0; k < 3; ++k)
			{
				outputfile.width(19); outputfile <<  "-";
				outputfile.width(3); outputfile << " | ";
			}
		}

		//dV magnitude
		outputfile.precision(5);
		outputfile.width(17); outputfile << dVmag;
		outputfile.width(3); outputfile << " | ";

		//thrust, Isp, power
		outputfile.precision(5);
		if (event_type == "coast" || event_type == "force-coast" || event_type == "upwr_flyby" || event_type == "intercept" || event_type == "interface" || event_type == "LT_rndzvs" || event_type == "departure" || event_type == "match_point" || event_type == "match-vinf" || event_type == "zeroflyby")
		{
			outputfile.width(14); outputfile << "-";
			outputfile.width(3); outputfile << " | ";
		}
		else if (Thrust > 1.0e-6)
		{
			outputfile.width(14); outputfile << Thrust;
			outputfile.width(3); outputfile << " | ";
		}
		else
		{
			outputfile.width(14); outputfile << "impulse";
			outputfile.width(3); outputfile << " | ";
		}

		outputfile.precision(0);
		if (event_type == "coast" || event_type == "force-coast" || event_type == "upwr_flyby" || event_type == "intercept" || event_type == "interface" || event_type == "LT_rndzvs" || event_type == "departure" || event_type == "match_point" || event_type == "match-vinf" || event_type == "zeroflyby")
		{
			outputfile.width(14); outputfile << "-";
			outputfile.width(3); outputfile << " | ";
		}
		else if (Isp > 0)
		{
			outputfile.width(14); outputfile << Isp;
			outputfile.width(3); outputfile << " | ";
		} 
		else if (event_type == "launch")
		{
			outputfile.width(14); outputfile << "LV-supplied";
			outputfile.width(3); outputfile << " | ";
		}
		else
		{
			outputfile.width(14); outputfile << "UNHANDLED EVENT TYPE";
			outputfile.width(3); outputfile << " | ";
		}


		outputfile.precision(5);
		if (AvailPower > 0)
		{
			outputfile.width(14); outputfile << AvailPower;
			outputfile.width(3); outputfile << " | ";
		}
		else
		{
			outputfile.width(14); outputfile << "-";
			outputfile.width(3); outputfile << " | ";
		}

		if (event_type == "SFthrust" || event_type == "FBLTthrust" || event_type == "PSBIthrust" || event_type == "begin_spiral" || event_type == "end_spiral")
		{
			outputfile.precision(8);
			outputfile.width(19); outputfile << scientific << mdot << fixed;
		
		}
		else
		{
			outputfile.width(19); outputfile << "-";
		}
		outputfile.width(3); outputfile << " | ";

		//mass
		outputfile.width(14); outputfile.precision(4); outputfile << state[6];

		outputfile.width(3); outputfile << " | ";
		
		//number of active engines
		outputfile.width(14);
		if (event_type == "SFthrust" || event_type == "FBLTthrust" || event_type == "PSBIthrust" || event_type == "begin_spiral" || event_type == "end_spiral")
			outputfile << number_of_active_engines;
		else
			outputfile << "-";
		
		outputfile.width(3); outputfile << " | ";
		outputfile.precision(5);
		if (event_type == "SFthrust" || event_type == "FBLTthrust" || event_type == "PSBIthrust" || event_type == "begin_spiral" || event_type == "end_spiral")
		{
			outputfile.width(14); outputfile << active_power;
			outputfile.width(3); outputfile << " | ";
		}
		else
		{
			outputfile.width(14); outputfile << "-";
			outputfile.width(3); outputfile << " | ";
		}

		outputfile << endl;
	}

	

	void phase::calcbounds_left_boundary(const string& prefix, int first_X_entry_in_phase, vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, int j, int p,  EMTG::Astrodynamics::universe* Universe, missionoptions* options)
	{
		//if applicable, vary the journey initial mass increment
		if (options->journey_starting_mass_increment[j] > 0.0 && options->journey_variable_mass_increment[j] && p == 0)
		{
			Xlowerbounds->push_back(0.0);
			Xupperbounds->push_back(1.0);
			Xdescriptions->push_back(prefix + "journey initial mass scale factor");
		}

		//next, we need to know if we are the first phase in the journey and the journey does not start with a flyby
		if (p == 0 && !(options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 6))
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
				//If not, then the first decision variable is the stay time at the first body in the journey (i.e. at the asteroid for sample return)
				Xlowerbounds->push_back(options->journey_wait_time_bounds[j][0]);
				Xupperbounds->push_back(options->journey_wait_time_bounds[j][1]);
				Xdescriptions->push_back(prefix + "stay time (days)");

			}

			//if this boundary point is at a free point in space, with the various elements either fixed or free
			//this is only relevant for the first journey - succcessive journeys will start from the right hand boundary of the previous journey
			if  (boundary1_location_code == -1 && j == 0)
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
			}\

			//then we have up to four variables to parameterize the departure
			if ((j == 0 || !(options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4 || options->journey_departure_type[j] == 6)))
			{
				if (!(options->journey_departure_type[j] == 5 || options->journey_departure_type[j] == 2)) //for journeys that do not start from a spiral or a free direct departure
				{
					//First, we have outgoing velocity
					Xlowerbounds->push_back(options->journey_initial_impulse_bounds[j][0]);
					Xupperbounds->push_back(options->journey_initial_impulse_bounds[j][1]);

					//make sure that if we are launching on a rocket that the maximum C3 does not exceed that of the rocket
					if (j == 0 && options->journey_departure_type[0] == 0)
					{
						double C3max = 100.0;
						double C3max_LV[] = {35,50,40,60,60,60,60,60,60,60,60,60,200,0,100,0};
						if (options->LV_type > 0)
							C3max = C3max_LV[options->LV_type - 1];
						else if (options->LV_type == -2)
							C3max = options->custom_LV_C3_bounds[1];

						if (Xupperbounds->back() > sqrt(C3max))
							Xupperbounds->back() = sqrt(C3max);
					}
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
				else if (options->journey_departure_type[j] == 5)//for journeys which start from a spiral - encode the spiral Isp if appropriate
				{
					if (options->engine_type == 4 || options->engine_type == 13) //generic VSI or Xenon Hall thruster
					{
						Xlowerbounds->push_back(options->IspLT_minimum);
						Xupperbounds->push_back(options->IspLT);
						Xdescriptions->push_back(prefix + "Escape spiral Isp");
					}
					else if (options->engine_type == 12) //VASIMR
					{
						Xlowerbounds->push_back(3000.0);
						Xupperbounds->push_back(5000.0);
						Xdescriptions->push_back(prefix + "Escape spiral Isp");
					}
				}

				if (j == 0 && options->allow_initial_mass_to_vary)
				{
					//if we have enabled varying the initial mass, then pass through a mass multiplier
					Xlowerbounds->push_back(0.1);
					Xupperbounds->push_back(1.0);
					Xdescriptions->push_back(prefix + "initial mass multiplier (0-1)");
				}
			}
		}
		else
		{
			//if this phase starts with a flyby

			//we need to encode the initial velocity increment in inertial coordinates
			//special case: if this is the first phase of a journey that starts with a fixed v-infinity flyby, then we do not need to encode anything
			if (!(p == 0 && options->journey_departure_type[j] == 4))
			{
				Xlowerbounds->push_back(-25.0);
				Xupperbounds->push_back(25.0);
				Xdescriptions->push_back(prefix + "initial velocity increment x");
				Xlowerbounds->push_back(-25.0);
				Xupperbounds->push_back(25.0);
				Xdescriptions->push_back(prefix + "initial velocity increment y");
				Xlowerbounds->push_back(-25.0);
				Xupperbounds->push_back(25.0);
				Xdescriptions->push_back(prefix + "initial velocity increment z");
			}

			if (p == 0 && options->journey_departure_type[j] == 6) //constraints for a zero-turn flyby
			{
				Flowerbounds->push_back(-math::SMALL);
				Fupperbounds->push_back(math::SMALL);
				Fdescriptions->push_back(prefix + "Incoming and outcoming v_infinity_x must match at the flyby");

				//Jacobian entry (nonlinear) for the flyby X velocity match constraint
				//the velocity match constraint is dependent on this phase's initial velocity vector components and the previous phase's terminal velocity vector components
				//
				//Initial velocity:
				for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
				{
					if ((*Xdescriptions)[entry].find("initial velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby X velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_constraints_X_indices.push_back(entry);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
						break;
					}
				}
				//Terminal velocity:
				for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
				{
					if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry-2);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby X velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry-2 << "]: " << (*Xdescriptions)[entry-2];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_constraints_X_indices.push_back(entry-2);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry-2] - (*Xlowerbounds)[entry-2]);
						break;
					}
				}

				Flowerbounds->push_back(-math::SMALL);
				Fupperbounds->push_back(math::SMALL);
				Fdescriptions->push_back(prefix + "Incoming and outcoming v_infinity_y must match at the flyby");

				//Jacobian entry (nonlinear) for the flyby Y velocity match constraint
				//the velocity match constraint is dependent on this phase's initial velocity vector components and the previous phase's terminal velocity vector components
				//
				//Initial velocity:
				for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
				{
					if ((*Xdescriptions)[entry].find("initial velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry+1);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby Y velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry+1 << "]: " << (*Xdescriptions)[entry+1];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_constraints_X_indices.push_back(entry+1);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry+1] - (*Xlowerbounds)[entry+1]);
						break;
					}
				}
				//Terminal velocity:
				for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
				{
					if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry-1);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby Y velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry-1 << "]: " << (*Xdescriptions)[entry-1];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_constraints_X_indices.push_back(entry-1);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry-1] - (*Xlowerbounds)[entry-1]);
						break;
					}
				}

				Flowerbounds->push_back(-math::SMALL);
				Fupperbounds->push_back(math::SMALL);
				Fdescriptions->push_back(prefix + "Incoming and outcoming v_infinity_z must match at the flyby");

				//Jacobian entry (nonlinear) for the flyby Z velocity match constraint
				//the velocity match constraint is dependent on this phase's initial velocity vector components and the previous phase's terminal velocity vector components
				//
				//Initial velocity:
				for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
				{
					if ((*Xdescriptions)[entry].find("initial velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry+2);
						stringstream entryNameStream;
						entryNameStream << "Derivative of " << prefix << " flyby Z velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry+2 << "]: " << (*Xdescriptions)[entry+2];
						Gdescriptions->push_back(entryNameStream.str());
						flyby_constraints_X_indices.push_back(entry+2);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry+2] - (*Xlowerbounds)[entry+2]);
						break;
					}
				}
				//Terminal velocity:
				for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
				{
					if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby Z velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_constraints_X_indices.push_back(entry);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
						break;
					}
				}
			}
			else //constraints for other types of flyby
			{
				//we must encode a flyby velocity match constraint
				Flowerbounds->push_back(-math::SMALL);
				Fupperbounds->push_back(math::SMALL);
				Fdescriptions->push_back(prefix + "Incoming and outcoming v_infinity must match at the flyby");

				//Jacobian entry (nonlinear) for the flyby velocity match constraint
				//the velocity match constraint is dependent on this phase's initial velocity vector components and the previous phase's terminal velocity vector components
				//
				//Initial velocity:
				for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
				{
					if ((*Xdescriptions)[entry].find("initial velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_constraints_X_indices.push_back(entry);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
					}
				}
				//Terminal velocity:
				int foundcount = 0;
				for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
				{
					if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby velocity magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_constraints_X_indices.push_back(entry);
						flyby_velocity_magnitude_constraint_G_indices.push_back(iGfun->size() - 1);
						flyby_constraints_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
						++foundcount;
					}
					if (foundcount >= 3)
						break;
				}

				//we also need to encode a no-collision constraint
				if (Universe->bodies[boundary1_location_code - 1].mass < 1.0e+25)
					Flowerbounds->push_back(-10.0);
				else
					Flowerbounds->push_back(-300.0);
				Fupperbounds->push_back(0.0);
				Fdescriptions->push_back(prefix + "flyby altitude constraint (above minimum altitude but below [10x/300x] altitude for [rocky/gas] planets");

				//Jacobian entry (nonlinear) for the flyby no-collision constraint
				//the no-collision constraint is dependent on this phase's initial velocity vector components and the previous phase's terminal velocity vector components
				//
				//Initial velocity:
				for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
				{
					if ((*Xdescriptions)[entry].find("initial velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby no-collision constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_altitude_constraint_G_indices.push_back(iGfun->size() - 1);
					}
				}
				//Terminal velocity:
				foundcount = 0;
				for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
				{
					if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
					{
						iGfun->push_back(Fdescriptions->size() - 1);
						jGvar->push_back(entry);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << prefix << " flyby no-collision constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
						Gdescriptions->push_back(EntryNameStream.str());
						flyby_altitude_constraint_G_indices.push_back(iGfun->size() - 1);
						++foundcount;
					}
					if (foundcount >= 3)
						break;
				}
			}
		}
	}

	void phase::calcbounds_flight_time(const string& prefix, int first_X_entry_in_phase, vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p,  EMTG::Astrodynamics::universe* Universe, missionoptions* options)
	{
		if (boundary1_location_code > 0)
		{
			this->a1 = Universe->bodies[boundary1_location_code-1].SMA;
            this->e1 = Universe->bodies[boundary1_location_code - 1].ECC;
		}
		else if (boundary1_location_code == -1) //end at fixed or free point in space
		{
			if (options->journey_departure_elements_type[j] == 0) //orbit defined in inertial coordinates
			{
				double temp_coordinates[6];
				double temp_elements[6];
				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_departure_elements_vary_flag[j][k])
						temp_coordinates[k] = options->journey_departure_elements_bounds[j][k][1];
					else
						temp_coordinates[k] = options->journey_departure_elements[j][k];
				}

				Astrodynamics::inertial2COE(temp_coordinates, Universe->mu, temp_elements);
                this->a1 = temp_elements[0];
                this->e1 = temp_elements[1];
			}
			else if (options->journey_departure_elements_type[j] == 1) //orbit defined in classical orbital elements
			{
				double temp_elements[6];
				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_departure_elements_vary_flag[j][k])
						temp_elements[k] = options->journey_departure_elements_bounds[j][k][1];
					else
						temp_elements[k] = options->journey_departure_elements[j][k];
				}

                this->a1 = temp_elements[0];
                this->e1 = temp_elements[1];
			}
		}

		if (boundary2_location_code > 0)
		{
            this->a2 = Universe->bodies[boundary2_location_code - 1].SMA;
            this->e2 = Universe->bodies[boundary2_location_code - 1].ECC;
		}
		else if (boundary2_location_code == -1) //end at fixed or free point in space
		{
			if (options->journey_arrival_elements_type[j] == 0) //orbit defined in inertial coordinates
			{
				double temp_coordinates[6];
				double temp_elements[6];
				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_arrival_elements_vary_flag[j][k])
						temp_coordinates[k] = options->journey_arrival_elements_bounds[j][k][1];
					else
						temp_coordinates[k] = options->journey_arrival_elements[j][k];
				}

				Astrodynamics::inertial2COE(temp_coordinates, Universe->mu, temp_elements);
                this->a2 = temp_elements[0];
                this->e2 = temp_elements[1];
			}
			else if (options->journey_arrival_elements_type[j] == 1) //orbit defined in classical orbital elements
			{
				double temp_elements[6];
				for (int k = 0; k < 6; ++k)
				{
					if (options->journey_arrival_elements_vary_flag[j][k])
						temp_elements[k] = options->journey_arrival_elements_bounds[j][k][1];
					else
						temp_elements[k] = options->journey_arrival_elements[j][k];
				}

                this->a2 = temp_elements[0];
                this->e2 = temp_elements[1];
			}
		}

        if (this->e1 < 1.0)
            this->pseudoa1 = this->a1 * (1 + this->e1);
		else
            this->pseudoa1 = Universe->r_SOI / 5.0;

        if (this->e2 < 1.0)
            this->pseudoa2 = this->a2 * (1 + e2);
		else
            this->pseudoa2 = Universe->r_SOI / 5.0;

        this->T1 = 2 * math::PI*sqrt(this->pseudoa1*this->pseudoa1*this->pseudoa1 / Universe->mu);// pseudo-period of body 1 in days
        this->T2 = 2 * math::PI*sqrt(this->pseudoa2*this->pseudoa2*this->pseudoa2 / Universe->mu);// pseudo-period of body 2 in days
	
		double forced_coast_this_phase = 0.0;
		if (p == 0 && j == 0)
			forced_coast_this_phase += options->forced_post_launch_coast;
		else if (p > 0 || (p == 0 && (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4)) )
			forced_coast_this_phase += options->forced_flyby_coast;
		if (p < options->number_of_phases[j] - 1 || (p == options->number_of_phases[j] - 1 && (options->journey_arrival_type[j] == 2 || options->journey_arrival_type[j] == 0)) )
			forced_coast_this_phase += options->forced_flyby_coast;

		if (boundary1_location_code == boundary2_location_code && boundary1_location_code > 0) //if this transfer is a repeat of the same planet, we have special rules
		{
            double lowerbound_temp = this->T1 * 0.5;
			Xlowerbounds->push_back(lowerbound_temp > forced_coast_this_phase ? lowerbound_temp : forced_coast_this_phase);
            Xupperbounds->push_back(this->T1 * 20.0);
		}
		else if (boundary1_location_code == -1 && boundary2_location_code == -1) //for transfers between two free or fixed orbits
		{
			double lowerbound_temp = 1.0;
			Xlowerbounds->push_back(lowerbound_temp > forced_coast_this_phase ? lowerbound_temp : forced_coast_this_phase);
            Xupperbounds->push_back(max(this->T1, this->T2) * 20.0);
		}
		else
		{
			//lower bound is the same for all non-resonant phases
            double lowerbound_temp = 0.1 * min(this->T1, this->T2);

			lowerbound_temp = lowerbound_temp > forced_coast_this_phase ? lowerbound_temp : forced_coast_this_phase;
			Xlowerbounds->push_back(lowerbound_temp > 10.0 *  Universe->TU ? 10.0 *  Universe->TU : lowerbound_temp);

			if (max(pseudoa1,pseudoa2)/Universe->LU < 2.0) //outermost body is an inner body with a < 2 LU
                Xupperbounds->push_back(2.0 * max(this->T1, this->T2) < 45.0 * Universe->TU ? 45.0 * Universe->TU : 2.0 * max(this->T1, this->T2));
			 
			else //outermost body is an outer body
                Xupperbounds->push_back(1.0 * max(this->T1, this->T2) < 45.0 * Universe->TU ? 45.0 * Universe->TU : 1.0 * max(this->T1, this->T2));
		}

		Xdescriptions->push_back(prefix + "phase flight time");

		//compute the synodic period of the boundary points, for use in the MBH synodic period perturbation
		//these are "true" periods, not the pseudo-periods used for computing the bounds
        this->T1 = 2 * math::PI*sqrt(this->a1*this->a1*this->a1 / Universe->mu);
        this->T2 = 2 * math::PI*sqrt(this->a2*this->a2*this->a2 / Universe->mu);
        synodic_periods->push_back(1.0 / (fabs(1.0 / this->T1 - 1.0 / this->T2)));
	}

	void phase::calcbounds_right_boundary(const string& prefix, int first_X_entry_in_phase, vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, int j, int p,  EMTG::Astrodynamics::universe* Universe, missionoptions* options)
	{
		//if we are the last phase in the journey, then encode any variables necessary for the right hand boundary condition
		if (p == (options->number_of_phases[j] - 1))
		{
			//first, variables necessary to choose the boundary point
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
		}

		//if this phase is a terminal intercept, then encode the terminal velocity increment
		if (p == options->number_of_phases[j] - 1 && options->journey_arrival_type[j] == 2)
		{
			Xlowerbounds->push_back(-options->journey_final_velocity[j][1]);
			Xupperbounds->push_back(options->journey_final_velocity[j][1]);
			Xdescriptions->push_back(prefix + "terminal velocity increment x");
			Xlowerbounds->push_back(-options->journey_final_velocity[j][1]);
			Xupperbounds->push_back(options->journey_final_velocity[j][1]);
			Xdescriptions->push_back(prefix + "terminal velocity increment y");
			Xlowerbounds->push_back(-options->journey_final_velocity[j][1]);
			Xupperbounds->push_back(options->journey_final_velocity[j][1]);
			Xdescriptions->push_back(prefix + "terminal velocity increment z");
		}
		//if this phase is NOT a terminal rendezvous OR a final velocity match, then encode the terminal velocity increment
		//also do not encode this if the escape condition (arrival type 6) is set or if the journey ends in a spiral
		else if (!(p == options->number_of_phases[j] - 1 && (options->journey_arrival_type[j] == 4 || options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 5 || options->journey_arrival_type[j] == 6 || options->journey_arrival_type[j] == 7)))
		{
			Xlowerbounds->push_back(-25.0);
			Xupperbounds->push_back(25.0);
			Xdescriptions->push_back(prefix + "terminal velocity increment x");
			Xlowerbounds->push_back(-25.0);
			Xupperbounds->push_back(25.0);
			Xdescriptions->push_back(prefix + "terminal velocity increment y");
			Xlowerbounds->push_back(-25.0);
			Xupperbounds->push_back(25.0);
			Xdescriptions->push_back(prefix + "terminal velocity increment z");
		}
		

		//arrival mass
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

		//capture spirals need to know the phase end-mass before the spiral, so the spiral Isp is the last item to be encoded
		if (options->journey_arrival_type[j] == 7)//for journeys which end with a spiral - encode the spiral Isp if appropriate
		{
			if (options->engine_type == 4 || options->engine_type == 13) //generic VSI or Xenon Hall thruster
			{
				Xlowerbounds->push_back(options->IspLT_minimum);
				Xupperbounds->push_back(options->IspLT);
				Xdescriptions->push_back(prefix + "Arrival spiral Isp");
			}
			else if (options->engine_type == 12) //VASIMR
			{
				Xlowerbounds->push_back(3000.0);
				Xupperbounds->push_back(5000.0);
				Xdescriptions->push_back(prefix + "Arrival spiral Isp");
			}
		}
	}

	//function to find dependencies of a constraint on other variables due to a spiral at the beginning of any preceeding journey
	void phase::find_dependencies_due_to_escape_spiral(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, int j, int p, missionoptions* options)
	{
		//loop over journeys
		for (int jj = 0; jj <= j; ++jj)
		{
			//the first phase of the journey already has the journey-beginning variables "baked in" to it and so does not need any additional entries
			//also we don't have any dependencies if there is no escape spiral in that journey
			if ( !( jj == j && p == 0) && options->journey_departure_type[jj] == 5)
			{
				//loop over all variables in the decision vector and determine the first entry in the journey of interest
				int first_entry_in_jj;
				stringstream pjprefix_stream;
				pjprefix_stream << "j" << jj << "p";
				string pjprefix = pjprefix_stream.str();
				for (size_t Xentry = 0; Xentry < Xdescriptions->size() - 1; ++Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find(pjprefix) < 1024)
					{
						first_entry_in_jj = Xentry;
						break;
					}
				}//end loop over all variables in the decision vector

				//this constraint has a derivative with respect to the mass at the beginning of any journey that has an escape spiral
				//therefore we have derivatives with respect to:
				//1. if the first journey, the initial mass scale factor if enabled
				if (jj == 0 && options->allow_initial_mass_to_vary)
				{
					for (size_t Xentry = 0; Xentry < Xdescriptions->size() - 1; ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("initial mass multiplier (0-1)") < 1024)
						{
							iGfun->push_back(Fdescriptions->size() - 1);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
							break;
						}
					}
				}
				//2. if NOT the first journey, the arrival mass at the end of the previous journey
				if (jj > 0)
				{
					for (int Xentry = first_entry_in_jj; Xentry > 0; --Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("arrival mass") < 1024)
						{
							//first check for duplicates
							bool duplicateflag = false;
							stringstream entry_tag_stream;
							entry_tag_stream << "F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
							for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
							{
								if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}
				//3. if present, the journey initial mass increment scale factor
				//we must first check for duplicates associated with the current constraint, because this can occur
				
				for (size_t Xentry = first_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find("journey initial mass scale factor") < 1024)
					{
						bool duplicateflag = false;
						for (int XXentry = Xdescriptions->size() - 1; XXentry > 0; --XXentry)
						{
							stringstream tempstream;
							tempstream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
							if ( (*Xdescriptions)[XXentry].find("journey initial mass scale factor") < 1024)
							{
								duplicateflag = true;
								break;
							}
						}
						if (!duplicateflag)
						{
							iGfun->push_back(Fdescriptions->size() - 1);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
							break;
						}
					}

				}

				//this constraint has a derivative with respect to the Isp at the beginning of any journey that has an escape spiral
				if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
				{
					for (size_t Xentry = first_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("Escape spiral Isp") < 1024)
						{
							iGfun->push_back(Fdescriptions->size() - 1);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
							break;
						}
					}
				}

				//all spirals have a dependency on the BOL power if it is a variable
				if (options->objective_type == 13)
				{
					for (size_t Xentry = first_entry_in_jj; Xentry < Xdescriptions->size(); ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("engine input power (kW)") < 1024 )
						{
							bool duplicateflag = false;
							stringstream entry_tag_stream;
							entry_tag_stream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";

							for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
							{
								if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}
			}
		}//end loop over journeys
	}

	//function to find dependencies of a constraint on other variables due to a spiral at the end of any preceeding journey
	void phase::find_dependencies_due_to_capture_spiral(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, int j, int p, missionoptions* options)
	{
		//loop over journeys
		for (int jj = 0; jj < j; ++jj)
		{
			//we have a dependency on the arrival mass and the capture spiral Isp (if applicable) for each preceding journey that has a capture spiral
			if (options->journey_arrival_type[jj] == 7)
			{
				//loop over all variables in the decision vector and determine the first entry in the journey of interest
				int last_entry_in_jj;
				stringstream pjprefix_stream;
				pjprefix_stream << "j" << jj << "p";
				string pjprefix = pjprefix_stream.str();
				for (int Xentry = Xdescriptions->size() - 1; Xentry > 0 ; --Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find(pjprefix) < 1024)
					{
						last_entry_in_jj = Xentry;
						break;
					}
				}//end loop over all variables in the decision vector


				//this constraint has a derivative with respect to the arrival mass at the end of any journey that has a capture spiral
				//the exception is that the first phase of a journey already includes an entry for the arrival mass at the end of the previous journey
				//so we don't need to count it twice
				if (!( jj == (j - 1) && p == 0))
				{
				
					for (int Xentry = last_entry_in_jj; Xentry > 0; --Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("arrival mass") < 1024)
						{
							//first check for duplicates
							bool duplicateflag = false;
							stringstream entry_tag_stream;
							entry_tag_stream << "F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
							for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
							{
								if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}

				//this constraint has a derivative with respect to the Isp at the beginning of any journey that has a capture spiral
				if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
				{
					for (size_t Xentry = last_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("Capture spiral Isp") < 1024)
						{
							iGfun->push_back(Fdescriptions->size() - 1);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
							break;
						}
					}
				}	

				//all spirals have a dependency on the BOL power if it is a variable
				if (options->objective_type == 13)
				{
					for (size_t Xentry = last_entry_in_jj; Xentry < Xdescriptions->size() - 1; --Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("engine input power (kW)") < 1024 )
						{
							bool duplicateflag = false;
							for (int XXentry = Xdescriptions->size() - 1; XXentry > 0; --XXentry)
							{
								stringstream tempstream;
								tempstream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
								if ( (*Xdescriptions)[XXentry].find("engine input power (kW)") < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}
			}

			
		}//end loop over journeys
	}

    void phase::calcbounds_phase_thruster_parameters(const string& prefix, int first_X_entry_in_phase, vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
    {
        if (options->engine_type == 1 && (j == 0 && p == 0))
        {
            //constant Isp, efficiency, EMTG computes input power
            Xlowerbounds->push_back(options->engine_input_power_bounds[0]);
            Xupperbounds->push_back(options->engine_input_power_bounds[1]);
            Xdescriptions->push_back(prefix + "engine input power (kW)");
        }
        else if (options->engine_type == 2 && (j == 0 && p == 0))
        {
            //constant power, EMTG chooses Isp
            Xlowerbounds->push_back(options->IspLT_minimum);
            Xupperbounds->push_back(options->IspLT);
            Xdescriptions->push_back(prefix + "engine Isp (s)");
        }
        else if (options->objective_type == 13 && j == 0 && p == 0)
        {
            //EMTG varies input power for whatever engine/power model you have
            Xlowerbounds->push_back(math::SMALL);
            Xupperbounds->push_back(options->power_at_1_AU);
            Xdescriptions->push_back(prefix + "engine input power (kW)");
        }
    }
} /* namespace EMTG */
