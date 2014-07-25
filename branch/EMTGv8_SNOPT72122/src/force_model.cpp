//force model for EMTGv8
//Jacob Englander 4/4/2013

#include "Astrodynamics.h"
#include "EMTG_math.h"
#include "EMTG_Matrix.h"

namespace EMTG { namespace Astrodynamics {
	int force_model(EMTG::missionoptions* options,
					EMTG::Astrodynamics::universe* Universe,
					double* spacecraft_state,
					double* epoch,
					double* launch_epoch,
					double* control,
					double* max_thrust,
					double* max_mass_flow_rate,
					double* Isp,
					double* power,
					double* active_power,
					int* number_of_active_engines,
					double* force_vector,
					const bool& generate_derivatives,
					double* dTdP,
					double* dmdotdP,
					double* dTdIsp,
					double* dmdotdIsp,
					double* dPdr,
					double* dPdt,
					double* dFSRPdr,
					vector<double>& dagravdRvec,
					vector<double>& dagravdtvec)
	{
		//Note: all thrusts are returned in Newtons

		//we need to know how far we are from the Sun
		double distance_from_sun_in_AU;
		double position_relative_to_sun_in_AU[3];

		if (!(Universe->central_body_SPICE_ID == 10))
		{
			double central_body_state_in_km[6];	
			Universe->locate_central_body(*epoch, central_body_state_in_km, options);

			position_relative_to_sun_in_AU[0] = (central_body_state_in_km[0] + spacecraft_state[0]) / options->AU;
			position_relative_to_sun_in_AU[1] = (central_body_state_in_km[1] + spacecraft_state[1]) / options->AU;
			position_relative_to_sun_in_AU[2] = (central_body_state_in_km[2] + spacecraft_state[2]) / options->AU;

			distance_from_sun_in_AU = math::norm(position_relative_to_sun_in_AU, 3);
		}
		else //if we are orbiting the sun, it's quite silly to look up the position of the sun relative to itself, so don't bother
			distance_from_sun_in_AU = math::norm(spacecraft_state, 3)/Universe->LU;

		//compute the thrust available from the engine
		EMTG::Astrodynamics::find_engine_parameters(options,
													distance_from_sun_in_AU,
													*epoch - *launch_epoch,
													max_thrust,
													max_mass_flow_rate,
													Isp, 
													power,
													active_power,
													number_of_active_engines,
													generate_derivatives,
													dTdP,
													dmdotdP,
													dTdIsp,
													dmdotdIsp,
													dPdr,
													dPdt);

		//modify the thrust by the duty cycle of the engine
		double applied_thrust = *max_thrust * options->engine_duty_cycle;

		//compute thrust in Newtons
		force_vector[0] = control[0] * applied_thrust;
		force_vector[1] = control[1] * applied_thrust;
		force_vector[2] = control[2] * applied_thrust;

		//compute Solar Radiation Pressure if applicable
		if (options->perturb_SRP)
		{
			double F_SRP = -(4.60689e-6 * options->spacecraft_area * options->coefficient_of_reflectivity) / (distance_from_sun_in_AU*distance_from_sun_in_AU * 1000000);
			
			force_vector[0] += F_SRP * position_relative_to_sun_in_AU[0] / distance_from_sun_in_AU;
			force_vector[1] += F_SRP * position_relative_to_sun_in_AU[1] / distance_from_sun_in_AU;
			force_vector[2] += F_SRP * position_relative_to_sun_in_AU[2] / distance_from_sun_in_AU;

			*dFSRPdr = 2.0 * (4.60689e-6 * options->spacecraft_area * options->coefficient_of_reflectivity) / (distance_from_sun_in_AU*distance_from_sun_in_AU*distance_from_sun_in_AU * 1000000);
		}

		//compute third body perturbations if applicable
		if (options->perturb_thirdbody)
		{
			double Fgravity = 0;
			dagravdRvec[0] = 0.0;
			dagravdRvec[1] = 0.0;
			dagravdRvec[2] = 0.0;
			dagravdtvec[0] = 0.0;
			dagravdtvec[1] = 0.0;
			dagravdtvec[2] = 0.0;

			//perturbations due to the sun
			//don't bother calculating these if we are orbiting the sun; that would be silly
			if (!(Universe->central_body_SPICE_ID == 10))
			{
				double distance_from_sun_in_km = distance_from_sun_in_AU * options->AU;
				Fgravity = -1.32712428e11 / (distance_from_sun_in_km * distance_from_sun_in_km) * 1000;

				force_vector[0] += Fgravity * position_relative_to_sun_in_AU[0] / distance_from_sun_in_AU;
				force_vector[1] += Fgravity * position_relative_to_sun_in_AU[1] / distance_from_sun_in_AU;
				force_vector[2] += Fgravity * position_relative_to_sun_in_AU[2] / distance_from_sun_in_AU;

				if (generate_derivatives)
				{
					//3rd body derivatives
					double dadr_coeff = 3 * 1.32712428e11 / (distance_from_sun_in_km*distance_from_sun_in_km*distance_from_sun_in_km*distance_from_sun_in_km);
					double dadxyz = -1.32712428e11 / (distance_from_sun_in_km*distance_from_sun_in_km*distance_from_sun_in_km);

					//note: "thingy" is whatever variable is of interest, be it position, velocity, or time - they are all the same. See Donald's notes 8/29/2013
					double dadthingy_x = (dadr_coeff * position_relative_to_sun_in_AU[0]*position_relative_to_sun_in_AU[0] / distance_from_sun_in_AU)/options->AU + dadxyz;
					double dadthingy_y = (dadr_coeff * position_relative_to_sun_in_AU[1]*position_relative_to_sun_in_AU[2] / distance_from_sun_in_AU)/options->AU + dadxyz;
					double dadthingy_z = (dadr_coeff * position_relative_to_sun_in_AU[1]*position_relative_to_sun_in_AU[2] / distance_from_sun_in_AU)/options->AU + dadxyz;

					dagravdRvec[0] += dadthingy_x;
					dagravdRvec[1] += dadthingy_y;
					dagravdRvec[2] += dadthingy_z;

					dagravdtvec[0] -= dadthingy_x;
					dagravdtvec[1] -= dadthingy_y;
					dagravdtvec[2] -= dadthingy_z;
				}
			}

			//perturbations due to bodies in the flyby menu
			for (size_t b = 0; b < Universe->perturbation_menu.size(); ++b)
			{
				double body_state_in_km[9];
				double spacecraft_position_relative_to_body_in_km[3];
				double distance_from_body_in_km;

				Universe->bodies[Universe->perturbation_menu[b]].locate_body(*epoch,
																			body_state_in_km,
																			generate_derivatives && options->derivative_type > 2,
																			options);

				spacecraft_position_relative_to_body_in_km[0] = spacecraft_state[0] - body_state_in_km[0];
				spacecraft_position_relative_to_body_in_km[1] = spacecraft_state[1] - body_state_in_km[1];
				spacecraft_position_relative_to_body_in_km[2] = spacecraft_state[2] - body_state_in_km[2];

				distance_from_body_in_km = math::norm(spacecraft_position_relative_to_body_in_km, 3);

				//to avoid singularities and to avoid screwing up the flyby model, we will only include a third body perturbations when we are not too close to the body
				if (distance_from_body_in_km > Universe->bodies[Universe->perturbation_menu[b]].radius + Universe->bodies[Universe->perturbation_menu[b]].minimum_safe_flyby_altitude)
				{
					Fgravity = -Universe->bodies[Universe->perturbation_menu[b]].mu / (distance_from_body_in_km*distance_from_body_in_km) * 1000;

					force_vector[0] += Fgravity * spacecraft_position_relative_to_body_in_km[0] / distance_from_body_in_km;
					force_vector[1] += Fgravity * spacecraft_position_relative_to_body_in_km[1] / distance_from_body_in_km;
					force_vector[2] += Fgravity * spacecraft_position_relative_to_body_in_km[2] / distance_from_body_in_km;

					if (generate_derivatives)
					{
						//3rd body derivatives
						double dadr_coeff = 3 * Universe->bodies[Universe->perturbation_menu[b]].mu / (distance_from_body_in_km*distance_from_body_in_km*distance_from_body_in_km*distance_from_body_in_km);
						double dadxyz = -Universe->bodies[Universe->perturbation_menu[b]].mu / (distance_from_body_in_km*distance_from_body_in_km*distance_from_body_in_km);

						//note: "thingy" is whatever variable is of interest, be it position, velocity, or time - they are all the same. See Donald's notes 8/29/2013
						double dadthingy_x = (dadr_coeff * spacecraft_position_relative_to_body_in_km[0]*spacecraft_position_relative_to_body_in_km[0] / distance_from_body_in_km) + dadxyz;
						double dadthingy_y = (dadr_coeff * spacecraft_position_relative_to_body_in_km[1]*spacecraft_position_relative_to_body_in_km[2] / distance_from_body_in_km) + dadxyz;
						double dadthingy_z = (dadr_coeff * spacecraft_position_relative_to_body_in_km[1]*spacecraft_position_relative_to_body_in_km[2] / distance_from_body_in_km) + dadxyz;

						dagravdRvec[0] = dadthingy_x;
						dagravdRvec[1] = dadthingy_y;
						dagravdRvec[2] = dadthingy_z;

						dagravdtvec[0] -= dadthingy_x;
						dagravdtvec[1] -= dadthingy_y;
						dagravdtvec[2] -= dadthingy_z;
					}
				}
			}
		}

		return 0;
	}
}} //close namespace
