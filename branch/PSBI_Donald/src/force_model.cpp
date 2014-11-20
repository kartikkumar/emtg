//force model for EMTGv8
//Jacob Englander 4/4/2013

#include <iomanip>
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
					vector<double>& dagravdtvec,
                    vector<double>& central_body_state_mks)
	{
		//Note: all thrusts are returned in Newtons

		//we need to know how far we are from the Sun
		double distance_from_sun_in_AU;
		double position_relative_to_sun_in_AU[3];

		if (!(Universe->central_body_SPICE_ID == 10))
		{
            Universe->locate_central_body(*epoch, central_body_state_mks.data(), options, generate_derivatives);

            position_relative_to_sun_in_AU[0] = (central_body_state_mks[0] + spacecraft_state[0]) / options->AU;
            position_relative_to_sun_in_AU[1] = (central_body_state_mks[1] + spacecraft_state[1]) / options->AU;
            position_relative_to_sun_in_AU[2] = (central_body_state_mks[2] + spacecraft_state[2]) / options->AU;

			distance_from_sun_in_AU = math::norm(position_relative_to_sun_in_AU, 3);
		}
        else //if we are orbiting the sun, it's quite silly to look up the position of the sun relative to itself, so don't bother
        {
            distance_from_sun_in_AU = math::norm(spacecraft_state, 3) / Universe->LU;
            for (size_t k = 0; k < (generate_derivatives && options->derivative_type > 2? 12 : 6); ++k)
                central_body_state_mks[k] = 0.0;
        }

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

	int FBLT_force_model(EMTG::missionoptions* options,
		EMTG::Astrodynamics::universe* Universe,
		double* spacecraft_state_relative_to_central_body_in_km,
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
		EMTG::math::Matrix<double> & A,
		vector<double>& dagravdRvec,
		vector<double>& dagravdtvec,
        vector<double>& central_body_state_mks)
	{
		//Note: all thrusts are returned in Newtons

		bool normalized_STMs = true;

		//we need to know how far we are from the Sun
		double spacecraft_distance_from_sun_in_LU;
		double spacecraft_distance_from_sun_in_AU;
		double spacecraft_distance_from_sun_in_km;
		double spacecraft_position_relative_to_sun_in_LU[3];
		double spacecraft_position_relative_to_sun_in_AU[3];
		double spacecraft_position_relative_to_sun_in_km[3];

		//we need the spacecraft's state in LU/TU
		double spacecraft_state_relative_to_central_body_in_LU[7];
		double spacecraft_distance_from_central_body_in_LU;
		double spacecraft_distance_from_central_body_in_km;

		for (size_t k = 0; k < 3; ++k)
		{
			spacecraft_state_relative_to_central_body_in_LU[k] = spacecraft_state_relative_to_central_body_in_km[k] / Universe->LU;
			spacecraft_state_relative_to_central_body_in_LU[k + 3] = spacecraft_state_relative_to_central_body_in_km[k + 3] / Universe->LU * Universe->TU;
		}
		spacecraft_distance_from_central_body_in_LU = math::norm(spacecraft_state_relative_to_central_body_in_LU, 3);
		spacecraft_distance_from_central_body_in_km = math::norm(spacecraft_state_relative_to_central_body_in_km, 3);

		if (!(Universe->central_body_SPICE_ID == 10))
		{
			//locate the central body's position w.r.t. the Sun
            Universe->locate_central_body(*epoch, central_body_state_mks.data(), options, generate_derivatives);

            spacecraft_position_relative_to_sun_in_km[0] = (central_body_state_mks[0] + spacecraft_state_relative_to_central_body_in_km[0]);
            spacecraft_position_relative_to_sun_in_km[1] = (central_body_state_mks[1] + spacecraft_state_relative_to_central_body_in_km[1]);
            spacecraft_position_relative_to_sun_in_km[2] = (central_body_state_mks[2] + spacecraft_state_relative_to_central_body_in_km[2]);

			spacecraft_position_relative_to_sun_in_LU[0] = spacecraft_position_relative_to_sun_in_km[0] / Universe->LU;
			spacecraft_position_relative_to_sun_in_LU[1] = spacecraft_position_relative_to_sun_in_km[1] / Universe->LU;
			spacecraft_position_relative_to_sun_in_LU[2] = spacecraft_position_relative_to_sun_in_km[2] / Universe->LU;

			spacecraft_distance_from_sun_in_km = math::norm(spacecraft_position_relative_to_sun_in_km, 3);
			spacecraft_distance_from_sun_in_AU = spacecraft_distance_from_sun_in_km / options->AU;
			spacecraft_distance_from_sun_in_LU = spacecraft_distance_from_sun_in_km / Universe->LU;
		}
		else //if we are orbiting the sun, it's quite silly to look up the position of the sun relative to itself, so don't bother
		{		
			spacecraft_distance_from_sun_in_km = math::norm(spacecraft_state_relative_to_central_body_in_km, 3);
			spacecraft_distance_from_sun_in_AU = spacecraft_distance_from_sun_in_km / options->AU;
		}

		double mass_kg = spacecraft_state_relative_to_central_body_in_km[6];
		double mass_normalized = spacecraft_state_relative_to_central_body_in_km[6] / options->maximum_mass;

        /*
		double r_pert = 1.0e-6;
		double Pforward;
		double Pbackward;
		double rtemp = spacecraft_distance_from_sun_in_AU;

		for (size_t loopCount = 0; loopCount < 3; ++loopCount)
		{

			if (loopCount == 0)
				spacecraft_distance_from_sun_in_AU += r_pert;
			else if (loopCount == 1)
			{
				spacecraft_distance_from_sun_in_AU -= 2.0*r_pert;
			}
			else
				spacecraft_distance_from_sun_in_AU = rtemp;
                
			//compute the maximum thrust available from the engines  
			EMTG::Astrodynamics::find_engine_parameters(options,
				spacecraft_distance_from_sun_in_AU,
				*epoch - *launch_epoch,
				max_thrust, //kN
				max_mass_flow_rate, // kg/s
				Isp, // seconds
				power, // kW
				active_power, // kW
				number_of_active_engines,
				generate_derivatives,
				dTdP, // kN/kW 
				dmdotdP, // kg/kW
				dTdIsp, // kN/s
				dmdotdIsp, // kg/s
				dPdr, // kW/AU
				dPdt // kW/s
				);
            
		if (loopCount == 0)
			Pforward = *power;
		else if (loopCount == 1)
			Pbackward = *power;

		} //end finite difference loop
        */
        //compute the maximum thrust available from the engines  
        EMTG::Astrodynamics::find_engine_parameters(options,
                                                    spacecraft_distance_from_sun_in_AU,
                                                    *epoch - *launch_epoch,
                                                    max_thrust, //kN
                                                    max_mass_flow_rate, // kg/s
                                                    Isp, // seconds
                                                    power, // kW
                                                    active_power, // kW
                                                    number_of_active_engines,
                                                    generate_derivatives,
                                                    dTdP, // kN/kW 
                                                    dmdotdP, // kg/kW
                                                    dTdIsp, // kN/s
                                                    dmdotdIsp, // kg/s
                                                    dPdr, // kW/AU
                                                    dPdt // kW/s
                                                    );
        /*
		double dPdr_FD = (Pforward - Pbackward) / (2.0 * r_pert);
		std::cout << setprecision(16);
	    std::cout << "Finite differenced dPdr: " << dPdr_FD << std::endl;
		std::cout << "Analytical dPdr: " << (*dPdr) << std::endl;
        std::cout << "error: " << dPdr_FD - (*dPdr) << std::endl;
		//getchar();
        */
        

		//depend. Theing on whether we want normalized or MKS STMs we must convert
		//the engine model derivative units
		if (normalized_STMs)
		{
			*max_thrust = (*max_thrust) / options->maximum_mass / Universe->LU * Universe->TU * Universe->TU;
			*max_mass_flow_rate = (*max_mass_flow_rate) / options->maximum_mass * Universe->TU;
			//we don't need to scale Isp, because none of the EOM's directly depend on it
			//we don't want to scale Isp, because it will mess up the output to file
			//*Isp = (*Isp) / Universe->TU;
			*dTdP = (*dTdP) / options->maximum_mass / Universe->LU * Universe->TU * Universe->TU;
			*dmdotdP = (*dmdotdP) / options->maximum_mass;
			*dTdIsp = (*dTdIsp) / options->maximum_mass / Universe->LU * Universe->TU * Universe->TU * Universe->TU;
			*dmdotdIsp = (*dmdotdIsp) / options->maximum_mass * Universe->TU;
			*dPdr = (*dPdr) / options->AU * Universe->LU;
			*dPdt = (*dPdt) * Universe->TU;
		}
		//if we don't want normalized units, we still need to fix one derivative that is in kW/AU
		else
		{
			*dPdr = (*dPdr) / options->AU;
		}

		//modify the thrust by the duty cycle of the engine
		double maximum_available_thrust = (*max_thrust) * options->engine_duty_cycle;

		//compute thrust...control vector consists of throttle decision parameters for this FBLT segment
		force_vector[0] = control[0] * maximum_available_thrust;
		force_vector[1] = control[1] * maximum_available_thrust;
		force_vector[2] = control[2] * maximum_available_thrust;

		double epsilon = 1.0e-10; //avoid divide by zero if thruster is off
		double control_norm = EMTG::math::norm(control, 3) + epsilon;

		//****************************************
		//
		//Auxilliary Derivatives of Interest
		//
		//****************************************
		

		double one_over_spacecraft_distance_from_central_body_in_km = 1.0 / spacecraft_distance_from_central_body_in_km;
		double one_over_spacecraft_distance_from_central_body_in_LU = 1.0 / spacecraft_distance_from_central_body_in_LU;

		double one_over_control_norm_plus_epsilon = 1.0 / control_norm;

		//gradient of the magnitude of the spacecraft position vector
		double drdx = spacecraft_state_relative_to_central_body_in_km[0] * one_over_spacecraft_distance_from_central_body_in_km;
		double drdy = spacecraft_state_relative_to_central_body_in_km[1] * one_over_spacecraft_distance_from_central_body_in_km;
		double drdz = spacecraft_state_relative_to_central_body_in_km[2] * one_over_spacecraft_distance_from_central_body_in_km;

		double drdx_normalized = spacecraft_state_relative_to_central_body_in_LU[0] * one_over_spacecraft_distance_from_central_body_in_LU;
		double drdy_normalized = spacecraft_state_relative_to_central_body_in_LU[1] * one_over_spacecraft_distance_from_central_body_in_LU;
		double drdz_normalized = spacecraft_state_relative_to_central_body_in_LU[2] * one_over_spacecraft_distance_from_central_body_in_LU;

		//gradient of the magnitude of the control vector
		double dcontrol_normdux = control[0] * one_over_control_norm_plus_epsilon;
		double dcontrol_normduy = control[1] * one_over_control_norm_plus_epsilon;
		double dcontrol_normduz = control[2] * one_over_control_norm_plus_epsilon;

		if (generate_derivatives)
		{
			//****************************************
			//
			//State Propagation (A) matrix calculation
			//
			//****************************************

			//Top Row Identity
			A(0, 3) = 1.0;
			A(1, 4) = 1.0;
			A(2, 5) = 1.0;

			//A21 dadr (everything except for the third body terms, they are handled below)
			
			if (normalized_STMs)
			{
				double three_muCB_over_spacecraft_distance_from_CB_in_LU5 = 3.0 * 1.0 / pow(spacecraft_distance_from_central_body_in_LU, 5.0);
				double muCB_over_spacecraft_distance_from_CB_in_LU3 = 1.0 / pow(spacecraft_distance_from_central_body_in_LU, 3.0);
				double D_dTdP_dPdr_over_msc = (options->engine_duty_cycle) * (*dTdP) * (*dPdr) / mass_normalized;
				double* State = spacecraft_state_relative_to_central_body_in_LU;

				A(3, 0) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[0] * State[0] - muCB_over_spacecraft_distance_from_CB_in_LU3 + control[0] * D_dTdP_dPdr_over_msc * drdx_normalized;
				A(3, 1) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[0] * State[1]												  + control[0] * D_dTdP_dPdr_over_msc * drdy_normalized;
				A(3, 2) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[0] * State[2]												  + control[0] * D_dTdP_dPdr_over_msc * drdz_normalized;
				A(4, 0) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[1] * State[0]												  + control[1] * D_dTdP_dPdr_over_msc * drdx_normalized;
				A(4, 1) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[1] * State[1] - muCB_over_spacecraft_distance_from_CB_in_LU3 + control[1] * D_dTdP_dPdr_over_msc * drdy_normalized;
				A(4, 2) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[1] * State[2]                                                + control[1] * D_dTdP_dPdr_over_msc * drdz_normalized;
				A(5, 0) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[2] * State[0]                                                + control[2] * D_dTdP_dPdr_over_msc * drdx_normalized;
				A(5, 1) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[2] * State[1]                                                + control[2] * D_dTdP_dPdr_over_msc * drdy_normalized;
				A(5, 2) = three_muCB_over_spacecraft_distance_from_CB_in_LU5 * State[2] * State[2] - muCB_over_spacecraft_distance_from_CB_in_LU3 + control[2] * D_dTdP_dPdr_over_msc * drdz_normalized;
			}
			else
			{
				double three_muCB_over_spacecraft_distance_from_CB_in_km5 = 3.0 * Universe->mu / pow(spacecraft_distance_from_central_body_in_km, 5.0);
				double muCB_over_spacecraft_distance_from_CB_in_km3 = Universe->mu / pow(spacecraft_distance_from_central_body_in_km, 3.0);
				double D_dTdP_dPdr_over_msc = (options->engine_duty_cycle) * (*dTdP) * (*dPdr) / mass_kg;
				double* State = spacecraft_state_relative_to_central_body_in_km;

				A(3, 0) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[0] * State[0] - muCB_over_spacecraft_distance_from_CB_in_km3 + control[0] * D_dTdP_dPdr_over_msc * drdx;
				A(3, 1) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[0] * State[1]												  + control[0] * D_dTdP_dPdr_over_msc * drdy;
				A(3, 2) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[0] * State[2]												  + control[0] * D_dTdP_dPdr_over_msc * drdz;
				A(4, 0) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[1] * State[0]												  + control[1] * D_dTdP_dPdr_over_msc * drdx;
				A(4, 1) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[1] * State[1] - muCB_over_spacecraft_distance_from_CB_in_km3 + control[1] * D_dTdP_dPdr_over_msc * drdy;
				A(4, 2) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[1] * State[2]                                                + control[1] * D_dTdP_dPdr_over_msc * drdz;
				A(5, 0) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[2] * State[0]                                                + control[2] * D_dTdP_dPdr_over_msc * drdx;
				A(5, 1) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[2] * State[1]                                                + control[2] * D_dTdP_dPdr_over_msc * drdy;
				A(5, 2) = three_muCB_over_spacecraft_distance_from_CB_in_km5 * State[2] * State[2] - muCB_over_spacecraft_distance_from_CB_in_km3 + control[2] * D_dTdP_dPdr_over_msc * drdz;
			}

			//A23 dadm
			if (normalized_STMs)
			{
				//maximum_available_thrust has already been normalized
				double D_thrust_over_msc2 = maximum_available_thrust / (mass_normalized * mass_normalized);

				A(3, 6) = -control[0] * D_thrust_over_msc2;
				A(4, 6) = -control[1] * D_thrust_over_msc2;
				A(5, 6) = -control[2] * D_thrust_over_msc2;
			}
			else
			{
				double D_thrust_over_msc2 = maximum_available_thrust / (mass_kg * mass_kg);

				A(3, 6) = -control[0] * D_thrust_over_msc2;
				A(4, 6) = -control[1] * D_thrust_over_msc2;
				A(5, 6) = -control[2] * D_thrust_over_msc2;
			}

			//A24 dadu
			if (normalized_STMs)
			{
				double D_thrust_over_msc = maximum_available_thrust / mass_normalized;

				A(3, 7) = D_thrust_over_msc;
				A(3, 8) = 0.0;
				A(3, 9) = 0.0;
				A(4, 7) = 0.0;
				A(4, 8) = D_thrust_over_msc;
				A(4, 9) = 0.0;
				A(5, 7) = 0.0;
				A(5, 8) = 0.0;
				A(5, 9) = D_thrust_over_msc;
			}
			else
			{
				double D_thrust_over_msc = maximum_available_thrust / mass_kg;

				A(3, 7) = D_thrust_over_msc;
				A(3, 8) = 0.0;
				A(3, 9) = 0.0;
				A(4, 7) = 0.0;
				A(4, 8) = D_thrust_over_msc;
				A(4, 9) = 0.0;
				A(5, 7) = 0.0;
				A(5, 8) = 0.0;
				A(5, 9) = D_thrust_over_msc;
			}

			//A31 dmdotdr
            double control_norm_D_dmdotdP_dPdr = control_norm * (options->engine_duty_cycle) * (*dmdotdP) * (*dPdr);

			if (normalized_STMs)
			{
				A(6, 0) = -control_norm_D_dmdotdP_dPdr * drdx_normalized;
				A(6, 1) = -control_norm_D_dmdotdP_dPdr * drdy_normalized;
				A(6, 2) = -control_norm_D_dmdotdP_dPdr * drdz_normalized;
			}
			else
			{
				A(6, 0) = -control_norm_D_dmdotdP_dPdr * drdx;
				A(6, 1) = -control_norm_D_dmdotdP_dPdr * drdy;
				A(6, 2) = -control_norm_D_dmdotdP_dPdr * drdz;
			}

			//A34 dmdotdu
            double D_mdot = (options->engine_duty_cycle) * (*max_mass_flow_rate);
            
			A(6, 7) = -dcontrol_normdux * D_mdot;
			A(6, 8) = -dcontrol_normduy * D_mdot;
			A(6, 9) = -dcontrol_normduz * D_mdot;
		}




		//compute third body perturbations if applicable
		if (options->perturb_thirdbody)
		{
			double Fgravity = 0.0;

			//perturbations due to the Sun
			//don't bother calculating these if we are orbiting the sun; that would be silly
			//if we are not orbiting the Sun, then we ALWAYS calculate Sun perturbations
			if (!(Universe->central_body_SPICE_ID == 10))
			{
				double spacecraft_distance_from_sun_in_km = spacecraft_distance_from_sun_in_AU * options->AU;

				// Force of gravity due to the Sun in Newtons
				double mu_sun_mks = 1.32712440018e11;

				Fgravity = -mu_sun_mks / (spacecraft_distance_from_sun_in_km * spacecraft_distance_from_sun_in_km) * 1000.0;


				force_vector[0] += Fgravity * spacecraft_position_relative_to_sun_in_km[0] / spacecraft_distance_from_sun_in_km;
				force_vector[1] += Fgravity * spacecraft_position_relative_to_sun_in_km[1] / spacecraft_distance_from_sun_in_km;
				force_vector[2] += Fgravity * spacecraft_position_relative_to_sun_in_km[2] / spacecraft_distance_from_sun_in_km;

				if (generate_derivatives)
				{
					if (normalized_STMs)
					{
						double mu_sun_normalized = 1.32712440018e11 / Universe->mu;
						double three_muSun_over_spacecraft_distance_from_sun_in_LU5 = 3.0 * mu_sun_normalized / pow(spacecraft_distance_from_sun_in_LU, 5.0);
						double muSun_over_spacecraft_distance_from_sun_in_LU3 = mu_sun_normalized / pow(spacecraft_distance_from_sun_in_LU, 3.0);

						A(3, 0) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[0] * spacecraft_position_relative_to_sun_in_LU[0] - muSun_over_spacecraft_distance_from_sun_in_LU3;
						A(3, 1) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[0] * spacecraft_position_relative_to_sun_in_LU[1];
						A(3, 2) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[0] * spacecraft_position_relative_to_sun_in_LU[2];
						A(4, 0) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[1] * spacecraft_position_relative_to_sun_in_LU[0];
						A(4, 1) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[1] * spacecraft_position_relative_to_sun_in_LU[1] - muSun_over_spacecraft_distance_from_sun_in_LU3;
						A(4, 2) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[1] * spacecraft_position_relative_to_sun_in_LU[2];
						A(5, 0) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[2] * spacecraft_position_relative_to_sun_in_LU[0];
						A(5, 1) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[2] * spacecraft_position_relative_to_sun_in_LU[1];
						A(5, 2) += three_muSun_over_spacecraft_distance_from_sun_in_LU5 * spacecraft_position_relative_to_sun_in_LU[2] * spacecraft_position_relative_to_sun_in_LU[2] - muSun_over_spacecraft_distance_from_sun_in_LU3;
					}
					else
					{
						double three_muSun_over_spacecraft_distance_from_sun_in_km5 = 3.0 * mu_sun_mks / pow(spacecraft_distance_from_sun_in_km, 5.0);
						double muSun_over_spacecraft_distance_from_sun_in_km3 = mu_sun_mks / pow(spacecraft_distance_from_sun_in_km, 3.0);

						A(3, 0) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[0] * spacecraft_position_relative_to_sun_in_km[0] - muSun_over_spacecraft_distance_from_sun_in_km3;
						A(3, 1) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[0] * spacecraft_position_relative_to_sun_in_km[1];
						A(3, 2) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[0] * spacecraft_position_relative_to_sun_in_km[2];
						A(4, 0) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[1] * spacecraft_position_relative_to_sun_in_km[0];
						A(4, 1) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[1] * spacecraft_position_relative_to_sun_in_km[1] - muSun_over_spacecraft_distance_from_sun_in_km3;
						A(4, 2) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[1] * spacecraft_position_relative_to_sun_in_km[2];
						A(5, 0) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[2] * spacecraft_position_relative_to_sun_in_km[0];
						A(5, 1) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[2] * spacecraft_position_relative_to_sun_in_km[1];
						A(5, 2) += three_muSun_over_spacecraft_distance_from_sun_in_km5 * spacecraft_position_relative_to_sun_in_km[2] * spacecraft_position_relative_to_sun_in_km[2] - muSun_over_spacecraft_distance_from_sun_in_km3;
					}

				}
			}

			//perturbations due to bodies in the perturbation list
			for (size_t b = 0; b < Universe->perturbation_menu.size(); ++b)
			{
				double third_body_state_relative_to_central_body_in_km[9];
				double spacecraft_position_relative_to_third_body_in_km[3];
				double spacecraft_distance_from_third_body_in_km;
				double spacecraft_position_relative_to_third_body_in_LU[3];
				double spacecraft_distance_from_third_body_in_LU;

				//extract the third body state vector from the ephemeris
				Universe->bodies[Universe->perturbation_menu[b]].locate_body(*epoch,
					third_body_state_relative_to_central_body_in_km,
					generate_derivatives && options->derivative_type > 2,
					options);

				//form the vector from the third body to the spacecraft in the central body reference frame
				spacecraft_position_relative_to_third_body_in_km[0] = spacecraft_state_relative_to_central_body_in_km[0] - third_body_state_relative_to_central_body_in_km[0];
				spacecraft_position_relative_to_third_body_in_km[1] = spacecraft_state_relative_to_central_body_in_km[1] - third_body_state_relative_to_central_body_in_km[1];
				spacecraft_position_relative_to_third_body_in_km[2] = spacecraft_state_relative_to_central_body_in_km[2] - third_body_state_relative_to_central_body_in_km[2];

				spacecraft_position_relative_to_third_body_in_LU[0] = spacecraft_position_relative_to_third_body_in_km[0] / Universe->LU;
				spacecraft_position_relative_to_third_body_in_LU[1] = spacecraft_position_relative_to_third_body_in_km[1] / Universe->LU;
				spacecraft_position_relative_to_third_body_in_LU[2] = spacecraft_position_relative_to_third_body_in_km[2] / Universe->LU;

				spacecraft_distance_from_third_body_in_km = math::norm(spacecraft_position_relative_to_third_body_in_km, 3);
				spacecraft_distance_from_third_body_in_LU = math::norm(spacecraft_position_relative_to_third_body_in_LU, 3);

				//to avoid singularities and to avoid screwing up the flyby model, we will only include a third body perturbations when we are not too close to the body
				if (spacecraft_distance_from_third_body_in_km > Universe->bodies[Universe->perturbation_menu[b]].radius + Universe->bodies[Universe->perturbation_menu[b]].minimum_safe_flyby_altitude)
				{
					Fgravity = -Universe->bodies[Universe->perturbation_menu[b]].mu / (spacecraft_distance_from_third_body_in_km*spacecraft_distance_from_third_body_in_km) * 1000.0;

					force_vector[0] += Fgravity * spacecraft_position_relative_to_third_body_in_km[0] / spacecraft_distance_from_third_body_in_km;
					force_vector[1] += Fgravity * spacecraft_position_relative_to_third_body_in_km[1] / spacecraft_distance_from_third_body_in_km;
					force_vector[2] += Fgravity * spacecraft_position_relative_to_third_body_in_km[2] / spacecraft_distance_from_third_body_in_km;

					if (generate_derivatives)
					{
						if (normalized_STMs)
						{
							double mu_3B_normalized = Universe->bodies[Universe->perturbation_menu[b]].mu / Universe->mu;
							double three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 = 3.0 * mu_3B_normalized / pow(spacecraft_distance_from_third_body_in_LU, 5.0);
							double mu3B_over_spacecraft_distance_from_third_body_in_LU3 = mu_3B_normalized / pow(spacecraft_distance_from_third_body_in_LU, 3.0);

							A(3, 0) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[0] * spacecraft_position_relative_to_third_body_in_LU[0] - mu3B_over_spacecraft_distance_from_third_body_in_LU3;
							A(3, 1) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[0] * spacecraft_position_relative_to_third_body_in_LU[1];
							A(3, 2) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[0] * spacecraft_position_relative_to_third_body_in_LU[2];
							A(4, 0) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[1] * spacecraft_position_relative_to_third_body_in_LU[0];
							A(4, 1) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[1] * spacecraft_position_relative_to_third_body_in_LU[1] - mu3B_over_spacecraft_distance_from_third_body_in_LU3;
							A(4, 2) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[1] * spacecraft_position_relative_to_third_body_in_LU[2];
							A(5, 0) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[2] * spacecraft_position_relative_to_third_body_in_LU[0];
							A(5, 1) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[2] * spacecraft_position_relative_to_third_body_in_LU[1];
							A(5, 2) += three_mu3B_over_spacecraft_distance_from_third_body_in_LU5 * spacecraft_position_relative_to_third_body_in_LU[2] * spacecraft_position_relative_to_third_body_in_LU[2] - mu3B_over_spacecraft_distance_from_third_body_in_LU3;
						}
						else
						{
							double mu_3B = Universe->bodies[Universe->perturbation_menu[b]].mu;
							double three_mu3B_over_spacecraft_distance_from_third_body_in_km5 = 3.0 * mu_3B / pow(spacecraft_distance_from_third_body_in_km, 5.0);
							double mu3B_over_spacecraft_distance_from_third_body_in_km3 = mu_3B / pow(spacecraft_distance_from_third_body_in_km, 3.0);

							A(3, 0) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[0] * spacecraft_position_relative_to_third_body_in_km[0] - mu3B_over_spacecraft_distance_from_third_body_in_km3;
							A(3, 1) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[0] * spacecraft_position_relative_to_third_body_in_km[1];
							A(3, 2) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[0] * spacecraft_position_relative_to_third_body_in_km[2];
							A(4, 0) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[1] * spacecraft_position_relative_to_third_body_in_km[0];
							A(4, 1) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[1] * spacecraft_position_relative_to_third_body_in_km[1] - mu3B_over_spacecraft_distance_from_third_body_in_km3;
							A(4, 2) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[1] * spacecraft_position_relative_to_third_body_in_km[2];
							A(5, 0) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[2] * spacecraft_position_relative_to_third_body_in_km[0];
							A(5, 1) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[2] * spacecraft_position_relative_to_third_body_in_km[1];
							A(5, 2) += three_mu3B_over_spacecraft_distance_from_third_body_in_km5 * spacecraft_position_relative_to_third_body_in_km[2] * spacecraft_position_relative_to_third_body_in_km[2] - mu3B_over_spacecraft_distance_from_third_body_in_km3;
						}
					}
				}
			}

			

			
		}

		return 0;
	}// end FBLT_force_model

}} //close namespace
