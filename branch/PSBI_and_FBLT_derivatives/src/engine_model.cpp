//code for engine model
//all data is public domain except for VASIMR and Xenon Hall thruster models which are kept in a separate sourcefile and are NOT TO BE DISTRIBUTED

#include "Astrodynamics.h"
#include "missionoptions.h"
#include "EMTG_math.h"

//if the proprietary switch is turned on, include the proprietary functions. This header file is NOT distributed with the open-source version of EMTG
#ifdef _EMTG_proprietary
#include "proprietary_functions.h"
#endif

namespace EMTG { namespace Astrodynamics {

int find_engine_parameters(	EMTG::missionoptions* options,
								const double& r,
								const double& total_flight_time,
								double* max_thrust,
								double* max_mass_flow_rate,
								double* Isp,
								double* power,
								double* active_power,
								int* number_of_active_engines,
								const bool& generate_derivatives,
								double* dTdP,
								double* dmdotdP,
								double* dTdIsp,
								double* dmdotdIsp,
								double* dPdr,
								double* dPdt)
{
	//placeholders for engine performance coefficients
	double at, bt, ct, dt, et, ht, gt, af, bf, cf, df, ef, hf, gf, minP, maxP, power_penalty = 1.0;

	if (options->engine_type == 0) //fixed thrust, Isp
	{
		*Isp = options->IspLT;
		*max_thrust = options->Thrust;
		*max_mass_flow_rate = *max_thrust / *Isp / options->g0;
		*active_power = 0.0;

		if (generate_derivatives)
		{
			*dTdP = 0.0;
			*dmdotdP = 0.0;
			*dTdIsp = 0.0;
			*dmdotdIsp = 0.0;
			*dPdr = 0.0;
			*dPdt = 0.0;
		}
	}
	else if (options->engine_type == 1) //power is independent variable, efficiency and Isp are specified
	{
		//extract Isp from the options structure
		*Isp = options->IspLT;

		//compute thrust
		*dTdP = 2000.0 * (options->user_defined_engine_efficiency / (*Isp * options->g0));
		*max_thrust = *dTdP * options->power_at_1_AU;
		*active_power = options->power_at_1_AU;

		//compute mass flow rate
		*dmdotdP = *dTdP /  *Isp / options->g0;
		*max_mass_flow_rate = *dmdotdP * options->power_at_1_AU;

		if (generate_derivatives)
		{
			*dTdIsp = 0.0;
			*dmdotdIsp = 0.0;
			*dPdr = 0.0;
			*dPdt = 0.0;
		}
	}
	else //all engine types where the available power is not fixed
	{
		//how much power is produced?
		double input_power, spacecraft_power;
		double r2 = r*r;
		if (options->power_source_type == 0) //solar power
		{
			double g0 = options->solar_power_gamma[0];
			double g1 = options->solar_power_gamma[1];
			double g2 = options->solar_power_gamma[2];
			double g3 = options->solar_power_gamma[3];
			double g4 = options->solar_power_gamma[4];
			
			input_power = options->power_at_1_AU / r2 * ( (g0 + g1/r + g2/r2) / (1.0 + g3 * r + g4 * r2) );

			if (generate_derivatives)
			{
				double r3 = r2*r;
				double r4 = r3*r;
				double r5 = r4*r;
				*dPdr = -options->power_at_1_AU * (4*g2 + 3*g1*r + g3*(3*g0*r3 + 4*g1*r2 + 5*g2*r) + 2*g0*r2 + g4*(4*g0*r4 + 5*g1*r3 + 6*g2*r2))/(r5*(g4*r2 + g3*r + 1)*(g4*r2 + g3*r + 1));

			}
		}
		else //constant power
		{
			input_power = options->power_at_1_AU;
			*dPdr = 0.0;
		}

		//if applicable, model the decay rate of the power system
		if (options->power_decay_rate > 1.0e-5)
		{
			double decay_coeff = pow( (1 - options->power_decay_rate), total_flight_time / (365.25 * 86400.0) );
			input_power *= decay_coeff;

			if (generate_derivatives)
			{
				*dPdr *= decay_coeff;
				*dPdt = log(1 - options->power_decay_rate) * input_power / (365.25 * 86400.0);
			}
		}
		else
			*dPdt = 0.0;

		//how much power is used by the spacecraft?
		if (options->spacecraft_power_model_type == 0)
		{
			spacecraft_power = options->spacecraft_power_coefficients[0] + options->spacecraft_power_coefficients[1]/r + options->spacecraft_power_coefficients[2]/r2;
			*dPdr -= -options->spacecraft_power_coefficients[1]/r2 - options->spacecraft_power_coefficients[2]/(r2*r);
		}
		else
		{
			if (input_power > options->spacecraft_power_coefficients[0])
			{
				spacecraft_power = options->spacecraft_power_coefficients[0];
				//no contribution to dPdr
			}
			else
			{
				spacecraft_power = options->spacecraft_power_coefficients[0] + options->spacecraft_power_coefficients[1] * (options->spacecraft_power_coefficients[2] - input_power);
				*dPdr -= -options->spacecraft_power_coefficients[1] * (*dPdr);
			}
		}

		//how much power is available to the engines?
		*power = (1.0 - options->power_margin) * (input_power - spacecraft_power);
		if (generate_derivatives)
			*dPdr *= (1.0 - options->power_margin);
		*power = *power < math::SMALL ? math::SMALL : *power;

		//now, what subtype of engine is this?
		if (options->engine_type == 2 || options->engine_type == 4) //Isp is independent variable, power uses power model, and efficiency is fixed
		{
			if (*power < options->engine_input_power_bounds[0])
			{
				//compute thrust
				*dTdP = 0;
				*max_thrust = 0;

				//compute mass flow rate
				*dmdotdP = 0;
				*max_mass_flow_rate = 0;

				*active_power = 0.0;
				*power = input_power;

				if (generate_derivatives)
				{
					*dTdIsp = 0;
					*dmdotdIsp = 0;
				}
			}
			else
			{
				double effective_power = *power > options->engine_input_power_bounds[1] ? options->engine_input_power_bounds[1] : *power;

				//compute thrust
				*dTdP = 2000 * (options->user_defined_engine_efficiency / (*Isp * options->g0));
				*max_thrust = *dTdP * effective_power;

				//compute mass flow rate
				*dmdotdP = *dTdP /  *Isp / options->g0;
				*max_mass_flow_rate = *dmdotdP * effective_power;

				*active_power = effective_power;
				*power = input_power;

				if (generate_derivatives)
				{
					*dTdIsp = -1.0 / (*Isp) * (*max_thrust);
					*dmdotdIsp = -2.0 / (*Isp) * (*max_mass_flow_rate);
				}
			}
		}
		else if (options->engine_type == 3) //constant Isp, efficiency, compute thrust based on available power
		{
			if (*power < options->engine_input_power_bounds[0])
			{
				//compute thrust
				*dTdP = 0;
				*max_thrust = 0;

				//compute mass flow rate
				*dmdotdP = 0;
				*max_mass_flow_rate = 0;

				*active_power = 0.0;
				*power = input_power;

				if (generate_derivatives)
				{
					*dTdIsp = 0;
					*dmdotdIsp = 0;
				}
			}
			else
			{
				double effective_power = *power > options->engine_input_power_bounds[1] ? options->engine_input_power_bounds[1] : *power;

				//extract Isp from the options structure
				*Isp = options->IspLT;

				//compute thrust
				*dTdP = 2000 * (options->user_defined_engine_efficiency / (*Isp * options->g0));
				*max_thrust = *dTdP * effective_power;

				//compute mass flow rate
				*dmdotdP = *dTdP /  *Isp / options->g0;
				*max_mass_flow_rate = *dmdotdP * effective_power;

				*active_power = effective_power;
				*power = input_power;
			}
		}
		else if (options->engine_type == 12) //VASIMR analytical model
			{
#ifdef _EMTG_proprietary
			EMTG::Proprietary::VASIMR_model(power, Isp, max_thrust, max_mass_flow_rate, dTdP, dmdotdP, options);
#else
				cout << "VASIMR model not included in open-source version" << endl;
				throw 1711;
#endif
			}
		else if (options->engine_type == 13) //Xenon hall thruster analytical model
			{
#ifdef _EMTG_proprietary
			EMTG::Proprietary::HallThrusterXenon_model(power, Isp, max_thrust, max_mass_flow_rate, dTdP, dmdotdP, options);
#else
				cout << "Xenon Hall thruster model not included in open-source version" << endl;
				throw 1711;
#endif
			}
		
		else//custom engine polynomials or standard engine
		{
			if (options->engine_type == 5) 
			{
				at = options->engine_input_thrust_coefficients[0];
				bt = options->engine_input_thrust_coefficients[1];
				ct = options->engine_input_thrust_coefficients[2];
				dt = options->engine_input_thrust_coefficients[3];
				et = options->engine_input_thrust_coefficients[4];
				gt = options->engine_input_thrust_coefficients[5];
				ht = options->engine_input_thrust_coefficients[6];

				af = options->engine_input_mass_flow_rate_coefficients[0];
				bf = options->engine_input_mass_flow_rate_coefficients[1];
				cf = options->engine_input_mass_flow_rate_coefficients[2];
				df = options->engine_input_mass_flow_rate_coefficients[3];
				ef = options->engine_input_mass_flow_rate_coefficients[4];
				gf = options->engine_input_mass_flow_rate_coefficients[5];
				hf = options->engine_input_mass_flow_rate_coefficients[6];

				minP = options->engine_input_power_bounds[0];
				maxP = options->engine_input_power_bounds[1];
			}
			else if (options->engine_type == 21) //13 kW STMD Hall high-Isp
			{
#ifdef _EMTG_proprietary
				EMTG::Proprietary::STMDHall13kWHIsp_model(&at,
														&bt,
														&ct,
														&dt,
														&et,
														&gt,
														&ht,
														&af,
														&bf,
														&cf,
														&df,
														&ef,
														&gf,
														&hf,
														&minP,
														&maxP);
#else
				cout << "STMD 13 kW Hall thruster model not included in open-source version" << endl;
				throw 1711;
#endif
			}
			else if (options->engine_type == 22) //13 kW STMD Hall high-thrust
			{
#ifdef _EMTG_proprietary
				EMTG::Proprietary::STMDHall13kWHthrust_model(&at,
															&bt,
															&ct,
															&dt,
															&et,
															&gt,
															&ht,
															&af,
															&bf,
															&cf,
															&df,
															&ef,
															&gf,
															&hf,
															&minP,
															&maxP);
#else
				cout << "STMD 13 kW Hall thruster model not included in open-source version" << endl;
				throw 1711;
#endif
			}
			else if (options->engine_type >= 26 && options->engine_type <= 28) //13 kW STMD Hall 9-8-2014
			{
#ifdef _EMTG_proprietary
				EMTG::Proprietary::STMD_13kW_Hall_10_1_2014(options->engine_type - 26,
															&at,
															&bt,
															&ct,
															&dt,
															&et,
															&gt,
															&ht,
															&af,
															&bf,
															&cf,
															&df,
															&ef,
															&gf,
															&hf,
															&minP,
															&maxP);
#else
				cout << "STMD 13 kW Hall thruster model not included in open-source version" << endl;
				throw 1711;
#endif
			}
			else //standard engine from list - these models are available in the public literature
			{
				static double min_power[] = {0.525, 0.436, 0.302, 0.302, 0.302, 1.252, 0.638, 0.638, 1.15, 16.2, 7.0, 5.0, 0.354, 0.64, 0.64, 0.64};
				static double max_power[] = {2.6, 5.03, 4.839, 4.839, 4.839, 7.455, 7.266, 7.266, 4.91, 23.04, 12.0, 17.5, 3.821, 7.36, 7.36, 7.36};
                int k = options->engine_type < 12 ? options->engine_type - 6 : options->engine_type - 8;
				k = options->engine_type > 20 ? k - 2 : k;
				k = options->engine_type > 25 ? k - 3 : k;
                
	
				//1: NSTAR, 2: XIPS-25, 3: BPT-4000 High-Isp, 4: BPT-4000 High-Thrust, 5: BPT-4000 Ex-High-Isp, 6: NEXT high-Isp old, 7: NEXT high-Isp v10, 8: NEXT high-thrust v10, 9: BPT-4000 MALTO, 10: NEXIS, 11: H6MS, 12: BHT20K, 13: Aerojet HiVHAC EM
	
				//first, the coefficients for thrust
				static double Ht[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
				static double Gt[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
				static double Et[] = { 5.145602, 0.0367, -0.095437, 0.173870, 1.174296, 0.02334, -0.19082, 0.09591, -0.6574, 0.0, 0.0, 0.0, 0.0, 0.101855017, -0.111563126, 0.085120672 };
				static double Dt[] = { -36.720293, -0.4966, 1.637023, -1.150940, -10.102479, -0.6815, 2.96519, -1.98537, 6.2683, 0.0, 0.0, 0.0, 0.0, -2.04053417, 1.72548416, -1.42659172 };
				static double Ct[] = { 90.486509, 1.4111, -9.517167, -2.118891, 19.422224, 6.882, -14.41789, 11.47980, -14.2820, 0.119, -0.775176, 0.035143, -4.2897, 11.4181412, -7.91621814, 5.17797704 };
				static double Bt[] = { -51.694393, 35.3591, 72.030104, 77.342132, 47.927765, 3.7746, 54.05382, 15.06977, 37.751, 5.042, 58.399388, 51.393286, 54.696, 16.0989424, 40.543251, 37.1873936 };
				static double At[] = { 26.337459, -0.3984, -7.181341, -8.597025, 1.454064, 36.467, -1.92224e-6, 14.51552, 49.466, 293.9, 64.444882, -10.812857, 1.0202, 1.19388817E+01, 3.68945763, -0.804281458 };
	
				//next, the coefficients for mass flow rate
				static double Hf[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
				static double Gf[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
				static double Ef[] = { 0.36985, 0.0091, -0.008432, -0.011949, 0.086106, 0.002892, -0.004776, 0.01492, 0.0025, 0.0, 0.0, 0.0, -0.0271, 0.011021367, -0.002913991, 0.009019519 };
				static double Df[] = { -2.5372, -0.1219, 0.148511, 0.235144, -0.727280, -0.0718, 0.05717, -0.27539, -0.2629, 0.0, 0.0, 0.0, 0.2499, -0.207253445, 0.029887398, -0.138009326 };
				static double Cf[] = { 6.2539, 0.5402, -0.802790, -1.632373, 1.328508, 0.6470, -0.09956, 1.60966, 2.712, 0.01084, -0.019616, 0.0012, -0.9465, 1.21670237, 0.027771576, 0.535910391 };
				static double Bf[] = { -5.3568, 0.0426, 3.743362, 6.847936, 1.998082, -1.7266, 0.03211, -2.53056, -7.232, -0.49113, 1.5188, 2.0157, 2.4711, -1.71102132, -0.180919262, 0.534442545 };
				static double Af[] = { 2.5060, 0.8211, 1.244345, 0.352444, 1.653105, 3.630, 2.13781, 3.22089, 12.204, 11.979, 1.4558, 0.9977, 1.251, 2.75956482, 2.22052155, 1.40535083 };

				minP = min_power[k];
				maxP = max_power[k];
				at = At[k];
				bt = Bt[k];
				ct = Ct[k];
				dt = Dt[k];
				et = Et[k];
				gt = Gt[k];
				ht = Ht[k];

				af = Af[k];
				bf = Bf[k];
				cf = Cf[k];
				df = Df[k];
				ef = Ef[k];
				gf = Gf[k];
				hf = Hf[k];
			}

			//how many engines will we operate and at what power level?
			double P_eff = 0.0;
			double P_eff2, P_eff3, P_eff4, P_eff5, P_eff6;
			*number_of_active_engines = 0;
			if (*power >= minP)
			{
				if (options->throttle_logic_mode == 0)
				{ //maximum power use
					*number_of_active_engines = options->number_of_engines;
				
					for (int n = options->number_of_engines; n > 0; --n)
					{
						if (*power > n*minP && *power > (n-1)*maxP)
						{
							P_eff = min(*power / n, maxP);
							*number_of_active_engines = n;

							break;
						}
					}
				}
				else if (options->throttle_logic_mode == 1)
				{//maximum thrust
					//first, how many thrusters can be fired at once?
					int max_thrusters = min( (int)(*power / minP), options->number_of_engines );
					double BestThrust = 0.0;

					//now compute the total thrust for each number of thrusters from n to 1
					//assume equal distribution of power among the active thrusters
					for (int n = max_thrusters; n > 0; --n)
					{
						if (*power > n*minP)
						{
							P_eff = min(*power / n, maxP);
							P_eff2 = P_eff*P_eff;
							P_eff3 = P_eff*P_eff2;
							P_eff4 = P_eff*P_eff3;
							P_eff5 = P_eff*P_eff4;
							P_eff6 = P_eff*P_eff5;
							double CurrentThrust = n*(ht * P_eff6 + gt * P_eff5 + et*P_eff4 + dt*P_eff3 + ct*P_eff2 + bt*P_eff + at);
							if (CurrentThrust > BestThrust)
							{
								BestThrust = CurrentThrust;
								*number_of_active_engines = n;
							}
						}
					}
					P_eff = min(*power / *number_of_active_engines, maxP);
				}
				else if (options->throttle_logic_mode == 2)
				{//maximum Isp
					//first, how many thrusters can be fired at once?
					int max_thrusters = min( (int)(*power / minP), options->number_of_engines );
					double BestIsp = 0.0;

					//now compute the total flow rate for each number of thrusters from n to 1
					//assume equal distribution of power among the active thrusters
					for (int n = max_thrusters; n > 0; --n)
					{
						if (*power > n*minP)
						{
							P_eff = min(*power / n, maxP);
							P_eff2 = P_eff*P_eff;
							P_eff3 = P_eff*P_eff2;
							P_eff4 = P_eff*P_eff3;
							P_eff5 = P_eff*P_eff4;
							P_eff6 = P_eff*P_eff5;
							double CurrentFlowRate = n*(hf * P_eff6 + gf * P_eff5 + ef*P_eff4 + df*P_eff3 + cf*P_eff2 + bf*P_eff + af);
							double CurrentThrust = n*(ht * P_eff6 + gt * P_eff5 + et*P_eff4 + dt*P_eff3 + ct*P_eff2 + bt*P_eff + at);
							double CurrentIsp = CurrentThrust / CurrentFlowRate / options->g0 * 1000.0;
							if (CurrentIsp > BestIsp)
							{
								BestIsp = CurrentIsp;
								*number_of_active_engines = n;
							}
						}
					}
					P_eff = min(*power / *number_of_active_engines, maxP);
				}
				else if (options->throttle_logic_mode == 3)
				{//maximum efficiency
					//first, how many thrusters can be fired at once?
					int max_thrusters = min( (int)(*power / minP), options->number_of_engines );
					double BestEfficiency = 0.0;

					//now compute the total flow rate for each number of thrusters from n to 1
					//assume equal distribution of power among the active thrusters
					for (int n = max_thrusters; n > 0; --n)
					{
						if (*power > n*minP)
						{
							P_eff = min(*power / n, maxP);
							P_eff2 = P_eff*P_eff;
							P_eff3 = P_eff*P_eff2;
							P_eff4 = P_eff*P_eff3;
							P_eff5 = P_eff*P_eff4;
							P_eff6 = P_eff*P_eff5;
							double CurrentFlowRate = n*(hf * P_eff6 + gf * P_eff5 + ef*P_eff4 + df*P_eff3 + cf*P_eff2 + bf*P_eff + af);
							double CurrentThrust = n*(ht * P_eff6 + gt * P_eff5 + et*P_eff4 + dt*P_eff3 + ct*P_eff2 + bt*P_eff + at);
							double CurrentIsp = CurrentThrust / CurrentFlowRate / options->g0 * 1000.0;
							double CurrentEfficiency = CurrentThrust * CurrentIsp * options->g0 / (2000 * P_eff);
							if (CurrentEfficiency > BestEfficiency)
							{
								BestEfficiency = CurrentEfficiency;
								*number_of_active_engines = n;
							}
						}
					}
					P_eff = min(*power / *number_of_active_engines, maxP);
				}
			}
			else
				*number_of_active_engines = 0;

			//now either power the engines, or, if insufficient power, switch off the engines
			double T, F;
			
			double dP = 0.0;
			double dX = 0.0;
			double ddXdP = 0.0;

			if (options->throttle_sharpness < 100.0)
			{
				double smoother_width = 1.0 / (options->throttle_sharpness + 0.01);
				if (*power - *number_of_active_engines*minP <= smoother_width)
				{
					dP = *power - *number_of_active_engines*minP - smoother_width;
					dX = dP*dP*dP*(dP*(dP*6 - 15) + 10);
					ddXdP = 30*dP*dP * (dP * (dP - 2) + 1);
				}
				else if (*power - (*number_of_active_engines - 1) * maxP <= smoother_width)
				{
					dP = *power - (*number_of_active_engines - 1) * maxP - smoother_width;
					dX = dP*dP*dP*(dP*(dP*6 - 15) + 10);
					ddXdP = 30*dP*dP * (dP * (dP - 2) + 1);
				}
				else if (*number_of_active_engines * maxP - P_eff <= smoother_width)
				{
					dP = *number_of_active_engines * maxP - P_eff - smoother_width;
					dX = -dP*dP*dP*(dP*(dP*6 - 15) + 10);
					ddXdP = -30*dP*dP * (dP * (dP - 2) + 1);
				}
				else if (*power - minP <= smoother_width)
				{
					dP = *power - minP - P_eff - smoother_width;
					dX = dP*dP*dP*(dP*(dP*6 - 15) + 10);
					ddXdP = 30*dP*dP * (dP * (dP - 2) + 1);
				}
			}

			//next compute Thrust in mN
			P_eff2 = P_eff*P_eff;
			P_eff3 = P_eff*P_eff2;
			P_eff4 = P_eff*P_eff3;
			P_eff5 = P_eff*P_eff4;
			P_eff6 = P_eff*P_eff5;
			T = (1.0 - dX) * (ht * P_eff6 + gt * P_eff5 + et*P_eff4 + dt*P_eff3 + ct*P_eff2 + bt*P_eff + at);
			
			//next compute mass flow rate in mg/s
			F = (1.0 - dX) * (hf * P_eff6 + gf * P_eff5 + ef*P_eff4 + df*P_eff3 + cf*P_eff2 + bf*P_eff + af);
			

			if (*power < minP)
			{
				power_penalty = exp(-10* minP / *power);
				T *= power_penalty;
				F *= power_penalty;
			}

			//return Thrust in N (convert from mN) and mass flow rate in kg/s (convert from mg/s)
			*max_thrust = 1.0e-3 * T * *number_of_active_engines;
			*max_mass_flow_rate = 1.0e-6 * F * *number_of_active_engines;
			*Isp = *max_thrust / *max_mass_flow_rate / options->g0;
			*active_power = *number_of_active_engines * P_eff;
			
			

			if (generate_derivatives)
			{
				*dTdP = ((1.0 - dX)*(6*ht*P_eff5 + 5*gt*P_eff4 + 4*et*P_eff3 + 3*dt*P_eff2 + 2*ct*P_eff + bt) + ddXdP*T) * 1.0e-3 * *number_of_active_engines;
				*dmdotdP = ((1.0 - dX)*(6*hf*P_eff5 + 5*gf*P_eff4 + 4*ef*P_eff3 + 3*df*P_eff2 + 2*cf*P_eff + bf) + ddXdP*F)  * 1.0e-6 * *number_of_active_engines;
				*dTdIsp = 0.0;
				*dmdotdIsp = 0.0;


				if (*power < minP)
				{
					double dpenalty_dP = (10 * minP * power_penalty )/ ((*power) * (*power));
					*dTdP = *dTdP*power_penalty + T*dpenalty_dP;
					*dmdotdP = *dmdotdP*power_penalty + F*dpenalty_dP;
				}

				if (*power > maxP * *number_of_active_engines)
				{
					*dTdP = 0.0;
					*dmdotdP = 0.0;
				}
			}

			*power = input_power;
		}
	}

	*max_thrust /= 1000.0; //N to kN conversion
	*dTdIsp /= 1000.0; //N to kN conversion
	*dTdP /= 1000.0; //N to kN conversion

	return 0;
}

}}