//Battin spiral code

#include "missionoptions.h"
#include "Astrodynamics.h"
#include <math.h>

namespace EMTG { namespace Astrodynamics {
	void Battin_spiral	(const double& mass_before_spiral,
								const double& r0,
								const double& R0_planet,
								const double& time_since_launch,
								double& Isp,
								double& spiral_thrust,
								double& spiral_mdot,
								int& spiral_number_of_engines,
								double& spiral_power,
								double& spiral_active_power,
								const double& mu,
								const bool& generate_derivatives,
								missionoptions* options,
								double* t_spiral,
								double* mass_after_spiral,
								double* dv_spiral,
								double* dm_after_dm_before,
								double* dt_spiral_dm_before,
								double* dm_after_dIsp,
								double* dt_spiral_dIsp
								)
	{
		//compute the engine parameters at the beginning of the spiral, to be held constant
		double dTdP, dmdotdP, dTdIsp, dmdotdIsp, dPdr, dPdt;

		find_engine_parameters(	options,
								R0_planet,
								time_since_launch,
								&spiral_thrust,
								&spiral_mdot,
								&Isp,
								&spiral_power,
								&spiral_active_power,
								&spiral_number_of_engines,
								generate_derivatives,
								&dTdP,
								&dmdotdP,
								&dTdIsp,
								&dmdotdIsp,
								&dPdr,
								&dPdt);



		//compute tangential acceleration at the beginning of the spiral, to be held constant
		double a_t = spiral_thrust * options->engine_duty_cycle / (1000*mass_before_spiral);

		//compute circular orbit velocity
		double v_c = sqrt(mu / r0);

		//compute time for the spiral via Battin equation 8.110, p417
		double v_c2 = v_c * v_c;
		double v_c4 = v_c2 * v_c2;
		*t_spiral = -(v_c * (pow(( 20 * a_t*a_t * r0*r0)/v_c4, (1.0/8.0)) - 1))/a_t;

		//compute delta-v
		*dv_spiral = *t_spiral * a_t;

		//compute mass after spiral
		double expfun = exp(-*dv_spiral * 1000 / (Isp * options->g0));
		*mass_after_spiral = mass_before_spiral * expfun;

		//compute derivatives if desired
		if (generate_derivatives)
		{
			//derivative of spiral time with respect to initial mass
			double m_bs2 = mass_before_spiral * mass_before_spiral;
			double v_c3 = v_c * v_c2;
			double v_c7 = v_c3 * v_c4;
			double r02 = r0*r0;
			double r04 = r02*r02;
			double a_t2 = a_t*a_t;

			
			//derivative of spiral time with respect to acceleration
			double dt_spiral_da_t = 1.090661575 * a_t2 * r04 / (v_c7 * pow(a_t2 * r02 / v_c4, 15.0/8.0)) - v_c/a_t2;
			
			//derivative of acceleration with respect to mass before spiral
			double da_t_dm_bs = -a_t / mass_before_spiral;

			//derivative of spiral time with respect to mass before spiral
			*dt_spiral_dm_before = dt_spiral_da_t * da_t_dm_bs;
			
			//derivative of spiral delta-v with respect to mass before spiral
			double ddv_spiral_dm_bs = a_t * *dt_spiral_dm_before + *t_spiral * da_t_dm_bs;

			//derivative of mass after spiral with respect to mass before spiral
			*dm_after_dm_before = expfun * (1.0 - 1.0 / (Isp * options->g0) * ddv_spiral_dm_bs);

			//derivative of acceleration with respect to Isp
			double da_t_dIsp = options->engine_duty_cycle / 1000.0 / mass_before_spiral;

			//derivative of spiral time with respect to Isp
			*dt_spiral_dIsp = dt_spiral_da_t * da_t_dIsp;

			//derivative of spiral delta-v with respect to Isp
			double ddv_spiral_dIsp = a_t * *dt_spiral_dIsp + *t_spiral * da_t_dIsp;

			//derivative of mass after spiral with respect to Isp
			*dm_after_dIsp = *mass_after_spiral * 1.0 / (Isp * options->g0) * (*dv_spiral/Isp - ddv_spiral_dIsp);
		}
	}//close function Battin_spiral

	void Edelbaum_spiral(const double& mass_before_spiral,
						const double& r0,
						const double& R0_planet,
						const double& r_SOI,
						const double& time_since_launch,
						double& Isp,
						double& spiral_thrust,
						double& spiral_mdot,
						int& spiral_number_of_engines,
						double& spiral_power,
						double& spiral_active_power,
						const double& mu,
						const bool& generate_derivatives,
						missionoptions* options,
						double* t_spiral,
						double* mass_after_spiral,
						double* dv_spiral,
						double* dm_after_dm_before,
						double* dt_spiral_dm_before,
						double* dm_after_dIsp,
						double* dt_spiral_dIsp
						)
	{
		//compute the engine parameters at the beginning of the spiral, to be held constant
		double dTdP, dmdotdP, dTdIsp, dmdotdIsp, dPdr, dPdt;

		find_engine_parameters(	options,
								R0_planet,
								time_since_launch,
								&spiral_thrust,
								&spiral_mdot,
								&Isp,
								&spiral_power,
								&spiral_active_power,
								&spiral_number_of_engines,
								generate_derivatives,
								&dTdP,
								&dmdotdP,
								&dTdIsp,
								&dmdotdIsp,
								&dPdr,
								&dPdt);

		//compute delta-v
		*dv_spiral = sqrt(mu/r0) - sqrt(mu/r_SOI);

		//compute mass after spiral
		double expfun = exp(-*dv_spiral * 1000 / (Isp * options->g0));
		*mass_after_spiral = mass_before_spiral * expfun;

		//compute time for spiral
		*t_spiral = (mass_before_spiral - *mass_after_spiral) / spiral_mdot;

		//compute derivatives if desired
		if (generate_derivatives)
		{
			//propulsion system efficiency
			double P2 = spiral_active_power * spiral_active_power;
			double g02 = options->g0 * options->g0;
			double Isp2 = Isp * Isp;
			double eta = (Isp*options->g0*options->g0*spiral_mdot)/(2*P2);

			//derivative of spiral time with respect to mass before spiral
			*dt_spiral_dm_before = -(Isp2*g02 * (expfun - 1))/(2*P2*eta);
			
			//derivative of mass after spiral with respect to mass before spiral
			*dm_after_dm_before = expfun;

			//derivative of spiral time with respect to Isp
			*dt_spiral_dIsp = (Isp*g02*mass_before_spiral*(1 - expfun))/(P2*eta) - (options->g0*mass_before_spiral*expfun*(*dv_spiral))/(2*P2*eta);

			//derivative of mass after spiral with respect to Isp
			*dm_after_dIsp = (mass_before_spiral*expfun*(*dv_spiral))/(Isp2*options->g0);
		}
	}//close function Edelbaum_spiral



}}//close namespace