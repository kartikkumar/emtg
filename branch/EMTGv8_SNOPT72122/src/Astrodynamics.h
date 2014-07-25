//header file for astrodynamics functions

#include "missionoptions.h"
#include "universe.h"

#ifndef _ASTRODYNAMICS
#define _ASTRODYNAMICS

//#define _use_universal_kepler

namespace EMTG{ namespace Astrodynamics
{
	//******************************
	//insertion burn
	double insertion_burn(const double *Vin, const double *Vplanet, const double mu, const double R_SOI, const double destination_a, const double destination_e, double* v_inf_inbound);

	//******************************
	//unpowered flyby
	void unpowered_flyby (const double* v_initial, const double *planet_v, const double mu, //INPUT
	           const double rplanet, const double R_SOI, double rp, const double beta,//INPUT
	           double* v_final, double* outgoing_C3, double* flyby_delta, double *flyby_orbit_energy);//OUTPUT
	
	//******************************
	//powered flyby
	void powered_flyby (const double *v_initial, const double *v_final, const double *planet_v, const double mu, //INPUT
	           const double rplanet, const double R_SOI,//INPUT
	           double *dVmag, double *flyby_altitude, double* flyby_delta, double* outgoing_C3, double *flyby_orbit_energy); //output

	//functions to convert between inertial coordinates and classical orbit elements
	void COE2inertial(const double* E_COE, const double mu, double* state);
	void inertial2COE(const double* state, const double mu, double* E_COE);

	//functions to convert between classical orbit elements and modified equinoctial orbit elements
	void COE2MEE(const double* E_COE, double* E_MEE);
	void MEE2COE(const double* E_MEE, double* E_COE);

	//******************************
	//thruster code
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
								double* dPdt);

	//******************************
	//force model code
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
					double* dSRPdr,
					vector<double>& dagravdRvec,
					vector<double>& dagravdtvec);

	//******************************
	//launch vehicle code
	void find_mass_to_orbit(double C3, double DLA, int LV_type, double* mass, double* dmdC3, missionoptions* options);

	//******************************
	//spiral approximation code
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
								);

	void Edelbaum_spiral	(const double& mass_before_spiral,
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
							);
}}

#endif //_ASTRODYNAMICS
