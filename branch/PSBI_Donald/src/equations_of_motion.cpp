//differential equations of motion for EMTG v8

#include <cmath>

#include "missionoptions.h"
#include "EMTG_math.h"
#include "Astrodynamics.h"
#include "equations_of_motion.h"
#include "universe.h"
#include "SpiceUsr.h"


namespace EMTG { namespace Astrodynamics {namespace EOM
{	
	//equations of motion for an object moving in the heliocentric inertial frame with a thrust term
	void EOM_inertial_continuous_thrust(double* x,
										const double& t,
										const double& t0,
										double* u,
										double* f,
										double* thrust,
										double* mdot,
										double* Isp,
										double* power,
										double* active_power,
										int* number_of_active_engines,
										int &STMrows,
										int &STMcolumns,
										void* optionsvoidpointer, 
										void* Universepointer, 
										void* ControllerPointer)
	{
		missionoptions* options = (missionoptions*) optionsvoidpointer;
		EMTG::Astrodynamics::universe* Universe = (EMTG::Astrodynamics::universe*) Universepointer;
		
		double ForceVector[3];
		double spacecraft_state[7];
		double epoch = t * Universe->TU;
		double launch_epoch = t0;

		double dTdP,dmdotdP,dTdIsp,dmdotdIsp,dPdr,dPdt,dFSRPdr;
		
        static vector<double> dagravdRvec(3), dagravdtvec(3), central_body_state_mks((options->derivative_type > 2) ? 12 : 6);

		double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); //magnitude of position vector
		r = fabs(r) < EMTG::math::SMALL ? EMTG::math::sgn(r) * EMTG::math::SMALL : r;

		double mass = x[6] * options->maximum_mass;
		mass = mass < EMTG::math::SMALL ? EMTG::math::SMALL : mass;

		//determine the forces acting on the object 
		//and calculate the state propagation matrix A

		static EMTG::math::Matrix<double> A (STMrows, STMcolumns, 0.0);

		//The s/c state must be converted to MKS in order to interface with the force and engine models
		for (int k = 0; k < 3; ++k)
		{
			spacecraft_state[k] = x[k] * Universe->LU;
			spacecraft_state[k+3] = x[k+3] * Universe->LU / Universe->TU;
		}
		spacecraft_state[6] = mass;

		if (options->engine_type == 1)
			*Isp = options->IspLT;
		else if (options->engine_type == 2)
			*power = options->power_at_1_AU;

		if (options->mission_type == 3)
		{
			EMTG::Astrodynamics::FBLT_force_model(options,
				Universe,
				spacecraft_state,
				&epoch,
				&launch_epoch,
				u,
				thrust,
				mdot,
				Isp,
				power,
				active_power,
				number_of_active_engines,
				ForceVector,
				true,
				&dTdP,
				&dmdotdP,
				&dTdIsp,
				&dmdotdIsp,
				&dPdr,
				&dPdt,
				&dFSRPdr,
				A,
				dagravdRvec,
				dagravdtvec,
                central_body_state_mks);
		}

		else
		{
			EMTG::Astrodynamics::force_model(options,
				Universe,
				spacecraft_state,
				&epoch,
				&launch_epoch,
				u,
				thrust,
				mdot,
				Isp,
				power,
				active_power,
				number_of_active_engines,
				ForceVector,
				false,
				&dTdP,
				&dmdotdP,
				&dTdIsp,
				&dmdotdIsp,
				&dPdr,
				&dPdt,
				&dFSRPdr,
				dagravdRvec,
                dagravdtvec,
                central_body_state_mks);
		}

		//for (int k = 0; k < 3; ++k)
			//ForceVector[k] /= ( (Universe->LU / (Universe->TU * Universe->TU)) );

		//convert mdot from kg/s to SpacecraftMassUnits/TU
		//double working_mdot = *mdot * Universe->TU / options->maximum_mass;


		double mu_normalized = 1.0;

		//EOM

		//x
		f[0] = x[3];

		//y
		f[1] = x[4];

		//z
		f[2] = x[5];

		//xdot
		f[3] = -mu_normalized*x[0] / (r*r*r) + ForceVector[0] / x[6];

		//ydot
		f[4] = -mu_normalized*x[1] / (r*r*r) + ForceVector[1] / x[6];

		//zdot
		f[5] = -mu_normalized*x[2] / (r*r*r) + ForceVector[2] / x[6];

		//mass
		f[6] = -EMTG::math::norm(u, 3) * (*mdot) * options->engine_duty_cycle;



		//control vector time rates of change (these are decision variables and are not impacted by dynamics and are constant across an FBLT time step)
		//these are just a place holders in the STM and will never have a diffeqs associated with them
		f[7] = 0.0000000000000000;
		f[8] = 0.0000000000000000;
		f[9] = 0.0000000000000000;

		//TOF time rates of change (this is also a decision variable; it is not impacted by dynamics and is constant across an FBLT time step/phase)
		//this is just a place holder in the STM and will never have a diffeq associated with it
		f[10] = 0.0000000000000000;

		//*******************
		//
		//Phi dot calculation
		//
		//*******************

		//Form the STM
		//The STM is comprised of state variables

		EMTG::math::Matrix<double> STM(STMrows, STMcolumns, 0.0);
		EMTG::math::Matrix<double> STMdot(STMrows, STMcolumns, 0.0);

		int xcount = 11;

		for (size_t i = 0; i < STMrows; ++i)
		{
			for (size_t j = 0; j < STMcolumns; ++j)
			{
				STM(i, j) = x[xcount];
				++xcount;
			}
		}

		//differential equation for STM creation
		STMdot = A * STM;

		//Package the Phidot entries into the gradient vector f behind the spacecraft's state variables as well as the control and TOF gradients which are always zero
		int fcount = 11; //STM entry diffeq's are placed behind the augmented state vector
		for (size_t i = 0; i < STMrows; ++i)
		{
			for (size_t j = 0; j < STMcolumns; ++j)
			{
				f[fcount] = STMdot(i, j);
				++fcount;
			}
		}
		
		
	}

	//equations of motion for an object moving in the heliocentric inertial frame with a thrust term using the Sundman transformation
	void EOM_inertial_continuous_thrust_sundman(double* x,
												const double& t,
												const double& t0,
												double* u,
												double* f,
												double* thrust,
												double* mdot,
												double* Isp,
												double* power,
												double* active_power,
												int* number_of_active_engines,
												void* optionsvoidpointer, 
												void* Universepointer, 
												void* ControllerPointer)
	{
		missionoptions* options = (missionoptions*) optionsvoidpointer;
		EMTG::Astrodynamics::universe* Universe = (EMTG::Astrodynamics::universe*) Universepointer;
		double mu_normalized = 1.0;
		
		double ForceVector[3];
		double epoch = x[7]*Universe->TU;
		double launch_epoch = t0;

		double dTdP,dmdotdP,dTdIsp,dmdotdIsp,dPdr,dPdt,dFSRPdr;
		
        static vector<double> dagravdRvec(3), dagravdtvec(3), central_body_state_mks((options->derivative_type > 2) ? 12 : 6);

		double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); //magnitude of position vector
		r = fabs(r) < EMTG::math::SMALL ? EMTG::math::sgn(r) * EMTG::math::SMALL : r;

		double mass = x[6] * options->maximum_mass;
		mass = mass < EMTG::math::SMALL ? EMTG::math::SMALL : mass;

		//determine the forces acting on the object

		EMTG::Astrodynamics::force_model(options,
										Universe,
										x,
										&epoch,
										&launch_epoch,
										u,
										thrust,
										mdot,
										Isp,
										power,
										active_power,
										number_of_active_engines,
										ForceVector,
										false,
										&dTdP,
										&dmdotdP,
										&dTdIsp,
										&dmdotdIsp,
										&dPdr,
										&dPdt,
										&dFSRPdr,
										dagravdRvec,
                                        dagravdtvec,
                                        central_body_state_mks);

		for (int k = 0; k < 3; ++k)
			ForceVector[k] /= (Universe->LU / (Universe->TU * Universe->TU));

		//convert mdot from kg/s to SpacecraftMassUnits/TU
		double working_mdot = *mdot * Universe->TU / options->maximum_mass;

		//EOM

		//x
		f[0] = r * x[3];

		//y
		f[1] = r * x[4];

		//z
		f[2] = r * x[5];

		//xdot
        f[3] = -mu_normalized*x[0] / (r*r) + r * ForceVector[0] / mass;

		//ydot
        f[4] = -mu_normalized*x[1] / (r*r) + r * ForceVector[1] / mass;

		//zdot
        f[5] = -mu_normalized*x[2] / (r*r) + r * ForceVector[2] / mass;

		//mass
		f[6] = -r * EMTG::math::norm(u, 3) * working_mdot;

		//time
		f[7] = r;
	}

}}}