//differential equations of motion for EMTG v8


#include "missionoptions.h"
#include "EMTG_math.h"
#include "Astrodynamics.h"
#include "equations_of_motion.h"
#include "universe.h"
#include "SpiceUsr.h"

#include <math.h>

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
										void* optionsvoidpointer, 
										void* Universepointer, 
										void* ControllerPointer)
	{
		missionoptions* options = (missionoptions*) optionsvoidpointer;
		EMTG::Astrodynamics::universe* Universe = (EMTG::Astrodynamics::universe*) Universepointer;
		double mu = 1;
		
		double ForceVector[3];
		double spacecraft_state[6];
		double epoch = t / 86400 * Universe->TU;
		double launch_epoch = t0;

		double dTdP,dmdotdP,dTdIsp,dmdotdIsp,dPdr,dPdt,dFSRPdr;
		
		static vector<double> dagravdRvec(3), dagravdtvec(3);

		double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); //magnitude of position vector
		r = fabs(r) < EMTG::math::SMALL ? EMTG::math::sgn(r) * EMTG::math::SMALL : r;

		double mass = x[6] * options->maximum_mass;
		mass = mass < EMTG::math::SMALL ? EMTG::math::SMALL : mass;

		//determine the forces acting on the object
		for (int k = 0; k < 3; ++k)
		{
			spacecraft_state[k] = x[k] * Universe->LU;
			spacecraft_state[k+3] = x[k+3] * Universe->LU / Universe->TU;
		}

		if (options->engine_type == 1)
			*Isp = options->IspLT;
		else if (options->engine_type == 2)
			*power = options->power_at_1_AU;

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
										dagravdtvec);

		for (int k = 0; k < 3; ++k)
			ForceVector[k] /= ( 1000.0 * (Universe->LU / (Universe->TU * Universe->TU)) );

		//convert mdot from kg/s to SpacecraftMassUnits/TU
		*mdot *= Universe->TU / options->maximum_mass;

		//EOM

		//x
		f[0] = x[3];

		//y
		f[1] = x[4];

		//z
		f[2] = x[5];

		//xdot
		f[3] = -mu*x[0]/(r*r*r) + ForceVector[0] / mass;

		//ydot
		f[4] = -mu*x[1]/(r*r*r) + ForceVector[1] / mass;

		//zdot
		f[5] = -mu*x[2]/(r*r*r) + ForceVector[2] / mass;

		//mass
		f[6] = -EMTG::math::norm(u, 3) * *mdot * options->engine_duty_cycle;
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
		double mu = 1;
		
		double ForceVector[3];
		double epoch = x[7]*Universe->TU / 86400;
		double launch_epoch = t0;

		double dTdP,dmdotdP,dTdIsp,dmdotdIsp,dPdr,dPdt,dFSRPdr;
		
		static vector<double> dagravdRvec(3), dagravdtvec(3);

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
										dagravdtvec);

		for (int k = 0; k < 3; ++k)
			ForceVector[k] /= 1000.0 * (Universe->LU / (Universe->TU * Universe->TU));

		//convert mdot from kg/s to SpacecraftMassUnits/TU
		*mdot *= Universe->TU / options->maximum_mass;

		//EOM

		//x
		f[0] = r * x[3];

		//y
		f[1] = r * x[4];

		//z
		f[2] = r * x[5];

		//xdot
		f[3] = -mu*x[0]/(r*r) + r * ForceVector[0] / mass;

		//ydot
		f[4] = -mu*x[1]/(r*r) + r * ForceVector[1] / mass;

		//zdot
		f[5] = -mu*x[2]/(r*r) + r * ForceVector[2] / mass;

		//mass
		f[6] = -r * EMTG::math::norm(u, 3) * *mdot;

		//time
		f[7] = r;
	}

}}}