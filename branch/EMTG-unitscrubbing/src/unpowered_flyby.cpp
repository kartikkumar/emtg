//code for unpowered flyby
//Algorithm by Dario Izzo 2006?
//appears in many places in published literature
//implementation by Jacob Englander

#include <cmath>
#include "Astrodynamics.h"
#include "EMTG_math.h"

using namespace EMTG;

namespace EMTG{ namespace Astrodynamics {

void unpowered_flyby (const double* v_initial, const double *planet_v, const double mu, //INPUT
	           const double rplanet, const double R_SOI, double rp, const double beta,//INPUT
	           double* v_final, double* outgoing_C3, double* flyby_delta, double *flyby_orbit_energy)//OUTPUT
{
	double Vinf[3];
	double Vinf_mag, e;
	
	//wrap the input variables into the format for this function
	for (int k=0;k<3;++k)
	{
		Vinf[k] = v_initial[k] - planet_v[k];
	}

	//find the magnitude of Vinf_in
	Vinf_mag = sqrt(Vinf[0]*Vinf[0] + Vinf[1]*Vinf[1] + Vinf[2]*Vinf[2]);

	//convert rp from planetary radii to km
	rp *= rplanet;

	//find the eccentricity of the flyby
	e = 1 + (rp*Vinf_mag*Vinf_mag/mu);

	//find the flyby turning angle
	*flyby_delta = 2 * asin(1/e);

	//find the local unit vectors for the flyby
	double ihat[3], jhat[3], khat[3], jnorm;
	//first find ihat
	for (int k=0;k<3;++k)
		ihat[k] = Vinf[k]/Vinf_mag;
	//then find jhat
	math::cross(ihat, planet_v, jhat);
	jnorm = sqrt(jhat[0]*jhat[0] + jhat[1]*jhat[1] + jhat[2]*jhat[2]);
	for (int k=0;k<3;++k)
		jhat[k] /= jnorm;
	//and finally khat
	math::cross(ihat, jhat, khat);

	//now compute the unit vector in the direction of V_inf_out
	double cdelta = cos(*flyby_delta);
	double cbeta = cos(beta);
	double sdelta = sin(*flyby_delta);
	double sbeta = sin(beta);
	for (int k=0;k<3;++k)
	{
		Vinf[k] = cdelta*ihat[k] + cbeta*sdelta*jhat[k] + sbeta*sdelta*khat[k];
	}

	//find the actual V_inf_out by multiplying the unit vector by the magnitude
	for (int k=0;k<3;++k)
		Vinf[k] *= Vinf_mag;

	//and finally find the heliocentric state_final
	for (int k=0;k<3;++k)
	{
		v_final[k] = Vinf[k] + planet_v[k];
	}

	//compute the outgoing C3
	*outgoing_C3 = Vinf_mag*Vinf_mag;

	//apply a penalty to low-velocity flybys
	//we need the following quantity, the energy of the spacecraft at the flyby with a 10% fudge factor margin, to be positive
	*flyby_orbit_energy = (Vinf_mag*Vinf_mag*0.9*0.9)/2 - mu/R_SOI;

	return;
}

}}