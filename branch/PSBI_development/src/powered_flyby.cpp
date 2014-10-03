//function to compute the deltaV required for a powered flyby
//appears in many places in published literature
//for use with HOCP or as a stand-alone C++ file

#include <cmath>
#include "Astrodynamics.h"
#include "EMTG_math.h"

using namespace std;

namespace EMTG{ namespace Astrodynamics {

void powered_flyby (const double *v_initial, const double *v_final, const double *planet_v, const double mu, //INPUT
	           const double rplanet, const double R_SOI,//INPUT
	           double *dVmag, double *flyby_altitude, double* flyby_delta, double* outgoing_C3, double *flyby_orbit_energy)//OUTPUT
{
	double Vinf_in[3], Vinf_out[3];
	double Vinf_in_mag, Vinf_out_mag, rp_new, e_in, e_out, f, df, rp;
	int iter;
	
	//wrap the input variables into the format for this function
	for (int k = 0; k < 3; ++k)
	{
		Vinf_in[k] = v_initial[k] - planet_v[k];
		Vinf_out[k] = v_final[k] - planet_v[k];
	}

	double denom = sqrt(EMTG::math::dot(Vinf_in, Vinf_in, 3)*EMTG::math::dot(Vinf_out, Vinf_out, 3));
	*flyby_delta = acos(EMTG::math::dot(Vinf_in, Vinf_out, 3)/((denom!=0) ? denom : 1e-6)); //swing-by turning angle

	//convert Vinf_in and Vinf_out to km/s
	Vinf_in_mag = sqrt(EMTG::math::dot(Vinf_in, Vinf_in, 3));
	Vinf_out_mag = sqrt(EMTG::math::dot(Vinf_out, Vinf_out, 3));

	//find flyby_rp via Newton iteration
	rp_new = rplanet;//initial guess for the flyby periapse distance
	rp = 0;
	iter = 0;
	while (abs(rp - rp_new) > 1e-8 && iter <  30)
	{
		rp = rp_new;
		
		//compute the incoming and outgoing eccentricities
		e_in = 1 + Vinf_in_mag*Vinf_in_mag*rp/mu;
		e_in = (e_in!=0) ? e_in : 1e-4;
		e_out = 1 + Vinf_out_mag*Vinf_out_mag*rp/mu;
		e_out = (e_out!=0) ? e_out : 1e-4;
		
		//compute functions for Newton iteration
		f = asin(1/e_in) + asin(1/e_out) - *flyby_delta;
		df = -Vinf_in_mag*Vinf_in_mag / (mu * e_in*e_in * sqrt(1 - 1/(e_in*e_in))) - Vinf_out_mag*Vinf_out_mag / (mu * e_out*e_out * sqrt(1 - 1/(e_out*e_out)));
		
		//compute the new flyby_rp
		rp_new = rp - f/df;
		
		rp_new = rp_new < 0.0 ? 0.01 : rp_new;

		iter = iter + 1;
	}

	*flyby_altitude = rp - rplanet;
	
	//compute the dV in km/s
	*dVmag = abs(sqrt(Vinf_in_mag*Vinf_in_mag + 2*mu/rp) - sqrt(Vinf_out_mag*Vinf_out_mag + 2*mu/rp));
	
	//compute outgoing C3
	*outgoing_C3 = Vinf_out_mag * Vinf_out_mag;

	//apply a penalty to low-velocity flybys
	//we need the following quantity, the energy of the spacecraft at the flyby with a 10% fudge factor margin, to be positive
	*flyby_orbit_energy = (Vinf_in_mag*Vinf_in_mag*0.9*0.9)/2 - mu/R_SOI;

	return;
}

}}