//function to compute the insertion burn

#include "EMTG_math.h"
#include <math.h>

using namespace EMTG;

namespace EMTG { namespace Astrodynamics {

double insertion_burn(const double *Vin, const double *Vplanet, const double mu, const double R_SOI, const double destination_a, const double destination_e, double* v_inf_inbound)
{
	//function to determine the magnitude of the insertion burn using the
	//patched conic method
	//insertion occurs at (t-epoch)

	double dV_heliocentric[3];
	double r_p, v_p_ellipse, v_p_hyperbola, dV;

	//get the inbound velocity in km/s
	for (int i=0;i<3;++i)
		dV_heliocentric[i] = Vin[i] - Vplanet[i];
	
	*v_inf_inbound = math::norm(dV_heliocentric,3) + 1.0e-10;

	//find the periapse distance
	r_p = destination_a * (1 - destination_e);

	//find the velocity of the spacecraft at periapse of the inbound hyperbola
	v_p_hyperbola = sqrt(*v_inf_inbound * *v_inf_inbound + 2*mu/r_p);

	//find the velocity of the spacecraft at periapse of the initial elliptical
	//orbit
	v_p_ellipse = sqrt(mu*(2/r_p - 1/destination_a));

	//calculate the final deltaV in km/s
	dV = fabs(v_p_hyperbola - v_p_ellipse);

	//add a penalty to ensure that the incoming orbit is indeed a hyperbola - otherwise the above model does not work
	//TODO insertion_burn needs to apply an explicit constraint rather than a penalty, since penalties are no longer used
	if (((*v_inf_inbound* *v_inf_inbound*0.9*0.9)/2 - mu/R_SOI) < 0) dV += 1/(*v_inf_inbound);

	return dV;
}

}}