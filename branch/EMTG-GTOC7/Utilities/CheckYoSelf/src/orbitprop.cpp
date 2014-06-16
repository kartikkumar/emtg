#include "GTOC7_solution_check.h"


void orbitprop(Body * body, const double & delta_t, const double & mu_sun)
{
	body->M = body->M + sqrt(mu_sun/((body->a)*(body->a)*(body->a)))*delta_t;

	//Correct the final value of M for over 2*pi rotations
	while(body->M>2*PI)
	{
		body->M = body->M-2*PI;
	}
	
	
	//Calculate new E and tru values as well
	body->E = laguerreConway(body->ecc, body->M);
	body->tru = 2.0*atan(sqrt((1.0+body->ecc)/(1.0-body->ecc))*tan(body->E/2.0));

	//Update Cartesian coordinates also
	coe2cartesian(body, mu_sun);
}

void orbitprop(Spacecraft * probe, const double & delta_t, const double & phase, const double & timestep, const double & mu_sun)
{
	probe->M[phase][timestep] = probe->M[phase][timestep] + sqrt(mu_sun / ((probe->a[phase][timestep])*(probe->a[phase][timestep])*(probe->a[phase][timestep])))*delta_t;

	//Correct the final value of M for over 2*pi rotations
	while (probe->M[phase][timestep] > 2.0 * PI)
	{
		probe->M[phase][timestep] = probe->M[phase][timestep] - 2.0 * PI;
	}


	//Calculate new E and tru values as well
	probe->E[phase][timestep] = laguerreConway(probe->ecc[phase][timestep], probe->M[phase][timestep]);
	probe->tru[phase][timestep] = 2.0*atan(sqrt((1.0 + probe->ecc[phase][timestep]) / (1.0 - probe->ecc[phase][timestep]))*tan(probe->E[phase][timestep] / 2.0));

	//Update Cartesian coordinates also
	coe2cartesian(probe, phase, timestep, mu_sun);
}