#include "GTOC7_solution_check.h"


void coe2cartesian(Body * body, const double & mu_sun)
{
	
	double a = body->a;
	double ecc = body->ecc;
	double inc = body->inc;
	double LAN = body->LAN;
	double omega = body->omega;
	double M = body->M;
	double E = body->E;
	double tru = body->tru;
	
	

	body->r = a*(1.0-ecc*ecc)/(1.0+ecc*cos(tru));
	body->v = sqrt(2.0*mu_sun/(body->r)-mu_sun/a);
	body->h = sqrt(mu_sun*a*(1.0-ecc*ecc));

	body->x = (body->r)*(cos(tru+omega)*cos(LAN)-sin(tru+omega)*cos(inc)*sin(LAN));
	body->y = (body->r)*(cos(tru+omega)*sin(LAN)+sin(tru+omega)*cos(inc)*cos(LAN));
    body->z = (body->r)*(sin(tru+omega)*sin(inc));

    body->vx = -(mu_sun/(body->h))*(cos(LAN)*(sin(tru+omega)+ecc*sin(omega))+sin(LAN)*(cos(tru+omega)+ecc*cos(omega))*cos(inc));
    body->vy = -(mu_sun/(body->h))*(sin(LAN)*(sin(tru+omega)+ecc*sin(omega))-cos(LAN)*(cos(tru+omega)+ecc*cos(omega))*cos(inc));
    body->vz =  (mu_sun/(body->h))*(cos(tru+omega)+ecc*cos(omega))*sin(inc);
}

void coe2cartesian(Spacecraft * probe, const int & phase, const int & timestep, const double & mu_sun)
{
	
	double a = probe->a[phase][timestep];
	double ecc = probe->ecc[phase][timestep];
	double inc = probe->inc[phase][timestep];
	double LAN = probe->LAN[phase][timestep];
	double omega = probe->omega[phase][timestep];
	double M = probe->M[phase][timestep];
	double E = probe->E[phase][timestep];
	double tru = probe->tru[phase][timestep];



	probe->r[phase][timestep] = a*(1.0 - ecc*ecc) / (1.0 + ecc*cos(tru));
	probe->v[phase][timestep] = sqrt(2.0*mu_sun / (probe->r[phase][timestep]) - mu_sun / a);
	probe->h[phase][timestep] = sqrt(mu_sun*a*(1.0 - ecc*ecc));

	probe->x[phase][timestep] = (probe->r[phase][timestep])*(cos(tru + omega)*cos(LAN) - sin(tru + omega)*cos(inc)*sin(LAN));
	probe->y[phase][timestep] = (probe->r[phase][timestep])*(cos(tru + omega)*sin(LAN) + sin(tru + omega)*cos(inc)*cos(LAN));
	probe->z[phase][timestep] = (probe->r[phase][timestep])*(sin(tru + omega)*sin(inc));

	probe->vx[phase][timestep] = -(mu_sun / (probe->h[phase][timestep]))*(cos(LAN)*(sin(tru + omega) + ecc*sin(omega)) + sin(LAN)*(cos(tru + omega) + ecc*cos(omega))*cos(inc));
	probe->vy[phase][timestep] = -(mu_sun / (probe->h[phase][timestep]))*(sin(LAN)*(sin(tru + omega) + ecc*sin(omega)) - cos(LAN)*(cos(tru + omega) + ecc*cos(omega))*cos(inc));
	probe->vz[phase][timestep] = (mu_sun / (probe->h[phase][timestep]))*(cos(tru + omega) + ecc*cos(omega))*sin(inc);
}