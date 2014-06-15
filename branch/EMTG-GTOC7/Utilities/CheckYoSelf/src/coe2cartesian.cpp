#include "GTOC7_solution_check.h"


void coe2cartesian(Body * body)
{
	double mu_sun = 132712440018.0;
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