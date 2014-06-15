#include "GTOC7_solution_check.h"


void cartesian2coe(Body * body, double & mu_sun)
{
	
	double x = body->x;
	double y = body->y;
	double z = body->z;
	double vx = body->vx;
	double vy = body->vy;
	double vz = body->vz;
	double r = sqrt(x*x + y*y + z*z);
	double v = sqrt(vx*vx + vy*vy + vz*vz);
	double h;
	double n;

	std::vector <double> I = { 1.0, 0.0, 0.0 };
	std::vector <double> J = { 0.0, 1.0, 0.0 };
	std::vector <double> K = { 0.0, 0.0, 1.0 };
	std::vector <double> rvec = { x, y, z };
	std::vector <double> vvec = {vx, vy, vz};
	std::vector <double> evec (3, 0.0);
	std::vector <double> hvec (3, 0.0);
	std::vector <double> hvecHat (3, 0.0);
	std::vector <double> nvec (3, 0.0);
	

	for(size_t i = 0; i < 3; ++i)
		evec[i] = (1.0/mu_sun)*((v*v-mu_sun/r)*rvec[i]-dot(rvec,vvec)*vvec[i]);
	body->ecc =  sqrt(evec[0]*evec[0] + evec[1]*evec[1] + evec[2]*evec[2]);

	body->a = r/(2.0 - (r*v*v)/mu_sun);
	
	hvec = cross(rvec,vvec);
	h = sqrt(hvec[0]*hvec[0] + hvec[1]*hvec[1] + hvec[2]*hvec[2]);
	for(size_t i = 0; i < 3; ++i)
		hvecHat[i] = hvec[i]/h;


	nvec = cross(K,hvecHat);
	n = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);

	body->inc = acos(dot(hvec, K)/h);

	body->omega = acos(dot(nvec, evec)/(n*(body->ecc)));
	if(dot(evec, K) < 0)
		body->omega = 2.0 * PI - body->omega;

	body->LAN = acos(dot(nvec,I)/n);
	if(dot(nvec, J) < 0)
		body->LAN = 2.0*PI-body->LAN;

	body->tru = acos(dot(evec,rvec)/(body->ecc*body->r));
	if(dot(rvec,vvec) < 0)
		body->tru = 2.0*PI-body->tru;

	body->E = 2.0*atan(sqrt((1.0-body->ecc)/(1.0+body->ecc))*tan((body->tru)/2));
	body->M = (body->E)-(body->ecc)*sin(body->E);

	body->r = r;
	body->v = v;
	body->h = h;
}

void cartesian2coe_probe(Spacecraft * probe, const int & phase, const int & timestep, double & mu_sun)
{
	
	double x = probe->x[phase][timestep];
	double y = probe->y[phase][timestep];
	double z = probe->z[phase][timestep];
	double vx = probe->vx[phase][timestep];
	double vy = probe->vy[phase][timestep];
	double vz = probe->vz[phase][timestep];
	double r = sqrt(x*x + y*y + z*z);
	double v = sqrt(vx*vx + vy*vy + vz*vz);
	double h;
	double n;

	std::vector <double> I = { 1, 0, 0 };
	std::vector <double> J = { 0, 1, 0 };
	std::vector <double> K = { 0, 0, 1 };
	std::vector <double> rvec = { x, y, z };
	std::vector <double> vvec = { vx, vy, vz };
	std::vector <double> evec(3, 0.0);
	std::vector <double> hvec(3, 0.0);
	std::vector <double> hvecHat(3, 0.0);
	std::vector <double> nvec(3, 0.0);

	probe->r[phase].push_back(r);
	probe->v[phase].push_back(v);

	for (size_t i = 0; i < 3; ++i)
		evec[i] = (1.0 / mu_sun)*((v*v - mu_sun / r)*rvec[i] - dot(rvec, vvec)*vvec[i]);
	probe->ecc[phase].push_back( sqrt(evec[0] * evec[0] + evec[1] * evec[1] + evec[2] * evec[2]) );

	probe->a[phase].push_back( r / (2.0 - (r*v*v) / mu_sun) );

	hvec = cross(rvec, vvec);
	h = sqrt(hvec[0] * hvec[0] + hvec[1] * hvec[1] + hvec[2] * hvec[2]);
	for (size_t i = 0; i < 3; ++i)
		hvecHat[i] = hvec[i] / h;


	nvec = cross(K,hvecHat);
	n = sqrt(nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]);

	probe->inc[phase].push_back( acos(dot(hvec, K) / h) );

	probe->omega[phase].push_back ( (dot(nvec, evec) / (n*(probe->ecc[phase][timestep]))) );
	if (dot(evec, K) < 0)
		probe->omega[phase][timestep] = 2.0 * PI - probe->omega[phase][timestep];

	probe->LAN[phase].push_back( acos(dot(nvec, I) / n) );
	if (dot(nvec, J) < 0)
		probe->LAN[phase][timestep] = 2.0*PI - probe->LAN[phase][timestep];

	probe->tru[phase].push_back( (dot(evec, rvec) / (probe->ecc[phase][timestep] * probe->r[phase][timestep])) );
	if (dot(rvec, vvec) < 0)
		probe->tru[phase][timestep] = 2.0*PI - probe->tru[phase][timestep];

	probe->E[phase].push_back( 2.0*atan(sqrt((1.0 - probe->ecc[phase][timestep]) / (1.0 + probe->ecc[phase][timestep]))*tan((probe->tru[phase][timestep]) / 2)) );
	probe->M[phase].push_back( (probe->E[phase][timestep]) - (probe->ecc[phase][timestep])*sin(probe->E[phase][timestep]) );

	
	probe->h[phase].push_back( h );
}