//source file for orbit element conversions
#include <cmath>
#include <iostream>

#include "EMTG_math.h"




namespace EMTG{ namespace Astrodynamics {

//convert classical orbit elements to modified equinoctial elements
//see http://www.cdeagle.com/pdf/mee.pdf

void COE2MEE(const double* E_COE, double* E_MEE)
{
	double e = E_COE[1];
	double RAAN = E_COE[3];
	double w_plus_RAAN = E_COE[4] + RAAN;
	double ti = tan(E_COE[2]/2);

	//p
	E_MEE[0] = E_COE[0] * (1 - e*e);

	//f
	E_MEE[1] = e * cos(w_plus_RAAN);

	//g
	E_MEE[2] = e * sin(w_plus_RAAN);

	//h
	E_MEE[3] = ti * cos(RAAN);

	//k
	E_MEE[4] = ti * sin(RAAN);

	//L
	E_MEE[5] = w_plus_RAAN + E_COE[5];
}

//convert modified equinoctial elements to classical elements
void MEE2COE(const double* E_MEE, double* E_COE)
{
	double p = E_MEE[0];
	double f = E_MEE[1];
	double g = E_MEE[2];
	double h = E_MEE[3];
	double k = E_MEE[4];

	//a
	E_COE[0] = p / (1 - f*f - g*g);

	//e
	E_COE[1] = sqrt(f*f + g*g);

	//i
	E_COE[2] = atan2(2*sqrt(h*h + k*k), 1 - h*h - k*k);

	//RAAN
	E_COE[3] = atan2(k, h);

	//w
	E_COE[4] = atan2(g*h - f*k, f*h + g*k);

	//f
	E_COE[5] = E_MEE[5] - E_COE[3] - E_COE[4];
}

void inertial2COE(const double* state, const double mu, double* E_COE)
{
	double r, v, h, n, rdotv, s;
	double evec[3], nvec[3], hvec[3];
	double ihat[] = {1, 0, 0};
	double jhat[] = {0, 1, 0};
	double khat[] = {0, 0, 1};

	r = EMTG::math::norm(state, 3);
	v = EMTG::math::norm(state+3, 3);

	//a
	E_COE[0] = r / (2 - r*v*v/mu);

	//angular momentum vector and scalar
	EMTG::math::cross(state, state+3, hvec);
	h = EMTG::math::norm(hvec, 3);

	//eccentricity vector
	rdotv = EMTG::math::dot(state, state+3, 3);
	s = (v*v - mu/r);

	evec[0] = 1/mu * (s*state[0] - rdotv*state[3]);
	evec[1] = 1/mu * (s*state[1] - rdotv*state[4]);
	evec[2] = 1/mu * (s*state[2] - rdotv*state[5]);

	//eccentricity scalar
	E_COE[1] = EMTG::math::norm(evec, 3);

	//nodal vector
	EMTG::math::cross(khat, hvec, nvec);
	for (size_t k=0;k<3;++k)
		nvec[k] /= h;

	n = EMTG::math::norm(nvec, 3);

	//inclination
	E_COE[2] = acos(hvec[2] / h);

	//RAAN

	if (nvec[1] >= 0)
		E_COE[3] = acos(nvec[0] / n);
	else
		E_COE[3] = 2*EMTG::math::PI - acos(nvec[0] / n);

	if (n == 0)
		E_COE[3] = 0;

	//w
	if (evec[2] >= 0)
		E_COE[4] = acos(EMTG::math::dot(nvec, evec, 3) / (n * E_COE[1]));
	else
		E_COE[4] = 2*EMTG::math::PI - acos(EMTG::math::dot(nvec, evec, 3) / (n * E_COE[1]));

	if (n == 0) //if no inclination, then eccentricity vector points to the periapse
		E_COE[4] = atan2(evec[1], evec[0]);

	//f
	if (rdotv >= 0)
		E_COE[5] = acos(EMTG::math::dot(evec, state, 3) / (r * E_COE[1]));
	else
		E_COE[5] = 2*EMTG::math::PI - acos(EMTG::math::dot(evec, state, 3) / (r * E_COE[1]));
}

void COE2inertial(const double* E_COE, const double mu, double* state)
{
	double a = E_COE[0];
	double e = E_COE[1];
	double i = E_COE[2];
	double W = E_COE[3];
	double w = E_COE[4];
	double f = E_COE[5];
	double theta = w + f;

	double cosi = cos(i);
	double sini = sin(i);
	double cosW = cos(W);
	double sinW = sin(W);
	double cosw = cos(w);
	double sinw = sin(w);
	double ctheta = cos(theta);
	double stheta = sin(theta);

	double p = a * (1 - e*e);

	if (p < 1.0e-30)
	{
		std::cout << "Error converting parabolic orbit to Cartesian coordinates" << std::endl;
		throw 25;
	}

	double r = p / (1 + e * cos(f));

	double h = sqrt(mu * a * (1 - e*e));

	//position
	state[0] = r * (cosW*ctheta - sinW*stheta*cosi);
	state[1] = r * (sinW*ctheta + cosW*stheta*cosi);
	state[2] = r * stheta * sini;

	//velocity
	state[3] = -mu/h * (cosW * (stheta + e*sinw) + sinW * (ctheta + e*cosw) * cosi);
	state[4] = -mu/h * (sinW * (stheta + e*sinw) - cosW * (ctheta + e*cosw) * cosi);
	state[5] =  mu/h * (ctheta + e*cosw) * sini;

	
}

}} // namespace EMTG