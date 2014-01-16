//EMTG universal Kepler propagator
//Jacob Englander 3/31/2013
//adaptation of Shepperd's method, originally written in MATLAB by David Eagle

#include "UniversalKeplerPropagator.h"
#include "EMTG_math.h"

namespace EMTG { namespace Astrodynamics {

//constructor
UniversalKeplerPropagator::UniversalKeplerPropagator()
{
	rv0.resize(3, 1);
	vv0.resize(3, 1);
	rvf.resize(3, 1);
	vvf.resize(3, 1);

	hv.resize(3, 1);

	M.resize(3, 3);
	STM.resize(6, 6);
}

//destructor
UniversalKeplerPropagator::~UniversalKeplerPropagator(){}


//methods

int UniversalKeplerPropagator::propagate_and_compute_STM(const double* state0, double* statef, const double& deltat_seconds, const double& mu_seconds, const double& LU, const double& TU, bool compute_STM)
{
	// two body state transition matrix and propagation
	//
	// Shepperd's method
	//
	// input
	//
	//  tau = propagation time interval (seconds)
	//  ri  = initial eci position vector (kilometers)
	//  vi  = initial eci velocity vector (km/sec)
	//
	// output
	//
	//  rf  = final eci position vector (kilometers)
	//  vf  = final eci velocity vector (km/sec)
	//  stm = state transition matrix
	//
	// Adapted from David Eagle code
	// Orbital Mechanics with Matlab
	//
	//***************************************************

	// convergence criterion

	double tol = 1.0e-12;

	//decompose the initial and final state vectors
	rv0.assign_all(state0);
	vv0.assign_all(state0 + 3);

	rv0 /= LU;
	vv0 *= TU/LU;
	double mu = 1.0;
	double deltat = deltat_seconds / TU;

	double r0 = rv0.norm();
	double v0 = vv0.norm();
	double rf;

	//compute the initial angular momentum vector
	rv0.cross_in_place(vv0, hv);

	//compute the momentum
	double n0 = rv0.dot(vv0);
   
   
	double beta = (2 * mu / r0) - v0*v0;

	double u = 0;

	double umax = beta == 0 ? 1.0e+24 : 1 / sqrt(fabs(beta));
    
	double umin = -umax;

	double delu;

	if (beta <= 0)
	   delu = 0;
	else
	{
		double p = math::TwoPI * mu * sqrt(beta*beta*beta);
		double n = floor((1.0 / p) * (deltat + 0.5 * p - 2 * n0 / beta));
		delu = math::TwoPI * n * sqrt(beta*beta*beta*beta*beta);
	}

	double tsav = 1.0e99;
	int niter = 0;

	//declare the Universal functions
	double U0, U1, U2, U3;

	//declare other variables that will be used throughout
	double uu;
	int n, l, d, k;
	double a, b, g;

	// Kepler iteration loop

	while (1)
	{
		++niter;
      
		double q = beta * u * u;
		q /= (1 + q);
		U0 = 1 - 2 * q;
		U1 = 2 * u * (1 - q);
      
		// continued fraction iteration

		n = 0;
		l = 3;
		d = 15;
		k = -9;
		a = 1.0;
		b = 1.0;
		g = 1.0;
		
		int liter = 0;
		while (1)
		{
			++liter;
			double gsav = g;
			k *= -1;
			l += 2;
			d += 4 * l;
			n += (1 + k) * l;
			a = d / (d - n * a * q);
			b *= (a - 1);
			g += b;
    
			if (fabs(g - gsav) < tol)
				break; 
		}
  
		uu = (16.0 / 15.0) * U1 * U1 * U1 * U1 * U1 * g + delu;
       
		U2 = 2.0 * U1 * U1;
		U1 = 2.0 * U0 * U1;
		U0 = 2.0 * U0 * U0 - 1;
		U3 = beta * uu + U1 * U2 / 3.0;
      
		rf = r0 * U0 + n0 * U1 + mu * U2;
		double t = r0 * U1 + n0 * U2 + mu * U3;
       
		double  dtdu = 4 * rf * (1 - q);
     
		//check for time convergence

		if (fabs(t - tsav) < tol)
			break;

		double usav = u;
		tsav = t;
		double terr = deltat - t;

		if (fabs(terr) < fabs(deltat) * tol)
			break;  

		double du = terr / dtdu;
      
		if (du < 0)
		{
			umax = u;
			u = u + du;
			if (u < umin)
			u = 0.5 * (umin + umax);
		} 
		else
		{
			umin = u;
			u = u + du;
			if (u > umax)
			u = 0.5 * (umin + umax);
		}  

		// check for independent variable convergence

		if (fabs(u - usav) < tol)
			break;  
  
		// check for more than 50 iterations
  
		if (niter > 50)
		{
			std::cout << "Universal Kepler solver exceeded 50 iterations" << endl;
			throw 1000;   
		}
	}

	double uc = u;
	double fm = -mu * U2 / r0;
	double ggm = -mu * U2 / rf;
	double F = 1 + fm;
	double G = r0 * U1 + n0 * U2;
	double FF = -mu * U1 / (r0 * rf);
	double GG = 1 + ggm;

	// compute final state vector
	for (int k = 0; k < 3; ++k)
	{
		rvf(k) = rv0(k) * F + vv0(k) * G;
		vvf(k) = rv0(k) * FF + vv0(k) * GG;
	}

	//compute the state transition matrix using Shepperd's method

	if (compute_STM)
	{
		uu = G * U2 + 3 * mu * uu;

		double a0 = mu / (r0*r0*r0);
		double a1 = mu / (rf*rf*rf);

		//create the M matrix
		M.assign_zeros();

		M(0,0) = FF * (U0 / (r0 * rf) + 1 / (r0 * r0) + 1 / (rf * rf)) - a0 * a1 * uu;
		M(0,1) = (FF * U1 + (ggm / rf)) / rf;
		M(0,2) = ggm * U1 / rf - a1 * uu;
		M(1,0) = -(FF * U1 + (fm / r0)) / r0;
		M(1,1) = -FF * U2;
		M(1,2) = -ggm * U2;
		M(2,0) = fm * U1 / r0 - a0 * uu;
		M(2,1) = fm * U2;
		M(2,2) = G * U2 - uu;

		//create the state transition matrix
		for (int i = 0; i < 2; ++i)
		{
			int x = 2 * (i+1) - 3;
			int ib = 3 * (2 - (i+1));
			for (int j = 0; j < 2; ++j)
			{
				int jb = 3 * (j);
				for (int ii = 0; ii < 3; ++ii)
				{
					double t1 = rvf(ii) * M(i,j) + vvf(ii) * M(i+1,j);
					double t2 = rvf(ii) * M(i,j+1) + vvf(ii) * M(i+1,j+1);
					for (int jj = 0; jj < 3; ++jj)
					{
						STM(ii+ib, jj+jb) = x * (t1 * rv0(jj) + t2 * vv0(jj));
					}
				}
			}
		}

		for (int i = 0; i < 3; ++i)
		{
			int j = i + 3;
			STM(i,i) = STM(i,i) + F;
			STM(i,j) = STM(i,j) + G;
			STM(j,i) = STM(j,i) + FF;
			STM(j,j) = STM(j,j) + GG;
		}

		/*  For now, we will leave the STM scaled in the interest of numerical stability
		//unscale
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 3; j < 6; ++j)
			{
				STM(i,j) *= TU;
			}
		}

		for (int i = 3; i < 6; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				STM(i,j) /= TU;
			}
		}
		*/
	}

	//unscale the state vectors
	rvf *= LU;
	vvf *= LU/TU;

	for (int k = 0; k < 3; ++k)
	{
		statef[k] = rvf(k);
		statef[k+3] = vvf(k);
	}

	return 0;
}


}} //close namespace
