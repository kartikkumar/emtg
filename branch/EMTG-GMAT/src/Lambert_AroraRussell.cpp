//Lambert solver by Arora and Russell
//A FAST AND ROBUST MULTIPLE REVOLUTION LAMBERT ALGORITHM USING A COSINE TRANSFORMATION
//AAS Hilton Head 2013
//implemented by Jacob Englander in C++ for EMTG, 3-5-2014

#include "Lambert_AroraRussell.h"
#include "EMTG_math.h"

#include <cmath>

namespace EMTG { namespace Astrodynamics {
	//assume input and output units are consistent (i.e. designed for km and s but should work for anything)
	void Lambert_AroraRussell(const double* R1,
							  const double* R2, 
							  const double& TOF,
							  const double& mu,
							  const int& Nrev,
							  const bool& LongWay,
							  const double& tolerance,
							  const int& max_iterations,
							  double* V1,
							  double* V2,
							  double& error,
							  int& iterations)
	{
		//Step 0: declare some variables
		double k; //iteration variable
		double deltak = 1.0e+10; //error in current solution for k
		double r1 = math::norm(R1, 3); //magnitude of initial position
		double r2 = math::norm(R2, 3); //magnitude of final position
		double ctheta = math::dot(R1, R2, 3) / (r1 * r2); //cosine of the transfer angle
		double S = sqrt((r1 + r2) * (r1 + r2) * (r1 + r2) / mu); //semi-perimeter
		double theta = acos(ctheta); //transfer angle
		iterations = 0; //iterations counter
		double sq2 = sqrt(2.0);
		double eps = 2.0e-2; //to prevent singularity at k = sqrt(2)
		double d = theta <= math::PI ? 1.0 : -1.0;
		double tau = d * sqrt(r1 * r2 * (1 + ctheta)) / (r1 + r2); //lambert geometry parameter

		//Step 1: generate appropriate initial guess
		//declare some variables that will be used in the initial guess

		//Step 1.1 compare the desired time of flight to the parabolic time of flight
		double T_parabolic = S * sqrt(1 - sq2*tau) * (1 + sq2*tau) / 3.0;
		
		if (TOF <= T_parabolic) //a hyperbolic trajectory is desired
		{
			double k_n, k_m, k_i, Z, alpha, F_0, F_1, F_i, F_star;
			double TOF20 = S * sqrt(1.0 - 20.0 * tau) * (tau + 0.04940968903 * (1.0 - 20.0 * tau));
			double TOF100 = S * sqrt(1.0 - 100.0 * tau) * (tau + 0.00999209404 * (1.0 - 100.0 * tau));
			if (d > 0)
			{
				k_n = sq2;
				k_m = 1.0/tau;
				k_i = (k_n + k_m) / 2.0;
				Z = 1.0 / sq2;
				alpha = 0.5;
				F_0 = T_parabolic;
				F_1 = 0.0;
				double m_i = 2 - k_i * k_i;
				double W = compute_W(k_i, m_i, Nrev);
				F_i = S * sqrt(1 - k_i*tau) * (tau + (1 - k_i*tau) * W);
				F_star = TOF;

				double x_star = pow( (Z * (F_0 - F_star)*(F_1 - F_i) ) / ( (F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star) ), 1.0/alpha);
				k = k_n + (k_m - k_n) * x_star;
			}
			else if (TOF < TOF20) //Arora's "H1" region
			{
				k_n = sq2;
				k_m = 20.0;
				k_i = (2*k_n + k_m) / 3.0;
				Z = 1.0 / 3.0;
				alpha = 1.0;
				F_0 = T_parabolic;
				F_1 = TOF20;
				double m_i = 2 - k_i * k_i;
				double W = compute_W(k_i, m_i, Nrev);
				F_i = S * sqrt(1 - k_i*tau) * (tau + (1 - k_i*tau) * W);
				F_star = TOF;

				double x_star = pow( (Z * (F_0 - F_star)*(F_1 - F_i) ) / ( (F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star) ), 1.0/alpha);
				k = k_n + (k_m - k_n) * x_star;
			}
			else //Arora's "H2" region
			{
			}
		}
		else if (Nrev >= 1.0) //a multi-revolution elliptical orbit is desired
		{
		}
		else if (Nrev < 1.0) // a single-revolution elliptical orbit is desired
		{
		}

		//Step 2: iterate to find k
		while (fabs(deltak) > tolerance && iterations < max_iterations)
		{
			//Step 2.1 increment the iterations counter
			++iterations;

			//Step 2.2 compute W, dW, ddW
			double m = 2 - k*k;
			double sgnk = k >= 0 ? 1.0 : -1.0;
			double W, dW, ddW;

			W = compute_W(k, m, Nrev);

			//dW and ddW
			if (k < sq2 - eps)
			{
				dW = (-2.0 + 3.0 * W * k) / m;
				ddW = (5.0 * dW * k + 3.0 * W) / m;
			}
			else if (k > sq2 + eps)
			{
				dW = (-2.0 + 3.0 * W * k) / -m;
				ddW = (5.0 * dW * k + 3.0 * W) / -m;
			}
			else
			{
				double v = k - sq2;
				double v2 = v*v;
				double v3 = v*v2;
				double v4 = v3*v;
				double v5 = v4*v;
				double v6 = v5*v;
				double v7 = v6*v;
				dW = - 1 / 5.0
					+ sq2 * v * (4.0/35.0)
					- v2 * (6.0 / 63.0)
					+ sq2 * v3 * (8.0 / 231.0)
					- v4 * (10.0 / 429.0)
					+ sq2 * v5 * (48.0 / 6435.0)
					- v6 * (56.0 / 12155.0)
					+ sq2 * v7 * (64.0 / 46189.0);
				ddW = sq2 * (4.0/35.0)
					- v * (12.0 / 63.0)
					+ sq2 * v2 * (24.0 / 231.0)
					- v3 * (40.0 / 429.0)
					+ sq2 * v4 * (240.0 / 6435.0)
					- v5 * (336.0 / 12155.0)
					+ sq2 * v6 * (448.0 / 46189.0);
			}
			

			//Step 2.3 compute TOFc, dTOFc, ddTOFc
			double TOFc = S * sqrt(1 - k*tau) * (tau + (1 - k*tau) * W);
			double c = (1 - k * tau) / tau;
			double sqrtctau = sqrt(1 - k*tau);
			double dTOFc = -TOFc / (2.0 * c) + S * tau * sqrtctau * (dW * c - W);
			double ddTOFc = - TOFc / (4.0 * c*c) + S * tau * sqrtctau * (W/c + c*ddW - 3.0*dW);

			//Step 2.4 compute deltak
			deltak = -(TOFc - TOF) / (dTOFc - (TOFc - TOF) * ddTOFc / (2.0 * dTOFc));

			//Step 2.5 update k from deltak
			k += deltak;
		}

		//Step 3: compute f, g, gdot (we don't need fdot)
		double f = 1 - (1 - k * tau) / r1;
		double g = tau * (r1 + r2) * sqrt(1 - k * tau);
		double gdot = 1 - (1 - k * tau) / r2;

		//Step 4: compute V1 and V2
		V1[0] = (R2[0] - f * R1[0]) / g;
		V1[1] = (R2[1] - f * R1[1]) / g;
		V1[2] = (R2[2] - f * R1[2]) / g;
		V2[0] = (gdot * R2[0] - R1[0]) / g;
		V2[1] = (gdot * R2[1] - R1[1]) / g;
		V2[2] = (gdot * R2[2] - R1[2]) / g;
	}

	//fast computation of acosh
	double acoshAR(const double& b)
	{
		return log(b + sqrt(b*b - 1.0));
	}

	//fast, rough computation of acos
	double acosAR(const double& x)
	{
		double coeff;
		double fx = fabs(x);
		double sgnx = x >= 0.0 ? 1.0 : -1.0;

		if ( fx <= 0.6 )
			coeff = (0.000014773722 + (1.1782782 - 0.52020038 * fx) * fx) / (1.1793469 + (-0.53277664 - 0.14454764 * fx) * fx);
		else if ( fx <= 0.97 )
			coeff = (0.011101554 + (8.9810074 + (-14.816468 + 5.9249913 * fx) * fx) * fx) / (9.2299851 + (-16.001036 + 6.8381053 * fx) * fx);
		else if ( fx <= 0.99 )
			coeff = (-35.750586 + (107.24325 - 70.780244 * fx) * fx) / (27.105764 - 26.638535 * fx);

		else
			coeff = asin(fx);
		
		return math::PIover2 - sgnx * coeff;
	}

	//function to compute the parameter W
	double compute_W(const double& k,
					 const double& m,
					 const int& Nrev)
	{
		const double sq2 = sqrt(2.0);
		const double eps = 2.0e-2;
		int sgnk = k < 0.0 ? -1 : 1;
		if (-sq2 <= k && k < sq2 - eps) //elliptical orbit case
				return ( (1 - sgnk) * math::PI + sgnk*acos(1-m) + 2*math::PI*Nrev ) / sqrt(m*m*m) - k/m;
			else if (k > sq2 + eps) //hyperbolic orbits
				return acoshAR(1 - m) / sqrt(-m*m*m) - k/m;
			else if (sq2 - eps <= k && k <= sq2 + eps) //Nrev = 0 case
			{
				double v = k - sq2;
				double v2 = v*v;
				double v3 = v*v2;
				double v4 = v3*v;
				double v5 = v4*v;
				double v6 = v5*v;
				double v7 = v6*v;
				double v8 = v7*v;
				return sq2 / 3.0
					- v / 5.0
					+ sq2 * v2 * (2.0/35.0)
					- v3 * (2.0 / 63.0)
					+ sq2 * v4 * (2.0 / 231.0)
					- v5 * (2.0 / 429.0)
					+ sq2 * v6 * (8.0 / 6435.0)
					- v7 * (8.0 / 12155.0)
					+ sq2 * v8 * (8.0 / 46189.0);
			}
			else
			{
				throw 200000;
				return math::LARGE;
			}
	}
}} //close namespace