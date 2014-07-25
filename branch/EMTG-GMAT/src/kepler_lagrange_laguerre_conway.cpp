//Lagrange method Kepler solver

#include <math.h>
#include <iostream>

#include "kepler_lagrange_laguerre_conway.h"
#include "STM.h"
#include "EMTG_math.h"


namespace Kepler
{
	//Lagrange-Laguerre-Conway Kepler solver
	void KeplerLagrangeLaguerreConway(const double* state0, double* state, const double& mu, const double& propTime, double& F, double& G, double& Ft, double& Gt, double& Ftt, double& Gtt, STM& stm, const bool& compute_STM_flag)
	{
		double r;
		double deltaE = 0.0, deltaH = 0.0, cdeltaE = 0.0, sdeltaE = 0.0, chdeltaH = 0.0, shdeltaH = 0.0;
		double U0, U1, U2, U3, X;
		int iteration_count = 0;

		//scaling code
		/*
		double LU = 149597870.691;
		double TU = sqrt(LU*LU*LU/mu_in);
		double mu = 1.0;
		double state0[6];
		for (int k = 0; k < 3; ++k)
			state0[k] = state0_in[k] / LU;
		for (int k = 3; k < 6; ++k)
			state0[k] = state0_in[k] * TU/LU;
		double propTime = propTime_in / TU;
		*/
		//set Laguerre-Conway n
		const static int n = 5;

		//step 1: determine if we are an ellipse or a hyperbola/parabola
		double sqmu = sqrt(mu);
		double r0 = sqrt(state0[0] * state0[0] + state0[1] * state0[1] + state0[2] * state0[2]);
		double v0 = sqrt(state0[3] * state0[3] + state0[4] * state0[4] + state0[5] * state0[5]);

		//Battin equation 4.23
		double sigma0 = (state0[0] * state0[3] + state0[1] * state0[4] + state0[2] * state0[5]) / sqmu;
		double sqrtp = sqrt(2 * r0 * sigma0 * sigma0);
		double sigma;
		double a = 1.0 / (2.0 / r0 - v0*v0 / mu);
		double coeff = (1 - r0 / a);
		double Energy = v0*v0 / 2.0 - mu / r0;
		const static double alpha_tolerance = 1.0e-8;

		//alpha for the Universal Kepler solver - this takes different values for each type of conic section
		double alpha = fabs(a) > alpha_tolerance ? 1.0 / a : 0.0;

		if (alpha > alpha_tolerance) //ellipse
		{
			//step 2: solve Kepler's equation

			//Step 2.1: initialize ellipse-specific quantities
			double sqrta = sqrt(a);

			//Battin equation 4.35
			double deltaM = sqmu / (sqrta*sqrta*sqrta) * propTime;

			//Step 2.2: solve for deltaE using the Laguerre-Conway method
			//initial guess is deltaE = deltaM
			deltaE = 1.0e+100;
			double deltaE_new = deltaM;

			while (fabs(deltaE - deltaE_new) > 1.0e-12 && iteration_count < 100)
			{
				//iteration count update
				++iteration_count;

				//deltaE update
				deltaE = deltaE_new;

				//trig evaluations that only need to be done once
				cdeltaE = cos(deltaE);
				sdeltaE = sin(deltaE);

				//Kepler's equation for an elliptical orbit
				double F = deltaE - deltaM + sigma0 / sqrta * (1.0 - cdeltaE) - coeff * sdeltaE;

				//derivative with respect to deltaE
				double dF = 1 + sigma0 / sqrta * sdeltaE - coeff * cdeltaE;

				//second derivative with respect to deltaE
				double ddF = sigma0 / sqrta * cdeltaE + coeff * sdeltaE;

				//Laguerre-Conway update
				int sgn = dF >= 0 ? 1 : -1;
				deltaE_new = deltaE - n*F / (dF + sgn * sqrt(fabs( (n-1)*(n-1)*dF*dF - n * (n - 1) * F * ddF )));

				if ( !(deltaE_new == deltaE_new) )
				{
					std::cout << "deltaE_new is a NaN!" << std::endl;
					throw 1000000;
				}
			}

			//Step 3: find X, U0, U1, U2, U3
			double sqalpha = sqrt(alpha);
			X = sqrta * deltaE;
			double sqalphaX = sqalpha * X;
			U0 = cos(sqalphaX);
			U1 = sin(sqalphaX) / sqalpha;
			U2 = (1 - U0) / alpha;
			U3 = (sqalphaX / sqalpha - U1) / alpha;
		}
		else if (alpha < -alpha_tolerance) //hyperbola
		{
			//Step 2.1: initialize hyperbola-specific quantities
			double sqrtma = sqrt(-a);

			//Battin equation 4.51
			double deltaN = sqmu / sqrt(-a*a*a) * propTime;

			//Step 2.2: solve for deltaH using the Laguerre-Conway method
			//initial guess is deltaH = deltaN
			deltaH = 0;
			double deltaH_new = propTime < 0.0 ? -1 : 1;

			while (fabs(deltaH - deltaH_new) > 1.0e-12 && iteration_count < 100)
			{
				//iteration count update
				++iteration_count;

				//deltaE update
				deltaH = deltaH_new;

				//trig evaluations that only need to be done once
				chdeltaH = cosh(deltaH);
				shdeltaH = sinh(deltaH);

				//Kepler's equation for a hyperbolic orbit
				double F = -deltaN - deltaH + sigma0/sqrtma * (chdeltaH - 1.0) + coeff * shdeltaH;

				//derivative with respect to deltaH
				double dF = -1 + sigma0/sqrtma * shdeltaH + coeff * chdeltaH;

				//second derivative with respect to deltaH
				double ddF = sigma0/sqrtma * chdeltaH + coeff * shdeltaH;

				//Laguerre-Conway update
				int sgn = dF >= 0 ? 1 : -1;
				deltaH_new = deltaH - n*F / (dF + sgn * sqrt(fabs( (n-1)*(n-1)*dF*dF - n * (n - 1) * F * ddF )));

				if ( !(deltaH_new == deltaH_new) )
				{
					std::cout << "deltaH_new is a NaN!" << std::endl;
					throw 1000000;
				}
			}

			//Step 3: find X, U0, U1, U2
			double sqmalpha = sqrt(-alpha);
			X = sqrtma * deltaH;
			double sqmalphaX = sqmalpha * X;
			U0 = cosh(sqmalphaX);
			U1 = sinh(sqmalphaX) / sqmalpha;
			U2 = (U0 - 1) / -alpha;
			U3 = (U1 - sqmalphaX / sqmalpha) / -alpha;
		}
		else //parabola - solve via Universal method
		{
			//std::cout << "Oy! A parabola!" << std::endl;

			//initial guess for the universal variable X, from Prussing and Conway p39
			double rp;
			//angular momentum vector and scalar
			double hvec[3];
			EMTG::math::cross(state0, state0+3, hvec);
			double h = EMTG::math::norm(hvec, 3);

			if (fabs(alpha) < alpha_tolerance) //parabola
			{
				rp = h*h / mu;
			}
			else
			{
				//eccentricity vector
				double evec[3];
				double rdotv = EMTG::math::dot(state0, state0+3, 3);
				double s = (v0*v0 - mu/r0);

				evec[0] = 1.0/mu * (s*state0[0] - rdotv*state0[3]);
				evec[1] = 1.0/mu * (s*state0[1] - rdotv*state0[4]);
				evec[2] = 1.0/mu * (s*state0[2] - rdotv*state0[5]);

				//eccentricity scalar
				double e = EMTG::math::norm(evec, 3);

				//periapse distance
				rp = a * (1.0 - e);
			}
			double Xplus = sqmu * propTime / rp;
			//compute U0, U1, U2, U3 as per Battin Problem 4-21, p180
			U0 = 1.0;
			U1 = Xplus;
			U2 = Xplus*Xplus / 2.0;
			U3 = U2 * Xplus / 3.0;
		
			double FXplus = r0 * U1 + sigma0 * U2 + U3 - sqmu * propTime;
			double X_new = (mu * propTime * propTime) / (rp * FXplus + sqmu * propTime);

			//Step 2: solve Kepler's equation via a Universal Laguerre-Conway method
			X = 1.0e+100;
			while (fabs(X - X_new) > 1.0e-12 && iteration_count < 100)
			{
				++iteration_count;

				//Step 2.1
				//X update
				X = X_new;

				//Step 2.2
				//compute U0, U1, U2, U3 as per Battin Problem 4-21, p180
				U0 = 1.0;
				U1 = X;
				U2 = X*X / 2.0;
				U3 = U2 * X / 3.0;

				//Step 2.3 compute r and sigma via equations 4.82 and 4.83, p178
				r = r0 * U0 + sigma0 * U1 + U2;
				sigma = sigma0 * U0 + (1 - alpha * r0) * U1;

				//Step 2.4 Universal form of Kepler's equation and its derivatives
				double FX = r0 * U1 + sigma0 * U2 + U3 - sqmu * propTime;
				double dFX = r;
				double ddFX = sigma;

				//Laguerre-Conway update
				int sgn = dFX >= 0 ? 1 : -1;
				X_new = X - n*FX / (dFX + sgn * sqrt(fabs( (n-1)*(n-1)*dFX*dFX - n * (n - 1) * F * ddFX )));
			}

			if ( !(X == X) )
			{
				std::cout << "X is a NaN!" << std::endl;
				throw 1000000;
			}
		}

		//Step 3: find F, G, Ft, Gt
		r = r0 * U0 + sigma0 * U1 + U2;
		sigma = sigma0 * U0 + (1 - alpha * r0) * U1;
		F = 1.0 - U2 / r0;
		G = (r0 * U1 + sigma0 * U2) / sqmu;
		Ft = -sqmu / (r0 * r) * U1;
		Gt = 1.0 - U2 / r;

		//Step 4: compute the final state as functions of F, G, Ft, Gt
		//Battin equation 3.33
		state[0] = (F*state0[0] + G*state0[3]);// * LU;
		state[1] = (F*state0[1] + G*state0[4]);// * LU;
		state[2] = (F*state0[2] + G*state0[5]);// * LU;
		state[3] = (Ft*state0[0] + Gt*state0[3]);// * LU/TU;
		state[4] = (Ft*state0[1] + Gt*state0[4]);// * LU/TU;
		state[5] = (Ft*state0[2] + Gt*state0[5]);// * LU/TU;

		//Step 5: Compute the state transition matrix if requested
		if (compute_STM_flag)
		{
			//STM calculation code
			//Step 5.1 compute the universal functions Ui and their derivatives
			//from the recursion relation 4.76
			double U4 = a * (X*X / 2.0 - U2);
			double U5 = a * (X*X*X / 6.0 - U3);

			//derivatives
			double dXdt = sqmu / r;
			double U0dot, U1dot, U2dot;

			if (alpha > alpha_tolerance) //ellipse
			{
				double sqalpha = sqrt(alpha);
				double sqalphaX = sqalpha * X;
				U0dot =  -sqalpha * sin(sqalphaX) * dXdt;
				U1dot = U0 * dXdt;
				U2dot = U1 * dXdt;

				if (!(U0dot == U0dot && U1dot == U1dot && U2dot == U2dot))
				{
					std::cout << "Failure to calculate Udots" << std::endl;
				}
			}
			else if (alpha < -alpha_tolerance) //hyperbola
			{
				double sqmalpha = sqrt(-alpha);
				double sqmalphaX = sqmalpha * X;
				U0dot = sqmalpha * sinh(sqmalphaX);
				U1dot = U0 * dXdt;
				U2dot = U1 * dXdt;

				if (!(U0dot == U0dot && U1dot == U1dot && U2dot == U2dot))
				{
					std::cout << "Failure to calculate Udots" << std::endl;
				}
			}
			else //parabola
			{
				U0dot = 0;
				U1dot = dXdt;
				U2dot = X * dXdt;

				if (!(U0dot == U0dot && U1dot == U1dot && U2dot == U2dot))
				{
					std::cout << "Failure to calculate Udots" << std::endl;
				}
			}

			//Step 5.2 compute C, which along with F, G, Ft, and Gt can be used to compute the STM
			//Battin equation 9.74
			double C = 1.0/sqmu * (3*U5 - X*U4) - propTime * U2;

			//Step 5.3 compute R, V, R~ and V~, Battin equations 9.84 - 9.87
			
			//R~
			stm(0,0) = r/mu*(state[3] - state0[3])*(state[3] - state0[3]) + (r0*(1.0 - F)*(state[0]*state0[0]) + C*(state0[0]*state[3]))/(r0*r0*r0) + F;
			stm(0,1) = r/mu*(state[3] - state0[3])*(state[4] - state0[4]) + (r0*(1.0 - F)*(state[0]*state0[1]) + C*(state0[1]*state[3]))/(r0*r0*r0);
			stm(0,2) = r/mu*(state[3] - state0[3])*(state[5] - state0[5]) + (r0*(1.0 - F)*(state[0]*state0[2]) + C*(state0[2]*state[3]))/(r0*r0*r0);
			stm(1,0) = r/mu*(state[4] - state0[4])*(state[3] - state0[3]) + (r0*(1.0 - F)*(state[1]*state0[0]) + C*(state0[0]*state[4]))/(r0*r0*r0);
			stm(1,1) = r/mu*(state[4] - state0[4])*(state[4] - state0[4]) + (r0*(1.0 - F)*(state[1]*state0[1]) + C*(state0[1]*state[4]))/(r0*r0*r0) + F;
			stm(1,2) = r/mu*(state[4] - state0[4])*(state[5] - state0[5]) + (r0*(1.0 - F)*(state[1]*state0[2]) + C*(state0[2]*state[4]))/(r0*r0*r0);
			stm(2,0) = r/mu*(state[5] - state0[5])*(state[3] - state0[3]) + (r0*(1.0 - F)*(state[2]*state0[0]) + C*(state0[0]*state[5]))/(r0*r0*r0);
			stm(2,1) = r/mu*(state[5] - state0[5])*(state[4] - state0[4]) + (r0*(1.0 - F)*(state[2]*state0[1]) + C*(state0[1]*state[5]))/(r0*r0*r0);
			stm(2,2) = r/mu*(state[5] - state0[5])*(state[5] - state0[5]) + (r0*(1.0 - F)*(state[2]*state0[2]) + C*(state0[2]*state[5]))/(r0*r0*r0) + F;

			//R
			stm(0,3) = r0/mu*(1.0-F)*(state0[3]*(state[0] - state0[0]) - state0[0]*(state[3] - state0[3])) + C/mu*(state[3]*state0[3]) + G;
			stm(0,4) = r0/mu*(1.0-F)*(state0[4]*(state[0] - state0[0]) - state0[1]*(state[3] - state0[3])) + C/mu*(state[3]*state0[4]);
			stm(0,5) = r0/mu*(1.0-F)*(state0[5]*(state[0] - state0[0]) - state0[2]*(state[3] - state0[3])) + C/mu*(state[3]*state0[5]);
			stm(1,3) = r0/mu*(1.0-F)*(state0[3]*(state[1] - state0[1]) - state0[0]*(state[4] - state0[4])) + C/mu*(state[4]*state0[3]);
			stm(1,4) = r0/mu*(1.0-F)*(state0[4]*(state[1] - state0[1]) - state0[1]*(state[4] - state0[4])) + C/mu*(state[4]*state0[4]) + G;
			stm(1,5) = r0/mu*(1.0-F)*(state0[5]*(state[1] - state0[1]) - state0[2]*(state[4] - state0[4])) + C/mu*(state[4]*state0[5]);
			stm(2,3) = r0/mu*(1.0-F)*(state0[3]*(state[2] - state0[2]) - state0[0]*(state[5] - state0[5])) + C/mu*(state[5]*state0[3]);
			stm(2,4) = r0/mu*(1.0-F)*(state0[4]*(state[2] - state0[2]) - state0[1]*(state[5] - state0[5])) + C/mu*(state[5]*state0[4]);
			stm(2,5) = r0/mu*(1.0-F)*(state0[5]*(state[2] - state0[2]) - state0[2]*(state[5] - state0[5])) + C/mu*(state[5]*state0[5]) + G;

			//V~
			stm(3,0) = -(state0[0]*(state[3] - state0[3]))/(r0*r0) - (state[0]*(state[3] - state0[3]))/(r*r) + Ft*(1.0 - (state[0]*state[0])/(r*r) + ((state[1]*(state[0]*state[4] - state[3]*state[1]) + state[2]*(state[0]*state[5] - state[3]*state[2]))*(state[3] - state0[3]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[0]*state0[0]);
			stm(3,1) = -(state0[1]*(state[3] - state0[3]))/(r0*r0) - (state[0]*(state[4] - state0[4]))/(r*r) + Ft*(	   - (state[0]*state[1])/(r*r) + ((state[1]*(state[0]*state[4] - state[3]*state[1]) + state[2]*(state[0]*state[5] - state[3]*state[2]))*(state[4] - state0[4]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[0]*state0[1]);
			stm(3,2) = -(state0[2]*(state[3] - state0[3]))/(r0*r0) - (state[0]*(state[5] - state0[5]))/(r*r) + Ft*(    - (state[0]*state[2])/(r*r) + ((state[1]*(state[0]*state[4] - state[3]*state[1]) + state[2]*(state[0]*state[5] - state[3]*state[2]))*(state[5] - state0[5]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[0]*state0[2]);
			stm(4,0) = -(state0[0]*(state[4] - state0[4]))/(r0*r0) - (state[1]*(state[3] - state0[3]))/(r*r) + Ft*(    - (state[1]*state[0])/(r*r) + (-(state[0]*(state[0]*state[4] - state[3]*state[1]) + state[2]*(state[1]*state[5] - state[4]*state[2]))*(state[3] - state0[3]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[1]*state0[0]);
			stm(4,1) = -(state0[1]*(state[4] - state0[4]))/(r0*r0) - (state[1]*(state[4] - state0[4]))/(r*r) + Ft*(1.0 - (state[1]*state[1])/(r*r) + (-(state[0]*(state[0]*state[4] - state[3]*state[1]) + state[2]*(state[1]*state[5] - state[4]*state[2]))*(state[4] - state0[4]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[1]*state0[1]);
			stm(4,2) = -(state0[2]*(state[4] - state0[4]))/(r0*r0) - (state[1]*(state[5] - state0[5]))/(r*r) + Ft*(    - (state[1]*state[2])/(r*r) + (-(state[0]*(state[0]*state[4] - state[3]*state[1]) + state[2]*(state[1]*state[5] - state[4]*state[2]))*(state[5] - state0[5]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[1]*state0[2]);
			stm(5,0) = -(state0[0]*(state[5] - state0[5]))/(r0*r0) - (state[2]*(state[3] - state0[3]))/(r*r) + Ft*(    - (state[2]*state[0])/(r*r) + (-(state[0]*(state[0]*state[5] - state[3]*state[2]) - state[1]*(state[1]*state[5] - state[4]*state[2]))*(state[3] - state0[3]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[2]*state0[0]);
			stm(5,1) = -(state0[1]*(state[5] - state0[5]))/(r0*r0) - (state[2]*(state[4] - state0[4]))/(r*r) + Ft*(    - (state[2]*state[1])/(r*r) + (-(state[0]*(state[0]*state[5] - state[3]*state[2]) - state[1]*(state[1]*state[5] - state[4]*state[2]))*(state[4] - state0[4]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[2]*state0[1]);
			stm(5,2) = -(state0[2]*(state[5] - state0[5]))/(r0*r0) - (state[2]*(state[5] - state0[5]))/(r*r) + Ft*(1.0 - (state[2]*state[2])/(r*r) + (-(state[0]*(state[0]*state[5] - state[3]*state[2]) - state[1]*(state[1]*state[5] - state[4]*state[2]))*(state[5] - state0[5]))/(mu*r)) - mu*C/(r*r*r*r0*r0*r0)*(state[2]*state0[2]);

			//V
			stm(3,3) = r0/mu*(state[3] - state0[3])*(state[3] - state0[3]) + (r0*(1.0 - F)*(state[0]*state0[0]) - C*(state[0]*state0[3]))/(r*r*r) + Gt;
			stm(3,4) = r0/mu*(state[3] - state0[3])*(state[4] - state0[4]) + (r0*(1.0 - F)*(state[0]*state0[1]) - C*(state[0]*state0[4]))/(r*r*r);
			stm(3,5) = r0/mu*(state[3] - state0[3])*(state[5] - state0[5]) + (r0*(1.0 - F)*(state[0]*state0[2]) - C*(state[0]*state0[5]))/(r*r*r);
			stm(4,3) = r0/mu*(state[4] - state0[4])*(state[3] - state0[3]) + (r0*(1.0 - F)*(state[1]*state0[0]) - C*(state[1]*state0[3]))/(r*r*r);
			stm(4,4) = r0/mu*(state[4] - state0[4])*(state[4] - state0[4]) + (r0*(1.0 - F)*(state[1]*state0[1]) - C*(state[1]*state0[4]))/(r*r*r) + Gt;
			stm(4,5) = r0/mu*(state[4] - state0[4])*(state[5] - state0[5]) + (r0*(1.0 - F)*(state[1]*state0[2]) - C*(state[1]*state0[5]))/(r*r*r);
			stm(5,3) = r0/mu*(state[5] - state0[5])*(state[3] - state0[3]) + (r0*(1.0 - F)*(state[2]*state0[0]) - C*(state[2]*state0[3]))/(r*r*r);
			stm(5,4) = r0/mu*(state[5] - state0[5])*(state[4] - state0[4]) + (r0*(1.0 - F)*(state[2]*state0[1]) - C*(state[2]*state0[4]))/(r*r*r);
			stm(5,5) = r0/mu*(state[5] - state0[5])*(state[5] - state0[5]) + (r0*(1.0 - F)*(state[2]*state0[2]) - C*(state[2]*state0[5]))/(r*r*r) + Gt;

			//Step 5.4 compute Ftt and Gtt
			double rdot = r0 * U0dot + sigma0 * U1dot + U2dot;
			double r2 = r*r;
			Ftt = -sqmu / r0 * (U1dot / r - U1 * rdot / r2);

			Gtt = -(U2dot / r - U2 * rdot / r2);
		}
	}

	void KeplerLagrangeLaguerreConway(const double* state0, double* state, const double& mu, const double& propTime, double* F, double* G, double* Ft, double* Gt)
	{
		double r, sqrta, sqrtma, deltaE, deltaH;
		int iteration_count = 0;

		//set Laguerre-Conway n
		int n = 5;

		//step 1: determine if we are an ellipse or a hyperbola/parabola
		double sqmu = sqrt(mu);
		double r0 = sqrt(state0[0] * state0[0] + state0[1] * state0[1] + state0[2] * state0[2]);
		double v0 = sqrt(state0[3] * state0[3] + state0[4] * state0[4] + state0[5] * state0[5]);

		//Battin equation 4.23
		double sigma0 = (state0[0] * state0[3] + state0[1] * state0[4] + state0[2] * state0[5]) / sqmu;

		double a = 1.0 / (2.0 / r0 - v0*v0 / mu);
		double coeff = (1 - r0 / a);
		

		 
		if (a > 0.0) //ellipse
		{
			//step 2: solve Kepler's equation

			//Step 2.1: initialize ellipse-specific quantities
			sqrta = sqrt(a);

			//Battin equation 4.35
			double deltaM = sqmu / sqrt(a*a*a) * propTime;

			//Step 2.2: solve for deltaE using the Laguerre-Conway method
			//initial guess is deltaE = deltaM
			deltaE = 1000.0;
			double deltaE_new = deltaM;

			while (fabs(deltaE - deltaE_new) > 1.0e-8 && iteration_count < 100)
			{
				//iteration count update
				++iteration_count;

				//deltaE update
				deltaE = deltaE_new;

				//trig evaluations that only need to be done once
				double cdeltaE = cos(deltaE);
				double sdeltaE = sin(deltaE);

				//Kepler's equation for an elliptical orbit
				double F = deltaE - deltaM + sigma0 / sqrta * (1.0 - cdeltaE) - coeff * sdeltaE;

				//derivative with respect to deltaE
				double dF = 1 + sigma0 / sqrta * sdeltaE - coeff * cdeltaE;

				//second derivative with respect to deltaE
				double ddF = sigma0 / sqrta * cdeltaE + coeff * sdeltaE;

				//Laguerre-Conway update
				int sgn = dF >= 0 ? 1 : -1;
				deltaE_new = deltaE - n*F / (dF + sgn * sqrt(fabs( (n-1)*(n-1)*dF*dF - n * (n - 1) * F * ddF )));
			}

			//Step 3: find F, G, r Ft, Gt

			//Step 3.1 find F and G
			double cdeltaE = cos(deltaE);
			double sdeltaE = sin(deltaE);
			//Battin equations 4.41
			*F = 1 - a/r0 * (1.0 - cdeltaE);

			*G = a*sigma0/sqmu * (1.0 - cdeltaE) + r0*sqrta/sqmu * sdeltaE;

			//Step 3.2 find r
			//Battin equation 4.42
			r = a + (r0 - a) * cdeltaE + sigma0*sqrta * sdeltaE;

			//Step 3.3: find Ft, Gt
			//Battin equations 4.41
			*Ft = -sqrta*sqmu / (r*r0) * sdeltaE;

			*Gt = 1.0 - a/r * (1.0 -cdeltaE);
		}
		else //parabola or hyperbola
		{
			//Step 2.1: initialize hyperbola-specific quantities
			sqrtma = sqrt(-a);

			//Battin equation 4.51
			double deltaN = sqmu / sqrt(-a*a*a) * propTime;

			//Step 2.2: solve for deltaH using the Laguerre-Conway method
			//initial guess is deltaH = deltaN
			deltaH = 0;
			double deltaH_new = propTime < 0.0 ? -1 : 1;

			while (fabs(deltaH - deltaH_new) > 1.0e-8 && iteration_count < 100)
			{
				//iteration count update
				++iteration_count;

				//deltaE update
				deltaH = deltaH_new;

				//trig evaluations that only need to be done once
				double cdeltaH = cosh(deltaH);
				double sdeltaH = sinh(deltaH);

				//Kepler's equation for a hyperbolic orbit
				double F = -deltaN - deltaH + sigma0/sqrtma * (cdeltaH - 1.0) + coeff * sdeltaH;

				//derivative with respect to deltaH
				double dF = -1 + sigma0/sqrtma * sdeltaH + coeff * cdeltaH;

				//second derivative with respect to deltaH
				double ddF = sigma0/sqrtma * cdeltaH + coeff * sdeltaH;

				//Laguerre-Conway update
				int sgn = dF >= 0 ? 1 : -1;
				deltaH_new = deltaH - n*F / (dF + sgn * sqrt(fabs( (n-1)*(n-1)*dF*dF - n * (n - 1) * F * ddF )));
			}



			//Step 3: find F, G, r Ft, Gt

			//Step 3.1 find F and G
			double chdeltaH = cosh(deltaH);
			double shdeltaH = sinh(deltaH);
			//Battin equations 4.62
			*F = 1.0 - a/r0 * (1.0 - chdeltaH);

			*G = a*sigma0/sqmu * (1.0 - chdeltaH) + r0*sqrtma/sqmu * shdeltaH;

			//Step 3.2 find r
			//Battin equations 4.63
			//note: Battin writes r = -a + (r0 + a) * chdeltaH + sigma0*sqrtma * shdeltaH;
			//but this is incorrect - he left a sign change in the first and second terms
			r = a + (r0 - a) * chdeltaH + sigma0*sqrtma * shdeltaH;

			//Step 3.3: find Ft, Gt
			//Battin equations 4.62
			*Ft = -sqrtma*sqmu / (r*r0) * shdeltaH;

			*Gt = 1.0 - a/r * (1.0 - chdeltaH);
		}

		//Step 4: compute the final state as functions of F, G, Ft, Gt
		//Battin equation 3.33
		state[0] = *F*state0[0] + *G*state0[3];
		state[1] = *F*state0[1] + *G*state0[4];
		state[2] = *F*state0[2] + *G*state0[5];
		state[3] = *Ft*state0[0] + *Gt*state0[3];
		state[4] = *Ft*state0[1] + *Gt*state0[4];
		state[5] = *Ft*state0[2] + *Gt*state0[5];
	}

	//Lagrange-Laguerre-Conway Kepler solver
	double KeplerLaguerreConway(const double& ECC, const double& MN)
	{
		int iteration_count = 0;

		//set Laguerre-Conway n
		int n = 5;
		 
		if (ECC < 1.0) //ellipse
		{
			//Step 2: solve for E using the Laguerre-Conway method
			//initial guess is E = M
			double E = 1000.0;
			double E_new = MN;

			while (fabs(E - E_new) > 1.0e-12 && iteration_count < 100)
			{
				//iteration count update
				++iteration_count;

				//deltaE update
				E = E_new;

				//trig evaluations that only need to be done once
				double cE = cos(E);
				double sE = sin(E);

				//Kepler's equation for an elliptical orbit
				double F = E - ECC*sE - MN;

				//derivative with respect to E
				double dF = -ECC*cE + 1;

				//second derivative with respect to E
				double ddF = ECC*sE;

				//Laguerre-Conway update
				int sgn = dF >= 0 ? 1 : -1;
				E_new = E - n*F / (dF + sgn * sqrt(fabs( (n-1)*(n-1)*dF*dF - n * (n - 1) * F * ddF )));
			}

			return E;
		}
		else //parabola or hyperbola
		{
			//Step 2: solve for H using the Laguerre-Conway method
			double H = 1000;
			double H_new = 1;

			while (fabs(H - H_new) > 1.0e-12 && iteration_count < 100)
			{
				//iteration count update
				++iteration_count;

				//deltaE update
				H = H_new;

				//trig evaluations that only need to be done once
				double cH = cosh(H);
				double sH = sinh(H);

				//Kepler's equation for a hyperbolic orbit
				double F = H - ECC*sH - MN;

				//derivative with respect to deltaH
				double dF = -ECC*cH + 1;

				//second derivative with respect to deltaH
				double ddF = -ECC*sH;

				//Laguerre-Conway update
				int sgn = dF >= 0 ? 1 : -1;
				H_new = H - n*F / (dF + sgn * sqrt(fabs( (n-1)*(n-1)*dF*dF - n * (n - 1) * F * ddF )));
			}

			return H;
		}
	}
}//end namespace Kepler