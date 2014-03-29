//Laguerre-Conway Kepler solver using the iterative n method of Ger
//from "Orbital Mechanics" by Prussing and Conway, Chapter 2
//and "Classical and Advanced Kepler Algorithms" by Gim J. Der
//STMs by Ellison

#include "STM.h"
#include "Kepler_Lagrange_Laguerre_Conway_Der.h"

#include <math.h>
#include <iostream>

namespace Kepler
{
	void Kepler_Lagrange_Laguerre_Conway_Der(const double* state0_kms,
											 double* state_kms,
											 const double& mu_km2s3,
											 const double& LU,
											 const double& propTime_seconds,
											 double& F,
											 double& G,
											 double& Ft,
											 double& Gt,
											 double& Ftt,
											 double& Gtt,
											 STM& stm,
											 const bool& compute_STM_flag)
	{
		//Step 0: declare necessary variables
		double mu = 1.0;
		double sqmu = sqrt(mu);
		const double alphatol = 1.0e-12;
		const double Xtol = 1.0e-12;
		const int maximum_order = 10;
		const int maximum_iterations_per_order = 10;
		double r, sigma, U0, U1, U2, U3, X = 1.0e+100, X_new, dX;
		double sqalpha, sqmalpha, sqmalphaX;
		double state0[6], state[6];

		//Step 1: compute r0 and v0 in LUs and TUs, and propTime in TUs
		double TU = sqrt(LU*LU*LU/mu_km2s3);

		for (int k = 0; k < 3; ++k)
			state0[k] = state0_kms[k] / LU;
		for (int k = 3; k < 6; ++k)
			state0[k] = state0_kms[k] / LU * TU;


		double r0 = sqrt(state0[0] * state0[0] + state0[1] * state0[1] + state0[2] * state0[2]);
		double v0 = sqrt(state0[3] * state0[3] + state0[4] * state0[4] + state0[5] * state0[5]);
		double propTime = propTime_seconds / TU;

		//Step 2: compute alpha, which determines whether we are on an ellipse, parabola, or hyperbola
		//and also sigma0 which we need in the universal Kepler iteration
		double alpha = (2.0 / r0 - v0*v0);
		double sigma0 = (state0[0] * state0[3] + state0[1] * state0[4] + state0[2] * state0[5]);
		if (alpha > alphatol) //ellipse
			sqalpha = sqrt(alpha);
		else if (alpha < -alphatol)
			sqmalpha = sqrt(-alpha);

		//Step 3: compute an initial guess
		if (alpha > alphatol)
			X_new = alpha * mu * propTime;
		else
			X_new = 0.1 * mu * propTime / r0;

		//Step 4: Perform the Laguerre-Conway-Ger iteration
		int N = 2; //current Laguerre iteration order
		int iteration_this_N = 0; //number of iterations for the current order
		int total_iterations = 0;
		while (fabs(X - X_new) > Xtol && N < maximum_order)
		{
			//Step 4.1: increment the iteration count
			//if we have maxed out the number of iterations for this order, increase the order
			++iteration_this_N;
			++total_iterations;
			if (iteration_this_N >= maximum_iterations_per_order)
			{
				++N; //increase the order
				iteration_this_N = 0; //reset the number of iterations
			}

			//Step 4.2: update X
			X = X_new;

			//Step 4.3: compute U0, U1, U2, U3 for the candidate X
			if (alpha > alphatol) //ellipse
			{
				//compute U0, U1, U2, U3 via Stumpff functions
				double y = alpha * X * X;
				double C = (1.0 - cos(sqrt(y))) / y;
				double S = (sqrt(y) - sin(sqrt(y))) / sqrt(y*y*y);
				U1 = X * (1.0 - y * S);
				U2 = X * X * C;
				U3 = X * X * X * S;
				U0 = 1.0 - alpha * U2;
			}
			else if (alpha < -alphatol) //hyperbola
			{
				
				sqmalphaX = sqmalpha * X;
				if (sqmalphaX > 50.0)
					sqmalphaX = 50.0;
				else if (sqmalphaX < -50.0)
					sqmalphaX = -50.0;
				U0 = cosh(sqmalphaX);
				U1 = sinh(sqmalphaX) / sqmalpha;
				U2 = (1.0 - U0) / alpha;
				U3 = 1.0 / alpha * (X - U1);
				/*double y = alpha * X * X;
				double yqr = sqrt(-y);
				double C = (1.0 - cosh(yqr)) / y;
				double S = (sinh(yqr) - sqrt(-y)) / (-yqr*yqr*yqr);
				U1 = X * (1.0 - y * S);
				U2 = X * X * C;
				U3 = X * X * X * S;
				U0 = 1.0 - alpha * U2;*/
			}
			else //parabola
			{
				U0 = 1.0;
				U1 = X;
				U2 = U1*X / 2.0;
				U3 = U2*X / 3.0;
			}
			


			//Step 4.4: compute r and sigma
			r = r0 * U0 + sigma0 * U1 + U2;
			sigma = sigma0 * U0 + (1 - alpha * r0) * U1;

			//Step 4.5 compute F, dF, and ddF
			double FX = r0 * U1 + sigma0 * U2 + U3 - sqmu * propTime;
			double dFX = r;
			double ddFX = sigma;

			//Step 4.6 Laguerre-Conway or Newton iteration depending on the situation
			int sgn = dFX >= 0 ? 1 : -1;
			double denom = fabs( (N-1)*(N-1)*dFX*dFX - N * (N - 1) * FX * ddFX );
			if (denom > 0.0) //Laguerre-Conway iteration
			{
				dX = N*FX / (dFX + sgn * sqrt(denom));	
			}
			else //Newton iteration for very sensitive cases, i.e. when on a hyperbolic/parabolic asymptote
				dX = FX / dFX;

			X_new -= dX;
		}

		//Step 5: find F, G, Ft, Gt
		F = 1.0 - U2 / r0;
		G = (r0 * U1 + sigma0 * U2) / sqmu;
		Ft = -sqmu / (r0 * r) * U1;
		Gt = 1.0 - U2 / r;

		//Step 6: compute the final state as functions of F, G, Ft, Gt
		state[0] = (F*state0[0] + G*state0[3]);
		state[1] = (F*state0[1] + G*state0[4]);
		state[2] = (F*state0[2] + G*state0[5]);
		state[3] = (Ft*state0[0] + Gt*state0[3]);
		state[4] = (Ft*state0[1] + Gt*state0[4]);
		state[5] = (Ft*state0[2] + Gt*state0[5]);

		for (int k = 0; k < 3; ++k)
			state_kms[k] = state[k] * LU;
		for (int k = 3; k < 6; ++k)
			state_kms[k] = state[k] * LU / TU;

		//Step 7: Compute the state transition matrix if requested
		if (compute_STM_flag)
		{
			//STM calculation code
			//Step 7.1 compute the universal functions Ui and their derivatives
			//from the recursion relation 4.76
			double a;
			if (fabs(alpha) < alphatol)
				a = 1.0e+30;
			else
				a = 1.0 / alpha;

			double U4 = a * (X*X / 2.0 - U2);
			double U5 = a * (X*X*X / 6.0 - U3);

			//derivatives
			double dXdt = sqmu / r;
			double U0dot = -alpha * U1 * dXdt;
			double U1dot = U0 * dXdt;
			double U2dot = U1 * dXdt;

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

			//scale the STM
			//lower left
			for (int i = 3; i < 6; ++i)
				for (int j = 0; j < 3; ++j)
					stm(i,j) /= TU;
			//upper right
			for (int i = 0; i < 3; ++i)
				for (int j = 3; j < 6; ++j)
					stm(i,j) *= TU;

			//Step 5.4 compute Ftt and Gtt
			double rdot = r0 * U0dot + sigma0 * U1dot + U2dot;
			double r2 = r*r;
			Ftt = -sqmu / r0 * (U1dot / r - U1 * rdot / r2);

			Gtt = -(U2dot / r - U2 * rdot / r2);

			//Step 5.4 scale F, G, Ft, Gt, Ftt, Gtt
			//F converts km to km so there is no scalig
			//G converts km/s to km so it must be multiplied by TU
			G *= TU;
			//Ft converts km to km/s so it must be divided by TU
			Ft /= TU;
			//Gt converts km/s to km/s so there is no scaling
			//Ftt converts km to km/s^2 so it must be divided by TU twice
			Ftt /= (TU * TU);
			//Gtt converts km/s to km/s^2 so it must be divided by TU
			Gtt /= TU;
		}
	}
}