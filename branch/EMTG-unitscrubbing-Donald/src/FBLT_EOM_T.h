#include <iostream>
#include <iomanip>
#include "SpiceUsr.h"
#include "EMTG_math.h"
#include "EMTG_Matrix.h"

#include "GSAD.h"

#ifndef _FBLT_EOM_T
#define _FBLT_EOM_T

//Function prototypes
namespace EMTG
{
	namespace Astrodynamics
	{
		namespace EOM
		{
			template <typename T>
			void FBLT_EOM(std::vector<T> & X, //state vector
						std::vector<T> & dXdTOF, //derivative of state vector with respect to time of flight, where time of flight is a decision variable
						std::vector<T> & f, //gradient of the state vector (instantaneous time derivative of the state vector)
						std::vector<T> & dfdTOF, //instantaneous derivative of state vector time derivative with respect to time of flight, where time of flight is a decision variable
						T & t_left_step, //epoch of the left hand side of the current RK sub-step
						T & dt_left_stepdTOF, //derivative of the left hand side epoch w.r.t. TOF
						const double & c, //current RK8(7)13M node
						T & h, //integrator step size
						T & dhdTOF, //partial derivative of the integrator step size w.r.t. TOF
						const double & launch_epoch,
						const double & mu_3B,
						const double & mu_cb, 
						const double & DU,
						const double & TU)
			{

				std::cout << std::setprecision(16);
				//STM size
				int rows = 11;
				int columns = 11;

				//cbtosc: spacecraft w.r.t. central body
				//cbto3B: 3rd body w.r.t. central body
				//3Btosc: spacecraft w.r.t. 3rd body
				T rvec_cbtosc[3] = { X[0], X[1], X[2] };
				T rvec_cbto3B[3];
				T rvec_3Btosc[3];
				T vvec_cbtosc[3] = { X[3], X[4], X[5] };
				T vvec_cbto3B[3];
				T vvec_3Btosc[3];

				//Locate a 3rd Body 599:Jupiter
				int spice_ID = 599;
				int central_body_spice_ID = 10; //the Sun
				double LT_dump;
				double state3B[6];
				double current_epoch;

				
				current_epoch = t_left_step + c*h;


				spkez_c(spice_ID, current_epoch*TU - (51544.5 * 86400.0), "J2000", "NONE", central_body_spice_ID, state3B, &LT_dump);

				//break the SPICE output state into r and v vectors
				for (size_t i = 0; i < 3; ++i)
				{
					rvec_cbto3B[i] = state3B[i] / DU;
					vvec_cbto3B[i] = state3B[i + 3] / DU * TU;
				}

				//derivative of current epoch w.r.t. TOF
				T dcurrent_epochdTOF = dt_left_stepdTOF + c*dhdTOF;

				//Third body position sensitivity w.r.t. TOF
				//Just the derivative of the SPICE position Chebychev polynomial
				//If we are using SPICE to locate small bodies or satellites we will need to find a way to obtain the position derivatives.....Ryan Russell's FIRE ephemeris reader would fix this
				T dx3BdTOF = vvec_cbto3B[0] * dcurrent_epochdTOF;
				T dy3BdTOF = vvec_cbto3B[1] * dcurrent_epochdTOF;
				T dz3BdTOF = vvec_cbto3B[2] * dcurrent_epochdTOF;

				
				//Form the spacecraft position and velocity vectors w.r.t. the third body
				rvec_3Btosc[0] = rvec_cbtosc[0] - rvec_cbto3B[0];
				rvec_3Btosc[1] = rvec_cbtosc[1] - rvec_cbto3B[1];
				rvec_3Btosc[2] = rvec_cbtosc[2] - rvec_cbto3B[2];
				vvec_3Btosc[0] = vvec_cbtosc[0] - vvec_cbto3B[0];
				vvec_3Btosc[1] = vvec_cbtosc[1] - vvec_cbto3B[1];
				vvec_3Btosc[2] = vvec_cbtosc[2] - vvec_cbto3B[2];

				//define vector magnitudes for the three positions and velocities of interest in the restricted three-body problem
				T r_cbtosc = sqrt(rvec_cbtosc[0] * rvec_cbtosc[0] + rvec_cbtosc[1] * rvec_cbtosc[1] + rvec_cbtosc[2] * rvec_cbtosc[2]);
				T r_cbto3B = sqrt(rvec_cbto3B[0] * rvec_cbto3B[0] + rvec_cbto3B[1] * rvec_cbto3B[1] + rvec_cbto3B[2] * rvec_cbto3B[2]);
				T r_3Btosc = sqrt(rvec_3Btosc[0] * rvec_3Btosc[0] + rvec_3Btosc[1] * rvec_3Btosc[1] + rvec_3Btosc[2] * rvec_3Btosc[2]);
				T v_cbtosc = sqrt(vvec_cbtosc[0] * vvec_cbtosc[0] + vvec_cbtosc[1] * vvec_cbtosc[1] + vvec_cbtosc[2] * vvec_cbtosc[2]);
				T v_cbto3B = sqrt(vvec_cbto3B[0] * vvec_cbto3B[0] + vvec_cbto3B[1] * vvec_cbto3B[1] + vvec_cbto3B[2] * vvec_cbto3B[2]);
				T v_3Btosc = sqrt(vvec_3Btosc[0] * vvec_3Btosc[0] + vvec_3Btosc[1] * vvec_3Btosc[1] + vvec_3Btosc[2] * vvec_3Btosc[2]);

				//create some additional local variables from the state vector and the magnitude of the control vector
				T m_sc = X[6];
				T u[3] = { X[7], X[8], X[9] };
				T unorm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

				//Initialize STM equation components -- Phidot = A*Phi
				EMTG::math::Matrix<double> A(rows, columns, 0.0);
				EMTG::math::Matrix<double> Phi(rows, columns, 0.0);
				EMTG::math::Matrix<double> Phidot(rows, columns, 0.0);


				//*************************************
				//
				//Spacecraft Parameters and Power Model
				//
				//*************************************
				double g0 = 9.80665 / 1000.0 / DU * TU * TU; //acceleration due to Earth's gravity
				double D = 0.95; //duty cycle
				double thrust = 0.472 / 1000.0 / DU * TU * TU; //thrust
				double Isp = 4336.0 / TU; //specific impulse
				double mdotmax = thrust / Isp / g0; //max mass flow rate 

				//*******************
				//
				//Equations of Motion
				//
				//*******************

				//Spacecraft velocity
				f[0] = X[3];
				f[1] = X[4];
				f[2] = X[5];

				//Effects of the central body on acceleration
				T mu_over_r_cbtosc_cubed = mu_cb / (r_cbtosc*r_cbtosc*r_cbtosc);

				f[3] = -mu_over_r_cbtosc_cubed*rvec_cbtosc[0];
				f[4] = -mu_over_r_cbtosc_cubed*rvec_cbtosc[1];
				f[5] = -mu_over_r_cbtosc_cubed*rvec_cbtosc[2];

				//Effects of the third body on acceleration
				T one_over_r3Btosc_cubed = 1.0 / (r_3Btosc*r_3Btosc*r_3Btosc);
				T one_over_rcbto3B_cubed = 1.0 / (r_cbto3B*r_cbto3B*r_cbto3B);

				f[3] -= mu_3B * (rvec_3Btosc[0] * one_over_r3Btosc_cubed + rvec_cbto3B[0] * one_over_rcbto3B_cubed);
				f[4] -= mu_3B * (rvec_3Btosc[1] * one_over_r3Btosc_cubed + rvec_cbto3B[1] * one_over_rcbto3B_cubed);
				f[5] -= mu_3B * (rvec_3Btosc[2] * one_over_r3Btosc_cubed + rvec_cbto3B[2] * one_over_rcbto3B_cubed);

				//Effects of the thruster on acceleration
				T D_thrust_over_msc = D * thrust / m_sc;

				f[3] += u[0] * D_thrust_over_msc;
				f[4] += u[1] * D_thrust_over_msc;
				f[5] += u[2] * D_thrust_over_msc;


				//mass
				f[6] = -unorm*D*mdotmax;

				//control vector time rates of change (these are decision variables and are not impacted by dynamics and are constant across an FBLT time step)
				//these are just a place holders in the STM and will never have a diffeqs associated with them
				f[7] = 0.0000000000000000;
				f[8] = 0.0000000000000000;
				f[9] = 0.0000000000000000;

				//TOF time rates of change (this is also a decision variable; it is not impacted by dynamics and is constant across an FBLT time step/phase)
				//this is just a place holder in the STM and will never have a diffeq associated with it
				f[10] = 0.0000000000000000;

				

				//***********************
				//
				//Power Model Derivatives
				//
				//***********************
				T dTdP = 0.0000000000000000;
				T dPdr = 0.0000000000000000;
				T dmdotmaxdP = 0.0000000000000000;
				T dmdotmaxdTOF = 0.0000000000000000;
				T dthrustdTOF = 0.0000000000000000;

				//****************************************
				//
				//Other Auxilliary Derivatives of Interest
				//
				//****************************************
				double epsilon = 1.0e-10; //avoid divide by zero if thruster is off

				T one_over_rcbtosc = 1.0 / r_cbtosc;
				T one_over_unorm_plus_epsilon = 1.0 / (unorm + epsilon);

				//gradient of the magnitude of the spacecraft position vector
				T drdx = rvec_cbtosc[0] * one_over_rcbtosc;
				T drdy = rvec_cbtosc[1] * one_over_rcbtosc;
				T drdz = rvec_cbtosc[2] * one_over_rcbtosc;

				//gradient of the magnitude of the control vector
				T dunormdux = u[0] * one_over_unorm_plus_epsilon;
				T dunormduy = u[1] * one_over_unorm_plus_epsilon;
				T dunormduz = u[2] * one_over_unorm_plus_epsilon;


				//********************
				//
				//A matrix calculation
				//
				//********************

				//Top Row Identity
				A(0, 3) = 1.0;
				A(1, 4) = 1.0;
				A(2, 5) = 1.0;

				//A21 dadr
				T three_mucb_over_rcbtosc5 = 3.0 * mu_cb / (r_cbtosc * r_cbtosc * r_cbtosc * r_cbtosc * r_cbtosc);
				T mucb_over_rcbtosc3 = mu_cb / (r_cbtosc * r_cbtosc * r_cbtosc);
				T three_mu3B_over_r3Btosc5 = 3.0 * mu_3B / (r_3Btosc * r_3Btosc * r_3Btosc * r_3Btosc * r_3Btosc);
				T mu3B_over_r3Btosc3 = mu_3B / (r_3Btosc * r_3Btosc * r_3Btosc);
				T D_dTdP_dPdr_over_msc = D * dTdP * dPdr / m_sc;

				A(3, 0) = three_mucb_over_rcbtosc5 * rvec_cbtosc[0] * rvec_cbtosc[0] - mucb_over_rcbtosc3 + three_mu3B_over_r3Btosc5 * rvec_3Btosc[0] * rvec_3Btosc[0] - mu3B_over_r3Btosc3 + u[0] * D_dTdP_dPdr_over_msc * drdx;
				A(3, 1) = three_mucb_over_rcbtosc5 * rvec_cbtosc[0] * rvec_cbtosc[1]                      + three_mu3B_over_r3Btosc5 * rvec_3Btosc[0] * rvec_3Btosc[1]                      + u[0] * D_dTdP_dPdr_over_msc * drdy;
				A(3, 2) = three_mucb_over_rcbtosc5 * rvec_cbtosc[0] * rvec_cbtosc[2]                      + three_mu3B_over_r3Btosc5 * rvec_3Btosc[0] * rvec_3Btosc[2]                      + u[0] * D_dTdP_dPdr_over_msc * drdz;
				A(4, 0) = three_mucb_over_rcbtosc5 * rvec_cbtosc[1] * rvec_cbtosc[0]                      + three_mu3B_over_r3Btosc5 * rvec_3Btosc[1] * rvec_3Btosc[0]                      + u[1] * D_dTdP_dPdr_over_msc * drdx;
				A(4, 1) = three_mucb_over_rcbtosc5 * rvec_cbtosc[1] * rvec_cbtosc[1] - mucb_over_rcbtosc3 + three_mu3B_over_r3Btosc5 * rvec_3Btosc[1] * rvec_3Btosc[1] - mu3B_over_r3Btosc3 + u[1] * D_dTdP_dPdr_over_msc * drdy;
				A(4, 2) = three_mucb_over_rcbtosc5 * rvec_cbtosc[1] * rvec_cbtosc[2]                      + three_mu3B_over_r3Btosc5 * rvec_3Btosc[1] * rvec_3Btosc[2]                      + u[1] * D_dTdP_dPdr_over_msc * drdz;
				A(5, 0) = three_mucb_over_rcbtosc5 * rvec_cbtosc[2] * rvec_cbtosc[0]                      + three_mu3B_over_r3Btosc5 * rvec_3Btosc[2] * rvec_3Btosc[0]                      + u[2] * D_dTdP_dPdr_over_msc * drdx;
				A(5, 1) = three_mucb_over_rcbtosc5 * rvec_cbtosc[2] * rvec_cbtosc[1]                      + three_mu3B_over_r3Btosc5 * rvec_3Btosc[2] * rvec_3Btosc[1]                      + u[2] * D_dTdP_dPdr_over_msc * drdy;
				A(5, 2) = three_mucb_over_rcbtosc5 * rvec_cbtosc[2] * rvec_cbtosc[2] - mucb_over_rcbtosc3 + three_mu3B_over_r3Btosc5 * rvec_3Btosc[2] * rvec_3Btosc[2] - mu3B_over_r3Btosc3 + u[2] * D_dTdP_dPdr_over_msc * drdz;
				
				//A23 dadm
				T D_thrust_over_msc2 = D * thrust / (m_sc*m_sc);

				A(3, 6) = -u[0] * D_thrust_over_msc2;
				A(4, 6) = -u[1] * D_thrust_over_msc2;
				A(5, 6) = -u[2] * D_thrust_over_msc2;
				
				//A24 dadu
				A(3, 7) = D_thrust_over_msc;
				A(3, 8) = 0.0;
				A(3, 9) = 0.0;
				A(4, 7) = 0.0;
				A(4, 8) = D_thrust_over_msc;
				A(4, 9) = 0.0;
				A(5, 7) = 0.0;
				A(5, 8) = 0.0;
				A(5, 9) = D_thrust_over_msc;

				//A31 dmdotdr
				T unorm_D_dmdotmaxdP_dPdr = unorm * D * dmdotmaxdP;

				A(6, 0) = -unorm_D_dmdotmaxdP_dPdr * drdx;
				A(6, 1) = -unorm_D_dmdotmaxdP_dPdr * drdy;
				A(6, 2) = -unorm_D_dmdotmaxdP_dPdr * drdz;

				//A34 dmdotdu
				T D_mdotmax = D * mdotmax;

				A(6, 7) = -dunormdux * D_mdotmax;
				A(6, 8) = -dunormduy * D_mdotmax;
				A(6, 9) = -dunormduz * D_mdotmax;


				//****************************************************
				//
				//TOF derivative calculation for this DOPRI 8(7) stage
				//
				//****************************************************

				//dadr3B -- sensitivity of s/c acceleration to third body position
				T dadr3B[9];

				//daxdr3B
				dadr3B[0] = -mu_3B * (3.0 * rvec_3Btosc[0] * rvec_3Btosc[0] / pow(r_3Btosc, 5.0) - 1.0 / pow(r_3Btosc, 3.0) - 3.0 * rvec_cbto3B[0] * rvec_cbto3B[0] / pow(r_cbto3B, 5.0) + 1.0 / pow(r_cbto3B, 3.0));
				dadr3B[1] = -mu_3B * (3.0 * rvec_3Btosc[0] * rvec_3Btosc[1] / pow(r_3Btosc, 5.0)                            - 3.0 * rvec_cbto3B[0] * rvec_cbto3B[1] / pow(r_cbto3B, 5.0));
				dadr3B[2] = -mu_3B * (3.0 * rvec_3Btosc[0] * rvec_3Btosc[2] / pow(r_3Btosc, 5.0)                            - 3.0 * rvec_cbto3B[0] * rvec_cbto3B[2] / pow(r_cbto3B, 5.0));

				//daydr3B
				dadr3B[3] = -mu_3B * (3.0 * rvec_3Btosc[1] * rvec_3Btosc[0] / pow(r_3Btosc, 5.0)                            - 3.0 * rvec_cbto3B[1] * rvec_cbto3B[0] / pow(r_cbto3B, 5.0));
				dadr3B[4] = -mu_3B * (3.0 * rvec_3Btosc[1] * rvec_3Btosc[1] / pow(r_3Btosc, 5.0) - 1.0 / pow(r_3Btosc, 3.0) - 3.0 * rvec_cbto3B[1] * rvec_cbto3B[1] / pow(r_cbto3B, 5.0) + 1.0 / pow(r_cbto3B, 3.0));
				dadr3B[5] = -mu_3B * (3.0 * rvec_3Btosc[1] * rvec_3Btosc[2] / pow(r_3Btosc, 5.0)                            - 3.0 * rvec_cbto3B[1] * rvec_cbto3B[2] / pow(r_cbto3B, 5.0));

				//dazdr3B
				dadr3B[6] = -mu_3B * (3.0 * rvec_3Btosc[2] * rvec_3Btosc[0] / pow(r_3Btosc, 5.0)                            - 3.0 * rvec_cbto3B[2] * rvec_cbto3B[0] / pow(r_cbto3B, 5.0));
				dadr3B[7] = -mu_3B * (3.0 * rvec_3Btosc[2] * rvec_3Btosc[1] / pow(r_3Btosc, 5.0)                            - 3.0 * rvec_cbto3B[2] * rvec_cbto3B[1] / pow(r_cbto3B, 5.0));
				dadr3B[8] = -mu_3B * (3.0 * rvec_3Btosc[2] * rvec_3Btosc[2] / pow(r_3Btosc, 5.0) - 1.0 / pow(r_3Btosc, 3.0) - 3.0 * rvec_cbto3B[2] * rvec_cbto3B[2] / pow(r_cbto3B, 5.0) + 1.0 / pow(r_cbto3B, 3.0));

				//Construct the partial derivatives w.r.t. TOF for the state gradient
				T dadthrust[3];

				dadthrust[0] = u[0] * D / m_sc;
				dadthrust[1] = u[1] * D / m_sc;
				dadthrust[2] = u[2] * D / m_sc;

				dfdTOF[0] = dXdTOF[3];
				dfdTOF[1] = dXdTOF[4];
				dfdTOF[2] = dXdTOF[5];
				//dxddotdTOF = dxddotdxsc*dxscdTOF + dxddotdysc*dyscdTOF +  dxddotdzsc*dzscdTOF + dxddotdx3B*dx3BdTOF + dxddotdy3B*dy3BdTOF + dxddotdz3B*dx3BdTOF + dxddotdmsc*dmscdTOF + dxddotdTx*dTxdTOF
				dfdTOF[3] = A(3, 0)*dXdTOF[0] + A(3, 1)*dXdTOF[1] + A(3, 2)*dXdTOF[2] + dadr3B[0] * dx3BdTOF + dadr3B[1] * dy3BdTOF + dadr3B[2] * dz3BdTOF + A(3, 6)*dXdTOF[6] + dadthrust[0] * dthrustdTOF;
				dfdTOF[4] = A(4, 0)*dXdTOF[0] + A(4, 1)*dXdTOF[1] + A(4, 2)*dXdTOF[2] + dadr3B[3] * dx3BdTOF + dadr3B[4] * dy3BdTOF + dadr3B[5] * dz3BdTOF + A(4, 6)*dXdTOF[6] + dadthrust[1] * dthrustdTOF;
				dfdTOF[5] = A(5, 0)*dXdTOF[0] + A(5, 1)*dXdTOF[1] + A(5, 2)*dXdTOF[2] + dadr3B[6] * dx3BdTOF + dadr3B[7] * dy3BdTOF + dadr3B[8] * dz3BdTOF + A(5, 6)*dXdTOF[6] + dadthrust[2] * dthrustdTOF;
				dfdTOF[6] = -unorm * D * dmdotmaxdTOF;

				//*******************
				//
				//Phi dot calculation
				//
				//*******************

				//Form the Phi matrix
				//The Phi matrix is comprised of state variables
				int xcount = 11;

				for (size_t i = 0; i < rows; ++i)
				{
					for (size_t j = 0; j < columns; ++j)
					{
						Phi(i, j) = X[xcount];
						++xcount;
					}
				}

				//differential equation for STM creation
				Phidot = A*Phi;
				
				//Package the Phidot entries into the gradient vector f behind the spacecraft's state variables as well as the control and TOF gradients which are always zero
				int fcount = 11; //STM entry diffeq's are placed behind the augmented state vector
				for (size_t i = 0; i < rows; ++i)
				{
					for (size_t j = 0; j < columns; ++j)
					{
						f[fcount] = Phidot(i, j);
						++fcount;
					}
				}

			}//end FBLT_EOM

		}
	}
}

#endif