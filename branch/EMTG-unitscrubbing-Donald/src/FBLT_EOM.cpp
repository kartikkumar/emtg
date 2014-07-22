#include <iostream>
#include <iomanip>
#include "GSAD.h"
#include "SpiceUsr.h"
#include "EMTG_math.h"
#include "EMTG_Matrix.h"



namespace EMTG {namespace Astrodynamics {namespace EOM {

	void FBLT_EOM(GSAD::adouble * x, GSAD::adouble * f, const double & current_epoch, double & mu_3B, double & mu_cb, double & DU, double & TU)
	{
		bool DEBUG = false;
		std::cout << std::setprecision(16);
		//STM size
		int rows = 10;
		int columns = 10;

		//cbtosc: spacecraft w.r.t. central body
		//cbto3B: 3rd body w.r.t. central body
		//3Btosc: spacecraft w.r.t. 3rd body
		GSAD::adouble rvec_cbtosc[3] = { x[0], x[1], x[2] };
		GSAD::adouble rvec_cbto3B[3];
		GSAD::adouble rvec_3Btosc[3];
		GSAD::adouble vvec_cbtosc[3] = { x[3], x[4], x[5] };
		GSAD::adouble vvec_cbto3B[3];
		GSAD::adouble vvec_3Btosc[3];

		//Locate a 3rd Body 599:Jupiter
		int spice_ID = 599;
		int central_body_spice_ID = 10; //the Sun
		double LT_dump;
		double state3B[6];

		spkez_c(spice_ID, current_epoch*TU - (51544.5 * 86400.0), "J2000", "NONE", central_body_spice_ID, state3B, &LT_dump);

		for (size_t i = 0; i < 3; ++i)
		{
			rvec_cbto3B[i] = state3B[i] / DU;
			vvec_cbto3B[i] = state3B[i + 3] / DU * TU;
		}

		//DEBUG
		if (DEBUG)
			std::cout << "Jupiter state: " << state3B[0] << " " << state3B[1] << " " << state3B[2] << " " << state3B[3] << " " << state3B[4] << " " << state3B[5] << std::endl << std::endl;

		//Form the spacecraft position and velocity vectors w.r.t. the third body
		rvec_3Btosc[0] = rvec_cbtosc[0] - rvec_cbto3B[0];
		rvec_3Btosc[1] = rvec_cbtosc[1] - rvec_cbto3B[1];
		rvec_3Btosc[2] = rvec_cbtosc[2] - rvec_cbto3B[2];
		vvec_3Btosc[0] = vvec_cbtosc[0] - vvec_cbto3B[0];
		vvec_3Btosc[1] = vvec_cbtosc[1] - vvec_cbto3B[1];
		vvec_3Btosc[2] = vvec_cbtosc[2] - vvec_cbto3B[2];

		//define some vector magnitudes
		GSAD::adouble r_cbtosc = sqrt(rvec_cbtosc[0] * rvec_cbtosc[0] + rvec_cbtosc[1] * rvec_cbtosc[1] + rvec_cbtosc[2] * rvec_cbtosc[2]);
		GSAD::adouble r_cbto3B = sqrt(rvec_cbto3B[0] * rvec_cbto3B[0] + rvec_cbto3B[1] * rvec_cbto3B[1] + rvec_cbto3B[2] * rvec_cbto3B[2]);
		GSAD::adouble r_3Btosc = sqrt(rvec_3Btosc[0] * rvec_3Btosc[0] + rvec_3Btosc[1] * rvec_3Btosc[1] + rvec_3Btosc[2] * rvec_3Btosc[2]);
		GSAD::adouble v_cbtosc = sqrt(vvec_cbtosc[0] * vvec_cbtosc[0] + vvec_cbtosc[1] * vvec_cbtosc[1] + vvec_cbtosc[2] * vvec_cbtosc[2]);
		GSAD::adouble v_cbto3B = sqrt(vvec_cbto3B[0] * vvec_cbto3B[0] + vvec_cbto3B[1] * vvec_cbto3B[1] + vvec_cbto3B[2] * vvec_cbto3B[2]);
		GSAD::adouble v_3Btosc = sqrt(vvec_3Btosc[0] * vvec_3Btosc[0] + vvec_3Btosc[1] * vvec_3Btosc[1] + vvec_3Btosc[2] * vvec_3Btosc[2]);

		GSAD::adouble m_sc = x[6];
		GSAD::adouble u[3] = { x[7], x[8], x[9] };
		GSAD::adouble unorm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

		//Initialize STM equation components Phidot = A*Phi
		EMTG::math::Matrix<GSAD::adouble> A(rows, columns, 0.0);
		EMTG::math::Matrix<GSAD::adouble> Phi(rows, columns, 0.0);
		EMTG::math::Matrix<GSAD::adouble> Phidot(rows, columns, 0.0);

		//*************************************
		//
		//Spacecraft Parameters and Power Model
		//
		//*************************************
		double g0 = 9.80665 / 1000.0 / DU * TU * TU; //acceleration due to Earth's gravity
		double D = 0.95; //duty cycle
		double T = 0.472 / 1000.0 / DU * TU * TU; //thrust
		double Isp = 4336.0 / TU; //specific impulse
		double mdotmax = T / Isp / g0; //max mass flow rate 

		//*******************
		//
		//Equations of Motion
		//
		//*******************

		//DEBUG
		if (DEBUG)
		{
			for (size_t i = 0; i < 107; ++i)
				std::cout << x[i] << std::endl;
			getchar();
		}

		//Spacecraft velocity
		f[0] = x[3];
		f[1] = x[4];
		f[2] = x[5];

		//Effects of the central body on acceleration
		f[3] = -mu_cb / (r_cbtosc*r_cbtosc*r_cbtosc)*rvec_cbtosc[0];
		f[4] = -mu_cb / (r_cbtosc*r_cbtosc*r_cbtosc)*rvec_cbtosc[1];
		f[5] = -mu_cb / (r_cbtosc*r_cbtosc*r_cbtosc)*rvec_cbtosc[2];

		//Effects of the third body on acceleration
		f[3] -= mu_3B * (rvec_3Btosc[0] / (r_3Btosc*r_3Btosc*r_3Btosc) + rvec_cbto3B[0] / (r_cbto3B*r_cbto3B*r_cbto3B));
		f[4] -= mu_3B * (rvec_3Btosc[1] / (r_3Btosc*r_3Btosc*r_3Btosc) + rvec_cbto3B[1] / (r_cbto3B*r_cbto3B*r_cbto3B));
		f[5] -= mu_3B * (rvec_3Btosc[2] / (r_3Btosc*r_3Btosc*r_3Btosc) + rvec_cbto3B[2] / (r_cbto3B*r_cbto3B*r_cbto3B));

		//Effects of the thruster on acceleration
		f[3] += u[0] * D * T / m_sc;
		f[4] += u[1] * D * T / m_sc;
		f[5] += u[2] * D * T / m_sc;


		//mass
		f[6] = -unorm*D*mdotmax;

		//control vector time rates of change (these are decision variables and are not impacted by dynamics and are constant across an FBLT time step)
		f[7] = 0.0000000000000000;
		f[8] = 0.0000000000000000;
		f[9] = 0.0000000000000000;

		//***********************
		//
		//Power Model Derivatives
		//
		//***********************
		GSAD::adouble dTdP = 0.0000000000000000;
		GSAD::adouble dPdr = 0.0000000000000000;
		GSAD::adouble dmdotmaxdP = 0.0000000000000000;

		//*****************
		//
		//Other Derivatives
		//
		//*****************
		double epsilon = 1.0e-10; //avoid divide by zero if thruster is off
		GSAD::adouble drdx = rvec_cbtosc[0] / r_cbtosc;
		GSAD::adouble drdy = rvec_cbtosc[1] / r_cbtosc;
		GSAD::adouble drdz = rvec_cbtosc[2] / r_cbtosc;
		GSAD::adouble dunormdux = u[0] / (unorm + epsilon);
		GSAD::adouble dunormduy = u[1] / (unorm + epsilon);
		GSAD::adouble dunormduz = u[2] / (unorm + epsilon);


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
		A(3, 0) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[0] * rvec_cbtosc[0] - mu_cb / pow(r_cbtosc, 3.0) + 3.0 * mu_3B*rvec_3Btosc[0] * rvec_3Btosc[0] / pow(r_3Btosc, 5.0) - mu_3B / pow(r_3Btosc, 3.0) + u[0] * D / m_sc*dTdP*dPdr*drdx;
		A(3, 1) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[0] * rvec_cbtosc[1]                              + 3.0 * mu_3B*rvec_3Btosc[0] * rvec_3Btosc[1] / pow(r_3Btosc, 5.0)								+ u[0] * D / m_sc*dTdP*dPdr*drdy;
		A(3, 2) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[0] * rvec_cbtosc[2]                              + 3.0 * mu_3B*rvec_3Btosc[0] * rvec_3Btosc[2] / pow(r_3Btosc, 5.0)								+ u[0] * D / m_sc*dTdP*dPdr*drdz;
		A(4, 0) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[1] * rvec_cbtosc[0]                              + 3.0 * mu_3B*rvec_3Btosc[1] * rvec_3Btosc[0] / pow(r_3Btosc, 5.0)								+ u[1] * D / m_sc*dTdP*dPdr*drdx;
		A(4, 1) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[1] * rvec_cbtosc[1] - mu_cb / pow(r_cbtosc, 3.0) + 3.0 * mu_3B*rvec_3Btosc[1] * rvec_3Btosc[1] / pow(r_3Btosc, 5.0) - mu_3B / pow(r_3Btosc, 3.0) + u[1] * D / m_sc*dTdP*dPdr*drdy;
		A(4, 2) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[1] * rvec_cbtosc[2]								+ 3.0 * mu_3B*rvec_3Btosc[1] * rvec_3Btosc[2] / pow(r_3Btosc, 5.0)								+ u[1] * D / m_sc*dTdP*dPdr*drdz;
		A(5, 0) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[2] * rvec_cbtosc[0]								+ 3.0 * mu_3B*rvec_3Btosc[2] * rvec_3Btosc[0] / pow(r_3Btosc, 5.0)								+ u[2] * D / m_sc*dTdP*dPdr*drdx;
		A(5, 1) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[2] * rvec_cbtosc[1]								+ 3.0 * mu_3B*rvec_3Btosc[2] * rvec_3Btosc[1] / pow(r_3Btosc, 5.0)								+ u[2] * D / m_sc*dTdP*dPdr*drdy;
		A(5, 2) = 3.0 * mu_cb / pow(r_cbtosc, 5.0)*rvec_cbtosc[2] * rvec_cbtosc[2] - mu_cb / pow(r_cbtosc, 3.0) + 3.0 * mu_3B*rvec_3Btosc[2] * rvec_3Btosc[2] / pow(r_3Btosc, 5.0) - mu_3B / pow(r_3Btosc, 3.0) + u[2] * D / m_sc*dTdP*dPdr*drdz;

		//A23 dadm
		A(3, 6) = -u[0] * D * T / (m_sc*m_sc);
		A(4, 6) = -u[1] * D * T / (m_sc*m_sc);
		A(5, 6) = -u[2] * D * T / (m_sc*m_sc);
		//
		//A24 dadu
		A(3, 7) = D * T / m_sc;
		A(3, 8) = 0.0;
		A(3, 9) = 0.0;
		A(4, 7) = 0.0;
		A(4, 8) = D * T / m_sc;
		A(4, 9) = 0.0;
		A(5, 7) = 0.0;
		A(5, 8) = 0.0;
		A(5, 9) = D * T / m_sc;

		//A31 dmdotdr
		A(6, 0) = -unorm * D * dmdotmaxdP * dPdr * drdx;
		A(6, 1) = -unorm * D * dmdotmaxdP * dPdr * drdy;
		A(6, 2) = -unorm * D * dmdotmaxdP * dPdr * drdz;

		//A34 dmdotdu
		A(6, 7) = -dunormdux * D * mdotmax;
		A(6, 8) = -dunormduy * D * mdotmax;
		A(6, 9) = -dunormduz * D * mdotmax;


		//DEBUG
		if (DEBUG)
		{
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
					std::cout << A(i, j) << " ";

				std::cout << std::endl;
			}
			getchar();
		}

		//*******************
		//
		//Phi dot calculation
		//
		//*******************

		//Form the Phi matrix
		//The Phi matrix is comprised of decision variables
		int xcount = 10;

		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < columns; ++j)
			{
				Phi(i, j) = x[xcount];
				++xcount;
			}
		}

		//DEBUG
		if (DEBUG)
		{
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
					std::cout << Phi(i, j) << " ";

				std::cout << std::endl;
			}
			getchar();
		}

		//Magical STM creation equation
		Phidot = A*Phi;
		int fcount = 10;

		//DEBUG
		if (DEBUG)
		{
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
					std::cout << Phidot(i, j) << " ";

				std::cout << std::endl;
			}
			getchar();
		}

		//Package the Phidot entries into the gradient vector f
		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < columns; ++j)
			{
				f[fcount] = Phidot(i, j);
				++fcount;
			}
		}


	}//end FBLT_EOM

}}} //close off the namespaces