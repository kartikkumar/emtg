// Dormand-Prince (DOPRI) 8th(7th) order 13 step algorithm
// DOPRI constants are from Numerical Recipies
// Jacob Englander 9/10/2012

#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cstring>

#include "rk8713M.h"
#include "equations_of_motion.h"

using namespace std;

namespace EMTG {
	namespace integration
	{
		rk8713M::rk8713M(){} //the default constructor will never be called

		rk8713M::rk8713M(int ns_in)
		{
			ns = ns_in;
			
			x_left = new double[ns];
			f1.resize(ns);
			f2.resize(ns);
			f3.resize(ns);
			f4.resize(ns);
			f5.resize(ns);
			f6.resize(ns);
			f7.resize(ns);
			f8.resize(ns);
			f9.resize(ns);
			f10.resize(ns);
			f11.resize(ns);
			f12.resize(ns);
			f13.resize(ns);
			df1dTOF.resize(ns);
			df2dTOF.resize(ns);
			df3dTOF.resize(ns);
			df4dTOF.resize(ns);
			df5dTOF.resize(ns);
			df6dTOF.resize(ns);
			df7dTOF.resize(ns);
			df8dTOF.resize(ns);
			df9dTOF.resize(ns);
			df10dTOF.resize(ns);
			df11dTOF.resize(ns);
			df12dTOF.resize(ns);
			df13dTOF.resize(ns);
			y.resize(ns);
			dydTOF.resize(ns);
			dX_leftdTOF.resize(ns);
			dX_rightdTOF.resize(ns);


			//x_left = new double[ns];
			//f1 = new double[ns];
			//f2 = new double[ns];
			//f3 = new double[ns];
			//f4 = new double[ns];
			//f5 = new double[ns];
			//f6 = new double[ns];
			//f7 = new double[ns];
			//f8 = new double[ns];
			//f9 = new double[ns];
			//f10 = new double[ns];
			//f11 = new double[ns];
			//f12 = new double[ns];
			//f13 = new double[ns];
			//y = new double[ns];
		}

		//destructor
		rk8713M::~rk8713M()
		{
			//clean everything up
			//delete[] y;
			//delete[] f13;
			//delete[] f12;
			//delete[] f11;
			//delete[] f10;
			//delete[] f9;
			//delete[] f8;
			//delete[] f7;
			//delete[] f6;
			//delete[] f5;
			//delete[] f4;
			//delete[] f3;
			//delete[] f2;
			//delete[] f1;
			delete[] x_left;
		}

		void rk8713M::rk8713M_step(double * x_left, // spacecraft state at the LHS of the current RK sub-step 
			double * x_right, // pointer to store the spacecraft's state at the RHS of the current sub-step
			double * u, // control 3 vector
			std::vector <double> & f1, // gradient vector from the first RK stage
			std::vector <double> & df1dTOF, // TOF derivative of gradient vector from the first RK stage
			const double & t_left_step, // epoch at the LHS of the current RK sub-step
			const double & dt_left_stepdTOF,
			const double & t_0, // launch epoch
			double & h, // RK sub-step size
			double & dhdTOF, // TOF derivative of RK sub-step size
			double * error, // pointer to store error between 7th order and 8th order solutions
			
			void(*EOM)(double * x, // spacecraft's current state at left hand side of the current RK sub-step
			double * dx_dTOF,
			const double & t_left_step, // current epoch in TU's
			const double & dt_left_stepdTOF,
			const double & c2,
			const double & h,
			const double & dhdTOF,
			const double & t0, // launch epoch in seconds
			double * u, // throttle parameter vector
			double * f, // EOM gradient vector
			double * dfdTOF,
			double * thrust, // pointer that will extract info from the engine model (for storage in an archive)
			double * mdot, // pointer that will extract info from the engine model (for storage in an archive)
			double * Isp, // pointer that will extract info from the engine model (for storage in an archive)
			double * power, // pointer that will extract info from the engine model (for storage in an archive)
			double * active_power, // pointer that will extract info from the engine model (for storage in an archive)
			int * number_of_active_engines, // pointer that will extract info from the engine model (for storage in an archive)
			int & STMrows,
			int & STMcolumns,
			void * optionsvoidpointer, // passes the options structure through
			void * Universepointer, // passes the universe structure through
			void * ControllerPointer), // I assume Jacob was experimenting with feedback controllers?

			double * thrust, // pointer to extract thrust from engine model
			double * mdot, // pointer to extract mass flow rate from engine model
			double * Isp, // pointer to extract Isp from engine model
			double * power, // pointer to extract total power available to the spacecraft
			double * active_power, // total power used by ALL engines (if there are multiple)
			int * number_of_active_engines, // pointer to extract the number of engines that can fire
			int & STMrows,
			int & STMcolumns,
			void * optionspointer, // pass in options
			void * Universepointer, // pass in universe
			void * ControllerPointer // pass in controller
			)
		{

			*error = 0.0;//sufficiently high

			//RK node constants
			//these encode the positions within the RK step of each stage
			const static double c2 = 1.0 / 18.0;
			const static double c3 = 1.0 / 12.0;
			const static double c4 = 1.0 / 8.0;
			const static double c5 = 5.0 / 16.0;
			const static double c6 = 3.0 / 8.0;
			const static double c7 = 59.0 / 400.0;
			const static double c8 = 93.0 / 200.0;
			const static double c9 = 5490023248.0 / 9719169821.0;
			const static double c10 = 13.0 / 20.0;
			const static double c11 = 1201146811.0 / 1299019798.0;
			const static double c12 = 1.0;
			const static double c13 = 1.0;

			const static double a21 = 1.0 / 18.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + (h*a21)*f1[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a21*f1[i] + h*a21*df1dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c2,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f2,
				   df2dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer); //This is the call to the EOM function, replace this call with the call to your EOM function of choice.

			const static double a31 = 1.0 / 48.0;
			const static double a32 = 1.0 / 16.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a31*f1[i] + h*a32*f2[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a31*f1[i] + h*a31*df1dTOF[i] + dhdTOF*a32*f2[i] + h*a32*df2dTOF[i];
			}


			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c3,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f3,
				   df3dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a41 = 1.0 / 32.0;
			//const static double a42 = 0;
			const static double a43 = 3.0 / 32.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + (h*a41)*f1[i] + h*a43*f3[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a41*f1[i] + h*a41*df1dTOF[i] + dhdTOF*a43*f3[i] + h*a43*df3dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c4,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f4,
				   df4dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a51 = 5.0 / 16.0;
			//const static double a52 = 0;
			const static double a53 = -75.0 / 64.0;
			const static double a54 = 75.0 / 64.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a51*f1[i] + h*a53*f3[i] + h*a54*f4[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a51*f1[i] + h*a51*df1dTOF[i] + dhdTOF*a53*f3[i] + h*a53*df3dTOF[i] + dhdTOF*a54*f4[i] + h*a54*df4dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c5,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f5,
				   df5dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a61 = 3.0 / 80.0;
			//const static double a62 = 0;
			//const static double a63 = 0;
			const static double a64 = 3.0 / 16.0;
			const static double a65 = 3.0 / 20.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a61*f1[i] + h*a64*f4[i] + h*a65*f5[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a61*f1[i] + h*a61*df1dTOF[i] + dhdTOF*a64*f4[i] + h*a64*df4dTOF[i] + dhdTOF*a65*f5[i] + h*a65*df5dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c6,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f6,
				   df6dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a71 = 29443841.0 / 614563906.0;
			//const static double a72 = 0;
			//const static double a73 = 0;
			const static double a74 = 77736538.0 / 692538347.0;
			const static double a75 = -28693883.0 / 1125000000.0;
			const static double a76 = 23124283.0 / 1800000000.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a71*f1[i] + h*a74*f4[i] + h*a75*f5[i] + h*a76*f6[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a71*f1[i] + h*a71*df1dTOF[i] + dhdTOF*a74*f4[i] + h*a74*df4dTOF[i] + dhdTOF*a75*f5[i] + h*a75*df5dTOF[i] + dhdTOF*a76*f6[i] + h*a76*df6dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c7,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f7,
				   df7dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a81 = 16016141.0 / 946692911.0;
			//const static double a82 = 0;
			//const static double a83 = 0;
			const static double a84 = 61564180.0 / 158732637.0;
			const static double a85 = 22789713.0 / 633445777.0;
			const static double a86 = 545815736.0 / 2771057229.0;
			const static double a87 = -180193667.0 / 1043307555.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a81*f1[i] + h*a84*f4[i] + h*a85*f5[i] + h*a86*f6[i] + h*a87*f7[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a81*f1[i] + h*a81*df1dTOF[i] + dhdTOF*a84*f4[i] + h*a84*df4dTOF[i] + dhdTOF*a85*f5[i] + h*a85*df5dTOF[i] + dhdTOF*a86*f6[i] + h*a86*df6dTOF[i] + dhdTOF*a87*f7[i] + h*a87*df7dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c8,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f8,
				   df8dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a91 = 39632708.0 / 573591083.0;
			//const static double a92 = 0;
			//const static double a93 = 0;
			const static double a94 = -433636366.0 / 683701615.0;
			const static double a95 = -421739975.0 / 2616292301.0;
			const static double a96 = 100302831.0 / 723423059.0;
			const static double a97 = 790204164.0 / 839813087.0;
			const static double a98 = 800635310.0 / 3783071287.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a91*f1[i] + h*a94*f4[i] + h*a95*f5[i] + h*a96*f6[i] + h*a97*f7[i] + h*a98*f8[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a91*f1[i] + h*a91*df1dTOF[i] + dhdTOF*a94*f4[i] + h*a94*df4dTOF[i] + dhdTOF*a95*f5[i] + h*a95*df5dTOF[i] + dhdTOF*a96*f6[i] + h*a96*df6dTOF[i] + dhdTOF*a97*f7[i] + h*a97*df7dTOF[i] + dhdTOF*a98*f8[i] + h*a98*df8dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c9,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f9,
				   df9dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a10_1 = 246121993.0 / 1340847787.0;
			//const static double a10_2 = 0;
			//const static double a10_3 = 0;
			const static double a10_4 = -37695042795.0 / 15268766246.0;
			const static double a10_5 = -309121744.0 / 1061227803.0;
			const static double a10_6 = -12992083.0 / 490766935.0;
			const static double a10_7 = 6005943493.0 / 2108947869.0;
			const static double a10_8 = 393006217.0 / 1396673457.0;
			const static double a10_9 = 123872331.0 / 1001029789.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a10_1*f1[i] + h*a10_4*f4[i] + h*a10_5*f5[i] + h*a10_6*f6[i] + h*a10_7*f7[i] + h*a10_8*f8[i] + h*a10_9*f9[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a10_1*f1[i] + h*a10_1*df1dTOF[i] + dhdTOF*a10_4*f4[i] + h*a10_4*df4dTOF[i] + dhdTOF*a10_5*f5[i] + h*a10_5*df5dTOF[i] + dhdTOF*a10_6*f6[i] + h*a10_6*df6dTOF[i] + dhdTOF*a10_7*f7[i] + h*a10_7*df7dTOF[i] + dhdTOF*a10_8*f8[i] + h*a10_8*df8dTOF[i] + dhdTOF*a10_9*f9[i] + h*a10_9*df9dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c10,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f10,
				   df10dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a11_1 = -1028468189.0 / 846180014.0;
			//const static double a11_2 = 0;
			//const static double a11_3 = 0;
			const static double a11_4 = 8478235783.0 / 508512852.0;
			const static double a11_5 = 1311729495.0 / 1432422823.0;
			const static double a11_6 = -10304129995.0 / 1701304382.0;
			const static double a11_7 = -48777925059.0 / 3047939560.0;
			const static double a11_8 = 15336726248.0 / 1032824649.0;
			const static double a11_9 = -45442868181.0 / 3398467696.0;
			const static double a11_10 = 3065993473.0 / 597172653.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a11_1*f1[i] + h*a11_4*f4[i] + h*a11_5*f5[i] + h*a11_6*f6[i] + h*a11_7*f7[i] + h*a11_8*f8[i] + h*a11_9*f9[i] + h*a11_10*f10[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a11_1*f1[i] + h*a11_1*df1dTOF[i] + dhdTOF*a11_4*f4[i] + h*a11_4*df4dTOF[i] + dhdTOF*a11_5*f5[i] + h*a11_5*df5dTOF[i] + dhdTOF*a11_6*f6[i] + h*a11_6*df6dTOF[i] + dhdTOF*a11_7*f7[i] + h*a11_7*df7dTOF[i] + dhdTOF*a11_8*f8[i] + h*a11_8*df8dTOF[i] + dhdTOF*a11_9*f9[i] + h*a11_9*df9dTOF[i] + dhdTOF*a11_10*f10[i] + h*a11_10*df10dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c11,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f11,
				   df11dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a12_1 = 185892177.0 / 718116043.0;
			//const static double a12_2 = 0;
			//const static double a12_3 = 0;
			const static double a12_4 = -3185094517.0 / 667107341.0;
			const static double a12_5 = -477755414.0 / 1098053517.0;
			const static double a12_6 = -703635378.0 / 230739211.0;
			const static double a12_7 = 5731566787.0 / 1027545527.0;
			const static double a12_8 = 5232866602.0 / 850066563.0;
			const static double a12_9 = -4093664535.0 / 808688257.0;
			const static double a12_10 = 3962137247.0 / 1805957418.0;
			const static double a12_11 = 65686358.0 / 487910083.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a12_1*f1[i] + h*a12_4*f4[i] + h*a12_5*f5[i] + h*a12_6*f6[i] + h*a12_7*f7[i] + h*a12_8*f8[i] + h*a12_9*f9[i] + h*a12_10*f10[i] + h*a12_11*f11[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a12_1*f1[i] + h*a12_1*df1dTOF[i] + dhdTOF*a12_4*f4[i] + h*a12_4*df4dTOF[i] + dhdTOF*a12_5*f5[i] + h*a12_5*df5dTOF[i] + dhdTOF*a12_6*f6[i] + h*a12_6*df6dTOF[i] + dhdTOF*a12_7*f7[i] + h*a12_7*df7dTOF[i] + dhdTOF*a12_8*f8[i] + h*a12_8*df8dTOF[i] + dhdTOF*a12_9*f9[i] + h*a12_9*df9dTOF[i] + dhdTOF*a12_10*f10[i] + h*a12_10*df10dTOF[i] + dhdTOF*a12_11*f11[i] + h*a12_11*df11dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c12,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f12,
				   df12dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);

			const static double a13_1 = 403863854.0 / 491063109.0;
			//const static double a13_2 = 0;
			//const static double a13_3 = 0;
			const static double a13_4 = -5068492393.0 / 434740067.0;
			const static double a13_5 = -411421997.0 / 543043805.0;
			const static double a13_6 = 652783627.0 / 914296604.0;
			const static double a13_7 = 11173962825.0 / 925320556.0;
			const static double a13_8 = -13158990841.0 / 6184727034.0;
			const static double a13_9 = 3936647629.0 / 1978049680.0;
			const static double a13_10 = -160528059.0 / 685178525.0;
			const static double a13_11 = 248638103.0 / 1413531060.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*a13_1*f1[i] + h*a13_4*f4[i] + h*a13_5*f5[i] + h*a13_6*f6[i] + h*a13_7*f7[i] + h*a13_8*f8[i] + h*a13_9*f9[i] + h*a13_10*f10[i] + h*a13_11*f11[i];
				dydTOF[i] = dx_leftdTOF[i] + dhdTOF*a13_1*f1[i] + h*a13_1*df1dTOF[i] + dhdTOF*a13_4*f4[i] + h*a13_4*df4dTOF[i] + dhdTOF*a13_5*f5[i] + h*a13_5*df5dTOF[i] + dhdTOF*a13_6*f6[i] + h*a13_6*df6dTOF[i] + dhdTOF*a13_7*f7[i] + h*a13_7*df7dTOF[i] + dhdTOF*a13_8*f8[i] + h*a13_8*df8dTOF[i] + dhdTOF*a13_9*f9[i] + h*a13_9*df9dTOF[i] + dhdTOF*a13_10*f10[i] + h*a13_10*df10dTOF[i] + dhdTOF*a13_11*f11[i] + h*a13_11*df11dTOF[i];
			}

			(*EOM)(y,
				   dydTOF,
				   t_left_step,
				   dt_left_stepdTOF,
				   c13,
				   h,
				   dhdTOF,
				   t_0,
				   u,
				   f13,
				   df13dTOF,
				   thrust,
				   mdot,
				   Isp,
				   power,
				   active_power,
				   number_of_active_engines,
				   STMrows,
				   STMcolumns,
				   optionspointer,
				   Universepointer,
				   ControllerPointer);


			//8th order solution -- do this first because we need x_left unchanged; when I did 5th order first I was then accumulating on the value of x_left twice for y
			const static double b1upper = 14005451.0 / 335480064.0;
			//const static double b2upper = 0;
			//const static double b3upper = 0;
			//const static double b4upper = 0;
			//const static double b5upper = 0;
			const static double b6upper = -59238493.0 / 1068277825.0;
			const static double b7upper = 181606767.0 / 758867731.0;
			const static double b8upper = 561292985.0 / 797845732.0;
			const static double b9upper = -1041891430.0 / 1371343529.0;
			const static double b10upper = 760417239.0 / 1151165299.0;
			const static double b11upper = 118820643.0 / 751138087.0;
			const static double b12upper = -528747749.0 / 2220607170.0;
			const static double b13upper = 1.0 / 4.0;

			for (int i = 0; i < ns; ++i)
			{
				y[i] = x_left[i] + h*b1upper*f1[i] + h*b6upper*f6[i] + h*b7upper*f7[i] + h*b8upper*f8[i] + h*b9upper*f9[i] + h*b10upper*f10[i] + h*b11upper*f11[i] + h*b12upper*f12[i] + h*b13upper*f13[i];
			}

			//7th order solution (compared against the 8th order solution to determine the error due to one RK step)
			const static double b1lower = 13451932.0 / 455176623.0;
			//const static double b2lower = 0; 
			//const static double b3lower = 0; 
			//const static double b4lower = 0; 
			//const static double b5lower = 0; 
			const static double b6lower = -808719846.0 / 976000145.0;
			const static double b7lower = 1757004468.0 / 5645159321.0;
			const static double b8lower = 656045339.0 / 265891186.0;
			const static double b9lower = -3867574721.0 / 1518517206.0;
			const static double b10lower = 465885868.0 / 322736535.0;
			const static double b11lower = 53011238.0 / 667516719.0;
			const static double b12lower = 2.0 / 45.0;
			//const static double b13lower = 0;


			for (int i = 0; i < ns; ++i)
			{
				//Determine right hand values of states to 7th order
				x_right[i] = x_left[i] + h*b1lower*f1[i] + h*b6lower*f6[i] + h*b7lower*f7[i] + h*b8lower*f8[i] + h*b9lower*f9[i] + h*b10lower*f10[i] + h*b11lower*f11[i] + h*b12lower*f12[i];

				//these are the TOF derivatives of the right-hand state, they become the derivatives of the left hand state for the next step (or are the final TOF derivatives)
				dx_rightdTOF[i] = dx_leftdTOF[i] + dhdTOF*b1lower*f1[i] + h*b1lower*df1dTOF[i] + dhdTOF*b6lower*f6[i] + h*b6lower*df6dTOF[i] + dhdTOF*b7lower*f7[i] + h*b7lower*df7dTOF[i] + dhdTOF*b8lower*f8[i] + h*b8lower*df8dTOF[i] + dhdTOF*b9lower*f9[i] + h*b9lower*df9dTOF[i] + dhdTOF*b10lower*f10[i] + h*b10lower*df10dTOF[i] + dhdTOF*b11lower*f11[i] + h*b11lower*df11dTOF[i] + dhdTOF*b12lower*f12[i] + h*b12lower*df12dTOF[i];


				//take the 8th order solution as truth and compare it with the 7th order solution to quantify the error for this substep
				*error = max(*error, fabs(x_right[i] - y[i]));
			}
		}

		void rk8713M::adaptive_step_int(double * x_left_in, // spacecraft's state at the left boundary of the current segment (FBLT "step")
			double * x_right, // pointer to the spacecraft state at the RHS of the segment (for state data archive)
			double * uleft, // 3 vector encoding the three throttle parameters for this FBLT segment
			const double & t_left, // current epoch in TU's
			const double & t_0, // launch epoch (NOT time at beginning of phase, unless this is a 1 - phase mission)
			double const & local_step, // rough guess at how big the integration step size should be (FBLT segment time / 2) in TU's
			double * resumeH, // NO LONGER USED...in the original implementation we forcibly discretized a priori to the integration
			double * resumeError, // NO LONGER USED...same reason, the outer for loop around the step do-while is no longer present
			double const & PRECISION_TARGET, // integration error tolerance, currently 1.0e-8 in EMTG, this is set in FBLTphase.cpp (hard-coded)
			
			void(*EOM)(double * x, // spacecraft's current state at left hand side of the current RK sub-step
			double * dx_dTOF,
			const double & t_left_step, // current epoch in TU's
			const double & dt_left_stepdTOF, 
			const double & c2,
			const double & h,
			const double & dhdTOF,
			const double & t0, // launch epoch in seconds
			double * u, // throttle parameter vector
			double * f, // EOM gradient vector
			double * dfdTOF,
			double * thrust, // pointer that will extract info from the engine model (for storage in an archive)
			double * mdot, // pointer that will extract info from the engine model (for storage in an archive)
			double * Isp, // pointer that will extract info from the engine model (for storage in an archive)
			double * power, // pointer that will extract info from the engine model (for storage in an archive)
			double * active_power, // pointer that will extract info from the engine model (for storage in an archive)
			int * number_of_active_engines, // pointer that will extract info from the engine model (for storage in an archive)
			int & STMrows,
			int & STMcolumns,
			void * optionsvoidpointer, // passes the options structure through
			void * Universepointer, // passes the universe structure through
			void * ControllerPointer), // I assume Jacob was experimenting with feedback controllers?

			double* thrust,
			double* mdot,
			double* Isp,
			double* power,
			double* active_power,
			int* number_of_active_engines,
			int &STMrows,
			int &STMcolumns,
			void* optionspointer, void* Universepointer, void* ControllerPointer)
		{

			double accumulatedH = 0.0;
			double effectiveH = *resumeH;
			double precision_error = *resumeError;
			bool last_substep;

			if (*resumeH > local_step)
				effectiveH = local_step;

			for (int k = 0; k < ns; ++k)
				x_left[k] = x_left_in[k];

			last_substep = false; //at the beginning of a segment we are on the FIRST substep

			//Forward Integration of the states
			do
			{ //loop until we get all the way through a full h (a full segment)

				//-> INSERT EOM 1st call and STORE f1 values. Replace this with your EOM call of choice.	
				(*EOM)(x_left, 
					   t_left, 
					   t_0, 
					   uleft, 
					   f1, 
					   thrust, 
					   mdot, 
					   Isp, 
					   power, 
					   active_power, 
					   number_of_active_engines, 
					   STMrows, 
					   STMcolumns, 
					   optionspointer, 
					   Universepointer, 
					   ControllerPointer);


				//take a trial substep
				do
				{ //cycle until trial substep is done with sufficient accuracy
					//copy in for the current step our left point

					if (!last_substep)
					{
						//no error!  give it a real value so we don't divide by zero.
						if (precision_error == 0.0)
						{
							precision_error = 1e-15; //Almost zero!
						}

						if (precision_error >= PRECISION_TARGET)
							effectiveH = 0.98*effectiveH*pow(PRECISION_TARGET / precision_error, 0.17);

						else
						{
							effectiveH = 1.01*effectiveH*pow(PRECISION_TARGET / precision_error, 0.18);

							//if our increased step kicks us too long, make it shorter anyway and just run it
							if (fabs(local_step - accumulatedH) < fabs(effectiveH) && !last_substep)
							{
								effectiveH = local_step - accumulatedH;
							}

						}

						if (fabs(effectiveH) < 1e-14)
						{//H is too small
							cout << "rk8713M: H Got too Small. The integrator has Alexed. Aborting." << endl;
							throw 13;
						}

					}

					else if (precision_error > PRECISION_TARGET && last_substep)
					{
						//we got here because we thought it was the last substep, and upon calculation it 
						//was too big and not precise enough so we have to shrink it
						last_substep = false; //not last step after all
						effectiveH = 0.98*effectiveH*pow(PRECISION_TARGET / precision_error, 0.17);
					}


					//rk takes in x1 as the left-hand side, but returns it as the answer (right hand side)
					//rk8713M (x_left, x_right, uleft, f1, effectiveH, precision_error, segment, ns, nc, nseg, sat_num, *design);
					rk8713M_step(x_left, x_right, uleft, f1, t_left + accumulatedH, t_0, effectiveH, &precision_error, EOM, thrust, mdot, Isp, power, active_power, number_of_active_engines, STMrows, STMcolumns, optionspointer, Universepointer, ControllerPointer);

				} while (precision_error > PRECISION_TARGET);

				//trial substep was accurate enough; it becomes the new left
				for (int statenum = 0; statenum < ns; ++statenum)
				{
					x_left[statenum] = x_right[statenum];
				}


				accumulatedH += effectiveH;
				//if our next step will push us over, reduce it down to be as small as necessary to hit target exactly
				if (fabs(local_step - accumulatedH) < fabs(effectiveH) && fabs(local_step - accumulatedH) > 0 && !last_substep)
				{
					*resumeH = effectiveH;
					*resumeError = precision_error;
					effectiveH = local_step - accumulatedH;
					last_substep = true; //assume that the next substep will be the last substep now
				}


			} while (fabs(accumulatedH) < fabs(local_step));

			//end of integration routine
		}


	}
} //namespace EMTG::integration