//code to find mass delivered to orbit by the launch vehicle as a function of C3 and DLA
//all LV functions are polynomial fit to performance curves from KSC launch vehicle website
//all data is public domain except FH and SLS which are kept in a separate source file, NOT TO BE DISTRIBUTED

#include <math.h>
#include <iostream>
#include "missionoptions.h"

//if the proprietary switch is turned on, include the proprietary functions. This header file is NOT distributed with the open-source version of EMTG
#ifdef _EMTG_proprietary
#include "proprietary_functions.h"
#endif


namespace EMTG { namespace Astrodynamics {

void find_mass_to_orbit(double C3, double DLA, int LV_type, double* mass, double* dmdC3, missionoptions* options)
{
	//-2: custom launch vehicle
	//1: Falcon 9 (v1.0)			NLSII
	//2: Atlas V (501)				NLSII
	//3: Falcon 9 (v1.1)			NLSII
	//4: Atlas V (401)				NLSII
	//5: Atlas V (511)				NLSII
	//6: Atlas V (411)				NLSII
	//7: Atlas V (521)				NLSII
	//8: Atlas V (421)				NLSII
	//9: Atlas V (531)				NLSII
	//10: Atlas V (431)				NLSII
	//11: Atlas V (541)				NLSII
	//12: Atlas V (551)				NLSII
	//13: Atlas V (551) w/Star 48	NLSI
	//14: Falcon 9 Heavy			
	//15: Delta IV Heavy			NLSI
	//16: SLS Block 1

	if (LV_type > 0)
	{
		if (LV_type == 14)
		{
#ifdef _EMTG_proprietary
			EMTG::Proprietary::FalconHeavy_C3_curve(C3, mass, dmdC3);
#else
			std::cout << "Falcon Heavy model not included in open-source version" << std::endl;
			throw 1711;
#endif
			return;
		}

		if (LV_type == 16)
		{
#ifdef _EMTG_proprietary
			EMTG::Proprietary::SLS_C3_curve(C3, mass, dmdC3);
#else
			std::cout << "SLS model not included in open-source version" << std::endl;
			throw 1711;
#endif
			return;
		}

		double a1[] = {0,0,-8.00E-05,0,0,0,0,0,0,0,0,0,-3.5649E-08, 0, 8.0e-7, 0,-3.5649E-08, 0, 8.0e-7, 0};
		double a2[] = {2.81E-05,5.62E-06,-3.88E-05,-2.48E-05,-0.00009382,0.00006198,0.00012334,-0.00010901,0.00006465,0.00007576,-0.00138695,-0.00004196,2.31305E-05,0, -0.0002, 0};
		double a3[] = {-0.001517208,-0.00067462,0.00190765,0.00211196,0.00403555,-0.00295026,-0.00712978,0.00596291,-0.00276635,-0.00405132,0.03050117,-0.0035373,-0.006458248, 0, 0.0145, 0};
		double a4[] = {0.357956243,0.44268134,0.43698409,0.47075449,0.26854604,0.41398944,0.61102598,0.37650144,0.60963424,0.70733066,0.30110723,0.91375291,1.065903806, 0, 0.5825, 0};
		double a5[] = {-64.48375826,-78.8465652,-88.38856438,-98.4962944,-55.39915501,-69.35547443,-83.52984026,-91.90752777,-102.9890546,-111.2601399,-69.96585082,-110.0132867,-114.5157292, 0, -183.69, 0};
		double a6[] = {3034.683258,3930.871041,4655.158371,5235.260181,2094.89,3265.42,4195.86,4939.98,5595.09,6106.14,1974.88,3634.59,6595.738063, 0, 12170, 0};
		double C3min[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		double C3max[] = {35,50,40,60,60,60,60,60,60,60,60,60,200,0,100,0};
		int k = LV_type - 1;
		double C32 = C3*C3;
		double C33 = C3*C32;
		double C34 = C3*C33;
		double C35 = C3*C34;

		if (C3 > C3max[k] || C3 < C3min[k])
		{
			double expfun = exp(-C3/C3max[k]);
			*mass = expfun; //make sure that a nonzero, differentiable value is returned so that SNOPT can define the point
			*dmdC3 = -expfun / C3max[k];
		}
		else
		{
			*mass = a1[k]*C35 + a2[k]*C34 + a3[k]*C33 + a4[k]*C32 + a5[k]*C3 + a6[k] - options->LV_adapter_mass;
			*dmdC3 = 5*a1[k]*C34 + 4*a2[k]*C33 + 3*a3[k]*C32 + 2*a4[k]*C3 + a5[k];
		}
	}
	else if (LV_type == -2)
	{
		double C32 = C3*C3;
		double C33 = C3*C32;
		double C34 = C3*C33;
		double C35 = C3*C34;

		if (C3 > options->custom_LV_C3_bounds[1] || C3 < options->custom_LV_C3_bounds[0])
		{
			double expfun = exp(-C3/options->custom_LV_C3_bounds[1]);
			*mass = expfun - options->LV_adapter_mass; //make sure that a nonzero, differentiable value is returned so that SNOPT can define the point
			*dmdC3 = -expfun / options->custom_LV_C3_bounds[1];
		}
		else
		{
			double* a = options->custom_LV_coefficients;
			*mass = a[0]*C35 + a[1]*C34 + a[2]*C33 + a[3]*C32 + a[4]*C3 + a[5] - options->LV_adapter_mass;
			*dmdC3 = 5*a[0]*C34 + 4*a[1]*C33 + 3*a[2]*C32 + 2*a[3]*C3 + a[4];
		}
	}


	return;
}

}}