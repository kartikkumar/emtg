//differential equations of motion for EMTG v8

#include "missionoptions.h"

#ifndef _EOM
#define _EOM

namespace EMTG {namespace Astrodynamics { namespace EOM
{
	//equations of motion for an object moving in the heliocentric inertial frame with a thrust term
	void EOM_inertial_continuous_thrust(std::vector <double> & x,
		                                EMTG::math::Matrix <double> & dxdTOF,
										const double & epoch_step_left,
										std::vector <double> & depoch_step_leftdTOF,
										const double & c,
										const double & h,
										const double & dhdTOF,
										const double & t0,
										const std::vector <double> & u,
										std::vector <double> & f,
										EMTG::math::Matrix <double> & dfdTOF,
										const int & phase_num,
										double * thrust,
										double * mdot,
										double * Isp,
										double * power,
										double * active_power,
                                        int* number_of_active_engines,
                                        int & STMrows,
                                        int & STMcolumns,
										void * optionsvoidpointer, 
										void * Universepointer, 
										void * ControllerPointer);


	//equations of motion for an object moving in the heliocentric inertial frame with a thrust term, using the Sundman transformation for non-uniform spacing of control points
	void EOM_inertial_continuous_thrust_sundman(double* x,
												const double& t,
												const double& t0,
												double* u,
												double* f,
												double* thrust,
												double* mdot,
												double* Isp,
												double* power,
												double* active_power,
												int* number_of_active_engines,
												void* optionsvoidpointer, 
												void* Universepointer, 
												void* ControllerPointer);
}}}

#endif //_EOM