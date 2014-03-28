//header file for Laguerre-Conway Kepler solver using the iterative n method of Ger
//from "Orbital Mechanics" by Prussing and Conway, Chapter 2
//and "Classical and Advanced Kepler Algorithms" by Gim J. Der
//STMs by Ellison

#include "STM.h"


#ifndef KEPLERLAGRANGELAGUERRECONWAYGER
#define KEPLERLAGRANGELAGUERRECONWAYGER

namespace Kepler
{
	void Kepler_Lagrange_Laguerre_Conway_Ger(const double* state0_kms,
											 double* state_kms,
											 const double& mu,
											 const double& LU,
											 const double& propTime,
											 double& F,
											 double& G,
											 double& Ft,
											 double& Gt,
											 double& Ftt,
											 double& Gtt,
											 STM& stm,
											 const bool& compute_STM_flag);
}

#endif