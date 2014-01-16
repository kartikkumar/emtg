//header file for Lagrangian Kepler

#include "STM.h"
#include <math.h>


namespace Kepler
{
	//Lagrange-Laguerre-Conway Kepler propagator
	void KeplerLagrangeLaguerreConway(const double* state0, double* state, const double& mu, const double& propTime, STM& stm, const bool& compute_STM_flag);
	void KeplerLagrangeLaguerreConway(const double* state0, double* state, const double& mu, const double& propTime, double* F, double* G, double* Ft, double* Gt);


	//Laguerre-Conway solver for Kepler's equation
	double KeplerLaguerreConway(const double& ECC, const double& MN);

}//end namespace Kepler