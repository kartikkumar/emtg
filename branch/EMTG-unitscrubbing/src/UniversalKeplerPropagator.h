//EMTG state transition matrix class
//Jacob Englander 3/25/2013

#ifndef EMTG_STATETRANSIT
#define EMTG_STATETRANSIT

#include "EMTG_Matrix.h"

namespace EMTG { namespace Astrodynamics {

class UniversalKeplerPropagator
{
public:
	//constructor
	UniversalKeplerPropagator();

	//destructor
	virtual ~UniversalKeplerPropagator();

	//methods
	int propagate_and_compute_STM(const double* state0, double* statef, const double& deltat_seconds, const double& mu_seconds, const double& LU, const double& TU, bool computeSTM);

	//fields
	
	//state vectors
	math::Matrix<double> rv0;
	math::Matrix<double> vv0;
	math::Matrix<double> rvf;
	math::Matrix<double> vvf;

	//specific angular momentum
	math::Matrix<double> hv;

	//STM-related matrices
	math::Matrix<double> M;
	math::Matrix<double> STM;
};

}}

#endif //EMTG_STATETRANSIT