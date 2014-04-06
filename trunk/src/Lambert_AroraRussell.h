//header file for Lambert solver by Arora and Russell Placeholder Stub
//this is the expected interface of the future Lambert solver to be included in a upcoming EMTG release
//In the meantime, if you want a lambert solver, supply your own, and make sure its interface matches the stuff below.


#ifndef LAMBERT_ARORA_RUSSELL
#define LAMBERT_ARORA_RUSSELL

namespace EMTG { namespace Astrodynamics {
	void Lambert_AroraRussell(const double* R1,
							  const double* R2, 
							  const double& TOF,
							  const double& mu,
							  const int& Nrev,
							  const bool& LongWay,
							  const double& tolerance,
							  const int& max_iterations,
							  double* V1,
							  double* V2,
							  double& error,
							  int& iterations);

	double acoshAR(const double& b);

	double acosAR(const double& x);

	double compute_W(const double& k,
					 const double& m,
					 const int& Nrev);
}} //close namespace
#endif