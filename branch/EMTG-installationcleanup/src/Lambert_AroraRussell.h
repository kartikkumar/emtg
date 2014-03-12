//header file for Lambert solver by Arora and Russell
//A FAST AND ROBUST MULTIPLE REVOLUTION LAMBERT ALGORITHM USING A COSINE TRANSFORMATION
//AAS Hilton Head 2013
//implemented by Jacob Englander in C++ for EMTG, 3-5-2014

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