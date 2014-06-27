//header file for Lambert solver by Arora and Russell
//A FAST AND ROBUST MULTIPLE REVOLUTION LAMBERT ALGORITHM USING A COSINE TRANSFORMATION
//AAS Hilton Head 2013
//implemented by Jacob Englander and Matthew A. Vavrina in C++ on 6/11/2014

#ifndef LAMBERT_ARORA_RUSSELL
#define LAMBERT_ARORA_RUSSELL

#include <cmath>
#include <vector>
#include <math.h> 
#include <complex>

namespace EMTG { namespace Astrodynamics {
	void Lambert_AroraRussell(const double* R1,
							  const double* R2, 
							  const double& TOFin,
							  const double& mu,
							  const int& Nrev,
							  const bool& LongWay,
							  const bool& ShortPeriod,
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

	double compute_TOF(const double &k, 
					   const double &S, 
					   const double &tau,
					   const double &W);

	double compute_kb(const double &k_bGuess, 
					const double &tau, 
					const double &S,
					const double &Nrev, 
					const double &tolerance, 
					const double &max_iterations, 
					const double &sq2, 
					const double &eps);

}
} //close namespace
#endif