//header file for Lambert solver by Arora and Russell Placeholder Stub
//this is the expected interface of the future Lambert solver to be included in a upcoming EMTG release
//In the meantime, if you want a lambert solver, supply your own, and make sure its interface matches the stuff below.

#include "Lambert_AroraRussell.h"
#include "EMTG_math.h"

#include <cmath>

namespace EMTG { namespace Astrodynamics {
	//assume input and output units are consistent (i.e. designed for km and s but should work for anything)
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
							  int& iterations)
	{
		throw 42428642; //This throw is because the code is not present.  Once ready, it will be pulled from a branch and this file will be replaced.  If you are supplying your own, add Lambert code here.
	}
}} //close namespace