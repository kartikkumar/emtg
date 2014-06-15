#include "GTOC7_solution_check.h"

double mag(std::vector <double> & vector)
{
	double magnitude;
	double temp = 0.0;

	for (size_t i=0; i < vector.size(); ++i)
	     temp += vector[i]*vector[i];

	if ( fabs(temp) >= 1.0e-16 )
	{
		magnitude=sqrt(temp);
	}
	else
	{
		magnitude=0;
	}
	return magnitude;
}