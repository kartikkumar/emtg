#include "GTOC7_solution_check.h"

double dot(std::vector <double> & vectorA, std::vector <double> & vectorB)
{
	double ScalarProduct=0;
	int i;
		
	for( i = 0; i < vectorA.size(); ++i )
		ScalarProduct=ScalarProduct+vectorA[i]*vectorB[i];
	
	return ScalarProduct;
}