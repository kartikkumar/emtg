#include "GTOC7_solution_check.h"

std::vector <double> cross(std::vector <double> & A, std::vector <double> & B)
{  

	std::vector <double> answer (3, 0.0);

    answer[0]=A[1]*B[2]-A[2]*B[1];
    answer[1]=A[2]*B[0]-A[0]*B[2];
    answer[2]=A[0]*B[1]-A[1]*B[0];

	return answer;
}