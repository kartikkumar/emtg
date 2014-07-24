//math functions for EMTG

#ifndef _EMTG_MATH
#define _EMTG_MATH

#include <cmath>
#include <cstdlib>

namespace EMTG { namespace math
{
	//constants
	const double PI = 3.14159265358979323;
	const double PIover2 = 3.14159265358979323 / 2.0;
	const double TwoPI = 3.14159265358979323 * 2.0;
	const double SMALL  = 1.0e-13;
	const double LARGE  = 1.0e+30;

	template <typename T> int sgn(T val) 
	{
		return (T(0) < val) - (val < T(0));
	}

	inline void unitv(const double *V_in, double *Ver_out)
	{
		double v_mod = 0;
		int i;

		for (i = 0;i < 3;i++)
		{
			v_mod += V_in[i]*V_in[i];
		}
	
		double sqrtv_mod = sqrt(v_mod);

		for (i = 0;i < 3;i++)
		{
			Ver_out[i] = (sqrtv_mod!=0) ? V_in[i]/sqrtv_mod : 0;
		}
	}

	inline void cross (const double A[], const double B[], double C[])
	{
		C[0]=A[1]*B[2]-A[2]*B[1];
		C[1]=A[2]*B[0]-A[0]*B[2];
		C[2]=A[0]*B[1]-A[1]*B[0];   
		return;
	}

	inline double acosh(double x) {return 2*log(sqrt((x+1.0)/2.0) + sqrt((x-1.0)/2.0));}

	inline double asinh(double x) {	return log (x+sqrt(1+x*x));}

	inline double norm(const double A[], size_t dim)
	{
		double norm2 = 0.0;
		for (size_t k=0; k<dim; ++k)
			norm2 += A[k]*A[k];

		return sqrt(norm2);
	}

	inline double dot(const double A[], const double B[], size_t dim)
	{
		double dot = 0.0;

		for (size_t k=0;k<dim;++k)
			dot += A[k]*B[k];

		return dot;
	}

	inline double maxabsvec(const double A[], size_t dim)
	{
		double maxval = 0.0;

		for (size_t k=0;k<dim;++k)
		{
			double absval = fabs(A[k]);

			if (absval > maxval)
				maxval = absval;
		}

		return maxval;
	}

	inline double mod(const double A, const double B)
	{
		//note B must be positive
		double principle_value = fmod(A,B);

		if (principle_value < 0)
			principle_value += B;

		return principle_value;
	}

	inline double angle_vector(const double A[], const double B[], size_t dim)
	{
		//computes the angle between two vectors
		double a = math::norm(A, dim);
		double b = math::norm(B, dim);

		return acos(math::dot(A, B, dim) / (a * b));
	}


	inline double EMTG_C_FUNC(double y)  //transcendental C
	{
			return (y > 0) ? (1 - cos(sqrt(y)))/y : (cosh(sqrt(-y)) - y)/(-y);
	}

	inline double EMTG_S_FUNC(double y)  //transcendental S
	{
		(y > 0) ? (sqrt(y) - sin(sqrt(y)))/sqrt(y*y*y) : (sinh(sqrt(-y)) - sqrt(-y)/sqrt(-y*y*y));
	}

	inline void matrix_vector_multiply(double* matrix, double* vector_in, double* vector_out, const int& dimension)
	{
		for (int i = 0; i < dimension; ++i)
		{
			vector_out[i] = 0.0;

			for (int j = 0; j < dimension; ++j)
			{
				vector_out[i] += matrix[i*dimension + j] * vector_in[j];
			}
		}
	}

}}
#endif //_EMTG_MATH
