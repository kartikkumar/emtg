#include "GTOC7_solution_check.h"

/*
Author: Donald Ellison
Date: March 7th 2012

This function solves Kepler's equation using the Laguerre-Conway
numerical algorithm.

INPUTS: ecc: satellite eccentricity
        M: satellite mean anomaly

LOCAL VARIABLES: E: the calculated mean anomaly 
                 f: Keper's equation trancendental constraint
				 fprime: first derivative of Kepler's equation w.r.t. E
				 fdoubleprime: second derivative of Kepler's equation w.r.t E
				 tolerance: allowable tolerance for f constraint
				 n: integer variable appearing in Laguerre's algorithm

OUTPUTS: E: eccentric anomaly

FUNCTIONS CALLED: None
*/

double laguerreConway(double & ecc, double & M)
{
	double E;
	double f = 0.0;
	double fprime = 0.0;
	double fdoubleprime = 0.0;
	double tolerance = 1.0e-13;
	
	int n = 5;

	//Initial guess for E.
	E = (M*(1.0-sin(M+ecc))+(M+ecc)*sin(M))/(1.0+sin(M)-sin(M+ecc));

	//Set up constraint to zero out.
	f = M-(E-ecc*sin(E));
	fprime = -1.0 + ecc*cos(E);
	fdoubleprime = -ecc*sin(E);

	while(fabs(f)>tolerance)
	{
		while(fabs(f)>tolerance && fprime>0.0)
		{
			E = E-(n*f)/(fprime+(sqrt(fabs((n-1)*(n-1)*(fprime)*(fprime)-n*(n-1)*f*fdoubleprime))));

			f = M-(E-ecc*sin(E));
			fprime = -1+ecc*cos(E);
			fdoubleprime = -ecc*sin(E);
		}

		while(fabs(f)>tolerance && fprime<0.0)
		{
			E = E-(n*f)/(fprime-(sqrt(fabs((n-1)*(n-1)*(fprime)*(fprime)-n*(n-1)*f*fdoubleprime))));

			f = M-(E-ecc*sin(E));
			fprime = -1.0 + ecc*cos(E);
			fdoubleprime = -ecc*sin(E);
		}
	}

	//Correct the final value of E for over 2*pi rotations

	
	while(E>2.0*PI)
	{
		E = E-2.0*PI;
	}
	
	return E;
}