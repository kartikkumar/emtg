//Amsterdam root finder

#include <math.h> // required for fabs()		
#include "AmsterdamRootFinder.h"
#include "EMTG_math.h"

namespace EMTG { namespace Astrodynamics
{

	double f(double x, const double p1, const double p2)
	{
		return p1*tan(x) - log(tan(0.5*x + 0.25*EMTG::math::PI)) - p2;
	}

	double Amsterdam_Method(double a, double c, const double p1, const double p2) 
	{
		int max_iterations = 500;
		double tolerance = 1e-15;

		double fa = f(a,p1,p2), b = 0.5 * ( a + c ), fc = f(c,p1,p2), fb = fa * fc;
		double delta, dab, dcb;
		int i;

		  // If the initial estimates do not bracket a root, set the err flag and //
		  // return.  If an initial estimate is a root, then return the root.     //

		if ( fb >= 0.0 )
			if ( fb > 0.0 )  
			{
				return 0.0; 
			}
			else return ( fa == 0.0 ) ? a : c;

		  // Insure that the initial estimate a < c. //

		if ( a > c ) 
		{
			delta = a; 
			a = c; 
			c = delta; 
			delta = fa; 
			fa = fc; 
			fc = delta;
		}

		  // If the function at the left endpoint is positive, and the function //
		  // at the right endpoint is negative.  Iterate reducing the length    //
		  // of the interval by either bisection or quadratic inverse           //
		  // interpolation.  Note that the function at the left endpoint        //
		  // remains nonnegative and the function at the right endpoint remains //
		  // nonpositive.                                                       //

		if ( fa > 0.0 )
			for ( i = 0; i < max_iterations; ++i)
			{

				// Are the two endpoints within the user specified tolerance ?

				if ( ( c - a ) < tolerance ) return 0.5 * ( a + c);

					// No, Continue iteration.

				fb = f(b,p1,p2);

				// Check that we are converging or that we have converged near //
				// the left endpoint of the inverval.  If it appears that the  //
				// interval is not decreasing fast enough, use bisection.      //
				if ( ( c - a ) < tolerance ) return 0.5 * ( a + c);
				if ( ( b - a ) < tolerance )
					if ( fb > 0 )
					{
						a = b;
						fa = fb;
						b = 0.5 * ( a + c );
						continue;
					}
					else return b;

				// Check that we are converging or that we have converged near //
				// the right endpoint of the inverval.  If it appears that the //
				// interval is not decreasing fast enough, use bisection.      //

				if ( ( c - b ) < tolerance )
					if ( fb < 0 ) 
					{
						c = b;
						fc = fb; 
						b = 0.5 * ( a + c );
						continue;
					}
					else return b;

				// If quadratic inverse interpolation is feasible, try it. //

				if (  ( fa > fb ) && ( fb > fc ) )
				{
					delta = AMSTERDAM_DENOMINATOR(fa,fb,fc);
					if ( delta != 0.0 )
					{
						dab = a - b;
						dcb = c - b;
						delta = AMSTERDAM_NUMERATOR(dab,dcb,fa,fb,fc) / delta;

						// Will the new estimate of the root be within the   //
						// interval?  If yes, use it and update interval.    //
						// If no, use the bisection method.                  //

						if ( delta > dab && delta < dcb ) 
						{
							if ( fb > 0.0 ) 
							{ 
								a = b; 
								fa = fb;
							}
							else if ( fb < 0.0 ) 
							{
								c = b;
								fc = fb;
							}
							else return b;

							b += delta;

							continue;
						}
					}
				}

				// If not, use the bisection method. //

				fb > 0 ? ( a = b, fa = fb ) : ( c = b, fc = fb );
				b = 0.5 * ( a + c );
			}
			else

			// If the function at the left endpoint is negative, and the function //
			// at the right endpoint is positive.  Iterate reducing the length    //
			// of the interval by either bisection or quadratic inverse           //
			// interpolation.  Note that the function at the left endpoint        //
			// remains nonpositive and the function at the right endpoint remains //
			// nonnegative.                                                       //

			for ( i = 0; i < max_iterations; ++i)
			{
				if ( ( c - a ) < tolerance ) return 0.5 * ( a + c);
				fb = f(b,p1,p2);

				if ( ( b - a ) < tolerance )
				if ( fb < 0 )
				{
					a = b;
					fa = fb;
					b = 0.5 * ( a + c );
					continue;
				}
				else return b;

				if ( ( c - b ) < tolerance )
					if ( fb > 0 )
					{
						c = b;
						fc = fb;
						b = 0.5 * ( a + c ); 
						continue;
					}
					else return b;

				if (  ( fa < fb ) && ( fb < fc ) ) 
				{
					delta = AMSTERDAM_DENOMINATOR(fa,fb,fc);
					if ( delta != 0.0 ) {
						dab = a - b;
						dcb = c - b;
						delta = AMSTERDAM_NUMERATOR(dab,dcb,fa,fb,fc) / delta;
						if ( delta > dab && delta < dcb )
						{
							if ( fb < 0.0 ) 
							{ 
								a = b;
								fa = fb;
							}
							else if ( fb > 0.0 ) 
							{ 
								c = b;
								fc = fb;
							}
							else return b;

							b += delta;

							continue;
						}
					}
				}
				fb < 0 ? ( a = b, fa = fb ) : ( c = b, fc = fb );
				b = 0.5 * ( a + c );
			}
			return  b;
	}
}}