//Amsterdam root finder


#ifndef _AMSTERDAM
#define _AMSTERDAM

#define AMSTERDAM_NUMERATOR(dab,dcb,fa,fb,fc) fb*(dab*fc*(fc-fb)-fa*dcb*(fa-fb))
#define AMSTERDAM_DENOMINATOR(fa,fb,fc) (fc-fb)*(fa-fb)*(fa-fc)

namespace EMTG { namespace Astrodynamics {

	double f(double x, const double p1, const double p2);

	double Amsterdam_Method(double a, double c, const double p1, const double p2);


}};

#endif //_AMSTERDAM