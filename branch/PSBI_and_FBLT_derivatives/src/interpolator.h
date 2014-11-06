//EMTG interpolation class header file
//Jacob Englander 6/14/2013

#ifndef _EMTG_INTERPOLATOR
#define _EMTG_INTERPOLATOR

#include <vector>
#include <utility>

#include "EMTG_math.h"




using namespace std;


namespace EMTG { namespace math {

	class interpolator
	{
	public:
		//default constructor
		interpolator(void);

		//constructor for a known data table
		interpolator ( vector < pair<double, double> > InputTable);

		//destructor
		virtual ~interpolator(void);

		//fields
		vector< pair<double, double> > DataTable;

		//methods
		double interpolate(double x);
	};

}}//close namespace EMTG::math

#endif