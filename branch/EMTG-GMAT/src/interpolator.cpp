//EMTG interpolation class header file
//Jacob Englander 6/14/2013

#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>

#include "interpolator.h"
#include "EMTG_math.h"



using namespace std;


namespace EMTG { namespace math {

	//default constructor
	interpolator::interpolator(void)
	{
	}

	//constructor for a known data table
	interpolator::interpolator( vector < pair<double, double> > InputTable)
	{
		DataTable = InputTable;
	}

	//destructor
	interpolator::~interpolator(void)
	{
	}

	//function to interpolate the data table
	double interpolator::interpolate(double x)
	{
		// Check if x is out of bound
		if ( x > DataTable.back().first )
			return -LARGE;
		if ( x < DataTable.front().first )
			return LARGE;

		vector<pair<double, double> >::iterator it, it2;

		it = lower_bound(DataTable.begin(), DataTable.end(), make_pair(x, -LARGE));

		// Corner case

		if (it == DataTable.begin()) 
			return it->second;

		it2 = it;
		--it2;

		return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
	}

}} //close namespace EMTG::math