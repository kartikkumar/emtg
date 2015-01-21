//useful debug functions for EMTG


#ifndef _EMTG_DEBUG_FUNC
#define _EMTG_DEBUG_FUNC

#include <iostream>

using namespace std;


namespace EMTG { namespace debug_functions {

	void print_1D_array(double* printarray, int length)
	{
		for (int k = 0; k < length; ++k)
			cout << printarray[k] << " ";
		cout << endl;
	}
}}//close namespace

#endif