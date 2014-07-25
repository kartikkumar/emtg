//SNOPT user function

#include <ctime>
#include <iostream>

#include "problem.h"

#include "snopt.hh"




namespace EMTG { namespace Solvers {

	int SNOPT_user_function(integer    *Status, integer *n,    doublereal x[],
							integer    *needF,  integer *neF,  doublereal F[],
							integer    *needG,  integer *neG,  doublereal G[],
							char       *cu,     integer *lencu,
							integer    iu[],    integer *leniu,
							doublereal ru[],    integer *lenru )
	{
		//Step 1: create a pointer to the Problem object
		EMTG::problem* Problem = (EMTG::problem*) ru;

		//Step 2: unscale the decision vector
		Problem->unscale(x);

		//Step 3: call the fitness function
		try
		{
			Problem->evaluate(&(Problem->X[0]), F, G, *needG, Problem->iGfun, Problem->jGvar);
		}
		catch (int errorcode)
		{
			if (errorcode == 13)  //integration step error
				*Status = -1;
			if (errorcode == 1000000) //Kepler solver error
				*Status = -1;
		}

		//Step 4: If we have exceeded the allotted time, stop SNOPT by setting Status = -2 (Status < -1 causes SNOPT to quit)
		time_t now = time(NULL);
		if (now - (time_t) *iu > Problem->options.snopt_max_run_time)
		{
			std::cout << "Exceeded SNOPT time limit of " << Problem->options.snopt_max_run_time << " seconds. Aborting SNOPT run." << std::endl;
			*Status = -2;
		}

		return 0;
	}
}} //close namespace