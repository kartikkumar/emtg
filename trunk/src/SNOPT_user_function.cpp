//SNOPT user function

#include <ctime>
#include <iostream>

#include "problem.h"

#include "snoptProblemExtension.h"


namespace EMTG { namespace Solvers {

#ifdef Heritage_SNOPT7
	int
#else
	void
#endif
		SNOPT_user_function(SNOPT_INT_TYPE    *Status, SNOPT_INT_TYPE *n, SNOPT_DOUBLE_TYPE x[],
							SNOPT_INT_TYPE    *needF, SNOPT_INT_TYPE *neF, SNOPT_DOUBLE_TYPE F[],
							SNOPT_INT_TYPE    *needG, SNOPT_INT_TYPE *neG, SNOPT_DOUBLE_TYPE G[],
							char       *cu, SNOPT_INT_TYPE *lencu,
							SNOPT_INT_TYPE    iu[], SNOPT_INT_TYPE *leniu,
							SNOPT_DOUBLE_TYPE ru[], SNOPT_INT_TYPE *lenru)
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

#ifdef Heritage_SNOPT7
		return 0;
#else
		return;
#endif
	}
}} //close namespace