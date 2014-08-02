//SNOPT user function

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
							SNOPT_DOUBLE_TYPE ru[], SNOPT_INT_TYPE *lenru);
}}