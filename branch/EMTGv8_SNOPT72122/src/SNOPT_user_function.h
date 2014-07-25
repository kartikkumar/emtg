//SNOPT user function

#include "snopt.hh"


namespace EMTG { namespace Solvers {

	int SNOPT_user_function(integer    *Status, integer *n,    doublereal x[],
								  integer    *needF,  integer *neF,  doublereal F[],
								  integer    *needG,  integer *neG,  doublereal G[],
								  char       *cu,     integer *lencu,
								  integer    iu[],    integer *leniu,
								  doublereal ru[],    integer *lenru );
}}