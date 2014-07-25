#include "EMTG.h"
#include "options.h"
#include "snopt.h"

#ifndef _EMTG_SOLVERS
#define _EMTG_SOLVERS

int EMTG_usrf_(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
	     integer    iu[],    integer *leniu,
	     doublereal ru[],    integer *lenru );

int EMTG_pseudo_SF_usrf_(integer    *Status, integer *n,    doublereal x[],
						  integer    *needF,  integer *neF,  doublereal F[],
						  integer    *needG,  integer *neG,  doublereal G[],
						  char       *cu,     integer *lencu,
						  integer    iu[],    integer *leniu,
						  doublereal ru[],    integer *lenru );

int EMTG_SF_usrf_(integer    *Status, integer *n,    doublereal x[],
						  integer    *needF,  integer *neF,  doublereal F[],
						  integer    *needG,  integer *neG,  doublereal G[],
						  char       *cu,     integer *lencu,
						  integer    iu[],    integer *leniu,
						  doublereal ru[],    integer *lenru );

#endif //EMTG_SOLVERS