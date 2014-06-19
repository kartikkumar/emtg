//code to convert from MJD to month, day, year, hours, minutes, seconds
//Brent Barbee
//added to EMTG 7-24-2012

namespace EMTG
{
#ifndef _MJD_TO_MDYHMS
#define _MJD_TO_MDYHMS


void mjd_to_mdyhms(double MJD, int *month, int *day, int *year, int *hrs, int *mins, double *secs);


#endif

}