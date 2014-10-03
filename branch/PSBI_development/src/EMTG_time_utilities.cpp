//time utilities functions

#include "EMTG_time_utilities.h"
#include <cmath>

namespace EMTG 
{
    namespace time_utilities 
    {
        //convert from ET to TDB
        //adapted from http://sourcecodebrowser.com/astronomical-almanac/5.6/tdb_8c.html
        double time_utilities::convert_ET_to_TDB(const double& ETepoch)
        {
            double M, T;

            /* Find time T in Julian centuries from J2000.  */
            T = (ETepoch - 2451545.0) / 36525.0;

            /* Mean anomaly of sun = l' (J. Laskar) */
            M = 129596581.038354 * T + 1287104.76154;

            /* Reduce arc seconds mod 360 degrees.  */
            M = M - 1296000.0 * floor(M / 1296000.0);

            M += ((((((((
                1.62e-20 * T
                - 1.0390e-17) * T
                - 3.83508e-15) * T
                + 4.237343e-13) * T
                + 8.8555011e-11) * T
                - 4.77258489e-8) * T
                - 1.1297037031e-5) * T
                + 1.4732069041e-4) * T
                - 0.552891801772) * T * T;

            M *= 4.8481368110953599359e-6;
            /* TDB - TDT, in seconds.  */
            T = 0.001658 * sin(M) + 0.000014 * sin(M + M);

            T = ETepoch + T / 86400.0;

            return T;
        }

        //function to convert ET double-precision value to UTC string (in STK-happy format)
        std::string convert_ET_to_UTC_string(const double& ETepoch)
        {
            SpiceChar UTC_char[25];
            timout_c(   ETepoch, 
                        "DD MON YYYY HR:MN:SC.### :UTC",
                        25,
                        UTC_char);
            
            std::string UTC_string(UTC_char);

            return UTC_string;
        }
    }//close namespace time_utilities
}//close namespace EMTG