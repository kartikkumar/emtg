//header for time utilities functions
#ifndef EMTG_TIME_UTILITIES
#define EMTG_TIME_UTILITIES

#include "SpiceUsr.h"
#include <string>

namespace EMTG 
{
    namespace time_utilities 
    {
        //ET to TDB conversion
        double convert_ET_to_TDB(const double& ETepoch);

        //ET to UTC conversion
        std::string convert_ET_to_UTC_string(const double& ETepoch);
    }
}

#endif //EMTG_TIME_UTILITIES