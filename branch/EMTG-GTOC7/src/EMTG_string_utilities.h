#ifndef EMTG_STRING_H_
#define EMTG_STRING_H_

#include <string>
#include <sstream>

namespace EMTG {
	namespace string_utilities {

		std::string convert_number_to_formatted_string(const double& number, int expSize);
	}
}

#endif //EMTG_STRING_H_