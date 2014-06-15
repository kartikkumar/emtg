#include "EMTG_string_utilities.h"

namespace EMTG {
	namespace string_utilities {

		//adapted from from http://stackoverflow.com/questions/7132957/c-scientific-notation-format-number
		std::string convert_number_to_formatted_string(const double& number, int expSize)
		{
			std::ostringstream oss;
			oss.precision(14);
			std::string output;
			oss << std::scientific << number;
			std::string numstring = oss.str();
			unsigned int ePos = numstring.find("e");
			numstring.replace(ePos, 1, "E");
			unsigned int dPos = numstring.find(".");
			if (ePos == 0)
			{
				//no exponent
				return numstring;
			}
			else if (dPos == 0)
			{
				//not decimal
				return numstring;
			}
			else
			{
				if (number < 0.0)
					output = "";
				else
					output = " ";
				output += numstring.substr(0, ePos);

				while (output.size() < 17)
					output += "0";

				output += numstring.substr(ePos, 2);
				if (numstring.size() - expSize > ePos + 1)
					output += numstring.substr(numstring.size() - expSize, numstring.size());
				else
				{
					//expSize too big (or bug -> e used but no exponent?)
				}
				return output;
			}
		}
	}
}
