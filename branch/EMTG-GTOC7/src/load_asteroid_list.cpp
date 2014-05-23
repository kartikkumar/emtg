//Donald Ellison May 23rd 2014

#include "lazy_race_tree_search.h"

namespace EMTG{
	void load_asteroid_list(std::string & asteroid_filename, std::vector <int> & asteroid_list)
	{

		std::ifstream inputfile(asteroid_filename.c_str());

		int number;

		if (!inputfile.is_open())
		{
			std::cout << "Cannot find asteroid file for lazy race tree search: " + asteroid_filename << std::endl;
		}

		while (!inputfile.eof())
		{
			inputfile >> number;
			asteroid_list.push_back(number);
		}

	}
}