//Donald Ellison May 23rd 2014

#include "lazy_race_tree_search.h"

namespace EMTG{

	void write_branch_summary(EMTG::mission & branch_mission, EMTG::missionoptions & branch_options, std::string & tree_summary_file_location, int & tree_level, const int & branch)
	{
		//open the summary file for appending
		std::ofstream outputfile(tree_summary_file_location.c_str(), std::ios::app);

		//Write things of interest for this branch
		outputfile.width(15); outputfile << left << tree_level;
		outputfile.width(15); outputfile << branch;
		outputfile.width(35); outputfile << branch_options.destination_list[0][0];
		outputfile.width(35); outputfile << branch_options.destination_list[0][1];
		outputfile.width(25); outputfile.precision(10); outputfile << branch_options.launch_window_open_date;
		outputfile.width(20); outputfile.precision(10); outputfile << branch_mission.Xopt[0];
		outputfile.width(20); outputfile.precision(10); outputfile << branch_mission.Xopt[1];
		outputfile.width(20); outputfile.precision(10); outputfile << branch_mission.Xopt[0] + branch_mission.Xopt[1];
		outputfile.width(20); outputfile << branch_options.maximum_mass;
		outputfile.width(20); outputfile << branch_mission.Xopt[2];

		
		outputfile << std::endl;
		//close the tree summary file
		outputfile.close();
	}

}