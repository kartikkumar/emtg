//Donald Ellison
//May 28th 2014
//Automatic .emtgopt and .qsub file generator for EMTG cluster runs

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>

#include "boost/filesystem.hpp"
//#include "boost/filesystem/fstream.hpp"
//#include "boost/date_time.hpp"
//#include "boost/date_time/local_time/local_date_time.hpp"
//#include "boost/lexical_cast.hpp"
//#include "boost/program_options.hpp"


int main(int argc, char* argv[])
{
	
	std::cout << "EMTG automatic options file generation starting" << std::endl << std::endl;

	//parse the options file
	std::string starting_body_ID_file_name = "ERROR";
	
	std::string options_and_qsub_files_directory = "..//Generated_Options_Files";
	
	std::vector <int> starting_body_ID_list; //NOTE, if we are making mothership launch options files, this is the list of destinations
	std::vector <double> epoch_list;

	if (argc == 1)
		starting_body_ID_file_name = "default.seed";
	else if (argc == 2)
		starting_body_ID_file_name.assign(argv[1]);

	
	if (starting_body_ID_file_name.compare("ERROR") == 0)
	{
		std::cout << "Aborting program run." << std::endl;
		std::cin.ignore();
		return 0;
	}

	std::cout << starting_body_ID_file_name << " read successfully" << std::endl << std::endl;


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//
	//
	// User specified stuff
	//
	//
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	double engine_duty_cycle;
	int objective_type;
	int MotherShipOpt;
	double MotherShipLaunchEpoch;
	int single_point_drop_off;

	std::cout << std::endl << std::endl << "EMTG Options File Entries:" << std::endl << std::endl << std::endl;

	std::cout << std::endl << "Enter an engine duty cycle [0.0, 1.0]:" << std::endl;
	std::cin >> engine_duty_cycle;

	std::cout << std::endl << "Enter an objective function type (0: Minimum deltaV, 1: Minimum TOF, 2: Min. propellant usage/adaptive TOF + min. prop. for lazy race tree)" << std::endl;
	std::cin >> objective_type;

	std::cout << std::endl << "Are you generating MotherShip launch files? (0: no, 1: yes)" << std::endl;
	std::cin >> MotherShipOpt;

	if (MotherShipOpt == 1)
	{
		std::cout << std::endl << "What is the mothership's launch window opening epoch (MJD)?" << std::endl;
		std::cin >> MotherShipLaunchEpoch; 
	}

	std::cout << std::endl << "Are you performing a single-point mothership drop-off? (0: no, 1: yes)" << std::endl;
	std::cin >> single_point_drop_off;

	std::cout << std::endl << std::endl << "Qsub File Entries:" << std::endl << std::endl << std::endl;

	std::string email;
	std::string queue;
	int nodes, ppn, walltime, cluster;
	bool BlueWaters = false;

	std::cout << std::endl << "Which cluster are you using? (1:Blue Waters or 2:Taub)" << std::endl;
	std::cin >> cluster;

	if (cluster == 1)
		BlueWaters = true;

	std::cout << std::endl << "What queue are you requesting (Blue Waters: low, normal or high Taub: cse or secondary)?" << std::endl;
	std::cin >> queue;

	std::cout << std::endl << "What is your email address? (Blue Waters only, enter anything for Taub)" << std::endl;
	std::cin >> email;

	std::cout << std::endl << "What is the requested walltime (integer number of minutes)?" << std::endl;
	std::cin >> walltime;

	std::cout << std::endl << "How many nodes are you requesting?" << std::endl;
	std::cin >> nodes;

	std::cout << std::endl << "How many cores per node? (Blue Waters: 16 max, taub: 12 max)" << std::endl;
	std::cin >> ppn;

	std::cout << std::endl;
	std::cin.ignore();





	//load the text file of starting_body_ID's and epochs if applicable
	bool read_epoch = true;
	int number;
	double epoch = 59215.0;

	if (MotherShipOpt)
		read_epoch = false;


	std::ifstream starting_body_ID_file(starting_body_ID_file_name.c_str());

	if (!starting_body_ID_file.is_open())
	{
		std::cout << "Cannot find starting body ID file for automatic EMTG option file generation: " + starting_body_ID_file_name << std::endl;
	}

	std::string ProbeName; //first line in seed file is probe #....1,2 or 3
	int ProbeID;

	if (single_point_drop_off)
	{
		starting_body_ID_file >> ProbeName;

		if (ProbeName == "Probe1")
			ProbeID = 1;
		else if (ProbeName == "Probe2")
			ProbeID = 2;
		else if (ProbeName == "Probe3")
			ProbeID = 3;
		else
		{
			std::cout << "Invalid probe name in the .seed file" << std::endl;
			return 0;
		}
	}
		
	while (!starting_body_ID_file.eof())
	{
		if (read_epoch == true)
		{
			starting_body_ID_file >> number; starting_body_ID_file >> epoch;
			starting_body_ID_list.push_back(number);
			epoch_list.push_back(epoch);
		}

		else
		{
			starting_body_ID_file >> number;
			starting_body_ID_list.push_back(number);
			epoch_list.push_back(epoch);
		}
	}
	
	
	std::string choice;
	std::string temp_line;
	std::string sub_temp_line;

	//create directory to dump the options files and qsub subdirectories to
	boost::filesystem::path p(options_and_qsub_files_directory);

	try
	{

		boost::filesystem::create_directories(p);
	}
	catch (std::exception &e)
	{
		std::cerr << "Error " << e.what() << ": Directory creation failed" << std::endl;
	}
	
	

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//
	//
	// EMTG OPTIONS FILE GENERATION
	//
	//
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	std::string options_file_name;
	
	if (MotherShipOpt)
		options_file_name = "Mommy_model.emtgopt";
	else
		options_file_name = "LRTS_model.emtgopt";

	std::cout << "Populating directories with options files" << std::endl;

	clock_t tStart = clock();

	std::string the_whole_part;
    std::string the_decimal_part;
	std::string epoch_as_string;
	std::ostringstream convert;
	int end_position;
	convert.precision(10);



	for (size_t index = 0; index < starting_body_ID_list.size(); ++index)
	{
		std::ifstream options_file(options_file_name.c_str(), std::ios::in);
		
		//std::string new_options_file_subdirectory = options_and_qsub_files_directory + "//" + "LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + std::to_string(int(epoch)) + "//";
		////create a subdirectory to dump the options file to
		//boost::filesystem::path p(new_options_file_subdirectory);
		//try
		//{

		//	boost::filesystem::create_directories(p);
		//}
		//catch (std::exception &e)
		//{
		//	std::cerr << "Error " << e.what() << ": Directory creation failed" << std::endl;
		//}
		convert.str("");
		the_whole_part.clear();
		the_decimal_part.clear();
		epoch_as_string.clear();
		convert << epoch_list[index];
		epoch_as_string = convert.str();
		
		the_whole_part.assign(epoch_as_string.substr(0, epoch_as_string.find('.')));

		end_position = epoch_as_string.end() - epoch_as_string.begin();
		the_decimal_part.assign(epoch_as_string.substr(epoch_as_string.find('.')+1, end_position));

		std::string new_options_file_name;

		if (MotherShipOpt)
			new_options_file_name = options_and_qsub_files_directory + "//" + std::to_string(starting_body_ID_list[index]) + "_Mommy.emtgopt";
		else if (single_point_drop_off)
			new_options_file_name = options_and_qsub_files_directory + "//" + std::to_string(index) + "_" + std::to_string(ProbeID) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".emtgopt";
		else
			new_options_file_name = options_and_qsub_files_directory + "//" + std::to_string(index) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".emtgopt";
		
		
		std::ofstream new_options_file(new_options_file_name.c_str(), std::ios::out);

		new_options_file.precision(10);
		while (!options_file.eof())
		{
			//read the current line from file and store it temporarily
			std::getline(options_file, temp_line); //default delimiter is newline

			//grab the variable name only
			sub_temp_line.assign(temp_line.substr(0, temp_line.find(' ')));

			if (sub_temp_line.compare("lazy_race_tree_start_location_ID") == 0)
			{
				//we have found the line of interest, modify it and write it to the new options file
				new_options_file << std::left << sub_temp_line << ' ' << std::to_string(starting_body_ID_list[index]) << std::endl;
			}
			else if (single_point_drop_off == 1 && sub_temp_line.compare("lazy_race_tree_target_list_file") == 0)
			{
				if (ProbeID == 3)
					new_options_file << std::left << sub_temp_line << ' ' << "../Koronis_high.asteroidlist" << std::endl;
				else if (ProbeID == 2)
					new_options_file << std::left << sub_temp_line << ' ' << "../Koronis_mid.asteroidlist" << std::endl;
				else if (ProbeID == 1)
					new_options_file << std::left << sub_temp_line << ' ' << "../Koronis_low.asteroidlist" << std::endl;
			}
			else if (sub_temp_line.compare("mission_name") == 0 && MotherShipOpt == 1)
			{
				new_options_file << std::left << sub_temp_line << ' ' << std::to_string(starting_body_ID_list[index]) + "_Mommy" << std::endl;
			}
			else if (sub_temp_line.compare("mission_name") == 0 && MotherShipOpt == 0)
			{
				new_options_file << std::left << sub_temp_line << ' ' << std::to_string(index) + "_" + std::to_string(ProbeID) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part << std::endl;
			}
			else if (sub_temp_line.compare("objective_type") == 0)
			{
				new_options_file << std::left << sub_temp_line << ' ' << objective_type << std::endl;
			}
			else if (sub_temp_line.compare("launch_window_open_date") == 0)
			{
				if (MotherShipOpt == 0)
					new_options_file << std::left << sub_temp_line << ' ' << std::to_string(epoch_list[index]) << std::endl;
				else if (MotherShipOpt == 1)
					new_options_file << std::left << sub_temp_line << ' ' << std::to_string(MotherShipLaunchEpoch) << std::endl;
			}
			else if (sub_temp_line.compare("engine_duty_cycle") == 0)
			{
				new_options_file << std::left << sub_temp_line << ' ' << std::to_string(engine_duty_cycle) << std::endl;
			}
			else if (sub_temp_line.compare("global_timebounded") == 0)
			{
				if (MotherShipOpt)
					new_options_file << std::left << sub_temp_line << ' ' << '1' << std::endl;
				else
				{
					if (objective_type == 1)
						new_options_file << std::left << sub_temp_line << ' ' << '0' << std::endl;
					else if (objective_type == 2)
						new_options_file << std::left << sub_temp_line << ' ' << '1' << std::endl;
				}
			}
			else if (sub_temp_line.compare("destination_list") == 0)
			{
				new_options_file << std::left << sub_temp_line << " 1 " << std::to_string(starting_body_ID_list[index]) << std::endl;
			}
			else
			{
				//write the temp line to the new options file unaltered
				new_options_file << std::left << temp_line << std::endl;
			}

		}

		options_file.close();
		new_options_file.close();
 	}

	double elapsed_time = double((clock() - tStart)) / CLOCKS_PER_SEC;
	std::cout.precision(4);
	std::cout << "Options file generation completed in: " << std::to_string(elapsed_time) << " s" << std::endl << std::endl;
	

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//
	//
	// QSUB FILE GENERATION
	//
	//
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//Now we want to modify and write a qsub to the same directory that we just wrote the options file to
	std::cout << "Populating directories with qsub files" << std::endl;
	std::string qsub_file_name;

	if (BlueWaters)
		qsub_file_name = "LRTS_model_BlueWaters.qsub";
	else if (MotherShipOpt)
		qsub_file_name = "Mommy_model_Taub.qsub";
	else
		qsub_file_name = "LRTS_model_Taub.qsub";

	tStart = clock();

	for (size_t index = 0; index < starting_body_ID_list.size(); ++index)
	{ 
		std::ifstream qsub_file(qsub_file_name.c_str(), std::ios::in);
		//std::string new_qsub_file_subdirectory = options_and_qsub_files_directory + "//" + "LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + std::to_string(int(epoch)) + "//";

		convert.str("");
		the_whole_part.clear();
		the_decimal_part.clear();
		epoch_as_string.clear();
		convert << epoch_list[index];
		epoch_as_string = convert.str();

		the_whole_part.assign(epoch_as_string.substr(0, epoch_as_string.find('.')));

		end_position = epoch_as_string.end() - epoch_as_string.begin();
		the_decimal_part.assign(epoch_as_string.substr(epoch_as_string.find('.') + 1, end_position));
		
		std::string new_qsub_file_name;

		if (MotherShipOpt)
			new_qsub_file_name = options_and_qsub_files_directory + "//" + std::to_string(starting_body_ID_list[index]) + "_Mommy.qsub";
		else if (single_point_drop_off)
			new_qsub_file_name = options_and_qsub_files_directory + "//" + std::to_string(index) + "_" + std::to_string(ProbeID) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".qsub";
		else
			new_qsub_file_name = options_and_qsub_files_directory + "//" + std::to_string(index) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".qsub";
		

		std::ofstream new_qsub_file(new_qsub_file_name.c_str(), std::ios::out);

		while (!qsub_file.eof())
		{
			//read the current line from file and store it temporarily
			std::getline(qsub_file, temp_line); //default delimiter is newline

			//grab the variable name only
			sub_temp_line.assign(temp_line.substr(0, temp_line.find(' ')));

			if (BlueWaters)
			{
				if (temp_line.compare("#PBS -M email") == 0)
				{
					new_qsub_file << std::left << "#PBS -M " + email << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -q queue") == 0)
				{
					new_qsub_file << std::left << "#PBS -q " + queue << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -l nodes=nodenum:ppn=corenum:xe") == 0)
				{
					//we have found the line of interest, modify it and write it to the new options file
					new_qsub_file << std::left << "#PBS -l nodes=" + std::to_string(nodes) + ":ppn=" + std::to_string(ppn) + ":xe" << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -l walltime=00:min:00") == 0)
				{
					new_qsub_file << std::left << "#PBS -l walltime=00:" + std::to_string(walltime) + ":00" << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -N jobname") == 0)
				{
					if (MotherShipOpt)
					{
						new_qsub_file << std::left << "#PBS -N " + std::to_string(starting_body_ID_list[index]) + "_Mommy" << std::endl;
					}
					else if (single_point_drop_off)
					{
						new_qsub_file << std::left << "#PBS -N " + std::to_string(index) + "_" + std::to_string(ProbeID) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part << std::endl;
					}
					else
					{
						new_qsub_file << std::left << "#PBS -N " + std::to_string(index) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part << std::endl;
					}
					continue;
				}
				else if (sub_temp_line.compare("aprun") == 0)
				{
					//we have found the line of interest, modify it and write it to the new options file
					if (MotherShipOpt)
						new_qsub_file << std::left << sub_temp_line << " -n " + std::to_string(nodes*ppn) + " ../emtg " << std::to_string(starting_body_ID_list[index]) + "_Mommy.emtgopt" << std::endl;
					else if (single_point_drop_off)
						new_qsub_file << std::left << sub_temp_line << " -n " + std::to_string(nodes*ppn) + " ../emtg " << std::to_string(index) + "_" + std::to_string(ProbeID) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".emtgopt" << std::endl;
					else
						new_qsub_file << std::left << sub_temp_line << " -n " + std::to_string(nodes*ppn) + " ../emtg " << std::to_string(index) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".emtgopt" << std::endl;
					continue;
				}

			}
			else if (MotherShipOpt)
			{
				if (sub_temp_line.compare("../emtg") == 0)
				{
					//we have found the line of interest, modify it and write it to the new options file
					new_qsub_file << std::left << sub_temp_line << ' ' <<  std::to_string(starting_body_ID_list[index]) + "_Mommy.emtgopt" << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -q queue") == 0)
				{
					new_qsub_file << std::left << "#PBS -q " + queue << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -l nodes=nodenum:ppn=corenum") == 0)
				{
					//we have found the line of interest, modify it and write it to the new options file
					new_qsub_file << std::left << "#PBS -l nodes=" + std::to_string(nodes) + ":ppn=" + std::to_string(ppn) << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -N jobname") == 0)
				{
					new_qsub_file << std::left << "#PBS -N " + std::to_string(starting_body_ID_list[index]) + "_Mommy" << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -l walltime=00:min:00") == 0)
				{
					new_qsub_file << std::left << "#PBS -l walltime=00:" + std::to_string(walltime) + ":00" << std::endl;
					continue;
				}
			}
			else
			{
				if (sub_temp_line.compare("mpiexec") == 0)
				{
					//we have found the line of interest, modify it and write it to the new options file
					if (single_point_drop_off)
						new_qsub_file << std::left << sub_temp_line << " ../emtg " << std::to_string(index) + "_" + std::to_string(ProbeID) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".emtgopt" << std::endl;
					else
						new_qsub_file << std::left << sub_temp_line << " ../emtg " << std::to_string(index) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part + ".emtgopt" << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -q queue") == 0)
				{
					new_qsub_file << std::left << "#PBS -q " + queue << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -l nodes=nodenum:ppn=corenum") == 0)
				{
					//we have found the line of interest, modify it and write it to the new options file
					new_qsub_file << std::left << "#PBS -l nodes=" + std::to_string(nodes) + ":ppn=" + std::to_string(ppn) << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -N jobname") == 0)
				{
					if (single_point_drop_off)
						new_qsub_file << std::left << "#PBS -N " + std::to_string(index) + "_" + std::to_string(ProbeID) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part << std::endl;
					else	
						new_qsub_file << std::left << "#PBS -N " + std::to_string(index) + "_LRTS_" + std::to_string(starting_body_ID_list[index]) + '_' + the_whole_part + 'd' + the_decimal_part << std::endl;
					continue;
				}
				else if (temp_line.compare("#PBS -l walltime=00:min:00") == 0)
				{
					new_qsub_file << std::left << "#PBS -l walltime=00:" + std::to_string(walltime) + ":00" << std::endl;
					continue;
				}
			}
			
			//write the temp line to the new options file unaltered
			new_qsub_file << std::left << temp_line << std::endl;
			

		}

		qsub_file.close();
		new_qsub_file.close();
	}

	elapsed_time = double((clock() - tStart)) / CLOCKS_PER_SEC;
	std::cout.precision(4);
	std::cout << "Qsub file generation completed in: " << std::to_string(elapsed_time) << " s" << std::endl << std::endl << std::endl;

	std::cout << "Generation complete, hit enter to exit." << std::endl;
	getchar();

	return 0;

}