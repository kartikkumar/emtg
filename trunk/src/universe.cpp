//header file for EMTG body class

#include <string>
#include <vector>
#include <fstream>

#include "missionoptions.h"
#include "body.h"
#include "universe.h"
#include "EMTG_math.h"

#include "SpiceUsr.h"



#include "boost/algorithm/string.hpp"

using namespace std;

namespace EMTG {namespace Astrodynamics {
	//default constructor does nothing
	universe::universe(){}

	//constructor to load a data file
	universe::universe(const int& j, string universe_file, missionoptions* options)
	{
		load_universe_data(j, universe_file, options);
	}

	//destructor
	universe::~universe(){}

	//function to load a data file
	int universe::load_universe_data(const int& j, string universefile, missionoptions* options)
	{
		ifstream inputfile(universefile.c_str());
		int linenumber = 0;
		string choice;
		string peek;
		double value;
		char dump_buffer[1024];

		if (!inputfile.is_open())
		{
			cout << "Failure to read " << universefile << endl;
			return 2;
		}

		while (!inputfile.eof())
		{
			peek = inputfile.peek();
			if (peek == "#" || peek == "\r" || peek == "\n") {
				//comment or blank line, do not parse
				inputfile.getline(dump_buffer, 1024);
				++linenumber;
			}
			else 
			{
				inputfile >> choice;

				if (choice == "central_body_name")
					inputfile >> central_body_name;
				else if (choice == "central_body_SPICE_ID")
					inputfile >> central_body_SPICE_ID;
				else if (choice == "mu")
				{
					inputfile >> value;
					mu = value;
				}
				else if (choice == "central_body_radius")
				{
					inputfile >> value;
					central_body_radius = value;
				}
				else if (choice == "LU")
				{
					inputfile >> value;
					LU = value;
				}
				else if (choice == "reference_angles")
				{
					double reference_angles[6];
					inputfile >> reference_angles[0];
					inputfile >> reference_angles[1];
					inputfile >> reference_angles[2];
					inputfile >> reference_angles[3];
					inputfile >> reference_angles[4];
					inputfile >> reference_angles[5];
					
					if (central_body_name == "EARTH")
						LocalFrame.initialize();
					else
						LocalFrame.initialize(reference_angles[0] * math::PI / 180.0, reference_angles[1] * math::PI / 180.0, reference_angles[2] * math::PI / 180.0, reference_angles[3] * math::PI / 180.0, reference_angles[4] * math::PI / 180.0, reference_angles[5] * math::PI / 180.0);
				}
				else if (choice == "r_SOI")
				{
					inputfile >> value;
					r_SOI = value;
				}
				else if (choice == "minimum_safe_distance")
				{
					inputfile >> value;
					minimum_safe_distance = value;
				}		
				else if (choice == "begin_body_list")
				{
					string tempname;
					string tempshortname;
					int temp_bodycode;
					int temp_SPICEnum;
					double temp_minimum_altitude;
					double temp_mass;
					double temp_radius;
					double temp_epoch;
					vector<double> temp_reference_angles(6);
					vector<double> temp_elements(6);

					do
					{
						tempname.clear();

						//if not yet reached end of list:
						inputfile >> tempname;
						inputfile >> tempshortname;
						inputfile >> temp_bodycode;
						if (tempname == "end_body_list")
							break;
						inputfile >> temp_SPICEnum;
						inputfile >> temp_minimum_altitude;
						inputfile >> temp_mass;
						inputfile >> temp_radius;
						inputfile >> temp_epoch;
						temp_epoch *= 86400.0;

						for (int k = 0; k < 6; ++k)
						{
							inputfile >> temp_reference_angles[k];
							temp_reference_angles[k] *= math::PI / 180.0;
						}

						for (int k = 0; k < 6; ++k)
							inputfile >> temp_elements[k];

						bodies.push_back(new body(temp_bodycode, tempname, tempshortname, temp_SPICEnum, temp_minimum_altitude, temp_mass, temp_radius, temp_epoch, temp_reference_angles, temp_elements, mu, central_body_SPICE_ID, central_body_name, central_body_radius, options));

					} while(true);
				}

			}
		}

		//compute derived quantities
		TU = sqrt(LU*LU*LU/mu);

		//create the flyby menu
		create_flyby_and_perturbation_menus(j, options);

		return 0;
	}

	//function to locate the central body relative to the sun
	int universe::locate_central_body(const double& epoch, double* state, missionoptions* options)
	{
		if (!(boost::to_upper_copy(this->central_body_name) == "SUN"))
		{
			double LT_dump;
			spkez_c (central_body_SPICE_ID, epoch - (51544.5 * 86400.0), "J2000", "NONE", 10, state, &LT_dump);
		}
		else
		{
			state[0] = 0.0;
			state[1] = 0.0;
			state[2] = 0.0;
			state[3] = 0.0;
			state[4] = 0.0;
			state[5] = 0.0;
		}

		return 0;
	}

	//function to create the flyby menu - creates a list of bodies, by SPICE ID, which are flyby capable
	void universe::create_flyby_and_perturbation_menus(const int& j, missionoptions* options)
	{
		for (size_t k = 0; k < bodies.size(); ++k)
		{
			//flyby menu
			if (bodies[k].minimum_safe_flyby_altitude > 0.0)
				flyby_menu.push_back(k); //store the position in the universe list for flyby-capable bodies
			
			//perturbation menu
			if (options->perturb_thirdbody)
			{
				for (size_t b = 0; b < options->journey_perturbation_bodies[j].size(); ++b)
				{
					if (bodies[k].body_code == options->journey_perturbation_bodies[j][b])
						perturbation_menu.push_back(k); //store the position in the universe list for perturbation bodies
				}
			}
		}

		size_of_flyby_menu = flyby_menu.size() * 2;
	}

	//function to print the flyby menu
	void universe::print_flyby_and_perturbation_menus(string filename, missionoptions* options)
	{
		ofstream outputfile (filename.c_str(), ios::app);

		outputfile << endl;
		outputfile << "Flyby menu:" << endl;
		outputfile << "Position in flyby list | Position in Universe List | Name         | SPICE ID" << endl;
		outputfile << "----------------------------------------------------------------------------" << endl;
		for (int k = 0; k < size_of_flyby_menu/2; ++k)
		{
			outputfile.width(29);
			outputfile << left << k+1 << " | " ;
			outputfile.width(25);
			outputfile << left << flyby_menu[k] + 1 << " | ";
			outputfile.width(11);
			outputfile << bodies[flyby_menu[k]].name;
			outputfile.width(3);
			outputfile << " | " << bodies[flyby_menu[k]].spice_ID << endl;
		}

		if (options->perturb_thirdbody)
		{
			outputfile << endl;
			outputfile << "Perturbation menu:" << endl;
			outputfile << "Position in perturbation list | Position in Universe List | Name         | SPICE ID" << endl;
			outputfile << "-----------------------------------------------------------------------------------" << endl;
			for (size_t k = 0; k < perturbation_menu.size(); ++k)
			{
				outputfile.width(29);
				outputfile << left << k+1 << " | " ;
				outputfile.width(25);
				outputfile << left << perturbation_menu[k] + 1 << " | ";
				outputfile.width(11);
				outputfile << bodies[perturbation_menu[k]].name;
				outputfile.width(3);
				outputfile << " | " << bodies[perturbation_menu[k]].spice_ID << endl;
			}
		}

		outputfile.close();
	}

	//function to print the universe to a file
	void universe::print_universe(string filename, missionoptions* options)
	{
		ofstream outputfile(filename.c_str(), ios::trunc);

		outputfile << "Central body name: " << central_body_name << endl;
		outputfile << "Central body SPICE ID: " << central_body_SPICE_ID << endl;
		outputfile << "mu (km^3/s^2) = " << mu << endl;
		outputfile << "LU (km) = " << LU << endl;
		outputfile << "TU (s) = " << TU << endl;
		outputfile << "alpha0 (degrees) = " << LocalFrame.alpha0 * 180.0 / math::PI << endl;
		outputfile << "alphadot (degrees/century) = " << LocalFrame.alphadot * 180.0 / math::PI << endl;
		outputfile << "delta0 (degrees) = " << LocalFrame.delta0 * 180.0 / math::PI << endl;
		outputfile << "deltadot (degrees/century) = " << LocalFrame.deltadot * 180.0 / math::PI << endl;
		outputfile << "Sphere of influence radius (km) = " << r_SOI << endl;
		outputfile << "Minimum safe distance (km) = " << minimum_safe_distance << endl;
		outputfile << endl;
		outputfile << "Bodies:" << endl;
		outputfile << endl;

		outputfile.close();

		for (size_t k = 0; k < bodies.size(); ++k)
			bodies[k].print_body_to_screen(filename);

		print_flyby_and_perturbation_menus(filename, options);
	}

	
}}//close namespace