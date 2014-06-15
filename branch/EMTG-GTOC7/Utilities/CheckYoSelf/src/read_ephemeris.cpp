#include "GTOC7_solution_check.h"

std::vector <Body> read_ephemeris(std::string & GTOC7_ephemeris_file_name)
{
	
	std::ifstream ephemeris_file(GTOC7_ephemeris_file_name.c_str(), std::ios::in);
	std::vector <Body> asteroids;

	if (!ephemeris_file.is_open())
		std::cout << "Error in read_ephemeris. Cannot find GTOC7 ephemeris file: " + GTOC7_ephemeris_file_name << std::endl;
	
	//Grab the first two GTOC7 asteroid file header lines and toss them
	char dump_buffer[1024];
	ephemeris_file.getline(dump_buffer, 1024);
	ephemeris_file.getline(dump_buffer, 1024);

	int index = 0;
	char blank;

	while (!ephemeris_file.eof())
	{
		//add a new blank asteroid with the default constructor
		asteroids.push_back(Body());

		ephemeris_file >> asteroids[index].GTOC7_num;
		ephemeris_file >> asteroids[index].time_stamp;
		ephemeris_file >> asteroids[index].a;
		ephemeris_file >> asteroids[index].ecc;
		ephemeris_file >> asteroids[index].inc;
		ephemeris_file >> asteroids[index].omega;
		ephemeris_file >> asteroids[index].LAN;
		ephemeris_file >> asteroids[index].M;

		//convert to km and radians
		asteroids[index].a *= 149597870.691;
		asteroids[index].inc *= PI / 180.0;
		asteroids[index].omega *= PI / 180.0;
		asteroids[index].LAN *= PI / 180.0;
		asteroids[index].M *= PI / 180.0;

		//convert MJD to seconds
		asteroids[index].time_stamp *= 86400.0;

		//read the asteroid name, then strip the leading tab character
		std::getline(ephemeris_file, asteroids[index].name);
		asteroids[index].name.erase(asteroids[index].name.begin());

		//now calculate E, f and the cartesian state for the asteroid
		asteroids[index].E = laguerreConway(asteroids[index].ecc, asteroids[index].M);
		asteroids[index].tru = 2.0*atan(sqrt((1.0 + asteroids[index].ecc) / (1.0 - asteroids[index].ecc))*tan(asteroids[index].E / 2.0));
		coe2cartesian(&asteroids[index]);

		++index;
	}
	
	
	ephemeris_file.close();

	return asteroids;
}
