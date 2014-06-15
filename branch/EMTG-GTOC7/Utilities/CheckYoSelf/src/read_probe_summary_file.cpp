#include "GTOC7_solution_check.h"


Spacecraft read_probe_summary_file(std::string & probe_summary_file_name)
{
	Spacecraft probe;
	std::ifstream probe_summary_file(probe_summary_file_name.c_str(), std::ios::in);

	if (!probe_summary_file.is_open())
		std::cout << "Error in probe summary file reader. Cannot find probe summary file: " + probe_summary_file_name << std::endl;

	std::vector <int> asteroid_sequence;
	std::string peek;
	char dump_buffer[1024];
	int linenumber = 0;
	int end_position;
	int probe_num;
	int asteroid_num;

	std::string temp_line;
	std::string sub_temp_line;
	std::stringstream convert;
	convert.clear();
	convert.str(std::string());

	//Probe number
	std::getline(probe_summary_file, temp_line);
	end_position = temp_line.end() - temp_line.begin();
	sub_temp_line.assign(temp_line.substr(end_position - 1, end_position));
	convert.str(sub_temp_line);
	convert >> probe_num;
	probe.name = "Probe " + std::to_string(probe_num);

	//First phase number -- don't care
	probe_summary_file.getline(dump_buffer, 1024);

	//First asteroid number
	convert.clear();
	convert.str(std::string());
	std::getline(probe_summary_file, temp_line);
	end_position = temp_line.end() - temp_line.begin();
	sub_temp_line.assign(temp_line.substr(temp_line.find("Asteroid")+9, end_position));
	convert.str(sub_temp_line);
	convert >> asteroid_num;
	asteroid_sequence.push_back(asteroid_num);

	//Header line -- don't care
	probe_summary_file.getline(dump_buffer, 1024);

	int phase = 0;
	int timestep = 0;
	double MJD, x, y, z, vx, vy, vz, mass, Tx, Ty, Tz;
	std::vector <double> newphase;

	//add a new row (phase) to the state history
	probe.time_stamp.push_back(newphase);
	probe.x.push_back(newphase);
	probe.y.push_back(newphase);
	probe.z.push_back(newphase);
	probe.vx.push_back(newphase);
	probe.vy.push_back(newphase);
	probe.vz.push_back(newphase);
	probe.mass.push_back(newphase);
	probe.Tx.push_back(newphase);
	probe.Ty.push_back(newphase);
	probe.Tz.push_back(newphase);
	probe.a.push_back(newphase);
	probe.ecc.push_back(newphase);
	probe.inc.push_back(newphase);
	probe.LAN.push_back(newphase);
	probe.omega.push_back(newphase);
	probe.M.push_back(newphase);
	probe.tru.push_back(newphase);
	probe.E.push_back(newphase);
	probe.h.push_back(newphase);
	probe.r.push_back(newphase);
	probe.v.push_back(newphase);
	

	while (!probe_summary_file.eof())
	{
		

		peek = probe_summary_file.peek();
		if (peek == "\n")
			dump_buffer[0] = probe_summary_file.get();

		peek = probe_summary_file.peek();
		if (peek == "#")
		{
			//Phase line -- don't care
			probe_summary_file.getline(dump_buffer, 1024);
			++linenumber;

			//Starting asteroid number for the phase
			convert.clear();
			convert.str(std::string());
			std::getline(probe_summary_file, temp_line);
			end_position = temp_line.end() - temp_line.begin();

			sub_temp_line.assign(temp_line.substr(temp_line.find("A") + 9, end_position));
			sub_temp_line.assign(sub_temp_line.substr(0, sub_temp_line.find(' ')));
			convert.str(sub_temp_line);
			convert >> asteroid_num;
			asteroid_sequence.push_back(asteroid_num);

			sub_temp_line.assign(temp_line.substr(temp_line.find_last_of("A") + 9, end_position));
			convert.clear();
			convert.str(std::string());
			convert.str(sub_temp_line);
			convert >> asteroid_num;
			asteroid_sequence.push_back(asteroid_num);

			//Header line -- don't care
			probe_summary_file.getline(dump_buffer, 1024);
			++linenumber;

			//add a new row (phase) to the state history
			probe.time_stamp.push_back(newphase);
			probe.x.push_back(newphase);
			probe.y.push_back(newphase);
			probe.z.push_back(newphase);
			probe.vx.push_back(newphase);
			probe.vy.push_back(newphase);
			probe.vz.push_back(newphase);
			probe.mass.push_back(newphase);
			probe.Tx.push_back(newphase);
			probe.Ty.push_back(newphase);
			probe.Tz.push_back(newphase);
			probe.a.push_back(newphase);
			probe.ecc.push_back(newphase);
			probe.inc.push_back(newphase);
			probe.LAN.push_back(newphase);
			probe.omega.push_back(newphase);
			probe.M.push_back(newphase);
			probe.tru.push_back(newphase);
			probe.E.push_back(newphase);
			probe.h.push_back(newphase);
			probe.r.push_back(newphase);
			probe.v.push_back(newphase);

			++phase;

			timestep = 0;
		}


		//read current trajectory data line from probe summary file
		probe_summary_file >> MJD;
		probe_summary_file >> x;
		probe_summary_file >> y;
		probe_summary_file >> z;
		probe_summary_file >> vx;
		probe_summary_file >> vy;
		probe_summary_file >> vz;
		probe_summary_file >> mass;
		probe_summary_file >> Tx;
		probe_summary_file >> Ty;
		probe_summary_file >> Tz;

		//load the current state into the probe structure
		probe.time_stamp[phase].push_back(MJD*86400.0);
		probe.x[phase].push_back(x);
		probe.y[phase].push_back(y);
		probe.z[phase].push_back(z);
		probe.vx[phase].push_back(vx);
		probe.vy[phase].push_back(vy);
		probe.vz[phase].push_back(vz);
		probe.mass[phase].push_back(mass);
		probe.Tx[phase].push_back(Tx);
		probe.Ty[phase].push_back(Ty);
		probe.Tz[phase].push_back(Tz);

		

		cartesian2coe_probe(&probe, phase, timestep);

		++timestep;

	}

	return probe;

}