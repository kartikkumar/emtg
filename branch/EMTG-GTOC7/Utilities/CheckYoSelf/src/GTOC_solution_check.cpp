//Donald Ellison
//June 14th 2014
//This program verifies GTOC7 output files

#include "GTOC7_solution_check.h"


int main(int argc, char *argv[])
{

	std::string probe_summary_file_name = "ERROR";

	if (argc == 1)
		probe_summary_file_name = "default.GTOC7";
	else if (argc == 2)
		probe_summary_file_name.assign(argv[1]);


	if (probe_summary_file_name.compare("ERROR") == 0)
	{
		std::cout << "Input file is not valid, aborting program run." << std::endl;
		std::cin.ignore();
		return 0;
	}

	//Parse the input probe summary file
	Spacecraft probe;
	std::cout << "Reading probe summary file..." << std::endl;
	probe = read_probe_summary_file(probe_summary_file_name);
	std::cout << "Done reading probe summary file." << std::endl << std::endl;

	std::string GTOC7_ephemeris_file_name = "GTOC7ASTEROIDS.TXT";

	std::vector <Body> asteroids;
	Body Earth;

	std::cout << "Reading ephemeris..." << std::endl;

	//scan in GTOC7 asteroid ephemeris data
	asteroids = read_ephemeris(GTOC7_ephemeris_file_name);

	//Earth's ephemeris data
	Earth.a = 0.999988049532578 * 149597870.691;
	Earth.ecc = 0.01671681163160;
	Earth.inc = 0.0008854353079654 * PI / 180.0;
	Earth.omega = 287.61577546182 * PI / 180.0;
	Earth.LAN = 175.40647696473 * PI / 180.0;
	Earth.M = 257.60683707535 * PI / 180.0;
	Earth.time_stamp = 54000.0 * 86400.0;

	Earth.name = "Earth";
	Earth.GTOC7_num = 0;

	Earth.E = laguerreConway(Earth.ecc, Earth.M);
	Earth.tru = 2.0*atan(sqrt((1.0 + Earth.ecc) / (1.0 - Earth.ecc))*tan(Earth.E / 2.0));
	coe2cartesian(&Earth);

	std::cout << "Done reading ephemeris." << std::endl << std::endl;

	//Propagate all bodies to Jan. 1st 2021 (MJD 59215)
	double delta_t = (59215.0 - 56800.0)*86400.0;

	for (size_t i = 0; i < asteroids.size(); ++i)
		orbitprop(&asteroids[i], delta_t);

	delta_t = (59215.0 - 54000.0)*86400.0;
	orbitprop(&Earth, delta_t);

	//Start with the first phase
	int phase = 0;

	//if we are the probe that stays at the first asteroid, then move
	//on to check the first "real" phase
	if (probe.x[0].size() == 1)
		phase = 1;


	Spacecraft myprobe = probe; //make a copy of the probe

	double h;
	double error = 1.0e+20;
	int ns = 7;
	std::vector <double> X_left, X_right;
	std::vector <double> dX(7, 0.0);
	std::vector <double> Tvec(3, 0.0);


	double mu_sun = 132712440018.0;
	bool normalized_integrator = true;
	bool adaptive_step = true;
	double precisionTarget;
	double DU, TU;


	if (normalized_integrator)
	{
		DU = 149597870.691;
		TU = sqrt(DU * DU * DU / mu_sun);
		mu_sun = 1.0;
		precisionTarget = 1.0e-13;
	}
	else
	{
		DU = 1.0;
		TU = 1.0;
		precisionTarget = 1.0e-7;
	}

	if (adaptive_step)
	{

		//for phases -- add this loop in to check all phases once we are confident in the code for a single phase

		X_left = { myprobe.x[phase][0] / DU, myprobe.y[phase][0] / DU, myprobe.z[phase][0] / DU, myprobe.vx[phase][0] * TU / DU, myprobe.vy[phase][0] * TU / DU, myprobe.vz[phase][0] * TU / DU, myprobe.mass[phase][0] };
		//X_left = { 3.82858625335622E+08, - 3.05100026716626E+05,  2.89100959348077E+06, - 2.60587333616076E-01 , 1.97262246721978E+01,  8.74345055541848E-01,  1.85287741774152E+03 };
		//X_left = { Earth.x + 90.0, Earth.y + 90.0, Earth.z + 90.0, Earth.vx + 0.000009, Earth.vy + 0.000009, Earth.vz + 0.000009, 1.0 };


		for (size_t timestep = 0; timestep < 1; ++timestep)
		{
			h = 86400.0 / TU;
			//convert from N to kN
			Tvec[0] = probe.Tx[phase][timestep] / 1000.0 * TU * TU / DU;
			Tvec[1] = probe.Ty[phase][timestep] / 1000.0 * TU * TU / DU;
			Tvec[2] = probe.Tz[phase][timestep] / 1000.0 * TU * TU / DU;
			//Tvec[0] = 0.0; Tvec[1] = 0.0; Tvec[2] = 0.0;
			//num_sub_steps = 86400;

			X_right = adaptive_step_int(X_left, Tvec, h, ns, precisionTarget, DU, TU, mu_sun);
			X_left = X_right;
		}

		std::cout << std::setprecision(14) << ' ' << X_left[0] * DU << ' ' << X_left[1] * DU << ' ' << X_left[2] * DU << ' ' << X_left[3] * DU / TU << ' ' << X_left[4] * DU / TU << ' ' << X_left[5] * DU / TU << ' ' << X_left[6] << std::endl;

		getchar();

	}


	else
	{
		h = 0.01;
		int num_sub_steps;
		//Analytical propagation of Earth
		//std::cout << std::setprecision(16) << "Initial " << Earth.a / 149597870.691 << ' ' << Earth.ecc << ' ' << Earth.inc*180.0 / PI << ' ' << Earth.omega*180.0 / PI << ' ' << Earth.LAN*180.0 / PI << ' ' << Earth.M*180.0 / PI << std::endl;
		//delta_t = 365.0*86400.0;
		//orbitprop(&Earth, delta_t);

		//std::cout << std::setprecision(16) << "Analytical " << Earth.a / 149597870.691 << ' ' << Earth.ecc << ' ' << Earth.inc*180.0 / PI << ' ' << Earth.omega*180.0 / PI << ' ' << Earth.LAN*180.0 / PI << ' ' << Earth.M*180.0 / PI << std::endl;

		//for (size_t timestep = probe.x[phase].size() - 2; timestep < probe.x[phase].size() - 1; ++timestep)
		//for (size_t timestep = 0; timestep < 365; ++timestep)
		for (size_t timestep = 0; timestep < 10; ++timestep)
		{
			num_sub_steps = int(probe.time_stamp[phase][timestep + 1] - probe.time_stamp[phase][timestep]) / h;

			//convert from N to kN
			Tvec[0] = probe.Tx[phase][timestep] / 1000.0 * TU * TU / DU;
			Tvec[1] = probe.Ty[phase][timestep] / 1000.0 * TU * TU / DU;
			Tvec[2] = probe.Tz[phase][timestep] / 1000.0 * TU * TU / DU;
			//Tvec[0] = 0.0; Tvec[1] = 0.0; Tvec[2] = 0.0;
			//num_sub_steps = 86400;

			for (size_t sub_step = 0; sub_step < num_sub_steps; ++sub_step)
			{
				dX = GTOC7EOM(X_left, Tvec, DU, TU, mu_sun);
				X_right = rk8713M(X_left, Tvec, dX, h, ns, error, DU, TU, mu_sun);
				X_left = X_right;

				//std::cout << std::setprecision(10) << sub_step << ' ' << X_left[0] << ' ' << X_left[1] << ' ' << X_left[2] << ' ' << X_left[3] << ' ' << X_left[4] << ' ' << X_left[5] << ' ' << X_left[6] << std::endl;
			}
			//std::cout << std::setprecision(10) << timestep << ' ' << X_left[0] << ' ' << X_left[1] << ' ' << X_left[2] << ' ' << X_left[3] << ' ' << X_left[4] << ' ' << X_left[5] << ' ' << X_left[6] << std::endl;
			//std::cout << timestep << std::endl;
		}
		std::cout << std::setprecision(14) << ' ' << X_left[0] << ' ' << X_left[1] << ' ' << X_left[2] << ' ' << X_left[3] << ' ' << X_left[4] << ' ' << X_left[5] << ' ' << X_left[6] << std::endl;
		//Earth.x = X_left[0]; Earth.y = X_left[1]; Earth.z = X_left[2]; Earth.vx = X_left[3]; Earth.vy = X_left[4]; Earth.vz = X_left[5];

		//cartesian2coe(&Earth);

		//std::cout << std::setprecision(16) << "Integrated " << Earth.a / 149597870.691 << ' ' << Earth.ecc << ' ' << Earth.inc*180.0 / PI << ' ' << Earth.omega*180.0 / PI << ' ' << Earth.LAN*180.0/PI << ' ' << Earth.M*180.0/PI << std::endl;

		//std::cout << std::setprecision(10) << X_left[0] << ' ' << X_left[1] << ' ' << X_left[2] << ' ' << X_left[3] << ' ' << X_left[4] << ' ' << X_left[5] << ' ' << X_left[6] << std::endl;
		getchar();
		//advance probe to next departure
	}



	return 0;
}