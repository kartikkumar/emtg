//EMTG body class

#include "missionoptions.h"
#include "body.h"
#include "EMTG_math.h"
#include "Astrodynamics.h"
#include "kepler_lagrange_laguerre_conway.h"

#include "SpiceUsr.h"

#include <string>
#include <vector>
#include <math.h>
#include <fstream>

namespace EMTG {namespace Astrodynamics {

	//default constructor doesn't do anything
	body::body(){}

	//constructor takes in data and calls load_body_data()
	body::body(const int& ibody_code, const string& iname, const string& ishortname,  const int& ispice_ID, const double& imininum_altitude, const double& imass, const double& iradius, const double& iepoch, vector<double>& ireference_angles, vector<double>& iclassical_orbit_elements, const double& iuniverse_mu, const int& icentral_body_SPICE_ID, const string& icentral_body_name, const double& icentral_body_radius, missionoptions* options)
	{
		load_body_data(ibody_code, iname, ishortname, ispice_ID, imininum_altitude, imass, iradius, iepoch, ireference_angles, iclassical_orbit_elements, iuniverse_mu, icentral_body_SPICE_ID, icentral_body_name, icentral_body_radius, options);
	}

	//destructor
	body::~body(){}

	//function to load new data into the body
	void body::load_body_data(const int& ibody_code, const string& iname, const string& ishortname, const int& ispice_ID, const double& imininum_altitude, const double& imass, const double& iradius, const double& iepoch, vector<double>& ireference_angles, vector<double>& iclassical_orbit_elements, const double& iuniverse_mu, const int& icentral_body_SPICE_ID, const string& icentral_body_name, const double& icentral_body_radius, missionoptions* options)
	{
		//copy information from the inputs into the body
		name = iname;
		short_name = ishortname;
		universe_mu = iuniverse_mu;
		body_code = ibody_code;
		central_body_spice_ID = icentral_body_SPICE_ID;
		central_body_name = icentral_body_name;
		central_body_radius = icentral_body_radius;

		spice_ID = ispice_ID;
		minimum_safe_flyby_altitude = imininum_altitude;
		mass = imass;
		radius = iradius;
		reference_epoch = iepoch;
		SMA = iclassical_orbit_elements[0];
		ECC = iclassical_orbit_elements[1];
		INC = iclassical_orbit_elements[2] * EMTG::math::PI / 180.0;
		RAAN = iclassical_orbit_elements[3] * EMTG::math::PI / 180.0;
		AOP = iclassical_orbit_elements[4] * EMTG::math::PI / 180.0;
		MA = iclassical_orbit_elements[5] * EMTG::math::PI / 180.0;
		


		//compute additional values
		mu = options->G * mass;
		if (ECC < 0.2)
			r_SOI = SMA * pow(mu / universe_mu, 0.4);
		else
			r_SOI = SMA * (1 - ECC) * pow(mu / (3.0 * universe_mu), 0.333333333333333333333333);

		//determine which ephemeris to draw from
		if (options->ephemeris_source == 0)
		{
			body_ephemeris_source = 0; //use static ephemeris
			ephemeris_start_date = -0;
			ephemeris_end_date = 1e+10;
		}
		else if (options->ephemeris_source == 1)
		{
			//first, check to see if the body exists in the currently loaded SPICE kernels
			double temp_state[6];
			double LT_dump;
			spkez_c (spice_ID, unitim_c(reference_epoch + 2400000.5, "JDTDB", "ET"), "J2000", "NONE", central_body_spice_ID, temp_state, &LT_dump);
			if (failed_c())
				reset_c();


			if (fabs(temp_state[0]) > 1.0e-6 && fabs(temp_state[0]) < 1.0e+50)
			{
				body_ephemeris_source = 1; //body can be located using SPICE
			}
			else
			{
				cout << "Warning, body " << name << " does not have a SPICE ephemeris file." << endl;
				body_ephemeris_source = 0; //use static ephemeris
				ephemeris_start_date = 0;
				ephemeris_end_date = 1e+10;
			}
		}

		J2000_body_equatorial_frame.initialize(ireference_angles[0], ireference_angles[1], ireference_angles[2], ireference_angles[3], ireference_angles[4], ireference_angles[5]);
	}

	//function to find the body state vector at epoch
	int body::locate_body(const double& epoch, double* state, const bool& need_deriv, missionoptions* options)
	{
		double DT, n, M, E, V[6];

		switch (body_ephemeris_source)
		{
			case 1: //SPICE
				double LT_dump;
				spkez_c (spice_ID, unitim_c(epoch + 2400000.5, "JDTDB", "ET"), "J2000", "NONE", central_body_spice_ID, state, &LT_dump);

				if (need_deriv)
				{
					double statepert[6];
					spkez_c (spice_ID, unitim_c(epoch + 2400000.5, "JDTDB", "ET") + (1.0 / 86400.0), "J2000", "NONE", central_body_spice_ID, statepert, &LT_dump);
					state[6] = (statepert[3] - state[3]) / (1.0 / 86400.0);
					state[7] = (statepert[4] - state[4]) / (1.0 / 86400.0);
					state[8] = (statepert[5] - state[5]) / (1.0 / 86400.0);
				}

				break;
			case 0: //static ephemeris
					//TODO static ephemeris is not ready!
					//note, always should give in Earth equatorial J2000 coordinates for internal processing

					DT = ( epoch - reference_epoch ) * 86400;
					

					if (SMA > 0.0)
						n = sqrt( universe_mu / (SMA*SMA*SMA) );
					else
						n = sqrt( universe_mu / (-SMA*SMA*SMA) );
					
					M = MA + n*DT;
					M = fmod(M, 2 * EMTG::math::PI);

					E = Kepler::KeplerLaguerreConway(ECC, M);
					V[0] = SMA; 
					V[1] = ECC; 
					V[2] = INC;
					V[3] = RAAN;
					V[4] = AOP; 
					V[5] = E;

					COE2inertial(V, universe_mu, state);

					if (need_deriv)
					{
						double r = sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
						double r3 = r*r*r;
						state[6] = -universe_mu/r3 * state[0];
						state[7] = -universe_mu/r3 * state[1];
						state[8] = -universe_mu/r3 * state[2];
					}
				break;
			default:
				cout << "Invalid ephemeris source " << body_ephemeris_source << " for object " << name << endl;
				cout << "Program halted. Press enter to quit." << endl;
#ifndef BACKGROUND_MODE
				cin.ignore();
#endif
		}

		return 0;
	}
	
	//function to locate a point on the sphere of influence in cartesian coordinates (Earth Equatorial J2000, measured from central body of current universe)
	int body::locate_point_on_SOI(const double& theta, const double& phi, double* point_relative_to_body)
	{
		point_relative_to_body[0] = r_SOI * cos(theta)*cos(phi);
		point_relative_to_body[1] = r_SOI * sin(theta)*cos(phi);
		point_relative_to_body[2] = r_SOI * sin(phi);

		return 0;
	}

	//function to print body to screen (for debug purposes)
	void body::print_body_to_screen(string filename)
	{
		ofstream outputfile(filename.c_str(), ios::app);
		outputfile << "Body name: " << name << endl;
		outputfile << "Short name: " << short_name << endl;
		outputfile << "Body position in menu: " << body_code << endl;
		outputfile << "SPICE ID: " << spice_ID << endl;
		outputfile << "Valid flyby target? " << (minimum_safe_flyby_altitude > 0.0 ? "True" : "False") << endl;
		if (minimum_safe_flyby_altitude > 0.0)
			outputfile << "Minimum safe flyby altitude (km) " << minimum_safe_flyby_altitude << endl;
		outputfile << "Mass (kg): " << mass << endl;
		outputfile << "Radius (km): " << radius << endl;
		outputfile << "Ephemeris source: " << body_ephemeris_source << endl;
		outputfile << "R_SOI: " << r_SOI << endl;
		outputfile << "Reference Epoch (MJD): " << reference_epoch << endl;
		outputfile << "SMA (km): " << SMA << endl;
		outputfile << "ECC: " << ECC << endl;
		outputfile << "INC (deg): " << INC * 180.0 / EMTG::math::PI << endl;
		outputfile << "RAAN (deg): " << RAAN * 180.0 / EMTG::math::PI << endl;
		outputfile << "AOP (deg): " << AOP * 180.0 / EMTG::math::PI << endl;
		outputfile << "MA (deg): " << MA * 180.0 / EMTG::math::PI << endl;
		outputfile << endl;

		outputfile.close();
	}
	
	//comparator
	bool body::operator== (const body& OtherBody)
	{
		//compare three fields for accuracy
		if (this->name == OtherBody.name && this->spice_ID == OtherBody.spice_ID && this->mass == OtherBody.mass)
		{
			return true;
		}
		return false;
	}

	bool body::operator!= (const body& OtherBody)
	{
		return !(*this == OtherBody);
	}
	
}}//close namespace