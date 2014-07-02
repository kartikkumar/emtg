/* 
 * GMATScripter.cpp
 * 
 * Creation Date: july 1st 2014
 * Author: ryne beeson
 *
 * Purpose: methods for the GMATScripter Class, which generates a GMAT script file 
 *			following an EMTG run. Methods build on work conducted by Max Schadegg (summer 2013)
 *			and Jacob Englander
 */


#include "GMATScripter.h"
#include "mission.h"

//#include "boost/algorithm/string.hpp"
//#include "boost/algorithm/string/predicate.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"

#include "SpiceUsr.h"

#include <string>  //string
#include <fstream> //ofstream
//#include <sstream> //stringstream


namespace EMTG {

// default constructor; not intended for use
gmatscripter::gmatscripter(){
}

// constructor
gmatscripter::gmatscripter(mission* mission_in){
	this->ptr_gmatmission = mission_in;
}

// destructor
gmatscripter::~gmatscripter(){
}


// method to write out the GMAT script
void gmatscripter::write_GMAT_script(){

	// open a file for writing
	this->create_file();
	
	// collect the bodies used during the GMAT mission
	this->parse_bodies_list();

	// write out the preamble
	this->create_preamble();

	// write out the spacecraft
	this->write_GMAT_spacecraft();

	// write out the spacecraft hardware (i.e. thrusters and fuel tanks)
	this->write_GMAT_hardware();
	
	// write out the forcemodels
	this->write_GMAT_forcemodels();

	// write out the propagators
	this->write_GMAT_propagators();

	// write out the burn objects
	this->write_GMAT_burns();

	// write out arrays and variables
	this->write_GMAT_variables();
		
	// write out coordinate systems
	this->write_GMAT_coordinatesystems();

	// write out the solvers
	this->write_GMAT_solvers();

	// write out subscribers (i.e. plot views and reports)
	this->write_GMAT_subscribers();

}


// method to create a GMAT script file
void gmatscripter::create_file(){

	//create a filename
	string filename = this->ptr_gmatmission->options.working_directory + "//" + 
					  this->ptr_gmatmission->options.mission_name + "_" + 
					  this->ptr_gmatmission->options.description + "_GMAT_Test.script";
	//open the ofstream object called GMATfile
	this->GMATfile.open(filename.c_str(), ios::trunc);
	//set floating point decimal precision
	this->GMATfile.precision(25);

}


// method to parse the bodies used in the mission
void gmatscripter::parse_bodies_list(){

	//add bodies of first journey to the mission body list
	//TODO:: add other arrival/departure types
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		if (j == 0)
		{
			for (int p = 0; p < this->ptr_gmatmission->options.number_of_phases[0] + 1; ++p)
				missionbodies.push_back(this->ptr_gmatmission->TheUniverse[0].bodies[this->ptr_gmatmission->options.sequence[0][p] - 1]);
		}
		//add all other bodies to mission body list
		else
		{
			for (int p = 1; p < this->ptr_gmatmission->options.number_of_phases[j] + 1; ++p)
				missionbodies.push_back(this->ptr_gmatmission->TheUniverse[j].bodies[this->ptr_gmatmission->options.sequence[j][p] - 1]);
		}
	}

	int body_flag = 0;
	//remove duplicate entries in visited bodies list
	for (int index_body_visited = 0; index_body_visited < missionbodies.size(); ++index_body_visited)
	{
		body_flag = 0;
		for (int index_body_visited_unique = 0; index_body_visited_unique < missionbodies_unique.size(); ++index_body_visited_unique)
		{
			if (missionbodies[index_body_visited].spice_ID == missionbodies_unique[index_body_visited_unique].spice_ID)
			{
				body_flag = 1;
			}
		}
		//if body flag not switch 'on', then the body is unique; add it to the missionbodies_unique vector
		if (body_flag == 0)
		{
			missionbodies_unique.push_back(missionbodies[index_body_visited]);
		}
	}

	//make sure that each body doesnt start with a number for GMAT's sake
	for (int index_body_visited = 0; index_body_visited < missionbodies.size(); ++index_body_visited)
	{
		if ((missionbodies[index_body_visited].name[0] == 0) || (missionbodies[index_body_visited].name[0] == 1) || (missionbodies[index_body_visited].name[0] == 2) || (missionbodies[index_body_visited].name[0] == 3) || (missionbodies[index_body_visited].name[0] == 4) || (missionbodies[index_body_visited].name[0] == 5) || (missionbodies[index_body_visited].name[0] == 6) || (missionbodies[index_body_visited].name[0] == 7) || (missionbodies[index_body_visited].name[0] == 8) || (missionbodies[index_body_visited].name[0] == 9))
		{
			missionbodies[index_body_visited].name = "A" + missionbodies[index_body_visited].name;
		}
	}
	for (int index_body_visited = 0; index_body_visited < missionbodies_unique.size(); ++index_body_visited)
	{
		if ((missionbodies_unique[index_body_visited].name[0] == 0) || (missionbodies_unique[index_body_visited].name[0] == 1) || (missionbodies_unique[index_body_visited].name[0] == 2) || (missionbodies_unique[index_body_visited].name[0] == 3) || (missionbodies_unique[index_body_visited].name[0] == 4) || (missionbodies_unique[index_body_visited].name[0] == 5) || (missionbodies_unique[index_body_visited].name[0] == 6) || (missionbodies_unique[index_body_visited].name[0] == 7) || (missionbodies_unique[index_body_visited].name[0] == 8) || (missionbodies_unique[index_body_visited].name[0] == 9))
		{
			missionbodies_unique[index_body_visited].name = "A" + missionbodies_unique[index_body_visited].name;
		}
	}

}//end of parse_bodies_list() method


// method to create the script preamble
void gmatscripter::create_preamble(){

	//get the current timestamp from boost and assign it to now
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::stringstream timestream;
	timestream << static_cast<int>(now.date().month()) << "/" << now.date().day() << "/" << now.date().year() << " " << now.time_of_day().hours() << ":" << now.time_of_day().minutes() << ":" << now.time_of_day().seconds();

	//preamble
	GMATfile << "%--------------------------------------------------------------------------------" << endl;
	GMATfile << "%GMAT script created by EMTGv8" << endl;
	GMATfile << "%EMTG options file: " << this->ptr_gmatmission->options.working_directory + "/" + 
										  this->ptr_gmatmission->options.mission_name + ".emtgopt" << endl;
	GMATfile << "%EMTG output file: " <<  this->ptr_gmatmission->options.outputfile << endl;
	GMATfile << "%EMTG output written on: " << timestream.str() << endl;
	GMATfile << "%Author(s): Ryne Beeson     (2014_07_01)" << endl;
	GMATfile << "%           Max Schadegg    (2013_06_03)" << endl;
	GMATfile << "%           Jacob Englander (2013_08_09)" << endl;
	GMATfile << "%--------------------------------------------------------------------------------" << endl;
	GMATfile << endl;
	GMATfile << endl;

}


// method to create spacecraft information
void gmatscripter::write_GMAT_spacecraft(){

	//declarations
	int body_index = 0;
	double epoch = 0.0;

	//spacecraft header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Spacecraft" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//create and write out two s/c for each phase in each journey (forward + backward)
	//for-statement for counting through every journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for-statement for counting through each phase in a journey
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//create a name for the forward propagating spacecraft
			std::ostringstream spacecraft_name;
			spacecraft_name << "SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward";
			//add the name to the vector 'spacecraft_names'
			spacecraft_names.push_back(spacecraft_name.str());

			epoch = (this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch / 86400.0) + 2400000.5 - 2430000;
			//write out the spacecraft information
			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << ", forward propagated s/c" << endl;
			GMATfile << "Create Spacecraft " << spacecraft_name.str() << ";" << endl;
			GMATfile << spacecraft_name.str() << ".DateFormat = TAIModJulian;" << endl;
			GMATfile << spacecraft_name.str() << ".Epoch = " << epoch << ";" << endl;
			GMATfile << spacecraft_name.str() << ".DryMass = 0" << endl;
			GMATfile << spacecraft_name.str() << ".CoordinateSystem = " << missionbodies[body_index].name << "J2000Eq;" << endl;
			GMATfile << spacecraft_name.str() << ".Tanks = {FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward};" << endl;
			GMATfile << spacecraft_name.str() << ".Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward};" << endl;
			GMATfile << endl;

			//create a name for the backward propagating spacecraft
			spacecraft_name.str("");
			spacecraft_name << "SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward";
			spacecraft_names.push_back(spacecraft_name.str());

			epoch = (this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch + 
					 this->ptr_gmatmission->journeys[j].phases[p].TOF) / 86400.0 + 2400000.5 - 2430000;
			//write out the spacecraft information
			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << ", backward propagated s/c" << endl;
			GMATfile << "Create Spacecraft " << spacecraft_name.str() << ";" << endl;
			GMATfile << spacecraft_name.str() << ".DateFormat = TAIModJulian;" << endl;
			GMATfile << spacecraft_name.str() << ".Epoch = " << epoch << ";" << endl;
			GMATfile << spacecraft_name.str() << ".DryMass = 0" << endl;
			GMATfile << spacecraft_name.str() << ".CoordinateSystem = " << missionbodies[body_index + 1].name << "J2000Eq;" << endl;
			GMATfile << spacecraft_name.str() << ".Tanks = {FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward};" << endl;
			GMATfile << spacecraft_name.str() << ".Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward};" << endl;

			//increment
			++body_index;

			//spacing
			GMATfile << endl;
			GMATfile << endl;

		}// end of phases for-statement
	}//end of journeys for-statement
}//end of write_GMAT_spacecraft() method


// method to create hardware information
void gmatscripter::write_GMAT_hardware(){

	//declarations
	int body_index = 0;

	//hardware header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Hardware components" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//create two tanks and two thrusters for each phase in each journey (forward + backward)
			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << ", forward s/c tanks and thrusters" << endl;
			GMATfile << "Create FuelTank FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
			GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.AllowNegativeFuelMass = false;" << endl;
			GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.Volume = 10;" << endl;
			GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.FuelMass = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[6] << endl;
			GMATfile << endl;

			GMATfile << "Create Thruster Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.CoordinateSystem = " << missionbodies[body_index].name << "J2000Eq" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1 = 1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2 = 0;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3 = 0;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.DutyCycle = 1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.Tank = FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustScaleFactor = 1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.DecrementMass = true;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1 = .1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.K1 = 3000;" << endl;
			GMATfile << endl;
			GMATfile << endl;

			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << ", backward s/c tanks and thrusters" << endl;
			GMATfile << "Create FuelTank FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
			GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.AllowNegativeFuelMass = false;" << endl;
			GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.Volume = 10;" << endl;
			GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[6] << endl;

			GMATfile << endl;
			GMATfile << "Create Thruster Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.CoordinateSystem = " << missionbodies[body_index + 1].name << "J2000Eq" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection1 = 1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection2 = 0;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection3 = 0;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.DutyCycle = 1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.Tank = FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustScaleFactor = 1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.DecrementMass = true;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.C1 = .1;" << endl;
			GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.K1 = 3000;" << endl;
			GMATfile << endl;
			GMATfile << endl;
			++body_index;
		}//end of journeys for-statement
	}//end of phases for-statement
	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;
}//end of write_GMAT_hardware() method


// method to create forcemodel information
void gmatscripter::write_GMAT_forcemodels(){

	//declarations
	string filestring;
	vector <fs::path> SPICE_files;
	SPICEDOUBLE_CELL(spice_coverage, 10000);
	SpiceInt number_of_windows = 0;

	//force model header 
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Force Models" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//create central force model (TODO:: assume for now that first body's central body is the principle central body)
	GMATfile << "%Create 2-body force model for central body" << endl;
	GMATfile << "Create ForceModel " << missionbodies_unique[0].central_body_name << "2Bod;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "2Bod.CentralBody = " << missionbodies_unique[0].central_body_name << ";" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "2Bod.PointMasses = {" << missionbodies_unique[0].central_body_name << "};" << endl; //, Mercury, Venus, Earth, Luna, Mars, Jupiter, Saturn, Uranus, Neptune};" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "2Bod.Drag = None;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "2Bod.SRP = Off;" << endl;
	GMATfile << endl;

	GMATfile << "%Create 2-body force models for all bodies visited" << endl;

	//for each body in mission
	for (int index_body_visited = 0; index_body_visited < missionbodies_unique.size(); ++index_body_visited)
	{
		//must create any bodies that are not already defined in GMAT
		if ((missionbodies[index_body_visited].name != "Sun") && (missionbodies[index_body_visited].name != "Mercury") && (missionbodies[index_body_visited].name != "Venus") && (missionbodies[index_body_visited].name != "Earth") && (missionbodies[index_body_visited].name != "Mars") && (missionbodies[index_body_visited].name != "Jupiter") && (missionbodies[index_body_visited].name != "Saturn") && (missionbodies[index_body_visited].name != "Uranus") && (missionbodies[index_body_visited].name != "Neptune") && (missionbodies[index_body_visited].name != "Pluto"))
		{
			GMATfile << "%Must create model for body visited" << endl;
			GMATfile << "Create Planet " << missionbodies_unique[index_body_visited].name << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".NAIFId = " << missionbodies_unique[index_body_visited].spice_ID << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".EquatorialRadius = " << missionbodies_unique[index_body_visited].radius << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".Mu = " << missionbodies_unique[index_body_visited].mu << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".PosVelSource = 'SPICE'" << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".CentralBody = '" << missionbodies_unique[index_body_visited].central_body_name << "'" << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".RotationDataSource = 'IAUSimplified'" << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".OrientationEpoch = 21545" << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".SpinAxisRAConstant = " << missionbodies_unique[index_body_visited].J2000_body_equatorial_frame.alpha0 << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".SpinAxisRARate = " << missionbodies_unique[index_body_visited].J2000_body_equatorial_frame.alphadot << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".SpinAxisDECConstant = " << missionbodies_unique[index_body_visited].J2000_body_equatorial_frame.delta0 << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".SpinAxisDECRate = " << missionbodies_unique[index_body_visited].J2000_body_equatorial_frame.deltadot << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".RotationConstant = " << missionbodies_unique[index_body_visited].J2000_body_equatorial_frame.W << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".RotationRate = " << missionbodies_unique[index_body_visited].J2000_body_equatorial_frame.Wdot << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << ".OrbitSpiceKernelName = {";

			//find which spice file the body is located
			EMTG::filesystem::get_all_files_with_extension(fs::path(this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

			for (size_t k = 0; k < SPICE_files.size(); ++k)
			{
				filestring = this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

				//check if body is located in spice file
				scard_c(0, &spice_coverage);
				spkcov_c(filestring.c_str(), missionbodies_unique[index_body_visited].spice_ID, &spice_coverage);
				number_of_windows = wncard_c(&spice_coverage);
				if (number_of_windows > 0)
				{
					GMATfile << "'" << filestring << "', ";
				}
				number_of_windows = 0;
			}
			GMATfile << "'" << filestring << "'};" << endl;
			GMATfile << endl;
		}

		//must create any central bodies that are not already defined in GMAT
		if ((missionbodies_unique[index_body_visited].central_body_name != "Sun") && (missionbodies_unique[index_body_visited].central_body_name != "Mercury") && (missionbodies_unique[index_body_visited].central_body_name != "Venus") && (missionbodies_unique[index_body_visited].central_body_name != "Earth") && (missionbodies_unique[index_body_visited].central_body_name != "Mars") && (missionbodies_unique[index_body_visited].central_body_name != "Jupiter") && (missionbodies_unique[index_body_visited].central_body_name != "Saturn") && (missionbodies_unique[index_body_visited].central_body_name != "Uranus") && (missionbodies_unique[index_body_visited].central_body_name != "Neptune") && (missionbodies_unique[index_body_visited].central_body_name != "Pluto"))
		{
			GMATfile << "%Must create model for central body" << endl;
			GMATfile << "Create Planet " << missionbodies_unique[index_body_visited].central_body_name << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".NAIFId = " << missionbodies_unique[index_body_visited].central_body_spice_ID << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".EquatorialRadius = " << missionbodies_unique[index_body_visited].central_body_radius << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".Mu = " << missionbodies_unique[index_body_visited].universe_mu << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".PosVelSource = 'SPICE'" << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".CentralBody = 'Sun'" << endl; //assume Sun for now
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".RotationDataSource = 'IAUSimplified'" << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".OrientationEpoch = 21545" << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisRAConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alpha0 << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisRARate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alphadot << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisDECConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.delta0 << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisDECRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.deltadot << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".RotationConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.W << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".RotationRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.Wdot << endl;
			GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".OrbitSpiceKernelName = {";

			//find which spice file the body is located
			EMTG::filesystem::get_all_files_with_extension(fs::path(this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

			for (size_t k = 0; k < SPICE_files.size(); ++k)
			{
				filestring = this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

				//check if body is located in spice file
				scard_c(0, &spice_coverage);
				spkcov_c(filestring.c_str(), missionbodies_unique[index_body_visited].central_body_spice_ID, &spice_coverage);
				number_of_windows = wncard_c(&spice_coverage);
				if (number_of_windows > 0)
				{
					GMATfile << "'" << filestring << "', ";
				}
				number_of_windows = 0;
			}
			GMATfile << "'" << filestring << "'};" << endl;
			GMATfile << endl;
		}

		//create 2body force model for each body visited
		GMATfile << "Create ForceModel " << missionbodies_unique[index_body_visited].name << "2Bod;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "2Bod.CentralBody = " << missionbodies_unique[index_body_visited].name << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "2Bod.PointMasses = {" << missionbodies_unique[index_body_visited].name << "};" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "2Bod.Drag = None;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "2Bod.SRP = Off;" << endl;
		GMATfile << endl;
	}
	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;
}//end of write_GMAT_forcemodels() method


// method to create propagator information
void gmatscripter::write_GMAT_propagators(){

	//propagator header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Propagators" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//far planet propagator, larger time steps
	//TODO:: for now it is assumed that the principle central body is the first central body
	GMATfile << "%Create propagation model for central body" << endl;
	GMATfile << "Create Propagator " << missionbodies_unique[0].central_body_name << "Prop;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "Prop.FM = " << missionbodies_unique[0].central_body_name << "2Bod;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "Prop.Type                     = PrinceDormand78;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "Prop.InitialStepSize          = 60;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "Prop.Accuracy                 = 1e-11;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "Prop.MinStep                  = 0.0;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "Prop.MaxStep                  = 86400;" << endl;
	GMATfile << endl;

	//create propagators for close body approaches, smaller time steps
	GMATfile << "%Create propagation models for other bodies" << endl;
	for (int index_body_visited = 0; index_body_visited < missionbodies_unique.size(); ++index_body_visited)
	{
		GMATfile << "Create Propagator " << missionbodies_unique[index_body_visited].name << "Prop;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "Prop.FM = " << missionbodies_unique[index_body_visited].name << "2Bod;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "Prop.Type                     = PrinceDormand78;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "Prop.InitialStepSize          = 30;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "Prop.Accuracy                 = 1e-11;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "Prop.MinStep                  = 0.0;" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "Prop.MaxStep                  = 8640;" << endl;
		GMATfile << endl;
	}
	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_propagators() method


// method to create burn information
void gmatscripter::write_GMAT_burns(){

	//burn header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Burns" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//impulsive or finitie burn
	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//create two finite burn objects for each phase in each journey (forward + backward)
			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << ", forward and backward finite burn objects" << endl;
			GMATfile << "Create FiniteBurn FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
			GMATfile << "FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward.Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward};" << endl;
			GMATfile << "Create FiniteBurn FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
			GMATfile << "FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward.Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward};" << endl;
			GMATfile << endl;
		}//end of journeys for-statement
	}//end of phases for-statement
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_burns() method


// method to create array and variable information
void gmatscripter::write_GMAT_variables(){

	//arrays and variables header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Arrays, Variables, Strings" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//create mission variables
	GMATfile << "Create Variable ObjectiveFunction FinalEpoch LaunchEpoch_Scaled" << endl;
	//DEBUG: Check that these are the only impulsive-thrust phase types
	if (this->ptr_gmatmission->options.mission_type < 2 || this->ptr_gmatmission->options.mission_type == 5) //impulsive-thrust phase types
	{
		GMATfile << "Create Variable FinalMass_Scaled ThrusterISP" << endl;
	}
	else //low-thrust phase types
	{
		GMATfile << "Create Variable FinalMass_Scaled ThrusterISP ThrusterDutyCycle ThrusterMaxThrust" << endl;
	}


	//create variable for each s/c create rdotv temp variable
	for (int index_SC = 0; index_SC < spacecraft_names.size(); ++index_SC)
	{
		GMATfile << "Create Variable " << spacecraft_names[index_SC] << "_RdotV_Scaled " << spacecraft_names[index_SC] << "_PeriapseRadius_Scaled" << endl;
	}
	GMATfile << endl;

	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{

			//create the time and state variables
			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << " time variables" << endl;
			//create journey wait time variable
			if (p == 0)
			{
				GMATfile << "Create Variable Journey" << j + 1 << "_WaitTime_Scaled" << endl;
			}
			//create time of flight variables
			GMATfile << "Create Variable TOF_Journey" << j + 1 << "Phase" << p + 1 << "Scaled" << endl;
			GMATfile << "Create Variable TOF_Journey" << j + 1 << "Phase" << p + 1 << endl;
			GMATfile << "Create Variable TimeStepLength_Journey" << j + 1 << "Phase" << p + 1 << endl;
			GMATfile << "Create Variable LeftBoundarySOITime_Journey" << j + 1 << "Phase" << p + 1 << endl;
			GMATfile << "Create Variable RightBoundarySOITime_Journey" << j + 1 << "Phase" << p + 1 << endl;

			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << " forward variables" << endl;
			//created scaled state variables
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_X_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Y_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Z_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_VX_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_VY_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_VZ_Scaled" << endl;
			GMATfile << "Create Variable FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled" << endl;
			GMATfile << endl;
			GMATfile << "% Journey #" << j + 1 << ", Phase #" << p + 1 << " backward variables" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_X_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Y_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Z_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_VX_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_VY_Scaled" << endl;
			GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_VZ_Scaled" << endl;
			GMATfile << "Create Variable FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled" << endl;
			GMATfile << endl;

			//create the control variables
			//create arrays to hold thrust direction vector and unit vectors for each timestep
			GMATfile << "Create Array ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "[" << this->ptr_gmatmission->options.num_timesteps << ", 3];" << endl;
			GMATfile << "Create Array ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "[" << this->ptr_gmatmission->options.num_timesteps << ", 1];" << endl;

			GMATfile << endl;
		}
	}
	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_variables() method


// method to create coordinate systems information
void gmatscripter::write_GMAT_coordinatesystems(){

	//coordinate systems header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Coordinate systems" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//always in J2000 Eq. Syst.
	GMATfile << "%Create coordinate systems for plotting/viewing" << endl;
	//TODO:: as of now it is assumed that the first central body is the principle central body
	GMATfile << "Create CoordinateSystem " << missionbodies_unique[0].central_body_name << "J2000Eq;" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "J2000Eq.Origin = " << missionbodies_unique[0].central_body_name << ";" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "J2000Eq.Axes = MJ2000Eq;" << endl;
	GMATfile << endl;

	for (int b = 0; b < missionbodies_unique.size(); ++b)
	{
		GMATfile << "Create CoordinateSystem " << missionbodies_unique[b].name << "J2000Eq;" << endl;
		GMATfile << missionbodies_unique[b].name << "J2000Eq.Origin = " << missionbodies_unique[b].name << ";" << endl;
		GMATfile << missionbodies_unique[b].name << "J2000Eq.Axes = MJ2000Eq;" << endl;
		GMATfile << endl;
	}
	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_coordinatesystems() method


// method to create solver information
void gmatscripter::write_GMAT_solvers(){

	//solvers header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Solvers" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//fmincon or vf13ad
	//default optimizer is VF13
	GMATfile << "Create VF13ad NLPObject;" << endl;
	GMATfile << "NLPObject.ShowProgress = true;" << endl;
	GMATfile << "NLPObject.ReportStyle = Normal;" << endl;
	GMATfile << "NLPObject.ReportFile = 'VF13adNLPObject.data';" << endl;
	GMATfile << "NLPObject.MaximumIterations = 100;" << endl;
	GMATfile << "NLPObject.Tolerance = 1e-004;" << endl;
	GMATfile << "NLPObject.UseCentralDifferences = false;" << endl;
	GMATfile << "NLPObject.FeasibilityTolerance = 0.1;" << endl;
	GMATfile << endl;

	//uncomment for use with fmincon
	//DEBUG: should have option for using VF13ad or Fmincon
	GMATfile << "%Uncomment for use with fmincon" << endl;
	GMATfile << "%Create FminconOptimizer NLPObject;" << endl;
	GMATfile << "%NLPObject.DiffMaxChange = '0.1000';" << endl;
	GMATfile << "%NLPObject.DiffMinChange = '1.0000e-08';" << endl;
	GMATfile << "%NLPObject.MaxFunEvals = '1000';" << endl;
	GMATfile << "%NLPObject.TolX = '1.0000e-06';" << endl;
	GMATfile << "%NLPObject.TolFun = '1.0000e-06';" << endl;
	GMATfile << "%NLPObject.TolCon = '1.0000e-06';" << endl;
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_solvers() method


// method to create subscriber information
void gmatscripter::write_GMAT_subscribers(){

	//subscribers header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Subscribers" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//includes orbit views, plots, etc...
	//add central universe body view
	GMATfile << "%Create subscriber for central body view" << endl;
	GMATfile << "Create OrbitView " << missionbodies_unique[0].central_body_name << "View" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.ShowPlot =		true" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.SolverIterations =	 All" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.RelativeZOrder =	501" << endl;

	//add which bodies and s/c to plot 
	GMATfile << missionbodies_unique[0].central_body_name << "View.Add =	{";
	for (int index_SC = 0; index_SC < spacecraft_names.size(); ++index_SC)
	{
		GMATfile << spacecraft_names[index_SC] << ", ";
	}
	for (int index_body_visited = 0; index_body_visited < missionbodies_unique.size(); ++index_body_visited)
	{
		GMATfile << missionbodies_unique[index_body_visited].name;
		if (index_body_visited < missionbodies_unique.size() - 1)
		{
			GMATfile << ", ";
		}
	}
	GMATfile << ", " << missionbodies_unique[0].central_body_name << "}" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.CoordinateSystem =		" << missionbodies_unique[0].central_body_name << "J2000Eq" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.DrawObject = [";

	//create flag parameters for plotting
	for (int index_plot = 0; index_plot < (missionbodies_unique.size() + spacecraft_names.size()); ++index_plot)
	{
		GMATfile << "true ";
	}
	GMATfile << " true]" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.DataCollectFrequency   = 1" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.UpdatePlotFrequency    = 50" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.NumPointsToRedraw      = 300" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.ViewScaleFactor        = 35" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.ViewPointReference	  = " << missionbodies_unique[0].central_body_name << ";" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.ViewDirection		  = " << missionbodies_unique[0].central_body_name << ";" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.ViewPointVector		  = [ 0 0 30000000 ];" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.ViewUpAxis             = X" << endl;
	GMATfile << missionbodies_unique[0].central_body_name << "View.UseInitialView         = On" << endl;
	GMATfile << endl;

	//create orbit views for all bodies visited
	GMATfile << "%Create subscribers for other body views" << endl;
	for (int index_body_visited = 0; index_body_visited < missionbodies_unique.size(); ++index_body_visited)
	{
		GMATfile << "Create OrbitView " << missionbodies_unique[index_body_visited].name << "View" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.ShowPlot               = true" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.SolverIterations       = All" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.RelativeZOrder         = 501" << endl;

		//add which bodies and s/c to plot
		GMATfile << missionbodies_unique[index_body_visited].name << "View.Add                    = {";
		for (int index_SC = 0; index_SC < spacecraft_names.size(); ++index_SC)
		{
			GMATfile << spacecraft_names[index_SC] << ", ";
		}
		for (int index_body_visited = 0; index_body_visited < missionbodies_unique.size(); ++index_body_visited)
		{
			GMATfile << missionbodies_unique[index_body_visited].name;
			if (index_body_visited < missionbodies_unique.size() - 1)
			{
				GMATfile << ", ";
			}
		}
		GMATfile << "}" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.CoordinateSystem       = " << missionbodies_unique[index_body_visited].name << "J2000Eq" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.DrawObject             = [";

		//create flag parameters for plotting
		for (int index_plot = 0; index_plot < (missionbodies_unique.size() + spacecraft_names.size()); ++index_plot)
		{
			GMATfile << "true ";
		}
		GMATfile << "]" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.DataCollectFrequency   = 1" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.UpdatePlotFrequency    = 50" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.NumPointsToRedraw      = 300" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.ViewScaleFactor        = 35" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.ViewUpAxis             = X" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.UseInitialView         = On" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.ViewPointReference	  = " << missionbodies_unique[index_body_visited].name << ";" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.ViewDirection		  = " << missionbodies_unique[index_body_visited].name << ";" << endl;
		GMATfile << missionbodies_unique[index_body_visited].name << "View.ViewPointVector		  = [ 0 0 3000000 ];" << endl;
		GMATfile << endl;
	}
	GMATfile << endl;
	GMATfile << endl;

	//create reports for debugging purposes
	GMATfile << "%Create reports for debugging purposes" << endl;
	GMATfile << "Create ReportFile Report_SpacecraftState;" << endl;
	GMATfile << "Report_SpacecraftState.SolverIterations = Current;" << endl;
	GMATfile << "Report_SpacecraftState.UpperLeft = [ 0 0 ];" << endl;
	GMATfile << "Report_SpacecraftState.Size = [ 0 0 ];" << endl;
	GMATfile << "Report_SpacecraftState.RelativeZOrder = 0;" << endl;
	GMATfile << "Report_SpacecraftState.Maximized = false;" << endl;
	GMATfile << "Report_SpacecraftState.Filename = 'ReportFile1.txt';" << endl;
	GMATfile << "Report_SpacecraftState.Precision = 16;" << endl;
	GMATfile << "Report_SpacecraftState.WriteHeaders = true;" << endl;
	GMATfile << "Report_SpacecraftState.LeftJustify = On;" << endl;
	GMATfile << "Report_SpacecraftState.ZeroFill = Off;" << endl;
	GMATfile << "Report_SpacecraftState.ColumnWidth = 20;" << endl;
	GMATfile << "Report_SpacecraftState.WriteReport = true;" << endl;

	GMATfile << "Create ReportFile Report_SpacecraftControl;" << endl;
	GMATfile << "Report_SpacecraftControl.SolverIterations = Current;" << endl;
	GMATfile << "Report_SpacecraftControl.UpperLeft = [ 0 0 ];" << endl;
	GMATfile << "Report_SpacecraftControl.Size = [ 0 0 ];" << endl;
	GMATfile << "Report_SpacecraftControl.RelativeZOrder = 0;" << endl;
	GMATfile << "Report_SpacecraftControl.Maximized = false;" << endl;
	GMATfile << "Report_SpacecraftControl.Filename = 'ReportFile2.txt';" << endl;
	GMATfile << "Report_SpacecraftControl.Precision = 16;" << endl;
	GMATfile << "Report_SpacecraftControl.WriteHeaders = true;" << endl;
	GMATfile << "Report_SpacecraftControl.LeftJustify = On;" << endl;
	GMATfile << "Report_SpacecraftControl.ZeroFill = Off;" << endl;
	GMATfile << "Report_SpacecraftControl.ColumnWidth = 20;" << endl;
	GMATfile << "Report_SpacecraftControl.WriteReport = true;" << endl;
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_subscribers() method








} // end of EMTG namespace