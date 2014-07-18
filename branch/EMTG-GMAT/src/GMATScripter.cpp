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
#include "EMTG_math.h"

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
	this->create_GMAT_file();
	
	// collect the bodies used during the GMAT mission
	this->get_GMAT_bodieslist();

	// collect the launchdate used for the GMAT mission
	this->get_GMAT_missionlevelparameters();

	// write out the preamble
	this->write_GMAT_preamble();

	// write out the spacecraft
	this->write_GMAT_spacecraft();

	// write out the spacecraft hardware (i.e. thrusters and fuel tanks)
	this->write_GMAT_hardware();
	
	// write out additional models for bodies to be visited that are not planets
	this->write_GMAT_nonstandardbody();

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

	// write out the state inital guess
	this->write_GMAT_initialguess();

	// write out the mission initial boundary conditions
	this->write_GMAT_initialboundaryconditions();

	// write out the mission propagation
	this->write_GMAT_missionpropagate();

	// write out the mission final boundary conditions
	this->write_GMAT_finalboundaryconditions();

	// write out the objective function for the optimizer
	this->write_GMAT_objectivefunction();

}


// method to create a GMAT script file
void gmatscripter::create_GMAT_file(){

	//create a filename
	string filename = this->ptr_gmatmission->options.working_directory + "//" + 
					  this->ptr_gmatmission->options.mission_name + "_" + 
					  this->ptr_gmatmission->options.description + "_GMAT_Test.script";
	//open the ofstream object called GMATfile
	this->GMATfile.open(filename.c_str(), ios::trunc);
	//set floating point decimal precision
	this->GMATfile.precision(25);

	//create a debug file
	filename.erase(filename.begin(), filename.end());
	filename = this->ptr_gmatmission->options.working_directory + "//" +
   			   this->ptr_gmatmission->options.mission_name + "_" +
		       this->ptr_gmatmission->options.description + "_GMATDebug.script";
	//open the ofstream object called GMATfile
	this->GMATDebug.open(filename.c_str(), ios::trunc);
	//set floating point decimal precision
	this->GMATDebug.precision(10);

}


// method to parse the bodies used in the mission
void gmatscripter::get_GMAT_bodieslist(){

	//declaration
	int start_int = 0;
	double boundary_state[6];
	double state_difference[3];
	double state_difference_norm;

	//a loop structure over journeys and phases
	//	1.) collect all the bodies for the mission 
	//	2.) flag each body for ''
	for (int j = 0; j < this->ptr_gmatmission->number_of_journeys; ++j)
	{
		if (j > 0) { start_int = 1; }
		//add bodies of each phase for each journey to the missionbodies vector
		for (int p = start_int; p < this->ptr_gmatmission->journeys[j].number_of_phases + 1; ++p) {
			//push back body onto the vector
			missionbodies.push_back(this->ptr_gmatmission->TheUniverse[j].bodies[this->ptr_gmatmission->options.sequence[j][p] - 1]);
			//figure out if the spacecraft is within the SOI of the body
			if (p < this->ptr_gmatmission->journeys[j].number_of_phases) {
				//populate the boundary_state vector with the states of the body
				missionbodies.back().locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);
				//calculate whether the spacecraft for this phase will start within the SOI of the body
				for (int index = 0; index < 3; ++index) {
					state_difference[index] = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[index] - boundary_state[index];
				}
				state_difference_norm = math::norm(state_difference, 3);
				//compare the norm of the state to the bodies SOI radius
				if (state_difference_norm < missionbodies.back().r_SOI) { isSpaceCraftInSOI.push_back(true); }
				else { isSpaceCraftInSOI.push_back(false); }
			}
			else {
				//populate the boundary_state vector with the states of the body
				missionbodies.back().locate_body(this->ptr_gmatmission->journeys[j].phases[p - 1].phase_end_epoch, boundary_state, false, &this->ptr_gmatmission->options);
				//calculate whether the spacecraft for this phase will end within the SOI of the body
				for (int index = 0; index < 3; ++index) {
					state_difference[index] = this->ptr_gmatmission->journeys[j].phases[p - 1].state_at_end_of_phase[index] - boundary_state[index];
				}
				state_difference_norm = math::norm(state_difference, 3);
				//compare the norm of the state to the bodies SOI radius
				if (state_difference_norm < missionbodies.back().r_SOI) { isSpaceCraftInSOI.push_back(true); }
				else { isSpaceCraftInSOI.push_back(false); }
			}
			// for a departure phase, look at the journey conditions and flag the body for 'useCentralBodyInSOI'
			if (p == 0) {
				//    Journey Departure Types
				switch (this->ptr_gmatmission->options.journey_departure_type[j]) {
				case 0: //launch or direct insertion
					useCentralBodyInSOI.push_back(true);
					break;
				case 1: //depart from parking orbit(you can use this one in place of a launch vehicle model, and the departure burn will be done with the EDS motor)
					useCentralBodyInSOI.push_back(false);
					break;
				case 2: //free direct departure, i.e. do not burn to get the departure v_infinity(used for when operations about a small body are not modeled but the departure velocity is known)
					useCentralBodyInSOI.push_back(true);
					break;
				case 3: //flyby(only valid for successive journeys)
					useCentralBodyInSOI.push_back(false);
					break;
				case 4: //flyby with fixed v - infinity - out(only valid for successive journeys)
					useCentralBodyInSOI.push_back(false);
					break;
				case 5: //spiral - out from circular orbit(low - thrust missions only)
					useCentralBodyInSOI.push_back(false);
					break;
				case 6: //zero - turn flyby(for small bodies)
					useCentralBodyInSOI.push_back(false);
					break;
				}//end of switch
			}//end of if(phase == 0)
			// for an arrival phase, look at the journey conditions and flag the body for 'useCentralBodyInSOI'
			else if (p == this->ptr_gmatmission->journeys[j].number_of_phases) {
				//    Journey Arrival Types
				switch (this->ptr_gmatmission->options.journey_arrival_type[j]) {
				case 0: // insertion into parking orbit(use chemical Isp)
					useCentralBodyInSOI.push_back(false);
					break;
				case 1: // rendezvous(use chemical Isp)
					useCentralBodyInSOI.push_back(true);
					break;
				case 2: // intercept with bounded V_infinity
					useCentralBodyInSOI.push_back(true);
					break;
				case 3: // low - thrust rendezvous(does not work if terminal phase is not low - thrust)
					useCentralBodyInSOI.push_back(true);
					break;
				case 4: // match final v - infinity vector
					useCentralBodyInSOI.push_back(true);
					break;
				case 5: // match final v - infinity vector(low - thrust)
					useCentralBodyInSOI.push_back(true);
					break;
				case 6: // escape(E = 0)
					useCentralBodyInSOI.push_back(true);
					break;
				case 7: // capture spiral
					useCentralBodyInSOI.push_back(true);
					break;
				}//end of switch
			}//end of else if(phase == n - 1)
			// else the phase body is a flyby, so include 
			else {
				useCentralBodyInSOI.push_back(false);
			}
		}//end of phase for-statement

		//MyTODO:: A perturb list of bodies for each journey to use in the .PointMasses field of Force Models
		//string some_string = "blah blah";
		//forcemodel_pointmass_list.pushback( some_string );

	}//end of journeys for-statement


	//TODO:: add other arrival/departure types
	//for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	//{
	//	//add bodies of first journey to the mission body list
	//	if (j == 0)
	//	{
	//		for (int p = 0; p < this->ptr_gmatmission->options.number_of_phases[0] + 1; ++p)
	//			missionbodies.push_back(this->ptr_gmatmission->TheUniverse[0].bodies[this->ptr_gmatmission->options.sequence[0][p] - 1]);
	//	}
	//	//add all other bodies to mission body list
	//	else
	//	{
	//		for (int p = 1; p < this->ptr_gmatmission->options.number_of_phases[j] + 1; ++p)
	//			missionbodies.push_back(this->ptr_gmatmission->TheUniverse[j].bodies[this->ptr_gmatmission->options.sequence[j][p] - 1]);
	//	}
	//}


	int body_flag = 0;
	//remove duplicate entries in visited bodies list
	for (int body_index = 0; body_index < missionbodies.size(); ++body_index)
	{
		body_flag = 0;
		for (int body_index_unique = 0; body_index_unique < missionbodies_unique.size(); ++body_index_unique)
		{
			if (missionbodies[body_index].spice_ID == missionbodies_unique[body_index_unique].spice_ID)
			{
				body_flag = 1;
			}
		}
		//if body flag not switch 'on', then the body is unique; add it to the missionbodies_unique vector
		if (body_flag == 0)
		{
			missionbodies_unique.push_back(missionbodies[body_index]);

			//push back the 'isbodySOI' and 'isbodyFlyby' to appropriate vectors
			//use .size() then .erase( .begin(), .begin() + .size()) after push_back() is complete

		}
	}

	//make sure that each body doesnt start with a number for GMAT's sake
	for (int body_index = 0; body_index < missionbodies.size(); ++body_index)
	{
		if ((missionbodies[body_index].name[0] == 0) || (missionbodies[body_index].name[0] == 1) || (missionbodies[body_index].name[0] == 2) || (missionbodies[body_index].name[0] == 3) || (missionbodies[body_index].name[0] == 4) || (missionbodies[body_index].name[0] == 5) || (missionbodies[body_index].name[0] == 6) || (missionbodies[body_index].name[0] == 7) || (missionbodies[body_index].name[0] == 8) || (missionbodies[body_index].name[0] == 9))
		{
			missionbodies[body_index].name = "A" + missionbodies[body_index].name;
		}
	}
	for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index)
	{
		if ((missionbodies_unique[body_index].name[0] == 0) || (missionbodies_unique[body_index].name[0] == 1) || (missionbodies_unique[body_index].name[0] == 2) || (missionbodies_unique[body_index].name[0] == 3) || (missionbodies_unique[body_index].name[0] == 4) || (missionbodies_unique[body_index].name[0] == 5) || (missionbodies_unique[body_index].name[0] == 6) || (missionbodies_unique[body_index].name[0] == 7) || (missionbodies_unique[body_index].name[0] == 8) || (missionbodies_unique[body_index].name[0] == 9))
		{
			missionbodies_unique[body_index].name = "A" + missionbodies_unique[body_index].name;
		}
	}

	GMATDebug << "Central Bodies: ";
	for (int b = 0; b < missionbodies.size(); b++){
		GMATDebug << missionbodies[b].central_body_name << ", ";
	}
	GMATDebug << endl;
	GMATDebug << endl;

	GMATDebug << "Mission Bodies: ";
	for (int b = 0; b < missionbodies.size(); b++){
		GMATDebug << missionbodies[b].name << ", ";
	}
	GMATDebug << endl;
	GMATDebug << endl;

	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		GMATDebug << "Perturb ThirdBody: " << this->ptr_gmatmission->options.perturb_thirdbody << endl;
		int numperb = this->ptr_gmatmission->options.journey_number_of_perturbation_bodies[j];
		GMATDebug << "j" << j << " number of perturbation bodies: " << numperb << endl;
		GMATDebug << "bodies: ";
		for (int num = 0; num < numperb; ++num) {
			GMATDebug << this->ptr_gmatmission->options.journey_perturbation_bodies[j][num] << ", ";
		}
		GMATDebug << endl;
		GMATDebug << endl;
	}

}//end of get_GMAT_bodieslist() method


// method to get the mission level parameters  
void gmatscripter::get_GMAT_missionlevelparameters(){

	//declarations
	stringstream prefixstream;
	string prefix;
	int body_index = 0;

	//find index in Xdescriptions where the time bounds for the mission are located
	//it should always be the first entry, but for safety we will check the decision vector description until we 
	//find the appropriate identifier 'j0p0: launch epoch (MJD)'
	for (int iX = 0; iX < this->ptr_gmatmission->Xdescriptions.size(); ++iX)
	{
		if (this->ptr_gmatmission->Xdescriptions[iX] == "j0p0: launch epoch (MJD)")
		{
			//define these bounds temporarily
			LaunchDate_LowerBounds = this->ptr_gmatmission->Xlowerbounds[iX];
			LaunchDate_UpperBounds = this->ptr_gmatmission->Xupperbounds[iX];
		}
	}

	//find the lower and upper bounds for both distance and velocity during flybys
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{

			//construct a string that identifies the j#p#
			prefix.clear();
			prefixstream << "j" << j << "p" << p << ": ";
			prefix = prefixstream.str();

			if ((p == 0) || (p == this->ptr_gmatmission->journeys[j].number_of_phases - 1))
			{
				//TODO:: hardcoded in for journey arrivals and departures :-/
				if (missionbodies[body_index].mass < 1.0e25)
				{
					Forward_Flyby_Distance_LowerBound.push_back(-10 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
					Forward_Flyby_Distance_UpperBound.push_back( 10 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
				}
				else
				{
					Forward_Flyby_Distance_LowerBound.push_back(-300 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
					Forward_Flyby_Distance_UpperBound.push_back( 300 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
				}
				if (missionbodies[body_index + 1].mass < 1.0e25)
				{
					Backward_Flyby_Distance_LowerBound.push_back(-10 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
					Backward_Flyby_Distance_UpperBound.push_back( 10 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
				}
				else
				{
					Backward_Flyby_Distance_LowerBound.push_back(-300 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
					Backward_Flyby_Distance_UpperBound.push_back( 300 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
				}
				Forward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + (-25) * (-25)));
				Forward_Flyby_Velocity_UpperBound.push_back( sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + 25 * 25));
				Backward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + (-25) * (-25)));
				Backward_Flyby_Velocity_UpperBound.push_back( sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + 25 * 25));
			}
			else
			{
				for (int iX = 0; iX < this->ptr_gmatmission->Xdescriptions.size(); ++iX)
				{
					//find index in Xdescriptions where the flyby altitude bounds for that phase are located
					if (this->ptr_gmatmission->Xdescriptions[iX] == prefix + "flyby altitude constraint (above minimum altitude but below [100x/300x] altitude for [rocky/gas] planets")
					{
						Forward_Flyby_Distance_LowerBound.push_back( this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
						Forward_Flyby_Distance_UpperBound.push_back(-this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
						Backward_Flyby_Distance_LowerBound.push_back( this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
						Backward_Flyby_Distance_UpperBound.push_back(-this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
					}
					//find index in Xdescriptions where the flyby velocity bounds for that phase are located
					if (this->ptr_gmatmission->Xdescriptions[iX] == prefix + "initial velocity increment x")
					{
						Forward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xlowerbounds[iX] * this->ptr_gmatmission->Xlowerbounds[iX]));
						Forward_Flyby_Velocity_UpperBound.push_back( sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xupperbounds[iX] * this->ptr_gmatmission->Xupperbounds[iX]));
						Backward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xlowerbounds[iX] * this->ptr_gmatmission->Xlowerbounds[iX]));
						Backward_Flyby_Velocity_UpperBound.push_back( sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xupperbounds[iX] * this->ptr_gmatmission->Xupperbounds[iX]));
					}
				}
			}
			//increment
			body_index++;
		}//end of phases for-statement
	}//end of journeys for-statement

}//end of get_GMAT_missionlevelparameters() method


// method to get the phase level parameters
void gmatscripter::get_GMAT_phaselevelparameters(){

	//record whether each of the 





	//for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	//{
	//	for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
	//	{
	//	}
	//}

}


// method to create the script preamble
void gmatscripter::write_GMAT_preamble(){

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
			spacecraft_forward_names.push_back(spacecraft_name.str());

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
			spacecraft_backward_names.push_back(spacecraft_name.str());

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


// method to create models for bodies that are non-standard (i.e. not a planet nor the moon)
void gmatscripter::write_GMAT_nonstandardbody(){

	//force model header 
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- NonStandard Body Models" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//declarations
	string filestring;
	vector <fs::path> SPICE_files;
	SPICEDOUBLE_CELL(spice_coverage, 10000);
	SpiceInt number_of_windows = 0;

	//for each body in mission
	for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index)
	{
		//must create any bodies that are not already defined in GMAT
		if ((missionbodies[body_index].name != "Sun") && (missionbodies[body_index].name != "Mercury") && (missionbodies[body_index].name != "Venus") && (missionbodies[body_index].name != "Earth") && (missionbodies[body_index].name != "Mars") && (missionbodies[body_index].name != "Jupiter") && (missionbodies[body_index].name != "Saturn") && (missionbodies[body_index].name != "Uranus") && (missionbodies[body_index].name != "Neptune") && (missionbodies[body_index].name != "Pluto"))
		{
			GMATfile << "%Must create model for body visited" << endl;
			GMATfile << "Create Planet " << missionbodies_unique[body_index].name << endl;
			GMATfile << missionbodies_unique[body_index].name << ".NAIFId = " << missionbodies_unique[body_index].spice_ID << endl;
			GMATfile << missionbodies_unique[body_index].name << ".EquatorialRadius = " << missionbodies_unique[body_index].radius << endl;
			GMATfile << missionbodies_unique[body_index].name << ".Mu = " << missionbodies_unique[body_index].mu << endl;
			GMATfile << missionbodies_unique[body_index].name << ".PosVelSource = 'SPICE'" << endl;
			GMATfile << missionbodies_unique[body_index].name << ".CentralBody = '" << missionbodies_unique[body_index].central_body_name << "'" << endl;
			GMATfile << missionbodies_unique[body_index].name << ".RotationDataSource = 'IAUSimplified'" << endl;
			GMATfile << missionbodies_unique[body_index].name << ".OrientationEpoch = 21545" << endl;
			GMATfile << missionbodies_unique[body_index].name << ".SpinAxisRAConstant = " << missionbodies_unique[body_index].J2000_body_equatorial_frame.alpha0 << endl;
			GMATfile << missionbodies_unique[body_index].name << ".SpinAxisRARate = " << missionbodies_unique[body_index].J2000_body_equatorial_frame.alphadot << endl;
			GMATfile << missionbodies_unique[body_index].name << ".SpinAxisDECConstant = " << missionbodies_unique[body_index].J2000_body_equatorial_frame.delta0 << endl;
			GMATfile << missionbodies_unique[body_index].name << ".SpinAxisDECRate = " << missionbodies_unique[body_index].J2000_body_equatorial_frame.deltadot << endl;
			GMATfile << missionbodies_unique[body_index].name << ".RotationConstant = " << missionbodies_unique[body_index].J2000_body_equatorial_frame.W << endl;
			GMATfile << missionbodies_unique[body_index].name << ".RotationRate = " << missionbodies_unique[body_index].J2000_body_equatorial_frame.Wdot << endl;
			GMATfile << missionbodies_unique[body_index].name << ".OrbitSpiceKernelName = {";

			//find which spice file the body is located
			EMTG::filesystem::get_all_files_with_extension(fs::path(this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

			for (size_t k = 0; k < SPICE_files.size(); ++k)
			{
				filestring = this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

				//check if body is located in spice file
				scard_c(0, &spice_coverage);
				spkcov_c(filestring.c_str(), missionbodies_unique[body_index].spice_ID, &spice_coverage);
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
		if ((missionbodies_unique[body_index].central_body_name != "Sun") && (missionbodies_unique[body_index].central_body_name != "Mercury") && (missionbodies_unique[body_index].central_body_name != "Venus") && (missionbodies_unique[body_index].central_body_name != "Earth") && (missionbodies_unique[body_index].central_body_name != "Mars") && (missionbodies_unique[body_index].central_body_name != "Jupiter") && (missionbodies_unique[body_index].central_body_name != "Saturn") && (missionbodies_unique[body_index].central_body_name != "Uranus") && (missionbodies_unique[body_index].central_body_name != "Neptune") && (missionbodies_unique[body_index].central_body_name != "Pluto"))
		{
			GMATfile << "%Must create model for central body" << endl;
			GMATfile << "Create Planet " << missionbodies_unique[body_index].central_body_name << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".NAIFId = " << missionbodies_unique[body_index].central_body_spice_ID << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".EquatorialRadius = " << missionbodies_unique[body_index].central_body_radius << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".Mu = " << missionbodies_unique[body_index].universe_mu << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".PosVelSource = 'SPICE'" << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".CentralBody = 'Sun'" << endl; //assume Sun for now
			GMATfile << missionbodies_unique[body_index].central_body_name << ".RotationDataSource = 'IAUSimplified'" << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".OrientationEpoch = 21545" << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".SpinAxisRAConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alpha0 << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".SpinAxisRARate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alphadot << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".SpinAxisDECConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.delta0 << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".SpinAxisDECRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.deltadot << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".RotationConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.W << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".RotationRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.Wdot << endl;
			GMATfile << missionbodies_unique[body_index].central_body_name << ".OrbitSpiceKernelName = {";

			//find which spice file the body is located
			EMTG::filesystem::get_all_files_with_extension(fs::path(this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

			for (size_t k = 0; k < SPICE_files.size(); ++k)
			{
				filestring = this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

				//check if body is located in spice file
				scard_c(0, &spice_coverage);
				spkcov_c(filestring.c_str(), missionbodies_unique[body_index].central_body_spice_ID, &spice_coverage);
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
	}
}//end of write_GMAT_nonstandardbody() method


// method to create forcemodel information
void gmatscripter::write_GMAT_forcemodels(){

	//force model header 
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Force Models" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//create central force model (TODO:: assume for now that first body's central body is the principle central body)
	GMATfile << "%Create 2-body force model for central body" << endl;
	this->aux_GMAT_forcemodel(missionbodies_unique[0].central_body_name + "2Bod", 
							  missionbodies_unique[0].central_body_name, 
							  missionbodies_unique[0].central_body_name);

	GMATfile << "%Create 2-body force models for all bodies visited" << endl;
	for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++) {
		this->aux_GMAT_forcemodel(missionbodies_unique[body_index].name + "2Bod",
								  missionbodies_unique[body_index].name,
								  missionbodies_unique[body_index].name);
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
	this->aux_GMAT_propagate(missionbodies_unique[0].central_body_name + "Prop", 
							 missionbodies_unique[0].central_body_name + "2Bod", 
							 false);

	//create propagators for close body approaches, smaller time steps
	GMATfile << "%Create propagation models for other bodies" << endl;
	for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index) {
		this->aux_GMAT_propagate(missionbodies_unique[body_index].name + "Prop",
			missionbodies_unique[body_index].name + "2Bod",
			true);
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
			GMATfile << endl;
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

			//create temporary strings to identify time steps during phases
			for (int time_step = 0; time_step < this->ptr_gmatmission->options.num_timesteps; time_step++){
				GMATfile << "Create String tempString_j" << j << "p" << p << "s" << time_step << "_tminus" << endl;
				GMATfile << "tempString_j" << j << "p" << p << "s" << time_step << "_tminus = 'j" << j << "p" << p << "s" << time_step << "_tminus: '" << endl;
				GMATfile << "Create String tempString_j" << j << "p" << p << "s" << time_step << "_tplus" << endl;
				GMATfile << "tempString_j" << j << "p" << p << "s" << time_step << "_tplus = 'j" << j << "p" << p << "s" << time_step << "_tplus: '" << endl;
			}

		}//end of phase for-statement
	}//end of journeys for-statement

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
	for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index)
	{
		GMATfile << missionbodies_unique[body_index].name;
		if (body_index < missionbodies_unique.size() - 1)
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
	for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index)
	{
		GMATfile << "Create OrbitView " << missionbodies_unique[body_index].name << "View" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.ShowPlot               = true" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.SolverIterations       = All" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.RelativeZOrder         = 501" << endl;

		//add which bodies and s/c to plot
		GMATfile << missionbodies_unique[body_index].name << "View.Add                    = {";
		for (int index_SC = 0; index_SC < spacecraft_names.size(); ++index_SC)
		{
			GMATfile << spacecraft_names[index_SC] << ", ";
		}
		for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index)
		{
			GMATfile << missionbodies_unique[body_index].name;
			if (body_index < missionbodies_unique.size() - 1)
			{
				GMATfile << ", ";
			}
		}
		GMATfile << "}" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.CoordinateSystem       = " << missionbodies_unique[body_index].name << "J2000Eq" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.DrawObject             = [";

		//create flag parameters for plotting
		for (int index_plot = 0; index_plot < (missionbodies_unique.size() + spacecraft_names.size()); ++index_plot)
		{
			GMATfile << "true ";
		}
		GMATfile << "]" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.DataCollectFrequency   = 1" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.UpdatePlotFrequency    = 50" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.NumPointsToRedraw      = 300" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.ViewScaleFactor        = 35" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.ViewUpAxis             = X" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.UseInitialView         = On" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.ViewPointReference	  = " << missionbodies_unique[body_index].name << ";" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.ViewDirection		  = " << missionbodies_unique[body_index].name << ";" << endl;
		GMATfile << missionbodies_unique[body_index].name << "View.ViewPointVector		  = [ 0 0 3000000 ];" << endl;
		GMATfile << endl;
	}
	GMATfile << endl;
	GMATfile << endl;

	//declarations
	string report_name;
	stringstream report_name_stream;

	//create reports for debugging purposes
	GMATfile << "%Create reports for debugging purposes" << endl;
	GMATfile << "%Create a report for the central body and each unique body of the mission" << endl;

	report_name_stream << "Report_Spacecraft_" << missionbodies[0].central_body_name << "_States";
	report_name = report_name_stream.str();
	GMATfile << "Create ReportFile " << report_name << endl;
	GMATfile << report_name << ".SolverIterations = Current;" << endl;
	GMATfile << report_name << ".UpperLeft = [ 0 0 ];" << endl;
	GMATfile << report_name << ".Size = [ 0 0 ];" << endl;
	GMATfile << report_name << ".RelativeZOrder = 0;" << endl;
	GMATfile << report_name << ".Maximized = false;" << endl;
	GMATfile << report_name << ".Filename = 'Report_" << missionbodies[0].central_body_name << "Centered_States.txt';" << endl;
	GMATfile << report_name << ".Precision = 16;" << endl;
	GMATfile << report_name << ".WriteHeaders = true;" << endl;
	GMATfile << report_name << ".LeftJustify = On;" << endl;
	GMATfile << report_name << ".ZeroFill = Off;" << endl;
	GMATfile << report_name << ".ColumnWidth = 20;" << endl;
	GMATfile << report_name << ".WriteReport = true;" << endl;
	GMATfile << endl;

	for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++){
		// new report name
		report_name.erase (report_name.begin(), report_name.end());
		//report_name_stream << "Report_Spacecraft_" << missionbodies_unique[body_index].name << "_States";
		//report_name = report_name_stream.str();
		report_name = "Report_Spacecraft_" + missionbodies_unique[body_index].name + "_States";
		// GMAT script printing
		GMATfile << "Create ReportFile " << report_name << endl;
		GMATfile << report_name << ".SolverIterations = Current;" << endl;
		GMATfile << report_name << ".UpperLeft = [ 0 0 ];" << endl;
		GMATfile << report_name << ".Size = [ 0 0 ];" << endl;
		GMATfile << report_name << ".RelativeZOrder = 0;" << endl;
		GMATfile << report_name << ".Maximized = false;" << endl;
		GMATfile << report_name << ".Filename = 'Report_" << missionbodies_unique[body_index].name << "Centered_States.txt';" << endl;
		GMATfile << report_name << ".Precision = 16;" << endl;
		GMATfile << report_name << ".WriteHeaders = true;" << endl;
		GMATfile << report_name << ".LeftJustify = On;" << endl;
		GMATfile << report_name << ".ZeroFill = Off;" << endl;
		GMATfile << report_name << ".ColumnWidth = 20;" << endl;
		GMATfile << report_name << ".WriteReport = true;" << endl;
		GMATfile << endl;
	}

	GMATfile << "Create ReportFile Report_SpacecraftControl;" << endl;
	GMATfile << "Report_SpacecraftControl.SolverIterations = Current;" << endl;
	GMATfile << "Report_SpacecraftControl.UpperLeft = [ 0 0 ];" << endl;
	GMATfile << "Report_SpacecraftControl.Size = [ 0 0 ];" << endl;
	GMATfile << "Report_SpacecraftControl.RelativeZOrder = 0;" << endl;
	GMATfile << "Report_SpacecraftControl.Maximized = false;" << endl;
	GMATfile << "Report_SpacecraftControl.Filename = 'Report_SpaceCraftControlHistory.txt';" << endl;
	GMATfile << "Report_SpacecraftControl.Precision = 16;" << endl;
	GMATfile << "Report_SpacecraftControl.WriteHeaders = true;" << endl;
	GMATfile << "Report_SpacecraftControl.LeftJustify = On;" << endl;
	GMATfile << "Report_SpacecraftControl.ZeroFill = Off;" << endl;
	GMATfile << "Report_SpacecraftControl.ColumnWidth = 20;" << endl;
	GMATfile << "Report_SpacecraftControl.WriteReport = true;" << endl;
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_subscribers() method


// method to create the initial state guesses for the mission sequence
void gmatscripter::write_GMAT_initialguess(){

	//declarations
	double ThrustUnitVector_lowerbounds;
	double ThrustUnitVector_upperbounds;
	double boundary_state[6];
	math::Matrix<double> Vinf_in(3, 1);
	math::Matrix<double> Vinf_out(3, 1);
	math::Matrix<double> periapse_state_vector(6, 1);
	math::Matrix<double> periapse_position_vector(3, 1);
	math::Matrix<double> periapse_velocity_vector(3, 1);
	int body_index = 0;
	int name_index = 0;

	//begin mission sequence
	GMATfile << "BeginMissionSequence" << endl;
	GMATfile << endl;
	GMATfile << endl;

	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Initial State Guesses" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	GMATfile << "BeginScript 'Initial Guess Values' " << endl;
	GMATfile << endl;

	//TODO:: hook up engine inputs
	GMATfile << "	%Engine model parameters" << endl;
	if (this->ptr_gmatmission->options.mission_type < 2 || this->ptr_gmatmission->options.mission_type == 5) //impulsive-thrust phase types
	{
		GMATfile << "	ThrusterISP = " << this->ptr_gmatmission->options.IspChem << endl;

		GMATfile << endl;

		//for each journey
		for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
			{
				GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".K1 = ThrusterISP;" << endl;
				GMATfile << endl;
			}
		}
	}
	else //for all low-thrust phase types
	{
		switch (this->ptr_gmatmission->options.engine_type)
		{
		case 0: //fixed thrust and ISP
			GMATfile << "	ThrusterISP = "       << this->ptr_gmatmission->options.IspLT             << endl;
			GMATfile << "	ThrusterMaxThrust = " << this->ptr_gmatmission->options.Thrust            << endl;
			GMATfile << "	ThrusterDutyCycle = " << this->ptr_gmatmission->options.engine_duty_cycle << endl;
			break;
		case 1: //constant ISP, efficiency, EMTG chooses power
			break;
		case 2: //choice of power model, constant efficiency, EMTG chooses ISP
			break;
		case 3: //choice of power model, constant efficiency and ISP
			break;
		case 4: //continuously-varying specific impulse (constant efficiency)
			break;
		case 5: //custom thrust and mass flow rate
			break;
		case 6: //NSTAR
			break;
		case 7: //XIPS-25
			break;
		case 8: //BPT-4000 high ISP
			break;
		case 9: //BPT-4000 high thrust
			break;
		case 10: //BPT-4000 ex-high ISP
			break;
		case 11: //NEXT high ISP
			break;
		case 12: //VASIMR
			break;
		case 13: //Hall thruster
			break;
		case 14: //NEXT v10 high Isp
			break;
		case 15: //NEXT v10 high Thrust
			break;
		case 16: //BPT-4000 MALTO curve
			break;
		}

		GMATfile << endl;

		//for each journey
		for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
			{
				GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.DutyCycle  = ThrusterDutyCycle;"	<< endl;
				GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1         = ThrusterMaxThrust;" << endl;
				GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.K1         = ThrusterISP;"	    << endl;
				GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.DutyCycle = ThrusterDutyCycle;" << endl;
				GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.C1        = ThrusterMaxThrust;" << endl;
				GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.K1        = ThrusterISP;"	    << endl;
				GMATfile << endl;
			}
		}
	} //end code for low-thrust models

	//EMTG-GMAT needs to reference the upper and lower bounds vectors, so we must recalculate them
	//according to Jacob, this shouldn't be necessary anymore: 2014_07_15
	//this->ptr_gmatmission->calcbounds();

	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//FORWARD PROPAGATED SPACECRAFT
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Forward s/c initial conditions" << endl;
			GMATDebug << " I should be printing right now." << endl;
			//interpret state from beginning of each journey based on departure type (//TODO::)
			if (p == 0)
			{
				GMATDebug << "p: " << p << endl;
				//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
				missionbodies[body_index].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);

				switch (this->ptr_gmatmission->options.journey_departure_type[j])
				{
				case 0: //launch or direct insertion // (EMTG default mission)
					GMATDebug << "case: 0" << endl;
					//calculate v_infinity vectors
					for (int k = 0; k < 3; ++k)
					{
						Vinf_out(k) = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
					}
					//calculate inc from vinfinity, then make guess at min altitude at planet
					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[0] - boundary_state[0] << ";" << endl;
					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[1] - boundary_state[1] << ";" << endl;
					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[2] - boundary_state[2] << ";" << endl;
					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[3] - boundary_state[3] << ";" << endl;
					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[4] - boundary_state[4] << ";" << endl;
					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[5] - boundary_state[5] << ";" << endl;
					GMATDebug << "Initial States: ";
					for (int counter = 0; counter < 6; counter++){
						GMATDebug << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[counter] << ", ";
					}
					GMATDebug << endl;
					break;
				case 1: //depart from parking orbit
					GMATDebug << "case: 1" << endl;
					//TODO::
					//journey_departure_elements_type == 1 means given in COE, which is only relevant if 
					//departure type is 'depart from parking orbit' or 'free - direct departure' OR
					//arrival type is 'insertion into parking orbit (use chemical Isp)' or ...

					//I dont think John's code works correctly
					//periapse_state_vector = journeys[j].phases[p].calculate_periapse_state_from_asymptote_and_parking_orbit(Vinf_out, options.journey_departure_elements[j][2], options.journey_departure_elements[j][0] - missionbodies[body_index].radius, journeys[j].phases[p].phase_start_epoch, &TheUniverse[options.number_of_journeys - 1], body_index + 1);

					//can convert to inertial and then to appropriate frame
					//COE2inertial(const double* E_COE, const double mu, double* state)
					//reference frame transformation here if needed
					if (this->ptr_gmatmission->options.journey_departure_elements_type[j] == 1)
					{
						for (int k = 0; k < 3; ++k)
						{
							periapse_position_vector(k) = periapse_state_vector(k);
							periapse_velocity_vector(k) = periapse_state_vector(k + 3);
						}
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = " << periapse_position_vector(0) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = " << periapse_position_vector(1) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = " << periapse_position_vector(2) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << periapse_velocity_vector(0) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << periapse_velocity_vector(1) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << periapse_velocity_vector(2) << ";" << endl;
					}
					break;
				case 2: //free direct departure
					GMATDebug << "case: 2" << endl;
					//TODO::
					//journey_departure_elements_type == 1 means given in COE, which is only relevant if 
					//departure type is 'depart from parking orbit' or 'free - direct departure' OR
					//arrival type is 'insertion into parking orbit (use chemical Isp)' or ...

					//I dont think John's code works correctly
					//periapse_state_vector = journeys[j].phases[p].calculate_periapse_state_from_asymptote_and_parking_orbit(Vinf_out, options.journey_departure_elements[j][2], options.journey_departure_elements[j][0] - missionbodies[body_index].radius, journeys[j].phases[p].phase_start_epoch, &TheUniverse[options.number_of_journeys - 1], body_index + 1);

					//can convert to inertial and then to appropriate frame
					//COE2inertial(const double* E_COE, const double mu, double* state)
					//reference frame transformation here if needed
					if (this->ptr_gmatmission->options.journey_departure_elements_type[j] == 1)
					{
						for (int k = 0; k < 3; ++k)
						{
							periapse_position_vector(k) = periapse_state_vector(k);
							periapse_velocity_vector(k) = periapse_state_vector(k + 3);
						}
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = " << periapse_position_vector(0) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = " << periapse_position_vector(1) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = " << periapse_position_vector(2) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << periapse_velocity_vector(0) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << periapse_velocity_vector(1) << ";" << endl;
						GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << periapse_velocity_vector(2) << ";" << endl;
					}
					break;
				case 3: //depart from flyby
					GMATDebug << "case: 3" << endl;
					break;
				}
			}//end of if(p == 0)

			//for all other phases, interpret state from beginning of each phase and then add the minimum flyby altitude
			else
			{
				//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
				missionbodies[body_index].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);

				//calculate v_infinity vectors
				for (int k = 0; k < 3; ++k)
				{
					Vinf_in(k)  = this->ptr_gmatmission->journeys[j].phases[p - 1].state_at_end_of_phase[k + 3]   - boundary_state[k + 3];
					Vinf_out(k) = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
				}

				periapse_state_vector = this->ptr_gmatmission->journeys[j].phases[p].calculate_flyby_periapse_state(Vinf_in, Vinf_out, this->ptr_gmatmission->journeys[j].phases[p - 1].flyby_altitude, missionbodies[body_index]);

				//add this position vector to state's initial guess
				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = "  << periapse_state_vector(0) << ";" << endl;
				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = "  << periapse_state_vector(1) << ";" << endl;
				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = "  << periapse_state_vector(2) << ";" << endl;
				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << periapse_state_vector(3) << ";" << endl;
				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << periapse_state_vector(4) << ";" << endl;
				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << periapse_state_vector(5) << ";" << endl;
			}
			GMATfile << "	" << spacecraft_forward_names[name_index] << ".CoordinateSystem = " << missionbodies[body_index].name << "J2000Eq;" << endl;
			GMATfile << endl;


			//BACKWARD PROPAGATED SPACECRAFT
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward s/c initial conditions" << endl;

			//interpret state from end of each journey based on arrival type (//TODO::)
			if (p == this->ptr_gmatmission->journeys[j].number_of_phases - 1)
			{
				//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
				missionbodies[body_index + 1].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_end_epoch, boundary_state, false, &this->ptr_gmatmission->options);

				switch (this->ptr_gmatmission->options.journey_arrival_type[j])
				{
				case 0: //parking orbit insertion
					break;
				case 1: //rendezvous (chemical)
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5] << ";" << endl;
					break;
				case 2: //flyby with bounded VHP
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5] << ";" << endl;
					break;
				case 3: //rendezvous (LT) // (EMTG default mission)
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4] << ";" << endl;
					GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5] << ";" << endl;
					break;
				case 4: //match Vinf vector (chemical)
					break;
				case 5: //match Vinf vector (LT)
					break;
				}
			}

			//for all other phases, interpret state from end of each phase and then add the minimum flyby altitude
			else
			{

				//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
				missionbodies[body_index + 1].locate_body(this->ptr_gmatmission->journeys[j].phases[p + 1].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);

				//calculate v_infinity vectors
				for (int k = 0; k < 3; ++k)
				{
					Vinf_in(k)  = this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[k + 3]           - boundary_state[k + 3];
					Vinf_out(k) = this->ptr_gmatmission->journeys[j].phases[p + 1].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
				}

				periapse_state_vector = this->ptr_gmatmission->journeys[j].phases[p].calculate_flyby_periapse_state(Vinf_in, Vinf_out, this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude, missionbodies[body_index + 1]);

				//add this position vector to state's initial guess
				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << periapse_state_vector(0) << ";" << endl;
				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << periapse_state_vector(1) << ";" << endl;
				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << periapse_state_vector(2) << ";" << endl;
				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << periapse_state_vector(3) << ";" << endl;
				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << periapse_state_vector(4) << ";" << endl;
				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << periapse_state_vector(5) << ";" << endl;
			}
			GMATfile << "	" << spacecraft_backward_names[name_index] << ".CoordinateSystem = " << missionbodies[body_index + 1].name << "J2000Eq;" << endl;
			GMATfile << endl;

			//insert scaled launch epoch and wait times for each subsequent journey
			if ((j > 0) && (p == 0))
			{
				GMATfile << "	%Guess for scaled journey wait times" << endl;
				GMATfile << "	Journey" << j + 1 << "_WaitTime_Scaled = " << ((this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch - this->ptr_gmatmission->journeys[j - 1].phases[this->ptr_gmatmission->journeys[j - 1].number_of_phases - 1].phase_end_epoch) - this->ptr_gmatmission->options.journey_wait_time_bounds[j][0]) / (this->ptr_gmatmission->options.journey_wait_time_bounds[j][1] - this->ptr_gmatmission->options.journey_wait_time_bounds[j][0]) << ";" << endl;
			}
			if (j == 0)
			{
				if (p == 0)
				{
					GMATfile << "	%Guess for scaled spacecraft launch epoch" << endl;
					GMATfile << "	LaunchEpoch_Scaled = " << (this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch - LaunchDate_LowerBounds) / (LaunchDate_UpperBounds - LaunchDate_LowerBounds) << ";" << endl;
				}
			}
			GMATfile << endl;

			//insert scaled states for forward and backward s/c (doing the math in GMAT bc of all the if statements above!)
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Forward s/c scaled initial conditions" << endl;
			GMATfile << "	" << spacecraft_forward_names[name_index] << "_X_Scaled = ("  << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X - "  << Forward_Flyby_Distance_LowerBound[p] << ") / (" << Forward_Flyby_Distance_UpperBound[p] << " - " << Forward_Flyby_Distance_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_forward_names[name_index] << "_Y_Scaled = ("  << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y - "  << Forward_Flyby_Distance_LowerBound[p] << ") / (" << Forward_Flyby_Distance_UpperBound[p] << " - " << Forward_Flyby_Distance_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_forward_names[name_index] << "_Z_Scaled = ("  << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z - "  << Forward_Flyby_Distance_LowerBound[p] << ") / (" << Forward_Flyby_Distance_UpperBound[p] << " - " << Forward_Flyby_Distance_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_forward_names[name_index] << "_VX_Scaled = (" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX - " << Forward_Flyby_Velocity_LowerBound[p] << ") / (" << Forward_Flyby_Velocity_UpperBound[p] << " - " << Forward_Flyby_Velocity_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_forward_names[name_index] << "_VY_Scaled = (" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY - " << Forward_Flyby_Velocity_LowerBound[p] << ") / (" << Forward_Flyby_Velocity_UpperBound[p] << " - " << Forward_Flyby_Velocity_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_forward_names[name_index] << "_VZ_Scaled = (" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ - " << Forward_Flyby_Velocity_LowerBound[p] << ") / (" << Forward_Flyby_Velocity_UpperBound[p] << " - " << Forward_Flyby_Velocity_LowerBound[p] << ")" << endl;
			GMATfile << "	FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[6] / (this->ptr_gmatmission->options.maximum_mass) << endl;
			GMATfile << endl;

			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward s/c scaled initial conditions" << endl;
			GMATfile << "	" << spacecraft_backward_names[name_index] << "_X_Scaled = ("  << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X - "  << Backward_Flyby_Distance_LowerBound[p] << ") / (" << Backward_Flyby_Distance_UpperBound[p] << " - " << Backward_Flyby_Distance_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_backward_names[name_index] << "_Y_Scaled = ("  << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y - "  << Backward_Flyby_Distance_LowerBound[p] << ") / (" << Backward_Flyby_Distance_UpperBound[p] << " - " << Backward_Flyby_Distance_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_backward_names[name_index] << "_Z_Scaled = ("  << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z - "  << Backward_Flyby_Distance_LowerBound[p] << ") / (" << Backward_Flyby_Distance_UpperBound[p] << " - " << Backward_Flyby_Distance_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_backward_names[name_index] << "_VX_Scaled = (" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX - " << Backward_Flyby_Velocity_LowerBound[p] << ") / (" << Backward_Flyby_Velocity_UpperBound[p] << " - " << Backward_Flyby_Velocity_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_backward_names[name_index] << "_VY_Scaled = (" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY - " << Backward_Flyby_Velocity_LowerBound[p] << ") / (" << Backward_Flyby_Velocity_UpperBound[p] << " - " << Backward_Flyby_Velocity_LowerBound[p] << ")" << endl;
			GMATfile << "	" << spacecraft_backward_names[name_index] << "_VZ_Scaled = (" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ - " << Backward_Flyby_Velocity_LowerBound[p] << ") / (" << Backward_Flyby_Velocity_UpperBound[p] << " - " << Backward_Flyby_Velocity_LowerBound[p] << ")" << endl;
			GMATfile << "	FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[6] / (this->ptr_gmatmission->options.maximum_mass) << endl;
			GMATfile << endl;
			++name_index;
			++body_index;

			//initial guess for inter-phase control
			//this means thrust vectors (and sometimes Isp) for low-thrust phases
			//or burn index for MGADSM and MGANDSM phases. No control for MGA phases
			//initialize thrust vector directions
			//define thrust unit vector bounds to scale to between 0 and 1
			double ThrustUnitVector_lowerbounds = -1;
			double ThrustUnitVector_upperbounds =  1;

			//propagate forward s/c using finite burns (must scale unit vectors to between 0 and 1)
			GMATfile << endl;
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Forward s/c scaled thrust vectors" << endl;
			for (int step = 0; step < (this->ptr_gmatmission->options.num_timesteps / 2); ++step)
			{
				GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][0] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
				GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][1] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
				GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][2] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
				GMATfile << endl;
			}
			//propagate backward s/c using finite burns (must scale unit vectors to between 0 and 1)
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward s/c scaled thrust vectors" << endl;
			for (int step = this->ptr_gmatmission->options.num_timesteps - 1; (step >= this->ptr_gmatmission->options.num_timesteps / 2); --step)
			{
				GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << ", Time Step #" << step + 1 << ", Backward Propagation" << endl;
				GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][0] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
				GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][1] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
				GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][2] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
				GMATfile << endl;
			}
		}//end of phase for-statement
	}//end of journey for-statement
	GMATfile << "EndScript" << endl;
	GMATfile << endl;

}//end of write_GMAT_initialguess() method


// method to create the initial boundary conditions for the mission sequence
void gmatscripter::write_GMAT_initialboundaryconditions(){

	//declarations
	double CumulatedTOF = 0;
	int body_index = 0;
	int name_index = 0;

	//begin optimization loop
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%------------ Initial Boundary Constraints" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;
	GMATfile << "Optimize 'OptimizeSequence' NLPObject {SolveMode = Solve, ExitMode = DiscardAndContinue}" << endl;
	GMATfile << endl;

	//vary the launch date and journey wait times
	GMATfile << "	%Vary launch epoch and journey wait times" << endl;
	GMATfile << "	Vary 'VaryLaunchEpoch_Scaled' NLPObject(LaunchEpoch_Scaled = LaunchEpoch_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;

	//for each journey vary each subsequent journey wait time
	for (int j = 1; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		GMATfile << "	Vary 'VaryJourney" << j + 1 << "_WaitTime_Scaled' NLPObject(Journey" << j + 1 << "_WaitTime_Scaled = Journey" << j + 1 << "_WaitTime_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	}

	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//insert unscaled epoch values
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " s/c unscaled epochs" << endl;

			//if subsequent journey, add unscaled wait times
			if (j > 0)
			{
				GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << ".Epoch.TAIModJulian' " << spacecraft_forward_names[name_index] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_UpperBounds - LaunchDate_LowerBounds) << " + " << LaunchDate_LowerBounds << ") / 86400.0";
				for (int jj = 0; jj < j + 1; ++jj)
				{
					GMATfile << " + (Journey" << jj + 1 << "_WaitTime_Scaled * " << (this->ptr_gmatmission->options.journey_wait_time_bounds[jj][1] - this->ptr_gmatmission->options.journey_wait_time_bounds[jj][0]) << " + " << this->ptr_gmatmission->options.journey_wait_time_bounds[jj][0] << ")";
				}
				GMATfile << " + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
				//add time of flight for backward s/c (end of phase)
				CumulatedTOF = CumulatedTOF + this->ptr_gmatmission->journeys[j].phases[p].TOF;
				GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << ".Epoch.TAIModJulian' " << spacecraft_backward_names[name_index] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_UpperBounds - LaunchDate_LowerBounds) << " + " << LaunchDate_LowerBounds << ") / 86400.0";
				for (int jj = 0; jj < j + 1; ++jj)
				{
					GMATfile << " + (Journey" << jj + 1 << "_WaitTime_Scaled * " << (this->ptr_gmatmission->options.journey_wait_time_bounds[jj][1] - this->ptr_gmatmission->options.journey_wait_time_bounds[jj][0]) << " + " << this->ptr_gmatmission->options.journey_wait_time_bounds[jj][0] << ")";
				}
				GMATfile << " + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
			}
			//for the first journey
			else
			{
				GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << ".Epoch.TAIModJulian' " << spacecraft_forward_names[name_index] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_UpperBounds - LaunchDate_LowerBounds) << " + " << LaunchDate_LowerBounds << ") / 86400.0 + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
				//add time of flight for backward s/c (end of phase)
				CumulatedTOF = CumulatedTOF + this->ptr_gmatmission->journeys[j].phases[p].TOF;
				GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << ".Epoch.TAIModJulian' " << spacecraft_backward_names[name_index] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_UpperBounds - LaunchDate_LowerBounds) << " + " << LaunchDate_LowerBounds << ") / 86400.0 + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
			}
			name_index++;
			GMATfile << endl;
		}
	}
	//store final epoch value
	GMATfile << "	'CalcFinalEpoch' FinalEpoch = " << spacecraft_backward_names[name_index - 1] << ".Epoch.TAIModJulian" << endl;

	//reset spacecraft name index to zero
	name_index = 0;

	//add inequality constraints for min flyby altitude and ensure that each s/c begins at periapse
	GMATfile << "	%Add periapse and flyby altitude constraints for each phase" << endl;
	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			if ((p < this->ptr_gmatmission->journeys[j].number_of_phases) && (p != 0))
			{
				GMATfile << "	'Constraint" << spacecraft_forward_names[name_index] << "_PeriapseRadius_Scaled' " << spacecraft_forward_names[name_index] << "_PeriapseRadius_Scaled = " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << ".RadPer / " << missionbodies[body_index].minimum_safe_flyby_altitude << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_PeriapseRadius_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_PeriapseRadius_Scaled >= 1)" << endl;
				GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "_RdotV_Scaled' " << spacecraft_forward_names[name_index] << "_RdotV_Scaled = (" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X * " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX + " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y * " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY + " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z * " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ + " << sqrt(2 * missionbodies[body_index].mu * missionbodies[body_index].radius) << ") / " << 2 * sqrt(2 * missionbodies[body_index].mu * missionbodies[body_index].radius) << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_RdotV_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_RdotV_Scaled = 0.5);" << endl;
				GMATfile << endl;
			}
			//increment our spacecraft name and body indices
			name_index++;
			body_index++;
		}
	}

	//reset spacecraft name index to zero
	name_index = 0;
	body_index = 0;

	//vary the initial states of s/c
	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//vary all forward s/c
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " vary forward s/c scaled states" << endl;
			GMATfile << "	Vary 'Vary" << spacecraft_forward_names[name_index] << "_X_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_X_Scaled = " << spacecraft_forward_names[name_index] << "_X_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'Vary" << spacecraft_forward_names[name_index] << "_Y_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_Y_Scaled = " << spacecraft_forward_names[name_index] << "_Y_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'Vary" << spacecraft_forward_names[name_index] << "_Z_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_Z_Scaled = " << spacecraft_forward_names[name_index] << "_Z_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'Vary" << spacecraft_forward_names[name_index] << "_VX_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_VX_Scaled = " << spacecraft_forward_names[name_index] << "_VX_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'Vary" << spacecraft_forward_names[name_index] << "_VY_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_VY_Scaled = " << spacecraft_forward_names[name_index] << "_VY_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'Vary" << spacecraft_forward_names[name_index] << "_VZ_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_VZ_Scaled = " << spacecraft_forward_names[name_index] << "_VZ_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'VaryFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled' NLPObject(FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << endl;
			name_index++;
		}
	}
	//vary last s/c state for mission (backward from last body)
	GMATfile << "	%Journey #" << this->ptr_gmatmission->options.number_of_journeys << ", Phase #" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << " vary backward s/c scaled states" << endl;
	GMATfile << "	Vary 'Vary" << spacecraft_backward_names[name_index - 1] << "_X_Scaled' NLPObject(" << spacecraft_backward_names[name_index - 1] << "_X_Scaled = " << spacecraft_backward_names[name_index - 1] << "_X_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	GMATfile << "	Vary 'Vary" << spacecraft_backward_names[name_index - 1] << "_Y_Scaled' NLPObject(" << spacecraft_backward_names[name_index - 1] << "_Y_Scaled = " << spacecraft_backward_names[name_index - 1] << "_Y_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	GMATfile << "	Vary 'Vary" << spacecraft_backward_names[name_index - 1] << "_Z_Scaled' NLPObject(" << spacecraft_backward_names[name_index - 1] << "_Z_Scaled = " << spacecraft_backward_names[name_index - 1] << "_Z_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	GMATfile << "	Vary 'Vary" << spacecraft_backward_names[name_index - 1] << "_VX_Scaled' NLPObject(" << spacecraft_backward_names[name_index - 1] << "_VX_Scaled = " << spacecraft_backward_names[name_index - 1] << "_VX_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	GMATfile << "	Vary 'Vary" << spacecraft_backward_names[name_index - 1] << "_VY_Scaled' NLPObject(" << spacecraft_backward_names[name_index - 1] << "_VY_Scaled = " << spacecraft_backward_names[name_index - 1] << "_VY_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	GMATfile << "	Vary 'Vary" << spacecraft_backward_names[name_index - 1] << "_VZ_Scaled' NLPObject(" << spacecraft_backward_names[name_index - 1] << "_VZ_Scaled = " << spacecraft_backward_names[name_index - 1] << "_VZ_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	GMATfile << "	Vary 'VaryFuelMass_Journey" << this->ptr_gmatmission->options.number_of_journeys << "Phase" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << "Backward_Scaled' NLPObject(FuelMass_Journey" << this->ptr_gmatmission->options.number_of_journeys << "Phase" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << "Backward_Scaled = FuelMass_Journey" << this->ptr_gmatmission->options.number_of_journeys << "Phase" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << "Backward_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	GMATfile << "	'CalcFinalMass_Scaled' FinalMass_Scaled = FuelMass_Journey" << this->ptr_gmatmission->options.number_of_journeys << "Phase" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << "Backward_Scaled" << endl;
	GMATfile << endl;

	//reset spacecraft name index to zero
	name_index = 0;

	//now we must define states to their unscaled values
	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " set forward s/c unscaled states" << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X' " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = " << spacecraft_forward_names[name_index] << "_X_Scaled * " << (Forward_Flyby_Distance_UpperBound[p] - Forward_Flyby_Distance_LowerBound[p]) << " + " << Forward_Flyby_Distance_LowerBound[p] << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y' " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = " << spacecraft_forward_names[name_index] << "_Y_Scaled * " << (Forward_Flyby_Distance_UpperBound[p] - Forward_Flyby_Distance_LowerBound[p]) << " + " << Forward_Flyby_Distance_LowerBound[p] << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z' " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = " << spacecraft_forward_names[name_index] << "_Z_Scaled * " << (Forward_Flyby_Distance_UpperBound[p] - Forward_Flyby_Distance_LowerBound[p]) << " + " << Forward_Flyby_Distance_LowerBound[p] << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX' " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << spacecraft_forward_names[name_index] << "_VX_Scaled * " << (Forward_Flyby_Velocity_UpperBound[p] - Forward_Flyby_Velocity_LowerBound[p]) << " + " << Forward_Flyby_Velocity_LowerBound[p] << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY' " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << spacecraft_forward_names[name_index] << "_VY_Scaled * " << (Forward_Flyby_Velocity_UpperBound[p] - Forward_Flyby_Velocity_LowerBound[p]) << " + " << Forward_Flyby_Velocity_LowerBound[p] << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ' " << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << spacecraft_forward_names[name_index] << "_VZ_Scaled * " << (Forward_Flyby_Velocity_UpperBound[p] - Forward_Flyby_Velocity_LowerBound[p]) << " + " << Forward_Flyby_Velocity_LowerBound[p] << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.FuelMass' " << spacecraft_forward_names[name_index] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.FuelMass = FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled * " << (this->ptr_gmatmission->options.maximum_mass) << endl;
			GMATfile << endl;

			//match forward and backward s/c states and epochs		
			if (name_index < (spacecraft_forward_names.size() - 1))
			{
				if (p < this->ptr_gmatmission->journeys[j].number_of_phases - 1)
				{
					//match only odd phase backward s/c == subsequent even phase forward s/c
					GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " match forward and backward s/c unscaled states and epochs" << endl;
					GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X' " << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = " << spacecraft_forward_names[name_index + 1] << "." << missionbodies[body_index + 1].name << "J2000Eq.X;" << endl;
					GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y' " << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = " << spacecraft_forward_names[name_index + 1] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y;" << endl;
					GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z' " << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = " << spacecraft_forward_names[name_index + 1] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z;" << endl;
					GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX' " << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << spacecraft_forward_names[name_index + 1] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX;" << endl;
					GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY' " << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << spacecraft_forward_names[name_index + 1] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY;" << endl;
					GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ' " << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << spacecraft_forward_names[name_index + 1] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ;" << endl;
					GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << ".Epoch.TAIModJulian' " << spacecraft_backward_names[name_index] << ".Epoch.TAIModJulian = " << spacecraft_forward_names[name_index + 1] << ".Epoch.TAIModJulian;" << endl;
					if (p < this->ptr_gmatmission->journeys[j].number_of_phases - 1)
					{
						GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass' " << spacecraft_backward_names[name_index] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass = " << spacecraft_forward_names[name_index + 1] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 2 << "Forward.FuelMass" << endl;
					}
					else
					{
						GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass' " << spacecraft_backward_names[name_index] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass = " << spacecraft_forward_names[name_index + 1] << ".FuelTank_Journey" << j + 2 << "Phase" << 1 << "Forward.FuelMass" << endl;
					}
					GMATfile << endl;
				}
			}
			name_index++;
			body_index++;
		}//end of phase for-statement
	}//end of journey for-statement

	//set last backward states to unscaled values
	GMATfile << "	%Journey #" << this->ptr_gmatmission->options.number_of_journeys << ", Phase #" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << " set backward s/c unscaled states" << endl;
	GMATfile << "	'Calc" << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.X' "  << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.X = "  << spacecraft_backward_names[name_index - 1] << "_X_Scaled * "  << (Backward_Flyby_Distance_UpperBound.back() - Backward_Flyby_Distance_LowerBound.back()) << " + " << Backward_Flyby_Distance_LowerBound.back() << endl;
	GMATfile << "	'Calc" << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.Y' "  << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.Y = "  << spacecraft_backward_names[name_index - 1] << "_Y_Scaled * "  << (Backward_Flyby_Distance_UpperBound.back() - Backward_Flyby_Distance_LowerBound.back()) << " + " << Backward_Flyby_Distance_LowerBound.back() << endl;
	GMATfile << "	'Calc" << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.Z' "  << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.Z = "  << spacecraft_backward_names[name_index - 1] << "_Z_Scaled * "  << (Backward_Flyby_Distance_UpperBound.back() - Backward_Flyby_Distance_LowerBound.back()) << " + " << Backward_Flyby_Distance_LowerBound.back() << endl;
	GMATfile << "	'Calc" << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.VX' " << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << spacecraft_backward_names[name_index - 1] << "_VX_Scaled * " << (Backward_Flyby_Velocity_UpperBound.back() - Backward_Flyby_Velocity_LowerBound.back()) << " + " << Backward_Flyby_Velocity_LowerBound.back() << endl;
	GMATfile << "	'Calc" << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.VY' " << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << spacecraft_backward_names[name_index - 1] << "_VY_Scaled * " << (Backward_Flyby_Velocity_UpperBound.back() - Backward_Flyby_Velocity_LowerBound.back()) << " + " << Backward_Flyby_Velocity_LowerBound.back() << endl;
	GMATfile << "	'Calc" << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.VZ' " << spacecraft_backward_names[name_index - 1] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << spacecraft_backward_names[name_index - 1] << "_VZ_Scaled * " << (Backward_Flyby_Velocity_UpperBound.back() - Backward_Flyby_Velocity_LowerBound.back()) << " + " << Backward_Flyby_Velocity_LowerBound.back() << endl;
	GMATfile << "	'Calc" << spacecraft_backward_names[name_index - 1] << ".FuelTank_Journey" << this->ptr_gmatmission->options.number_of_journeys << "Phase" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << "Backward.FuelMass' " << spacecraft_backward_names[name_index - 1] << ".FuelTank_Journey" << this->ptr_gmatmission->options.number_of_journeys << "Phase" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << "Backward.FuelMass = FuelMass_Journey" << this->ptr_gmatmission->options.number_of_journeys << "Phase" << this->ptr_gmatmission->journeys[this->ptr_gmatmission->options.number_of_journeys - 1].number_of_phases << "Backward_Scaled * " << (this->ptr_gmatmission->options.maximum_mass) << endl;
	GMATfile << endl;
	GMATfile << endl;

}// end of write_GMAT_initialboundaryconditions() method


// method to create the mission propagation lines
void gmatscripter::write_GMAT_missionpropagate(){

	//declarations
	int name_index = 0;
	int body_index = 0;
	double elapsed_secs;

	GMATfile << "	%-------------------------------------------------------------------------" << endl;
	GMATfile << "	%---------- Propagation" << endl;
	GMATfile << "	%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//declarations
			double periapse_velocity_magnitude;
			math::Matrix<double> Vinf_out(3, 1), Vinf_in(3, 1);
			double delta_t;
			double boundary_state[6];
			int index_delta_t;

			//declare thrust unit vector bounds to scale between 0 and 1
			double ThrustUnitVector_lowerbounds = -1;
			double ThrustUnitVector_upperbounds =  1;

			//GMAT PenUp()
			this->aux_GMAT_penUp(); GMATfile << endl;

			//propagate forward s/c using finite burns (must scale unit vectors to between 0 and 1)
			for (int step = 0; step < (this->ptr_gmatmission->options.num_timesteps / 2); ++step)
			{
				//header
				GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << ", Time Step #" << step + 1 << ", Forward Propagation";
				GMATfile << endl;
				//vary thrust vector for each step
				GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				//calculate the thrust unit vector magnitude of the step
				GMATfile << "	'CalcThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << 
							"' ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = sqrt(" <<
							"(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " << 
						    "(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " << 
							"(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2);"  << endl;
				//create a nonlinear constraint for the thrust unit vector magnitude
				GMATfile << "	NonlinearConstraint 'ConstraintThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) <= 1)" << endl;
				//assign thrust unit directions to thruster
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
				//assign value to thrust coefficient
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1 = ThrusterMaxThrust * ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;

				//locate the position of the bodies at beginning and end of phase
				missionbodies[body_index].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);

				GMATDebug << "j" << j << "p" << p << "s" << step << " ";
				//calculate time spent in SOI of body and during which timestep the boundary occurs
				for (int k = 0; k < 3; ++k)
				{
					//calculate the vinf_out vector
					Vinf_out(k) = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
					GMATDebug << "  |  Vinf_out(" << k << "): " << Vinf_out(k);
				}

				//velocity magnitude at the periapse
				periapse_velocity_magnitude = sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + Vinf_out.dot(Vinf_out));
				GMATDebug << "  |  velocity at periapse: " << periapse_velocity_magnitude;
				//delta time in (days)
				delta_t = (missionbodies[body_index].r_SOI / periapse_velocity_magnitude) / 86400;
				GMATDebug << "  |  delta time: " << delta_t;
				//the step index at which the spacecraft will leave the SOI
				index_delta_t = floor(delta_t / (this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps));
				GMATDebug << "  |  step index: " << index_delta_t;

				//report things for debugging
				this->write_GMAT_report(j, p, step, spacecraft_forward_names[name_index], missionbodies[0].central_body_name, true, true, true);
				for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++){
					this->write_GMAT_report(j, p, step, spacecraft_forward_names[name_index], missionbodies_unique[body_index].name, true, true, false);
				}

				//GMAT BeginFiniteBurn()
				this->aux_GMAT_beginburn(j, p, step, spacecraft_forward_names[name_index], "Forward");
				//GMAT PenDown()
				this->aux_GMAT_penDown();

				//check if inside SOI
				//TODO:: remove second part of if-statement after fixing arrival/departures
				if ((step < index_delta_t) && (p != 0))
				{
					GMATDebug << "  |  if()" << endl;
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps;
					this->aux_GMAT_propagate(j, p, step, spacecraft_forward_names[name_index], "Forward", missionbodies[body_index].name, elapsed_secs);
				}
				//check if it leaves SOI during time step  
				else if ((step == index_delta_t) && (p != 0))
				{
					GMATDebug << "  |  else if()" << endl;
					//if so, propagate only until outside estimated SOI of body
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = delta_t - (this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps) * index_delta_t;
					this->aux_GMAT_propagate(j, p, step, spacecraft_forward_names[name_index], "Forward", missionbodies[body_index].name, elapsed_secs);

					//then propagate the rest of the timestep
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = (this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps) * (index_delta_t + 1) - delta_t;
					this->aux_GMAT_propagate(j, p, step, spacecraft_forward_names[name_index], "Forward", missionbodies[0].central_body_name, elapsed_secs);
				}
				//otherwise propagate the full timestep
				else
				{
					GMATDebug << "  |  else()" << endl;
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps;
					this->aux_GMAT_propagate(j, p, step, spacecraft_forward_names[name_index], "Forward", missionbodies[0].central_body_name, elapsed_secs);
				}

				//GMAT EndFiniteBurn()
				aux_GMAT_endburn(j, p, step, spacecraft_forward_names[name_index], "Forward");

				//report things for debugging
				this->write_GMAT_report(j, p, step, spacecraft_forward_names[name_index], missionbodies[0].central_body_name, true, false, false);
				for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++){
					this->write_GMAT_report(j, p, step, spacecraft_forward_names[name_index], missionbodies_unique[body_index].name, true, false, false);
				}

				//GMAT PenUp()
				this->aux_GMAT_penUp(); GMATfile << endl;
			}


			//propagate backward s/c using finite burns (must scale unit vectors to between 0 and 1)
			for (int step = this->ptr_gmatmission->options.num_timesteps - 1; (step >= this->ptr_gmatmission->options.num_timesteps / 2); --step)
			{
				GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << ", Time Step #" << step + 1 << ", Backward Propagation";
				GMATfile << endl;
				GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	'CalcThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "' ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = sqrt((ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)* " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " << "(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " << "(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2);" << endl;

				//assign thrust unit directions to thruster
				GMATfile << "	NonlinearConstraint 'ConstraintThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) <= 1)" << endl;
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection1' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection1 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection2' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection2 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection3' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection3 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
				GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.C1' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.C1 = ThrusterMaxThrust * ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;

				//begin propagation of finite burn

				//locate the position of the bodies at beginning and end of phase
				missionbodies[body_index + 1].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_end_epoch, boundary_state, false, &this->ptr_gmatmission->options);

				//calculate time spent in SOI of body and during which timestep the boundary occurs
				for (int k = 0; k < 3; ++k)
				{
					Vinf_in(k) = this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[k + 3] - boundary_state[k + 3];
				}
				periapse_velocity_magnitude = sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + Vinf_in.dot(Vinf_in));
				delta_t = (missionbodies[body_index + 1].r_SOI / periapse_velocity_magnitude) / 86400;
				index_delta_t = (this->ptr_gmatmission->options.num_timesteps - 1) - floor(delta_t / (this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps));

				//report things for debugging
				this->write_GMAT_report(j, p, step, spacecraft_backward_names[name_index], missionbodies[0].central_body_name, false, true, true);
				for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++){
					this->write_GMAT_report(j, p, step, spacecraft_backward_names[name_index], missionbodies_unique[body_index].name, false, true, false);
				}

				//GMAT BeginFiniteBurn()
				this->aux_GMAT_beginburn(j, p, step, spacecraft_backward_names[name_index], "Backward");
				//GMAT PenDown()
				this->aux_GMAT_penDown();

				//check if inside SOI
				//TODO:: remove second part of if-statement after fixing arrival/departures
				if ((step > index_delta_t) && (p != this->ptr_gmatmission->options.number_of_phases[j] - 1))
				{
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = -this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps;
					this->aux_GMAT_propagate(j, p, step, spacecraft_backward_names[name_index], "Backward", missionbodies[body_index + 1].name, elapsed_secs);
				}

				//check if it leaves SOI during time step  
				else if ((step == index_delta_t) && (p != this->ptr_gmatmission->options.number_of_phases[j] - 1))
				{
					//is so, propagate only until outside estimated SOI of body
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = -(delta_t - (this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps) * ((this->ptr_gmatmission->options.num_timesteps - 1) - index_delta_t));
					this->aux_GMAT_propagate(j, p, step, spacecraft_backward_names[name_index], "Backward", missionbodies[body_index + 1].name, elapsed_secs);

					//then propagate the rest of the timestep
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = -((this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps) * (this->ptr_gmatmission->options.num_timesteps - index_delta_t) - delta_t);
					this->aux_GMAT_propagate(j, p, step, spacecraft_backward_names[name_index], "Backward", missionbodies[0].central_body_name, elapsed_secs);
				}

				//otherwise propagate the full timestep
				else
				{
					//calculate elapsed time for propagation and then propagate
					elapsed_secs = -(this->ptr_gmatmission->journeys[j].phases[p].TOF / this->ptr_gmatmission->options.num_timesteps);
					this->aux_GMAT_propagate(j, p, step, spacecraft_backward_names[name_index], "Backward", missionbodies[0].central_body_name, elapsed_secs);
				}

				//GMAT EndFiniteBurn()
				aux_GMAT_endburn(j, p, step, spacecraft_backward_names[name_index], "Backward");

				//report things for debugging
				this->write_GMAT_report(j, p, step, spacecraft_backward_names[name_index], missionbodies[0].central_body_name, false, false, false);
				for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++){
					this->write_GMAT_report(j, p, step, spacecraft_backward_names[name_index], missionbodies_unique[body_index].name, false, false, false);
				}

				//GMAT PenUp()
				this->aux_GMAT_penUp(); GMATfile << endl;
			}
			name_index++;
			body_index++;

			this->aux_GMAT_penDown();
			
		}// end of phase for-statement
	}// end of journey for-statement
	GMATfile << endl;
	GMATfile << endl;

}// end of write_GMAT_missionpropagate() method


// method to create the final boundary conditions for the mission sequence
void gmatscripter::write_GMAT_finalboundaryconditions(){

	//declarations
	int name_index = 0;
	int body_index = 0;

	//apply all boundary constraints at match points
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Final Boundary Constraints" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//for each journey
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	{
		//for each phase
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
		{
			//scale final forward states for constraints
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Forward, scale final states for constraints" << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "_X_Scaled' " << spacecraft_forward_names[name_index] << "_X_Scaled = " << spacecraft_forward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.X / " << this->ptr_gmatmission->TheUniverse[j].LU << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "_Y_Scaled' " << spacecraft_forward_names[name_index] << "_Y_Scaled = " << spacecraft_forward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.Y / " << this->ptr_gmatmission->TheUniverse[j].LU << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "_Z_Scaled' " << spacecraft_forward_names[name_index] << "_Z_Scaled = " << spacecraft_forward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.Z / " << this->ptr_gmatmission->TheUniverse[j].LU << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "_VX_Scaled' " << spacecraft_forward_names[name_index] << "_VX_Scaled = " << spacecraft_forward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.VX / " << this->ptr_gmatmission->TheUniverse[j].LU << " * " << this->ptr_gmatmission->TheUniverse[j].TU << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "_VY_Scaled' " << spacecraft_forward_names[name_index] << "_VY_Scaled = " << spacecraft_forward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.VY / " << this->ptr_gmatmission->TheUniverse[j].LU << " * " << this->ptr_gmatmission->TheUniverse[j].TU << endl;
			GMATfile << "	'Calc" << spacecraft_forward_names[name_index] << "_VZ_Scaled' " << spacecraft_forward_names[name_index] << "_VZ_Scaled = " << spacecraft_forward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.VZ / " << this->ptr_gmatmission->TheUniverse[j].LU << " * " << this->ptr_gmatmission->TheUniverse[j].TU << endl;
			GMATfile << endl;

			//scale final backward states for constraints
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward, scale final states for constraints" << endl;
			GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "_X_Scaled' " << spacecraft_backward_names[name_index] << "_X_Scaled = " << spacecraft_backward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.X / " << this->ptr_gmatmission->TheUniverse[j].LU << endl;
			GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "_Y_Scaled' " << spacecraft_backward_names[name_index] << "_Y_Scaled = " << spacecraft_backward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.Y / " << this->ptr_gmatmission->TheUniverse[j].LU << endl;
			GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "_Z_Scaled' " << spacecraft_backward_names[name_index] << "_Z_Scaled = " << spacecraft_backward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.Z / " << this->ptr_gmatmission->TheUniverse[j].LU << endl;
			GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "_VX_Scaled' " << spacecraft_backward_names[name_index] << "_VX_Scaled = " << spacecraft_backward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.VX / " << this->ptr_gmatmission->TheUniverse[j].LU << " * " << this->ptr_gmatmission->TheUniverse[j].TU << endl;
			GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "_VY_Scaled' " << spacecraft_backward_names[name_index] << "_VY_Scaled = " << spacecraft_backward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.VY / " << this->ptr_gmatmission->TheUniverse[j].LU << " * " << this->ptr_gmatmission->TheUniverse[j].TU << endl;
			GMATfile << "	'Calc" << spacecraft_backward_names[name_index] << "_VZ_Scaled' " << spacecraft_backward_names[name_index] << "_VZ_Scaled = " << spacecraft_backward_names[name_index] << "." << missionbodies[0].central_body_name << "J2000Eq.VZ / " << this->ptr_gmatmission->TheUniverse[j].LU << " * " << this->ptr_gmatmission->TheUniverse[j].TU << endl;
			GMATfile << endl;

			//apply match point constraints
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward, apply match point constraints" << endl;
			GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_X_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_X_Scaled = " << spacecraft_backward_names[name_index] << "_X_Scaled)" << endl;
			GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_Y_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_Y_Scaled = " << spacecraft_backward_names[name_index] << "_Y_Scaled)" << endl;
			GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_Z_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_Z_Scaled = " << spacecraft_backward_names[name_index] << "_Z_Scaled)" << endl;
			GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_VX_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_VX_Scaled = " << spacecraft_backward_names[name_index] << "_VX_Scaled)" << endl;
			GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_VY_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_VY_Scaled = " << spacecraft_backward_names[name_index] << "_VY_Scaled)" << endl;
			GMATfile << "	NonlinearConstraint 'Constraint" << spacecraft_forward_names[name_index] << "_VZ_Scaled' NLPObject(" << spacecraft_forward_names[name_index] << "_VZ_Scaled = " << spacecraft_backward_names[name_index] << "_VZ_Scaled)" << endl;
			GMATfile << endl;

			//apply fuel mass constraints
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " apply fuel mass constraint at match point" << endl;
			GMATfile << "	'CalcFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled' FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = " << spacecraft_forward_names[name_index] << ".TotalMass / " << (2 * this->ptr_gmatmission->options.maximum_mass) << endl;
			GMATfile << "	'CalcFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled' FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled = " << spacecraft_forward_names[name_index] << ".TotalMass / " << (2 * this->ptr_gmatmission->options.maximum_mass) << endl;
			GMATfile << "	NonlinearConstraint 'ConstraintFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "' NLPObject(FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled)" << endl;
			GMATfile << endl;
			GMATfile << endl;
			
			//increment
			name_index++;
			body_index++;
		}//end of phase for-statement
	}//end of journey for-statement

}//end of write_GMAT_finalboundaryconditions() method


// method to create the objective function
void gmatscripter::write_GMAT_objectivefunction(){

	//declarations
	int TotalNumberTimeSteps = 0;

	//optimize the user-defined objective function
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Objective Function" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	switch (this->ptr_gmatmission->options.objective_type)
	{
	case 0: // minimize delta-V
		GMATfile << "	%Objective function is to minimize delta-V" << endl;
		GMATfile << "	'CalcObjectiveFunction' ObjectiveFunction = (0";
		for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
		{
			for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
			{
				for (int step = 0; step < this->ptr_gmatmission->options.num_timesteps; ++step)
				{
					GMATfile << " + ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)";
					++TotalNumberTimeSteps;
				}
			}
		}
		GMATfile << ") / " << TotalNumberTimeSteps << endl;
		break;
	case 1: //minimize time
		GMATfile << "	%Objective function is to minimize time" << endl;
		GMATfile << "	'CalcObjectiveFunction' ObjectiveFunction = (FinalEpoch - LaunchEpoch_Scaled * " << (LaunchDate_UpperBounds - LaunchDate_LowerBounds) << " + " << LaunchDate_LowerBounds << ") / " << this->ptr_gmatmission->options.total_flight_time_bounds[1] << endl;
		break;
	case 2: //maximize final mass
		GMATfile << "	%Objective function is to maximize final mass" << endl;
		GMATfile << "	'CalcObjectiveFunction' ObjectiveFunction = -FinalMass_Scaled" << endl;
		break;
	case 3: //GTOC 1 asteroid deflection function
		break;
	case 4: //launch as late as possible in window
		break;
	case 5: //launch as early as possible in window
		break;
	case 6: //maximize orbit energy
		break;
	case 7: //minimize launch mass
		break;
	case 8: //arrive as early as possible
		break;
	case 9: //arrive as late as possible
		break;
	case 10: //minimum propellant
		break;
	}

	//currently vf13 does not like to obey set bounds if an objective function is given
	GMATfile << "	%Minimize objective function" << endl;
	GMATfile << "	%Minimize NLPObject(ObjectiveFunction)" << endl;
	GMATfile << endl;
	GMATfile << "EndOptimize" << endl;

}// end of write_GMAT_objectivefunction() method


// method to report results for debugging
void gmatscripter::write_GMAT_report(int j, int p, int s, string spacecraft_name, string body_name, 
									 bool isforwardspacecraft, bool isbeforemaneuver, bool writecontrolhistory){

	// declaration
	string prefix;
	string tempString;
	stringstream tempStringStream;

	// prefix to be added if Forward or Backward Spacecraft
	if (isforwardspacecraft) { prefix = "Forward";}
	else { prefix = "Backward"; }

	// prefix to be added if before a finite burn maneuver
	if (isbeforemaneuver) { 
		tempStringStream << "tempString_j" << j << "p" << p << "s" << s << "_tminus ";
		tempString = tempStringStream.str();
	}
	else {
		tempStringStream << "tempString_j" << j << "p" << p << "s" << s << "_tplus ";
		tempString = tempStringStream.str();
	}

	// write the spacecraft controls
	if (writecontrolhistory) {
		GMATfile << "	Report 'Report_SpacecraftControl' Report_SpacecraftControl " << tempString;
		GMATfile << spacecraft_name << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << prefix << ".ThrustDirection1 ";
		GMATfile << spacecraft_name << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << prefix << ".ThrustDirection2 ";
		GMATfile << spacecraft_name << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << prefix << ".ThrustDirection3 ";
		GMATfile << spacecraft_name << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << prefix << ".C1" << endl;
	}
	// write the spacecraft states in the central body frame
	GMATfile << "	Report 'Report_SpacecraftState' Report_Spacecraft_" << body_name << "_States " << tempString;
	GMATfile << spacecraft_name << "." << body_name << "J2000Eq.X ";
	GMATfile << spacecraft_name << "." << body_name << "J2000Eq.Y ";
	GMATfile << spacecraft_name << "." << body_name << "J2000Eq.Z ";
	GMATfile << spacecraft_name << "." << body_name << "J2000Eq.VX ";
	GMATfile << spacecraft_name << "." << body_name << "J2000Eq.VY ";
	GMATfile << spacecraft_name << "." << body_name << "J2000Eq.VZ ";
	GMATfile << spacecraft_name << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << prefix << ".FuelMass;" << endl;
	//GMATfile << "	Report 'Report_SpacecraftState' Report_SpacecraftState " << tempString;
	//GMATfile << spacecraft_name << "." << missionbodies[0].central_body_name << "J2000Eq.X ";
	//GMATfile << spacecraft_name << "." << missionbodies[0].central_body_name << "J2000Eq.Y ";
	//GMATfile << spacecraft_name << "." << missionbodies[0].central_body_name << "J2000Eq.Z ";
	//GMATfile << spacecraft_name << "." << missionbodies[0].central_body_name << "J2000Eq.VX ";
	//GMATfile << spacecraft_name << "." << missionbodies[0].central_body_name << "J2000Eq.VY ";
	//GMATfile << spacecraft_name << "." << missionbodies[0].central_body_name << "J2000Eq.VZ ";
	//GMATfile << spacecraft_name << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << prefix << ".FuelMass;" << endl;
	// write the spacecraft states in the unique body frames
	//for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++) {
	//	GMATfile << spacecraft_name << "." << missionbodies_unique[body_index].name << "J2000Eq.X ";
	//	GMATfile << spacecraft_name << "." << missionbodies_unique[body_index].name << "J2000Eq.Y ";
	//	GMATfile << spacecraft_name << "." << missionbodies_unique[body_index].name << "J2000Eq.Z ";
	//	GMATfile << spacecraft_name << "." << missionbodies_unique[body_index].name << "J2000Eq.VX ";
	//	GMATfile << spacecraft_name << "." << missionbodies_unique[body_index].name << "J2000Eq.VY ";
	//	GMATfile << spacecraft_name << "." << missionbodies_unique[body_index].name << "J2000Eq.VZ ";
	//}
	

}


// method to write a GMAT BeginFiniteBurn Command
void gmatscripter::aux_GMAT_beginburn(int j, int p, int s, string spacecraft_name, string prefix){

	GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << s + 1 
						  << prefix << "' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << prefix << "(" << spacecraft_name << ");" << endl;

}


// method to write a GMAT EndFiniteBurn Command
void gmatscripter::aux_GMAT_endburn(int j, int p, int s, string spacecraft_name, string prefix){

	GMATfile << "	EndFiniteBurn 'EndFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << s + 1 
					  << prefix << "' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << prefix << "(" << spacecraft_name << ");" << endl;

}


// method to write a GMAT Propagate Command
void gmatscripter::aux_GMAT_propagate(int j, int p, int s, string spacecraft_name, string prefix, string body_name, double elapsed_secs){

	//declarations
	string str;

	if (elapsed_secs < 0.0){ str = "BackProp "; }
	else { str = ""; }

	if (elapsed_secs == 0.0){
		GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << s + 1 << prefix << "' " << str
			<< body_name << "Prop(" << spacecraft_name << ");" << endl;
	}
	else {
		GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << s + 1 << prefix << "' " << str
			<< body_name << "Prop(" << spacecraft_name << ") {" << spacecraft_name << ".ElapsedSecs = " << elapsed_secs << "};" << endl;
	}

}


// method to write a GMAT PenUp Command
void gmatscripter::aux_GMAT_penUp(){

	GMATfile << "	PenUp 'PenUp' " << missionbodies[0].central_body_name << "View";
	for (int body_index = 0; body_index < missionbodies.size(); ++body_index)
	{
		GMATfile << " " << missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << endl;

}


// method to write a GMAT PenDown Command
void gmatscripter::aux_GMAT_penDown(){

	GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
	for (int body_index = 0; body_index < missionbodies.size(); ++body_index)
	{
		GMATfile << " " << missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << endl;

}


// method to write a GMAT ForceModel Resource
void gmatscripter::aux_GMAT_forcemodel(string forcemodelname, string centralbody, string pointmasses){

	GMATfile << "Create ForceModel " << forcemodelname << endl;
	GMATfile << forcemodelname << ".CentralBody = " << centralbody << endl;
	GMATfile << forcemodelname << ".PointMasses = {" << pointmasses << "};" << endl;
	GMATfile << forcemodelname << ".Drag = None;" << endl;
	GMATfile << forcemodelname << ".SRP = Off;" << endl;
	GMATfile << endl;

}


// method to write a GMAT Propagator Resource
void gmatscripter::aux_GMAT_propagate(string propagatorname, string forcemodelname, bool isCloseApproach) {

	//declarations
	double initialstepsize = 60.0;
	double maxstepsize = 86400.0;

	if (isCloseApproach) {
		initialstepsize = 30.0;
		maxstepsize = 8640.0;
	}

	GMATfile << "Create Propagator " << propagatorname << endl;
	GMATfile << propagatorname << ".FM = " << forcemodelname << endl;
	GMATfile << propagatorname << ".Type = PrinceDormand78; " << endl;
	GMATfile << propagatorname << ".InitialStepSize = " << initialstepsize << endl;
	GMATfile << propagatorname << ".Accuracy = 1e-11; " << endl;
	GMATfile << propagatorname << ".MinStep = 0.0; " << endl;
	GMATfile << propagatorname << ".MaxStep = " << maxstepsize << endl;
	GMATfile << endl;

}


} // end of EMTG namespace