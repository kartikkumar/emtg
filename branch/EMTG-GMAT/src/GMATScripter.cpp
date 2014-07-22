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

	// collect the step level parameters for the GMAT mission
	this->get_GMAT_steplevelparameters();

	// collect the mission level parameters for the GMAT mission
	this->get_GMAT_missionlevelparameters();

	// collect the journey level parameters for the GMAT mission
	this->get_GMAT_journeylevelparameters();

	// collect the phase level parameters for the GMAT mission
	this->get_GMAT_phaselevelparameters();

	// write out the preamble
	this->write_GMAT_preamble();

	// write out the spacecraft hardware (i.e. thrusters and fuel tanks)
	this->write_GMAT_hardware();

	// write out the spacecraft
	this->write_GMAT_spacecraft();

	// write out additional models for bodies to be visited that are not planets
	this->write_GMAT_nonstandardbody();

	// write out the forcemodels
	this->write_GMAT_forcemodels();

	// write out the propagators
	this->write_GMAT_propagators();

	// write out the burn objects
	this->write_GMAT_burns();
	
	// write out coordinate systems
	this->write_GMAT_coordinatesystems();

	// write out the solvers
	this->write_GMAT_solvers();

	// write out subscribers (i.e. plot views and reports)
	this->write_GMAT_subscribers();

	// write out arrays and variables
	this->write_GMAT_variables();

	// write the beginmissionsequence command
	this->write_GMAT_beginmissionsequence();

	// write out the state inital guess
	this->write_GMAT_initialguess();

	// write out the optimization phase of the mission sequence
	this->write_GMAT_optimization();

	// write out the mission initial boundary conditions
	//this->write_GMAT_initialboundaryconditions();

	// write out the mission propagation
	//this->write_GMAT_missionpropagate();

	// write out the mission final boundary conditions
	//this->write_GMAT_finalboundaryconditions();

	// write out the objective function for the optimizer
	//this->write_GMAT_objectivefunction();

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
void gmatscripter::get_GMAT_missionlevelparameters() {

	//declarations
	stringstream tempstream;
	vector <string> tempvector;

	//what type of mission is being solved
	//TODO: What about mission type 5? and confirm mission type 6-9 results in a mission type of 0-4.
	if (this->ptr_gmatmission->options.mission_type == 0 || this->ptr_gmatmission->options.mission_type == 1 || this->ptr_gmatmission->options.mission_type == 4) {
		isLT   = false;
		isFBLT = false;
	}
	else if (this->ptr_gmatmission->options.mission_type == 2) {
		isLT   = true;
		isFBLT = false;
	}
	else if (this->ptr_gmatmission->options.mission_type == 3) {
		isLT   = true;
		isFBLT = true;
	}

	////find index in Xdescriptions where the time bounds for the mission are located
	////it should always be the first entry, but for safety we will check the decision vector description until we 
	////find the appropriate identifier 'j0p0: launch epoch (MJD)'
	//for (int iX = 0; iX < this->ptr_gmatmission->Xdescriptions.size(); ++iX) {
	//	if (this->ptr_gmatmission->Xdescriptions[iX] == "j0p0: launch epoch (MJD)") {
	//		//define these bounds temporarily
	//		LaunchDate_LowerBounds = this->ptr_gmatmission->Xlowerbounds[iX];
	//		LaunchDate_UpperBounds = this->ptr_gmatmission->Xupperbounds[iX];
	//		LaunchDate = this->ptr_gmatmission->Xopt[iX];
	//		GMATDebug << "ix: " << iX << endl;
	//	}
	//}

	LaunchDate_LowerBounds = this->ptr_gmatmission->Xlowerbounds[0];
	LaunchDate_UpperBounds = this->ptr_gmatmission->Xupperbounds[0];
	LaunchWindow = LaunchDate_UpperBounds - LaunchDate_LowerBounds;
	LaunchDate = this->ptr_gmatmission->Xopt[0];

	GMATDebug << " ------------ " << endl;
	GMATDebug << "LaunchDate_LowerBounds " << LaunchDate_LowerBounds << endl;
	GMATDebug << "LaunchDate_UpperBounds " << LaunchDate_UpperBounds << endl;
	GMATDebug << "LaunchDate " << LaunchDate << endl;
	//GMATDebug << this->ptr_gmatmission->X[0] << endl;
	//GMATDebug << this->ptr_gmatmission->X0[0] << endl;
	//GMATDebug << this->ptr_gmatmission->Xopt[0] << endl;
	GMATDebug << " ------------ " << endl;

	// ---------------------------
	//      VARIABLE CREATION
	// ---------------------------
	//store mission level variables that will be needed during the optimization sequence in GMAT
	tempvector.push_back("ObjectiveFunction");
	mission_level_variables.push_back(tempvector);
	tempvector.clear();
	//launch window date
	tempvector.push_back("LaunchWindowOpenDate");
	tempstream << " = " << (LaunchDate_LowerBounds / 86400) + TAIModJOffset;
	tempvector.push_back(tempstream.str());
	mission_level_variables.push_back(tempvector);
	tempstream.str("");
	tempvector.clear();
	//launch window delta
	tempvector.push_back("LaunchWindow");
	tempstream << " = " << (LaunchWindow) / 86400.0;
	tempvector.push_back(tempstream.str());
	mission_level_variables.push_back(tempvector);
	tempstream.str("");
	tempvector.clear();
	//launch window scaling
	tempvector.push_back("LaunchWindowScaling");
	tempstream << " = " << (LaunchDate - LaunchDate_LowerBounds) / (LaunchWindow);
	tempvector.push_back(tempstream.str());
	mission_level_variables.push_back(tempvector);
	tempstream.str("");
	tempvector.clear();
	
	
	//this->ptr_gmatmission->options.total_flight_time_bounds;


	//variables depenedent on the objective function
	//if (this->ptr_gmatmission->options.objective_type == 2) {
	//	mission_level_variables.push_back("FinalMass_Scaled");
	//}

	//  Sometimes ENGINE Variables are Mission Type (i.e. fixed thrust Isp or Impulsive)
	//+ OTHERWISE SEE "get_GMAT_missionlevelparameters" for creation of 'phase_level_parameters'
	//if() low thrust & fixed thrust-Isp, else() impulsive
	if (isLT && this->ptr_gmatmission->options.engine_type == 0) {
		//Isp
		tempvector.push_back("ThrusterISP");
		tempstream << " = " << this->ptr_gmatmission->options.IspLT;
		tempvector.push_back(tempstream.str());
		mission_level_variables.push_back(tempvector);
		tempstream.str("");
		tempvector.clear();
		//DutyCycle
		tempvector.push_back("ThrusterDutyCycle");
		tempstream << " = " << this->ptr_gmatmission->options.Thrust;
		tempvector.push_back(tempstream.str());
		mission_level_variables.push_back(tempvector);
		tempstream.str("");
		tempvector.clear();
		//MaxThrust
		tempvector.push_back("ThrusterMaxThrust");
		tempstream << " = " << this->ptr_gmatmission->options.engine_duty_cycle;
		tempvector.push_back(tempstream.str());
		mission_level_variables.push_back(tempvector);
		tempstream.str("");
		tempvector.clear();
	}
	else if (!isLT) {
		tempvector.push_back("ThrusterISP");
		tempstream << " = " << this->ptr_gmatmission->options.IspChem;
		tempvector.push_back(tempstream.str());
		mission_level_variables.push_back(tempvector);
		tempstream.str("");
		tempvector.clear();
	}

}//end of get_GMAT_missionlevelparameters() method


// method to get the journey level parameters
void gmatscripter::get_GMAT_journeylevelparameters() {

	for (int index = 0; index < this->ptr_gmatmission->options.number_of_journeys; ++index) {
		
		// ---------------------------
		//      VARIABLE CREATION
		// ---------------------------
		// journey time bounds (0: unbounded, 1: bounded flight time, 2: bounded arrival date)
		if (this->ptr_gmatmission->options.journey_timebounded[index] == 1) {
			//journey date, delta, and scaling
			//GMATfile << "Create Variable Journey" << j + 1 << "_WaitTime_Scaled" << endl;
		}
	}
}


// method to get the phase level parameters
void gmatscripter::get_GMAT_phaselevelparameters() {

	//declarations
	stringstream tempstream;
	int body_index = 0;

	//find the lower and upper bounds for both distance and velocity during flybys
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {

			//construct a string that identifies the j#p#
			tempstream << "j" << j << "p" << p << ": ";
			

			if ((p == 0) || (p == this->ptr_gmatmission->journeys[j].number_of_phases - 1)) {
				//TODO:: hardcoded in for journey arrivals and departures :-/
				if (missionbodies[body_index].mass < 1.0e25) {
					Forward_Flyby_Distance_LowerBound.push_back(-10 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
					Forward_Flyby_Distance_UpperBound.push_back(10 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
				}
				else {
					Forward_Flyby_Distance_LowerBound.push_back(-300 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
					Forward_Flyby_Distance_UpperBound.push_back(300 * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
				}
				if (missionbodies[body_index + 1].mass < 1.0e25) {
					Backward_Flyby_Distance_LowerBound.push_back(-10 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
					Backward_Flyby_Distance_UpperBound.push_back(10 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
				}
				else {
					Backward_Flyby_Distance_LowerBound.push_back(-300 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
					Backward_Flyby_Distance_UpperBound.push_back(300 * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
				}
				Forward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + (-25) * (-25)));
				Forward_Flyby_Velocity_UpperBound.push_back(sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + 25 * 25));
				Backward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + (-25) * (-25)));
				Backward_Flyby_Velocity_UpperBound.push_back(sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + 25 * 25));
			}
			else {
				for (int iX = 0; iX < this->ptr_gmatmission->Xdescriptions.size(); ++iX) {
					//find index in Xdescriptions where the flyby altitude bounds for that phase are located
					if (this->ptr_gmatmission->Xdescriptions[iX] == tempstream.str() + "flyby altitude constraint (above minimum altitude but below [100x/300x] altitude for [rocky/gas] planets")
					{
						Forward_Flyby_Distance_LowerBound.push_back(this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
						Forward_Flyby_Distance_UpperBound.push_back(-this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index].minimum_safe_flyby_altitude + missionbodies[body_index].radius));
						Backward_Flyby_Distance_LowerBound.push_back(this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
						Backward_Flyby_Distance_UpperBound.push_back(-this->ptr_gmatmission->Xlowerbounds[iX] * (missionbodies[body_index + 1].minimum_safe_flyby_altitude + missionbodies[body_index + 1].radius));
					}
					//find index in Xdescriptions where the flyby velocity bounds for that phase are located
					if (this->ptr_gmatmission->Xdescriptions[iX] == tempstream.str() + "initial velocity increment x") {
						Forward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xlowerbounds[iX] * this->ptr_gmatmission->Xlowerbounds[iX]));
						Forward_Flyby_Velocity_UpperBound.push_back(sqrt(2 * missionbodies[body_index].mu / (missionbodies[body_index].radius + missionbodies[body_index].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xupperbounds[iX] * this->ptr_gmatmission->Xupperbounds[iX]));
						Backward_Flyby_Velocity_LowerBound.push_back(-sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xlowerbounds[iX] * this->ptr_gmatmission->Xlowerbounds[iX]));
						Backward_Flyby_Velocity_UpperBound.push_back(sqrt(2 * missionbodies[body_index + 1].mu / (missionbodies[body_index + 1].radius + missionbodies[body_index + 1].minimum_safe_flyby_altitude) + this->ptr_gmatmission->Xupperbounds[iX] * this->ptr_gmatmission->Xupperbounds[iX]));
					}
				}
			}
			//increment
			body_index++;
			//clear 'tempstream'
			tempstream.str("");

			// ---------------------------
			//      VARIABLE CREATION
			// ---------------------------
			//tempstream << "% Journey #" << j << ", Phase #" << p << " Time Variables: ";
			//phase_level_variables.push_back(tempstream.str());
			//tempstream.str("");
			//tempstream << "j" << j << "p" << p;
			//phase_level_variables.push_back(tempstream.str() + "_TOF");
			//phase_level_variables.push_back(tempstream.str() + "_TOF_Scaled");

			////create time of flight variables
			//GMATfile << "Create Variable TOF_Journey" << j + 1 << "Phase" << p + 1 << "Scaled" << endl;
			//GMATfile << "Create Variable TOF_Journey" << j + 1 << "Phase" << p + 1 << endl;
			//GMATfile << "Create Variable TimeStepLength_Journey" << j + 1 << "Phase" << p + 1 << endl;
			//GMATfile << "Create Variable LeftBoundarySOITime_Journey" << j + 1 << "Phase" << p + 1 << endl;
			//GMATfile << "Create Variable RightBoundarySOITime_Journey" << j + 1 << "Phase" << p + 1 << endl;
			//GMATfile << endl;
			//GMATfile << endl;

		}//end of phases for-statement
	}//end of journeys for-statement

}


// method to get the step level parameters
void gmatscripter::get_GMAT_steplevelparameters(){

	//declarations
	bool isForwardStep;
	int f2b_index;
	bool usePushBack;
	bool isStep_MatchPointAdjacent;
	double body_istates[6];
	double body_fstates[6];
	double iepoch;
	double fepoch;
	double stepsize;
	int body_index;
	bool isInSOI_at_Start;
	bool isInSOI_at_End;
	double initial_position_diff[3];
	double final_position_diff[3];
	double initial_velocity_diff[3];
	double final_velocity_diff[3];
	double periapse_velocity_magnitude;
	double approximate_time_in_SOI;
	stringstream tempstream;
	vector <string> forcemodel_storage;
	vector <string> propagator_storage;
	int gmatsteps_thus_far;
	int gmatsteps_thus_far_in_phase;
	int number_of_gmatsteps_added;
	vector <double> thrustvector;
	
	//new architecture (2014_07_21)
	//ClassMethod()
	//gmatstep aGMATstep(0, 0, 0);
	//GMATDebug << endl;
	//GMATDebug << " ----- class() gmatstep TESTING ----- " << endl;
	//GMATDebug << "j:  " << aGMATstep.j << endl;
	//GMATDebug << "p:  " << aGMATstep.p << endl;
	//GMATDebug << "gs: " << aGMATstep.gs << endl;
	//GMATDebug << "identifier: " << aGMATstep.identifier << endl;
	//GMATDebug << " ------------------------------------ " << endl;


	//loop over all the timesteps used for each journey, each phase
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {

			//reset 'isForwardStep'
			isForwardStep = true;
			//reset 'useInsert'
			usePushBack = true;
			//set 'f2b_index', which indicates where in the vector storage variables
			//items should be inserted when saving backward spacecraft information
			f2b_index = gmat_step_forcemodel_storage.size();
			//initialize counter of gmatsteps
			gmatsteps_thus_far_in_phase = 0;

			GMATDebug << "f2b_index: " << f2b_index << endl;

			for (int s = 0; s < this->ptr_gmatmission->options.num_timesteps; ++s) {
				
				GMATDebug << "Step: " << s;

				//save how many total gmat steps have been created thus far
				gmatsteps_thus_far = gmat_step_timesteps.size();

				//are we looking at a 'Forward' or 'Backward' step 's'
				if (s >= (this->ptr_gmatmission->options.num_timesteps / 2)) { isForwardStep = false; }

				//select the missionbody index to use in analysis
				if (isForwardStep) { body_index = p; }
				else { body_index = p + 1; }

				//is the current step adjacent to the matchpoint
				//if yes, then we will need to add an additional gmat step
				if (s == ((this->ptr_gmatmission->options.num_timesteps / 2) - 1) || s == (this->ptr_gmatmission->options.num_timesteps / 2)) { isStep_MatchPointAdjacent = true; }
				else { isStep_MatchPointAdjacent = false; }

				//additional gmat step for half-step after matchpoint
				if (!isForwardStep && isStep_MatchPointAdjacent) {

					GMATDebug << " MP+ ";

					//find the stepsize
					stepsize = 0.5*this->ptr_gmatmission->journeys[j].phases[p].time_step_sizes[s];
					//advance 'iepoch' and 'fepoch'
					iepoch = fepoch;
					fepoch += stepsize;

					//ClassMethod()
					gmatstep aGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					aGMATstep.setNames(isForwardStep);

					//NOTE: could eliminate if (usePushBack) here because it "should" always be true before this if-statement.
					//we assume that we are NOT exiting NOR entering bodies during the match point constraint, THERFORE,

					//name of the forcemodel
					tempstream << "FM_" << missionbodies[body_index].central_body_name;
					tempstream << "_3rdBodies_" << missionbodies[body_index - 1].name << "_" << missionbodies[body_index].name;
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					aGMATstep.setFMandProp(tempstream.str(), false);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str());
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");

					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					tempstream << missionbodies[body_index - 1].name << ", " << missionbodies[body_index].name;
					forcemodel_storage.push_back(tempstream.str());

					//clear the tempstream
					tempstream.str("");
					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(stepsize); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, stepsize); }
					//ClassMethod()
					aGMATstep.initial_timestep = stepsize;

					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(false); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, false); }
					//store the forcemodel
					if (usePushBack) { gmat_step_forcemodel_storage.push_back(forcemodel_storage); }
					else { gmat_step_forcemodel_storage.insert(gmat_step_forcemodel_storage.begin() + f2b_index, forcemodel_storage); }
					//store the propagator
					if (usePushBack) { gmat_step_propagator_storage.push_back(propagator_storage); }
					else { gmat_step_propagator_storage.insert(gmat_step_propagator_storage.begin() + f2b_index, propagator_storage); }

					//add the thrust vector 
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					gmat_step_thrust_vectors.push_back(thrustvector);
					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					aGMATstep.setThrust(thrustvector);

					//thruster names
					tempstream << "Thruster_j" << j << "p" << p << "_Backward";
					gmat_step_thruster_names.push_back(tempstream.str());
					tempstream.str("");

					
					//ClassMethod()
					gmat_steps.push_back(aGMATstep);


					//clear the temporary forcemodel_storage vector
					forcemodel_storage.clear();
					//clear the temporary propagator_storage vector
					propagator_storage.clear();



					//increment the steps taken thus far
					gmatsteps_thus_far_in_phase++;
					gmatsteps_thus_far++;

					//by definition, we only enter this for-statement if we are NOT on a forward step
					//we need to decrement 'f2b_index' to setup the correct index for using the insert() method for vectors
					f2b_index--;
					//bool flag for whether we should start using inserts instead of push_back
					usePushBack = false;

				}

				//get the current step size, which by definition can be a variable timestep although we always use uniform in practice
				if (s == 0 || s == this->ptr_gmatmission->options.num_timesteps - 1) { 
					stepsize = 0.5*this->ptr_gmatmission->journeys[j].phases[p].time_step_sizes[s]; 
				}
				else { 
					if (isForwardStep) { stepsize = 0.5*this->ptr_gmatmission->journeys[j].phases[p].time_step_sizes[s - 1] + 0.5*this->ptr_gmatmission->journeys[j].phases[p].time_step_sizes[s]; }
					else { stepsize = 0.5*this->ptr_gmatmission->journeys[j].phases[p].time_step_sizes[s] + 0.5*this->ptr_gmatmission->journeys[j].phases[p].time_step_sizes[s + 1]; }
				}

				//get the starting and ending epoch of the current step for analysis
				//note: we always have a half-step at the start and end of the phase as well as just before and after the match point
				if (s == 0) {
					iepoch = this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch;
					fepoch = iepoch + stepsize;
				}
				else {
					iepoch = fepoch;
					fepoch += stepsize;
				}

				//get body's states at 'iepoch' and 'fepoch'
				missionbodies[body_index].locate_body(iepoch, body_istates, false, &this->ptr_gmatmission->options);
				missionbodies[body_index].locate_body(fepoch, body_fstates, false, &this->ptr_gmatmission->options);
				//there is an odd case where we could analyze the 'body_istates' of the first body and 'body_fstates' of the
				//second body across the matchpoint, but this appears unnecessary at this time.

				//compute the state relative difference of the spacecraft and body at 'iepoch'
				if (isForwardStep) {
					if (s == 0) {
						for (int index = 0; index < 3; ++index) {
							initial_position_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[index] - body_istates[index];
							initial_velocity_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[index + 3] - body_istates[index + 3];
						}
					}
					else {
						for (int index = 0; index < 3; ++index) {
							initial_position_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s - 1][index] - body_istates[index];
							initial_velocity_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s - 1][index + 3] - body_istates[index + 3];
						}
					}
				}
				else {
					for (int index = 0; index < 3; ++index) {
						initial_position_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s][index] - body_istates[index];
						initial_velocity_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s][index + 3] - body_istates[index + 3];
					}
				}
				//compute the state relative difference of the spacecraft and body at 'fepoch'
				if (isForwardStep) {
					for (int index = 0; index < 3; ++index) {
						final_position_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s][index] - body_fstates[index];
						final_velocity_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s][index + 3] - body_fstates[index + 3];
					}
				}
				else {
					if (s == (this->ptr_gmatmission->options.num_timesteps - 1)) {
						for (int index = 0; index < 3; ++index) {
							final_position_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[index] - body_fstates[index];
							final_velocity_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[index + 3] - body_fstates[index + 3];
						}
					}
					else {
						for (int index = 0; index < 3; ++index) {
							final_position_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s + 1][index] - body_fstates[index];
							final_velocity_diff[index] = this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s + 1][index + 3] - body_fstates[index + 3];
						}
					}
				}

				//evaluate whether the spacecraft was in the body's SOI at 'iepoch' and 'fepoch'
				if (math::norm(initial_position_diff, 3) < missionbodies[body_index].r_SOI) { isInSOI_at_Start = true; }
				else { isInSOI_at_Start = false; }
				if (math::norm(final_position_diff, 3) < missionbodies[body_index].r_SOI) { isInSOI_at_End = true; }
				else { isInSOI_at_End = false; }

				//is this the step where we exit or enter the SOI?
				//if true, then we need to cut another GMAT step here to use a different Propagator (i.e. FM)
				tempstream << "FM_" << missionbodies[body_index].central_body_name;
				if (isInSOI_at_Start && isInSOI_at_End) {

					GMATDebug << " SOI-SOI ";

					//ClassMethod()
					gmatstep aGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					aGMATstep.setNames(isForwardStep);

					//name of the forcemodel
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					aGMATstep.setFMandProp(tempstream.str(), true);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str() + "_CloseApproach");
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					forcemodel_storage.push_back("");


					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(stepsize); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, stepsize); }
					//ClassMethod()
					aGMATstep.initial_timestep = stepsize;

					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(true); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, true); }

					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					aGMATstep.setThrust(thrustvector);
					//ClassMethod()
					if (usePushBack) { gmat_steps.push_back(aGMATstep); }
					else { gmat_steps.insert(gmat_steps.begin() + f2b_index, aGMATstep); }

					//increment the 'f2b_index' variable if we are still on a forward step
					if (isForwardStep) { f2b_index++; }

					//increment the variable counting gmatsteps in the current phase
					gmatsteps_thus_far_in_phase++;

				}
				else if (isInSOI_at_Start && !isInSOI_at_End) {

					GMATDebug << " SOI-> ";

					//ClassMethod()
					gmatstep aGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					aGMATstep.setNames(isForwardStep);

					//for this case we must figure out how long we are in the SOI before exiting
					periapse_velocity_magnitude = sqrt(2.0 * missionbodies[body_index].mu / (missionbodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + math::norm(initial_velocity_diff, 3)*math::norm(initial_velocity_diff, 3));
					approximate_time_in_SOI = (missionbodies[body_index].r_SOI / periapse_velocity_magnitude);

					//approximate_time_in_SOI should be shorter than stepsize
					
					//Creation of the First Model
					//name of the forcemodel
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					aGMATstep.setFMandProp(tempstream.str(), true);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str() + "_CloseApproach");
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					forcemodel_storage.push_back("");
					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(approximate_time_in_SOI); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, approximate_time_in_SOI); }
					//ClassMethod()
					aGMATstep.initial_timestep = approximate_time_in_SOI;

					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(true); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, true); }
					//store the forcemodel
					if (usePushBack) { gmat_step_forcemodel_storage.push_back(forcemodel_storage); }
					else { gmat_step_forcemodel_storage.insert(gmat_step_forcemodel_storage.begin() + f2b_index, forcemodel_storage); }
					//clear the temporary forcemodel_storage vector
					forcemodel_storage.clear();
					//store the propagator
					if (usePushBack) { gmat_step_propagator_storage.push_back(propagator_storage); }
					else { gmat_step_propagator_storage.insert(gmat_step_propagator_storage.begin() + f2b_index, propagator_storage); }
					//clear the temporary propagator_storage vector
					propagator_storage.clear();

					//increment the variable counting gmatsteps in the current phase
					gmatsteps_thus_far_in_phase++;

					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					aGMATstep.setThrust(thrustvector);
					//ClassMethod()
					if (usePushBack) { gmat_steps.push_back(aGMATstep); }
					else { gmat_steps.insert(gmat_steps.begin() + f2b_index, aGMATstep); }
					//ClassMethod()
					gmatstep bGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					bGMATstep.setNames(isForwardStep);


					//Creation of the Second Model
					//name of the forcemodel
					tempstream << "FM_" << missionbodies[body_index].central_body_name;
					if (isForwardStep) { tempstream << "_3rdBodies_" << missionbodies[body_index].name << "_" << missionbodies[body_index + 1].name; }
					else { tempstream << "_3rdBodies_" << missionbodies[body_index - 1].name << "_" << missionbodies[body_index].name; }
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					bGMATstep.setFMandProp(tempstream.str(), false);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str());
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					if (isForwardStep) { tempstream << missionbodies[body_index].name << ", " << missionbodies[body_index + 1].name; }
					else { tempstream << missionbodies[body_index - 1].name << ", " << missionbodies[body_index].name; }
					forcemodel_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(stepsize - approximate_time_in_SOI); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, stepsize - approximate_time_in_SOI); }
					//ClassMethod()
					bGMATstep.initial_timestep = stepsize - approximate_time_in_SOI;


					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(false); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, false); }

					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					bGMATstep.setThrust(thrustvector);
					//ClassMethod()
					if (usePushBack) { gmat_steps.push_back(bGMATstep); }
					else { gmat_steps.insert(gmat_steps.begin() + f2b_index, bGMATstep); }

					//increment the 'f2b_index' variable if we are still on a forward step
					if (isForwardStep) { f2b_index += 2; }

					//increment the variable counting gmatsteps in the current phase
					gmatsteps_thus_far_in_phase++;

				}
				else if (!isInSOI_at_Start && isInSOI_at_End) {

					GMATDebug << " ->SOI ";

					//ClassMethod()
					gmatstep aGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					aGMATstep.setNames(isForwardStep);

					//for this case we must figure out how long we are in the SOI after entering
					periapse_velocity_magnitude = sqrt(2.0 * missionbodies[body_index].mu / (missionbodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + math::norm(final_velocity_diff, 3)*math::norm(final_velocity_diff, 3));
					approximate_time_in_SOI = (missionbodies[body_index].r_SOI / periapse_velocity_magnitude);

					//approximate_time_in_SOI should be shorter than stepsize

					//Creation of the First Model
					//name of the forcemodel
					if (isForwardStep) { tempstream << "_3rdBodies_" << missionbodies[body_index].name << "_" << missionbodies[body_index + 1].name; }
					else { tempstream << "_3rdBodies_" << missionbodies[body_index - 1].name << "_" << missionbodies[body_index].name; }
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					aGMATstep.setFMandProp(tempstream.str(), false);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str());
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					if (isForwardStep) { tempstream << missionbodies[body_index].name << ", " << missionbodies[body_index + 1].name; }
					else { tempstream << missionbodies[body_index - 1].name << ", " << missionbodies[body_index].name; }
					forcemodel_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(stepsize - approximate_time_in_SOI); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, stepsize - approximate_time_in_SOI); }
					//ClassMethod()
					aGMATstep.initial_timestep = stepsize - approximate_time_in_SOI;

					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(false); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, false); }
					//store the forcemodel
					if (usePushBack) { gmat_step_forcemodel_storage.push_back(forcemodel_storage); }
					else { gmat_step_forcemodel_storage.insert(gmat_step_forcemodel_storage.begin() + f2b_index, forcemodel_storage); }
					//clear the temporary forcemodel_storage vector
					forcemodel_storage.clear();
					//store the propagator
					if (usePushBack) { gmat_step_propagator_storage.push_back(propagator_storage); }
					else { gmat_step_propagator_storage.insert(gmat_step_propagator_storage.begin() + f2b_index, propagator_storage); }
					//clear the temporary propagator_storage vector
					propagator_storage.clear();

					//increment the variable counting gmatsteps in the current phase
					gmatsteps_thus_far_in_phase++;


					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					aGMATstep.setThrust(thrustvector);
					//ClassMethod()
					if (usePushBack) { gmat_steps.push_back(aGMATstep); }
					else { gmat_steps.insert(gmat_steps.begin() + f2b_index, aGMATstep); }
					//ClassMethod()
					gmatstep bGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					bGMATstep.setNames(isForwardStep);


					//Creation of the Second Model
					//name of the forcemodel
					tempstream << "FM_" << missionbodies[body_index].central_body_name;
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					bGMATstep.setFMandProp(tempstream.str(), true);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str() + "_CloseApproach");
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					forcemodel_storage.push_back("");
					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(approximate_time_in_SOI); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, approximate_time_in_SOI); }
					//ClassMethod()
					aGMATstep.initial_timestep = approximate_time_in_SOI;


					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(true); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, true); }

					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					bGMATstep.setThrust(thrustvector);
					//ClassMethod()
					if (usePushBack) { gmat_steps.push_back(bGMATstep); }
					else { gmat_steps.insert(gmat_steps.begin() + f2b_index, bGMATstep); }

					//increment the 'f2b_index' variable if we are still on a forward step
					if (isForwardStep) { f2b_index += 2; }

					//increment the variable counting gmatsteps in the current phase
					gmatsteps_thus_far_in_phase++;

				}
				else {

					GMATDebug << " ---> ";

					//ClassMethod()
					gmatstep aGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					aGMATstep.setNames(isForwardStep);

					//name of the forcemodel
					if (isForwardStep) { tempstream << "_3rdBodies_" << missionbodies[body_index].name     << "_" << missionbodies[body_index + 1].name; }
					else {               tempstream << "_3rdBodies_" << missionbodies[body_index - 1].name << "_" << missionbodies[body_index].name; }
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					aGMATstep.setFMandProp(tempstream.str(), false);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str());
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					if (isForwardStep) { tempstream << missionbodies[body_index].name     << ", " << missionbodies[body_index + 1].name; }
					else {               tempstream << missionbodies[body_index - 1].name << ", " << missionbodies[body_index].name; }
					forcemodel_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(stepsize); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, stepsize); }
					//ClassMethod()
					aGMATstep.initial_timestep = stepsize;


					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(false); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, false); }
					
					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					aGMATstep.setThrust(thrustvector);
					//ClassMethod()
					if (usePushBack) { gmat_steps.push_back(aGMATstep); }
					else { gmat_steps.insert(gmat_steps.begin() + f2b_index, aGMATstep); }

					//increment the 'f2b_index' variable if we are still on a forward step
					if (isForwardStep) { f2b_index++; }

					//increment the variable counting gmatsteps in the current phase
					gmatsteps_thus_far_in_phase++;

				}

				//store the forcemodel
				if (usePushBack) { gmat_step_forcemodel_storage.push_back(forcemodel_storage); }
				else { gmat_step_forcemodel_storage.insert(gmat_step_forcemodel_storage.begin() + f2b_index, forcemodel_storage); }
				//clear the temporary forcemodel_storage vector
				forcemodel_storage.clear();
				//store the propagator
				if (usePushBack) { gmat_step_propagator_storage.push_back(propagator_storage); }
				else { gmat_step_propagator_storage.insert(gmat_step_propagator_storage.begin() + f2b_index, propagator_storage); }
				//clear the temporary propagator_storage vector
				propagator_storage.clear();

				//additional gmat step for half-step before matchpoint
				if (isForwardStep && isStep_MatchPointAdjacent) {

					GMATDebug << " -MP ";

					//ClassMethod()
					gmatstep aGMATstep(j, p, gmatsteps_thus_far_in_phase, missionbodies);
					aGMATstep.setNames(isForwardStep);

					//find the stepsize
					stepsize = 0.5*this->ptr_gmatmission->journeys[j].phases[p].time_step_sizes[s];
					//advance 'iepoch' and 'fepoch'
					iepoch = fepoch;
					fepoch += stepsize;
					
					//we assume that we are NOT exiting NOR entering bodies during the match point constraint, THERFORE,

					//name of the forcemodel
					tempstream << "FM_" << missionbodies[body_index].central_body_name;
					tempstream << "_3rdBodies_" << missionbodies[body_index].name << "_" << missionbodies[body_index + 1].name;
					forcemodel_storage.push_back(tempstream.str());
					//ClassMethod()
					aGMATstep.setFMandProp(tempstream.str(), false);

					//name of the propagator
					propagator_storage.push_back("Propagator_" + tempstream.str());
					propagator_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//the forcemodel central body
					forcemodel_storage.push_back(missionbodies[body_index].central_body_name);
					//the forcemodel 3rd body point masses
					tempstream << missionbodies[body_index].name << ", " << missionbodies[body_index + 1].name;
					forcemodel_storage.push_back(tempstream.str());
					//clear the tempstream
					tempstream.str("");
					//save the timestep length
					if (usePushBack) { gmat_step_timesteps.push_back(stepsize); }
					else { gmat_step_timesteps.insert(gmat_step_timesteps.begin() + f2b_index, stepsize); }
					//ClassMethod()
					aGMATstep.initial_timestep = stepsize;


					//save whether we are conducting a close approach
					if (usePushBack) { gmat_step_propagator_isCloseApproach.push_back(false); }
					else { gmat_step_propagator_isCloseApproach.insert(gmat_step_propagator_isCloseApproach.begin() + f2b_index, false); }

					//store the forcemodel
					if (usePushBack) { gmat_step_forcemodel_storage.push_back(forcemodel_storage); }
					else { gmat_step_forcemodel_storage.insert(gmat_step_forcemodel_storage.begin() + f2b_index, forcemodel_storage); }
					//clear the temporary forcemodel_storage vector
					forcemodel_storage.clear();
					//store the propagator
					if (usePushBack) { gmat_step_propagator_storage.push_back(propagator_storage); }
					else { gmat_step_propagator_storage.insert(gmat_step_propagator_storage.begin() + f2b_index, propagator_storage); }
					//clear the temporary propagator_storage vector
					propagator_storage.clear();

					//ClassMethod()
					aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
					aGMATstep.setThrust(thrustvector);
					//ClassMethod()
					if (usePushBack) { gmat_steps.push_back(aGMATstep); }
					else { gmat_steps.insert(gmat_steps.begin() + f2b_index, aGMATstep); }

					//increment the 'f2b_index' variable if we are still on a forward step
					if (isForwardStep) { f2b_index += 2; }

					//increment the variable counting gmatsteps in the current phase
					gmatsteps_thus_far_in_phase++;

				}



				//store the thrust vector initialize information for the gmat steps created above
				number_of_gmatsteps_added = gmat_step_timesteps.size() - gmatsteps_thus_far;
				aux_GMAT_populate_thrustvector(j, p, s, thrustvector, -1.0, 1.0);
				for (int index = 0; index < number_of_gmatsteps_added; ++index) {
					if (usePushBack) { gmat_step_thrust_vectors.push_back(thrustvector); }
					else { gmat_step_thrust_vectors.insert(gmat_step_thrust_vectors.begin() + f2b_index, thrustvector); }
				}
				thrustvector.clear();

				//store the thruster names for the gmat steps
				if (isForwardStep) {
					tempstream << "Thruster_j" << j << "p" << p << "_Forward";
					for (int index = 0; index < number_of_gmatsteps_added; ++index) {
						if (usePushBack) { gmat_step_thruster_names.push_back(tempstream.str()); }
						else { gmat_step_thruster_names.insert(gmat_step_thruster_names.begin() + f2b_index, tempstream.str()); }
					}
				}
				else {
					tempstream << "Thruster_j" << j << "p" << p << "_Backward";
					for (int index = 0; index < number_of_gmatsteps_added; ++index) {
						if (usePushBack) { gmat_step_thruster_names.push_back(tempstream.str()); }
						else { gmat_step_thruster_names.insert(gmat_step_thruster_names.begin() + f2b_index, tempstream.str()); }
					}
				}
				tempstream.str("");



				//on the last step, figure out how many gmat steps have been created for this phase
				if (s == (this->ptr_gmatmission->options.num_timesteps - 1)) {
					if (p == 0) { gmat_steps_per_phase.push_back(gmat_step_timesteps.size()); }
					else { gmat_steps_per_phase.push_back(gmat_step_timesteps.size() - gmat_steps_per_phase.back()); }
				}



				GMATDebug << endl;
				GMATDebug << endl;

				//Debugging 
				GMATDebug << "options.num_timesteps: " << this->ptr_gmatmission->options.num_timesteps << endl;
				GMATDebug << "timestep number: " << s << endl;
				GMATDebug << "-----------------------" << endl;
				GMATDebug << this->ptr_gmatmission->journeys[0].phases[0].spacecraft_state.size() << endl;
				for (int states = 0; states < this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s].size(); ++states) {
					GMATDebug << this->ptr_gmatmission->journeys[j].phases[p].spacecraft_state[s][states] << endl;
				}
				GMATDebug << "-----------------------" << endl;

				// for-statement for each step
				// save information of forcemodeling etc.. into a structure with an index counter 
				// that indicates upto which step some value is appropriate

			}//end step for-statement
		}//end phase for-statement
	}//end journey for-statement

	//debug check to make sure all the information was saved correctly
	int number_of_total_steps = 0;
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {
			number_of_total_steps += this->ptr_gmatmission->options.num_timesteps;
		}
	}
	GMATDebug << endl;
	GMATDebug << "/////////////////////////////////////////" << endl;
	GMATDebug << "number_of_total_steps: " << number_of_total_steps << endl;
	GMATDebug << "gmat_step_forcemodel_storage.size(): " << gmat_step_forcemodel_storage.size() << endl;
	GMATDebug << "gmat_step_propagator_storage.size(): " << gmat_step_propagator_storage.size() << endl;
	GMATDebug << "gmat_step_propagator_isCloseApproach.size(): " << gmat_step_propagator_isCloseApproach.size() << endl;
	GMATDebug << "gmat_step_timesteps.size(): " << gmat_step_timesteps.size() << endl;
	GMATDebug << "gmat_step_thrust_vectors.size(): " << gmat_step_thrust_vectors.size() << endl;
	GMATDebug << "gmat_step_thruster_name.size(): " << gmat_step_thruster_names.size() << endl;
	GMATDebug << "gmat_steps.size(): " << gmat_steps.size() << endl;
	for (int index = 0; index < gmat_step_forcemodel_storage.size(); ++index) {
		GMATDebug << endl;
		GMATDebug << " -----  GMAT STEP: " << index << " ----- " << endl;
		GMATDebug << endl;
		// Print ForceModel Info
		GMATDebug << "ForceModel: ";
		for (int i = 0; i < 3; ++i) {
			GMATDebug << gmat_step_forcemodel_storage[index][i] << "   ";
		}
		GMATDebug << endl;
		// Print Propagator Info
		GMATDebug << "Propagator: ";
		for (int i = 0; i < 2; ++i) {
			GMATDebug << gmat_step_propagator_storage[index][i] << "   ";
		}
		GMATDebug << gmat_step_propagator_isCloseApproach[index] << endl;
		// Print Control Info
		GMATDebug << "Control: ";
		for (int i = 0; i < 3; ++i) {
			GMATDebug << gmat_step_thrust_vectors[index][i] << "   ";
		}
		GMATDebug << endl;
		// Print the gmat timestep
		GMATDebug << "Timestep: " << gmat_step_timesteps[index] << endl;
		// Print the gmat thruster name
		GMATDebug << "Thruster Name: " << gmat_step_thruster_names[index] << endl;
	}
	GMATDebug << "/////////////////////////////////////////" << endl;

}//end method


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


// method to create hardware information
void gmatscripter::write_GMAT_hardware(){

	//declarations
	stringstream fueltankname_stream;
	stringstream thrustername_stream;
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
			//create a fuel tank for the forward propagating spacecraft of the phase
			fueltankname_stream << "FuelTank_j" << j << "p" << p << "_Forward";
			fueltank_forward_names.push_back(fueltankname_stream.str());
			this->create_GMAT_fueltank(j, p, fueltank_forward_names.back(), "Forward");
			//create a thruster for the foward propagating spacecraft of the phase
			thrustername_stream << "Thruster_j" << j << "p" << p << "_Forward";
			thruster_forward_names.push_back(thrustername_stream.str());
			this->create_GMAT_thruster(j, p, thruster_forward_names.back(), missionbodies[body_index].name, fueltank_forward_names.back());
			//clear the stringstreams
			fueltankname_stream.str("");
			thrustername_stream.str("");

			//create some whitespace
			GMATfile << endl;

			//create a fuel tank for the backward propagating spacecraft of the phase
			fueltankname_stream << "FuelTank_j" << j << "p" << p << "_Backward";
			fueltank_backward_names.push_back(fueltankname_stream.str());
			this->create_GMAT_fueltank(j, p, fueltank_backward_names.back(), "Backward");
			//create a thruster for the foward propagating spacecraft of the phase
			thrustername_stream << "Thruster_j" << j << "p" << p << "_Backward";
			thruster_backward_names.push_back(thrustername_stream.str());
			this->create_GMAT_thruster(j, p, thruster_backward_names.back(), missionbodies[body_index].name, fueltank_backward_names.back());
			//clear the stringstreams
			fueltankname_stream.str("");
			thrustername_stream.str("");

			//create some whitespace
			GMATfile << endl;
			//increment
			++body_index;
		}//end of journeys for-statement
	}//end of phases for-statement
	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_hardware() method


// method to create spacecraft information
void gmatscripter::write_GMAT_spacecraft(){

	//declarations
	stringstream spacecraftname_stream;
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
			spacecraftname_stream << "SpaceCraft_j" << j << "p" << p << "_Forward";
			//add the name to the vector 'spacecraft_names'
			spacecraft_names.push_back(spacecraftname_stream.str());
			spacecraft_forward_names.push_back(spacecraftname_stream.str());
			//create the GMAT spacecraft
			this->create_GMAT_spacecraft(j, p, spacecraft_forward_names.back(), "Forward", missionbodies[body_index].name);
			//clear the stringstream
			spacecraftname_stream.str("");

			//create a name for the backward propagating spacecraft
			spacecraftname_stream << "SpaceCraft_j" << j << "p" << p << "_Backward";
			spacecraft_names.push_back(spacecraftname_stream.str());
			spacecraft_backward_names.push_back(spacecraftname_stream.str());
			//create the GMAT spacecraft
			this->create_GMAT_spacecraft(j, p, spacecraft_backward_names.back(), "Backward", missionbodies[body_index + 1].name);
			//clear the stringstream
			spacecraftname_stream.str("");

			//increment
			++body_index;

			//spacing
			GMATfile << endl;
		}// end of phases for-statement
	}//end of journeys for-statement

}//end of write_GMAT_spacecraft() method


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

	//declaration
	vector <string> tempstrings;
	bool hasNotBeenCopied;

	//force model header 
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Force Models" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//save the first force model and copy its name to a temporary vector of strings
	//we will check against the temporary vector of strings to make sure we are not 
	//printing duplicates
	this->create_GMAT_forcemodel(gmat_step_forcemodel_storage[0][0],
								 gmat_step_forcemodel_storage[0][1],
								 gmat_step_forcemodel_storage[0][2]);
	tempstrings.push_back(gmat_step_forcemodel_storage[0][0]);

	//generate the rest of the forcemodels
	for (int i = 1; i < gmat_step_forcemodel_storage.size(); ++i) {
		//set bool to true
		hasNotBeenCopied = true;
		//loop through any saved forcemodel names
		for (int j = 0; j < tempstrings.size(); ++j) {
			//if we have already saved out the current forcemodel, 
			//then turn the bool flag to false and do not print the current forcemodel
			if (tempstrings[j].compare(gmat_step_forcemodel_storage[i][0]) == 0) { hasNotBeenCopied = false; }
		}
		//if it hasn't been printed yet, then print it!
		if (hasNotBeenCopied) {
			this->create_GMAT_forcemodel(gmat_step_forcemodel_storage[i][0],
				gmat_step_forcemodel_storage[i][1],
				gmat_step_forcemodel_storage[i][2]);
			tempstrings.push_back(gmat_step_forcemodel_storage[i][0]);
		}
	}

	//create central force model (TODO:: assume for now that first body's central body is the principle central body)
	//GMATfile << "%Create 2-body force model for central body" << endl;
	//this->create_GMAT_forcemodel(missionbodies_unique[0].central_body_name + "2Bod", 
	//						     missionbodies_unique[0].central_body_name, 
	//						     missionbodies_unique[0].central_body_name);

	//GMATfile << "%Create 2-body force models for all bodies visited" << endl;
	//for (int body_index = 0; body_index < missionbodies_unique.size(); body_index++) {
	//	this->create_GMAT_forcemodel(missionbodies_unique[body_index].name + "2Bod",
	//							     missionbodies_unique[body_index].name,
	//							     missionbodies_unique[body_index].name);
	//}

	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;
}//end of write_GMAT_forcemodels() method


// method to create propagator information
void gmatscripter::write_GMAT_propagators(){

	//declaration
	vector <string> tempstrings;
	bool hasNotBeenCopied;

	//propagator header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Propagators" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//save the first propagator and copy its name to a temporary vector of strings
	//we will check against the temporary vector of strings to make sure we are not 
	//printing duplicates
	this->create_GMAT_propagator(gmat_step_propagator_storage[0][0],
								 gmat_step_propagator_storage[0][1],
								 gmat_step_propagator_isCloseApproach[0]);
	tempstrings.push_back(gmat_step_propagator_storage[0][0]);

	//generate the rest of the propagators
	for (int i = 1; i < gmat_step_propagator_storage.size(); ++i) {
		//set bool to true
		hasNotBeenCopied = true;
		//loop through any saved propagator names
		for (int j = 0; j < tempstrings.size(); ++j) {
			//if we have already saved out the current propagator, 
			//then turn the bool flag to false and do not print the current propagator
			if (tempstrings[j].compare(gmat_step_propagator_storage[i][0]) == 0) { hasNotBeenCopied = false; }
		}

		//if it hasn't been printed yet, then print it!
		if (hasNotBeenCopied) {
			this->create_GMAT_propagator(gmat_step_propagator_storage[i][0],
										 gmat_step_propagator_storage[i][1],
										 gmat_step_propagator_isCloseApproach[i]);
			tempstrings.push_back(gmat_step_propagator_storage[i][0]);
		}
	}



	////far planet propagator, larger time steps
	////TODO:: for now it is assumed that the principle central body is the first central body
	//GMATfile << "%Create propagation model for central body" << endl;
	//this->create_GMAT_propagator(missionbodies_unique[0].central_body_name + "Prop", 
	//						     missionbodies_unique[0].central_body_name + "2Bod", 
	//						     false);

	////create propagators for close body approaches, smaller time steps
	//GMATfile << "%Create propagation models for other bodies" << endl;
	//for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index) {
	//	this->create_GMAT_propagator(missionbodies_unique[body_index].name + "Prop",
	//							     missionbodies_unique[body_index].name + "2Bod",
	//								 true);
	//}



	//add some vertical whitespace
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_propagators() method


// method to create burn information
void gmatscripter::write_GMAT_burns(){

	//declaration
	stringstream finiteburnname_stream;

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
			GMATfile << "% Journey #" << j << ", Phase #" << p << ", Forward and Backward FiniteBurn Objects" << endl;
			finiteburnname_stream << "FiniteBurn_j" << j << "p" << p << "_Forward";
			finiteburn_forward_names.push_back(finiteburnname_stream.str());
			this->create_GMAT_finiteburn(finiteburn_forward_names.back(), thruster_forward_names[p]);
			finiteburnname_stream.str("");
			finiteburnname_stream << "FiniteBurn_j" << j << "p" << p << "_Backward";
			finiteburn_backward_names.push_back(finiteburnname_stream.str());
			this->create_GMAT_finiteburn(finiteburn_backward_names.back(), thruster_backward_names[p]);
			//create some whitespace
			GMATfile << endl;



			//GMATfile << "Create FiniteBurn FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
			//GMATfile << "FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward.Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward};" << endl;
			//GMATfile << "Create FiniteBurn FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
			//GMATfile << "FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward.Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward};" << endl;
			//GMATfile << endl;
		}//end of journeys for-statement
	}//end of phases for-statement
	GMATfile << endl;

}//end of write_GMAT_burns() method


// method to create array and variable information
void gmatscripter::write_GMAT_variables(){

	//arrays and variables header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Arrays, Variables, Strings" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//write out mission level variables
	GMATfile << "% --- Mission Level: Arrays, Variables, and Strings" << endl;
	for (int index = 0; index < mission_level_variables.size(); ++index) {
		GMATfile << "Create Variable " << mission_level_variables[index][0] << endl;
		if (mission_level_variables[index].size() == 2) {
			GMATfile << mission_level_variables[index][0] << mission_level_variables[index][1] << endl;
		}
	}

	//write out journey level variables
	GMATfile << "% --- Journey Level: Arrays, Variables, and Strings" << endl;
	for (int index = 0; index < journey_level_variables.size(); ++index) {
		GMATfile << "Create Variable " << journey_level_variables[index][0] << endl;
		if (journey_level_variables[index].size() == 2) {
			GMATfile << journey_level_variables[index][0] << journey_level_variables[index][1] << endl;
		}
	}

	//write out phase level variables
	GMATfile << "% --- Phase Level: Arrays, Variables, and Strings" << endl;
	for (int index = 0; index < phase_level_variables.size(); ++index) {
		GMATfile << "Create Variable " << phase_level_variables[index][0] << endl;
		if (phase_level_variables[index].size() == 2) {
			GMATfile << phase_level_variables[index][0] << phase_level_variables[index][1] << endl;
		}
	}
	//create arrays for spacecraft thrust vectors
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {
			GMATfile << "Create Array ThrustVector_j" << j << "p" << p << "[" << gmat_steps_per_phase[p] << ", 3];" << endl;
			GMATfile << "Create Array ThrustUnitVectorMagnitude_j" << j << "p" << p << "[" << gmat_steps_per_phase[p] << ", 1];" << endl;
		}
	}
	GMATfile << endl;

	//create time related variables

	//write out spacecraft variables
	//create scaled states for each spacecraft
	for (int index = 0; index < spacecraft_names.size(); ++index) {
		GMATfile << "%Scaled States of " << spacecraft_names[index] << endl;
		GMATfile << "Create Variable " << spacecraft_names[index] << "_X_Scaled" << endl;
		GMATfile << "Create Variable " << spacecraft_names[index] << "_Y_Scaled" << endl;
		GMATfile << "Create Variable " << spacecraft_names[index] << "_Z_Scaled" << endl;
		GMATfile << "Create Variable " << spacecraft_names[index] << "_VX_Scaled" << endl;
		GMATfile << "Create Variable " << spacecraft_names[index] << "_VY_Scaled" << endl;
		GMATfile << "Create Variable " << spacecraft_names[index] << "_VZ_Scaled" << endl;
		GMATfile << "Create Variable " << spacecraft_names[index] << "_FuelMass_Scaled" << endl;
		GMATfile << endl;
	}
	GMATfile << endl;


	////write out variables without a current home
	//GMATfile << "% --- Homeless: Arrays, Variables, and Strings" << endl;
	//GMATfile << "Create Variable FinalEpoch" << endl;
	//GMATfile << "Create Variable LaunchEpoch_Scaled" << endl;
	////create variable for each s/c create rdotv temp variable
	//for (int index = 0; index < spacecraft_names.size(); ++index) {
	//	GMATfile << "Create Variable " << spacecraft_names[index] << "_RdotV_Scaled " << endl;
	//	GMATfile << "Create Variable " << spacecraft_names[index] << "_PeriapseRadius_Scaled" << endl;
	//}
	//GMATfile << endl;


	//write out report variables
	//create temporary strings to identify time steps during phases
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {
			for (int gs = 0; gs < gmat_steps_per_phase[p]; gs++){
				GMATfile << "Create String tempString_j" << j << "p" << p << "gs" << gs << "_tminus" << endl;
				GMATfile << "tempString_j" << j << "p" << p << "gs" << gs << "_tminus = 'j" << j << "p" << p << "gs" << gs << "_tminus: '" << endl;
				GMATfile << "Create String tempString_j" << j << "p" << p << "gs" << gs << "_tplus" << endl;
				GMATfile << "tempString_j" << j << "p" << p << "gs" << gs << "_tplus = 'j" << j << "p" << p << "gs" << gs << "_tplus: '" << endl;
			}
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
	//TODO:: as of now it is assumed that the first central body is the principle central body
	GMATfile << "%Create coordinate systems for plotting/viewing" << endl;
	this->create_GMAT_coordinatesystem(missionbodies_unique[0].central_body_name);

	for (int b = 0; b < missionbodies_unique.size(); ++b) {
		this->create_GMAT_coordinatesystem(missionbodies_unique[b].name);
	}
	//add some vertical whitespace
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
	GMATfile << "NLPObject.MaximumIterations = 5;" << endl;
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
	double boundary_state[6];
	math::Matrix<double> Vinf_in(3, 1);
	math::Matrix<double> Vinf_out(3, 1);
	math::Matrix<double> periapse_state_vector(6, 1);
	math::Matrix<double> periapse_position_vector(3, 1);
	math::Matrix<double> periapse_velocity_vector(3, 1);
	int body_index = 0;
	int name_index = 0;

	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Initial State Guesses" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	GMATfile << "BeginScript 'Initial Guess Values' " << endl;
	GMATfile << endl;

	// write out the initial guess for the spacecraft engine parameters
	for (int index = 0; index < spacecraft_names.size(); ++index) {
		if (!isLT) {
			GMATfile << "   " << spacecraft_names[index] << ".K1 = ThrusterISP" << endl;
		}
		else if (isLT && this->ptr_gmatmission->options.engine_type == 0) {
			GMATfile << "   " << spacecraft_names[index] << ".DutyCycle = ThrusterDutyCycle" << endl;
			GMATfile << "   " << spacecraft_names[index] << ".C1        = ThrusterMaxThrust" << endl;
			GMATfile << "   " << spacecraft_names[index] << ".K1        = ThrusterISP" << endl;
			GMATfile << endl;
		}
	}

	//write out the initial guess values for the thrust vector history
	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {
			if (p == 0) {
				for (int index = 0; index < gmat_steps_per_phase[p]; ++index) {
					for (int subindex = 0; subindex < 3; ++subindex) { GMATfile << "	ThrustVector_j" << j << "p" << p << "(" << index << ", " << subindex << ") = " << gmat_step_thrust_vectors[index][subindex] << endl; }
				}
			}
			else {
				for (int index = gmat_steps_per_phase[p - 1] - 1; index < gmat_steps_per_phase[p]; ++index) {
					for (int subindex = 0; subindex < 3; ++subindex) { GMATfile << "	ThrustVector_j" << j << "p" << p << "(" << index << ", " << subindex << ") = " << gmat_step_thrust_vectors[index][subindex] << endl; }
				}
			}
			GMATfile << endl;
		}//end of phase for-statement
	}//end of journey for-statement

	GMATfile << "EndScript" << endl;
	GMATfile << endl;

	////TODO:: hook up engine inputs
	//GMATfile << "	%Engine model parameters" << endl;
	//if (this->ptr_gmatmission->options.mission_type < 2 || this->ptr_gmatmission->options.mission_type == 5) //impulsive-thrust phase types
	//{
	//	GMATfile << "	ThrusterISP = " << this->ptr_gmatmission->options.IspChem << endl;

	//	GMATfile << endl;

	//	//for each journey
	//	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	//	{
	//		//for each phase
	//		for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
	//		{
	//			GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".K1 = ThrusterISP;" << endl;
	//			GMATfile << endl;
	//		}
	//	}
	//}
	//else //for all low-thrust phase types
	//{
	//	switch (this->ptr_gmatmission->options.engine_type)
	//	{
	//	case 0: //fixed thrust and ISP
	//		GMATfile << "	ThrusterISP = "       << this->ptr_gmatmission->options.IspLT             << endl;
	//		GMATfile << "	ThrusterMaxThrust = " << this->ptr_gmatmission->options.Thrust            << endl;
	//		GMATfile << "	ThrusterDutyCycle = " << this->ptr_gmatmission->options.engine_duty_cycle << endl;
	//		break;
	//	case 1: //constant ISP, efficiency, EMTG chooses power
	//		break;
	//	case 2: //choice of power model, constant efficiency, EMTG chooses ISP
	//		break;
	//	case 3: //choice of power model, constant efficiency and ISP
	//		break;
	//	case 4: //continuously-varying specific impulse (constant efficiency)
	//		break;
	//	case 5: //custom thrust and mass flow rate
	//		break;
	//	case 6: //NSTAR
	//		break;
	//	case 7: //XIPS-25
	//		break;
	//	case 8: //BPT-4000 high ISP
	//		break;
	//	case 9: //BPT-4000 high thrust
	//		break;
	//	case 10: //BPT-4000 ex-high ISP
	//		break;
	//	case 11: //NEXT high ISP
	//		break;
	//	case 12: //VASIMR
	//		break;
	//	case 13: //Hall thruster
	//		break;
	//	case 14: //NEXT v10 high Isp
	//		break;
	//	case 15: //NEXT v10 high Thrust
	//		break;
	//	case 16: //BPT-4000 MALTO curve
	//		break;
	//	}

	//	GMATfile << endl;

	//	//for each forward spacecraft, set the dutycycle, maxthrust and Isp
	//	for (int index = 0; index < thruster_forward_names.size(); ++index) {
	//		GMATfile << "	" << thruster_forward_names[index] << ".DutyCycle  = ThrusterDutyCycle;" << endl;
	//		GMATfile << "	" << thruster_forward_names[index] << ".C1         = ThrusterMaxThrust;" << endl;
	//		GMATfile << "	" << thruster_forward_names[index] << ".K1         = ThrusterISP;" << endl;
	//		GMATfile << endl;
	//	}
	//	//for each backward spacecraft, set the dutycycle, maxthrust and Isp
	//	for (int index = 0; index < thruster_backward_names.size(); ++index) {
	//		GMATfile << "	" << thruster_backward_names[index] << ".DutyCycle  = ThrusterDutyCycle;" << endl;
	//		GMATfile << "	" << thruster_backward_names[index] << ".C1         = ThrusterMaxThrust;" << endl;
	//		GMATfile << "	" << thruster_backward_names[index] << ".K1         = ThrusterISP;" << endl;
	//		GMATfile << endl;
	//	}

	//	//for each journey
	//	//for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	//	//{
	//	//	//for each phase
	//	//	for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
	//	//	{
	//	//		GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.DutyCycle  = ThrusterDutyCycle;"	<< endl;
	//	//		GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1         = ThrusterMaxThrust;" << endl;
	//	//		GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.K1         = ThrusterISP;"	    << endl;
	//	//		GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.DutyCycle = ThrusterDutyCycle;" << endl;
	//	//		GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.C1        = ThrusterMaxThrust;" << endl;
	//	//		GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.K1        = ThrusterISP;"	    << endl;
	//	//		GMATfile << endl;
	//	//	}
	//	//}

	//} //end code for low-thrust models


	


	////for each journey
	//for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j)
	//{
	//	//for each phase
	//	for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p)
	//	{
	//		//FORWARD PROPAGATED SPACECRAFT
	//		GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Forward s/c initial conditions" << endl;
	//		GMATDebug << " I should be printing right now." << endl;
	//		//interpret state from beginning of each journey based on departure type (//TODO::)
	//		if (p == 0)
	//		{
	//			GMATDebug << "p: " << p << endl;
	//			//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
	//			missionbodies[body_index].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);

	//			switch (this->ptr_gmatmission->options.journey_departure_type[j])
	//			{
	//			case 0: //launch or direct insertion // (EMTG default mission)
	//				GMATDebug << "case: 0" << endl;
	//				//calculate v_infinity vectors
	//				for (int k = 0; k < 3; ++k)
	//				{
	//					Vinf_out(k) = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
	//				}
	//				//calculate inc from vinfinity, then make guess at min altitude at planet
	//				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[0] - boundary_state[0] << ";" << endl;
	//				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[1] - boundary_state[1] << ";" << endl;
	//				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[2] - boundary_state[2] << ";" << endl;
	//				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[3] - boundary_state[3] << ";" << endl;
	//				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[4] - boundary_state[4] << ";" << endl;
	//				GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[5] - boundary_state[5] << ";" << endl;
	//				GMATDebug << "Initial States: ";
	//				for (int counter = 0; counter < 6; counter++){
	//					GMATDebug << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[counter] << ", ";
	//				}
	//				GMATDebug << endl;
	//				break;
	//			case 1: //depart from parking orbit
	//				GMATDebug << "case: 1" << endl;
	//				//TODO::
	//				//journey_departure_elements_type == 1 means given in COE, which is only relevant if 
	//				//departure type is 'depart from parking orbit' or 'free - direct departure' OR
	//				//arrival type is 'insertion into parking orbit (use chemical Isp)' or ...

	//				//I dont think John's code works correctly
	//				//periapse_state_vector = journeys[j].phases[p].calculate_periapse_state_from_asymptote_and_parking_orbit(Vinf_out, options.journey_departure_elements[j][2], options.journey_departure_elements[j][0] - missionbodies[body_index].radius, journeys[j].phases[p].phase_start_epoch, &TheUniverse[options.number_of_journeys - 1], body_index + 1);

	//				//can convert to inertial and then to appropriate frame
	//				//COE2inertial(const double* E_COE, const double mu, double* state)
	//				//reference frame transformation here if needed
	//				if (this->ptr_gmatmission->options.journey_departure_elements_type[j] == 1)
	//				{
	//					for (int k = 0; k < 3; ++k)
	//					{
	//						periapse_position_vector(k) = periapse_state_vector(k);
	//						periapse_velocity_vector(k) = periapse_state_vector(k + 3);
	//					}
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = " << periapse_position_vector(0) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = " << periapse_position_vector(1) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = " << periapse_position_vector(2) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << periapse_velocity_vector(0) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << periapse_velocity_vector(1) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << periapse_velocity_vector(2) << ";" << endl;
	//				}
	//				break;
	//			case 2: //free direct departure
	//				GMATDebug << "case: 2" << endl;
	//				//TODO::
	//				//journey_departure_elements_type == 1 means given in COE, which is only relevant if 
	//				//departure type is 'depart from parking orbit' or 'free - direct departure' OR
	//				//arrival type is 'insertion into parking orbit (use chemical Isp)' or ...

	//				//I dont think John's code works correctly
	//				//periapse_state_vector = journeys[j].phases[p].calculate_periapse_state_from_asymptote_and_parking_orbit(Vinf_out, options.journey_departure_elements[j][2], options.journey_departure_elements[j][0] - missionbodies[body_index].radius, journeys[j].phases[p].phase_start_epoch, &TheUniverse[options.number_of_journeys - 1], body_index + 1);

	//				//can convert to inertial and then to appropriate frame
	//				//COE2inertial(const double* E_COE, const double mu, double* state)
	//				//reference frame transformation here if needed
	//				if (this->ptr_gmatmission->options.journey_departure_elements_type[j] == 1)
	//				{
	//					for (int k = 0; k < 3; ++k)
	//					{
	//						periapse_position_vector(k) = periapse_state_vector(k);
	//						periapse_velocity_vector(k) = periapse_state_vector(k + 3);
	//					}
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = " << periapse_position_vector(0) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = " << periapse_position_vector(1) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = " << periapse_position_vector(2) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << periapse_velocity_vector(0) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << periapse_velocity_vector(1) << ";" << endl;
	//					GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << periapse_velocity_vector(2) << ";" << endl;
	//				}
	//				break;
	//			case 3: //depart from flyby
	//				GMATDebug << "case: 3" << endl;
	//				break;
	//			}
	//		}//end of if(p == 0)

	//		//for all other phases, interpret state from beginning of each phase and then add the minimum flyby altitude
	//		else
	//		{
	//			//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
	//			missionbodies[body_index].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);

	//			//calculate v_infinity vectors
	//			for (int k = 0; k < 3; ++k)
	//			{
	//				Vinf_in(k)  = this->ptr_gmatmission->journeys[j].phases[p - 1].state_at_end_of_phase[k + 3]   - boundary_state[k + 3];
	//				Vinf_out(k) = this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
	//			}

	//			periapse_state_vector = this->ptr_gmatmission->journeys[j].phases[p].calculate_flyby_periapse_state(Vinf_in, Vinf_out, this->ptr_gmatmission->journeys[j].phases[p - 1].flyby_altitude, missionbodies[body_index]);

	//			//add this position vector to state's initial guess
	//			GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X = "  << periapse_state_vector(0) << ";" << endl;
	//			GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y = "  << periapse_state_vector(1) << ";" << endl;
	//			GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z = "  << periapse_state_vector(2) << ";" << endl;
	//			GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX = " << periapse_state_vector(3) << ";" << endl;
	//			GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY = " << periapse_state_vector(4) << ";" << endl;
	//			GMATfile << "	" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ = " << periapse_state_vector(5) << ";" << endl;
	//		}
	//		GMATfile << "	" << spacecraft_forward_names[name_index] << ".CoordinateSystem = " << missionbodies[body_index].name << "J2000Eq;" << endl;
	//		GMATfile << endl;


	//		//BACKWARD PROPAGATED SPACECRAFT
	//		GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward s/c initial conditions" << endl;

	//		//interpret state from end of each journey based on arrival type (//TODO::)
	//		if (p == this->ptr_gmatmission->journeys[j].number_of_phases - 1)
	//		{
	//			//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
	//			missionbodies[body_index + 1].locate_body(this->ptr_gmatmission->journeys[j].phases[p].phase_end_epoch, boundary_state, false, &this->ptr_gmatmission->options);

	//			switch (this->ptr_gmatmission->options.journey_arrival_type[j])
	//			{
	//			case 0: //parking orbit insertion
	//				break;
	//			case 1: //rendezvous (chemical)
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5] << ";" << endl;
	//				break;
	//			case 2: //flyby with bounded VHP
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5] << ";" << endl;
	//				break;
	//			case 3: //rendezvous (LT) // (EMTG default mission)
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4] << ";" << endl;
	//				GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5] << ";" << endl;
	//				break;
	//			case 4: //match Vinf vector (chemical)
	//				break;
	//			case 5: //match Vinf vector (LT)
	//				break;
	//			}
	//		}

	//		//for all other phases, interpret state from end of each phase and then add the minimum flyby altitude
	//		else
	//		{

	//			//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
	//			missionbodies[body_index + 1].locate_body(this->ptr_gmatmission->journeys[j].phases[p + 1].phase_start_epoch, boundary_state, false, &this->ptr_gmatmission->options);

	//			//calculate v_infinity vectors
	//			for (int k = 0; k < 3; ++k)
	//			{
	//				Vinf_in(k)  = this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[k + 3]           - boundary_state[k + 3];
	//				Vinf_out(k) = this->ptr_gmatmission->journeys[j].phases[p + 1].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
	//			}

	//			periapse_state_vector = this->ptr_gmatmission->journeys[j].phases[p].calculate_flyby_periapse_state(Vinf_in, Vinf_out, this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude, missionbodies[body_index + 1]);

	//			//add this position vector to state's initial guess
	//			GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X = "  << periapse_state_vector(0) << ";" << endl;
	//			GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y = "  << periapse_state_vector(1) << ";" << endl;
	//			GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z = "  << periapse_state_vector(2) << ";" << endl;
	//			GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX = " << periapse_state_vector(3) << ";" << endl;
	//			GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY = " << periapse_state_vector(4) << ";" << endl;
	//			GMATfile << "	" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ = " << periapse_state_vector(5) << ";" << endl;
	//		}
	//		GMATfile << "	" << spacecraft_backward_names[name_index] << ".CoordinateSystem = " << missionbodies[body_index + 1].name << "J2000Eq;" << endl;
	//		GMATfile << endl;

	//		//insert scaled launch epoch and wait times for each subsequent journey
	//		if ((j > 0) && (p == 0))
	//		{
	//			GMATfile << "	%Guess for scaled journey wait times" << endl;
	//			GMATfile << "	Journey" << j + 1 << "_WaitTime_Scaled = " << ((this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch - this->ptr_gmatmission->journeys[j - 1].phases[this->ptr_gmatmission->journeys[j - 1].number_of_phases - 1].phase_end_epoch) - this->ptr_gmatmission->options.journey_wait_time_bounds[j][0]) / (this->ptr_gmatmission->options.journey_wait_time_bounds[j][1] - this->ptr_gmatmission->options.journey_wait_time_bounds[j][0]) << ";" << endl;
	//		}
	//		if (j == 0)
	//		{
	//			if (p == 0)
	//			{
	//				GMATfile << "	%Guess for scaled spacecraft launch epoch" << endl;
	//				GMATfile << "	LaunchEpoch_Scaled = " << (this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch - LaunchDate_LowerBounds) / (LaunchDate_UpperBounds - LaunchDate_LowerBounds) << ";" << endl;
	//			}
	//		}
	//		GMATfile << endl;

	//		//insert scaled states for forward and backward s/c (doing the math in GMAT bc of all the if statements above!)
	//		GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Forward s/c scaled initial conditions" << endl;
	//		GMATfile << "	" << spacecraft_forward_names[name_index] << "_X_Scaled = ("  << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.X - "  << Forward_Flyby_Distance_LowerBound[p] << ") / (" << Forward_Flyby_Distance_UpperBound[p] << " - " << Forward_Flyby_Distance_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_forward_names[name_index] << "_Y_Scaled = ("  << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Y - "  << Forward_Flyby_Distance_LowerBound[p] << ") / (" << Forward_Flyby_Distance_UpperBound[p] << " - " << Forward_Flyby_Distance_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_forward_names[name_index] << "_Z_Scaled = ("  << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.Z - "  << Forward_Flyby_Distance_LowerBound[p] << ") / (" << Forward_Flyby_Distance_UpperBound[p] << " - " << Forward_Flyby_Distance_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_forward_names[name_index] << "_VX_Scaled = (" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VX - " << Forward_Flyby_Velocity_LowerBound[p] << ") / (" << Forward_Flyby_Velocity_UpperBound[p] << " - " << Forward_Flyby_Velocity_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_forward_names[name_index] << "_VY_Scaled = (" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VY - " << Forward_Flyby_Velocity_LowerBound[p] << ") / (" << Forward_Flyby_Velocity_UpperBound[p] << " - " << Forward_Flyby_Velocity_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_forward_names[name_index] << "_VZ_Scaled = (" << spacecraft_forward_names[name_index] << "." << missionbodies[body_index].name << "J2000Eq.VZ - " << Forward_Flyby_Velocity_LowerBound[p] << ") / (" << Forward_Flyby_Velocity_UpperBound[p] << " - " << Forward_Flyby_Velocity_LowerBound[p] << ")" << endl;
	//		GMATfile << "   " << spacecraft_forward_names[name_index] << "_FuelMass_Scaled = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[6] / (this->ptr_gmatmission->options.maximum_mass) << endl;
	//		//GMATfile << "	FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[6] / (this->ptr_gmatmission->options.maximum_mass) << endl;
	//		GMATfile << endl;

	//		GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward s/c scaled initial conditions" << endl;
	//		GMATfile << "	" << spacecraft_backward_names[name_index] << "_X_Scaled = ("  << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.X - "  << Backward_Flyby_Distance_LowerBound[p] << ") / (" << Backward_Flyby_Distance_UpperBound[p] << " - " << Backward_Flyby_Distance_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_backward_names[name_index] << "_Y_Scaled = ("  << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Y - "  << Backward_Flyby_Distance_LowerBound[p] << ") / (" << Backward_Flyby_Distance_UpperBound[p] << " - " << Backward_Flyby_Distance_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_backward_names[name_index] << "_Z_Scaled = ("  << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.Z - "  << Backward_Flyby_Distance_LowerBound[p] << ") / (" << Backward_Flyby_Distance_UpperBound[p] << " - " << Backward_Flyby_Distance_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_backward_names[name_index] << "_VX_Scaled = (" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VX - " << Backward_Flyby_Velocity_LowerBound[p] << ") / (" << Backward_Flyby_Velocity_UpperBound[p] << " - " << Backward_Flyby_Velocity_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_backward_names[name_index] << "_VY_Scaled = (" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VY - " << Backward_Flyby_Velocity_LowerBound[p] << ") / (" << Backward_Flyby_Velocity_UpperBound[p] << " - " << Backward_Flyby_Velocity_LowerBound[p] << ")" << endl;
	//		GMATfile << "	" << spacecraft_backward_names[name_index] << "_VZ_Scaled = (" << spacecraft_backward_names[name_index] << "." << missionbodies[body_index + 1].name << "J2000Eq.VZ - " << Backward_Flyby_Velocity_LowerBound[p] << ") / (" << Backward_Flyby_Velocity_UpperBound[p] << " - " << Backward_Flyby_Velocity_LowerBound[p] << ")" << endl;
	//		GMATfile << "   " << spacecraft_backward_names[name_index] << "_FuelMass_Scaled = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[6] / (this->ptr_gmatmission->options.maximum_mass) << endl;
	//		//GMATfile << "	FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[6] / (this->ptr_gmatmission->options.maximum_mass) << endl;
	//		GMATfile << endl;
	//		++name_index;
	//		++body_index;



	//		////initial guess for inter-phase control
	//		////this means thrust vectors (and sometimes Isp) for low-thrust phases
	//		////or burn index for MGADSM and MGANDSM phases. No control for MGA phases
	//		////initialize thrust vector directions
	//		////define thrust unit vector bounds to scale to between 0 and 1
	//		//double ThrustUnitVector_lowerbounds = -1;
	//		//double ThrustUnitVector_upperbounds =  1;

	//		////propagate forward s/c using finite burns (must scale unit vectors to between 0 and 1)
	//		//GMATfile << endl;
	//		//GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Forward s/c scaled thrust vectors" << endl;
	//		//for (int step = 0; step < (this->ptr_gmatmission->options.num_timesteps / 2); ++step)
	//		//{
	//		//	GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][0] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
	//		//	GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][1] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
	//		//	GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][2] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
	//		//	GMATfile << endl;
	//		//}
	//		////propagate backward s/c using finite burns (must scale unit vectors to between 0 and 1)
	//		//GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << " Backward s/c scaled thrust vectors" << endl;
	//		//for (int step = this->ptr_gmatmission->options.num_timesteps - 1; (step >= this->ptr_gmatmission->options.num_timesteps / 2); --step)
	//		//{
	//		//	GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << ", Time Step #" << step + 1 << ", Backward Propagation" << endl;
	//		//	GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][0] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
	//		//	GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][1] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
	//		//	GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = " << (this->ptr_gmatmission->journeys[j].phases[p].control[step][2] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
	//		//	GMATfile << endl;
	//		//}



	//		




	//	}//end of phase for-statement
	//}//end of journey for-statement
	//GMATfile << "EndScript" << endl;
	//GMATfile << endl;

}//end of write_GMAT_initialguess() method


// method to create the beginmissionsequence statement
void gmatscripter::write_GMAT_beginmissionsequence(){

	GMATfile << "BeginMissionSequence" << endl;
	GMATfile << endl;
	GMATfile << endl;

}


// method to create the optimization sequence for the mission
void gmatscripter::write_GMAT_optimization() {

	//declarations
	stringstream tempstream0;
	stringstream tempstream1;
	stringstream tempstream2;
	stringstream tempstream3;
	stringstream tempstream;
	double epoch;

	GMATfile << "Optimize 'OptimizeSequence' NLPObject {SolveMode = Solve, ExitMode = DiscardAndContinue}" << endl;

	// vary launch epoch
	// WITH 
	// constraint on launch upper and lower bound
	GMATfile << "	" << "%Vary Launch Epoch" << endl;
	this->aux_GMAT_vary("LaunchWindowScaling");
	//set launch spacecraft epoch
	this->aux_GMAT_calculate(spacecraft[0].Name, "(LaunchWindowScaling * LaunchWindow + LaunchWindowOpenDate)");


	// vary journey waittime(s) for interphase spacecraft
	// WITH
	// constraint on journey times upper and lower bounds


	//vary the thrust vectors for each gmat step
	//gmat_step_thrust_vectors
	//for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
	//	for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {
	//		for (int gs = 0; gs < gmat_step_timesteps[j*p + p]; ++gs) {
	//			//vary thrust direction
	//			tempstream0 << "ThrustVector_j" << j << "p" << p << "(" << gs << ", 0)";
	//			tempstream1 << "ThrustVector_j" << j << "p" << p << "(" << gs << ", 1)";
	//			tempstream2 << "ThrustVector_j" << j << "p" << p << "(" << gs << ", 2)";
	//			this->aux_GMAT_vary(tempstream0.str());
	//			this->aux_GMAT_vary(tempstream1.str());
	//			this->aux_GMAT_vary(tempstream2.str());
	//			//calculate the thrust unit vector magnitude
	//			tempstream << "ThrustUnitVectorMagnitude_j" << j << "p" << p << "(" << gs << ", 1)";
	//			tempstream3 << "sqrt(( " << tempstream0.str() << " * 2 - 1) ^ 2 + ( " << tempstream1.str() << " * 2 - 1) ^ 2 + ( " << tempstream2.str() << " * 2 - 1) ^ 2)";
	//			this->aux_GMAT_calculate(tempstream.str(), tempstream3.str());
	//			//constrain the thrust unit vector magnitude
	//			this->aux_GMAT_nonlinearconstraint(tempstream.str(),"<=","1");


	//			GMATfile << endl;

	//			//clear the tempstreams
	//			tempstream0.str(""); tempstream1.str(""); tempstream2.str(""); tempstream3.str(""); tempstream.str("");


	//		}
	//	}
	//}





	//GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	//GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	//GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
	////calculate the thrust unit vector magnitude of the step
	//GMATfile << "	'CalcThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 <<
	//	"' ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = sqrt(" <<
	//	"(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " <<
	//	"(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " <<
	//	"(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2);" << endl;
	////create a nonlinear constraint for the thrust unit vector magnitude
	//GMATfile << "	NonlinearConstraint 'ConstraintThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) <= 1)" << endl;
	////assign thrust unit directions to thruster
	//GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
	//GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
	//GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
	////assign value to thrust coefficient
	//GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1' Thruster_Journey" << j + 1 <<

	





	



	

	// if launch or direct departure, can vary launch vector
	// WITH
	// constraint on max velocity magnitude

	// vary flyby locations for interphase spacecraft

	// vary spacecraft fuelmass 

	//optimize the user-defined objective function
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Objective Function" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//currently vf13 does not like to obey set bounds if an objective function is given
	GMATfile << "	%Minimize objective function" << endl;
	GMATfile << "	%Minimize NLPObject(ObjectiveFunction)" << endl;
	GMATfile << endl;
	GMATfile << "EndOptimize" << endl;

	//TEMPORARY DEBUG CODE
	GMATDebug << endl;
	GMATDebug << " *** spacecraft vector<struct> print out *** " << endl;
	for (int index = 0; index < spacecraft.size(); ++index) {
		GMATDebug << spacecraft[index].Name << " | " << "Epoch: " << spacecraft[index].Epoch << endl;
	}

}


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
				//this->aux_GMAT_beginburn(j, p, step, spacecraft_forward_names[name_index], "Forward");
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
				//aux_GMAT_endburn(j, p, step, spacecraft_forward_names[name_index], "Forward");

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
				//this->aux_GMAT_beginburn(j, p, step, spacecraft_backward_names[name_index], "Backward");
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
				//aux_GMAT_endburn(j, p, step, spacecraft_backward_names[name_index], "Backward");

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


// method to populate an array of type double that holds x,y,z thrust vector information
void gmatscripter::aux_GMAT_populate_thrustvector(int j, int p, int s, vector <double>& x, double lower_bound, double upper_bound) {
	for (int index = 0; index < 3; ++index) {
		x.push_back((this->ptr_gmatmission->journeys[j].phases[p].control[s][index] - lower_bound) / (upper_bound - lower_bound));
	}
}


// method to write a GMAT 'Vary' line
void gmatscripter::aux_GMAT_vary(string object2vary) {
	GMATfile << "	" << "Vary 'Vary " << object2vary << "' NLPObject(" << object2vary << " = " << object2vary << ", {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
}


// method to write a GMAT 'Calculate' line
void gmatscripter::aux_GMAT_calculate(string object2calculate, string rhs) {
	GMATfile << "   " << "'Calculate " << object2calculate << "' " << object2calculate << " = " << rhs << endl;
}


// method to write a GMAT 'NonlinearConstraint' line
void gmatscripter::aux_GMAT_nonlinearconstraint(string object2constrain, string relation, string rhs) {
	GMATfile << "	NonlinearConstraint '" << object2constrain << "' NLPObject( " << object2constrain << " " << relation << " " << rhs << " )" << endl;
}


// method to write a GMAT BeginFiniteBurn Command
void gmatscripter::aux_GMAT_beginburn(string finiteburnobject, string spacecraft_name){
	GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn " << finiteburnobject << "' " << finiteburnobject << "( " << spacecraft_name << " )" << endl;
	//previous aux_GMAT_beginburn resembled the following....
	//aux_GMAT_endburn(int j, int p, int s, string spacecraft_name, string prefix)
	//GMATfile << "	EndFiniteBurn 'EndFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << s + 1
	//	<< prefix << "' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << prefix << "(" << spacecraft_name << ");" << endl;
}


// method to write a GMAT EndFiniteBurn Command
void gmatscripter::aux_GMAT_endburn(string finiteburnobject, string spacecraft_name){
	GMATfile << "	EndFiniteBurn 'EndFiniteBurn " << finiteburnobject << "' " << finiteburnobject << "( " << spacecraft_name << " )" << endl;
	//GMATfile << "	EndFiniteBurn 'EndFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << s + 1
	//	<< prefix << "' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << prefix << "(" << spacecraft_name << ");" << endl;
}
//void gmatscripter::aux_GMAT_endburn(int j, int p, int s, string spacecraft_name, string prefix){
//
//	GMATfile << "	EndFiniteBurn 'EndFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << s + 1 
//					  << prefix << "' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << prefix << "(" << spacecraft_name << ");" << endl;
//
//}


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
	for (int body_index = 0; body_index < missionbodies.size(); ++body_index) {
		GMATfile << " " << missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << endl;

}


// method to write a GMAT PenDown Command
void gmatscripter::aux_GMAT_penDown(){

	GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
	for (int body_index = 0; body_index < missionbodies.size(); ++body_index) {
		GMATfile << " " << missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << endl;

}


// method to write a GMAT ForceModel Resource
void gmatscripter::create_GMAT_forcemodel(string forcemodelname, string centralbody, string pointmasses){

	GMATfile << "Create ForceModel " << forcemodelname << endl;
	GMATfile << forcemodelname << ".CentralBody = " << centralbody << endl;
	GMATfile << forcemodelname << ".PointMasses = {" << pointmasses << "};" << endl;
	GMATfile << forcemodelname << ".Drag = None;" << endl;
	GMATfile << forcemodelname << ".SRP = Off;" << endl;
	GMATfile << endl;

}


// method to write a GMAT Propagator Resource
void gmatscripter::create_GMAT_propagator(string propagatorname, string forcemodelname, bool isCloseApproach) {

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


// method to write a GMAT Spacecraft Resource
void gmatscripter::create_GMAT_spacecraft(int j, int p, string spacecraftname, string prefix, string thebody_coordinatesystem) {

	//declarations
	double epoch;
	string fueltankname;
	string thrustername;
	gmat_spacecraft aSpaceCraft;
	struct gmat_tank aTank;
	struct gmat_thruster aThruster;

	//decision branching
	if (prefix == "Forward") {
		epoch = (this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch / 86400.0) + 2400000.5 - 2430000;
		fueltankname = fueltank_forward_names[p];
		thrustername = thruster_forward_names[p];
	}
	else if (prefix == "Backward") {
		epoch = (this->ptr_gmatmission->journeys[j].phases[p].phase_end_epoch / 86400.0) + 2400000.5 - 2430000;
		fueltankname = fueltank_backward_names[p];
		thrustername = thruster_backward_names[p];
	}

	//write out the spacecraft information
	GMATfile << "% Journey #" << j << ", Phase #" << p << ", " << prefix << " Propagated Spacecraft" << endl;
	GMATfile << "Create Spacecraft " << spacecraftname << endl;
	GMATfile << spacecraftname << ".DateFormat = TAIModJulian;" << endl;
	GMATfile << spacecraftname << ".Epoch = " << epoch << ";" << endl;
	GMATfile << spacecraftname << ".DryMass = 0" << endl;
	GMATfile << spacecraftname << ".CoordinateSystem = " << thebody_coordinatesystem << "J2000Eq;" << endl;
	GMATfile << spacecraftname << ".Tanks = {" << fueltankname << "};" << endl;
	GMATfile << spacecraftname << ".Thrusters = {" << thrustername << "};" << endl;
	GMATfile << endl;

	//create a tank
	aTank.Name = "{ " + fueltankname + " }";
	aThruster.Name = "{ " + thrustername + " }";
	aThruster.Tank = aTank;
	//create and push_back a spacecraft struct type
	aSpaceCraft.Name = spacecraftname;
	aSpaceCraft.Epoch = epoch;
	aSpaceCraft.CoordinateSystem = thebody_coordinatesystem + "J2000Eq";
	aSpaceCraft.Tanks.push_back(aTank);
	aSpaceCraft.Thrusters.push_back(aThruster);
	spacecraft.push_back(aSpaceCraft);

}


// method to write a GMAT FuelTank Resource
void gmatscripter::create_GMAT_fueltank(int j, int p, string fueltankname, string prefix) {

	GMATfile << "% Journey #" << j << ", Phase #" << p << ", " << prefix << " Spacecraft FuelTank" << endl;
	GMATfile << "Create FuelTank " << fueltankname << endl;
	GMATfile << fueltankname << ".AllowNegativeFuelMass = false;" << endl;
	GMATfile << fueltankname << ".Volume = 10;" << endl;
	if (prefix == "Forward") {
		GMATfile << fueltankname << ".FuelMass = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_beginning_of_phase[6] << endl;
	}
	else if (prefix == "Backward") {
		GMATfile << fueltankname << ".FuelMass = " << this->ptr_gmatmission->journeys[j].phases[p].state_at_end_of_phase[6] << endl;
	}
	GMATfile << endl;

}


// method to write a GMAT Thruster Resource
void gmatscripter::create_GMAT_thruster(int j, int p, string thrustername, string thebody_coordinatesystem, string fueltankname) {

	GMATfile << "% Journey #" << j << ", Phase #" << p << ", Thruster Name: " << thrustername << endl;
	GMATfile << "Create Thruster " << thrustername << endl;
	GMATfile << thrustername << ".CoordinateSystem = " << thebody_coordinatesystem << "J2000Eq" << endl;
	GMATfile << thrustername << ".ThrustDirection1 = 1;" << endl;
	GMATfile << thrustername << ".ThrustDirection2 = 0;" << endl;
	GMATfile << thrustername << ".ThrustDirection3 = 0;" << endl;
	GMATfile << thrustername << ".DutyCycle = 1;" << endl;
	GMATfile << thrustername << ".Tank = " << fueltankname << endl;
	GMATfile << thrustername << ".ThrustScaleFactor = 1;" << endl;
	GMATfile << thrustername << ".DecrementMass = true;" << endl;
	GMATfile << thrustername << ".C1 = .1;" << endl;
	GMATfile << thrustername << ".K1 = 3000;" << endl;
	GMATfile << endl;

}


// method to write a GMAT FiniteBurn Resource
void gmatscripter::create_GMAT_finiteburn(string finiteburnname, string thrustername) {

	GMATfile << "Create FiniteBurn " << finiteburnname << endl;
	GMATfile << finiteburnname << ".Thrusters = {" << thrustername << "};" << endl;

}


// method to write a GMAT CoordinateSystem Resource
void gmatscripter::create_GMAT_coordinatesystem(string bodyname) {

	GMATfile << "Create CoordinateSystem " << bodyname << "J2000Eq;" << endl;
	GMATfile << bodyname << "J2000Eq.Origin = " << bodyname << ";" << endl;
	GMATfile << bodyname << "J2000Eq.Axes = MJ2000Eq;" << endl;
	GMATfile << endl;

}

} // end of EMTG namespace