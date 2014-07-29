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

#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"

#include "SpiceUsr.h"

#include <string>  //string
#include <fstream> //ofstream


namespace EMTG {

// default constructor; not intended for use
gmatscripter::gmatscripter(){}

// constructor
gmatscripter::gmatscripter(mission* mission_in){
	this->ptr_gmatmission = mission_in;
	//instantiate a 'gmatmission'
	gmatmission gmatMission(mission_in);
	GMATMission = gmatMission;
}

// destructor
gmatscripter::~gmatscripter(){}

// method to write out the GMAT script
void gmatscripter::write_GMAT_script(){

	//USE THE FOLLOWING to decided if we should use feasible or optimize mode
	//this->ptr_gmatmission->options.NLP_solver_mode

	// open a file for writing
	this->create_GMAT_file();

	// collect the mission level parameters for the GMAT mission
	this->create_GMAT_missions();

	// collect the journey level parameters for the GMAT mission
	this->create_GMAT_journeys();

	// collect the phase level parameters for the GMAT mission
	this->create_GMAT_phases();

	// collect the step level parameters for the GMAT mission
	this->create_GMAT_steps();

	// collect more phase level parameters for the GMAT mission
	this->postpass_GMAT_phases();

	// collect more journey level parameters for the GMAT mission
	this->postpass_GMAT_journeys();

	// collect more mission level parameters for the GMAT mission
	this->postpass_GMAT_missions();

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


// method to get the mission level parameters  
void gmatscripter::create_GMAT_missions() {

	//declarations
	double LaunchDate_LowerBounds;
	double LaunchDate_UpperBounds;
	double LaunchWindow;
	double LaunchDate;
	double LaunchScaling;
	//temporary descriptive parameters
	LaunchDate_LowerBounds = (ptr_gmatmission->Xlowerbounds[0] / 86400.0) + TAIModJOffset;
	LaunchDate_UpperBounds = (ptr_gmatmission->Xupperbounds[0] / 86400.0) + TAIModJOffset;
	LaunchWindow = LaunchDate_UpperBounds - LaunchDate_LowerBounds;
	LaunchDate = (ptr_gmatmission->Xopt[0] / 86400.0) + TAIModJOffset;
	LaunchScaling = (LaunchDate - LaunchDate_LowerBounds) / LaunchWindow;


	// ---------------------------
	//      VARIABLE CREATION
	// ---------------------------
	//store mission level variables that will be needed during the optimization sequence in GMAT
	GMATMission.setVariable("ObjectiveFunction");
	//launch window open date
	GMATMission.setVariable("LaunchWindowOpenDate", std::to_string(LaunchDate_LowerBounds));
	//launch window
	GMATMission.setVariable("LaunchWindow", std::to_string(LaunchWindow));
	//launch window scaling
	GMATMission.setVariable("LaunchWindowScaling", std::to_string(LaunchScaling));
	//  Sometimes ENGINE Variables are Mission Type (i.e. fixed thrust Isp or Impulsive)
	//+ OTHERWISE SEE "get_GMAT_missionlevelparameters" for creation of 'phase_level_parameters'
	//if() low thrust & fixed thrust-Isp, else() impulsive
	if (GMATMission.isLT && this->ptr_gmatmission->options.engine_type == 0) {
		//Isp
		GMATMission.setVariable("ThrusterIsp", std::to_string(this->ptr_gmatmission->options.IspLT));
		//DutyCycle
		GMATMission.setVariable("ThrusterDutyCycle", std::to_string(this->ptr_gmatmission->options.engine_duty_cycle));
		//MaxThrust
		GMATMission.setVariable("ThrusterMaxThrust", std::to_string(this->ptr_gmatmission->options.Thrust));
	}
	else if (!GMATMission.isLT) {
		//Isp
		GMATMission.setVariable("ThrusterIsp", std::to_string(this->ptr_gmatmission->options.IspChem));
	}


	// ------------------------
	//      VARY CREATION
	// ------------------------
	GMATMission.setVary("LaunchWindowScaling");



	// ----------------------------
	//      CALCULATE CREATION
	// ----------------------------
	


	// -----------------------------
	//      CONSTRAINT CREATION
	// -----------------------------
	

	
	//this->ptr_gmatmission->options.total_flight_time_bounds;


	//variables depenedent on the objective function
	//if (this->ptr_gmatmission->options.objective_type == 2) {
	//	mission_level_variables.push_back("FinalMass_Scaled");
	//}


	//DEBUG FILE
	GMATDebug << " ------------ " << endl;
	GMATDebug << "LaunchDate_LowerBounds: " << LaunchDate_LowerBounds << endl;
	GMATDebug << "LaunchDate_UpperBounds: " << LaunchDate_UpperBounds << endl;
	GMATDebug << "LaunchDate: " << LaunchDate << endl;
	GMATDebug << " ------------ " << endl;

}//end of get_GMAT_missionlevelparameters() method


// method to get the journey level parameters
void gmatscripter::create_GMAT_journeys() {

	GMATDebug << endl;
	GMATDebug << "create_GMAT_journeys()" << endl;
	GMATDebug << endl; 

	for (int index = 0; index < this->ptr_gmatmission->options.number_of_journeys; ++index) {
		
		// -----------------------------
		//         INSTANTIATION
		// -----------------------------
		gmatjourney agmatjourney(&GMATMission, index);

		// -----------------------------
		//       VARIABLE CREATION
		// -----------------------------
		agmatjourney.setVariable(agmatjourney.id + "_TimeScaling", 1.0);


		// -----------------------------
		//         VARY CREATION
		// -----------------------------
		agmatjourney.setVary(agmatjourney.id + "_TimeScaling", 1e-005, 0.0, 10e2, 1.0);


		// -----------------------------
		//      CALCULATE CREATION
		// -----------------------------
		


		// -----------------------------
		//      CONSTRAINT CREATION
		// -----------------------------
		// journey time bounds (0: unbounded, 1: bounded flight time, 2: bounded arrival date)
		if (this->ptr_gmatmission->options.journey_timebounded[index] == 0) {
			
		}



		// -----------------------------
		//          PUSH_BACK
		// -----------------------------
		GMATMission.myjourneys.push_back(agmatjourney);

	}
}


// method to get the phase level parameters
void gmatscripter::create_GMAT_phases() {

	GMATDebug << endl;
	GMATDebug << "create_GMAT_phases()" << endl;
	GMATDebug << endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].number_of_emtg_phases; ++p) {

			// -----------------------------
			//         INSTANTIATION
			// -----------------------------
			gmatphase agmatphase(&GMATMission.myjourneys[j], p);

			// -----------------------------
			//       VARIABLE CREATION
			// -----------------------------



			// -----------------------------
			//         VARY CREATION
			// -----------------------------



			// -----------------------------
			//      CALCULATE CREATION
			// -----------------------------
			if (agmatphase.p == 0) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + ".Epoch." + agmatphase.spacecraft_forward.DateFormat,
										"(LaunchWindowScaling * LaunchWindow + LaunchWindowOpenDate)");
			}



			// -----------------------------
			//      CONSTRAINT CREATION
			// -----------------------------
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " X", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.X", " =", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.X");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " Y", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.Y", " =", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.Y");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " Z", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.Z", " =", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.Z");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " VX", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.VX", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.VX");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " VY", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.VY", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.VY");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " VZ", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.VZ", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.VZ");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " Mass", agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.Thruster.Tank.Name + ".FuelMass", "=", agmatphase.spacecraft_backward.Name + "." + agmatphase.spacecraft_backward.Thruster.Tank.Name + ".FuelMass");


			// -----------------------------
			//          PUSH_BACK
			// -----------------------------
			GMATMission.myjourneys[j].myphases.push_back(agmatphase);

		}//end of phases for-statement
	}//end of journeys for-statement

}


// method to get the step level parameters
void gmatscripter::create_GMAT_steps() {

	GMATDebug << endl;
	GMATDebug << "create_GMAT_steps()" << endl;
	GMATDebug << endl;

	//declarations
	int body_index;
	double periapse_velocity_magnitude;
	double approximate_time_in_SOI;
	
	// ---------------------------
	//      VARIABLE CREATION
	// ---------------------------

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int s = 0; s < ptr_gmatmission->options.num_timesteps; ++s) {

				//instantiate a new 'gmatstep'
				gmatstep agmatstep(&GMATMission.myjourneys[j].myphases[p], s);

				//this is the step after the matchpoint; we will make this a half step
				if (!agmatstep.myspacecraft->isForward && agmatstep.isMatchPointStep) {
					GMATDebug << "MP+     ";
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(false);
					//scale the stepsize to be one-half
					agmatstep.scale_stepsize(0.5);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}

				//step is contained inside the SOI the entire time
				if (agmatstep.inSOIatStart && agmatstep.inSOIatEnd) {
					GMATDebug << "SOI-SOI ";
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(true);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
				//step starts in SOI, but exits it during emtg step; we should cut the step in two
				else if (agmatstep.inSOIatStart && !agmatstep.inSOIatEnd) {
					GMATDebug << "SOI --> ";
					//for this case we must figure out how long we are in the SOI before exiting
					if (agmatstep.myspacecraft->isForward) { body_index = 0; }
					else { body_index = 1; }
					periapse_velocity_magnitude = sqrt(2.0 * agmatstep.myphase->mybodies[body_index].mu / (agmatstep.myphase->mybodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + math::norm(agmatstep.initial_velocity_diff, 3)*math::norm(agmatstep.initial_velocity_diff, 3));
					approximate_time_in_SOI     = (agmatstep.myphase->mybodies[body_index].r_SOI / periapse_velocity_magnitude);
					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(approximate_time_in_SOI);
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(true);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();

					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(agmatstep.stepsize - approximate_time_in_SOI);
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(false);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
				//step starts outside SOI, but enters it during emtg step
				//we should cut the step in two
				else if (!agmatstep.inSOIatStart && agmatstep.inSOIatEnd) {
					GMATDebug << "-->SOI  ";
					//for this case we must figure out how long we are in the SOI after entering
					if (agmatstep.myspacecraft->isForward) { body_index = 0; }
					else { body_index = 1; }
					periapse_velocity_magnitude = sqrt(2.0 * agmatstep.myphase->mybodies[body_index].mu / (agmatstep.myphase->mybodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + math::norm(agmatstep.final_velocity_diff, 3)*math::norm(agmatstep.final_velocity_diff, 3));
					approximate_time_in_SOI     = (agmatstep.myphase->mybodies[body_index].r_SOI / periapse_velocity_magnitude);
					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(agmatstep.stepsize - approximate_time_in_SOI);
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(false);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();

					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(approximate_time_in_SOI);
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(true);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
				//step is outside SOI the entire time
				else {
					GMATDebug << " -----> ";
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(false);
					//set a variable for these time steps in "free-space"
					agmatstep.allowTheTimeStep2Vary = true;
					agmatstep.setVariable(agmatstep.id + "_TimeStep");
					agmatstep.setVariable(agmatstep.id + "_InitialTimeStep", agmatstep.stepsize);
					if (agmatstep.myspacecraft->isForward) { agmatstep.setCalculate(agmatstep.id + "_TimeStep", agmatstep.myphase->myjourney->id + "_TimeScaling * " + agmatstep.id + "_InitialTimeStep"); }
					else { agmatstep.setCalculate(agmatstep.id + "_TimeStep", "-" + agmatstep.myphase->myjourney->id + "_TimeScaling * " + agmatstep.id + "_InitialTimeStep"); }
					
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}

				//this is the step before the matchpoint; we will make this a half step
				if (agmatstep.myspacecraft->isForward && agmatstep.isMatchPointStep) {
					GMATDebug << "MP-     ";
					//set the ForceModel and Propagator Parameters
					agmatstep.setFMandProp(false);
					//scale the stepsize to be one-half
					agmatstep.scale_stepsize(0.5);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
			}//end of Step level
		}//end of Phase level
	}//end of Journey level

}//end method


// method to generate additional variables, vary, calculate and constraint commands post 'gmat' class() creations
void gmatscripter::postpass_GMAT_phases() {

	GMATDebug << endl;
	GMATDebug << "postpass_GMAT_phases()" << endl;
	GMATDebug << endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			// -----------------------------
			//      POST-DATA COLLECTION
			// -----------------------------
			//calculate the phase's allocated time that can/cannot be varied
			GMATMission.myjourneys[j].myphases[p].get_time_allocation();



			// -----------------------------
			//       VARIABLE CREATION
			// -----------------------------
			//a variable that represents the amount of time that can be scaled during each journey
			GMATMission.myjourneys[j].myphases[p].setVariable(GMATMission.myjourneys[j].myphases[p].id + "_EligableTime", GMATMission.myjourneys[j].myphases[p].eligabletime);
			GMATMission.myjourneys[j].myphases[p].setVariable(GMATMission.myjourneys[j].myphases[p].id + "_MatchPoint_PositionError");
			GMATMission.myjourneys[j].myphases[p].setVariable(GMATMission.myjourneys[j].myphases[p].id + "_MatchPoint_VelocityError");
			GMATMission.myjourneys[j].myphases[p].setVariable(GMATMission.myjourneys[j].myphases[p].id + "_MatchPoint_MassError");
			GMATMission.myjourneys[j].myphases[p].setVariable(GMATMission.myjourneys[j].myphases[p].id + "_TOF");



			// -----------------------------
			//         VARY CREATION
			// -----------------------------



			// -----------------------------
			//      CALCULATE CREATION
			// -----------------------------
			//backward spacecraft epoch
			GMATMission.myjourneys[j].myphases[p].setCalculate(GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".Epoch." + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.DateFormat,
				GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".Epoch." + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.DateFormat + " + ( " + std::to_string(GMATMission.myjourneys[j].myphases[p].ineligabletime) + " / 86400.0 )" +
				" + ( " + GMATMission.myjourneys[j].myphases[p].myjourney->id + "_TimeScaling * " + GMATMission.myjourneys[j].myphases[p].id + "_EligableTime ) / 86400.0");
			//TOF
			GMATMission.myjourneys[j].myphases[p].setCalculate(GMATMission.myjourneys[j].myphases[p].id + "_TOF",
				GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".Epoch." + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.DateFormat + " - " +
				GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".Epoch." + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.DateFormat);
			//position error at matchpoint
			GMATMission.myjourneys[j].myphases[p].setCalculate(GMATMission.myjourneys[j].myphases[p].id + "_MatchPoint_PositionError",
				"sqrt(( " + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".SunJ2000Eq.X - " + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".SunJ2000Eq.X) ^ 2 + " +
				"( " + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".SunJ2000Eq.Y - " + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".SunJ2000Eq.Y) ^ 2 + " +
				"( " + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".SunJ2000Eq.Z - " + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".SunJ2000Eq.Z) ^ 2 )", false);
			//velocity error at matchpoint
			GMATMission.myjourneys[j].myphases[p].setCalculate(GMATMission.myjourneys[j].myphases[p].id + "_MatchPoint_VelocityError",
				"sqrt(( " + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".SunJ2000Eq.VX - " + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".SunJ2000Eq.VX) ^ 2 + " +
				"( " + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".SunJ2000Eq.VY - " + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".SunJ2000Eq.VY) ^ 2 + " +
				"( " + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + ".SunJ2000Eq.VZ - " + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + ".SunJ2000Eq.VZ) ^ 2 )", false);
			//mass error at matchpoint
			GMATMission.myjourneys[j].myphases[p].setCalculate(GMATMission.myjourneys[j].myphases[p].id + "_MatchPoint_MassError",
				GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name + "." + GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Thruster.Tank.Name + ".FuelMass - " +
				GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name + "." + GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Thruster.Tank.Name + ".FuelMass", false);



			// -----------------------------
			//      CONSTRAINT CREATION
			// -----------------------------

		}//end of phases for-statement
	}//end of journeys for-statement

}


// method to generate additional variables, vary, calculate and constraint commands post 'gmat' class() creations
void gmatscripter::postpass_GMAT_journeys() {

	//declarations
	stringstream tempstream;
	
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		// -----------------------------
		//      POST-DATA COLLECTION
		// -----------------------------


		// -----------------------------
		//       VARIABLE CREATION
		// -----------------------------
		//if either mission or myjourney is time-constrained then we need to produce a TOF variable
		if (this->ptr_gmatmission->options.global_timebounded == 1 || this->ptr_gmatmission->options.journey_timebounded[GMATMission.myjourneys[j].j] == 1) {
			GMATMission.myjourneys[j].setVariable(GMATMission.myjourneys[j].id + "_TOF");
		}



		// -----------------------------
		//         VARY CREATION
		// -----------------------------



		// -----------------------------
		//      CALCULATE CREATION
		// -----------------------------
		//if either mission or myjourney is time-constrained then we need to calculate the TOF variable based on phase TOF
		if (this->ptr_gmatmission->options.global_timebounded == 1 || this->ptr_gmatmission->options.journey_timebounded[GMATMission.myjourneys[j].j] == 1) {
			tempstream << GMATMission.myjourneys[j].myphases[0].id << "_TOF";
			for (int p = 1; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
				tempstream << " + " << GMATMission.myjourneys[j].myphases[p].id << "_TOF";
			}
			GMATMission.myjourneys[j].setCalculate(GMATMission.myjourneys[j].id + "_TOF", tempstream.str(), false);
			tempstream.str("");
		}


		// -----------------------------
		//      CONSTRAINT CREATION
		// -----------------------------
		//if myjourney is time-constrained then we need constrain the TOF variable
		if (this->ptr_gmatmission->options.journey_timebounded[GMATMission.myjourneys[j].j] == 1) {
			GMATMission.myjourneys[j].setConstraint(GMATMission.myjourneys[j].id + "_TOF", ">=", this->ptr_gmatmission->options.journey_flight_time_bounds[GMATMission.myjourneys[j].j][0] / 86400.0 );
			GMATMission.myjourneys[j].setConstraint(GMATMission.myjourneys[j].id + "_TOF", "<=", this->ptr_gmatmission->options.journey_flight_time_bounds[GMATMission.myjourneys[j].j][1] / 86400.0);
		}

			
	}//end of journeys for-statement

}


// method to generate additional variables, vary, calculate and constraint commands post 'gmat' class() creations
void gmatscripter::postpass_GMAT_missions() {

	//declarations
	stringstream tempstream;

	// -----------------------------
	//      POST-DATA COLLECTION
	// -----------------------------

	// -----------------------------
	//       VARIABLE CREATION
	// -----------------------------
	//if a mission is time-constrained then we need to create a mission-level TOF to later constrain
	if (this->ptr_gmatmission->options.global_timebounded == 1) { GMATMission.setVariable("Mission_TOF"); }

	// -----------------------------
	//         VARY CREATION
	// -----------------------------


	// -----------------------------
	//      CALCULATE CREATION
	// -----------------------------
	//if a mission is time-constrained then we need to calculate the mission-level TOF
	if (this->ptr_gmatmission->options.global_timebounded == 1) {
		tempstream << GMATMission.myjourneys[0].id << "_TOF";
		for (int j = 1; j < GMATMission.myjourneys.size(); ++j) {
			tempstream << " + " << GMATMission.myjourneys[j].id << "_TOF";
		}
		GMATMission.setCalculate("Mission_TOF", tempstream.str(), false);
		tempstream.str("");
	}

	// -----------------------------
	//      CONSTRAINT CREATION
	// -----------------------------
	//if myjourney is time-constrained then we need constrain the TOF variable
	if (this->ptr_gmatmission->options.global_timebounded == 1) {
		GMATMission.setConstraint("Mission_TOF", ">=", this->ptr_gmatmission->options.total_flight_time_bounds[0] / 86400.0);
		GMATMission.setConstraint("Mission_TOF", "<=", this->ptr_gmatmission->options.total_flight_time_bounds[1] / 86400.0);
	}

}


// method to create the script preamble
void gmatscripter::write_GMAT_preamble() {

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
void gmatscripter::write_GMAT_spacecraft() {

	//spacecraft header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Spacecraft" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//write out the forward spacecraft
			this->create_GMAT_spacecraft(GMATMission.myjourneys[j].myphases[p].spacecraft_forward);
			//write out the backward spacecraft
			this->create_GMAT_spacecraft(GMATMission.myjourneys[j].myphases[p].spacecraft_backward);
		}
	}
	GMATfile << endl;

}//end of write_GMAT_spacecraft() method


// method to create hardware information
void gmatscripter::write_GMAT_hardware() {

	//hardware header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Hardware components" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//write out the forward spacecraft
			this->create_GMAT_fueltank(GMATMission.myjourneys[j].myphases[p].spacecraft_forward);
			this->create_GMAT_thruster(GMATMission.myjourneys[j].myphases[p].spacecraft_forward);
			//write out the backward spacecraft
			this->create_GMAT_fueltank(GMATMission.myjourneys[j].myphases[p].spacecraft_backward);
			this->create_GMAT_thruster(GMATMission.myjourneys[j].myphases[p].spacecraft_backward);
		}
	}
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
	for (int body_index = 0; body_index < GMATMission.missionbodies_unique.size(); ++body_index)
	{
		//must create any bodies that are not already defined in GMAT
		if ((GMATMission.missionbodies[body_index].name != "Sun") && (GMATMission.missionbodies[body_index].name != "Mercury") && (GMATMission.missionbodies[body_index].name != "Venus") && (GMATMission.missionbodies[body_index].name != "Earth") && (GMATMission.missionbodies[body_index].name != "Mars") && (GMATMission.missionbodies[body_index].name != "Jupiter") && (GMATMission.missionbodies[body_index].name != "Saturn") && (GMATMission.missionbodies[body_index].name != "Uranus") && (GMATMission.missionbodies[body_index].name != "Neptune") && (GMATMission.missionbodies[body_index].name != "Pluto"))
		{
			GMATfile << "%Must create model for body visited" << endl;
			GMATfile << "Create Planet " << GMATMission.missionbodies_unique[body_index].name << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".NAIFId = " << GMATMission.missionbodies_unique[body_index].spice_ID << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".EquatorialRadius = " << GMATMission.missionbodies_unique[body_index].radius << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".Mu = " << GMATMission.missionbodies_unique[body_index].mu << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".PosVelSource = 'SPICE'" << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".CentralBody = '" << GMATMission.missionbodies_unique[body_index].central_body_name << "'" << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".RotationDataSource = 'IAUSimplified'" << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".OrientationEpoch = 21545" << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisRAConstant = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.alpha0 << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisRARate = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.alphadot << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisDECConstant = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.delta0 << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisDECRate = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.deltadot << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".RotationConstant = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.W << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".RotationRate = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.Wdot << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".OrbitSpiceKernelName = {";

			//find which spice file the body is located
			EMTG::filesystem::get_all_files_with_extension(fs::path(this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

			for (size_t k = 0; k < SPICE_files.size(); ++k)
			{
				filestring = this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

				//check if body is located in spice file
				scard_c(0, &spice_coverage);
				spkcov_c(filestring.c_str(), GMATMission.missionbodies_unique[body_index].spice_ID, &spice_coverage);
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
		if ((GMATMission.missionbodies_unique[body_index].central_body_name != "Sun") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Mercury") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Venus") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Earth") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Mars") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Jupiter") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Saturn") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Uranus") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Neptune") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Pluto"))
		{
			GMATfile << "%Must create model for central body" << endl;
			GMATfile << "Create Planet " << GMATMission.missionbodies_unique[body_index].central_body_name << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".NAIFId = " << GMATMission.missionbodies_unique[body_index].central_body_spice_ID << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".EquatorialRadius = " << GMATMission.missionbodies_unique[body_index].central_body_radius << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".Mu = " << GMATMission.missionbodies_unique[body_index].universe_mu << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".PosVelSource = 'SPICE'" << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".CentralBody = 'Sun'" << endl; //assume Sun for now
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".RotationDataSource = 'IAUSimplified'" << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".OrientationEpoch = 21545" << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisRAConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alpha0 << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisRARate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alphadot << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisDECConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.delta0 << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisDECRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.deltadot << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".RotationConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.W << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".RotationRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.Wdot << endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".OrbitSpiceKernelName = {";

			//find which spice file the body is located
			EMTG::filesystem::get_all_files_with_extension(fs::path(this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

			for (size_t k = 0; k < SPICE_files.size(); ++k)
			{
				filestring = this->ptr_gmatmission->options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

				//check if body is located in spice file
				scard_c(0, &spice_coverage);
				spkcov_c(filestring.c_str(), GMATMission.missionbodies_unique[body_index].central_body_spice_ID, &spice_coverage);
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
	this->create_GMAT_forcemodel(GMATMission.myjourneys[0].myphases[0].mysteps[0].propagator.ForceModel);
	tempstrings.push_back(GMATMission.myjourneys[0].myphases[0].mysteps[0].propagator.ForceModel.Name);

	//generate the rest of the forcemodels
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				//set bool to true
				hasNotBeenCopied = true;
				//loop through any saved forcemodel names
				for (int k = 0; k < tempstrings.size(); ++k) {
					//if we have already saved out the current forcemodel, 
					//then turn the bool flag to false and do not print the current forcemodel
					if (tempstrings[k].compare(GMATMission.myjourneys[j].myphases[p].mysteps[gs].propagator.ForceModel.Name) == 0) { hasNotBeenCopied = false; }
				}
				//if it hasn't been printed yet, then print it!
				if (hasNotBeenCopied) {
					this->create_GMAT_forcemodel(GMATMission.myjourneys[j].myphases[p].mysteps[gs].propagator.ForceModel);
					tempstrings.push_back(GMATMission.myjourneys[j].myphases[p].mysteps[gs].propagator.ForceModel.Name);
				}
			}
		}
	}
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
	this->create_GMAT_propagator(GMATMission.myjourneys[0].myphases[0].mysteps[0].propagator);
	tempstrings.push_back(GMATMission.myjourneys[0].myphases[0].mysteps[0].propagator.Name);

	//generate the rest of the propagators
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				//set bool to true
				hasNotBeenCopied = true;
				//loop through any saved propagator names
				for (int k = 0; k < tempstrings.size(); ++k) {
					//if we have already saved out the current propagator, 
					//then turn the bool flag to false and do not print the current propagator
					if (tempstrings[k].compare(GMATMission.myjourneys[j].myphases[p].mysteps[gs].propagator.Name) == 0) { hasNotBeenCopied = false; }
				}
				//if it hasn't been printed yet, then print it!
				if (hasNotBeenCopied) {
					this->create_GMAT_propagator(GMATMission.myjourneys[j].myphases[p].mysteps[gs].propagator);
					tempstrings.push_back(GMATMission.myjourneys[j].myphases[p].mysteps[gs].propagator.Name);
				}
			}
		}
	}
	GMATfile << endl;

}//end of write_GMAT_propagators() method


// method to create burn information
void gmatscripter::write_GMAT_burns(){

	//burn header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Burns" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//write out the forward spacecraft
			this->create_GMAT_burn(GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Burn);
			//write out the backward spacecraft
			this->create_GMAT_burn(GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Burn);
		}
	}
	GMATfile << endl;

}//end of write_GMAT_burns() method


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
	this->create_GMAT_coordinatesystem(GMATMission.missionbodies_unique[0].central_body_name);

	for (int b = 0; b < GMATMission.missionbodies_unique.size(); ++b) {
		this->create_GMAT_coordinatesystem(GMATMission.missionbodies_unique[b].name);
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
	GMATfile << "Create OrbitView " << GMATMission.missionbodies_unique[0].central_body_name << "View" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ShowPlot =		false" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.SolverIterations =	 All" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.RelativeZOrder =	501" << endl;

	//add which bodies and s/c to plot 
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.Add =	{";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			GMATfile << GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name << ", ";
			GMATfile << GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name << ", ";
		}
	}
	for (int body_index = 0; body_index < GMATMission.missionbodies_unique.size(); ++body_index)
	{
		GMATfile << GMATMission.missionbodies_unique[body_index].name;
		if (body_index < GMATMission.missionbodies_unique.size() - 1)
		{
			GMATfile << ", ";
		}
	}
	GMATfile << ", " << GMATMission.missionbodies_unique[0].central_body_name << "}" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.CoordinateSystem =		" << GMATMission.missionbodies_unique[0].central_body_name << "J2000Eq" << endl;
	//bool flag parameters for drawing objects
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.DrawObject = [";
	for (int index_plot = 0; index_plot < (GMATMission.missionbodies_unique.size()); ++index_plot)
	{
		GMATfile << "true ";
	}
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			GMATfile << "true true ";
		}
	}
	GMATfile << "true]" << endl;
	//other parameters
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.DataCollectFrequency   = 1" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.UpdatePlotFrequency    = 50" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.NumPointsToRedraw      = 300" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewScaleFactor        = 35" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewPointReference	  = " << GMATMission.missionbodies_unique[0].central_body_name << ";" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewDirection		  = " << GMATMission.missionbodies_unique[0].central_body_name << ";" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewPointVector		  = [ 0 0 30000000 ];" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewUpAxis             = X" << endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.UseInitialView         = On" << endl;
	GMATfile << endl;

	//create orbit views for all bodies visited
	GMATfile << "%Create subscribers for other body views" << endl;
	for (int body_index = 0; body_index < GMATMission.missionbodies_unique.size(); ++body_index)
	{
		GMATfile << "Create OrbitView " << GMATMission.missionbodies_unique[body_index].name << "View" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ShowPlot               = false" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.SolverIterations       = All" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.RelativeZOrder         = 501" << endl;

		//add which bodies and s/c to plot
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.Add                    = {";
		for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
			for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
				GMATfile << GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name << ", ";
				GMATfile << GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name << ", ";
			}
		}
		for (int body_index = 0; body_index < GMATMission.missionbodies_unique.size(); ++body_index)
		{
			GMATfile << GMATMission.missionbodies_unique[body_index].name;
			if (body_index < GMATMission.missionbodies_unique.size() - 1)
			{
				GMATfile << ", ";
			}
		}
		GMATfile << "}" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.CoordinateSystem       = " << GMATMission.missionbodies_unique[body_index].name << "J2000Eq" << endl;
		//bool flag parameters for drawing objects
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.DrawObject             = [";
		for (int index_plot = 0; index_plot < (GMATMission.missionbodies_unique.size()); ++index_plot)
		{
			GMATfile << "true ";
		}
		for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
			for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
				GMATfile << "true true";
			}
		}
		GMATfile << "]" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.DataCollectFrequency   = 1" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.UpdatePlotFrequency    = 50" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.NumPointsToRedraw      = 300" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewScaleFactor        = 35" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewUpAxis             = X" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.UseInitialView         = On" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewPointReference	  = " << GMATMission.missionbodies_unique[body_index].name << ";" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewDirection		  = " << GMATMission.missionbodies_unique[body_index].name << ";" << endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewPointVector		  = [ 0 0 3000000 ];" << endl;
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

	report_name_stream << "Report_Spacecraft_" << GMATMission.missionbodies[0].central_body_name << "_States";
	report_name = report_name_stream.str();
	GMATfile << "Create ReportFile " << report_name << endl;
	GMATfile << report_name << ".SolverIterations = Current;" << endl;
	GMATfile << report_name << ".UpperLeft = [ 0 0 ];" << endl;
	GMATfile << report_name << ".Size = [ 0 0 ];" << endl;
	GMATfile << report_name << ".RelativeZOrder = 0;" << endl;
	GMATfile << report_name << ".Maximized = false;" << endl;
	GMATfile << report_name << ".Filename = 'Report_" << GMATMission.missionbodies[0].central_body_name << "Centered_States.txt';" << endl;
	GMATfile << report_name << ".Precision = 16;" << endl;
	GMATfile << report_name << ".WriteHeaders = true;" << endl;
	GMATfile << report_name << ".LeftJustify = On;" << endl;
	GMATfile << report_name << ".ZeroFill = Off;" << endl;
	GMATfile << report_name << ".ColumnWidth = 20;" << endl;
	GMATfile << report_name << ".WriteReport = true;" << endl;
	GMATfile << endl;

	for (int body_index = 0; body_index < GMATMission.missionbodies_unique.size(); body_index++){
		// new report name
		report_name.erase (report_name.begin(), report_name.end());
		//report_name_stream << "Report_Spacecraft_" << GMATMission.missionbodies_unique[body_index].name << "_States";
		//report_name = report_name_stream.str();
		report_name = "Report_Spacecraft_" + GMATMission.missionbodies_unique[body_index].name + "_States";
		// GMAT script printing
		GMATfile << "Create ReportFile " << report_name << endl;
		GMATfile << report_name << ".SolverIterations = Current;" << endl;
		GMATfile << report_name << ".UpperLeft = [ 0 0 ];" << endl;
		GMATfile << report_name << ".Size = [ 0 0 ];" << endl;
		GMATfile << report_name << ".RelativeZOrder = 0;" << endl;
		GMATfile << report_name << ".Maximized = false;" << endl;
		GMATfile << report_name << ".Filename = 'Report_" << GMATMission.missionbodies_unique[body_index].name << "Centered_States.txt';" << endl;
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

	//XY Plots for visually informing user of the optimization progress
	GMATfile << "%Creating some plots for to inform the user of the optimization progress" << endl;
	GMATfile << "Create XYPlot PositionErrorPlot" << endl;
	GMATfile << "PositionErrorPlot.SolverIterations = All" << endl;
	GMATfile << "PositionErrorPlot.UpperLeft = [0 0]" << endl;
	GMATfile << "PositionErrorPlot.Size = [0 0]" << endl;
	GMATfile << "PositionErrorPlot.RelativeZOrder = 0" << endl;
	GMATfile << "PositionErrorPlot.XVariable = Iterations" << endl;
	GMATfile << "PositionErrorPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_PositionError"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_PositionError"; }
		}
	}
	GMATfile << " }" << endl;
	GMATfile << "PositionErrorPlot.ShowGrid = true" << endl;
	GMATfile << "PositionErrorPlot.ShowPlot = true" << endl;
	GMATfile << endl;

	GMATfile << "Create XYPlot VelocityErrorPlot" << endl;
	GMATfile << "VelocityErrorPlot.SolverIterations = All" << endl;
	GMATfile << "VelocityErrorPlot.UpperLeft = [0 0]" << endl;
	GMATfile << "VelocityErrorPlot.Size = [0 0]" << endl;
	GMATfile << "VelocityErrorPlot.RelativeZOrder = 0" << endl;
	GMATfile << "VelocityErrorPlot.XVariable = Iterations" << endl;
	GMATfile << "VelocityErrorPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_VelocityError"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_VelocityError"; }
		}
	}
	GMATfile << " }" << endl;
	GMATfile << "VelocityErrorPlot.ShowGrid = true" << endl;
	GMATfile << "VelocityErrorPlot.ShowPlot = true" << endl;
	GMATfile << endl;

	GMATfile << "Create XYPlot MassErrorPlot" << endl;
	GMATfile << "MassErrorPlot.SolverIterations = All" << endl;
	GMATfile << "MassErrorPlot.UpperLeft = [0 0]" << endl;
	GMATfile << "MassErrorPlot.Size = [0 0]" << endl;
	GMATfile << "MassErrorPlot.RelativeZOrder = 0" << endl;
	GMATfile << "MassErrorPlot.XVariable = Iterations" << endl;
	GMATfile << "MassErrorPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_MassError"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_MassError"; }
		}
	}
	GMATfile << " }" << endl;
	GMATfile << "MassErrorPlot.ShowGrid = true" << endl;
	GMATfile << "MassErrorPlot.ShowPlot = true" << endl;
	GMATfile << endl;

	GMATfile << "Create XYPlot TOFPlot" << endl;
	GMATfile << "TOFPlot.SolverIterations = All" << endl;
	GMATfile << "TOFPlot.UpperLeft = [0 0]" << endl;
	GMATfile << "TOFPlot.Size = [0 0]" << endl;
	GMATfile << "TOFPlot.RelativeZOrder = 0" << endl;
	GMATfile << "TOFPlot.XVariable = Iterations" << endl;
	GMATfile << "TOFPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_TOF"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_TOF"; }
		}
	}
	GMATfile << " }" << endl;
	GMATfile << "TOFPlot.ShowGrid = true" << endl;
	GMATfile << "TOFPlot.ShowPlot = true" << endl;
	GMATfile << endl;

	//An iteration counter for the Optimization Sequence (and necessary for the XY Plots above)
	GMATMission.setVariable("Iterations", 0);
	GMATMission.setCalculate("Iterations", "Iterations + 1", false);
	GMATfile << endl;

}//end of write_GMAT_subscribers() method


// method to create array and variable information
void gmatscripter::write_GMAT_variables(){

	//arrays and variables header
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "%---------- Arrays, Variables, Strings" << endl;
	GMATfile << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	//write out mission level variables
	GMATfile << "% --- Mission Level: Arrays, Variables, and Strings" << endl;
	GMATMission.printVariable(GMATfile);
	GMATfile << endl;

	//write out journey level variables
	GMATfile << "% --- Journey Level: Arrays, Variables, and Strings" << endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) { GMATMission.myjourneys[j].printVariable(GMATfile); }
	GMATfile << endl;

	//write out the phase level variables
	GMATfile << "% --- Phase Level: Arrays, Variables, and Strings" << endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			GMATMission.myjourneys[j].myphases[p].printVariable(GMATfile);
		}
	}
	GMATfile << endl;

	//write out the step level variables
	GMATfile << "% --- Step Level: Arrays, Variables, and Strings" << endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printVariable(GMATfile);
			}	
		}
	}
	GMATfile << endl;

	//write out report variables
	//create temporary strings to identify time steps during phases
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				GMATfile << "Create String tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tminus" << endl;
				GMATfile << "tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tminus = '" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tminus: '" << endl;
				GMATfile << "Create String tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tplus" << endl;
				GMATfile << "tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tplus = '" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tplus: '" << endl;
			}
		}
	}
	GMATfile << endl;
	GMATfile << endl;

}//end of write_GMAT_variables() method


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

	//write out the phase level variables
	GMATfile << "% --- Phase Level: Arrays, Variables, and Strings" << endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//forward spacecraft
			this->create_GMAT_initialconditions(GMATMission.myjourneys[j].myphases[p].spacecraft_forward);
			//backward spacecraft
			this->create_GMAT_initialconditions(GMATMission.myjourneys[j].myphases[p].spacecraft_backward);
		}
	}

	GMATfile << endl;
	GMATfile << "EndScript" << endl;
	GMATfile << endl;


	////write out the initial guess values for the thrust vector history
	//for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
	//	for (int p = 0; p < this->ptr_gmatmission->journeys[j].number_of_phases; ++p) {
	//		if (p == 0) {
	//			for (int index = 0; index < gmat_steps_per_phase[p]; ++index) {
	//				for (int subindex = 0; subindex < 3; ++subindex) { GMATfile << "	ThrustVector_j" << j << "p" << p << "(" << index + 1 << ", " << subindex << ") = " << gmat_step_thrust_vectors[index][subindex] << endl; }
	//			}
	//		}
	//		else {
	//			for (int index = gmat_steps_per_phase[p - 1] - 1; index < gmat_steps_per_phase[p]; ++index) {
	//				for (int subindex = 0; subindex < 3; ++subindex) { GMATfile << "	ThrustVector_j" << j << "p" << p << "(" << index + 1 << ", " << subindex << ") = " << gmat_step_thrust_vectors[index][subindex] << endl; }
	//			}
	//		}
	//		GMATfile << endl;
	//	}//end of phase for-statement
	//}//end of journey for-statement

	//GMATfile << "EndScript" << endl;
	//GMATfile << endl;

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

	GMATfile << "	" << "%-------------------------------------------------------------------------" << endl;
	GMATfile << "	" << "%---------- Optimization Sequence" << endl;
	GMATfile << "	" << "%-------------------------------------------------------------------------" << endl;
	GMATfile << endl;

	GMATfile << "Optimize 'OptimizeSequence' NLPObject {SolveMode = Solve, ExitMode = DiscardAndContinue}" << endl;
	GMATfile << endl;

	// vary journey waittime(s) for interphase spacecraft
	// WITH
	// constraint on journey times upper and lower bounds

	//TODO:: 
	//include InitialConditions(SpaceCraft States)
	//include Vary(TOF, FuelMass)
	//include Calculate(TimeSteps)
	//include Constraint(TOF, MatchPoint)
	//include Active Subscribers(Reports, Figures)

	//Mission level vary and calculate commands
	GMATMission.printVary(GMATfile);
	GMATMission.printCalculate(GMATfile);
	GMATfile << endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		//Journey ID
		GMATfile << "   " << "% Journey ID: " << GMATMission.myjourneys[j].id << endl;
		//Jouney level vary and calculate commands
		GMATMission.myjourneys[j].printVary(GMATfile);
		GMATMission.myjourneys[j].printCalculate(GMATfile);
		GMATfile << endl;
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//Phase ID
			GMATfile << "   " << "% Phase ID: " << GMATMission.myjourneys[j].myphases[p].id << endl;
			//Phase level vary and calculate commands
			GMATMission.myjourneys[j].myphases[p].printVary(GMATfile);
			GMATMission.myjourneys[j].myphases[p].printCalculate(GMATfile);
			GMATfile << endl;
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				//Step ID
				GMATfile << "   " << "% Step ID: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << endl;
				//'PenUp' command
				PenUp();
				//Step level vary commands
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printVary(GMATfile);
				//Step level calculate commands
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printCalculate(GMATfile);
				//Step level constraint commands
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printConstraint(GMATfile);
				//'BeginBurn' command
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printBeginBurn(GMATfile);
				//'Propagate' command with '0.0' secs if the zeroPropagate flag is 'true' for the 'gmatstep'
				if (GMATMission.myjourneys[j].myphases[p].mysteps[gs].zeroPropagate) {
					this->aux_GMAT_propagate(GMATMission.myjourneys[j].myphases[p].mysteps[gs], true);
				}
				//'PenDown' command
				PenDown();
				//write out the 'Propagate' command
				GMATDebug << "spacecraft name: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].myspacecraft->Name << " isForward: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].myspacecraft->isForward << endl;
				GMATDebug << "id: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << " stepsize: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].stepsize << "  zeroPropagate: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].zeroPropagate;
				this->write_GMAT_report(GMATMission.myjourneys[j].myphases[p].mysteps[gs], true, true);
				this->aux_GMAT_propagate(GMATMission.myjourneys[j].myphases[p].mysteps[gs], false);
				this->write_GMAT_report(GMATMission.myjourneys[j].myphases[p].mysteps[gs], false, true);
				//'EndBurn' command
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printEndBurn(GMATfile);
				//Step level calculate commands
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printCalculate(GMATfile, 0);
				GMATfile << endl;
			}//end of Step level
			//Phase level constraint commands
			GMATMission.myjourneys[j].myphases[p].printCalculate(GMATfile, 0);
			GMATMission.myjourneys[j].myphases[p].printConstraint(GMATfile);
			GMATfile << endl;
		}//end of Phase level
		//Journey level calculate and constraint commands
		GMATMission.myjourneys[j].printCalculate(GMATfile, 0);
		GMATMission.myjourneys[j].printConstraint(GMATfile);
		GMATfile << endl;
	}//end of Journey level

	//Mission level calculate and constraint commands
	GMATMission.printCalculate(GMATfile, 0);
	GMATMission.printConstraint(GMATfile);
	GMATfile << endl;









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
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			GMATDebug << GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name << " | " << "Epoch: " << GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Epoch << endl;
			GMATDebug << GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name << " | " << "Epoch: " << GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Epoch << endl;
		}
	}

}





// method to report results for debugging
void gmatscripter::write_GMAT_report(class gmatstep& agmatstep, bool isbeforemaneuver, bool writecontrolhistory){

	// declaration
	int body_index;
	string tempString;
	stringstream tempStringStream;

	if (agmatstep.myspacecraft->isForward) { body_index = 0; }
	else { body_index = 1; }

	// prefix to be added if before a finite burn maneuver
	if (isbeforemaneuver) { 
		tempStringStream << "tempString_" << agmatstep.id << "_tminus ";
		tempString = tempStringStream.str();
	}
	else {
		tempStringStream << "tempString_" << agmatstep.id << "_tplus ";
		tempString = tempStringStream.str();
	}

	// write the spacecraft controls
	if (writecontrolhistory) {
		GMATfile << "	Report 'Report_SpacecraftControl' Report_SpacecraftControl " << tempString;
		GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->Thruster.Name << ".ThrustDirection1 ";
		GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->Thruster.Name << ".ThrustDirection2 ";
		GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->Thruster.Name << ".ThrustDirection3 ";
		GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->Thruster.Name << ".C1" << endl;
	}
	// write the spacecraft states in the central body frame
	GMATfile << "	Report 'Report_SpacecraftState' Report_Spacecraft_" << agmatstep.myphase->mybodies[body_index].name << "_States " << tempString;
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".X ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".Y ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".Z ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".VX ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".VY ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".VZ ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->Thruster.Tank.Name << ".FuelMass;" << endl;

}


// method to write a GMAT Propagate Command
void gmatscripter::aux_GMAT_propagate(class gmatstep& agmatstep, bool useZeroPropagate){
	if (useZeroPropagate && agmatstep.zeroPropagate) { GMATfile << "	Propagate 'Propagate " << agmatstep.myspacecraft->Name << "' " << agmatstep.propagator.Name << "( " << agmatstep.myspacecraft->Name << " ) { " << agmatstep.myspacecraft->Name << ".ElapsedSecs = " << 0.0 << " }" << endl; }
	else {
		//declarations
		string str;
		double elapsed_secs = agmatstep.stepsize;

		if (agmatstep.myspacecraft->isForward) {
			if (elapsed_secs < 0.0) { elapsed_secs = -1.0*elapsed_secs; }
			str = "";
		}
		else {
			if (elapsed_secs > 0.0) { elapsed_secs = -1.0*elapsed_secs; }
			str = "BackProp";
		}

		GMATDebug << " str: " << str << endl;
		//write out the propagate command line
		if (agmatstep.allowTheTimeStep2Vary) {
			GMATfile << "	Propagate 'Propagate " << agmatstep.myspacecraft->Name << "' " << str << " " 
				<< agmatstep.propagator.Name << "( " << agmatstep.myspacecraft->Name << " ) { " 
				<< agmatstep.myspacecraft->Name << ".ElapsedSecs = " << agmatstep.id << "_TimeStep" << " }" << endl;
		}
		else {
			GMATfile << "	Propagate 'Propagate " << agmatstep.myspacecraft->Name << "' " << str << " " 
				<< agmatstep.propagator.Name << "( " << agmatstep.myspacecraft->Name << " ) { " << agmatstep.myspacecraft->Name << ".ElapsedSecs = " << elapsed_secs << " }" << endl;
		}
	}
}


// method to write a GMAT PenUp Command
void gmatscripter::PenUp(){

	GMATfile << "	PenUp 'PenUp' " << GMATMission.missionbodies[0].central_body_name << "View";
	for (int body_index = 0; body_index < GMATMission.missionbodies.size(); ++body_index) {
		GMATfile << " " << GMATMission.missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << endl;

}


// method to write a GMAT PenDown Command
void gmatscripter::PenDown(){

	GMATfile << "	PenDown 'PenDown' " << GMATMission.missionbodies[0].central_body_name << "View";
	for (int body_index = 0; body_index < GMATMission.missionbodies.size(); ++body_index) {
		GMATfile << " " << GMATMission.missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << endl;

}


// method to write a GMAT ForceModel Resource
void gmatscripter::create_GMAT_forcemodel(struct gmat_forcemodel& forcemodel){

	GMATfile << "Create ForceModel " << forcemodel.Name << endl;
	GMATfile << forcemodel.Name << ".CentralBody = " << forcemodel.CentralBody << endl;
	GMATfile << forcemodel.Name << ".PointMasses = {";
	for (int item = 0; item < forcemodel.PointMasses.size(); ++item) {
		if (item == forcemodel.PointMasses.size() - 1) { GMATfile << forcemodel.PointMasses[item]; }
		else { GMATfile << forcemodel.PointMasses[item] << ", "; }
	}
	GMATfile << "};" << endl;
	GMATfile << forcemodel.Name << ".Drag = None;" << endl;
	GMATfile << forcemodel.Name << ".SRP = Off;" << endl;
	GMATfile << endl;

}


// method to write a GMAT Propagator Resource
void gmatscripter::create_GMAT_propagator(struct gmat_propagator& propagator) {

	//declarations
	double initialstepsize = 60.0;
	double maxstepsize = 86400.0;

	if (propagator.isCloseApproach) {
		initialstepsize = 30.0;
		maxstepsize = 8640.0;
	}

	GMATfile << "Create Propagator " << propagator.Name << endl;
	GMATfile << propagator.Name << ".FM = " << propagator.ForceModel.Name << endl;
	GMATfile << propagator.Name << ".Type = PrinceDormand78; " << endl;
	GMATfile << propagator.Name << ".InitialStepSize = " << initialstepsize << endl;
	GMATfile << propagator.Name << ".Accuracy = 1e-11; " << endl;
	GMATfile << propagator.Name << ".MinStep = 0.0; " << endl;
	GMATfile << propagator.Name << ".MaxStep = " << maxstepsize << endl;
	GMATfile << endl;

}


// method to write a GMAT Spacecraft Resource
void gmatscripter::create_GMAT_spacecraft(struct gmat_spacecraft& spacecraft) {

	//write out the spacecraft information
	GMATfile << "Create Spacecraft " << spacecraft.Name << endl;
	GMATfile << spacecraft.Name << ".DateFormat = " << spacecraft.DateFormat << endl;
	GMATfile << spacecraft.Name << ".Epoch      = " << spacecraft.Epoch << endl;
	GMATfile << spacecraft.Name << ".DryMass    = " << spacecraft.DryMass << endl;
	GMATfile << spacecraft.Name << ".CoordinateSystem = " << spacecraft.CoordinateSystem << endl;
	GMATfile << spacecraft.Name << ".Tanks     = {" << spacecraft.Thruster.Tank.Name << "}" << endl;
	GMATfile << spacecraft.Name << ".Thrusters = {" << spacecraft.Thruster.Name << "}" << endl;
	GMATfile << endl;

}


// method to write a GMAT FuelTank Resource
void gmatscripter::create_GMAT_fueltank(struct gmat_spacecraft& spacecraft) {

	GMATfile << "Create FuelTank " << spacecraft.Thruster.Tank.Name << endl;
	GMATfile << spacecraft.Thruster.Tank.Name << ".AllowNegativeFuelMass = false" << endl;
	GMATfile << spacecraft.Thruster.Tank.Name << ".Volume = 10" << endl;
	GMATfile << spacecraft.Thruster.Tank.Name << ".FuelMass = " << spacecraft.Thruster.Tank.FuelMass << endl;
	GMATfile << endl;
	
}


// method to write a GMAT Thruster Resource
void gmatscripter::create_GMAT_thruster(struct gmat_spacecraft& spacecraft) {

	//GMATfile << "% Journey #" << j << ", Phase #" << p << ", Thruster Name: " << thrustername << endl;
	GMATfile << "Create Thruster " << spacecraft.Thruster.Name << endl;
	GMATfile << spacecraft.Thruster.Name << ".CoordinateSystem = " << spacecraft.CoordinateSystem << endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustDirection1 = 1;" << endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustDirection2 = 0;" << endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustDirection3 = 0;" << endl;
	GMATfile << spacecraft.Thruster.Name << ".DutyCycle = 1;" << endl;
	GMATfile << spacecraft.Thruster.Name << ".Tank = " << spacecraft.Thruster.Tank.Name << endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustScaleFactor = 1;" << endl;
	GMATfile << spacecraft.Thruster.Name << ".DecrementMass = true;" << endl;
	GMATfile << spacecraft.Thruster.Name << ".C1 = .1;" << endl;
	GMATfile << spacecraft.Thruster.Name << ".K1 = 3000;" << endl;
	GMATfile << endl;

}


// method to write a GMAT FiniteBurn Resource
void gmatscripter::create_GMAT_burn(struct gmat_burn& burn) {

	if (burn.Type == "FiniteBurn") {
		GMATfile << "Create FiniteBurn " << burn.Name << endl;
		GMATfile << burn.Name << ".Thrusters = {" << burn.ThrusterName << "};" << endl;
	}

}


// method to write a GMAT CoordinateSystem Resource
void gmatscripter::create_GMAT_coordinatesystem(string bodyname) {

	GMATfile << "Create CoordinateSystem " << bodyname << "J2000Eq;" << endl;
	GMATfile << bodyname << "J2000Eq.Origin = " << bodyname << ";" << endl;
	GMATfile << bodyname << "J2000Eq.Axes = MJ2000Eq;" << endl;
	GMATfile << endl;

}


void gmatscripter::create_GMAT_initialconditions(struct gmat_spacecraft& spacecraft) {
	GMATfile << "   " << spacecraft.Name << "." << spacecraft.CoordinateSystem << ".X  = " << spacecraft.initialconditions[0] << endl;
	GMATfile << "   " << spacecraft.Name << "." << spacecraft.CoordinateSystem << ".Y  = " << spacecraft.initialconditions[1] << endl;
	GMATfile << "   " << spacecraft.Name << "." << spacecraft.CoordinateSystem << ".Z  = " << spacecraft.initialconditions[2] << endl;
	GMATfile << "   " << spacecraft.Name << "." << spacecraft.CoordinateSystem << ".VX = " << spacecraft.initialconditions[3] << endl;
	GMATfile << "   " << spacecraft.Name << "." << spacecraft.CoordinateSystem << ".VY = " << spacecraft.initialconditions[4] << endl;
	GMATfile << "   " << spacecraft.Name << "." << spacecraft.CoordinateSystem << ".VZ = " << spacecraft.initialconditions[5] << endl;
	GMATfile << endl;
}


// ------------------------- //
// ----- 'gmatmission' ----- //
// ------------------------- //

//constructor
gmatmission::gmatmission(){}

//constructor
gmatmission::gmatmission(mission* anemtgmission) {
	emtgmission = anemtgmission;
	this->get_mission_bodies();
	this->setMissionThrustType();
}

//destructor
gmatmission::~gmatmission(){}


// --------------------------- //
// ------ 'gmatjourney' ------ //
// --------------------------- //

//constructor
gmatjourney::gmatjourney(){}

//constructor
gmatjourney::gmatjourney(gmatmission* amission, int journey) {
	mymission = amission;
	j = journey;
	id = "j" + std::to_string(journey);
	number_of_emtg_phases = 1; // mymission->emtgmission->journeys[j].number_of_phases;
}

//destructor
gmatjourney::~gmatjourney(){}


// ------------------------- //
// ------ 'gmatphase' ------ //
// ------------------------- //

//method
void gmatphase::append_step(gmatstep agmatstep) {
	//push_back or insert the new 'gmatstep'
	if (agmatstep.usePushBack) { this->mysteps.push_back(agmatstep); }
	else { this->mysteps.insert(this->mysteps.begin() + agmatstep.theInsertIndex, agmatstep); }

	//set the zeroPropagate member of 'gmatstep' to false if it is not the first or last 'gmatstep' of this 'gmatphase'
	if (agmatstep.gs > 1) {
		//if I am at least the third element, then turn the last guy's 'zeroPropagate' to 'false' because
		//he isn't the last element anymore
		for (int index = 0; index < this->mysteps.size(); ++index) {
			if (this->mysteps[index].gs == agmatstep.gs - 1) {
				this->mysteps[index].zeroPropagate = false;
			}
		}
	}
}

//method
void gmatphase::get_time_allocation() {
	for (int index = 0; index < this->mysteps.size(); ++index) {
		if (this->mysteps[index].allowTheTimeStep2Vary) { this->eligabletime += this->mysteps[index].stepsize; }
		else { this->ineligabletime += this->mysteps[index].stepsize; }
		this->TOF += this->mysteps[index].stepsize;
	}
}


} // end of EMTG namespace