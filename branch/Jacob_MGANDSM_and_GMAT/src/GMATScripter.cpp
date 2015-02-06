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
#include "Astrodynamics.h"

#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"

#include "SpiceUsr.h"

#include <string>  //string
#include <fstream> //ofstream


namespace EMTG {

// default constructor; not intended for use
gmatscripter::gmatscripter(){}

// constructor
gmatscripter::gmatscripter(mission* mission_in) : GMATMission(mission_in){
	//instantiate a 'gmatmissison' and assign a point to the 'EMTG mission class'
	this->ptr_gmatmission = mission_in;
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
	this->write_GMAT_initialconditions();

	// write out the pre-optimization calculations
	this->write_GMAT_preoptimization_calculations();

	// write out the optimization phase of the mission sequence
	this->write_GMAT_optimization();

}


// method to create a GMAT script file
void gmatscripter::create_GMAT_file(){

	//create a filename
	string filename = this->ptr_gmatmission->options.working_directory + "//" + 
					  this->ptr_gmatmission->options.mission_name + "_" + 
					  this->ptr_gmatmission->options.description + "_GMAT.script";
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
	//find the Xdescription index related to the optimal launch date
	for (size_t iX = 0; iX < this->ptr_gmatmission->Xdescriptions.size(); ++iX) {
		if (this->ptr_gmatmission->Xdescriptions[iX] == "j0p0: launch epoch (MJD)") {
			LaunchDate_LowerBounds = (ptr_gmatmission->Xlowerbounds[iX] / 86400.0) + TAIModJOffset;
			LaunchDate_UpperBounds = (ptr_gmatmission->Xupperbounds[iX] / 86400.0) + TAIModJOffset;
			LaunchDate = (ptr_gmatmission->Xopt[iX] / 86400.0) + TAIModJOffset;
			break;
		}
	}
	//assignments
	LaunchWindow = LaunchDate_UpperBounds - LaunchDate_LowerBounds;
	LaunchScaling = (LaunchDate - LaunchDate_LowerBounds) / LaunchWindow;


	// ---------------------------
	//      VARIABLE CREATION
	// ---------------------------
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
	GMATDebug << " ------------ " << std::endl;
	GMATDebug << "LaunchDate_LowerBounds: " << LaunchDate_LowerBounds << std::endl;
	GMATDebug << "LaunchDate_UpperBounds: " << LaunchDate_UpperBounds << std::endl;
	GMATDebug << "LaunchDate: " << LaunchDate << std::endl;
	GMATDebug << " ------------ " << std::endl;

}//end of get_GMAT_missionlevelparameters() method


// method to get the journey level parameters
void gmatscripter::create_GMAT_journeys() {

	GMATDebug << std::endl;
	GMATDebug << "create_GMAT_journeys()" << std::endl;
	GMATDebug << std::endl; 

	for (int j = 0; j < this->ptr_gmatmission->options.number_of_journeys; ++j) {
		
		// -----------------------------
		//         INSTANTIATION
		// -----------------------------
		gmatjourney agmatjourney(&GMATMission, j);

		// -----------------------------
		//       VARIABLE CREATION
		// -----------------------------


		// -----------------------------
		//         VARY CREATION
		// -----------------------------


		// -----------------------------
		//      CALCULATE CREATION
		// -----------------------------


		// -----------------------------
		//      CONSTRAINT CREATION
		// -----------------------------
		// journey time bounds (0: unbounded, 1: bounded flight time, 2: bounded arrival date)
		if (this->ptr_gmatmission->options.journey_timebounded[j] == 0) {
			
		}


		// -----------------------------
		//          PUSH_BACK
		// -----------------------------
		GMATMission.myjourneys.push_back(agmatjourney);

	}
}


// method to get the phase level parameters
void gmatscripter::create_GMAT_phases() {

	GMATDebug << std::endl;
	GMATDebug << "create_GMAT_phases()" << std::endl;
	GMATDebug << std::endl;

	//  declarations
	double fuelwindow;
	double fuellowerbound = 0.0;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].number_of_emtg_phases; ++p) {
			// -----------------------------
			//        DATA COLLECTION
			// -----------------------------
			//  fuellowerbound
			fuellowerbound = 0.0;
			//  fuel window
			if (this->ptr_gmatmission->options.enable_maximum_propellant_mass_constraint == 1) { 
				fuelwindow = this->ptr_gmatmission->options.maximum_propellant_mass; 
			}
			else { 
                fuelwindow = this->ptr_gmatmission->options.maximum_mass - this->ptr_gmatmission->options.final_mass_constraint;
			}
			//  dla window and lower bound
			double dlawindow = (this->ptr_gmatmission->options.DLA_bounds[1] - this->ptr_gmatmission->options.DLA_bounds[0]) * math::PI / 180.0;
			double dlalowerbound = this->ptr_gmatmission->options.DLA_bounds[0] * math::PI / 180.0;
			//  departure type
			int depart_type = this->ptr_gmatmission->options.journey_departure_type[j];
			int lv_type = this->ptr_gmatmission->options.LV_type;
			//  arrival type
			int arrival_type = this->ptr_gmatmission->options.journey_arrival_type[j];


			// -----------------------------
			//         INSTANTIATION
			// -----------------------------
			gmatphase agmatphase(&GMATMission.myjourneys[j], p);



			// -----------------------------
			//       VARIABLE CREATION
			// -----------------------------
			//allow the initial mass to vary
			if (j == 0 && agmatphase.p == 0) {
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_FuelMass"); // , agmatphase.spacecraft_forward.Thruster.Tank.FuelMass);
				//  if an impulsive maneuver is being performed, then the initial fuel mass may be increased to account for the burn
				//+ only for forward impulsive (an annoying fact of doing impulsive with GMAT and backwards time. note: EDS forward doesn't count either).
				double Ffuelmass = agmatphase.spacecraft_forward.Thruster.Tank.FuelMass;
				if (agmatphase.spacecraft_forward.iBurn.UseImpulsive && !agmatphase.spacecraft_forward.iBurn.IsEDS) {
					Ffuelmass = Ffuelmass*std::exp(sqrt(this->ptr_gmatmission->journeys[agmatphase.myjourney->j].phases[0].C3_departure) / agmatphase.spacecraft_forward.iBurn.c);
					agmatphase.spacecraft_forward.Thruster.Tank.FuelMass = Ffuelmass;
				}
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_FuelScaling", Ffuelmass / fuelwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_FuelWindow", fuelwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_FuelLowerBound", fuellowerbound);
			}
			//set the backward spacecraft fuelmass as a variable to be varied by the gmat optimizer
			agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_FuelMass");
			//  if an impulsive maneuver is being performed, then the initial fuel mass may need to be increased to account for the burn
			//+ only for forward impulsive (an annoying fact of doing impulsive with GMAT and backwards time. note: EDS forward doesn't count either).
			double Bfuelmass = agmatphase.spacecraft_backward.Thruster.Tank.FuelMass;
			//if (agmatphase.spacecraft_backward.iBurn.UseImpulsive) {
			//	Bfuelmass = Bfuelmass*std::exp(sqrt(this->ptr_gmatmission->journeys[agmatphase.myjourney->j].phases[0].C3_arrival) / agmatphase.spacecraft_backward.iBurn.c);
			//	agmatphase.spacecraft_backward.Thruster.Tank.FuelMass = Bfuelmass;
			//}
			agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_FuelScaling", Bfuelmass / fuelwindow);
			agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_FuelWindow", fuelwindow);
			agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_FuelLowerBound", fuellowerbound);

			//  enforce left boundary conditions
			if (depart_type == 0 && j == 0 && p == 0) { agmatphase.reprint_forward_position_IC = true; }
			//  enforce right boundary conditions
			if ((arrival_type == 1 || arrival_type == 3) && agmatphase.isLastPhase) {
				agmatphase.reprint_backward_position_IC = true;
				agmatphase.reprint_backward_velocity_IC = true;
			}

			//  variation of launch vinf
			//  if using the EDS motor or chemical propulsion (i.e. not LV), to leave from a parking orbit
			if (agmatphase.isFirstPhase && agmatphase.spacecraft_forward.iBurn.UseImpulsive && depart_type == 1) {

				//  create a matrix for the vinf outbound vector
				double dvdeparture[3] = { this->ptr_gmatmission->journeys[0].phases[0].dVdeparture[0],
					this->ptr_gmatmission->journeys[0].phases[0].dVdeparture[1],
					this->ptr_gmatmission->journeys[0].phases[0].dVdeparture[2] };
				math::Matrix<double> Vinf_out(3, 1, dvdeparture);

				//  compute the r0, v0, dv vectors
				math::Matrix<double> periapse_state_vector;
				std::vector<double> v0(3, 0.0);
				periapse_state_vector = this->ptr_gmatmission->journeys[j].phases[p].calculate_periapse_state_from_asymptote_and_parking_orbit(Vinf_out, v0,
					this->ptr_gmatmission->options.parking_orbit_inclination * math::PI / 180.0,
					this->ptr_gmatmission->options.parking_orbit_altitude,
					this->ptr_gmatmission->journeys[j].phases[p].phase_start_epoch,
					&this->ptr_gmatmission->TheUniverse[this->ptr_gmatmission->options.number_of_journeys - 1],
					agmatphase.mybodies[0]);

				//  overwrite the spacecraft position initial conditions
				for (size_t i = 0; i < 3; i++) {
					agmatphase.spacecraft_forward.initialconditions[i] = periapse_state_vector(i);
				}

				//  overwrite the spacecraft velocity initial conditions
				for (size_t i = 3; i < 6; i++) {
					agmatphase.spacecraft_forward.initialconditions[i] = v0[i - 3];
				}

				//  set the impulsive burn for the spacecraft to the correct delta-v vector
				agmatphase.spacecraft_forward.iBurn.Element1 = periapse_state_vector(3) - v0[0];
				agmatphase.spacecraft_forward.iBurn.Element2 = periapse_state_vector(4) - v0[1];
				agmatphase.spacecraft_forward.iBurn.Element3 = periapse_state_vector(5) - v0[2];

				//  set the dv lowerbound and window
				double dvlowerbound, dvwindow;
				//  using EDS
				if (agmatphase.spacecraft_forward.iBurn.IsEDS) {
					dvlowerbound = -this->ptr_gmatmission->options.journey_initial_impulse_bounds[j][1];
					dvwindow = 2.0 * this->ptr_gmatmission->options.journey_initial_impulse_bounds[j][1];
				}
				//  using Chemical System
				else {
					//  window; slighty perturb the 'fuellowerbound' in log() so that we don't use 0.0. 
					//+ really, if the 'minimum_dry_mass' is set to non-zero for a realistic mission, then this won't be a problem
					double _perturb = 0.01 * this->ptr_gmatmission->options.maximum_mass;
                    if (this->ptr_gmatmission->options.final_mass_constraint > _perturb) { _perturb = 0.0; }
                    dvwindow = -2.0 * agmatphase.spacecraft_forward.iBurn.c * std::log((this->ptr_gmatmission->options.final_mass_constraint + fuellowerbound + _perturb) / (fuelwindow + fuellowerbound));
					//  lower bound
					dvlowerbound = -dvwindow / 2.0;
				}

				//  dVX arrival
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVX");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVXScaling", (agmatphase.spacecraft_forward.iBurn.Element1 - dvlowerbound) / dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVXWindow", dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVXLowerBound", dvlowerbound);
				//  dVY arrival
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVY");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVYScaling", (agmatphase.spacecraft_forward.iBurn.Element2 - dvlowerbound) / dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVYWindow", dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVYLowerBound", dvlowerbound);
				//  dVZ arrival
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVZ");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVZScaling", (agmatphase.spacecraft_forward.iBurn.Element3 - dvlowerbound) / dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVZWindow", dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_dVZLowerBound", dvlowerbound);

				if (agmatphase.spacecraft_forward.iBurn.IsEDS) {
					//  dV norm; will need to constrain within bounds
					agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_EDS_dV");
				}

				//  parameters for varying 'f' (true anomaly) and 'RAAN' (right ascension of the ascending node).
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_TA");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_TAScaling");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_TAWindow", 360.0);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_TALowerBound", 0.0);
				agmatphase.setPreOptimizationCalculation(agmatphase.spacecraft_forward.Name + "_TAScaling", agmatphase.spacecraft_forward.Name + ".TA / " + agmatphase.spacecraft_forward.Name + "_TAWindow");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RAAN");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RAANScaling");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RAANWindow", 360.0);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RAANLowerBound", 0.0);
				agmatphase.setPreOptimizationCalculation(agmatphase.spacecraft_forward.Name + "_RAANScaling", agmatphase.spacecraft_forward.Name + ".RAAN / " + agmatphase.spacecraft_forward.Name + "_RAANWindow");

			}
			//  use of a non-parking orbit departure
			else if (agmatphase.isFirstPhase && depart_type == 0) {

				//RLA
				double RLA = this->ptr_gmatmission->journeys[agmatphase.myjourney->j].phases[agmatphase.p].RA_departure;
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RLA");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RLAScaling", RLA / (4 * math::PI));
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RLAWindow", 4 * math::PI);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_RLALowerBound", -math::TwoPI);
				//DLA
				double DLA = this->ptr_gmatmission->journeys[agmatphase.myjourney->j].phases[agmatphase.p].DEC_departure;
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_DLA");
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_DLAScaling", (DLA - dlalowerbound) / dlawindow);
				if (agmatphase.myjourney->j == 0) {
					agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_DLAWindow", dlawindow);
					agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_DLALowerBound", dlalowerbound);
				}
				else { //  this is the bounds on the DLA?
					agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_DLAWindow", math::PI);
					agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_DLALowerBound", -math::PIover2);
				}
				//Vinf
				double vinf = sqrt(this->ptr_gmatmission->journeys[agmatphase.myjourney->j].phases[0].C3_departure);
				double vinfwindow = this->ptr_gmatmission->options.journey_initial_impulse_bounds[agmatphase.myjourney->j][1] - this->ptr_gmatmission->options.journey_initial_impulse_bounds[agmatphase.myjourney->j][0];

				//  "EDS", non-parking orbit
				if (agmatphase.spacecraft_forward.iBurn.UseImpulsive && depart_type == 0) {

					agmatphase.spacecraft_forward.iBurn.Element1 = this->ptr_gmatmission->journeys[0].phases[0].dVdeparture[0];
					agmatphase.spacecraft_forward.initialconditions[3] -= agmatphase.spacecraft_forward.iBurn.Element1;
					agmatphase.spacecraft_forward.iBurn.Element2 = this->ptr_gmatmission->journeys[0].phases[0].dVdeparture[1];
					agmatphase.spacecraft_forward.initialconditions[4] -= agmatphase.spacecraft_forward.iBurn.Element2;
					agmatphase.spacecraft_forward.iBurn.Element3 = this->ptr_gmatmission->journeys[0].phases[0].dVdeparture[2];
					agmatphase.spacecraft_forward.initialconditions[5] -= agmatphase.spacecraft_forward.iBurn.Element3;

				}
				//  compute the initial velocity magnitude. this is used to see if we need an extra boost at the start.
				double initial_velocity_magnitude = 0.0;
				for (size_t i = 3; i < 6; i++) { initial_velocity_magnitude += std::pow(agmatphase.spacecraft_forward.initialconditions[i], 2); }
				initial_velocity_magnitude = std::sqrt(initial_velocity_magnitude);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_Vinf");
				//  perturb the _vinf initial scaling if it is deemed to small with the current initial velocity
				//+ (i.e. the LT spacecraft may fall back into the planet after departure...)
				double C3;
				if (vinf < 0.1 && initial_velocity_magnitude < 0.1) { 
					C3 = 0.1*0.1;
					agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_VinfScaling", 0.1 / vinfwindow); 
				}
				else {
					C3 = vinf*vinf;
					agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_VinfScaling", vinf / vinfwindow); 
				}
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_VinfWindow", vinfwindow);
				agmatphase.setVariable(agmatphase.spacecraft_forward.Name + "_VinfLowerBound", this->ptr_gmatmission->options.journey_initial_impulse_bounds[agmatphase.myjourney->j][0]);
				//  if a launch vehicle is being used on the first journey and it isn't (-1 or 0) 'i.e. burn with EDS or Fixed Initial Mass'
				if (j == 0 && (lv_type != -1 && lv_type != 0)) {

					int k = lv_type - 1;

					double a1[] = { 0, 0, -8.00E-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.5649E-08, 0, 8.0e-7, 0, -3.5649E-08, 0, 8.0e-7, 0 };
					double a2[] = { 2.81E-05, 5.62E-06, -3.88E-05, -2.48E-05, -0.00009382, 0.00006198, 0.00012334, -0.00010901, 0.00006465, 0.00007576, -0.00138695, -0.00004196, 2.31305E-05, 0, -0.0002, 0 };
					double a3[] = { -0.001517208, -0.00067462, 0.00190765, 0.00211196, 0.00403555, -0.00295026, -0.00712978, 0.00596291, -0.00276635, -0.00405132, 0.03050117, -0.0035373, -0.006458248, 0, 0.0145, 0 };
					double a4[] = { 0.357956243, 0.44268134, 0.43698409, 0.47075449, 0.26854604, 0.41398944, 0.61102598, 0.37650144, 0.60963424, 0.70733066, 0.30110723, 0.91375291, 1.065903806, 0, 0.5825, 0 };
					double a5[] = { -64.48375826, -78.8465652, -88.38856438, -98.4962944, -55.39915501, -69.35547443, -83.52984026, -91.90752777, -102.9890546, -111.2601399, -69.96585082, -110.0132867, -114.5157292, 0, -183.69, 0 };
					double a6[] = { 3034.683258, 3930.871041, 4655.158371, 5235.260181, 2094.89, 3265.42, 4195.86, 4939.98, 5595.09, 6106.14, 1974.88, 3634.59, 6595.738063, 0, 12170, 0 };
					double C3min[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
					double C3max[] = { 60, 60, 60, 60, 50, 60, 60, 60, 60, 60, 35, 40, 200, 0, 100, 0 };

					//  if a custom LV type, then overwrite the local array entry that will be called at the end 
					//+ of this if() 
					if (lv_type == -2) {
						a1[k] = this->ptr_gmatmission->options.custom_LV_coefficients[0];
						a2[k] = this->ptr_gmatmission->options.custom_LV_coefficients[1];
						a3[k] = this->ptr_gmatmission->options.custom_LV_coefficients[2];
						a4[k] = this->ptr_gmatmission->options.custom_LV_coefficients[3];
						a5[k] = this->ptr_gmatmission->options.custom_LV_coefficients[4];
						a6[k] = this->ptr_gmatmission->options.custom_LV_coefficients[5];
						C3min[k] = this->ptr_gmatmission->options.custom_LV_C3_bounds[0];
						C3max[k] = this->ptr_gmatmission->options.custom_LV_C3_bounds[1];
					}
					
					double mass, dmdC3;
					EMTG::Astrodynamics::find_mass_to_orbit(C3, DLA, lv_type, &mass, &dmdC3, &this->ptr_gmatmission->options);

					agmatphase.setVariable("LV_MassCapability", this->ptr_gmatmission->options.LV_margin * mass);
					agmatphase.setVariable("LV_Margin", this->ptr_gmatmission->options.LV_margin);
					agmatphase.setVariable("LV_AdapterMass", this->ptr_gmatmission->options.LV_adapter_mass);
					agmatphase.setVariable("LV_a1", a1[k]);
					agmatphase.setVariable("LV_a2", a2[k]);
					agmatphase.setVariable("LV_a3", a3[k]);
					agmatphase.setVariable("LV_a4", a4[k]);
					agmatphase.setVariable("LV_a5", a5[k]);
					agmatphase.setVariable("LV_a6", a6[k]);
					agmatphase.setVariable("C31", C3);
					agmatphase.setVariable("C32", std::pow(C3, 2));
					agmatphase.setVariable("C33", std::pow(C3, 3));
					agmatphase.setVariable("C34", std::pow(C3, 4));
					agmatphase.setVariable("C35", std::pow(C3, 5));
					agmatphase.setVariable(agmatphase.id + "_InitialMass");
				}
			}

			//variation of chemical arrival
			if (agmatphase.isLastPhase && (arrival_type == 0 || arrival_type == 1)) {
				//  the dv arrival vector and magnitude
				double dvarrival[3] = { -this->ptr_gmatmission->journeys[0].phases[0].dVarrival[0],
										-this->ptr_gmatmission->journeys[0].phases[0].dVarrival[1],
										-this->ptr_gmatmission->journeys[0].phases[0].dVarrival[2] };
				//double dvarrival_mag = this->ptr_gmatmission->journeys[0].phases[0].dV_arrival_magnitude;
				//  set the iBurn and initial conditions for the spacecraft
				if (agmatphase.spacecraft_backward.iBurn.UseImpulsive) {
					agmatphase.spacecraft_backward.iBurn.Element1 = dvarrival[0];
					agmatphase.spacecraft_backward.initialconditions[3] -= agmatphase.spacecraft_backward.iBurn.Element1;
					agmatphase.spacecraft_backward.iBurn.Element2 = dvarrival[1];
					agmatphase.spacecraft_backward.initialconditions[4] -= agmatphase.spacecraft_backward.iBurn.Element2;
					agmatphase.spacecraft_backward.iBurn.Element3 = dvarrival[2];
					agmatphase.spacecraft_backward.initialconditions[5] -= agmatphase.spacecraft_backward.iBurn.Element3;
				}
				//  compute the initial velocity magnitude. this is used to see if we need an extra boost at the start.
				double _initial_velocity_magnitude = 0.0;
				for (size_t i = 3; i < 6; i++) { _initial_velocity_magnitude += std::pow(agmatphase.spacecraft_backward.initialconditions[i], 2); }
				_initial_velocity_magnitude = std::sqrt(_initial_velocity_magnitude);
				//  window; slighty perturb the 'fuellowerbound' in log() so that we don't use 0.0. 
				//+ really, if the 'minimum_dry_mass' is set to non-zero for a realistic mission, then this won't be a problem
				double _perturb = 0.01 * this->ptr_gmatmission->options.maximum_mass;
                if (this->ptr_gmatmission->options.final_mass_constraint > _perturb) { _perturb = 0.0; }
                double dvwindow = -2.0 * agmatphase.spacecraft_backward.iBurn.c * std::log((this->ptr_gmatmission->options.final_mass_constraint + fuellowerbound + _perturb) / (fuelwindow + fuellowerbound));
				//  lower bound
				double dvlowerbound = -dvwindow / 2.0;
				//  dVX arrival
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVX");
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVXScaling", (dvarrival[0] - dvlowerbound) / dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVXWindow", dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVXLowerBound", dvlowerbound);
				//  dVY arrival
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVY");
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVYScaling", (dvarrival[1] - dvlowerbound) / dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVYWindow", dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVYLowerBound", dvlowerbound);
				//  dVZ arrival
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVZ");
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVZScaling", (dvarrival[2] - dvlowerbound) / dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVZWindow", dvwindow);
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dVZLowerBound", dvlowerbound);
				//  dV norm for calculating fuel burn
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_dV");
				//  'c' for calculating fuel burn
				agmatphase.setVariable(agmatphase.spacecraft_backward.Name + "_c", agmatphase.spacecraft_backward.iBurn.c);
			}



			// -----------------------------
			//         VARY CREATION
			// -----------------------------
			//backward spacecraft fuel scaling
			agmatphase.setVary(agmatphase.spacecraft_backward.Name + "_FuelScaling");
			//forward spacecraft fuel scaling (1st journey only)
			if (j == 0 && agmatphase.p == 0 && this->ptr_gmatmission->options.allow_initial_mass_to_vary) {
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_FuelScaling");
			}
			//forward spacecraft at beginning of journey, launch scalings
			if (agmatphase.isFirstPhase && this->ptr_gmatmission->options.journey_departure_type[agmatphase.myjourney->j] == 0) {
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_RLAScaling");
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_DLAScaling");
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_VinfScaling");
			}
			if (agmatphase.isFirstPhase && agmatphase.spacecraft_forward.iBurn.UseImpulsive && depart_type == 1) {
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_dVXScaling");
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_dVYScaling");
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_dVZScaling");
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_TAScaling");
				agmatphase.setVary(agmatphase.spacecraft_forward.Name + "_RAANScaling");
			}
			//backward spacecraft chemical arrivals
			if (agmatphase.isLastPhase && (arrival_type == 0 || arrival_type == 1)) {
				agmatphase.setVary(agmatphase.spacecraft_backward.Name + "_dVXScaling");
				agmatphase.setVary(agmatphase.spacecraft_backward.Name + "_dVYScaling");
				agmatphase.setVary(agmatphase.spacecraft_backward.Name + "_dVZScaling");
			}



			// -----------------------------
			//      CALCULATE CREATION
			// -----------------------------
			//spacecraft drymass
			agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + ".DryMass", agmatphase.myjourney->id + "_DryMass");
			agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + ".DryMass", agmatphase.myjourney->id + "_DryMass");
			//forward spacecraft epoch 
			//	(first phase can vary using 'LaunchWindowScaling'). 
			//	otherwise phases use variable create from previous phase's backward spacecraft epoch [see postpass_GMAT_phases() method]
			if (agmatphase.myjourney->j == 0 && agmatphase.p == 0) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + ".Epoch." + agmatphase.spacecraft_forward.DateFormat,
										"(LaunchWindowScaling * LaunchWindow + LaunchWindowOpenDate)");
			}
			else if (agmatphase.myjourney->j != 0 && agmatphase.p == 0) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + ".Epoch." + agmatphase.spacecraft_forward.DateFormat,
					agmatphase.myjourney->id + "_WaitTimeScaling * " + agmatphase.myjourney->id + "_WaitTimeWindow + " + GMATMission.myjourneys[j - 1].myphases.back().spacecraft_backward.Name + "_Epoch + " + 
					agmatphase.myjourney->id + "_WaitTimeLowerBound");
			}
			else if (agmatphase.p > 0) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + ".Epoch." + agmatphase.spacecraft_forward.DateFormat, GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_Epoch");
			}
			//forward spacecraft fuel scaling
			//  for p == 0 and j == 0 the forward spacecraft can vary its fuel if the corresponding emtg option was flagged true.
			//	for p == 0 but j != 0 the forward spacecraft adopts the backward spacecraft from the previous journey and may have mass increments/decrements.
			//	for p > 0, the forward spacecraft adopts the backward spacecraft choice from the previous 'gmatstep'.
			if (j == 0 && agmatphase.p == 0 && this->ptr_gmatmission->options.allow_initial_mass_to_vary) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_FuelMass",
										agmatphase.spacecraft_forward.Name + "_FuelScaling * " + agmatphase.spacecraft_forward.Name + "_FuelWindow" + " + " + agmatphase.spacecraft_forward.Name + "_FuelLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.Thruster.Tank.Name + ".FuelMass", agmatphase.spacecraft_forward.Name + "_FuelMass");
			}
			else if (j != 0 && agmatphase.p == 0) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.Thruster.Tank.Name + ".FuelMass", GMATMission.myjourneys[j - 1].myphases.back().spacecraft_backward.Name + "_FuelMass");
			}
			else if (p > 0) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.Thruster.Tank.Name + ".FuelMass", GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_FuelMass");
			}
			//backward spacecraft fuel scaling
			agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "_FuelMass",
									agmatphase.spacecraft_backward.Name + "_FuelScaling * " + agmatphase.spacecraft_backward.Name + "_FuelWindow" + " + " + agmatphase.spacecraft_backward.Name + "_FuelLowerBound");
			agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "." + agmatphase.spacecraft_backward.Thruster.Tank.Name + ".FuelMass", agmatphase.spacecraft_backward.Name + "_FuelMass");
			//forward spacecraft X, Y, Z, VX, VY, VZ (after flyby phase)
			if (p > 0 && agmatphase.StartsWithFlyby) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".X", GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_X");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".Y", GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_Y");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".Z", GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_Z");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".VX", GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_VX");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".VY", GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_VY");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".VZ", GMATMission.myjourneys[j].myphases[p - 1].spacecraft_backward.Name + "_VZ");
			}
			//for propulsive departures
			else if (agmatphase.isFirstPhase && depart_type == 0) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_DLA", agmatphase.spacecraft_forward.Name + "_DLAScaling * " + agmatphase.spacecraft_forward.Name + "_DLAWindow + " + agmatphase.spacecraft_forward.Name + "_DLALowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_RLA", agmatphase.spacecraft_forward.Name + "_RLAScaling * " + agmatphase.spacecraft_forward.Name + "_RLAWindow + " + agmatphase.spacecraft_forward.Name + "_RLALowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_Vinf", agmatphase.spacecraft_forward.Name + "_VinfScaling * " + agmatphase.spacecraft_forward.Name + "_VinfWindow + " + agmatphase.spacecraft_forward.Name + "_VinfLowerBound");
				if (agmatphase.spacecraft_forward.iBurn.UseImpulsive) {
					agmatphase.setCalculate(agmatphase.spacecraft_forward.iBurn.Name + ".Element1", agmatphase.spacecraft_forward.Name + "_Vinf * cos( " + agmatphase.spacecraft_forward.Name + "_RLA ) * cos( " + agmatphase.spacecraft_forward.Name + "_DLA )");
					agmatphase.setCalculate(agmatphase.spacecraft_forward.iBurn.Name + ".Element2", agmatphase.spacecraft_forward.Name + "_Vinf * sin( " + agmatphase.spacecraft_forward.Name + "_RLA ) * cos( " + agmatphase.spacecraft_forward.Name + "_DLA )");
					agmatphase.setCalculate(agmatphase.spacecraft_forward.iBurn.Name + ".Element3", agmatphase.spacecraft_forward.Name + "_Vinf * sin( " + agmatphase.spacecraft_forward.Name + "_DLA )");
				}
				else {
					agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".VX", agmatphase.spacecraft_forward.Name + "_Vinf * cos( " + agmatphase.spacecraft_forward.Name + "_RLA ) * cos( " + agmatphase.spacecraft_forward.Name + "_DLA )");
					agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".VY", agmatphase.spacecraft_forward.Name + "_Vinf * sin( " + agmatphase.spacecraft_forward.Name + "_RLA ) * cos( " + agmatphase.spacecraft_forward.Name + "_DLA )");
					agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.CoordinateSystem + ".VZ", agmatphase.spacecraft_forward.Name + "_Vinf * sin( " + agmatphase.spacecraft_forward.Name + "_DLA )");
				}
			}
			//  departures from parking orbits
			else if (agmatphase.isFirstPhase && agmatphase.spacecraft_forward.iBurn.UseImpulsive && depart_type == 1) {
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_dVX", agmatphase.spacecraft_forward.Name + "_dVXScaling * " + agmatphase.spacecraft_forward.Name + "_dVXWindow + " + agmatphase.spacecraft_forward.Name + "_dVXLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_dVY", agmatphase.spacecraft_forward.Name + "_dVYScaling * " + agmatphase.spacecraft_forward.Name + "_dVYWindow + " + agmatphase.spacecraft_forward.Name + "_dVYLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_dVZ", agmatphase.spacecraft_forward.Name + "_dVZScaling * " + agmatphase.spacecraft_forward.Name + "_dVZWindow + " + agmatphase.spacecraft_forward.Name + "_dVZLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.iBurn.Name + ".Element1", agmatphase.spacecraft_forward.Name + "_dVX");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.iBurn.Name + ".Element2", agmatphase.spacecraft_forward.Name + "_dVY");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.iBurn.Name + ".Element3", agmatphase.spacecraft_forward.Name + "_dVZ");
				if (agmatphase.spacecraft_forward.iBurn.IsEDS) {
					agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_EDS_dV", "sqrt( " + agmatphase.spacecraft_forward.Name + "_dVX * " + agmatphase.spacecraft_forward.Name + "_dVX + " + agmatphase.spacecraft_forward.Name + "_dVY * " + agmatphase.spacecraft_forward.Name + "_dVY + " + agmatphase.spacecraft_forward.Name + "_dVZ * " + agmatphase.spacecraft_forward.Name + "_dVZ )");
				}
				//  'f' and 'RAAN'
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_TA", agmatphase.spacecraft_forward.Name + "_TAScaling * " + agmatphase.spacecraft_forward.Name + "_TAWindow + " + agmatphase.spacecraft_forward.Name + "_TALowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + ".TA", agmatphase.spacecraft_forward.Name + "_TA");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + "_RAAN", agmatphase.spacecraft_forward.Name + "_RAANScaling * " + agmatphase.spacecraft_forward.Name + "_RAANWindow + " + agmatphase.spacecraft_forward.Name + "_RAANLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_forward.Name + ".RAAN", agmatphase.spacecraft_forward.Name + "_RAAN");
			}
			//  C3 and Mass calculation for LV departures
			if (j == 0 && agmatphase.isFirstPhase && depart_type == 0 && (lv_type != -1 && lv_type != 0)) {
				//  calculate the C3 power terms
				agmatphase.setCalculate("C31", agmatphase.spacecraft_forward.Name + "_Vinf * " + agmatphase.spacecraft_forward.Name + "_Vinf");
				agmatphase.setCalculate("C32", "C31 * C31");
				agmatphase.setCalculate("C33", "C32 * C31");
				agmatphase.setCalculate("C34", "C33 * C31");
				agmatphase.setCalculate("C35", "C34 * C31");
				agmatphase.setCalculate("LV_MassCapability", "(1.0 - LV_Margin) * (LV_a1 * C35 + LV_a2 * C34 + LV_a3 * C33 + LV_a4 * C32 + LV_a5 * C31 + LV_a6 - LV_AdapterMass)");
				agmatphase.setCalculate(agmatphase.id + "_InitialMass", agmatphase.spacecraft_forward.Name + ".DryMass + " + agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.Thruster.Tank.Name + ".FuelMass");
			}

			//  chemical arrivals
			if (agmatphase.isLastPhase && (arrival_type == 0 || arrival_type == 1)) {
				agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "_dVX", agmatphase.spacecraft_backward.Name + "_dVXScaling * " + agmatphase.spacecraft_backward.Name + "_dVXWindow + " + agmatphase.spacecraft_backward.Name + "_dVXLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "_dVY", agmatphase.spacecraft_backward.Name + "_dVYScaling * " + agmatphase.spacecraft_backward.Name + "_dVYWindow + " + agmatphase.spacecraft_backward.Name + "_dVYLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "_dVZ", agmatphase.spacecraft_backward.Name + "_dVZScaling * " + agmatphase.spacecraft_backward.Name + "_dVZWindow + " + agmatphase.spacecraft_backward.Name + "_dVZLowerBound");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.iBurn.Name + ".Element1", agmatphase.spacecraft_backward.Name + "_dVX");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.iBurn.Name + ".Element2", agmatphase.spacecraft_backward.Name + "_dVY");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.iBurn.Name + ".Element3", agmatphase.spacecraft_backward.Name + "_dVZ");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "_dV", "sqrt( " + agmatphase.spacecraft_backward.Name + "_dVX * " + agmatphase.spacecraft_backward.Name + "_dVX + " + agmatphase.spacecraft_backward.Name + "_dVY * " + agmatphase.spacecraft_backward.Name + "_dVY + " + agmatphase.spacecraft_backward.Name + "_dVZ * " + agmatphase.spacecraft_backward.Name + "_dVZ )");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "_FuelMass", "( " + agmatphase.spacecraft_backward.Name + "_FuelMass + " + agmatphase.myjourney->id + "_DryMass ) * exp( " + agmatphase.spacecraft_backward.Name + "_dV / " + agmatphase.spacecraft_backward.Name + "_c ) - " + agmatphase.myjourney->id + "_DryMass");
				agmatphase.setCalculate(agmatphase.spacecraft_backward.Name + "." + agmatphase.spacecraft_backward.Thruster.Tank.Name + ".FuelMass", agmatphase.spacecraft_backward.Name + "_FuelMass");
			}


			// -----------------------------
			//      CONSTRAINT CREATION
			// -----------------------------
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " X", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.X", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.X");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " Y", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.Y", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.Y");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " Z", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.Z", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.Z");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " VX", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.VX", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.VX");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " VY", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.VY", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.VY");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " VZ", agmatphase.spacecraft_forward.Name + ".SunJ2000Eq.VZ", "=", agmatphase.spacecraft_backward.Name + ".SunJ2000Eq.VZ");
			agmatphase.setConstraint("MatchPoint " + agmatphase.id + " Mass", agmatphase.spacecraft_forward.Name + "." + agmatphase.spacecraft_forward.Thruster.Tank.Name + ".FuelMass", "=", agmatphase.spacecraft_backward.Name + "." + agmatphase.spacecraft_backward.Thruster.Tank.Name + ".FuelMass");
			//  constrain initial mass within LV capabilities
			if (j == 0 && agmatphase.isFirstPhase && depart_type == 0 && (lv_type != -1 && lv_type != 0)) {
				agmatphase.setConstraint("Launch Vehicle Mass Constraint", agmatphase.id + "_InitialMass", "<=", "LV_MassCapability");
			}
			//  constraint EDS burn within capabilities
			if (agmatphase.isFirstPhase && agmatphase.spacecraft_forward.iBurn.UseImpulsive && depart_type == 1 && agmatphase.spacecraft_forward.iBurn.IsEDS) {
				agmatphase.setConstraint("EDS Delta-V Constraint", agmatphase.spacecraft_forward.Name + "_EDS_dV", "<=", this->ptr_gmatmission->options.journey_initial_impulse_bounds[j][1]);
			}

			// -----------------------------
			//          PUSH_BACK
			// -----------------------------
			GMATMission.myjourneys[j].myphases.push_back(agmatphase);


		}//end of phases for-statement
	}//end of journeys for-statement

}


// method to get the step level parameters
void gmatscripter::create_GMAT_steps() {

	GMATDebug << std::endl;
	GMATDebug << "create_GMAT_steps()" << std::endl;
	GMATDebug << std::endl;

	//declarations
	int body_index;
	double periapse_velocity_magnitude;
	double initial_relative_velocity_norm;
	double final_relative_velocity_norm;
	double minimum_relative_velocity_norm;
	double approximate_time_in_SOI;
	
	// ---------------------------
	//      VARIABLE CREATION
	// ---------------------------

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int s = 0; s < ptr_gmatmission->options.num_timesteps; ++s) {

				//instantiate a new 'gmatstep'
				gmatstep agmatstep(&GMATMission.myjourneys[j].myphases[p], s);

				//Case: [ MP+ ] (i.e. the step just after the match-point); make this a half-step
				if (!agmatstep.myspacecraft->isForward && agmatstep.isMatchPointStep) {
					GMATDebug << "MP+     ";
					//set the ForceModel and Propagator Parameters (false: is not a CloseApproach)
					agmatstep.setFMandProp(false);
					//scale the stepsize to be one-half
					agmatstep.scale_stepsize(0.5);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}

				//Case: [ SOI-SOI ] (i.e. the step is complete contained within the SOI of its body for the entire step time)
				if (agmatstep.inSOIatStart && agmatstep.inSOIatEnd) {
					GMATDebug << "SOI-SOI ";
					//  set the ForceModel and Propagator Parameters (true: for CloseApproach, i.e. just use central body).
					//+ if leaving/arriving from/at a parking orbit, then use the larger propagator
					if (agmatstep.myspacecraft->iBurn.BCisParkingOrbit) { agmatstep.setFMandProp(false); } 
					else { agmatstep.setFMandProp(true); }
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
				//Case: [ SOI --> ] (i.e. the step starts within the SOI of its body, but exits it at somepoint during the step. this step may be cut in two)
				else if (agmatstep.inSOIatStart && !agmatstep.inSOIatEnd) {
					GMATDebug << "SOI --> ";
					//for this case we must figure out how long we are in the SOI before exiting
					//we start by calculating the norm of the initial and final relative velocity
					initial_relative_velocity_norm = math::norm(agmatstep.initial_velocity_diff, 3);
					final_relative_velocity_norm   = math::norm(agmatstep.final_velocity_diff,   3);
					//then for conservative reasons, select the minimum of the two norms just calculated
					if (initial_relative_velocity_norm <= final_relative_velocity_norm) { minimum_relative_velocity_norm = initial_relative_velocity_norm; }
					else { minimum_relative_velocity_norm = final_relative_velocity_norm; }
					//if the step is associated with the forward spacecraft
					if (agmatstep.myspacecraft->isForward) {
						body_index = 0;
						//if the step is part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force smaller max timesteps by the GMAT propagator)
						if (agmatstep.myphase->StartsWithFlyby) {
							periapse_velocity_magnitude = sqrt(2.0 * agmatstep.myphase->mybodies[body_index].mu / (agmatstep.myphase->mybodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + minimum_relative_velocity_norm*minimum_relative_velocity_norm);
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / periapse_velocity_magnitude);
						}
						//if the step is not part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force a longer 'gmatstep' timestep, such that the additional 3rd body perturbations don't get turned on to early.)
						else {
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / minimum_relative_velocity_norm);
						}
					}
					//the step is associated with the backward spacecraft (it is hard to imagine when this else block will run, but for completeness it has been added.)
					else { 
						body_index = 1; 
						//if the step is part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force smaller max timesteps by the GMAT propagator)
						if (agmatstep.myphase->EndsWithFlyby) {
							periapse_velocity_magnitude = sqrt(2.0 * agmatstep.myphase->mybodies[body_index].mu / (agmatstep.myphase->mybodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + minimum_relative_velocity_norm*minimum_relative_velocity_norm);
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / periapse_velocity_magnitude);
						}
						//if the step is not part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force a longer 'gmatstep' timestep, such that the additional 3rd body perturbations don't get turned on to early.)
						else {
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / minimum_relative_velocity_norm);
						}
					}
					//make sure the 'approximate_time_in_SOI' is not larger than the original 'stepsize'
					if (approximate_time_in_SOI > agmatstep.stepsize) { approximate_time_in_SOI = agmatstep.stepsize; }
					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(approximate_time_in_SOI);
					//  set the ForceModel and Propagator Parameters (true: for CloseApproach, i.e. just use central body).
					//+ if leaving/arriving from/at a parking orbit, then use the larger propagator
					if (agmatstep.myspacecraft->iBurn.BCisParkingOrbit) { agmatstep.setFMandProp(false); }
					else { agmatstep.setFMandProp(true); }
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();

					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(agmatstep.stepsize - approximate_time_in_SOI);
					//set the ForceModel and Propagator Parameters (false: is not a CloseApproach)
					agmatstep.setFMandProp(false);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
				//Case: [ --> SOI ] (i.e. the step starts outside of the SOI of its body, but eventually enters it. this step may be cut in two)
				else if (!agmatstep.inSOIatStart && agmatstep.inSOIatEnd) {
					GMATDebug << "-->SOI  ";
					//for this case we must figure out how long we are in the SOI before exiting
					//we start by calculating the norm of the initial and final relative velocity
					initial_relative_velocity_norm = math::norm(agmatstep.initial_velocity_diff, 3);
					final_relative_velocity_norm   = math::norm(agmatstep.final_velocity_diff,   3);
					//then for conservative reasons, select the minimum of the two norms just calculated
					if (initial_relative_velocity_norm <= final_relative_velocity_norm) { minimum_relative_velocity_norm = initial_relative_velocity_norm; }
					else { minimum_relative_velocity_norm = final_relative_velocity_norm; }
					//if the step is associated with the forward spacecraft (it is hard to imagine when this if block will run, but for completeness it has been added.)
					if (agmatstep.myspacecraft->isForward) {
						body_index = 0;
						//if the step is part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force smaller max timesteps by the GMAT propagator)
						if (agmatstep.myphase->StartsWithFlyby) {
							periapse_velocity_magnitude = sqrt(2.0 * agmatstep.myphase->mybodies[body_index].mu / (agmatstep.myphase->mybodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + minimum_relative_velocity_norm*minimum_relative_velocity_norm);
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / periapse_velocity_magnitude);
						}
						//if the step is not part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force a longer 'gmatstep' timestep, such that the additional 3rd body perturbations don't get turned on to early.)
						else {
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / minimum_relative_velocity_norm);
						}
					}
					//the step is associated with the backward spacecraft 
					else {
						body_index = 1;
						//if the step is part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force smaller max timesteps by the GMAT propagator)
						if (agmatstep.myphase->EndsWithFlyby) {
							periapse_velocity_magnitude = sqrt(2.0 * agmatstep.myphase->mybodies[body_index].mu / (agmatstep.myphase->mybodies[body_index].radius + this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude) + minimum_relative_velocity_norm*minimum_relative_velocity_norm);
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / periapse_velocity_magnitude);
						}
						//if the step is not part of a half-phase flyby. using minimum_relative_velocity_norm for conservative reasons (in this case it will force a longer 'gmatstep' timestep, such that the additional 3rd body perturbations don't get turned on to early.)
						else {
							approximate_time_in_SOI = (agmatstep.myphase->mybodies[body_index].r_SOI / minimum_relative_velocity_norm);
						}
					}
					//make sure the 'approximate_time_in_SOI' is not larger than the original 'stepsize'
					if (approximate_time_in_SOI > agmatstep.stepsize) { approximate_time_in_SOI = agmatstep.stepsize; }
					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(agmatstep.stepsize - approximate_time_in_SOI);
					//set the ForceModel and Propagator Parameters (false: is not a CloseApproach)
					agmatstep.setFMandProp(false);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();

					//reset the stepsize of the 'gmatstep'
					agmatstep.set_stepsize(approximate_time_in_SOI);
					//  set the ForceModel and Propagator Parameters (true: for CloseApproach, i.e. just use central body).
					//+ if leaving/arriving from/at a parking orbit, then use the larger propagator
					if (agmatstep.myspacecraft->iBurn.BCisParkingOrbit) { agmatstep.setFMandProp(false); }
					else { agmatstep.setFMandProp(true); }
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
				//Case: [ -----> ] (i.e. this step is traveling in "free-space" the entire time. said another way, it is never within a bodies SOI)
				else {
					GMATDebug << " -----> ";
					//set the ForceModel and Propagator Parameters (false: is not a CloseApproach)
					agmatstep.setFMandProp(false);
					//set a variable for these time steps in "free-space"
					agmatstep.allowTheTimeStep2Vary = true;
					agmatstep.setVariable(agmatstep.id + "_TimeStep");
					agmatstep.setVariable(agmatstep.id + "_InitialTimeStep", agmatstep.stepsize);
					if (agmatstep.myspacecraft->isForward) { agmatstep.setCalculate(agmatstep.id + "_TimeStep", agmatstep.myphase->id + "_TimeScaling * " + agmatstep.id + "_InitialTimeStep"); }
					else { agmatstep.setCalculate(agmatstep.id + "_TimeStep", "-" + agmatstep.myphase->id + "_TimeScaling * " + agmatstep.id + "_InitialTimeStep"); }
					
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}

				//Case: [ MP- ] (i.e. the step just before the match-point); make this a half-step
				if (agmatstep.myspacecraft->isForward && agmatstep.isMatchPointStep) {
					GMATDebug << "MP-     ";
					//set the ForceModel and Propagator Parameters (false: is not a CloseApproach)
					agmatstep.setFMandProp(false);
					//scale the stepsize to be one-half
					agmatstep.scale_stepsize(0.5);
					//append the 'gmatstep' to the 'gmatphase' collector
					GMATMission.myjourneys[j].myphases[p].append_step(agmatstep);
					GMATDebug << "gs: " << agmatstep.gs << " id: " << agmatstep.id << "   delta-t: " << agmatstep.stepsize / 86400.0;
					GMATDebug << "   iepoch: " << agmatstep.iepoch / 86400.0 << "   fepoch: " << agmatstep.fepoch / 86400.0 << "   usePushBack: " << agmatstep.usePushBack;
					GMATDebug << "   inSOI0: " << agmatstep.inSOIatStart << "   inSOI1: " << agmatstep.inSOIatEnd << std::endl;
					//reset the 'gmatstep'
					agmatstep.reset();
				}
			}//end of Step level
		}//end of Phase level
	}//end of Journey level

}//end method


// method to generate additional variables, vary, calculate and constraint commands post 'gmat' class() creations
void gmatscripter::postpass_GMAT_phases() {

	GMATDebug << std::endl;
	GMATDebug << "postpass_GMAT_phases()" << std::endl;
	GMATDebug << std::endl;

	//declaration
	gmatphase* agmatphase;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//assignment
			agmatphase = &GMATMission.myjourneys[j].myphases[p];

			// -----------------------------
			//      POST-DATA COLLECTION
			// -----------------------------
			//calculate the phase's allocated time that can/cannot be varied
			agmatphase->get_time_allocation();



			// -----------------------------
			//       VARIABLE CREATION
			// -----------------------------
			//a variable that represents the amount of time that can be scaled during each journey
			agmatphase->setVariable(agmatphase->id + "_TimeScaling", 1.0);
			agmatphase->setVariable(agmatphase->id + "_EligibleTime", agmatphase->eligibletime / 86400.0);
			agmatphase->setVariable(agmatphase->id + "_MatchPoint_PositionError");
			agmatphase->setVariable(agmatphase->id + "_MatchPoint_VelocityError");
			agmatphase->setVariable(agmatphase->id + "_MatchPoint_MassError");
			agmatphase->setVariable(agmatphase->id + "_TOF");
			agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_Epoch");
			if (agmatphase->EndsWithFlyby) {
				//X
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_XLowerBound", 0.0);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_XWindow", agmatphase->spacecraft_backward.flyby_distance_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_XScaling", agmatphase->spacecraft_backward.initialconditions[0] / agmatphase->spacecraft_backward.flyby_distance_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_X");
				//Y
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_YLowerBound", 0.0);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_YWindow", agmatphase->spacecraft_backward.flyby_distance_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_YScaling", agmatphase->spacecraft_backward.initialconditions[1] / agmatphase->spacecraft_backward.flyby_distance_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_Y");
				//Z
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_ZLowerBound", 0.0);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_ZWindow", agmatphase->spacecraft_backward.flyby_distance_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_ZScaling", agmatphase->spacecraft_backward.initialconditions[2] / agmatphase->spacecraft_backward.flyby_distance_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_Z");
				//VX
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VXLowerBound", 0.0);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VXWindow", agmatphase->spacecraft_backward.flyby_velocity_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VXScaling", agmatphase->spacecraft_backward.initialconditions[3] / agmatphase->spacecraft_backward.flyby_velocity_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VX");
				//VY
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VYLowerBound", 0.0);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VYWindow", agmatphase->spacecraft_backward.flyby_velocity_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VYScaling", agmatphase->spacecraft_backward.initialconditions[4] / agmatphase->spacecraft_backward.flyby_velocity_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VY");
				//VZ
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VZLowerBound", 0.0);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VZWindow", agmatphase->spacecraft_backward.flyby_velocity_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VZScaling", agmatphase->spacecraft_backward.initialconditions[5] / agmatphase->spacecraft_backward.flyby_velocity_upperbound);
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VZ");
				//Position and Velocity Norms
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_PositionNorm");
				agmatphase->setVariable(agmatphase->spacecraft_backward.Name + "_VelocityNorm");
			}



			// -----------------------------
			//         VARY CREATION
			// -----------------------------
			agmatphase->setVary(agmatphase->id + "_TimeScaling", 1e-005, 0.0, 20.0, 0.1);
			if (agmatphase->EndsWithFlyby) {
				agmatphase->setVary(agmatphase->spacecraft_backward.Name + "_XScaling");
				agmatphase->setVary(agmatphase->spacecraft_backward.Name + "_YScaling");
				agmatphase->setVary(agmatphase->spacecraft_backward.Name + "_ZScaling");
				agmatphase->setVary(agmatphase->spacecraft_backward.Name + "_VXScaling");
				agmatphase->setVary(agmatphase->spacecraft_backward.Name + "_VYScaling");
				agmatphase->setVary(agmatphase->spacecraft_backward.Name + "_VZScaling");
			}


			// -----------------------------
			//      CALCULATE CREATION
			// -----------------------------
			//backward spacecraft epoch (1. create variable [see above], 2. set variable, 3. assign variable to the spacecraft)
			agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_Epoch",
									 agmatphase->spacecraft_forward.Name + ".Epoch." + agmatphase->spacecraft_forward.DateFormat + " + " + std::to_string(agmatphase->ineligibletime / 86400.0) +
		  						     " + ( " + agmatphase->id + "_TimeScaling * " + agmatphase->id + "_EligibleTime ) ");
			agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + ".Epoch." + agmatphase->spacecraft_backward.DateFormat, agmatphase->spacecraft_backward.Name + "_Epoch");
			//TOF
			agmatphase->setCalculate(agmatphase->id + "_TOF",
									 agmatphase->spacecraft_backward.Name + ".Epoch." + agmatphase->spacecraft_backward.DateFormat + " - " +
		 							 agmatphase->spacecraft_forward.Name + ".Epoch." + agmatphase->spacecraft_forward.DateFormat);
			//backward spacecraft X, Y, Z, VX, VY, VZ scaling and setting
			if (agmatphase->EndsWithFlyby) {
				//X
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_X", agmatphase->spacecraft_backward.Name + "_XScaling * " + agmatphase->spacecraft_backward.Name + "_XWindow" + " + " + agmatphase->spacecraft_backward.Name + "_XLowerBound");
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "." + agmatphase->spacecraft_backward.CoordinateSystem + ".X", agmatphase->spacecraft_backward.Name + "_X");
				//Y
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_Y", agmatphase->spacecraft_backward.Name + "_YScaling * " + agmatphase->spacecraft_backward.Name + "_YWindow" + " + " + agmatphase->spacecraft_backward.Name + "_YLowerBound");
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "." + agmatphase->spacecraft_backward.CoordinateSystem + ".Y", agmatphase->spacecraft_backward.Name + "_Y");
				//Z
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_Z", agmatphase->spacecraft_backward.Name + "_ZScaling * " + agmatphase->spacecraft_backward.Name + "_ZWindow" + " + " + agmatphase->spacecraft_backward.Name + "_ZLowerBound");
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "." + agmatphase->spacecraft_backward.CoordinateSystem + ".Z", agmatphase->spacecraft_backward.Name + "_Z");
				//VX
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_VX", agmatphase->spacecraft_backward.Name + "_VXScaling * " + agmatphase->spacecraft_backward.Name + "_VXWindow" + " + " + agmatphase->spacecraft_backward.Name + "_VXLowerBound");
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "." + agmatphase->spacecraft_backward.CoordinateSystem + ".VX", agmatphase->spacecraft_backward.Name + "_VX");
				//VY
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_VY", agmatphase->spacecraft_backward.Name + "_VYScaling * " + agmatphase->spacecraft_backward.Name + "_VYWindow" + " + " + agmatphase->spacecraft_backward.Name + "_VYLowerBound");
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "." + agmatphase->spacecraft_backward.CoordinateSystem + ".VY", agmatphase->spacecraft_backward.Name + "_VY");
				//VZ
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_VZ", agmatphase->spacecraft_backward.Name + "_VZScaling * " + agmatphase->spacecraft_backward.Name + "_VZWindow" + " + " + agmatphase->spacecraft_backward.Name + "_VZLowerBound");
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "." + agmatphase->spacecraft_backward.CoordinateSystem + ".VZ", agmatphase->spacecraft_backward.Name + "_VZ");
				//Position and Velocity Norm
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_PositionNorm", "sqrt( " + agmatphase->spacecraft_backward.Name + "_X ^ 2 + " + agmatphase->spacecraft_backward.Name + "_Y ^ 2 + " + agmatphase->spacecraft_backward.Name + "_Z ^ 2 )");
				agmatphase->setCalculate(agmatphase->spacecraft_backward.Name + "_VelocityNorm", "sqrt( " + agmatphase->spacecraft_backward.Name + "_VX ^ 2 + " + agmatphase->spacecraft_backward.Name + "_VY ^ 2 + " + agmatphase->spacecraft_backward.Name + "_VZ ^ 2 )");
			}
			//position error at matchpoint
			agmatphase->setCalculate(agmatphase->id + "_MatchPoint_PositionError",
									"sqrt(( " + agmatphase->spacecraft_forward.Name + ".SunJ2000Eq.X - " + agmatphase->spacecraft_backward.Name + ".SunJ2000Eq.X) ^ 2 + " +
									"( " + agmatphase->spacecraft_forward.Name + ".SunJ2000Eq.Y - " + agmatphase->spacecraft_backward.Name + ".SunJ2000Eq.Y) ^ 2 + " +
									"( " + agmatphase->spacecraft_forward.Name + ".SunJ2000Eq.Z - " + agmatphase->spacecraft_backward.Name + ".SunJ2000Eq.Z) ^ 2 )", false);
			//velocity error at matchpoint
			agmatphase->setCalculate(agmatphase->id + "_MatchPoint_VelocityError",
									"sqrt(( " + agmatphase->spacecraft_forward.Name + ".SunJ2000Eq.VX - " + agmatphase->spacecraft_backward.Name + ".SunJ2000Eq.VX) ^ 2 + " +
									"( " + agmatphase->spacecraft_forward.Name + ".SunJ2000Eq.VY - " + agmatphase->spacecraft_backward.Name + ".SunJ2000Eq.VY) ^ 2 + " +
									"( " + agmatphase->spacecraft_forward.Name + ".SunJ2000Eq.VZ - " + agmatphase->spacecraft_backward.Name + ".SunJ2000Eq.VZ) ^ 2 )", false);
			//mass error at matchpoint
			agmatphase->setCalculate(agmatphase->id + "_MatchPoint_MassError",
									agmatphase->spacecraft_forward.Name + "." + agmatphase->spacecraft_forward.Thruster.Tank.Name + ".FuelMass - " +
									agmatphase->spacecraft_backward.Name + "." + agmatphase->spacecraft_backward.Thruster.Tank.Name + ".FuelMass", false);



			// -----------------------------
			//      CONSTRAINT CREATION
			// -----------------------------
			//backward spacecraft Position and Velocity norm constraints during flyby
			if (agmatphase->EndsWithFlyby) {
				agmatphase->setConstraint(agmatphase->spacecraft_backward.Name + "_PositionNorm", ">=", agmatphase->spacecraft_backward.flyby_distance_lowerbound);
				agmatphase->setConstraint(agmatphase->spacecraft_backward.Name + "_PositionNorm", "<=", agmatphase->spacecraft_backward.flyby_distance_upperbound);
				agmatphase->setConstraint(agmatphase->spacecraft_backward.Name + "_VelocityNorm", ">=", agmatphase->spacecraft_backward.flyby_velocity_lowerbound);
				agmatphase->setConstraint(agmatphase->spacecraft_backward.Name + "_VelocityNorm", "<=", agmatphase->spacecraft_backward.flyby_velocity_upperbound);
			}



			/////////////// TEMPORARY ///////////////////////
			if (p == 1) {
				math::Matrix<double> flyby_state;
				GMATDebug << std::endl;
				GMATDebug << "///////////////////       Calculate Flyby Periapse State       ///////////////////" << std::endl;
				flyby_state = this->ptr_gmatmission->journeys[j].phases[p].calculate_flyby_periapse_state(this->ptr_gmatmission->journeys[j].phases[p].V_infinity_in,
					this->ptr_gmatmission->journeys[j].phases[p].V_infinity_out, this->ptr_gmatmission->journeys[j].phases[p].flyby_altitude, agmatphase->mybodies[0]);
				for (int i = 0; i < 6; ++i) {
					GMATDebug << flyby_state(i, 0) << std::endl;
				}
				GMATDebug << "///////////////////////////////////////////////////////////////////////////////////";
				GMATDebug << std::endl;
			}
			/////////////// TEMPORARY ///////////////////////

		}//end of phases for-statement
	}//end of journeys for-statement

}


// method to generate additional variables, vary, calculate and constraint commands post 'gmat' class() creations
void gmatscripter::postpass_GMAT_journeys() {

	//declarations
	stringstream tempstream;
	gmatjourney* agmatjourney;
	double waittimewindow;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		//assignment
		agmatjourney = &GMATMission.myjourneys[j];

		// -----------------------------
		//      POST-DATA COLLECTION
		// -----------------------------
		//waittimewindow
		waittimewindow = (this->ptr_gmatmission->options.journey_wait_time_bounds[j][1] - this->ptr_gmatmission->options.journey_wait_time_bounds[j][0]) / 86400.0;



		// -----------------------------
		//       VARIABLE CREATION
		// -----------------------------
		//if either mission or myjourney is time-constrained then we need to produce a TOF variable OR if the objective function is min. time.
		if (this->ptr_gmatmission->options.global_timebounded == 1 || this->ptr_gmatmission->options.journey_timebounded[agmatjourney->j] == 1 ||
			this->ptr_gmatmission->options.objective_type == 1) {
			agmatjourney->setVariable(agmatjourney->id + "_TOF");
		}
		//wait times available for each journey
		if (agmatjourney->j > 0) {
			if (waittimewindow > 1e-6) { 
				agmatjourney->setVariable(agmatjourney->id + "_WaitTimeScaling", (agmatjourney->myphases.front().iepoch - GMATMission.myjourneys[j - 1].myphases.back().fepoch - this->ptr_gmatmission->options.journey_wait_time_bounds[j][0]) / 86400.0 / waittimewindow);
				agmatjourney->setVariable(agmatjourney->id + "_WaitTimeWindow", waittimewindow);
				agmatjourney->setVariable(agmatjourney->id + "_WaitTimeLowerBound", this->ptr_gmatmission->options.journey_wait_time_bounds[j][0] / 86400.0);
			}
			else{ 
				agmatjourney->setVariable(agmatjourney->id + "_WaitTimeScaling", 0.0); 
				agmatjourney->setVariable(agmatjourney->id + "_WaitTimeWindow",  0.0);
				agmatjourney->setVariable(agmatjourney->id + "_WaitTimeLowerBound", 0.0);
			}
		}
		//dry mass
		if (j == 0) { 
			agmatjourney->setVariable(agmatjourney->id + "_DryMass", agmatjourney->myphases[0].spacecraft_forward.DryMass); 
            agmatjourney->setVariable(agmatjourney->id + "_DryMassLowerBound", this->ptr_gmatmission->options.final_mass_constraint);
		}
		else { agmatjourney->setVariable(agmatjourney->id + "_DryMass");
		}
		//mass increment scaling and window
		if (this->ptr_gmatmission->options.journey_variable_mass_increment[j] == 1 &&
			this->ptr_gmatmission->options.journey_starting_mass_increment[j] > 0.0) {

			double drymass_scaling = 0.0;
			//find the optimal drymass scaling for journey 'j'
			for (size_t iX = 0; iX < this->ptr_gmatmission->Xdescriptions.size(); ++iX) {
				//find index in Xdescriptions where the launch time bounds are located
				if (this->ptr_gmatmission->Xdescriptions[iX] == "j" + std::to_string(j) + "p0: journey initial mass scale factor") {
					drymass_scaling = this->ptr_gmatmission->Xopt[iX];
					agmatjourney->setVariable(agmatjourney->id + "_MassIncrementScaling", this->ptr_gmatmission->Xopt[iX]);
					break;
				}
			}
			//mass increment window
			agmatjourney->setVariable(agmatjourney->id + "_MassIncrementWindow", this->ptr_gmatmission->options.journey_starting_mass_increment[j]);
		}


		// -----------------------------
		//         VARY CREATION
		// -----------------------------
		if (agmatjourney->j > 0) { agmatjourney->setVary(agmatjourney->id + "_WaitTimeScaling"); }
		//dry (pick-up) mass increment variation
		if (this->ptr_gmatmission->options.journey_variable_mass_increment[j] == 1 &&
			this->ptr_gmatmission->options.journey_starting_mass_increment[j] > 0.0) {
			agmatjourney->setVary(agmatjourney->id + "_MassIncrementScaling", 1e-005, 0.0, 1.0, 0.1);
		}


		// -----------------------------
		//      CALCULATE CREATION
		// -----------------------------
		//if either mission or myjourney is time-constrained then we need to calculate the TOF variable based on phase TOF, OR if the objective function is min. time.
		if (this->ptr_gmatmission->options.global_timebounded == 1 || this->ptr_gmatmission->options.journey_timebounded[agmatjourney->j] == 1 || 
			this->ptr_gmatmission->options.objective_type == 1) {
			tempstream << agmatjourney->myphases[0].id << "_TOF";
			for (int p = 1; p < agmatjourney->myphases.size(); ++p) {
				tempstream << " + " << agmatjourney->myphases[p].id << "_TOF";
			}
			agmatjourney->setCalculate(agmatjourney->id + "_TOF", tempstream.str(), false);
			tempstream.str("");
		}
		//dry mass calculation after possible (pick-up) mass increment and mass margin variation
		if (this->ptr_gmatmission->options.journey_variable_mass_increment[j] == 1 &&
			this->ptr_gmatmission->options.journey_starting_mass_increment[j] > 0.0) {
			if (j == 0) {
				agmatjourney->setCalculate(agmatjourney->id + "_DryMass", agmatjourney->id + "_DryMassLowerBound" + " + " + agmatjourney->id + "_MassIncrementScaling" + " * " + agmatjourney->id + "_MassIncrementWindow" + " + " + "MassMargin");
			}
			else {
				agmatjourney->setCalculate(agmatjourney->id + "_DryMass", GMATMission.myjourneys[j - 1].id + "_DryMass" + " + " + agmatjourney->id + "_MassIncrementScaling" + " * " + agmatjourney->id + "_MassIncrementWindow");
			}
		}
		else if (j == 0) {
			if (this->ptr_gmatmission->options.journey_starting_mass_increment[j] != 0.0) {
				agmatjourney->setCalculate(agmatjourney->id + "_DryMass", agmatjourney->id + "_DryMassLowerBound" + " + " + "MassMargin" + std::to_string(this->ptr_gmatmission->options.journey_starting_mass_increment[j]));
			}
			else {
				agmatjourney->setCalculate(agmatjourney->id + "_DryMass", agmatjourney->id + "_DryMassLowerBound" + " + " + "MassMargin");
			}
		}
		else {
			if (this->ptr_gmatmission->options.journey_starting_mass_increment[j] != 0) {
				agmatjourney->setCalculate(agmatjourney->id + "_DryMass", GMATMission.myjourneys[j - 1].id + "_DryMass" + " + " + std::to_string(this->ptr_gmatmission->options.journey_starting_mass_increment[j]));
			}
			else {
				agmatjourney->setCalculate(agmatjourney->id + "_DryMass", GMATMission.myjourneys[j - 1].id + "_DryMass");
			}
		}


		// -----------------------------
		//      CONSTRAINT CREATION
		// -----------------------------
		//if myjourney is time-constrained then we need constrain the TOF variable
		if (this->ptr_gmatmission->options.journey_timebounded[agmatjourney->j] == 1) {
			agmatjourney->setConstraint(agmatjourney->id + "_TOF", ">=", this->ptr_gmatmission->options.journey_flight_time_bounds[agmatjourney->j][0] / 86400.0 );
			agmatjourney->setConstraint(agmatjourney->id + "_TOF", "<=", this->ptr_gmatmission->options.journey_flight_time_bounds[agmatjourney->j][1] / 86400.0);
		}

			
	}//end of journeys for-statement

}


// method to generate additional variables, vary, calculate and constraint commands post 'gmat' class() creations
void gmatscripter::postpass_GMAT_missions() {

	//declarations
	stringstream tempstream;
	double massmargin = 0.0;
	double massmarginwindow = 0.0;

	// -----------------------------
	//      POST-DATA COLLECTION
	// -----------------------------

	// -----------------------------
	//       VARIABLE CREATION
	// -----------------------------
	GMATMission.setVariable("ObjectiveFunction");
	//if a mission is time-constrained then we need to create a mission-level TOF to later constrain, OR if the objective function is min. time.
	if (this->ptr_gmatmission->options.global_timebounded == 1 || this->ptr_gmatmission->options.objective_type == 1) { GMATMission.setVariable("Mission_TOF"); }
	//mass margin
	if (GMATMission.myjourneys[0].myphases[0].spacecraft_forward.mass_margin > 0.0) {
		GMATMission.setVariable("MassMargin", GMATMission.myjourneys[0].myphases[0].spacecraft_forward.mass_margin);
		massmargin = GMATMission.myjourneys[0].myphases[0].spacecraft_forward.mass_margin;
	}
	else {
		GMATMission.setVariable("MassMargin", 0.0);
	}
	if (GMATMission.emtgmission->options.allow_initial_mass_to_vary) {
        massmarginwindow = GMATMission.emtgmission->options.maximum_mass - GMATMission.emtgmission->options.final_mass_constraint;
		GMATMission.setVariable("MassMarginWindow", massmarginwindow);
		GMATMission.setVariable("MassMarginScaling", massmargin / massmarginwindow);
		GMATMission.setVariable("MaximumMass");
	}

	

	// -----------------------------
	//         VARY CREATION
	// -----------------------------
	//mass margin
	if (GMATMission.emtgmission->options.allow_initial_mass_to_vary) {
		GMATMission.setVary("MassMarginScaling");
	}

	// -----------------------------
	//      CALCULATE CREATION
	// -----------------------------
	//if a mission is time-constrained then we need to calculate the mission-level TOF OR if the objective function is min. time.
	if (this->ptr_gmatmission->options.global_timebounded == 1 || this->ptr_gmatmission->options.objective_type == 1) {
		tempstream << GMATMission.myjourneys[0].id << "_TOF";
		for (int j = 1; j < GMATMission.myjourneys.size(); ++j) {
			tempstream << " + " << GMATMission.myjourneys[j].id << "_TOF";
		}
		GMATMission.setCalculate("Mission_TOF", tempstream.str(), false);
		tempstream.str("");
	}
	//objective function types
	if (this->ptr_gmatmission->options.objective_type == 1) {
		GMATMission.setCalculate("ObjectiveFunction", "Mission_TOF", false);
	}
	else if (this->ptr_gmatmission->options.objective_type == 2) {
		GMATMission.setCalculate("ObjectiveFunction", "-" + GMATMission.myjourneys.back().myphases.back().spacecraft_backward.Name + "_FuelMass" + 
													  " - " + GMATMission.myjourneys.back().id + "_DryMass", false);
	}
	//mass margin
	if (GMATMission.emtgmission->options.allow_initial_mass_to_vary) {
		GMATMission.setCalculate("MassMargin", "MassMarginScaling * MassMarginWindow");
	}
	//maximum mass
	if (GMATMission.emtgmission->options.allow_initial_mass_to_vary) {
		GMATMission.setCalculate("MaximumMass", GMATMission.myjourneys[0].id + "_DryMassLowerBound" + " + " + "MassMargin" + " + " +
												GMATMission.myjourneys[0].myphases[0].spacecraft_forward.Name + "_FuelMass", false);
	}

	// -----------------------------
	//      CONSTRAINT CREATION
	// -----------------------------
	//if mymission is time-constrained then we need constrain the TOF variable
	if (this->ptr_gmatmission->options.global_timebounded == 1) {
		GMATMission.setConstraint("Mission_TOF", ">=", this->ptr_gmatmission->options.total_flight_time_bounds[0] / 86400.0);
		GMATMission.setConstraint("Mission_TOF", "<=", this->ptr_gmatmission->options.total_flight_time_bounds[1] / 86400.0);
	}
	//need to ensure that the sum of the dry mass, propellant, and margin at the start do not violate the maximum
	//this should only be necessary if we allow initial mass to vary because then mass margin and initial propellant may be more than the max.
	if (GMATMission.emtgmission->options.allow_initial_mass_to_vary) {
		GMATMission.setConstraint("MaximumMass", "<=", this->ptr_gmatmission->options.maximum_mass);
	}


}


// method to create the script preamble
void gmatscripter::write_GMAT_preamble() {

	//get the current timestamp from boost and assign it to now
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::stringstream timestream;
	timestream << static_cast<int>(now.date().month()) << "/" << now.date().day() << "/" << now.date().year() << " " << now.time_of_day().hours() << ":" << now.time_of_day().minutes() << ":" << now.time_of_day().seconds();

	//preamble
	GMATfile << "%--------------------------------------------------------------------------------" << std::endl;
	GMATfile << "%GMAT script created by EMTGv8" << std::endl;
	GMATfile << "%EMTG options file: " << this->ptr_gmatmission->options.working_directory + "/" + 
										  this->ptr_gmatmission->options.mission_name + ".emtgopt" << std::endl;
	GMATfile << "%EMTG output file: " <<  this->ptr_gmatmission->options.outputfile << std::endl;
	GMATfile << "%EMTG output written on: " << timestream.str() << std::endl;
	GMATfile << "%Author(s): Ryne Beeson     (2014_07_01)" << std::endl;
	GMATfile << "%           Max Schadegg    (2013_06_03)" << std::endl;
	GMATfile << "%           Jacob Englander (2013_08_09)" << std::endl;
	GMATfile << "%--------------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;
	GMATfile << std::endl;

}


// method to create spacecraft information
void gmatscripter::write_GMAT_spacecraft() {

	GMATDebug << "%-------------------------------------------------------------------------" << std::endl;
	GMATDebug << "%---------- Spacecraft" << std::endl;
	GMATDebug << "%-------------------------------------------------------------------------" << std::endl;
	GMATDebug << std::endl;

	//spacecraft header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Spacecraft" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//write out the forward spacecraft
			this->create_GMAT_spacecraft(GMATMission.myjourneys[j].myphases[p].spacecraft_forward);
			//write out the backward spacecraft
			this->create_GMAT_spacecraft(GMATMission.myjourneys[j].myphases[p].spacecraft_backward);

			GMATDebug << GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name << std::endl;
			GMATDebug << GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name << std::endl;

		}
	}
	GMATfile << std::endl;

	GMATDebug << "%-------------------------------------------------------------------------" << std::endl;

}//end of write_GMAT_spacecraft() method


// method to create hardware information
void gmatscripter::write_GMAT_hardware() {

	//hardware header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Hardware components" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

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
	GMATfile << std::endl;

}//end of write_GMAT_hardware() method


// method to create models for bodies that are non-standard (i.e. not a planet nor the moon)
void gmatscripter::write_GMAT_nonstandardbody(){

	//force model header 
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- NonStandard Body Models" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

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
			GMATfile << "%Must create model for body visited" << std::endl;
			GMATfile << "Create Planet " << GMATMission.missionbodies_unique[body_index].name << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".NAIFId = " << GMATMission.missionbodies_unique[body_index].spice_ID << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".EquatorialRadius = " << GMATMission.missionbodies_unique[body_index].radius << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".Mu = " << GMATMission.missionbodies_unique[body_index].mu << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".PosVelSource = 'SPICE'" << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".CentralBody = '" << GMATMission.missionbodies_unique[body_index].central_body_name << "'" << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".RotationDataSource = 'IAUSimplified'" << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".OrientationEpoch = 21545" << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisRAConstant = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.alpha0 << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisRARate = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.alphadot << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisDECConstant = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.delta0 << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".SpinAxisDECRate = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.deltadot << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".RotationConstant = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.W << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].name << ".RotationRate = " << GMATMission.missionbodies_unique[body_index].J2000_body_equatorial_frame.Wdot << std::endl;
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
			GMATfile << "'" << filestring << "'};" << std::endl;
			GMATfile << std::endl;
		}

		//must create any central bodies that are not already defined in GMAT
		if ((GMATMission.missionbodies_unique[body_index].central_body_name != "Sun") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Mercury") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Venus") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Earth") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Mars") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Jupiter") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Saturn") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Uranus") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Neptune") && (GMATMission.missionbodies_unique[body_index].central_body_name != "Pluto"))
		{
			GMATfile << "%Must create model for central body" << std::endl;
			GMATfile << "Create Planet " << GMATMission.missionbodies_unique[body_index].central_body_name << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".NAIFId = " << GMATMission.missionbodies_unique[body_index].central_body_spice_ID << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".EquatorialRadius = " << GMATMission.missionbodies_unique[body_index].central_body_radius << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".Mu = " << GMATMission.missionbodies_unique[body_index].universe_mu << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".PosVelSource = 'SPICE'" << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".CentralBody = 'Sun'" << std::endl; //assume Sun for now
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".RotationDataSource = 'IAUSimplified'" << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".OrientationEpoch = 21545" << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisRAConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alpha0 << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisRARate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.alphadot << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisDECConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.delta0 << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".SpinAxisDECRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.deltadot << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".RotationConstant = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.W << std::endl;
			GMATfile << GMATMission.missionbodies_unique[body_index].central_body_name << ".RotationRate = " << this->ptr_gmatmission->TheUniverse[0].LocalFrame.Wdot << std::endl;
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
			GMATfile << "'" << filestring << "'};" << std::endl;
			GMATfile << std::endl;
		}
	}
}//end of write_GMAT_nonstandardbody() method


// method to create forcemodel information
void gmatscripter::write_GMAT_forcemodels(){

	//declaration
	vector <string> tempstrings;
	bool hasNotBeenCopied;

	//force model header 
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Force Models" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

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
	GMATfile << std::endl;

}//end of write_GMAT_forcemodels() method


// method to create propagator information
void gmatscripter::write_GMAT_propagators(){

	//declaration
	vector <string> tempstrings;
	bool hasNotBeenCopied;

	//propagator header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Propagators" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

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
	GMATfile << std::endl;

}//end of write_GMAT_propagators() method


// method to create burn information
void gmatscripter::write_GMAT_burns(){

	//finite burn header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Finite Burns" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//write out the forward spacecraft
			this->create_GMAT_burn(GMATMission.myjourneys[j].myphases[p].spacecraft_forward.fBurn);
			GMATfile << std::endl;
			//write out the backward spacecraft
			this->create_GMAT_burn(GMATMission.myjourneys[j].myphases[p].spacecraft_backward.fBurn);
			GMATfile << std::endl;
		}
	}

	//impulsive burn header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Impulsive Burns" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//write out the forward spacecraft
			this->create_GMAT_burn(GMATMission.myjourneys[j].myphases[p].spacecraft_forward.iBurn);
			//write out the backward spacecraft
			this->create_GMAT_burn(GMATMission.myjourneys[j].myphases[p].spacecraft_backward.iBurn);
		}
	}
	GMATfile << std::endl;

}//end of write_GMAT_burns() method


// method to create coordinate systems information
void gmatscripter::write_GMAT_coordinatesystems(){

	//coordinate systems header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Coordinate systems" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	//always in J2000 Eq. Syst.
	//TODO:: as of now it is assumed that the first central body is the principle central body
	GMATfile << "%Create coordinate systems for plotting/viewing" << std::endl;
	this->create_GMAT_coordinatesystem(GMATMission.missionbodies_unique[0].central_body_name);

	for (int b = 0; b < GMATMission.missionbodies_unique.size(); ++b) {
		this->create_GMAT_coordinatesystem(GMATMission.missionbodies_unique[b].name);
	}
	//add some vertical whitespace
	GMATfile << std::endl;

}//end of write_GMAT_coordinatesystems() method


// method to create solver information
void gmatscripter::write_GMAT_solvers(){

	//solvers header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Solvers" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

    switch (this->ptr_gmatmission->options.GMAT_optimizer)
    {
        //fmincon or vf13ad
    case 0:
        //default optimizer is VF13
        GMATfile << "Create VF13ad NLPObject;" << std::endl;
        GMATfile << "NLPObject.ShowProgress = true;" << std::endl;
        GMATfile << "NLPObject.ReportStyle = Normal;" << std::endl;
        GMATfile << "NLPObject.ReportFile = 'VF13adNLPObject.data';" << std::endl;
        GMATfile << "NLPObject.MaximumIterations = 100;" << std::endl;
        GMATfile << "NLPObject.Tolerance = 1e-004;" << std::endl;
        GMATfile << "NLPObject.UseCentralDifferences = false;" << std::endl;
        GMATfile << "NLPObject.FeasibilityTolerance = 0.1;" << std::endl;
        GMATfile << std::endl;
        break;

    case 1:
        //SNOPT
        GMATfile << "Create SNOPT NLPObject;" << std::endl;
        GMATfile << "NLPObject.ShowProgress = true;" << std::endl;
        GMATfile << "NLPObject.ReportStyle = Normal;" << std::endl;
        GMATfile << "NLPObject.ReportFile = 'SNOPTObject.data';" << std::endl;
        GMATfile << "NLPObject.MaximumIterations = 100;" << std::endl;
        GMATfile << "NLPObject.Tolerance = 1e-004;" << std::endl;
        //GMATfile << "NLPObject.UseCentralDifferences = false;" << std::endl;
        //GMATfile << "NLPObject.FeasibilityTolerance = 0.1;" << std::endl;
        GMATfile << std::endl;
        break;

    case 2:
        //fmincon
        GMATfile << "Create FminconOptimizer NLPObject;" << std::endl;
        GMATfile << "NLPObject.DiffMaxChange = '0.1000';" << std::endl;
        GMATfile << "NLPObject.DiffMinChange = '1.0000e-08';" << std::endl;
        GMATfile << "NLPObject.MaxFunEvals = '1000';" << std::endl;
        GMATfile << "NLPObject.TolX = '1.0000e-06';" << std::endl;
        GMATfile << "NLPObject.TolFun = '1.0000e-06';" << std::endl;
        GMATfile << "NLPObject.TolCon = '1.0000e-06';" << std::endl;
        break;
    }
	GMATfile << std::endl;
	GMATfile << std::endl;

}//end of write_GMAT_solvers() method


// method to create subscriber information
void gmatscripter::write_GMAT_subscribers(){

	//subscribers header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Subscribers" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	//includes orbit views, plots, etc...
	//add central universe body view
	GMATfile << "%Create subscriber for central body view" << std::endl;
	GMATfile << "Create OrbitView " << GMATMission.missionbodies_unique[0].central_body_name << "View" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ShowPlot =		true" << std::endl;
    if (this->ptr_gmatmission->options.GMAT_plot_while_optimize)
	    GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.SolverIterations =	 All" << std::endl;
    else
        GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.SolverIterations =	 None" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.RelativeZOrder =	501" << std::endl;

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
	GMATfile << ", " << GMATMission.missionbodies_unique[0].central_body_name << "}" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.CoordinateSystem =		" << GMATMission.missionbodies_unique[0].central_body_name << "J2000Eq" << std::endl;
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
	GMATfile << "true]" << std::endl;
	//other parameters
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.DataCollectFrequency   = 1" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.UpdatePlotFrequency    = 50" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.NumPointsToRedraw      = 300" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewScaleFactor        = 35" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewPointReference	  = " << GMATMission.missionbodies_unique[0].central_body_name << ";" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewDirection		  = " << GMATMission.missionbodies_unique[0].central_body_name << ";" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewPointVector		  = [ 0 0 30000000 ];" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.ViewUpAxis             = X" << std::endl;
	GMATfile << GMATMission.missionbodies_unique[0].central_body_name << "View.UseInitialView         = On" << std::endl;
	GMATfile << std::endl;

	//create orbit views for all bodies visited
	GMATfile << "%Create subscribers for other body views" << std::endl;
	for (int body_index = 0; body_index < GMATMission.missionbodies_unique.size(); ++body_index)
	{
		GMATfile << "Create OrbitView " << GMATMission.missionbodies_unique[body_index].name << "View" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ShowPlot               = true" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.SolverIterations       = All" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.RelativeZOrder         = 501" << std::endl;

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
		GMATfile << "}" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.CoordinateSystem       = " << GMATMission.missionbodies_unique[body_index].name << "J2000Eq" << std::endl;
		//bool flag parameters for drawing objects
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.DrawObject             = [";
		for (int index_plot = 0; index_plot < (GMATMission.missionbodies_unique.size()); ++index_plot)
		{
			GMATfile << "true ";
		}
		for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
			for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
				GMATfile << "true true ";
			}
		}
		GMATfile << "]" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.DataCollectFrequency   = 1" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.UpdatePlotFrequency    = 50" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.NumPointsToRedraw      = 300" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewScaleFactor        = 35" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewUpAxis             = X" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.UseInitialView         = On" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewPointReference	  = " << GMATMission.missionbodies_unique[body_index].name << ";" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewDirection		  = " << GMATMission.missionbodies_unique[body_index].name << ";" << std::endl;
		GMATfile << GMATMission.missionbodies_unique[body_index].name << "View.ViewPointVector		  = [ 0 0 20000 ];" << std::endl;
		GMATfile << std::endl;
	}
	GMATfile << std::endl;
	GMATfile << std::endl;

	//declarations
	string report_name;
	stringstream report_name_stream;

	//create reports for debugging purposes
	GMATfile << "%Create reports for debugging purposes" << std::endl;
	GMATfile << "%Create a report for the central body and each unique body of the mission" << std::endl;

	report_name_stream << "Report_Spacecraft_" << GMATMission.missionbodies[0].central_body_name << "_States";
	report_name = report_name_stream.str();
	GMATfile << "Create ReportFile " << report_name << std::endl;
	GMATfile << report_name << ".SolverIterations = Current;" << std::endl;
	GMATfile << report_name << ".UpperLeft = [ 0 0 ];" << std::endl;
	GMATfile << report_name << ".Size = [ 0 0 ];" << std::endl;
	GMATfile << report_name << ".RelativeZOrder = 0;" << std::endl;
	GMATfile << report_name << ".Maximized = false;" << std::endl;
	GMATfile << report_name << ".Filename = 'Report_" << GMATMission.missionbodies[0].central_body_name << "Centered_States.txt';" << std::endl;
	GMATfile << report_name << ".Precision = 16;" << std::endl;
	GMATfile << report_name << ".WriteHeaders = true;" << std::endl;
	GMATfile << report_name << ".LeftJustify = On;" << std::endl;
	GMATfile << report_name << ".ZeroFill = Off;" << std::endl;
	GMATfile << report_name << ".ColumnWidth = 20;" << std::endl;
	GMATfile << report_name << ".WriteReport = true;" << std::endl;
	GMATfile << std::endl;

	for (int body_index = 0; body_index < GMATMission.missionbodies_unique.size(); body_index++){
		// new report name
		report_name.erase (report_name.begin(), report_name.end());
		//report_name_stream << "Report_Spacecraft_" << GMATMission.missionbodies_unique[body_index].name << "_States";
		//report_name = report_name_stream.str();
		report_name = "Report_Spacecraft_" + GMATMission.missionbodies_unique[body_index].name + "_States";
		// GMAT script printing
		GMATfile << "Create ReportFile " << report_name << std::endl;
		GMATfile << report_name << ".SolverIterations = Current;" << std::endl;
		GMATfile << report_name << ".UpperLeft = [ 0 0 ];" << std::endl;
		GMATfile << report_name << ".Size = [ 0 0 ];" << std::endl;
		GMATfile << report_name << ".RelativeZOrder = 0;" << std::endl;
		GMATfile << report_name << ".Maximized = false;" << std::endl;
		GMATfile << report_name << ".Filename = 'Report_" << GMATMission.missionbodies_unique[body_index].name << "Centered_States.txt';" << std::endl;
		GMATfile << report_name << ".Precision = 16;" << std::endl;
		GMATfile << report_name << ".WriteHeaders = true;" << std::endl;
		GMATfile << report_name << ".LeftJustify = On;" << std::endl;
		GMATfile << report_name << ".ZeroFill = Off;" << std::endl;
		GMATfile << report_name << ".ColumnWidth = 20;" << std::endl;
		GMATfile << report_name << ".WriteReport = true;" << std::endl;
		GMATfile << std::endl;
	}

	GMATfile << "Create ReportFile Report_SpacecraftControl;" << std::endl;
	GMATfile << "Report_SpacecraftControl.SolverIterations = Current;" << std::endl;
	GMATfile << "Report_SpacecraftControl.UpperLeft = [ 0 0 ];" << std::endl;
	GMATfile << "Report_SpacecraftControl.Size = [ 0 0 ];" << std::endl;
	GMATfile << "Report_SpacecraftControl.RelativeZOrder = 0;" << std::endl;
	GMATfile << "Report_SpacecraftControl.Maximized = false;" << std::endl;
	GMATfile << "Report_SpacecraftControl.Filename = 'Report_SpaceCraftControlHistory.txt';" << std::endl;
	GMATfile << "Report_SpacecraftControl.Precision = 16;" << std::endl;
	GMATfile << "Report_SpacecraftControl.WriteHeaders = true;" << std::endl;
	GMATfile << "Report_SpacecraftControl.LeftJustify = On;" << std::endl;
	GMATfile << "Report_SpacecraftControl.ZeroFill = Off;" << std::endl;
	GMATfile << "Report_SpacecraftControl.ColumnWidth = 20;" << std::endl;
	GMATfile << "Report_SpacecraftControl.WriteReport = true;" << std::endl;
	GMATfile << std::endl;
	GMATfile << std::endl;

	//XY Plots for visually informing user of the optimization progress
	GMATfile << "%Creating some plots for to inform the user of the optimization progress" << std::endl;
	GMATfile << "Create XYPlot PositionErrorPlot" << std::endl;
	GMATfile << "PositionErrorPlot.SolverIterations = All" << std::endl;
	GMATfile << "PositionErrorPlot.UpperLeft = [0 0]" << std::endl;
	GMATfile << "PositionErrorPlot.Size = [0 0]" << std::endl;
	GMATfile << "PositionErrorPlot.RelativeZOrder = 0" << std::endl;
	GMATfile << "PositionErrorPlot.XVariable = Iterations" << std::endl;
	GMATfile << "PositionErrorPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_PositionError"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_PositionError"; }
		}
	}
	GMATfile << " }" << std::endl;
	GMATfile << "PositionErrorPlot.ShowGrid = true" << std::endl;
	GMATfile << "PositionErrorPlot.ShowPlot = true" << std::endl;
	GMATfile << std::endl;

	GMATfile << "Create XYPlot VelocityErrorPlot" << std::endl;
	GMATfile << "VelocityErrorPlot.SolverIterations = All" << std::endl;
	GMATfile << "VelocityErrorPlot.UpperLeft = [0 0]" << std::endl;
	GMATfile << "VelocityErrorPlot.Size = [0 0]" << std::endl;
	GMATfile << "VelocityErrorPlot.RelativeZOrder = 0" << std::endl;
	GMATfile << "VelocityErrorPlot.XVariable = Iterations" << std::endl;
	GMATfile << "VelocityErrorPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_VelocityError"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_VelocityError"; }
		}
	}
	GMATfile << " }" << std::endl;
	GMATfile << "VelocityErrorPlot.ShowGrid = true" << std::endl;
	GMATfile << "VelocityErrorPlot.ShowPlot = true" << std::endl;
	GMATfile << std::endl;

	GMATfile << "Create XYPlot MassErrorPlot" << std::endl;
	GMATfile << "MassErrorPlot.SolverIterations = All" << std::endl;
	GMATfile << "MassErrorPlot.UpperLeft = [0 0]" << std::endl;
	GMATfile << "MassErrorPlot.Size = [0 0]" << std::endl;
	GMATfile << "MassErrorPlot.RelativeZOrder = 0" << std::endl;
	GMATfile << "MassErrorPlot.XVariable = Iterations" << std::endl;
	GMATfile << "MassErrorPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_MassError"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_MatchPoint_MassError"; }
		}
	}
	GMATfile << " }" << std::endl;
	GMATfile << "MassErrorPlot.ShowGrid = true" << std::endl;
	GMATfile << "MassErrorPlot.ShowPlot = true" << std::endl;
	GMATfile << std::endl;

	GMATfile << "Create XYPlot TOFPlot" << std::endl;
	GMATfile << "TOFPlot.SolverIterations = All" << std::endl;
	GMATfile << "TOFPlot.UpperLeft = [0 0]" << std::endl;
	GMATfile << "TOFPlot.Size = [0 0]" << std::endl;
	GMATfile << "TOFPlot.RelativeZOrder = 0" << std::endl;
	GMATfile << "TOFPlot.XVariable = Iterations" << std::endl;
	GMATfile << "TOFPlot.YVariables = { ";
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			if (j == 0 && p == 0) { GMATfile << GMATMission.myjourneys[j].myphases[p].id << "_TOF"; }
			else { GMATfile << " ," << GMATMission.myjourneys[j].myphases[p].id << "_TOF"; }
		}
	}
	GMATfile << " }" << std::endl;
	GMATfile << "TOFPlot.ShowGrid = true" << std::endl;
	GMATfile << "TOFPlot.ShowPlot = true" << std::endl;
	GMATfile << std::endl;

	GMATfile << "Create XYPlot ObjectiveFunctionPlot" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.SolverIterations = All" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.UpperLeft = [0 0]" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.Size = [0 0]" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.RelativeZOrder = 0" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.XVariable = Iterations" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.YVariables = { ObjectiveFunction }" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.ShowGrid = true" << std::endl;
	GMATfile << "ObjectiveFunctionPlot.ShowPlot = true" << std::endl;
	GMATfile << std::endl;

	//An iteration counter for the Optimization Sequence (and necessary for the XY Plots above)
	GMATMission.setVariable("Iterations", 0);
	GMATMission.setCalculate("Iterations", "Iterations + 1", false);
	GMATfile << std::endl;

}//end of write_GMAT_subscribers() method


// method to create array and variable information
void gmatscripter::write_GMAT_variables(){

	//arrays and variables header
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Arrays, Variables, Strings" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	//write out mission level variables
	GMATfile << "% --- Mission Level: Arrays, Variables, and Strings" << std::endl;
	GMATMission.printVariable(GMATfile);
	GMATfile << std::endl;

	//write out journey level variables
	GMATfile << "% --- Journey Level: Arrays, Variables, and Strings" << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) { GMATMission.myjourneys[j].printVariable(GMATfile); }
	GMATfile << std::endl;

	//write out the phase level variables
	GMATfile << "% --- Phase Level: Arrays, Variables, and Strings" << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			GMATMission.myjourneys[j].myphases[p].printVariable(GMATfile);
		}
		GMATfile << std::endl;
	}
	GMATfile << std::endl;

	//write out the step level variables
	GMATfile << "% --- Step Level: Arrays, Variables, and Strings" << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printVariable(GMATfile);
			}	
			GMATfile << std::endl;
		}
	}
	GMATfile << std::endl;

	//write out report variables
	//create temporary strings to identify time steps during phases
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				GMATfile << "Create String tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tminus" << std::endl;
				GMATfile << "tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tminus = '" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tminus: '" << std::endl;
				GMATfile << "Create String tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tplus" << std::endl;
				GMATfile << "tempString_" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tplus = '" << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << "_tplus: '" << std::endl;
			}
		}
	}
	GMATfile << std::endl;
	GMATfile << std::endl;

}//end of write_GMAT_variables() method


// method to create the beginmissionsequence statement
void gmatscripter::write_GMAT_beginmissionsequence(){

	GMATfile << "BeginMissionSequence" << std::endl;
	GMATfile << std::endl;
	GMATfile << std::endl;

}


// method to create the initial conditions for the mission sequence
void gmatscripter::write_GMAT_initialconditions(){

	//declarations
	double boundary_state[6];
	math::Matrix<double> Vinf_in(3, 1);
	math::Matrix<double> Vinf_out(3, 1);
	math::Matrix<double> periapse_state_vector(6, 1);
	math::Matrix<double> periapse_position_vector(3, 1);
	math::Matrix<double> periapse_velocity_vector(3, 1);
	int body_index = 0;
	int name_index = 0;

	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Initial State Guesses" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	GMATfile << "BeginScript 'Initial Guess Values' " << std::endl;
	GMATfile << std::endl;

	//write out the phase level variables
	GMATfile << "% --- Phase Level: Arrays, Variables, and Strings" << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//forward spacecraft
			this->create_GMAT_initialconditions(GMATMission.myjourneys[j].myphases[p].spacecraft_forward);
			GMATfile << std::endl;
			//backward spacecraft
			this->create_GMAT_initialconditions(GMATMission.myjourneys[j].myphases[p].spacecraft_backward);
			GMATfile << std::endl;
		}
	}

	GMATfile << std::endl;
	GMATfile << "% --- Additional / Updated Initial Guess" << std::endl;
	GMATfile << std::endl;

	GMATfile << std::endl;
	GMATfile << "EndScript" << std::endl;
	GMATfile << std::endl;

}//end of write_GMAT_initialguess() method


// method to create the initial conditions for the mission sequence
void gmatscripter::write_GMAT_preoptimization_calculations(){

	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Pre-Optimization Calculations" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	GMATfile << "BeginScript 'Pre Optimization Calculations' " << std::endl;
	GMATfile << std::endl;

	//write out the phase level variables
	GMATfile << "% --- Mission Level: " << std::endl;
	GMATMission.printPreOptCalculate(GMATfile);
	GMATfile << std::endl;

	//write out the phase level variables
	GMATfile << "% --- Journey Level: " << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		GMATMission.myjourneys[j].printPreOptCalculate(GMATfile);
	}
	GMATfile << std::endl;

	//write out the phase level variables
	GMATfile << "% --- Phase Level: " << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			GMATMission.myjourneys[j].myphases[p].printPreOptCalculate(GMATfile);
		}
	}
	GMATfile << std::endl;

	//write out the phase level variables
	GMATfile << "% --- Step Level: " << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printPreOptCalculate(GMATfile);
			}
		}
	}

	GMATfile << std::endl;
	GMATfile << "EndScript" << std::endl;
	GMATfile << std::endl;

}//end of write_GMAT_initialguess() method


// method to create the optimization sequence for the mission
void gmatscripter::write_GMAT_optimization() {

	GMATfile << "	" << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "	" << "%---------- Optimization Sequence" << std::endl;
	GMATfile << "	" << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	GMATfile << "Optimize 'OptimizeSequence' NLPObject {SolveMode = Solve, ExitMode = DiscardAndContinue}" << std::endl;
	GMATfile << std::endl;

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
	GMATfile << std::endl;

	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		//Journey ID
		GMATfile << "   " << "% Journey ID: " << GMATMission.myjourneys[j].id << std::endl;
		//Jouney level vary and calculate commands
		GMATMission.myjourneys[j].printVary(GMATfile);
		GMATMission.myjourneys[j].printCalculate(GMATfile);
		GMATfile << std::endl;
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			//Phase ID
			GMATfile << "   " << "% Phase ID: " << GMATMission.myjourneys[j].myphases[p].id << std::endl;
			//Phase level vary and calculate commands
			GMATMission.myjourneys[j].myphases[p].printVary(GMATfile);
			GMATMission.myjourneys[j].myphases[p].printCalculate(GMATfile);
			//  reprint phase boundary conditions if necessary (could have used the 'Calculate' command as well)
			this->rewrite_GMAT_initialconditions(GMATMission.myjourneys[j].myphases[p]);
			GMATfile << std::endl;
			for (int gs = 0; gs < GMATMission.myjourneys[j].myphases[p].mysteps.size(); ++gs) {
				//Step ID
				GMATfile << "   " << "% Step ID: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << std::endl;
				//Maneuver command if Impulsive Burn is used
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printManeuver(GMATfile);
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
				GMATDebug << "spacecraft name: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].myspacecraft->Name << " isForward: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].myspacecraft->isForward << std::endl;
				GMATDebug << "id: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].id << " stepsize: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].stepsize << "  zeroPropagate: " << GMATMission.myjourneys[j].myphases[p].mysteps[gs].zeroPropagate;
				this->write_GMAT_report(GMATMission.myjourneys[j].myphases[p].mysteps[gs], true, true);
				this->aux_GMAT_propagate(GMATMission.myjourneys[j].myphases[p].mysteps[gs], false);
				this->write_GMAT_report(GMATMission.myjourneys[j].myphases[p].mysteps[gs], false, true);
				//'EndBurn' command
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printEndBurn(GMATfile);
				//Step level calculate commands
				GMATMission.myjourneys[j].myphases[p].mysteps[gs].printCalculate(GMATfile, 0);
				GMATfile << std::endl;
			}//end of Step level
			//Phase level constraint commands
			GMATMission.myjourneys[j].myphases[p].printCalculate(GMATfile, 0);
			GMATMission.myjourneys[j].myphases[p].printConstraint(GMATfile);
			GMATfile << std::endl;
		}//end of Phase level
		//Journey level calculate and constraint commands
		GMATMission.myjourneys[j].printCalculate(GMATfile, 0);
		GMATMission.myjourneys[j].printConstraint(GMATfile);
		GMATfile << std::endl;
	}//end of Journey level

	//Mission level calculate and constraint commands
	GMATMission.printCalculate(GMATfile, 0);
	GMATMission.printConstraint(GMATfile);
	GMATfile << std::endl;









	//optimize the user-defined objective function
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << "%---------- Objective Function" << std::endl;
	GMATfile << "%-------------------------------------------------------------------------" << std::endl;
	GMATfile << std::endl;

	//currently vf13 does not like to obey set bounds if an objective function is given
	GMATfile << "	%Minimize objective function" << std::endl;
	GMATfile << "	Minimize NLPObject(ObjectiveFunction)" << std::endl;
	GMATfile << std::endl;
	GMATfile << "EndOptimize" << std::endl;

	//TEMPORARY DEBUG CODE
	GMATDebug << std::endl;
	GMATDebug << " *** spacecraft vector<struct> print out *** " << std::endl;
	for (int j = 0; j < GMATMission.myjourneys.size(); ++j) {
		for (int p = 0; p < GMATMission.myjourneys[j].myphases.size(); ++p) {
			GMATDebug << GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Name << " | " << "Epoch: " << GMATMission.myjourneys[j].myphases[p].spacecraft_forward.Epoch << std::endl;
			GMATDebug << GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Name << " | " << "Epoch: " << GMATMission.myjourneys[j].myphases[p].spacecraft_backward.Epoch << std::endl;
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
		GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->Thruster.Name << ".C1" << std::endl;
	}
	// write the spacecraft states in the central body frame
	GMATfile << "	Report 'Report_SpacecraftState' Report_Spacecraft_" << agmatstep.myphase->mybodies[body_index].name << "_States " << tempString;
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".X ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".Y ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".Z ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".VX ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".VY ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->CoordinateSystem << ".VZ ";
	GMATfile << agmatstep.myspacecraft->Name << "." << agmatstep.myspacecraft->Thruster.Tank.Name << ".FuelMass ";
	//if 'true' then write out the 'gmatsteps' phase variables
	if (true) {
		for (int i = 0; i < agmatstep.myphase->variables.size(); ++i) {
			GMATfile << agmatstep.myphase->variables[i][0] << " ";
		}
	}
	GMATfile << std::endl;

}


// method to write a GMAT Propagate Command
void gmatscripter::aux_GMAT_propagate(class gmatstep& agmatstep, bool useZeroPropagate){
	if (useZeroPropagate && agmatstep.zeroPropagate) { GMATfile << "	Propagate 'Propagate " << agmatstep.myspacecraft->Name << "' " << agmatstep.propagator.Name << "( " << agmatstep.myspacecraft->Name << " ) { " << agmatstep.myspacecraft->Name << ".ElapsedSecs = " << 0.0 << " }" << std::endl; }
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

		GMATDebug << " str: " << str << std::endl;
		//write out the propagate command line
		if (agmatstep.allowTheTimeStep2Vary) {
			GMATfile << "	Propagate 'Propagate " << agmatstep.myspacecraft->Name << "' " << str << " " 
				<< agmatstep.propagator.Name << "( " << agmatstep.myspacecraft->Name << " ) { " 
				<< agmatstep.myspacecraft->Name << ".ElapsedSecs = " << agmatstep.id << "_TimeStep" << " }" << std::endl;
		}
		else {
			GMATfile << "	Propagate 'Propagate " << agmatstep.myspacecraft->Name << "' " << str << " " 
				<< agmatstep.propagator.Name << "( " << agmatstep.myspacecraft->Name << " ) { " << agmatstep.myspacecraft->Name << ".ElapsedSecs = " << elapsed_secs << " }" << std::endl;
		}
	}
}


// method to write a GMAT PenUp Command
void gmatscripter::PenUp(){

	GMATfile << "	PenUp 'PenUp' " << GMATMission.missionbodies[0].central_body_name << "View";
	for (int body_index = 0; body_index < GMATMission.missionbodies.size(); ++body_index) {
		GMATfile << " " << GMATMission.missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << std::endl;

}


// method to write a GMAT PenDown Command
void gmatscripter::PenDown(){

	GMATfile << "	PenDown 'PenDown' " << GMATMission.missionbodies[0].central_body_name << "View";
	for (int body_index = 0; body_index < GMATMission.missionbodies.size(); ++body_index) {
		GMATfile << " " << GMATMission.missionbodies[body_index].name << "View";
	}
	GMATfile << ";" << std::endl;

}


// method to write a GMAT ForceModel Resource
void gmatscripter::create_GMAT_forcemodel(struct gmat_forcemodel& forcemodel){

	GMATfile << "Create ForceModel " << forcemodel.Name << std::endl;
	GMATfile << forcemodel.Name << ".CentralBody = " << forcemodel.CentralBody << std::endl;
	GMATfile << forcemodel.Name << ".PointMasses = {";
	for (int item = 0; item < forcemodel.PointMasses.size(); ++item) {
		if (item == forcemodel.PointMasses.size() - 1) { GMATfile << forcemodel.PointMasses[item]; }
		else { GMATfile << forcemodel.PointMasses[item] << ", "; }
	}
	GMATfile << "};" << std::endl;
	GMATfile << forcemodel.Name << ".Drag = None;" << std::endl;
	GMATfile << forcemodel.Name << ".SRP = Off;" << std::endl;
	GMATfile << std::endl;

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

	GMATfile << "Create Propagator " << propagator.Name << std::endl;
	GMATfile << propagator.Name << ".FM = " << propagator.ForceModel.Name << std::endl;
	GMATfile << propagator.Name << ".Type = PrinceDormand78; " << std::endl;
	GMATfile << propagator.Name << ".InitialStepSize = " << initialstepsize << std::endl;
	GMATfile << propagator.Name << ".Accuracy = 1e-11; " << std::endl;
	GMATfile << propagator.Name << ".MinStep = 0.0; " << std::endl;
	GMATfile << propagator.Name << ".MaxStep = " << maxstepsize << std::endl;
	GMATfile << std::endl;

}


// method to write a GMAT Spacecraft Resource
void gmatscripter::create_GMAT_spacecraft(struct gmat_spacecraft& spacecraft) {

	//write out the spacecraft information
	GMATfile << "Create Spacecraft " << spacecraft.Name << std::endl;
	GMATfile << spacecraft.Name << ".DateFormat = " << spacecraft.DateFormat << std::endl;
	GMATfile << spacecraft.Name << ".Epoch      = " << spacecraft.Epoch << std::endl;
	GMATfile << spacecraft.Name << ".DryMass    = " << spacecraft.DryMass << std::endl;
	GMATfile << spacecraft.Name << ".CoordinateSystem = " << spacecraft.CoordinateSystem << std::endl;
	GMATfile << spacecraft.Name << ".Tanks     = {" << spacecraft.Thruster.Tank.Name << "}" << std::endl;
	GMATfile << spacecraft.Name << ".Thrusters = {" << spacecraft.Thruster.Name << "}" << std::endl;
	//  initial conditions
	create_GMAT_initialconditions(spacecraft, "", false);
	GMATfile << std::endl;

}


// method to write a GMAT FuelTank Resource
void gmatscripter::create_GMAT_fueltank(struct gmat_spacecraft& spacecraft) {

	GMATfile << "Create FuelTank " << spacecraft.Thruster.Tank.Name << std::endl;
	GMATfile << spacecraft.Thruster.Tank.Name << ".AllowNegativeFuelMass = true" << std::endl;
	GMATfile << spacecraft.Thruster.Tank.Name << ".Volume = 10" << std::endl;
	GMATfile << spacecraft.Thruster.Tank.Name << ".FuelMass = " << spacecraft.Thruster.Tank.FuelMass << std::endl;
	GMATfile << std::endl;
	
}


// method to write a GMAT Thruster Resource
void gmatscripter::create_GMAT_thruster(struct gmat_spacecraft& spacecraft) {

	//GMATfile << "% Journey #" << j << ", Phase #" << p << ", Thruster Name: " << thrustername << std::endl;
	GMATfile << "% note: with the exception of .CoordinateSystem and .Tank, all values are 'dummy values'" << std::endl;
	GMATfile << "Create Thruster " << spacecraft.Thruster.Name << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".CoordinateSystem = " << spacecraft.CoordinateSystem << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustDirection1 = 1" << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustDirection2 = 0" << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustDirection3 = 0" << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".DutyCycle = 1" << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".Tank = " << spacecraft.Thruster.Tank.Name << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".ThrustScaleFactor = 1" << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".DecrementMass = true" << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".C1 = " << spacecraft.Thruster.C1 << std::endl;
	GMATfile << spacecraft.Thruster.Name << ".K1 = " << spacecraft.Thruster.K1 << std::endl;
	GMATfile << std::endl;

}


// method to write a GMAT FiniteBurn Resource
void gmatscripter::create_GMAT_burn(struct gmat_fburn& burn) {

	GMATfile << "Create FiniteBurn " << burn.Name << std::endl;
	GMATfile << burn.Name << ".Thrusters = {" << burn.ThrusterName << "};" << std::endl;

}

// method to write a GMAT ImpulsiveBurn Resource
void gmatscripter::create_GMAT_burn(struct gmat_iburn& burn) {

	if (burn.UseImpulsive) {
		GMATfile << "Create ImpulsiveBurn " << burn.Name << std::endl;
		GMATfile << burn.Name << ".CoordinateSystem = " << burn.CoordinateSystem << std::endl;
		GMATfile << "% " << burn.Name << ".Origin = " << burn.Origin << std::endl;
		GMATfile << "% " << burn.Name << ".Axes = " << burn.Axes << std::endl;
		GMATfile << burn.Name << ".Element1 = " << burn.Element1 << std::endl;
		GMATfile << burn.Name << ".Element2 = " << burn.Element2 << std::endl;
		GMATfile << burn.Name << ".Element3 = " << burn.Element3 << std::endl;

		GMATfile << burn.Name << ".DecrementMass = ";
		if (burn.DecrementMass) { GMATfile << "true" << std::endl; }
		else { GMATfile << "false" << std::endl; }

		GMATfile << burn.Name << ".Isp = " << burn.Isp << std::endl;
		GMATfile << burn.Name << ".GravitationalAccel = " << burn.g << std::endl;
		GMATfile << burn.Name << ".Tank = " << burn.TankName << std::endl;

		GMATfile << std::endl;
	}

}

// method to write a GMAT CoordinateSystem Resource
void gmatscripter::create_GMAT_coordinatesystem(string bodyname) {

	GMATfile << "Create CoordinateSystem " << bodyname << "J2000Eq;" << std::endl;
	GMATfile << bodyname << "J2000Eq.Origin = " << bodyname << ";" << std::endl;
	GMATfile << bodyname << "J2000Eq.Axes = MJ2000Eq;" << std::endl;
	GMATfile << std::endl;

}


// method to rewrite spacecraft initial conditions during the GMAT optimization segment
void gmatscripter::rewrite_GMAT_initialconditions(class gmatphase& agmatphase) {
	this->create_GMAT_initialconditions(agmatphase.spacecraft_forward, "   'Reset " + agmatphase.spacecraft_forward.Name + " Boundary Condition' ", true, agmatphase.reprint_forward_position_IC, agmatphase.reprint_forward_velocity_IC);
	this->create_GMAT_initialconditions(agmatphase.spacecraft_backward, "   'Reset " + agmatphase.spacecraft_backward.Name + " Boundary Condition' ", true, agmatphase.reprint_backward_position_IC, agmatphase.reprint_backward_velocity_IC);
}


// method to write out a spacecraft's initial conditions
void gmatscripter::create_GMAT_initialconditions(struct gmat_spacecraft& spacecraft, string specialchar, bool WithCoordinateSyst, 
												 bool PrintPosition, bool PrintVelocity) {
	
	std::string CoordinateSyst = "";
	if (WithCoordinateSyst) { CoordinateSyst = "." + spacecraft.CoordinateSystem; }

	if (PrintPosition) {
		GMATfile << specialchar << spacecraft.Name << CoordinateSyst << ".X  = " << spacecraft.initialconditions[0] << std::endl;
		GMATfile << specialchar << spacecraft.Name << CoordinateSyst << ".Y  = " << spacecraft.initialconditions[1] << std::endl;
		GMATfile << specialchar << spacecraft.Name << CoordinateSyst << ".Z  = " << spacecraft.initialconditions[2] << std::endl;
	}
	if (PrintVelocity) {
		GMATfile << specialchar << spacecraft.Name << CoordinateSyst << ".VX = " << spacecraft.initialconditions[3] << std::endl;
		GMATfile << specialchar << spacecraft.Name << CoordinateSyst << ".VY = " << spacecraft.initialconditions[4] << std::endl;
		GMATfile << specialchar << spacecraft.Name << CoordinateSyst << ".VZ = " << spacecraft.initialconditions[5] << std::endl;
	}

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
	number_of_emtg_phases = mymission->emtgmission->journeys[j].number_of_phases;
	if (j == 0) { isFirstJourney = true; }
	if (j == mymission->emtgmission->number_of_journeys - 1) { isLastJourney = true; }
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
		//summing the eligible / ineligible time for the 'gmatphase'
		if (this->mysteps[index].allowTheTimeStep2Vary) { this->eligibletime += this->mysteps[index].stepsize; }
		else { this->ineligibletime += this->mysteps[index].stepsize; }
		//summing the TOF for the 'gmatphase'
		this->TOF += this->mysteps[index].stepsize;
	}
}


} // end of EMTG namespace