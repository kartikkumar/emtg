/*
* GMATScripter.h
*
* Creation Date: july 1st 2014
* Author: ryne beeson
*
* Purpose: GMATScripter Class, which generates a GMAT script file
*		   following an EMTG run. Methods build on work conducted by Max Schadegg (summer 2013)
*		   and Jacob Englander
*/

#ifndef GMATSCRIPTER_H_
#define GMATSCRIPTER_H_

#include "mission.h"
#include "missionoptions.h"

#include <fstream> //ofstream
#include <vector>  
#include <string>


namespace EMTG {

class gmatscripter {

public:
	//constructors
	gmatscripter();
	gmatscripter(mission* mission_in);

	//destructor
	virtual ~gmatscripter();

	//methods
	//create file
	virtual void create_GMAT_file();
	//create data
	virtual void get_GMAT_bodieslist();
	virtual void get_GMAT_missionlevelparameters();
	virtual void get_GMAT_phaselevelparameters();
	//model setup
	virtual void write_GMAT_preamble();
	virtual void write_GMAT_spacecraft();
	virtual void write_GMAT_hardware();
	virtual void write_GMAT_nonstandardbody();
	virtual void write_GMAT_forcemodels();
	virtual void write_GMAT_propagators();
	virtual void write_GMAT_burns();
	virtual void write_GMAT_variables();
	virtual void write_GMAT_coordinatesystems();
	virtual void write_GMAT_solvers();
	virtual void write_GMAT_subscribers();
	//mission setup
	virtual void write_GMAT_initialguess();
	virtual void write_GMAT_initialboundaryconditions();
	virtual void write_GMAT_missionpropagate();
	virtual void write_GMAT_finalboundaryconditions();
	virtual void write_GMAT_objectivefunction();
	//reports
	virtual void write_GMAT_report(int j, int p, int s, string spacecraft_name, string body_name, 
								   bool isforwardspacecraft, bool isbeforemaneuver, bool writecontrolhistory);
	//GMAT Resource Methods
	virtual void aux_GMAT_forcemodel(string forcemodelname, string centralbody, string pointmasses);
	virtual void aux_GMAT_propagate(string propagatorname, string forcemodelname, bool isCloseApproach);
	//GMAT Command Methods
	virtual void aux_GMAT_beginburn(int j, int p, int s, string spacecraft_name, string prefix);
	virtual void aux_GMAT_endburn(int j, int p, int s, string spacecraft_name, string prefix);
	virtual void aux_GMAT_propagate(int j, int p, int s, string spacecraft_name, string prefix, string body_name, double elapsed_secs);
	virtual void aux_GMAT_penUp();
	virtual void aux_GMAT_penDown();
	


	//writeout the GMAT script
	virtual void write_GMAT_script();



	//members
	mission*  ptr_gmatmission;
	//missionoptions gmatoptions;
	std::ofstream GMATfile;
	//a temporary file for debugging purposes
	std::ofstream GMATDebug;
	
	vector <EMTG::Astrodynamics::body> missionbodies_unique;
	vector <EMTG::Astrodynamics::body> missionbodies;

	double LaunchDate_LowerBounds;
	double LaunchDate_UpperBounds;

	vector <double> Forward_Flyby_Distance_LowerBound;
	vector <double> Forward_Flyby_Distance_UpperBound;
	vector <double> Backward_Flyby_Distance_LowerBound;
	vector <double> Backward_Flyby_Distance_UpperBound;
	vector <double> Forward_Flyby_Velocity_LowerBound;
	vector <double> Forward_Flyby_Velocity_UpperBound;
	vector <double> Backward_Flyby_Velocity_LowerBound;
	vector <double> Backward_Flyby_Velocity_UpperBound;

	//vector of bool type for allowing simpler switching on/off of 
	//necessary 'force models', 'propagators', and mission sequence events in GMAT
	vector <bool> isSpaceCraftInSOI;
	vector <bool> useCentralBodyInSOI;

	//create a vector of strings for storing spacecraft names
	vector <string> spacecraft_names;
	vector <string> spacecraft_forward_names;
	vector <string> spacecraft_backward_names;


}; // end of class gmatscript

}  // end of EMTG namespace

#endif // end of GMATSCRIPTER_H_

