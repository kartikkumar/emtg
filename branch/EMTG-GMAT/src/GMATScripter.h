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
	virtual void get_GMAT_journeylevelparameters();
	virtual void get_GMAT_steplevelparameters();
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
	//optimization methods
	virtual void write_GMAT_beginmissionsequence();
	virtual void write_GMAT_initialguess();
	virtual void write_GMAT_optimization();

	//methods that "may" be deleted due to the new 'write_GMAT_optimization' method and approach - 2014_07_20
	virtual void write_GMAT_initialboundaryconditions();
	virtual void write_GMAT_missionpropagate();
	virtual void write_GMAT_finalboundaryconditions();
	virtual void write_GMAT_objectivefunction();
	
	//reports
	virtual void write_GMAT_report(int j, int p, int s, string spacecraft_name, string body_name, 
								   bool isforwardspacecraft, bool isbeforemaneuver, bool writecontrolhistory);
	//GMAT Resource Methods
	virtual void create_GMAT_forcemodel(string forcemodelname, string centralbody, string pointmasses);
	virtual void create_GMAT_propagator(string propagatorname, string forcemodelname, bool isCloseApproach);
	virtual void create_GMAT_spacecraft(int j, int p, string spacecraftname, string prefix, string thebody_coordinatesystem);
	virtual void create_GMAT_fueltank(int j, int p, string fueltankname, string prefix);
	virtual void create_GMAT_thruster(int j, int p, string thrustername, string thebody_coordinatesystem, string fueltankname);
	virtual void create_GMAT_finiteburn(string finiteburnname, string thrustername);
	virtual void create_GMAT_coordinatesystem(string bodyname);
	//GMAT Command Methods
	virtual void aux_GMAT_beginburn(string finiteburnobject, string spacecraft_name);
	//virtual void aux_GMAT_beginburn(int j, int p, int s, string spacecraft_name, string prefix);
	virtual void aux_GMAT_endburn(string finiteburnobject, string spacecraft_name);
	//virtual void aux_GMAT_endburn(int j, int p, int s, string spacecraft_name, string prefix);
	virtual void aux_GMAT_propagate(int j, int p, int s, string spacecraft_name, string prefix, string body_name, double elapsed_secs);
	virtual void aux_GMAT_penUp();
	virtual void aux_GMAT_penDown();

	//General Purpose Methods
	virtual void aux_GMAT_populate_thrustvector(int j, int p, int s, vector <double>& x, double lower_bound, double upper_bound);
	virtual void aux_GMAT_vary(string object2vary);
	virtual void aux_GMAT_calculate(string object2calculate, string rhs);
	virtual void aux_GMAT_nonlinearconstraint(string object2constrain, string relation, string rhs);
	


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


	//mission level parameters
	double LaunchDate_LowerBounds;
	double LaunchDate_UpperBounds;
	double LaunchWindow;
	double LaunchDate;
	double TAIModJOffset = -29999.5; // TAIModJOffset = 2400000.5 - 2430000;
	//bool for whether LT or Impulsive is being used
	bool isLT;
	//bool for whether FBLT or Sims-Flanagan is being used
	bool isFBLT;
	//variables to be used during mission sequence
	vector <vector <string>> mission_level_variables;

	//journey level parameters
	vector <vector <string>> journey_level_variables;

	//phase level parameters
	vector <double> Forward_Flyby_Distance_LowerBound;
	vector <double> Forward_Flyby_Distance_UpperBound;
	vector <double> Backward_Flyby_Distance_LowerBound;
	vector <double> Backward_Flyby_Distance_UpperBound;
	vector <double> Forward_Flyby_Velocity_LowerBound;
	vector <double> Forward_Flyby_Velocity_UpperBound;
	vector <double> Backward_Flyby_Velocity_LowerBound;
	vector <double> Backward_Flyby_Velocity_UpperBound;
	vector <vector <string>> phase_level_variables;

	//step level parameters
	vector < vector <string>> gmat_step_forcemodel_storage;
	vector < vector <string>> gmat_step_propagator_storage;
	vector < bool > gmat_step_propagator_isCloseApproach;
	vector <double> gmat_step_timesteps;
	vector < vector <double>> gmat_step_thrust_vectors;
	vector <string> gmat_step_thruster_names;
	vector <int> gmat_steps_per_phase;


	//vector of bool type for allowing simpler switching on/off of 
	//necessary 'force models', 'propagators', and mission sequence events in GMAT
	vector <bool> isSpaceCraftInSOI;
	vector <bool> useCentralBodyInSOI;

	//create a vector of strings for storing the spacecraft names for each phase x 2
	vector <string> spacecraft_names;
	vector <string> spacecraft_forward_names;
	vector <string> spacecraft_backward_names;
	//create a vector of strings for storing the fueltank names for each phase x 2
	vector <string> fueltank_forward_names;
	vector <string> fueltank_backward_names;
	//create a vector of strings for storing the thruster names for each phase x 2
	vector <string> thruster_forward_names;
	vector <string> thruster_backward_names;
	//create a vector of strings for storing the finiteburn names for each phase x 2
	vector <string> finiteburn_forward_names;
	vector <string> finiteburn_backward_names;

	//a struct type for gmat tank(s)
	struct gmat_tank {
		string Name;
		double FuelMass;
	};

	//a struct type for gmat thruster(s)
	struct gmat_thruster {
		string Name;
		struct gmat_tank Tank;
	};

	//a struct type for storing and accessing gmat spacecraft information
	struct gmat_spacecraft {
		//GMAT Specific Data
		string Name;
		string DateFormat = "TAIModJulian";
		double Epoch;
		double DryMass = 0.0;
		string CoordinateSystem;
		vector <struct gmat_tank> Tanks;
		vector <struct gmat_thruster> Thrusters;
		//Auxiliary Data
		bool isForward;
	};
	//TEMPORARY (2014_07_21) vector of gmat spacecraft(s)
	vector <struct gmat_spacecraft> spacecraft;

	//a struct type for forcemodel(s)
	struct gmat_forcemodel {
		string Name;
		string CentralBody;
		vector <string> PointMasses;
	};

	//a struct type for gmat propagator(s) 
	struct gmat_propagator {
		string Name;
		struct gmat_forcemodel ForceModel;
		bool isCloseApproach;
	};

	////a struct type for "gmat steps"
	//struct gmat_step {
	//	int j; //journey number
	//	int p; //phase number
	//	int gs; //gmat step number
	//	string identifier; //my identifier (e.g. 'j0p0gs0')
	//	struct gmat_spacecraft; //the spacecraft during my "gmat step"
	//	struct gmat_propagator; //the propagator to use for my "gmat step"
	//	string finiteburnname;
	//	vector <double> thrustvector;
	//	vector <string> vary; //collection of vary 'variables' to be performed in the optimize sequence. (n x 1)
	//	vector <vector <string>> calculate; //collection of calculate 'variables' to be performed in the optimize sequence. n x (2 x 1)
	//	vector <vector <string>> constraints; //collection of constraint 'variables' to be performed in the optimize sequence. n x (3 x 1)
	//};
	////vector of "gmat steps"
	vector <class gmatstep> gmat_steps;



}; // end of class gmatscript


class gmatstep {

public:
	//constructors
	gmatstep();
	gmatstep(int journey, int phase, int step, vector <EMTG::Astrodynamics::body>& bodies) {
		j = journey;
		p = phase;
		gs = step;
		longid = "j" + std::to_string(journey) + "p" + std::to_string(phase) + "gs" + std::to_string(step);
		shortid = "j" + std::to_string(journey) + "p" + std::to_string(phase);
		bodylist = bodies;
	};

	//destructor
	~gmatstep(){};

	//method
	void setNames(bool amIaForwardStep) {
		//declarations
		struct gmat_tank aTank;
		struct gmat_thruster aThruster;
		string append;

		//set the bool member of spacecraft, which indicates whether it is 'Forward' or 'Backward'
		this->spacecraft.isForward = amIaForwardStep;
		//assign 'append'
		if (amIaForwardStep) { append = "_Forward"; }
		else { append = "_BackWard"; }
		//name the spacecraft, tank, thruster, and finiteburnname
		this->spacecraft.Name = "SpaceCraft_" + this->shortid + append;
		aTank.Name = "FuelTank_" + this->shortid + append;
		aThruster.Tank = aTank;
		aThruster.Name = "Thruster_" + this->shortid + append;
		this->spacecraft.Thrusters.push_back(aThruster);
		this->finiteburnname = "FiniteBurn_" + this->shortid + append;
	}

	//method
	void setThrust(vector <double> tvector) {
		//declaration
		vector <string> tempvector;

		//store the thrust vector
		thrustvector = tvector;
		//create and store my thrust vector GMAT variable name
		thrustvectornames.push_back("ThrustVector_" + this->shortid + "(" + to_string(this->gs) + ", 0)");
		thrustvectornames.push_back("ThrustVector_" + this->shortid + "(" + to_string(this->gs) + ", 1)");
		thrustvectornames.push_back("ThrustVector_" + this->shortid + "(" + to_string(this->gs) + ", 2)");
		//store my thrust vector name, allowing my to vary during the optimization sequence
		vary.push_back(thrustvectornames[0]);
		vary.push_back(thrustvectornames[1]);
		vary.push_back(thrustvectornames[2]);
		//we will need to ensure that the thrust vector remains a unit vector during optimization, therefore 
		//we introduce a constraint. first do a calculation, then constrain it.
		tempvector.push_back("ThrustUnitVectorMagnitude_" + this->shortid + "(" + to_string(this->gs) + ", 1)");
		tempvector.push_back("sqrt(( " + thrustvectornames[0] + " * 2 - 1) ^ 2 + ( " + thrustvectornames[1] + " * 2 - 1) ^ 2 + ( " + thrustvectornames[2] + " * 2 - 1) ^ 2 )");
		calculate.push_back(tempvector);
		tempvector.pop_back();
		tempvector.push_back("<=");
		tempvector.push_back("1");
		constraints.push_back(tempvector);
	}

	//method
	void setFMandProp(string FMname, bool amIaCloseApproach) {
		//declarations
		int body_index;

		propagator.ForceModel.Name = FMname;
		propagator.Name = "Propagator_" + FMname;
		propagator.isCloseApproach = amIaCloseApproach;
		//based on 'isCloseApproach', my 'p', and whether 'spacecraft.isForward'
		//we can fill in the FM.PointMasses
		if (spacecraft.isForward) { body_index = this->p; }
		else { body_index = p + 1; }
		propagator.ForceModel.CentralBody = bodylist[body_index].central_body_name;
		if (amIaCloseApproach) {
			propagator.ForceModel.PointMasses.push_back("");
		}
		else {
			if (spacecraft.isForward) {
				propagator.ForceModel.PointMasses.push_back(bodylist[body_index].name);
				propagator.ForceModel.PointMasses.push_back(bodylist[body_index + 1].name);
			}
			else {
				propagator.ForceModel.PointMasses.push_back(bodylist[body_index - 1].name);
				propagator.ForceModel.PointMasses.push_back(bodylist[body_index].name);
			}
		}
	}

	//a struct type for gmat tank(s)
	struct gmat_tank {
		string Name;
		double FuelMass;
	};

	//a struct type for gmat thruster(s)
	struct gmat_thruster {
		string Name;
		struct gmat_tank Tank;
	};

	//a struct type for storing and accessing gmat spacecraft information
	struct gmat_spacecraft {
		//GMAT Specific Data
		string Name;
		string DateFormat = "TAIModJulian";
		double Epoch;
		double DryMass = 0.0;
		string CoordinateSystem;
		vector <struct gmat_thruster> Thrusters;
		//Auxiliary Data
		bool isForward;
	};

	//a struct type for forcemodel(s)
	struct gmat_forcemodel {
		string Name;
		string CentralBody;
		vector <string> PointMasses;
	};

	//a struct type for gmat propagator(s) 
	struct gmat_propagator {
		string Name;
		struct gmat_forcemodel ForceModel;
		bool isCloseApproach;
	};

	//members
	int j; //journey number
	int p; //phase number
	int gs; //gmat step number
	string longid; //my identifier (e.g. 'j0p0gs0')
	string shortid; //my shortened identifier (e.g. 'j0p0') 
	vector <EMTG::Astrodynamics::body> bodylist;
	struct gmat_spacecraft spacecraft; //the spacecraft during my "gmat step"
	struct gmat_propagator propagator; //the propagator to use for my "gmat step"
	string finiteburnname;
	vector <double> thrustvector;
	vector <string> thrustvectornames;
	double initial_timestep;
	vector <string> vary; //collection of vary 'variables' to be performed in the optimize sequence. (n x 1)
	vector <vector <string>> calculate; //collection of calculate 'variables' to be performed in the optimize sequence. n x (2 x 1)
	vector <vector <string>> constraints; //collection of constraint 'variables' to be performed in the optimize sequence. n x (3 x 1)

};

//void EMTG::gmatstep::setNames(bool amIaForwardStep) {
//	//declarations
//	struct gmat_tank aTank;
//	struct gmat_thruster aThruster;
//	string append;
//
//	//set the bool member of spacecraft, which indicates whether it is 'Forward' or 'Backward'
//	this->spacecraft.isForward = amIaForwardStep;
//	//assign 'append'
//	if (amIaForwardStep) { append = "_Forward";	}
//	else { append = "_BackWard"; }
//	//name the spacecraft, tank, thruster, and finiteburnname
//	this->spacecraft.Name = "SpaceCraft_" + this->shortid + append;
//	aTank.Name = "FuelTank_" + this->shortid + append;
//	aThruster.Tank = aTank;
//	aThruster.Name = "Thruster_" + this->shortid + append;
//	this->spacecraft.Thrusters.push_back(aThruster);
//	this->finiteburnname = "FiniteBurn_" + this->shortid + append;
//
//}

}  // end of EMTG namespace

#endif // end of GMATSCRIPTER_H_

