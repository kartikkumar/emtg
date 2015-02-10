/*
* GMATScripter.h
*
* Creation Date: july 1st 2014
* Author: ryne beeson
*
* Purpose: GMATScripter Class, which generates a GMAT script file
*		   following an EMTG run. Methods build on work conducted by Max Schadegg (summer 2013)
*		   and Jacob Englander
*
* Features that "should (i.e marginally tested)" be working as of 2014_08_01:
* - multiphase
* - multijourney (single-phase)
* - objective function (maximize final mass)
* - departures( launch direct-departure, free departure)
* - arrivals ( low-thrust rendezvous)
* - engine models (fixed isp / fixed thrust)
* - launch window and journey wait times enforced
* - bounded mission TOF enforced
* - bound journey TOF enforced
* - variable phase TOF
* - variable flyby position and velocity with bounds enabled
*
*/

#ifndef GMATSCRIPTER_H_
#define GMATSCRIPTER_H_

#include "mission.h"
#include "missionoptions.h"
#include "Astrodynamics.h"
#include "EMTG_math.h"

#include <fstream> //ofstream
#include <sstream> //stringstream
#include <vector>  
#include <string>


namespace EMTG {

	//declarations
	class gmatscripter;
	class gmatbaseclass;
	class gmatmission;
	class gmatjourney;
	class gmatphase;
	class gmatstep;


//a struct type for gmat tank(s)
struct gmat_electric_tank {
	string Name;
	double FuelMass;
};

//a struct type for gmat thruster(s)
struct gmat_electric_thruster {
	string Name;
	struct gmat_electric_tank Tank;
    double DutyCycle = 0.90;
    double ThrustScaleFactor = 1.0;
    double GravitationalAccel = 9.80665;
    string ThrustModel = "ThrustMassPolynomial";
    double MaximumUsablePower = 7.4;
    double MinimumUsablePower = 0.31;
    double Isp = 3219.12314;
    double ConstantThrust = 0.01243;
    double ThrustCoeff1 = 1.92e-006;
    double ThrustCoeff2 = 54.05382;
    double ThrustCoeff3 = -14.41789;
    double ThrustCoeff4 = 2.96519;
    double ThrustCoeff5 = -0.19082;
    double MassFlowCoeff1 = 2.13781;
    double MassFlowCoeff2 = 0.03211;
    double MassFlowCoeff3 = -0.09956;
    double MassFlowCoeff4 = 0.05717;
    double MassFlowCoeff5 = -0.004776;
    double FixedEfficiency = 0.654321;
};

//a struct type for gmat power system(s)
struct gmat_power_system {
    string Name;
    double InitialEpoch = 21547.00000039794;
    double InitialMaxPower = 1.2124;
    double AnnualDecayRate = 5.123;
    double Margin = 4.998;
    double BusCoeff1 = 0.32;
    double BusCoeff2 = 0.0001;
    double BusCoeff3 = 0.0001;
    int ShadowModel = 0; //switch to enable different GMAT shadow models
    double SolarCoeff1 = 1.33;
    double SolarCoeff2 = -0.11;
    double SolarCoeff3 = -0.12;
    double SolarCoeff4 = 0.11;
    double SolarCoeff5 = -0.02;
};

//a struct type for gmat finite burn
struct gmat_fburn {
	string Name;
	string Type = "FiniteBurn";
	string ThrusterName;
};

//a struct type for gmat impulsive burn
struct gmat_iburn {
	bool UseImpulsive = false;
	bool BCisParkingOrbit = false;
	bool IsEDS = false;
	string Name;
	string Type = "ImpulsiveBurn";
	string CoordinateSystem;
	string Origin;
	string Axes = "MJ2000Eq";
	double Element1 = 0.0;
	double Element2 = 0.0;
	double Element3 = 0.0;
	bool DecrementMass = true;
	double Isp;
    double g = 9.80665;
	double c;
	string TankName;
};

//a struct type for storing and accessing gmat spacecraft information
struct gmat_spacecraft {
	//GMAT Specific Data
	string Name;
	string DateFormat = "TAIModJulian";
	double Epoch;
	double DryMass = 0.0;
	string CoordinateSystem;
	struct gmat_electric_thruster Thruster;
    struct gmat_power_system PowerSystem;
	struct gmat_fburn fBurn;
	struct gmat_iburn iBurn;
	//Auxiliary Data
	bool isForward;
	double flyby_distance_lowerbound = 0;
	double flyby_distance_upperbound = 1e10;
	double flyby_velocity_lowerbound = 0;
	double flyby_velocity_upperbound = 1e10;
	math::Matrix<double> flyby_states;
	double initialconditions[6];
	double mass_margin = 0.0;
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

//a struct type for gmat calculate information
struct gmat_calculate {
	string message;
	string lhs;
	string rhs;
	bool writeAtTheFront = true;
};

//a struct type for gmat vary information
struct gmat_vary {
	string object2vary;
	double perturbation;
	double lowerbound;
	double upperbound;
	double maxstep;
};

//Class() 'gmatbaseclass'
class gmatbaseclass {

public:
	//constructor
	gmatbaseclass(){};

	//desstructor
	~gmatbaseclass(){};

	//collection members
	vector <vector<string> > variables;					//collection of 'variables' (n x (1 or 2))
	vector <struct gmat_vary> vary;						//collection of vary 'variables' to be performed in the optimize sequence.           
	vector <struct gmat_calculate> calculate;			//collection of calculate 'variables' to be performed in the optimize sequence.  
	vector <vector<string> > constraints;				//collection of constraint 'variables' to be performed in the optimize sequence. n x (3 x 1)
	vector <struct gmat_calculate> preopt_calculate;	//collection of calculate 'variables' to be performed prior to the optimization sequence.  

	//method
	void setVariable(string avariable) {
		//declarations
		vector <string> tempvector;
		//function
		tempvector.push_back(avariable);
		this->variables.push_back(tempvector);
	}
	
	//method
	void setVariable(string avariable, string anassignment) {
		//declarations
		vector <string> tempvector;
		//function
		tempvector.push_back(avariable);
		tempvector.push_back(" = " + anassignment);
		this->variables.push_back(tempvector);
	}

	//method 
	void setVariable(string avariable, double anassignment) {
		//declarations
		vector <string> tempvector;
		stringstream tempstream;
		//function
		tempstream.precision(25);
		tempstream << anassignment;
		tempvector.push_back(avariable);
		tempvector.push_back(" = " + tempstream.str());
		this->variables.push_back(tempvector);
	}

	//method
	void setVary(string object2vary, double perturbation = 0.00001, double lowerbound = 0.0, double upperbound = 1.0, double maxstep = 0.1) {
		//declaration
		struct gmat_vary astruct;
		//function
		astruct.object2vary = object2vary;
		astruct.perturbation = perturbation;
		astruct.lowerbound = lowerbound;
		astruct.upperbound = upperbound;
		astruct.maxstep = maxstep;
		this->vary.push_back(astruct);
		//  GMAT does not have control over a default lower and upper bound with VF13ad,
		//+ therefore we must create constraints here
		this->setConstraint(object2vary, ">=", lowerbound);
		this->setConstraint(object2vary, "<=", upperbound);
	}

	//method
	void setCalculate(string message, string lhs, string rhs, bool writeAtTheFront = true) {
		//declarations
		struct gmat_calculate astruct;
		//function
		astruct.message = message;
		astruct.lhs = lhs;
		astruct.rhs = rhs;
		astruct.writeAtTheFront = writeAtTheFront;
		this->calculate.push_back(astruct);
	}

	//method
	void setCalculate(string lhs, string rhs, bool writeAtTheFront = true) {
		//  call the other variant of setCalculate()
		this->setCalculate(lhs, lhs, rhs, writeAtTheFront);
	}

	//method
	void setConstraint(string lhs, string relation, double rhs) {
		//declarations
		vector <string> tempvector;
		stringstream tempstream;
		//function
		tempstream.precision(25);
		tempstream << rhs;
		tempvector.push_back(lhs);
		tempvector.push_back(lhs);
		tempvector.push_back(relation);
		tempvector.push_back(tempstream.str());
		this->constraints.push_back(tempvector);
	}

	//method
	void setConstraint(string lhs, string relation, string rhs) {
		//declarations
		vector <string> tempvector;
		//function
		tempvector.push_back(lhs);
		tempvector.push_back(lhs);
		tempvector.push_back(relation);
		tempvector.push_back(rhs);
		this->constraints.push_back(tempvector);
	}

	//method
	void setConstraint(string message, string lhs, string relation, double rhs) {
		//declarations
		vector <string> tempvector;
		stringstream tempstream;
		//function
		tempstream.precision(25);
		tempstream << rhs;
		tempvector.push_back(message);
		tempvector.push_back(lhs);
		tempvector.push_back(relation);
		tempvector.push_back(tempstream.str());
		this->constraints.push_back(tempvector);
	}
	
	//method
	void setConstraint(string message, string lhs, string relation, string rhs) {
		//declarations
		vector <string> tempvector;
		//function
		tempvector.push_back(message);
		tempvector.push_back(lhs);
		tempvector.push_back(relation);
		tempvector.push_back(rhs);
		this->constraints.push_back(tempvector);
	}

	//method
	void setPreOptimizationCalculation(string message, string lhs, string rhs, bool writeAtTheFront = true) {
		//declarations
		struct gmat_calculate astruct;
		//function
		astruct.message = message;
		astruct.lhs = lhs;
		astruct.rhs = rhs;
		astruct.writeAtTheFront = writeAtTheFront; // not currently being used for PreOpt (2014_11_25)
		this->preopt_calculate.push_back(astruct);
	}
	
	//method
	void setPreOptimizationCalculation(string lhs, string rhs, bool writeAtTheFront = true) {
		//  call the other variant of setCalculate()
		this->setPreOptimizationCalculation(lhs, lhs, rhs, writeAtTheFront);
	}
	
	//method to write a GMAT 'Variable' line
	void printVariable(std::ofstream& File) {
		for (int index = 0; index < this->variables.size(); ++index) {
			File << "Create Variable " << this->variables[index][0] << std::endl;
			if (this->variables[index].size() == 2) { File << variables[index][0] << variables[index][1] << std::endl; }
		}
	}

	//method to write a GMAT 'Vary' line
	void printVary(std::ofstream& File) {
		for (int index = 0; index < this->vary.size(); ++index) {
			File << "	" << "Vary 'Vary " << this->vary[index].object2vary << "' NLPObject(" << this->vary[index].object2vary << " = " << this->vary[index].object2vary << 
					", {Perturbation = " << this->vary[index].perturbation << 
					", Lower = " << this->vary[index].lowerbound << 
					", Upper = " << this->vary[index].upperbound << 
					", MaxStep = " << this->vary[index].maxstep << " })" << std::endl;
		}
	}

	//method to write a GMAT 'Calculate' line (int mode (1) if printing at "front" and (0) at "end") (see write_GMAT_optimization for e.g.)
	void printCalculate(std::ofstream& File, int mode = 1) {
		for (int index = 0; index < this->calculate.size(); ++index) {
			if (mode && this->calculate[index].writeAtTheFront) { File << "   " << "'Calculate " << this->calculate[index].message << "' " << this->calculate[index].lhs << " = " << this->calculate[index].rhs << std::endl; }
			else if (!mode && !this->calculate[index].writeAtTheFront) { File << "   " << "'Calculate " << this->calculate[index].message << "' " << this->calculate[index].lhs << " = " << this->calculate[index].rhs << std::endl; }
		}
	}

	//method to write a GMAT 'Calculate' line prior to the start of the optimization sequence
	void printPreOptCalculate(std::ofstream& File) {
		for (int index = 0; index < this->preopt_calculate.size(); ++index) {
			File << "   " << "'Calculate " << this->preopt_calculate[index].message << "' " << this->preopt_calculate[index].lhs << " = " << this->preopt_calculate[index].rhs << std::endl;
		}
	}

	// method to write a GMAT 'NonlinearConstraint' line
	void printConstraint(std::ofstream& File) {
		for (int index = 0; index < this->constraints.size(); ++index) {
			File << "   " << "NonlinearConstraint '" << this->constraints[index][0] << "' NLPObject( " << this->constraints[index][1] << " " << this->constraints[index][2] << " " << this->constraints[index][3] << " )" << std::endl;
		}
	}

};


//Class() 'gmatmission'
class gmatmission: public EMTG::gmatbaseclass {

public:
	//constructors
	gmatmission();
	gmatmission(mission* anemtgmission);

	//destructors
	~gmatmission();

	//members
	mission* emtgmission;
	vector <EMTG::Astrodynamics::body> missionbodies_unique;
	vector <EMTG::Astrodynamics::body> missionbodies;
	vector <int> starting_body_index_of_journey;
	bool isLT;    //bool for whether LT   or Impulsive is being used
	bool isFBLT;  //bool for whether FBLT or Sims-Flanagan is being used
	vector <gmatjourney> myjourneys;

	//method
	void get_mission_bodies() {

		//declaration
		int start_int = 0;
		int counter   = 0;

		//push_back()
		starting_body_index_of_journey.push_back(counter);

		//a loop structure over journeys and phases to collect all the bodies for the mission 
		for (int j = 0; j < this->emtgmission->number_of_journeys; ++j) {
			if (j > 0) { start_int = 1; }
			//add bodies of each phase for each journey to the missionbodies vector
			for (int p = start_int; p < this->emtgmission->journeys[j].number_of_phases + 1; ++p) {
				//push back body onto the vector
				missionbodies.push_back(this->emtgmission->TheUniverse[j].bodies[this->emtgmission->options.sequence[j][p] - 1]);
				counter++;
			}//end of phase for-statement
			starting_body_index_of_journey.push_back(counter);
		}//end of journeys for-statement

		int body_flag = 0;
		//remove duplicate entries in visited bodies list
		for (int body_index = 0; body_index < missionbodies.size(); ++body_index) {
			body_flag = 0;
			for (int body_index_unique = 0; body_index_unique < missionbodies_unique.size(); ++body_index_unique) {
				if (missionbodies[body_index].spice_ID == missionbodies_unique[body_index_unique].spice_ID) {
					body_flag = 1;
				}
			}
			//if body flag not switch 'on', then the body is unique; add it to the missionbodies_unique vector
			if (body_flag == 0) {
				missionbodies_unique.push_back(missionbodies[body_index]);
			}
		}

		//make sure that each body doesnt start with a number for GMAT's sake
		for (int body_index = 0; body_index < missionbodies.size(); ++body_index) {
			if ((missionbodies[body_index].name[0] == 0) || (missionbodies[body_index].name[0] == 1) || (missionbodies[body_index].name[0] == 2) || (missionbodies[body_index].name[0] == 3) || (missionbodies[body_index].name[0] == 4) || (missionbodies[body_index].name[0] == 5) || (missionbodies[body_index].name[0] == 6) || (missionbodies[body_index].name[0] == 7) || (missionbodies[body_index].name[0] == 8) || (missionbodies[body_index].name[0] == 9)) {
				missionbodies[body_index].name = "A" + missionbodies[body_index].name;
			}
		}
		for (int body_index = 0; body_index < missionbodies_unique.size(); ++body_index) {
			if ((missionbodies_unique[body_index].name[0] == 0) || (missionbodies_unique[body_index].name[0] == 1) || (missionbodies_unique[body_index].name[0] == 2) || (missionbodies_unique[body_index].name[0] == 3) || (missionbodies_unique[body_index].name[0] == 4) || (missionbodies_unique[body_index].name[0] == 5) || (missionbodies_unique[body_index].name[0] == 6) || (missionbodies_unique[body_index].name[0] == 7) || (missionbodies_unique[body_index].name[0] == 8) || (missionbodies_unique[body_index].name[0] == 9)) {
				missionbodies_unique[body_index].name = "A" + missionbodies_unique[body_index].name;
			}
		}
	}//end of get_mission_bodies() method

	//method
	void setMissionThrustType() {
		//what type of mission is being solved
		//TODO: What about mission type 5? and confirm mission type 6-9 results in a mission type of 0-4.
		if (this->emtgmission->options.mission_type == 0 || this->emtgmission->options.mission_type == 1 || this->emtgmission->options.mission_type == 4) {
			isLT = false;
			isFBLT = false;
		}
		else if (this->emtgmission->options.mission_type == 2) {
			isLT = true;
			isFBLT = false;
		}
		else if (this->emtgmission->options.mission_type == 3) {
			isLT = true;
			isFBLT = true;
		}
	}

};


//Class() 'gmatjourney'
class gmatjourney: public EMTG::gmatbaseclass {

public:
	//constructors
	gmatjourney();
	gmatjourney(gmatmission* amission, int journey);

	//destructors
	~gmatjourney();

	//members
	gmatmission* mymission;
	int j;
	string id;
	bool isFirstJourney = false;
	bool isLastJourney  = false;
	int number_of_emtg_phases;
	vector <gmatphase> myphases;

};


//Class() 'gmatphase'
class gmatphase: public EMTG::gmatbaseclass {

public:
	//constructors
	gmatphase();
	gmatphase(gmatjourney* ajourney, int phase) {
		myjourney = ajourney;
		p = phase;
		id = myjourney->id + "p" + std::to_string(p);
		if (p == 0) { isFirstPhase = true; }
		if (p == myjourney->number_of_emtg_phases - 1) { isLastPhase = true; }
		this->set_names();
		this->set_fuelmass_epoch_drymass();
		this->get_my_bodies();
        this->set_power_system();
		this->set_thruster();
		this->set_iBurn();
		this->get_flyby_data();
		this->set_initialconditions();
		this->set_epochs();
	}

	//destructor
	~gmatphase(){};

	//members
	gmatjourney* myjourney;
	int p;
	string id;
	struct gmat_spacecraft spacecraft_forward;
	struct gmat_spacecraft spacecraft_backward;
	vector <EMTG::Astrodynamics::body> mybodies;
	bool isFirstPhase = false;
	bool isLastPhase  = false;
	bool StartsWithFlyby = false;
	bool EndsWithFlyby   = false;
	double ineligibletime = 0.0;
	double eligibletime   = 0.0;
	double TOF = 0.0; 
	double iepoch;
	double fepoch;
	vector <gmatstep> mysteps;
	bool reprint_forward_position_IC = false;
	bool reprint_forward_velocity_IC = false;
	bool reprint_backward_position_IC = false;
	bool reprint_backward_velocity_IC = false;

	//method
	void set_names() {
		//assign names to my members
		spacecraft_forward.Name = "SpaceCraft_" + this->id + "_Forward";
        spacecraft_forward.PowerSystem.Name = "PowerSystem_" + this->id + "_Forward";
		spacecraft_forward.Thruster.Name = "Thruster_" + this->id + "_Forward";
		spacecraft_forward.Thruster.Tank.Name = "FuelTank_" + this->id + "_Forward";
		spacecraft_forward.fBurn.Name = "FiniteBurn_" + this->id + "_Forward";
		spacecraft_forward.fBurn.ThrusterName = spacecraft_forward.Thruster.Name;
		spacecraft_forward.iBurn.Name = "ImpulsiveBurn_" + this->id + "_Forward";
		spacecraft_forward.isForward = true;

		spacecraft_backward.Name = "SpaceCraft_" + this->id + "_Backward";
        spacecraft_backward.PowerSystem.Name = "PowerSystem_" + this->id + "_Backward";
		spacecraft_backward.Thruster.Name = "Thruster_" + this->id + "_Backward";
		spacecraft_backward.Thruster.Tank.Name = "FuelTank_" + this->id + "_Backward";
		spacecraft_backward.fBurn.Name = "FiniteBurn_" + this->id + "_Backward";
		spacecraft_backward.fBurn.ThrusterName = spacecraft_backward.Thruster.Name;
		spacecraft_backward.iBurn.Name = "ImpulsiveBurn_" + this->id + "_Backward";
		spacecraft_backward.isForward = false;
	}

	//method
	void set_fuelmass_epoch_drymass() {
		//assignment of relevant data from emtg mission
		double starting_mass = this->myjourney->mymission->emtgmission->journeys[myjourney->j].phases[p].state_at_beginning_of_phase[6];
		double final_mass = this->myjourney->mymission->emtgmission->journeys[myjourney->j].phases[p].state_at_end_of_phase[6];
		double drymass_min = this->myjourney->mymission->emtgmission->options.final_mass_constraint;
		double drymass_max_increment = this->myjourney->mymission->emtgmission->options.journey_starting_mass_increment[this->myjourney->j]; // this->myjourney->mymission->emtgmission->journeys[myjourney->j].phases[p].current_mass_increment;
		double drymass_increment_scaling = 1.0;
		//find the optimal drymass scaling for journey 'j'
		if (drymass_max_increment > 0.0) {
			for (size_t iX = 0; iX < this->myjourney->mymission->emtgmission->Xdescriptions.size(); ++iX) {
				//find index in Xdescriptions where the launch time bounds are located
				if (this->myjourney->mymission->emtgmission->Xdescriptions[iX] == "j" + std::to_string(this->myjourney->j) + "p0: journey initial mass scale factor") {
					drymass_increment_scaling = this->myjourney->mymission->emtgmission->Xopt[iX];
					break;
				}
			}
		}
		//forward spacecraft drymass
		if (this->myjourney->j == 0) {
			spacecraft_forward.DryMass = drymass_min + drymass_increment_scaling * drymass_max_increment;
		}
		else {
			spacecraft_forward.DryMass = this->myjourney->mymission->myjourneys[this->myjourney->j - 1].myphases.back().spacecraft_forward.DryMass +
										 drymass_increment_scaling * drymass_max_increment;
		}
		//forward spacecraft mass margin (i.e. difference between drymass and propellant maximum)
		if (this->myjourney->j == 0 && this->myjourney->mymission->emtgmission->options.enable_maximum_propellant_mass_constraint) {
			spacecraft_forward.mass_margin = this->myjourney->mymission->emtgmission->journeys[0].phases[0].state_at_beginning_of_phase[6]
											 - spacecraft_forward.DryMass
											 - this->myjourney->mymission->emtgmission->options.maximum_propellant_mass;
		}
		else if (this->myjourney->mymission->emtgmission->options.enable_maximum_propellant_mass_constraint) {
			spacecraft_forward.mass_margin = this->myjourney->mymission->myjourneys[this->myjourney->j - 1].myphases[0].spacecraft_forward.mass_margin;
		}
		//correct forward 'DryMass' when mass margin is positive
		if (this->myjourney->j == 0 && spacecraft_forward.mass_margin > 0.0) {
			spacecraft_forward.DryMass += spacecraft_forward.mass_margin;
		}
		//forward spacecraft fuelmass
		spacecraft_forward.Thruster.Tank.FuelMass = starting_mass - spacecraft_forward.DryMass;
		if (spacecraft_forward.Thruster.Tank.FuelMass < 0.0) { spacecraft_forward.Thruster.Tank.FuelMass = 0.0; }
		//forward spacecraft epoch
		spacecraft_forward.Epoch = (this->myjourney->mymission->emtgmission->journeys[myjourney->j].phases[p].phase_start_epoch / 86400.0) + 2400000.5 - 2430000;
		//backward spacecraft parameters
		spacecraft_backward.Epoch = (this->myjourney->mymission->emtgmission->journeys[myjourney->j].phases[p].phase_end_epoch / 86400.0) + 2400000.5 - 2430000;
		spacecraft_backward.DryMass = spacecraft_forward.DryMass;
		spacecraft_backward.Thruster.Tank.FuelMass = final_mass - spacecraft_backward.DryMass;
		if (spacecraft_backward.Thruster.Tank.FuelMass < 0.0) { spacecraft_backward.Thruster.Tank.FuelMass = 0.0; }
	}

	//method
	void get_my_bodies() {
		//declaration
		int start_index;

		if (myjourney->j > 0) { start_index = this->myjourney->mymission->starting_body_index_of_journey[myjourney->j] - 1; }
		else { start_index = 0; }

		//get the bodies
		mybodies.push_back(this->myjourney->mymission->missionbodies[p + start_index]);
		mybodies.push_back(this->myjourney->mymission->missionbodies[p + 1 + start_index]);
		//set the spacecraft coordinate systems
		spacecraft_forward.CoordinateSystem  = mybodies[0].name + "J2000Eq";
		spacecraft_backward.CoordinateSystem = mybodies[1].name + "J2000Eq";
	}

    //set the power system
    void set_power_system() {
        this->spacecraft_forward.PowerSystem.InitialMaxPower = this->myjourney->mymission->emtgmission->options.power_at_1_AU;
        this->spacecraft_forward.PowerSystem.InitialEpoch = 21547.00000039794;
        this->spacecraft_forward.PowerSystem.AnnualDecayRate = this->myjourney->mymission->emtgmission->options.power_decay_rate * 100.0;
        this->spacecraft_forward.PowerSystem.Margin = this->myjourney->mymission->emtgmission->options.power_margin * 100.0;
        this->spacecraft_forward.PowerSystem.BusCoeff1 = this->myjourney->mymission->emtgmission->options.spacecraft_power_coefficients[0];
        this->spacecraft_forward.PowerSystem.BusCoeff2 = this->myjourney->mymission->emtgmission->options.spacecraft_power_coefficients[1];
        this->spacecraft_forward.PowerSystem.BusCoeff3 = this->myjourney->mymission->emtgmission->options.spacecraft_power_coefficients[2];
        this->spacecraft_forward.PowerSystem.ShadowModel = 0; //switch to enable different GMAT shadow models
        this->spacecraft_forward.PowerSystem.SolarCoeff1 = this->myjourney->mymission->emtgmission->options.solar_power_gamma[0];
        this->spacecraft_forward.PowerSystem.SolarCoeff2 = this->myjourney->mymission->emtgmission->options.solar_power_gamma[1];
        this->spacecraft_forward.PowerSystem.SolarCoeff3 = this->myjourney->mymission->emtgmission->options.solar_power_gamma[2];
        this->spacecraft_forward.PowerSystem.SolarCoeff4 = this->myjourney->mymission->emtgmission->options.solar_power_gamma[3];
        this->spacecraft_forward.PowerSystem.SolarCoeff5 = this->myjourney->mymission->emtgmission->options.solar_power_gamma[4];
        //copy the forward power system to the backward power system
        this->spacecraft_backward.PowerSystem.InitialMaxPower = this->spacecraft_forward.PowerSystem.InitialMaxPower;
        this->spacecraft_backward.PowerSystem.InitialEpoch = this->spacecraft_forward.PowerSystem.InitialEpoch;
        this->spacecraft_backward.PowerSystem.AnnualDecayRate = this->spacecraft_forward.PowerSystem.AnnualDecayRate;
        this->spacecraft_backward.PowerSystem.Margin = this->spacecraft_forward.PowerSystem.Margin;
        this->spacecraft_backward.PowerSystem.BusCoeff1 = this->myjourney->mymission->emtgmission->options.spacecraft_power_coefficients[0];
        this->spacecraft_backward.PowerSystem.BusCoeff2 = this->myjourney->mymission->emtgmission->options.spacecraft_power_coefficients[1];
        this->spacecraft_backward.PowerSystem.BusCoeff3 = this->myjourney->mymission->emtgmission->options.spacecraft_power_coefficients[2];
        this->spacecraft_backward.PowerSystem.ShadowModel = this->spacecraft_forward.PowerSystem.ShadowModel;
        this->spacecraft_backward.PowerSystem.SolarCoeff1 = this->spacecraft_forward.PowerSystem.SolarCoeff1;
        this->spacecraft_backward.PowerSystem.SolarCoeff2 = this->spacecraft_forward.PowerSystem.SolarCoeff2;
        this->spacecraft_backward.PowerSystem.SolarCoeff3 = this->spacecraft_forward.PowerSystem.SolarCoeff3;
        this->spacecraft_backward.PowerSystem.SolarCoeff4 = this->spacecraft_forward.PowerSystem.SolarCoeff4;
        this->spacecraft_backward.PowerSystem.SolarCoeff5 = this->spacecraft_forward.PowerSystem.SolarCoeff5;
    }

	//  set the engine_type
	void set_thruster() {
		// fixed thrust/Isp
		if (this->myjourney->mymission->emtgmission->options.engine_type == 0) 
        {
            spacecraft_forward.Thruster.ThrustModel = "ConstantThrustAndIsp";
			spacecraft_forward.Thruster.ConstantThrust = this->myjourney->mymission->emtgmission->options.Thrust;
			spacecraft_forward.Thruster.Isp= this->myjourney->mymission->emtgmission->options.IspLT;
            spacecraft_backward.Thruster.ThrustModel = "ConstantThrustAndIsp";
            spacecraft_backward.Thruster.ConstantThrust = this->myjourney->mymission->emtgmission->options.Thrust;
            spacecraft_backward.Thruster.Isp = this->myjourney->mymission->emtgmission->options.IspLT;
		}
        //fixed Isp/efficiency
        else if (this->myjourney->mymission->emtgmission->options.engine_type == 3)
        {
            spacecraft_forward.Thruster.ThrustModel = "FixedEfficiency";
            spacecraft_forward.Thruster.FixedEfficiency = this->myjourney->mymission->emtgmission->options.user_defined_engine_efficiency;
            spacecraft_forward.Thruster.Isp = this->myjourney->mymission->emtgmission->options.IspLT;
            spacecraft_backward.Thruster.ThrustModel = "FixedEfficiency";
            spacecraft_backward.Thruster.FixedEfficiency = this->myjourney->mymission->emtgmission->options.user_defined_engine_efficiency;
            spacecraft_backward.Thruster.Isp = this->myjourney->mymission->emtgmission->options.IspLT;
        }
        //custom thruster coefficients
        else if (this->myjourney->mymission->emtgmission->options.engine_type == 5) 
        {
            spacecraft_forward.Thruster.ThrustModel = "ThrustMassPolynomial";
            spacecraft_forward.Thruster.MinimumUsablePower = this->myjourney->mymission->emtgmission->options.engine_input_power_bounds[0];
            spacecraft_forward.Thruster.MaximumUsablePower = this->myjourney->mymission->emtgmission->options.engine_input_power_bounds[1];
            spacecraft_forward.Thruster.ThrustCoeff1 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[0];
            spacecraft_forward.Thruster.ThrustCoeff2 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[1];
            spacecraft_forward.Thruster.ThrustCoeff3 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[2];
            spacecraft_forward.Thruster.ThrustCoeff4 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[3];
            spacecraft_forward.Thruster.ThrustCoeff5 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[4];
            spacecraft_forward.Thruster.MassFlowCoeff1 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[0];
            spacecraft_forward.Thruster.MassFlowCoeff2 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[1];
            spacecraft_forward.Thruster.MassFlowCoeff3 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[2];
            spacecraft_forward.Thruster.MassFlowCoeff4 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[3];
            spacecraft_forward.Thruster.MassFlowCoeff5 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[4];
            spacecraft_backward.Thruster.ThrustModel = "ThrustMassPolynomial";
            spacecraft_backward.Thruster.MinimumUsablePower = this->myjourney->mymission->emtgmission->options.engine_input_power_bounds[0];
            spacecraft_backward.Thruster.MaximumUsablePower = this->myjourney->mymission->emtgmission->options.engine_input_power_bounds[1];
            spacecraft_backward.Thruster.ThrustCoeff1 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[0];
            spacecraft_backward.Thruster.ThrustCoeff2 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[1];
            spacecraft_backward.Thruster.ThrustCoeff3 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[2];
            spacecraft_backward.Thruster.ThrustCoeff4 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[3];
            spacecraft_backward.Thruster.ThrustCoeff5 = this->myjourney->mymission->emtgmission->options.engine_input_thrust_coefficients[4];
            spacecraft_backward.Thruster.MassFlowCoeff1 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[0];
            spacecraft_backward.Thruster.MassFlowCoeff2 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[1];
            spacecraft_backward.Thruster.MassFlowCoeff3 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[2];
            spacecraft_backward.Thruster.MassFlowCoeff4 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[3];
            spacecraft_backward.Thruster.MassFlowCoeff5 = this->myjourney->mymission->emtgmission->options.engine_input_mass_flow_rate_coefficients[4];
        }
        //thruster coefficients from library
        else if (this->myjourney->mymission->emtgmission->options.engine_type > 5) 
        {
            double at, bt, ct, dt, et, ht, gt, af, bf, cf, df, ef, hf, gf, minP, maxP;
            EMTG::Astrodynamics::get_thruster_coefficients_from_library(&(this->myjourney->mymission->emtgmission->options),
                                                                        minP,
                                                                        maxP,
                                                                        at,
                                                                        bt,
                                                                        ct,
                                                                        dt,
                                                                        et,
                                                                        gt,
                                                                        ht,
                                                                        af,
                                                                        bf,
                                                                        cf,
                                                                        df,
                                                                        ef,
                                                                        gf,
                                                                        hf);
            spacecraft_forward.Thruster.ThrustModel = "ThrustMassPolynomial";
            spacecraft_forward.Thruster.MinimumUsablePower = minP;
            spacecraft_forward.Thruster.MaximumUsablePower = maxP;
            spacecraft_forward.Thruster.ThrustCoeff1 = at;
            spacecraft_forward.Thruster.ThrustCoeff2 = bt;
            spacecraft_forward.Thruster.ThrustCoeff3 = ct;
            spacecraft_forward.Thruster.ThrustCoeff4 = dt;
            spacecraft_forward.Thruster.ThrustCoeff5 = et;
            spacecraft_forward.Thruster.MassFlowCoeff1 = at;
            spacecraft_forward.Thruster.MassFlowCoeff2 = bt;
            spacecraft_forward.Thruster.MassFlowCoeff3 = ct;
            spacecraft_forward.Thruster.MassFlowCoeff4 = dt;
            spacecraft_forward.Thruster.MassFlowCoeff5 = et;
            spacecraft_forward.Thruster.ThrustModel = "ThrustMassPolynomial";
            spacecraft_forward.Thruster.MinimumUsablePower = minP;
            spacecraft_forward.Thruster.MaximumUsablePower = maxP;
            spacecraft_forward.Thruster.ThrustCoeff1 = at;
            spacecraft_forward.Thruster.ThrustCoeff2 = bt;
            spacecraft_forward.Thruster.ThrustCoeff3 = ct;
            spacecraft_forward.Thruster.ThrustCoeff4 = dt;
            spacecraft_forward.Thruster.ThrustCoeff5 = et;
            spacecraft_forward.Thruster.MassFlowCoeff1 = at;
            spacecraft_forward.Thruster.MassFlowCoeff2 = bt;
            spacecraft_forward.Thruster.MassFlowCoeff3 = ct;
            spacecraft_forward.Thruster.MassFlowCoeff4 = dt;
            spacecraft_forward.Thruster.MassFlowCoeff5 = et;
        }
	}

	//method
	void set_iBurn() {
		//  if this phase isn't the first nor last, then we do not currently use impulsive maneuvers
		if (!isFirstPhase && !isLastPhase) { return; }
		//  get the departure and arrival types for this phase
		int depart_type = this->myjourney->mymission->emtgmission->options.journey_departure_type[this->myjourney->j];
		int lv_type = this->myjourney->mymission->emtgmission->options.LV_type;
		int arrival_type = this->myjourney->mymission->emtgmission->options.journey_arrival_type[this->myjourney->j];
		double isp_chem = this->myjourney->mymission->emtgmission->options.IspChem;
		double isp_DS = this->myjourney->mymission->emtgmission->options.IspDS;
		 
		if (depart_type == 1) { spacecraft_forward.iBurn.BCisParkingOrbit = true; }
		if (arrival_type == 0) { spacecraft_backward.iBurn.BCisParkingOrbit = true; }

		//  if this is the first journey and first phase
		if ((this->myjourney->j == 0) && isFirstPhase) {
			//  if using launch/direct insertion with the EDS motor
			//+ OR
			//+ if using depart from parking orbit with the EDS motor
			if ((depart_type == 0 && lv_type == -1) || depart_type == 1){
				spacecraft_forward.iBurn.UseImpulsive = true;
				spacecraft_forward.iBurn.IsEDS = true;
				spacecraft_forward.iBurn.CoordinateSystem = spacecraft_forward.CoordinateSystem;
				spacecraft_forward.iBurn.Origin = mybodies[0].name;
				spacecraft_forward.iBurn.Isp = isp_DS;
				spacecraft_forward.iBurn.g = this->myjourney->mymission->emtgmission->options.g0;
				spacecraft_forward.iBurn.c = isp_DS * (spacecraft_forward.iBurn.g / 1000.0);
				spacecraft_forward.iBurn.DecrementMass = false;
				spacecraft_forward.iBurn.TankName = spacecraft_forward.Thruster.Tank.Name;
			}
		}
		//  if not the first journey, but is the first phase and using a direct insertion 
		//+ or departure from a parking orbit
		else if (isFirstPhase && (depart_type == 0 || depart_type == 1)) {
			spacecraft_forward.iBurn.UseImpulsive = true;
			spacecraft_forward.iBurn.CoordinateSystem = spacecraft_forward.CoordinateSystem;
			spacecraft_forward.iBurn.Origin = mybodies[0].name;
			spacecraft_forward.iBurn.Isp = isp_chem;
			spacecraft_forward.iBurn.g = this->myjourney->mymission->emtgmission->options.g0;
			spacecraft_forward.iBurn.c = isp_chem * (spacecraft_forward.iBurn.g / 1000.0);
			spacecraft_forward.iBurn.TankName = spacecraft_forward.Thruster.Tank.Name;
		}

		//  if the last phase and arrival is either insertion into a parking orbit or rendezvous with chemical
		if (isLastPhase && (arrival_type == 0 || arrival_type == 1)) {
			spacecraft_backward.iBurn.UseImpulsive = true;
			spacecraft_backward.iBurn.CoordinateSystem = spacecraft_backward.CoordinateSystem;
			spacecraft_backward.iBurn.Origin = mybodies[1].name;
			spacecraft_backward.iBurn.Isp = isp_chem;
			spacecraft_backward.iBurn.g = this->myjourney->mymission->emtgmission->options.g0;
			spacecraft_backward.iBurn.c = isp_chem * (spacecraft_forward.iBurn.g / 1000.0);
			spacecraft_backward.iBurn.DecrementMass = false;
			spacecraft_backward.iBurn.TankName = spacecraft_backward.Thruster.Tank.Name;
		}

		//  need to set the delta-v components
	
	}

	//method
	void get_flyby_data() {
		//first reset the spacecraft flyby initial conditions
		spacecraft_forward.flyby_distance_lowerbound = 0;
		spacecraft_forward.flyby_distance_upperbound = 1e10;
		spacecraft_forward.flyby_velocity_lowerbound = 0;
		spacecraft_forward.flyby_velocity_upperbound = 1e10;
		spacecraft_backward.flyby_distance_lowerbound = 0;
		spacecraft_backward.flyby_distance_upperbound = 1e10;
		spacecraft_backward.flyby_velocity_lowerbound = 0;
		spacecraft_backward.flyby_velocity_upperbound = 1e10;
		//a middle phase always starts with a flyby
		if (!isFirstPhase) { StartsWithFlyby = true; }
		//a middle phase always ends with a flyby
		if (!isLastPhase)  { EndsWithFlyby = true; }
		//a starting phase may start with a flyby if its journey has a flyby departure type (note: 3 and 4 only valid for successive journeys; therefore a previous journey should exist)
		if (this->myjourney->mymission->emtgmission->options.journey_departure_type[this->myjourney->j] == 3 || this->myjourney->mymission->emtgmission->options.journey_departure_type[this->myjourney->j] == 4 || this->myjourney->mymission->emtgmission->options.journey_departure_type[this->myjourney->j] == 6) {
			if (isFirstPhase) { StartsWithFlyby = true; }
		}
		//an ending phase may end with a flyby if the next journey is NOT the last journey and if it has a flyby departure type (note: 3 and 4 only valid for successive journeys; therefore a previous journey should exist)
		if (!this->myjourney->isLastJourney) {
			if (this->myjourney->mymission->emtgmission->options.journey_departure_type[this->myjourney->j + 1] == 3 || this->myjourney->mymission->emtgmission->options.journey_departure_type[this->myjourney->j + 1] == 4 || this->myjourney->mymission->emtgmission->options.journey_departure_type[this->myjourney->j + 1] == 6) {
				if (isLastPhase) { EndsWithFlyby = true; }
			}
		}
		//if we start with a flyby, then the forward spacecraft will need constraints applicable to a flyby at its body
		if (StartsWithFlyby) {
			//distance bounds
			spacecraft_forward.flyby_distance_lowerbound = mybodies[0].minimum_safe_flyby_altitude + mybodies[0].radius;
			if (mybodies[0].mass < 1.0e25) { spacecraft_forward.flyby_distance_upperbound = 10.0*spacecraft_forward.flyby_distance_lowerbound; }
			else { spacecraft_forward.flyby_distance_upperbound = 300.0*spacecraft_forward.flyby_distance_lowerbound; }
			//velocity bounds
			spacecraft_forward.flyby_velocity_lowerbound = sqrt(2.0 * mybodies[0].mu / (spacecraft_forward.flyby_distance_upperbound));
			spacecraft_forward.flyby_velocity_upperbound = sqrt(2.0 * mybodies[0].mu / (spacecraft_forward.flyby_distance_lowerbound) + 25.0*25.0);
			//based on the vinf_in, vinf_out, flyby_altitude, and body; calculate the spacecraft state at periapsis (calling a method of phase(). this method doesn't depend on phase() members, but still this could be done better to mitigate unforseen future implementation issues).
			if (isFirstPhase) {
				//fill in later when multiple journeys are addressed!!!!!!!!!!!
				//this may already be done [see Jacob's confirmation of my email from 2014_07_30 (i.e. members of phase() pertaining to flyby are with respect to left (in time) boundary point]
			}
			else {
				//puke
				spacecraft_forward.flyby_states = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].calculate_flyby_periapse_state(this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].V_infinity_in,
					this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].V_infinity_out, this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].flyby_altitude, mybodies[0]);
			}
		}
		if (EndsWithFlyby) {
			//distance bounds
			spacecraft_backward.flyby_distance_lowerbound = mybodies[1].minimum_safe_flyby_altitude + mybodies[1].radius;
			if (mybodies[1].mass < 1.0e25) { spacecraft_backward.flyby_distance_upperbound = 10.0*spacecraft_backward.flyby_distance_lowerbound; }
			else { spacecraft_backward.flyby_distance_upperbound = 300.0*spacecraft_backward.flyby_distance_lowerbound; }
			//velocity bounds
			spacecraft_backward.flyby_velocity_lowerbound = sqrt(2.0 * mybodies[1].mu / (spacecraft_backward.flyby_distance_upperbound));
			spacecraft_backward.flyby_velocity_upperbound = sqrt(2.0 * mybodies[1].mu / (spacecraft_backward.flyby_distance_lowerbound) + 25.0*25.0);
			//based on the vinf_in, vinf_out, flyby_altitude, and body; calculate the spacecraft state at periapsis (calling a method of phase(). this method doesn't depend on phase() members, but still this could be done better to mitigate unforseen future implementation issues).
			if (isLastPhase) {
				//fill in later when multiple journeys are addressed!!!!!!!!!!!
			}
			else {
				//puke
				spacecraft_backward.flyby_states = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].calculate_flyby_periapse_state(this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p + 1].V_infinity_in,
					this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p + 1].V_infinity_out, this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p + 1].flyby_altitude, mybodies[1]);
			}
		}
	}//end of method

	//method
	void set_initialconditions() {
		//declarations
		double epoch;
		double body_states[6];

		//get and set forward spacecraft initial conditions relative to its body
		if (StartsWithFlyby) {
			for (int i = 0; i < 6; ++i) {
				spacecraft_forward.initialconditions[i] = spacecraft_forward.flyby_states(i, 0);
			}
		}
		else {
			epoch = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].phase_start_epoch;
			this->mybodies[0].locate_body(epoch, body_states, false, &this->myjourney->mymission->emtgmission->options);
			for (int i = 0; i < 6; ++i) {
				spacecraft_forward.initialconditions[i] = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].state_at_beginning_of_phase[i] - body_states[i];
			}
		}

		//get and set backward spacecraft initial conditions relative to its body
		if (EndsWithFlyby) {
			for (int i = 0; i < 6; ++i) {
				spacecraft_backward.initialconditions[i] = spacecraft_backward.flyby_states(i, 0);
			}
		}
		else {
			epoch = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].phase_end_epoch;
			this->mybodies[1].locate_body(epoch, body_states, false, &this->myjourney->mymission->emtgmission->options);
			for (int i = 0; i < 6; ++i) {
				spacecraft_backward.initialconditions[i] = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].state_at_end_of_phase[i] - body_states[i];
			}
		}
	}

	//method
	void set_epochs() {
		this->iepoch = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].phase_start_epoch;
		this->fepoch = this->myjourney->mymission->emtgmission->journeys[this->myjourney->j].phases[this->p].phase_end_epoch;
	}

	//method
	virtual void get_time_allocation();

	//method
	virtual void append_step(gmatstep agmatstep);

};


//Class() 'gmatstep'
class gmatstep: public EMTG::gmatbaseclass {

public:
	//constructors
	gmatstep();
	gmatstep(gmatphase* aphase, int step) {
		myphase = aphase;
		s = step;
		gs = myphase->mysteps.size();
		id = myphase->id + "gs" + std::to_string(gs);
		this->find_my_spacecraft();
		this->set_names();
		this->isMatchPoint();
		this->usePushBackStill();
		this->set_stepsize();
		this->set_epochs();
		this->set_r_v_wrt_to_body();
		this->set_thrust();
	};

	//destructor
	~gmatstep(){};

	//members
	gmatphase* myphase;
	struct gmat_spacecraft* myspacecraft;
	int s;     //original emtg 's' step number
	int gs;    //gmat step number
	string id; //gmat step id
	struct gmat_propagator propagator; //the propagator to use for my "gmat step"
	vector <double> thrustvector;
	vector <string> thrustvectornames;
	double stepsize;
	double iepoch;
	double fepoch;
	double initial_position_diff[3];
	double initial_velocity_diff[3];
	double final_position_diff[3];
	double final_velocity_diff[3];
	bool inSOIatStart;
	bool inSOIatEnd;
	bool isMatchPointStep;
	bool usePushBack;
	int theInsertIndex;
	bool zeroPropagate = true;
	bool allowTheTimeStep2Vary = false;

	//method
	void find_my_spacecraft() {
		if (s < (this->myphase->myjourney->mymission->emtgmission->options.num_timesteps / 2)) {
			myspacecraft = &this->myphase->spacecraft_forward;
		}
		else {
			myspacecraft = &this->myphase->spacecraft_backward;
		}
	}

	//method
	void set_names() {
		//store my thrust vector name, allowing it to vary during the optimization sequence
		this->setVary("ThrustVector_" + this->id + "_Direction1");
		this->setVary("ThrustVector_" + this->id + "_Direction2");
		this->setVary("ThrustVector_" + this->id + "_Direction3");
		//we will need to ensure that the thrust vector remains a unit vector during optimization, therefore 
		//we introduce a constraint. first calculate the thrust vector magnitude, then constrain it
		this->setCalculate("ThrustUnitVectorMagnitude_" + this->id, "sqrt(( ThrustVector_" + this->id + "_Direction1" + " * 2 - 1) ^ 2 + ( ThrustVector_" + this->id + "_Direction2 * 2 - 1) ^ 2 + ( ThrustVector_" + this->id + "_Direction3 * 2 - 1) ^ 2 )");
		this->setConstraint("ThrustUnitVectorMagnitude_" + this->id, "<=", "1.0");
		//use the 'Equation' command (i.e. "calculate") in GMAT to assign the thrust directions
		this->setCalculate(myspacecraft->Thruster.Name + ".ThrustDirection1", "( ThrustVector_" + this->id + "_Direction1 * 2 - 1 ) / ThrustUnitVectorMagnitude_" + this->id);
		this->setCalculate(myspacecraft->Thruster.Name + ".ThrustDirection2", "( ThrustVector_" + this->id + "_Direction2 * 2 - 1 ) / ThrustUnitVectorMagnitude_" + this->id);
		this->setCalculate(myspacecraft->Thruster.Name + ".ThrustDirection3", "( ThrustVector_" + this->id + "_Direction3 * 2 - 1 ) / ThrustUnitVectorMagnitude_" + this->id);
		//use the 'Equation' command (i.e. "calculate") in GMAT to assign the thrust level
		//this->setCalculate(myspacecraft->Thruster.Name + ".C1", "( ThrusterMaxThrust * ThrustUnitVectorMagnitude_" + this->id + " )");
	}

	//method
	void isMatchPoint() {
		if (s == ((this->myphase->myjourney->mymission->emtgmission->options.num_timesteps / 2) - 1) || s == (this->myphase->myjourney->mymission->emtgmission->options.num_timesteps / 2)) { isMatchPointStep = true; }
		else { isMatchPointStep = false; }
	}

	//method
	void usePushBackStill() {
		//if the first step then use push_back
		if (gs == 0) { usePushBack = true; }
		//if this step is the first to use insert() instead of use push_back()
		else if (this->myphase->mysteps[gs - 1].isMatchPointStep && !this->myphase->mysteps[gs - 1].myspacecraft->isForward) {
			usePushBack = false;
			theInsertIndex = myphase->mysteps.size() - 1;
		}
		//else, continue to use the push_back method
		else { usePushBack = true; }

		//if the switch to using insert() has started, then continue to use insert()
		if (gs > 1) {
			if (!this->myphase->mysteps[gs - 2].usePushBack) {
				usePushBack = false;
				theInsertIndex = this->myphase->mysteps[gs - 2].theInsertIndex;
			}
		}
	}

	//method
	void set_stepsize() {
		//get the current step size, which by definition can be a variable timestep although we always use uniform in practice
		if (s == 0 || s == this->myphase->myjourney->mymission->emtgmission->options.num_timesteps - 1) {
			stepsize = 0.5*this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].time_step_sizes[s];
		}
		else {
			if (this->myspacecraft->isForward) {
				stepsize = 0.5*this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].time_step_sizes[s - 1]
					+ 0.5*this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].time_step_sizes[s];
			}
			else {
				stepsize = 0.5*this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].time_step_sizes[s]
					+ 0.5*this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].time_step_sizes[s + 1];
			}
		}
	}

	//method
	void set_stepsize(double astepsize) {
		//reset stepsize
		if (astepsize < 0.0) { stepsize = -astepsize; }
		else { stepsize = astepsize; }
		//set epochs and r,v w.r.t to body
		this->set_epochs();
		this->set_r_v_wrt_to_body();
	}

	//method
	void set_epochs() {
		//the first step
		if (myphase->mysteps.size() == 0)  {
			iepoch = this->myphase->iepoch;
			fepoch = iepoch + stepsize;
		}
		else {
			//still using push_back(), so we can use the back() method call
			if (usePushBack) {
				iepoch = this->myphase->mysteps.back().fepoch;
				fepoch = iepoch + stepsize;
			}
			//not using push_back()
			else {
				iepoch = this->myphase->mysteps[theInsertIndex].fepoch;
				fepoch = iepoch + stepsize;
			}
		}
	}

	//method
	void set_r_v_wrt_to_body() {
		//declarations
		double body_istates[6];
		double body_fstates[6];
		//get body's states at 'iepoch' and 'fepoch'
		if (myspacecraft->isForward) {
			myphase->mybodies[0].locate_body(iepoch, body_istates, false, &this->myphase->myjourney->mymission->emtgmission->options);
			myphase->mybodies[0].locate_body(fepoch, body_fstates, false, &this->myphase->myjourney->mymission->emtgmission->options);
		}
		else {
			myphase->mybodies[1].locate_body(iepoch, body_istates, false, &this->myphase->myjourney->mymission->emtgmission->options);
			myphase->mybodies[1].locate_body(fepoch, body_fstates, false, &this->myphase->myjourney->mymission->emtgmission->options);
		}
		//there is an odd case where we could analyze the 'body_istates' of the first body and 'body_fstates' of the
		//second body across the matchpoint, but this appears unnecessary at this time.
		//compute the state relative difference of the spacecraft and body at 'iepoch'
		if (myspacecraft->isForward) {
			if (s == 0) {
				for (int index = 0; index < 3; ++index) {
					initial_position_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].state_at_beginning_of_phase[index] - body_istates[index];
					initial_velocity_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].state_at_beginning_of_phase[index + 3] - body_istates[index + 3];
				}
			}
			else {
				for (int index = 0; index < 3; ++index) {
					initial_position_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s - 1][index] - body_istates[index];
					initial_velocity_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s - 1][index + 3] - body_istates[index + 3];
				}
			}
		}
		else {
			for (int index = 0; index < 3; ++index) {
				initial_position_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s][index] - body_istates[index];
				initial_velocity_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s][index + 3] - body_istates[index + 3];
			}
		}
		//compute the state relative difference of the spacecraft and body at 'fepoch'
		if (myspacecraft->isForward) {
			for (int index = 0; index < 3; ++index) {
				final_position_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s][index] - body_fstates[index];
				final_velocity_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s][index + 3] - body_fstates[index + 3];
			}
		}
		else {
			if (s == (this->myphase->myjourney->mymission->emtgmission->options.num_timesteps - 1)) {
				for (int index = 0; index < 3; ++index) {
					final_position_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].state_at_end_of_phase[index] - body_fstates[index];
					final_velocity_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].state_at_end_of_phase[index + 3] - body_fstates[index + 3];
				}
			}
			else {
				for (int index = 0; index < 3; ++index) {
					final_position_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s + 1][index] - body_fstates[index];
					final_velocity_diff[index] = this->myphase->myjourney->mymission->emtgmission->journeys[this->myphase->myjourney->j].phases[this->myphase->p].spacecraft_state[s + 1][index + 3] - body_fstates[index + 3];
				}
			}
		}
		//evaluate whether the spacecraft was in the body's SOI at 'iepoch' and 'fepoch'
		if (myspacecraft->isForward) {
			if (math::norm(initial_position_diff, 3) < myphase->mybodies[0].r_SOI) { inSOIatStart = true; }
			else { inSOIatStart = false; }
			if (math::norm(final_position_diff, 3) < myphase->mybodies[0].r_SOI) { inSOIatEnd = true; }
			else { inSOIatEnd = false; }
		}
		else {
			if (math::norm(initial_position_diff, 3) < myphase->mybodies[1].r_SOI) { inSOIatStart = true; }
			else { inSOIatStart = false; }
			if (math::norm(final_position_diff, 3) < myphase->mybodies[1].r_SOI) { inSOIatEnd = true; }
			else { inSOIatEnd = false; }
		}

	};

	//method
	void set_thrust() {
		//declarations
		int j = myphase->myjourney->j;
		int p = myphase->p;
		double lower_bound = -1.0;
		double upper_bound = 1.0;
		stringstream tempstream;
		tempstream.precision(25);
		//store control history
		for (int index = 0; index < 3; ++index) {
			thrustvector.push_back((this->myphase->myjourney->mymission->emtgmission->journeys[j].phases[p].control[s][index] - lower_bound) / (upper_bound - lower_bound));
			tempstream << thrustvector.back();
			setVariable("ThrustVector_" + this->id + "_Direction" + std::to_string(index + 1), tempstream.str());
			tempstream.str("");
		}
		setVariable("ThrustUnitVectorMagnitude_" + this->id);
	}

	//method
	void setFMandProp(bool amIaCloseApproach) {
		// ***********
		// Future Work: Include EMTG 3rd Body Perturbations in PointMasses member of Propagator
		//				Also, should techinically check to make sure the the central body is not also duplicated as a pointmass (e.g. a trajectory from a moon to an orbit about its central body)
		// ***********

		//declarations
		int body_index;
		bool isAFlyby;

		//select correct body index and whether we are performing a flyby
		if (myspacecraft->isForward) {
			body_index = 0; 
			if (myphase->StartsWithFlyby && amIaCloseApproach) { isAFlyby = true; }
			else { isAFlyby = false; }
		}
		else {
			body_index = 1; 
			if (myphase->EndsWithFlyby && amIaCloseApproach) { isAFlyby = true; }
			else { isAFlyby = false; }
		}

		//CloseApproach Flag
		propagator.isCloseApproach = amIaCloseApproach;
		//ForceModel Name
		if (amIaCloseApproach) {
			if (isAFlyby) { propagator.ForceModel.Name = "FM_" + myphase->mybodies[body_index].name; }
			else { propagator.ForceModel.Name = "FM_" + myphase->mybodies[body_index].central_body_name; }
		}
		else { 
			//use both bodies as 3rd body perturbation if they are different, otherwise we only need to name it once
			if (myphase->mybodies[0].name.compare(myphase->mybodies[1].name) == 0) { propagator.ForceModel.Name = "FM_" + myphase->mybodies[body_index].central_body_name + "_3rdBodies_" + myphase->mybodies[0].name; }
			else { propagator.ForceModel.Name = "FM_" + myphase->mybodies[body_index].central_body_name + "_3rdBodies_" + myphase->mybodies[0].name + "_" + myphase->mybodies[1].name; }
		}
		//Propagator Name
		if (amIaCloseApproach) { propagator.Name = "Propagator_" + propagator.ForceModel.Name + "_CloseApproach"; }
		else { propagator.Name = "Propagator_" + propagator.ForceModel.Name; }
		//ForceModel Central Body
		if (isAFlyby) { propagator.ForceModel.CentralBody = this->myphase->mybodies[body_index].name; }
		else { propagator.ForceModel.CentralBody = this->myphase->mybodies[body_index].central_body_name; }
		//ForeceModel PointMasses
		propagator.ForceModel.PointMasses.clear();
		propagator.ForceModel.PointMasses.push_back(propagator.ForceModel.CentralBody);
		if (amIaCloseApproach) {
			if (isAFlyby) { propagator.ForceModel.PointMasses.push_back(myphase->mybodies[body_index].central_body_name); }
		}
		else {
			//use both bodies as 3rd body perturbation if they are different, otherwise we only need to name it once
			if (myphase->mybodies[0].name.compare(myphase->mybodies[1].name) == 0) { propagator.ForceModel.PointMasses.push_back(this->myphase->mybodies[0].name); }
			else {
				propagator.ForceModel.PointMasses.push_back(this->myphase->mybodies[0].name);
				propagator.ForceModel.PointMasses.push_back(this->myphase->mybodies[1].name);
			}
		}
	}

	//method
	void reset() {
		//clear the names. if reusing the 'gmatstep' object but with a new id this is necessary
		thrustvectornames.clear();
		vary.clear();
		calculate.clear();
		constraints.clear();
		variables.clear();
		allowTheTimeStep2Vary = false;
		//update and reset
		gs++;
		id = myphase->id + "gs" + std::to_string(gs);
		this->set_names();
		this->isMatchPoint();
		this->usePushBackStill();
		this->set_stepsize();
		this->set_epochs();
		this->set_r_v_wrt_to_body();
		this->set_thrust();
	}

	//method
	void scale_stepsize(double scalar) {
		stepsize *= scalar;
		this->set_epochs();
		this->set_r_v_wrt_to_body();
	}

	//method
	void printManeuver(std::ofstream& File) {
		if (this->myspacecraft->iBurn.UseImpulsive && (this->gs == 0 || this->myphase->mysteps.size() - 1 == this->gs)) {
			File << "   " << "Maneuver '" << this->myspacecraft->iBurn.Type << " " << this->myspacecraft->iBurn.Name << "' " << this->myspacecraft->iBurn.Name << "( " << this->myspacecraft->Name << " )" << std::endl;
		}
	}
	
	//method
	void printBeginBurn(std::ofstream& File) {
		File << "   " << "Begin" << this->myspacecraft->fBurn.Type << " 'Begin" << this->myspacecraft->fBurn.Type << " " << this->myspacecraft->fBurn.Name << "' " << this->myspacecraft->fBurn.Name << "( " << this->myspacecraft->Name << " )" << std::endl;
	}

	//method
	void printEndBurn(std::ofstream& File) {
		File << "   " << "End" << this->myspacecraft->fBurn.Type << " 'End" << this->myspacecraft->fBurn.Type << " " << this->myspacecraft->fBurn.Name << "' " << this->myspacecraft->fBurn.Name << "( " << this->myspacecraft->Name << " )" << std::endl;
	}

};


//Class() 'gmatscripter'
class gmatscripter {

public:
	//constructors
	gmatscripter();
	gmatscripter(mission* mission_in);

	//destructor
	virtual ~gmatscripter();

	//members
	mission*  ptr_gmatmission;
	//missionoptions gmatoptions;
	std::ofstream GMATfile;
	//a temporary file for debugging purposes
	std::ofstream GMATDebug;
	//gmatclass objects
	gmatmission GMATMission;
	//TAIModJOffset = 2400000.5 - 2430000;
	double TAIModJOffset = -29999.5;

	//methods
	//create file
	virtual void create_GMAT_file();
	virtual void create_GMAT_missions();
	virtual void create_GMAT_journeys();
	virtual void create_GMAT_phases();
	virtual void create_GMAT_steps();
	//post creation pass for additional variable, vary, calculate, and constraint statements
	virtual void postpass_GMAT_phases();
	virtual void postpass_GMAT_journeys();
	virtual void postpass_GMAT_missions();
	//model setup
	virtual void write_GMAT_preamble();
	virtual void write_GMAT_spacecraft();
	virtual void write_GMAT_hardware();
	virtual void write_GMAT_nonstandardbody();
	virtual void write_GMAT_forcemodels();
	virtual void write_GMAT_propagators();
	virtual void write_GMAT_burns();
	virtual void write_GMAT_coordinatesystems();
	virtual void write_GMAT_solvers();
	virtual void write_GMAT_subscribers();
	virtual void write_GMAT_variables();
	//optimization methods
	virtual void write_GMAT_beginmissionsequence();
	virtual void write_GMAT_initialconditions();
	virtual void write_GMAT_preoptimization_calculations();
	virtual void write_GMAT_optimization();
	//write calls
	virtual void create_GMAT_spacecraft(struct gmat_spacecraft& spacecraft);
    virtual void create_GMAT_powersystem(struct gmat_spacecraft& spacecraft);
	virtual void create_GMAT_fueltank(struct gmat_spacecraft& spacecraft);
	virtual void create_GMAT_thruster(struct gmat_spacecraft& spacecraft);
	virtual void create_GMAT_forcemodel(struct gmat_forcemodel& forcemodel);
	virtual void create_GMAT_propagator(struct gmat_propagator& propagator);
	virtual void create_GMAT_burn(struct gmat_fburn& burn);
	virtual void create_GMAT_burn(struct gmat_iburn& burn);
	virtual void create_GMAT_coordinatesystem(string bodyname);
	virtual void rewrite_GMAT_initialconditions(class gmatphase& agmatphase);
	virtual void create_GMAT_initialconditions(struct gmat_spacecraft& spacecraft, string specialchar = "   ", bool WithCoordinateSyst = true,
											   bool PrintPosition = true, bool PrintVelocity = true);

	//GMAT Command Methods
	virtual void aux_GMAT_propagate(class gmatstep& agmatstep, bool useZeroPropagate);
	virtual void PenUp();
	virtual void PenDown();

	//reports
	virtual void write_GMAT_report(class gmatstep& agmatstep, bool isbeforemaneuver, bool writecontrolhistory);

	//writeout the GMAT script
	virtual void write_GMAT_script();

}; // end of class gmatscript


}  // end of EMTG namespace

#endif // end of GMATSCRIPTER_H_

