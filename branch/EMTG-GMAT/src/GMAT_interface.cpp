//EMTG GMAT interface
//collected code written by Max Schadegg
//also contains phase-specific code written by Jacob Englander

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

#include "mission.h"
#include "Astrodynamics.h"
#include "body.h"
#include "MGADSMphase.h"

#include "SpiceUsr.h"

#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"




using namespace std;

namespace EMTG
{
	//function to output .script for GMAT simulation (header preamble)
	void mission::output_GMAT_preamble()
	{
		//define body class objects
		vector <EMTG::Astrodynamics::body> missionbodies_unique;
		vector <EMTG::Astrodynamics::body> missionbodies;
		//create a vector of strings for storing spacecraft names
		vector <string> SC_created;
		vector <fs::path> SPICE_files;
		SPICEDOUBLE_CELL (spice_coverage, 10000);
		SpiceInt number_of_windows = 0;
		int index_SC = 0;
		int index_body_visited = 0;

		//create a filename, instantiate an ofstream object called GMATfile
		string filename = options.working_directory + "//" + options.mission_name + "_" + options.description + "_GMAT.script";
		string filestring;
		ofstream GMATfile(filename.c_str(), ios::trunc);
		//set floating point decimal precision
		GMATfile.precision(25);
	
		//get the current timestamp from boost and assign it to now
		boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
		std::stringstream timestream;
		timestream << static_cast<int>(now.date().month()) << "/" << now.date().day() << "/" << now.date().year() << " " << now.time_of_day().hours() << ":" << now.time_of_day().minutes() << ":" << now.time_of_day().seconds();

		//create GMAT header
		GMATfile << "%-------------------------------------------------------------------" << endl;
		GMATfile << "%Automated GMAT script created by EMTGv8" << endl;
		GMATfile << "%EMTG options file: " << options.working_directory + "/" + options.mission_name + ".emtgopt" << endl;
		GMATfile << "%EMTG output file: " << options.outputfile << endl;
		GMATfile << "%EMTG output written on: " << timestream.str() << endl;
		GMATfile << "%Script Author: Max Schadegg" << endl;
		GMATfile << "%Created on: 6/3/2013" << endl;
		GMATfile << "%Revised by: Jacob Englander" << endl;
		GMATfile << "%Last revision: 8/9/2013" << endl;
		GMATfile << "%-------------------------------------------------------------------" << endl;
		GMATfile << endl;
		GMATfile << endl;

		//create spacecraft models
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Spacecraft" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		//add bodies of first journey to the mission body list
		//TODO:: add other arrival/departure types
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			if (j == 0)
			{
				for (int p = 0; p < options.number_of_phases[0] + 1; ++p)
					missionbodies.push_back(TheUniverse[0].bodies[options.sequence[0][p] - 1]);
			}
			//add all other bodies to mission body list
			else
			{
				for (int p = 1; p < options.number_of_phases[j] + 1; ++p)
					missionbodies.push_back(TheUniverse[j].bodies[options.sequence[j][p] - 1]);
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

		//create and write out two s/c for each phase in each journey (forward + backward)
		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				//call method output_GMAT_spacecraft of phases class
				journeys[j].phases[p].output_GMAT_spacecraft(j, p, SC_created, index_SC, missionbodies, index_body_visited, GMATfile);
			}
		}
	
		//create hardware components for s/c
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Hardware components" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		index_body_visited = 0;
		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				journeys[j].phases[p].output_GMAT_fueltank_and_thruster(j, p, missionbodies, index_body_visited, GMATfile);
			}
		}
		GMATfile << endl;
	
		//create Force Models 
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
			if ((missionbodies[index_body_visited].name != "Sun") && (missionbodies[index_body_visited].name != "Mercury") && (missionbodies[index_body_visited].name != "Venus") && (missionbodies[index_body_visited].name != "Earth") && (missionbodies[index_body_visited].name != "Mars") && (missionbodies[index_body_visited].name != "Jupiter") && (missionbodies[index_body_visited].name != "Saturn") && (missionbodies[index_body_visited].name != "Uranus") && (missionbodies[index_body_visited].name != "Neptune") && (missionbodies[index_body_visited].name != "Pluto") )
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
				EMTG::filesystem::get_all_files_with_extension(fs::path(options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

				for (size_t k = 0; k < SPICE_files.size(); ++k)
				{
					filestring = options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

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
			if ((missionbodies_unique[index_body_visited].central_body_name != "Sun") && (missionbodies_unique[index_body_visited].central_body_name != "Mercury") && (missionbodies_unique[index_body_visited].central_body_name != "Venus") && (missionbodies_unique[index_body_visited].central_body_name != "Earth") && (missionbodies_unique[index_body_visited].central_body_name != "Mars") && (missionbodies_unique[index_body_visited].central_body_name != "Jupiter") && (missionbodies_unique[index_body_visited].central_body_name != "Saturn") && (missionbodies_unique[index_body_visited].central_body_name != "Uranus") && (missionbodies_unique[index_body_visited].central_body_name != "Neptune") && (missionbodies_unique[index_body_visited].central_body_name != "Pluto") )
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
				GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisRAConstant = " << TheUniverse[0].LocalFrame.alpha0 << endl;
				GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisRARate = " << TheUniverse[0].LocalFrame.alphadot << endl;
				GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisDECConstant = " << TheUniverse[0].LocalFrame.delta0 << endl;
				GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".SpinAxisDECRate = " << TheUniverse[0].LocalFrame.deltadot << endl;
				GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".RotationConstant = " << TheUniverse[0].LocalFrame.W << endl;
				GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".RotationRate = " << TheUniverse[0].LocalFrame.Wdot << endl;
				GMATfile << missionbodies_unique[index_body_visited].central_body_name << ".OrbitSpiceKernelName = {";
			
				//find which spice file the body is located
				EMTG::filesystem::get_all_files_with_extension(fs::path(options.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files);

				for (size_t k = 0; k < SPICE_files.size(); ++k)
				{
					filestring = options.universe_folder + "/ephemeris_files/" + SPICE_files[k].string();

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
	
		GMATfile << endl;
		GMATfile << endl;
		
		//create Propagators 
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
	
		GMATfile << endl;
		GMATfile << endl;

		//create (finite or impulsive) burn object to be redefined at each time step
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Burns" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				journeys[j].phases[p].output_GMAT_burn_objects(j, p, GMATfile);
			}
		}
		GMATfile << endl;
		GMATfile << endl;

		//create all arrays and variables
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Arrays, Variables, Strings" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		//create mission variables
		GMATfile << "Create Variable ObjectiveFunction FinalEpoch LaunchEpoch_Scaled" << endl;
		//DEBUG: Check that these are the only impulsive-thrust phase types
		if (options.mission_type < 2 || options.mission_type == 5) //impulsive-thrust phase types
		{
			GMATfile << "Create Variable FinalMass_Scaled ThrusterISP" << endl;
		}
		else //low-thrust phase types
		{
			GMATfile << "Create Variable FinalMass_Scaled ThrusterISP ThrusterDutyCycle ThrusterMaxThrust" << endl;
		}
		

		//create variable for each s/c create rdotv temp variable
		for (int index_SC = 0; index_SC < SC_created.size(); ++index_SC)
		{
			GMATfile << "Create Variable " << SC_created[index_SC] << "_RdotV_Scaled " << SC_created[index_SC] << "_PeriapseRadius_Scaled" << endl;
		}
		GMATfile << endl;

		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{			
				
				//create the time and state variables
				journeys[j].phases[p].output_GMAT_create_state_and_time_variables(j, p, GMATfile);

				//create the control variables
				journeys[j].phases[p].output_GMAT_create_interphase_control_variables(j, p, options, GMATfile);

				GMATfile << endl;
			}
		}
		GMATfile << endl;
		GMATfile << endl;

		//will always perform calculations in J2000 eq systems
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Coordinate systems" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		GMATfile << "%Create coordinate systems for plotting/viewing" << endl;
		//TODO:: as of now it is assumed that the first central body is the principle central body
		GMATfile << "Create CoordinateSystem " <<  missionbodies_unique[0].central_body_name << "J2000Eq;" << endl;
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
		GMATfile << endl;
		GMATfile << endl;

		//create solvers (can choose fmincon or VF13ad)
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Solvers" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

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

		//create subscribers (ie. orbit views, plots, etc)
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Subscribers" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		//add central universe body view
		GMATfile << "%Create subscriber for central body view" << endl;
		GMATfile << "Create OrbitView " << missionbodies_unique[0].central_body_name << "View" << endl;
		GMATfile << missionbodies_unique[0].central_body_name << "View.ShowPlot =		true" << endl;
		GMATfile << missionbodies_unique[0].central_body_name << "View.SolverIterations =	 All" << endl;
		GMATfile << missionbodies_unique[0].central_body_name << "View.RelativeZOrder =	501" << endl;
	
		//add which bodies and s/c to plot 
		GMATfile << missionbodies_unique[0].central_body_name << "View.Add =	{";
		for (int index_SC = 0; index_SC < SC_created.size(); ++index_SC)
		{
			GMATfile << SC_created[index_SC] << ", ";
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
		GMATfile << missionbodies_unique[0].central_body_name << "View.CoordinateSystem =		"<< missionbodies_unique[0].central_body_name << "J2000Eq" << endl;
		GMATfile << missionbodies_unique[0].central_body_name << "View.DrawObject = [";
	
		//create flag parameters for plotting
		for (int index_plot = 0; index_plot < (missionbodies_unique.size() + SC_created.size()); ++index_plot)
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
			for (int index_SC = 0; index_SC < SC_created.size(); ++index_SC)
			{
				GMATfile << SC_created[index_SC] << ", ";
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
			GMATfile << missionbodies_unique[index_body_visited].name << "View.CoordinateSystem       = "<< missionbodies_unique[index_body_visited].name << "J2000Eq" << endl;
			GMATfile << missionbodies_unique[index_body_visited].name << "View.DrawObject             = [";
		
			//create flag parameters for plotting
			for (int index_plot = 0; index_plot < (missionbodies_unique.size() + SC_created.size()); ++index_plot)
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
	}


	//function to output .script for GMAT simulation (mission sequence)
	void mission::output_GMAT_mission()
	{
		//declare variables
		string filename = options.working_directory + "//" + options.mission_name + "_" + options.description + "_GMAT.script";
		ofstream GMATfile(filename.c_str(), ios::app);
		GMATfile.precision(25);
		vector <EMTG::Astrodynamics::body> missionbodies;
		vector <string> SC_created;
		std::ostringstream SC_name_temp;
		stringstream prefixstream;
		string prefix;
		math::Matrix<double> Vinf_in(3, 1);
		math::Matrix<double> Vinf_out(3, 1);
		math::Matrix<double> periapse_state_vector(6, 1);
		math::Matrix<double> periapse_position_vector(3, 1);
		math::Matrix<double> periapse_velocity_vector(3, 1);
		math::Matrix<double> angular_momentum_vector(3, 1);
		double LaunchDate_lowerbounds;
		double LaunchDate_upperbounds;
		double Distance_lowerbounds_Forward;
		double Distance_upperbounds_Forward;
		double Distance_lowerbounds_Backward;
		double Distance_upperbounds_Backward;
		double Velocity_lowerbounds_Forward;
		double Velocity_upperbounds_Forward;
		double Velocity_lowerbounds_Backward;
		double Velocity_upperbounds_Backward;
		double ThrustUnitVector_lowerbounds;
		double ThrustUnitVector_upperbounds;
		double CumulatedTOF = 0;
		double boundary_state[6];
		double periapse_velocity_magnitude;
		double delta_t;
		phase* PreviousPhase;
		phase* CurrentPhase;
		phase* NextPhase;
		int index_SC = 0;
		int index_body_visited = 0;
		int index_delta_t;
		int TotalNumberTimeSteps = 0;

		//begin mission sequence
		GMATfile << "BeginMissionSequence" <<endl;
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
		if (options.mission_type < 2 || options.mission_type == 5) //impulsive-thrust phase types
		{
			GMATfile << "	ThrusterISP = " << options.IspChem << endl;

			GMATfile << endl;
	
			//for each journey
			for (int j = 0; j < options.number_of_journeys; ++j)
			{
				//for each phase
				for (int p = 0; p < journeys[j].number_of_phases; ++p)
				{
					GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".K1 = ThrusterISP;" << endl;
					GMATfile << endl;
				}
			}
		}
		else //for all low-thrust phase types
		{
			switch (options.engine_type)
			{
			case 0: //fixed thrust and ISP
				GMATfile << "	ThrusterISP = " << options.IspLT << endl;
				GMATfile << "	ThrusterMaxThrust = " << options.Thrust << endl;
				GMATfile << "	ThrusterDutyCycle = " << options.engine_duty_cycle << endl;
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
			for (int j = 0; j < options.number_of_journeys; ++j)
			{
				//for each phase
				for (int p = 0; p < journeys[j].number_of_phases; ++p)
				{
					GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.DutyCycle = ThrusterDutyCycle;" << endl;
					GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1 = ThrusterMaxThrust;" << endl;
					GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.K1 = ThrusterISP;" << endl;
					GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.DutyCycle = ThrusterDutyCycle;" << endl;
					GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.C1 = ThrusterMaxThrust;" << endl;
					GMATfile << "	Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.K1 = ThrusterISP;" << endl;
					GMATfile << endl;
				}
			}
		} //end code for low-thrust models


		//add bodies of first journey to the mission body list
		//TODO:: add other arrival/departure types
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			if (j == 0)
			{
				for (int p = 0; p < options.number_of_phases[0] + 1; ++p)
					missionbodies.push_back(TheUniverse[0].bodies[options.sequence[0][p] - 1]);
			}
			//add all other bodies to mission body list
			else 
			{
				for (int p = 1; p < options.number_of_phases[j] + 1; ++p)
					missionbodies.push_back(TheUniverse[j].bodies[options.sequence[j][p] - 1]);
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

		//EMTGMAT needs to reference the upper and lower bounds vectors, so we must recalculate them
		calcbounds();

		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				prefix.clear();
				prefixstream << "j" << j << "p" << p << ": ";
				prefix = prefixstream.str();
			
				//define bounds for scaling variables
				
				//flight time bounds
				//first, grab the flight time bounds from the phase

				
				for (int iX = 0; iX < Xdescriptions.size(); ++iX)
				{
					//find index in Xdescriptions where the launch time bounds are located
					if (Xdescriptions[iX] == prefix + "launch epoch (MJD)")
					{
						LaunchDate_lowerbounds = Xlowerbounds[iX];
						LaunchDate_upperbounds = Xupperbounds[iX];
					}
					//find index in Xdescriptions where the flyby altitude bounds for that phase are located
					if (Xdescriptions[iX] == prefix + "flyby altitude constraint (above minimum altitude but below [100x/300x] altitude for [rocky/gas] planets")
					{
						Distance_lowerbounds_Forward =  Xlowerbounds[iX] * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_upperbounds_Forward = -Xlowerbounds[iX] * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_lowerbounds_Backward =  Xlowerbounds[iX] * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
						Distance_upperbounds_Backward = -Xlowerbounds[iX] * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
					}
					//find index in Xdescriptions where the flyby velocity bounds for that phase are located
					if (Xdescriptions[iX] == prefix + "initial velocity increment x")
					{
						Velocity_lowerbounds_Forward = -sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + Xlowerbounds[iX] * Xlowerbounds[iX]);
						Velocity_upperbounds_Forward =  sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + Xupperbounds[iX] * Xupperbounds[iX]);
						Velocity_lowerbounds_Backward = -sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + Xlowerbounds[iX] * Xlowerbounds[iX]);
						Velocity_upperbounds_Backward =  sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + Xupperbounds[iX] * Xupperbounds[iX]);
					}			
				}
				if ((p == 0)||(p == journeys[j].number_of_phases - 1))
				{
					//TODO:: hardcoded in for journey arrivals and departures :-/
					if (missionbodies[index_body_visited].mass < 1.0e25)
					{
						Distance_lowerbounds_Forward = -10 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_upperbounds_Forward =  10 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
					}
					else
					{
						Distance_lowerbounds_Forward = -300 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_upperbounds_Forward =  300 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
					}
					if (missionbodies[index_body_visited + 1].mass < 1.0e25)
					{
						Distance_lowerbounds_Backward = -10 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
						Distance_upperbounds_Backward =  10 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
					}
					else
					{
						Distance_lowerbounds_Backward = -300 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
						Distance_upperbounds_Backward =  300 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
					}		
					Velocity_lowerbounds_Forward = -sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + (-25) * (-25));
					Velocity_upperbounds_Forward =  sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + 25 * 25);
					Velocity_lowerbounds_Backward = -sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + (-25) * (-25));
					Velocity_upperbounds_Backward =  sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + 25 * 25);
				}

				//FORWARD PROPAGATED SPACECRAFT
				//first create string for name of SC, add it to vector of all SC names
				SC_name_temp.str("");
				SC_name_temp << "SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward";
				SC_created.push_back(SC_name_temp.str());
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Forward s/c initial conditions" << endl;

				//interpret state from beginning of each journey based on departure type (//TODO::)
				if (p == 0)
				{
					//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
					missionbodies[index_body_visited].locate_body(journeys[j].phases[p].phase_start_epoch, boundary_state, false, &options);

					switch (options.journey_departure_type[j]) 
					{
					case 0: //launch or direct insertion

						//calculate v_infinity vectors
						for (int k = 0; k < 3; ++k)
						{
							Vinf_out(k) = journeys[j].phases[p].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
						}
						//calculate inc from vinfinity, then make guess at min altitude at planet

						//TODO::
						if (options.journey_departure_elements_type[j] == 200) //if given in orbital elements
						{
							//I dont think John's code works correctly
							//periapse_state_vector = journeys[j].phases[p].calculate_periapse_state_from_asymptote_and_parking_orbit(Vinf_out, options.journey_departure_elements[j][2], options.journey_departure_elements[j][0] - missionbodies[index_body_visited].radius, journeys[j].phases[p].phase_start_epoch, &TheUniverse[options.number_of_journeys - 1], index_body_visited + 1);
							for (int k = 0; k < 3; ++k)
							{
								periapse_position_vector(k) = periapse_state_vector(k);
								periapse_velocity_vector(k) = periapse_state_vector(k + 3);
							}						
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.X = " << periapse_position_vector(0) << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y = " << periapse_position_vector(1) << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z = " << periapse_position_vector(2) << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX = " << periapse_velocity_vector(0) << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY = " << periapse_velocity_vector(1) << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ = " << periapse_velocity_vector(2) << ";" << endl;
						
						}
						else 
						{
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.X = " << journeys[j].phases[p].state_at_beginning_of_phase[0] - boundary_state[0] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y = " << journeys[j].phases[p].state_at_beginning_of_phase[1] - boundary_state[1] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z = " << journeys[j].phases[p].state_at_beginning_of_phase[2] - boundary_state[2] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX = " << journeys[j].phases[p].state_at_beginning_of_phase[3] - boundary_state[3] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY = " << journeys[j].phases[p].state_at_beginning_of_phase[4] - boundary_state[4] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ = " << journeys[j].phases[p].state_at_beginning_of_phase[5] - boundary_state[5] << ";" << endl;
						}
						break;
					case 1: //depart from parking orbit

						break;
					case 2: //free direct departure
					
						break;
					case 3: //depart from flyby
						break;
					}
				}

				//for all other phases, interpret state from beginning of each phase and then add the minimum flyby altitude
				else
				{
					//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
					missionbodies[index_body_visited].locate_body(journeys[j].phases[p].phase_start_epoch, boundary_state, false, &options);

					//calculate v_infinity vectors
					for (int k = 0; k < 3; ++k)
					{
						Vinf_in(k) = journeys[j].phases[p - 1].state_at_end_of_phase[k + 3] - boundary_state[k + 3];
						Vinf_out(k) = journeys[j].phases[p].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
					}

					periapse_state_vector = journeys[j].phases[p].calculate_flyby_periapse_state(Vinf_in, Vinf_out, journeys[j].phases[p - 1].flyby_altitude, missionbodies[index_body_visited]);
			
					//add this position vector to state's initial guess
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.X = " << periapse_state_vector(0) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y = " << periapse_state_vector(1) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z = " << periapse_state_vector(2) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX = " << periapse_state_vector(3) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY = " << periapse_state_vector(4) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ = " << periapse_state_vector(5) << ";" << endl;
				}
				GMATfile << "	" << SC_created[index_SC] << ".CoordinateSystem = " << missionbodies[index_body_visited].name << "J2000Eq;" << endl;			
				GMATfile << endl;
				++index_SC;

	
				//BACKWARD PROPAGATED SPACECRAFT
				SC_name_temp.str("");
				SC_name_temp << "SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward";
				SC_created.push_back(SC_name_temp.str());
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Backward s/c initial conditions" << endl;

				//interpret state from end of each journey based on arrival type (//TODO::)
				if (p == journeys[j].number_of_phases - 1)
				{
					//first we must figure out where the initial position at each phase is, since EMTG goes through center of the body
					missionbodies[index_body_visited + 1].locate_body(journeys[j].phases[p].phase_end_epoch, boundary_state, false, &options);
				
					switch (options.journey_arrival_type[j]) 
					{
					case 0: //parking orbit insertion
							break;
					case 1: //rendezvous (chemical)
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X = " << journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y = " << journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z = " << journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX = " << journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY = " << journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ = " << journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5]  << ";" << endl;
						break;
					case 2: //flyby with bounded VHP
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X = " << journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y = " << journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z = " << journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX = " << journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY = " << journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ = " << journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5]  << ";" << endl;
						break;
					case 3: //rendezvous (LT)
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X = " << journeys[j].phases[p].state_at_end_of_phase[0] - boundary_state[0] << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y = " << journeys[j].phases[p].state_at_end_of_phase[1] - boundary_state[1]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z = " << journeys[j].phases[p].state_at_end_of_phase[2] - boundary_state[2]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX = " << journeys[j].phases[p].state_at_end_of_phase[3] - boundary_state[3]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY = " << journeys[j].phases[p].state_at_end_of_phase[4] - boundary_state[4]  << ";" << endl;
							GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ = " << journeys[j].phases[p].state_at_end_of_phase[5] - boundary_state[5]  << ";" << endl;
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
					missionbodies[index_body_visited + 1].locate_body(journeys[j].phases[p + 1].phase_start_epoch, boundary_state, false, &options);
				
					//calculate v_infinity vectors
					for (int k = 0; k < 3; ++k)
					{
						Vinf_in(k) = journeys[j].phases[p].state_at_end_of_phase[k + 3] - boundary_state[k + 3];
						Vinf_out(k) = journeys[j].phases[p + 1].state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
					}

					periapse_state_vector = journeys[j].phases[p].calculate_flyby_periapse_state(Vinf_in, Vinf_out, journeys[j].phases[p].flyby_altitude, missionbodies[index_body_visited + 1]);

					//add this position vector to state's initial guess
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X = " << periapse_state_vector(0) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y = " << periapse_state_vector(1) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z = " << periapse_state_vector(2) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX = " << periapse_state_vector(3) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY = " << periapse_state_vector(4) << ";" << endl;
					GMATfile << "	" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ = " << periapse_state_vector(5) << ";" << endl;
				}
				GMATfile << "	" << SC_created[index_SC] << ".CoordinateSystem = " << missionbodies[index_body_visited + 1].name << "J2000Eq;" << endl;			
				GMATfile << endl;
			
				//insert scaled launch epoch and wait times for each subsequent journey
				if ((j > 0) && (p == 0))
				{
					GMATfile <<"	%Guess for scaled journey wait times" << endl;
					GMATfile <<"	Journey" << j + 1 <<  "_WaitTime_Scaled = " << ((journeys[j].phases[p].phase_start_epoch - journeys[j - 1].phases[journeys[j - 1].number_of_phases - 1].phase_end_epoch) - options.journey_wait_time_bounds[j][0]) / (options.journey_wait_time_bounds[j][1] - options.journey_wait_time_bounds[j][0]) << ";" << endl;
				}
				if (j == 0)
				{
					if (p == 0)
					{
						GMATfile <<"	%Guess for scaled spacecraft launch epoch" << endl;
						GMATfile <<"	LaunchEpoch_Scaled = " << (journeys[j].phases[p].phase_start_epoch - LaunchDate_lowerbounds) / (LaunchDate_upperbounds - LaunchDate_lowerbounds) << ";" << endl;
					}
				}
				GMATfile << endl;

				//insert scaled states for forward and backward s/c (doing the math in GMAT bc of all the if statements above!)
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Forward s/c scaled initial conditions" << endl;
				GMATfile <<"	" << SC_created[index_SC - 1] << "_X_Scaled = (" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.X - " << Distance_lowerbounds_Forward << ") / (" << Distance_upperbounds_Forward << " - " << Distance_lowerbounds_Forward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC - 1] << "_Y_Scaled = (" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y - " << Distance_lowerbounds_Forward << ") / (" << Distance_upperbounds_Forward << " - " << Distance_lowerbounds_Forward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC - 1] << "_Z_Scaled = (" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z - " << Distance_lowerbounds_Forward << ") / (" << Distance_upperbounds_Forward << " - " << Distance_lowerbounds_Forward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC - 1] << "_VX_Scaled = (" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX - " << Velocity_lowerbounds_Forward << ") / (" << Velocity_upperbounds_Forward << " - " << Velocity_lowerbounds_Forward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC - 1] << "_VY_Scaled = (" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY - " << Velocity_lowerbounds_Forward << ") / (" << Velocity_upperbounds_Forward << " - " << Velocity_lowerbounds_Forward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC - 1] << "_VZ_Scaled = (" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ - " << Velocity_lowerbounds_Forward << ") / (" << Velocity_upperbounds_Forward << " - " << Velocity_lowerbounds_Forward << ")" << endl;
				GMATfile << "	FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = " << journeys[j].phases[p].state_at_beginning_of_phase[6] / (options.maximum_mass) << endl;
				GMATfile << endl;

				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Backward s/c scaled initial conditions" << endl;
				GMATfile <<"	" << SC_created[index_SC] << "_X_Scaled = (" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X - " << Distance_lowerbounds_Backward << ") / (" << Distance_upperbounds_Backward << " - " << Distance_lowerbounds_Backward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC] << "_Y_Scaled = (" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y - " << Distance_lowerbounds_Backward << ") / (" << Distance_upperbounds_Backward << " - " << Distance_lowerbounds_Backward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC] << "_Z_Scaled = (" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z - " << Distance_lowerbounds_Backward << ") / (" << Distance_upperbounds_Backward << " - " << Distance_lowerbounds_Backward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC] << "_VX_Scaled = (" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX - " << Velocity_lowerbounds_Backward << ") / (" << Velocity_upperbounds_Backward << " - " << Velocity_lowerbounds_Backward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC] << "_VY_Scaled = (" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY - " << Velocity_lowerbounds_Backward << ") / (" << Velocity_upperbounds_Backward << " - " << Velocity_lowerbounds_Backward << ")" << endl;
				GMATfile <<"	" << SC_created[index_SC] << "_VZ_Scaled = (" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ - " << Velocity_lowerbounds_Backward << ") / (" << Velocity_upperbounds_Backward << " - " << Velocity_lowerbounds_Backward << ")" << endl;
				GMATfile << "	FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled = " << journeys[j].phases[p].state_at_end_of_phase[6] / (options.maximum_mass) << endl;
				GMATfile << endl;	
				++index_SC;
				++index_body_visited;

				//initial guess for inter-phase control
				//this means thrust vectors (and sometimes Isp) for low-thrust phases
				//or burn index for MGADSM and MGANDSM phases. No control for MGA phases
				journeys[j].phases[p].output_GMAT_inter_phase_control_initial_guess(j, p, options, GMATfile);
			}
		}
		GMATfile << "EndScript" << endl;
		GMATfile << endl;
	
		//begin optimization loop
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%------------ Initial Boundary Constraints" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;
		GMATfile << "Optimize 'OptimizeSequence' NLPObject {SolveMode = Solve, ExitMode = DiscardAndContinue}" << endl;
		GMATfile << endl;
		index_SC = 0;
		index_body_visited = 0;

		//vary the launch date and journey wait times
		GMATfile << "	%Vary launch epoch and journey wait times" << endl;
		GMATfile << "	Vary 'VaryLaunchEpoch_Scaled' NLPObject(LaunchEpoch_Scaled = LaunchEpoch_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;

		//for each journey vary each subsequent journey wait time
		for (int j = 1; j < options.number_of_journeys; ++j)
		{			
			GMATfile << "	Vary 'VaryJourney" << j + 1 << "_WaitTime_Scaled' NLPObject(Journey" << j + 1 << "_WaitTime_Scaled = Journey" << j + 1 << "_WaitTime_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		}

		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{			
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				prefix.clear();
				prefixstream << "j" << j << "p" << p << ": ";
				prefix = prefixstream.str();
			
				//find index in Xdescriptions where the time bounds for that journey are located
				for (int iX = 0; iX < Xdescriptions.size(); ++iX)
				{
					if (Xdescriptions[iX] == prefix + "launch epoch (MJD)")
					{
						//define these bounds temporarily
						LaunchDate_lowerbounds = Xlowerbounds[iX];
						LaunchDate_upperbounds = Xupperbounds[iX];
					}
				}
				GMATfile << endl;

				//insert unscaled epoch values
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " s/c unscaled epochs" << endl;

				//if subsequent journey, add unscaled wait times
				if (j > 0)
				{
					GMATfile << "	'Calc" << SC_created[index_SC] << ".Epoch.TAIModJulian' " << SC_created[index_SC] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_upperbounds - LaunchDate_lowerbounds) << " + " << LaunchDate_lowerbounds << ") / 86400.0";
					for (int jj = 0; jj < j + 1; ++jj)
					{
						GMATfile << " + (Journey" << jj + 1 << "_WaitTime_Scaled * " << (options.journey_wait_time_bounds[jj][1] - options.journey_wait_time_bounds[jj][0]) << " + " << options.journey_wait_time_bounds[jj][0] << ")";
					}
					GMATfile << " + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
					//add time of flight for backward s/c (end of phase)
					CumulatedTOF = CumulatedTOF + journeys[j].phases[p].TOF;
					GMATfile << "	'Calc" << SC_created[index_SC + 1] << ".Epoch.TAIModJulian' " << SC_created[index_SC + 1] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_upperbounds - LaunchDate_lowerbounds) << " + " << LaunchDate_lowerbounds << ") / 86400.0";
					for (int jj = 0; jj < j + 1; ++jj)
					{
						GMATfile << " + (Journey" << jj + 1 << "_WaitTime_Scaled * " << (options.journey_wait_time_bounds[jj][1] - options.journey_wait_time_bounds[jj][0]) << " + " << options.journey_wait_time_bounds[jj][0] << ")";
					}
					GMATfile << " + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
				}
				//for the first journey
				else 
				{
					GMATfile << "	'Calc" << SC_created[index_SC] << ".Epoch.TAIModJulian' " << SC_created[index_SC] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_upperbounds - LaunchDate_lowerbounds) << " + " << LaunchDate_lowerbounds << ") / 86400.0 + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
					//add time of flight for backward s/c (end of phase)
					CumulatedTOF = CumulatedTOF + journeys[j].phases[p].TOF;
					GMATfile << "	'Calc" << SC_created[index_SC + 1] << ".Epoch.TAIModJulian' " << SC_created[index_SC + 1] << ".Epoch.TAIModJulian = (LaunchEpoch_Scaled * " << (LaunchDate_upperbounds - LaunchDate_lowerbounds) << " + " << LaunchDate_lowerbounds << ") / 86400.0 + " << CumulatedTOF / 86400.0 << " + 2400000.5 - 2430000" << endl;
				}
				index_SC = index_SC + 2;
				GMATfile << endl;
			} 
		}
		//store final epoch value
		GMATfile << "	'CalcFinalEpoch' FinalEpoch = " << SC_created[index_SC - 1] << ".Epoch.TAIModJulian" << endl;

		index_body_visited = 0;
		index_SC = 0;
		
		//add inequality constraints for min flyby altitude and ensure that each s/c begins at periapse
		GMATfile << "	%Add periapse and flyby altitude constraints for each phase" << endl;
		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				if ((index_SC + 1 % 2) && (p < journeys[j].number_of_phases) && (p != 0))
				{
					GMATfile << "	'Constraint" << SC_created[index_SC] << "_PeriapseRadius_Scaled' " << SC_created[index_SC] << "_PeriapseRadius_Scaled = " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << ".RadPer / " << missionbodies[index_body_visited].minimum_safe_flyby_altitude << endl;
					GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_PeriapseRadius_Scaled' NLPObject(" << SC_created[index_SC] << "_PeriapseRadius_Scaled >= 1)" << endl;
					GMATfile << "	'Calc" << SC_created[index_SC] << "_RdotV_Scaled' " << SC_created[index_SC] << "_RdotV_Scaled = (" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.X * " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX + " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y * " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY + " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z * " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ + " << sqrt(2 * missionbodies[index_body_visited].mu * missionbodies[index_body_visited].radius) << ") / " << 2 * sqrt(2 * missionbodies[index_body_visited].mu * missionbodies[index_body_visited].radius) << endl;
					GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_RdotV_Scaled' NLPObject(" << SC_created[index_SC] << "_RdotV_Scaled = 0.5);" << endl;
					GMATfile << endl;
				}
				index_SC = index_SC + 2;
				index_body_visited++;
			}
		}
		index_body_visited = 0;
		index_SC = 0;

		//vary the initial states of s/c
		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				//vary all forward s/c
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " vary forward s/c scaled states" << endl;
				GMATfile << "	Vary 'Vary" << SC_created[index_SC] << "_X_Scaled' NLPObject(" << SC_created[index_SC] << "_X_Scaled = " << SC_created[index_SC] << "_X_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'Vary" << SC_created[index_SC] << "_Y_Scaled' NLPObject(" << SC_created[index_SC] << "_Y_Scaled = " << SC_created[index_SC] << "_Y_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'Vary" << SC_created[index_SC] << "_Z_Scaled' NLPObject(" << SC_created[index_SC] << "_Z_Scaled = " << SC_created[index_SC] << "_Z_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'Vary" << SC_created[index_SC] << "_VX_Scaled' NLPObject(" << SC_created[index_SC] << "_VX_Scaled = " << SC_created[index_SC] << "_VX_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'Vary" << SC_created[index_SC] << "_VY_Scaled' NLPObject(" << SC_created[index_SC] << "_VY_Scaled = " << SC_created[index_SC] << "_VY_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'Vary" << SC_created[index_SC] << "_VZ_Scaled' NLPObject(" << SC_created[index_SC] << "_VZ_Scaled = " << SC_created[index_SC] << "_VZ_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << "	Vary 'VaryFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled' NLPObject(FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
				GMATfile << endl;
				index_SC = index_SC + 2;
			}
		}
		//vary last s/c state for mission (backward from last body)
		GMATfile <<"	%Journey #"<< options.number_of_journeys <<", Phase #"<< journeys[options.number_of_journeys - 1].number_of_phases << " vary backward s/c scaled states" << endl;
		GMATfile << "	Vary 'Vary" << SC_created[index_SC - 1] << "_X_Scaled' NLPObject(" << SC_created[index_SC - 1] << "_X_Scaled = " << SC_created[index_SC - 1] << "_X_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		GMATfile << "	Vary 'Vary" << SC_created[index_SC - 1] << "_Y_Scaled' NLPObject(" << SC_created[index_SC - 1] << "_Y_Scaled = " << SC_created[index_SC - 1] << "_Y_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		GMATfile << "	Vary 'Vary" << SC_created[index_SC - 1] << "_Z_Scaled' NLPObject(" << SC_created[index_SC - 1] << "_Z_Scaled = " << SC_created[index_SC - 1] << "_Z_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		GMATfile << "	Vary 'Vary" << SC_created[index_SC - 1] << "_VX_Scaled' NLPObject(" << SC_created[index_SC - 1] << "_VX_Scaled = " << SC_created[index_SC - 1] << "_VX_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		GMATfile << "	Vary 'Vary" << SC_created[index_SC - 1] << "_VY_Scaled' NLPObject(" << SC_created[index_SC - 1] << "_VY_Scaled = " << SC_created[index_SC - 1] << "_VY_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		GMATfile << "	Vary 'Vary" << SC_created[index_SC - 1] << "_VZ_Scaled' NLPObject(" << SC_created[index_SC - 1] << "_VZ_Scaled = " << SC_created[index_SC - 1] << "_VZ_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		GMATfile << "	Vary 'VaryFuelMass_Journey" << options.number_of_journeys << "Phase" << journeys[options.number_of_journeys - 1].number_of_phases << "Backward_Scaled' NLPObject(FuelMass_Journey" << options.number_of_journeys << "Phase" << journeys[options.number_of_journeys - 1].number_of_phases << "Backward_Scaled = FuelMass_Journey" << options.number_of_journeys << "Phase" << journeys[options.number_of_journeys - 1].number_of_phases << "Backward_Scaled, {Perturbation = 0.00001, Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
		GMATfile << "	'CalcFinalMass_Scaled' FinalMass_Scaled = FuelMass_Journey" << options.number_of_journeys << "Phase" << journeys[options.number_of_journeys - 1].number_of_phases << "Backward_Scaled" << endl;
		GMATfile << endl;

		index_SC = 0;
		//now we must define states to their unscaled values
		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				//define bounds for scaling variables
				for (int iX = 0; iX < Xdescriptions.size(); ++iX)
				{
					//find index in Xdescriptions where the flyby altitude bounds for that phase are located
					if (Xdescriptions[iX] == prefix + "flyby altitude constraint (above minimum altitude but below [100x/300x] altitude for [rocky/gas] planets")
					{
						Distance_lowerbounds_Forward =  Xlowerbounds[iX] * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_upperbounds_Forward = -Xlowerbounds[iX] * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_lowerbounds_Backward =  Xlowerbounds[iX] * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
						Distance_upperbounds_Backward = -Xlowerbounds[iX] * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
					}
					//find index in Xdescriptions where the flyby velocity bounds for that phase are located
					if (Xdescriptions[iX] == prefix + "initial velocity increment x")
					{
						Velocity_lowerbounds_Forward = -sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + Xlowerbounds[iX] * Xlowerbounds[iX]);
						Velocity_upperbounds_Forward =  sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + Xupperbounds[iX] * Xupperbounds[iX]);
						Velocity_lowerbounds_Backward = -sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + Xlowerbounds[iX] * Xlowerbounds[iX]);
						Velocity_upperbounds_Backward =  sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + Xupperbounds[iX] * Xupperbounds[iX]);
					}			
				}
				if ((p == 0)||(p == journeys[j].number_of_phases - 1))
				{
					//hardcoded in for journey arrivals and departures :-/
					if (missionbodies[index_body_visited].mass < 1.0e25)
					{
						Distance_lowerbounds_Forward = -10 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_upperbounds_Forward =  10 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
					}
					else
					{
						Distance_lowerbounds_Forward = -300 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
						Distance_upperbounds_Forward =  300 * (missionbodies[index_body_visited].minimum_safe_flyby_altitude + missionbodies[index_body_visited].radius);
					}
					if (missionbodies[index_body_visited + 1].mass < 1.0e25)
					{
						Distance_lowerbounds_Backward = -10 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
						Distance_upperbounds_Backward =  10 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
					}
					else
					{
						Distance_lowerbounds_Backward = -300 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
						Distance_upperbounds_Backward =  300 * (missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude + missionbodies[index_body_visited + 1].radius);
					}		
					Velocity_lowerbounds_Forward = -sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + (-25) * (-25));
					Velocity_upperbounds_Forward =  sqrt(2 * missionbodies[index_body_visited].mu / (missionbodies[index_body_visited].radius + missionbodies[index_body_visited].minimum_safe_flyby_altitude) + 25 * 25);
					Velocity_lowerbounds_Backward = -sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + (-25) * (-25));
					Velocity_upperbounds_Backward =  sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + missionbodies[index_body_visited + 1].minimum_safe_flyby_altitude) + 25 * 25);
				}

				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " set forward s/c unscaled states" << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.X' " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.X = " << SC_created[index_SC] << "_X_Scaled * " << (Distance_upperbounds_Forward - Distance_lowerbounds_Forward) << " + " << Distance_lowerbounds_Forward << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y' " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y = " << SC_created[index_SC] << "_Y_Scaled * " << (Distance_upperbounds_Forward - Distance_lowerbounds_Forward) << " + " << Distance_lowerbounds_Forward << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z' " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z = " << SC_created[index_SC] << "_Z_Scaled * " << (Distance_upperbounds_Forward - Distance_lowerbounds_Forward) << " + " << Distance_lowerbounds_Forward << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX' " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX = " << SC_created[index_SC] << "_VX_Scaled * " << (Velocity_upperbounds_Forward - Velocity_lowerbounds_Forward) << " + " << Velocity_lowerbounds_Forward << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY' " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY = " << SC_created[index_SC] << "_VY_Scaled * " << (Velocity_upperbounds_Forward - Velocity_lowerbounds_Forward) << " + " << Velocity_lowerbounds_Forward << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ' " << SC_created[index_SC] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ = " << SC_created[index_SC] << "_VZ_Scaled * " << (Velocity_upperbounds_Forward - Velocity_lowerbounds_Forward) << " + " << Velocity_lowerbounds_Forward << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.FuelMass' " << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.FuelMass = FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled * " << (options.maximum_mass) << endl;
				GMATfile << endl;
				index_SC++;

				//match forward and backward s/c states and epochs		
				if (index_SC + 1 < (SC_created.size() - 1))
				{			
					if (p < journeys[j].number_of_phases - 1)
					{
						//match only odd phase backward s/c == subsequent even phase forward s/c
						GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " match forward and backward s/c unscaled states and epochs" << endl;
						GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X' " << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X = " << SC_created[index_SC + 1] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.X;" << endl; 
						GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y' " << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y = " << SC_created[index_SC + 1] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Y;" << endl; 
						GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z' " << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z = " << SC_created[index_SC + 1] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.Z;" << endl; 
						GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX' " << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX = " << SC_created[index_SC + 1] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VX;" << endl; 
						GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY' " << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY = " << SC_created[index_SC + 1] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VY;" << endl; 
						GMATfile << "	'Calc" << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ' " << SC_created[index_SC] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ = " << SC_created[index_SC + 1] << "." << missionbodies[index_body_visited + 1].name << "J2000Eq.VZ;" << endl; 
						GMATfile << "	'Calc" << SC_created[index_SC] << ".Epoch.TAIModJulian' " << SC_created[index_SC] << ".Epoch.TAIModJulian = " << SC_created[index_SC + 1] << ".Epoch.TAIModJulian;" << endl; 
						if (p < journeys[j].number_of_phases - 1)
						{
							GMATfile << "	'Calc" << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass' " << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass = " << SC_created[index_SC + 1] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 2 << "Forward.FuelMass" << endl;
						}
						else
						{
							GMATfile << "	'Calc" << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass' " << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass = " << SC_created[index_SC + 1] << ".FuelTank_Journey" << j + 2 << "Phase" << 1 << "Forward.FuelMass" << endl;
						}			
						GMATfile << endl;
					}
				}
				index_SC++;
				index_body_visited++;
			}
		}
		//set last backward states to unscaled values
		GMATfile <<"	%Journey #"<< options.number_of_journeys <<", Phase #"<< journeys[options.number_of_journeys - 1].number_of_phases << " set backward s/c unscaled states" << endl;
		GMATfile << "	'Calc" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.X' " << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.X = " << SC_created[index_SC - 1] << "_X_Scaled * " << (Distance_upperbounds_Backward - Distance_lowerbounds_Backward) << " + " << Distance_lowerbounds_Backward << endl;
		GMATfile << "	'Calc" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y' " << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.Y = " << SC_created[index_SC - 1] << "_Y_Scaled * " << (Distance_upperbounds_Backward - Distance_lowerbounds_Backward) << " + " << Distance_lowerbounds_Backward << endl;
		GMATfile << "	'Calc" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z' " << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.Z = " << SC_created[index_SC - 1] << "_Z_Scaled * " << (Distance_upperbounds_Backward - Distance_lowerbounds_Backward) << " + " << Distance_lowerbounds_Backward << endl;
		GMATfile << "	'Calc" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX' " << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VX = " << SC_created[index_SC - 1] << "_VX_Scaled * " << (Velocity_upperbounds_Backward - Velocity_lowerbounds_Backward) << " + " << Velocity_lowerbounds_Backward << endl;
		GMATfile << "	'Calc" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY' " << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VY = " << SC_created[index_SC - 1] << "_VY_Scaled * " << (Velocity_upperbounds_Backward - Velocity_lowerbounds_Backward) << " + " << Velocity_lowerbounds_Backward << endl;
		GMATfile << "	'Calc" << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ' " << SC_created[index_SC - 1] << "." << missionbodies[index_body_visited].name << "J2000Eq.VZ = " << SC_created[index_SC - 1] << "_VZ_Scaled * " << (Velocity_upperbounds_Backward - Velocity_lowerbounds_Backward) << " + " << Velocity_lowerbounds_Backward << endl;
		GMATfile << "	'Calc" << SC_created[index_SC - 1] << ".FuelTank_Journey" <<  options.number_of_journeys << "Phase" << journeys[options.number_of_journeys - 1].number_of_phases << "Backward.FuelMass' " << SC_created[index_SC - 1] << ".FuelTank_Journey" <<  options.number_of_journeys << "Phase" << journeys[options.number_of_journeys - 1].number_of_phases << "Backward.FuelMass = FuelMass_Journey" << options.number_of_journeys << "Phase" << journeys[options.number_of_journeys - 1].number_of_phases << "Backward_Scaled * " << (options.maximum_mass) << endl;
		GMATfile << endl;
		GMATfile << endl;

		GMATfile << "	%-------------------------------------------------------------------------" << endl;
		GMATfile << "	%---------- Propagation" << endl;
		GMATfile << "	%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		index_SC = 0;
		index_body_visited = 0;
		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				journeys[j].phases[p].output_GMAT_propagate_phase_forward_and_back(j, p, missionbodies, index_body_visited, SC_created, index_SC, options, GMATfile);
			}
		}
		GMATfile << endl;
		GMATfile << endl;

		//apply all boundary constraints at match points
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Final Boundary Constraints" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		index_SC = 0;
		index_body_visited = 0;
		//for each journey
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			//for each phase
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				//scale final forward states for constraints
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Forward, scale final states for constraints" << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "_X_Scaled' " << SC_created[index_SC] << "_X_Scaled = " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.X / " << TheUniverse[j].LU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "_Y_Scaled' " << SC_created[index_SC] << "_Y_Scaled = " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.Y / " << TheUniverse[j].LU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "_Z_Scaled' " << SC_created[index_SC] << "_Z_Scaled = " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.Z / " << TheUniverse[j].LU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "_VX_Scaled' " << SC_created[index_SC] << "_VX_Scaled = " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VX / " << TheUniverse[j].LU << " * " << TheUniverse[j].TU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "_VY_Scaled' " << SC_created[index_SC] << "_VY_Scaled = " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VY / " << TheUniverse[j].LU << " * " << TheUniverse[j].TU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC] << "_VZ_Scaled' " << SC_created[index_SC] << "_VZ_Scaled = " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VZ / " << TheUniverse[j].LU << " * " << TheUniverse[j].TU << endl;
				GMATfile << endl;
						
				//scale final backward states for constraints
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Backward, scale final states for constraints" << endl;
				GMATfile << "	'Calc" << SC_created[index_SC + 1] << "_X_Scaled' " << SC_created[index_SC + 1] << "_X_Scaled = " << SC_created[index_SC + 1] << "." << missionbodies[0].central_body_name << "J2000Eq.X / " << TheUniverse[j].LU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC + 1] << "_Y_Scaled' " << SC_created[index_SC + 1] << "_Y_Scaled = " << SC_created[index_SC + 1] << "." << missionbodies[0].central_body_name << "J2000Eq.Y / " << TheUniverse[j].LU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC + 1] << "_Z_Scaled' " << SC_created[index_SC + 1] << "_Z_Scaled = " << SC_created[index_SC + 1] << "." << missionbodies[0].central_body_name << "J2000Eq.Z / " << TheUniverse[j].LU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC + 1] << "_VX_Scaled' " << SC_created[index_SC + 1] << "_VX_Scaled = " << SC_created[index_SC + 1] << "." << missionbodies[0].central_body_name << "J2000Eq.VX / " << TheUniverse[j].LU << " * " << TheUniverse[j].TU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC + 1] << "_VY_Scaled' " << SC_created[index_SC + 1] << "_VY_Scaled = " << SC_created[index_SC + 1] << "." << missionbodies[0].central_body_name << "J2000Eq.VY / " << TheUniverse[j].LU << " * " << TheUniverse[j].TU << endl;
				GMATfile << "	'Calc" << SC_created[index_SC + 1] << "_VZ_Scaled' " << SC_created[index_SC + 1] << "_VZ_Scaled = " << SC_created[index_SC + 1] << "." << missionbodies[0].central_body_name << "J2000Eq.VZ / " << TheUniverse[j].LU << " * " << TheUniverse[j].TU << endl;
				GMATfile << endl;

				//apply match point constraints
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Backward, apply match point constraints" << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_X_Scaled' NLPObject(" << SC_created[index_SC] << "_X_Scaled = " << SC_created[index_SC + 1] << "_X_Scaled)" << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_Y_Scaled' NLPObject(" << SC_created[index_SC] << "_Y_Scaled = " << SC_created[index_SC + 1] << "_Y_Scaled)" << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_Z_Scaled' NLPObject(" << SC_created[index_SC] << "_Z_Scaled = " << SC_created[index_SC + 1] << "_Z_Scaled)" << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_VX_Scaled' NLPObject(" << SC_created[index_SC] << "_VX_Scaled = " << SC_created[index_SC + 1] << "_VX_Scaled)" << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_VY_Scaled' NLPObject(" << SC_created[index_SC] << "_VY_Scaled = " << SC_created[index_SC + 1] << "_VY_Scaled)" << endl;
				GMATfile << "	NonlinearConstraint 'Constraint" << SC_created[index_SC] << "_VZ_Scaled' NLPObject(" << SC_created[index_SC] << "_VZ_Scaled = " << SC_created[index_SC + 1] << "_VZ_Scaled)" << endl;
				GMATfile << endl;

				//apply fuel mass constraints
				GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " apply fuel mass constraint at match point" << endl;
				GMATfile << "	'CalcFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled' FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = " << SC_created[index_SC] << ".TotalMass / " << (2* options.maximum_mass) << endl;
				GMATfile << "	'CalcFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled' FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled = " << SC_created[index_SC] << ".TotalMass / " << (2 * options.maximum_mass) << endl;
				GMATfile << "	NonlinearConstraint 'ConstraintFuelMass_Journey" << j + 1 << "Phase" << p + 1 << "' NLPObject(FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled = FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled)" << endl;
				GMATfile << endl;
				GMATfile << endl;
				index_SC = index_SC + 2;
				index_body_visited++;
			}
		}

		//optimize the user-defined objective function
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << "%---------- Objective Function" << endl;
		GMATfile << "%-------------------------------------------------------------------------" << endl;
		GMATfile << endl;

		switch (options.objective_type)
		{
		case 0: // minimize delta-V
			GMATfile << "	%Objective function is to minimize delta-V" << endl;
			GMATfile << "	'CalcObjectiveFunction' ObjectiveFunction = (0";
			for (int j = 0; j < options.number_of_journeys; ++j)
			{
				for (int p = 0; p < journeys[j].number_of_phases; ++p)
				{
					for (int step = 0; step < options.num_timesteps; ++step)
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
			GMATfile << "	'CalcObjectiveFunction' ObjectiveFunction = (FinalEpoch - LaunchEpoch_Scaled * " << (LaunchDate_upperbounds - LaunchDate_lowerbounds) << " + " << LaunchDate_lowerbounds << ") / " << options.total_flight_time_bounds[1] << endl;
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
	}//end method phase::output_GMAT_mission()


	void phase::output_GMAT_spacecraft(int& j, int& p, vector<string>& SC_created, int& index_SC, vector<EMTG::Astrodynamics::body>& missionbodies, int& index_body_visited, std::ofstream& GMATfile)
	{
		//first create string for name of SC, add it to vector of all SC names
		std::ostringstream SC_name_temp;
		SC_name_temp << "SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward";
		SC_created.push_back(SC_name_temp.str());

		//fill in all necessary parameters
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << ", forward propagated s/c" << endl;
		GMATfile << "Create Spacecraft " << SC_created[index_SC] << ";" << endl;
		GMATfile << SC_created[index_SC] << ".DateFormat = TAIModJulian;" << endl;
		GMATfile << SC_created[index_SC] << ".Epoch = " << (phase_start_epoch / 86400.0) + 2400000.5 - 2430000 << ";" << endl;
		GMATfile << SC_created[index_SC] << ".DryMass = 0" << endl;
		GMATfile << SC_created[index_SC] << ".CoordinateSystem = " << missionbodies[index_body_visited].name << "J2000Eq;" << endl;			
		GMATfile << SC_created[index_SC] << ".Tanks = {FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward};" << endl;
		GMATfile << SC_created[index_SC] << ".Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward};" << endl;
		++index_SC;
		GMATfile << endl;

		SC_name_temp.str("");
		SC_name_temp << "SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward";
		SC_created.push_back(SC_name_temp.str());

		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << ", backward propagated s/c" << endl;
		GMATfile << "Create Spacecraft " << SC_created[index_SC] << ";" << endl;
		GMATfile << SC_created[index_SC] << ".DateFormat = TAIModJulian;" << endl;
		GMATfile << SC_created[index_SC] << ".Epoch = " << (phase_start_epoch + TOF) / 86400.0 + 2400000.5 - 2430000 << ";" << endl;
		GMATfile << SC_created[index_SC] << ".DryMass = 0" << endl;
		GMATfile << SC_created[index_SC] << ".CoordinateSystem = " << missionbodies[index_body_visited + 1].name << "J2000Eq;" << endl;			
		GMATfile << SC_created[index_SC] << ".Tanks = {FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward};" << endl;
		GMATfile << SC_created[index_SC] << ".Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward};" << endl;
		++index_SC;
		++index_body_visited;

		GMATfile << endl;
		GMATfile << endl;
	}//end method phase::output_GMAT_spacecraft()

	void phase::output_GMAT_fueltank_and_thruster(int& j, int& p, vector<EMTG::Astrodynamics::body>& missionbodies, int& index_body_visited, std::ofstream& GMATfile)
	{
		//create two tanks and two thrusters for each phase in each journey (forward + backward)
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << ", forward s/c tanks and thrusters" << endl;
		GMATfile << "Create FuelTank FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.AllowNegativeFuelMass = false;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.Volume = 10;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.FuelMass = " << state_at_beginning_of_phase[6] << endl;
		GMATfile << endl;

		GMATfile << "Create Thruster Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.CoordinateSystem = " << missionbodies[index_body_visited].name << "J2000Eq" << endl;
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

		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << ", backward s/c tanks and thrusters" << endl;
		GMATfile << "Create FuelTank FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.AllowNegativeFuelMass = false;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.Volume = 10;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass = " << state_at_end_of_phase[6] << endl;
			
		GMATfile << endl;
		GMATfile << "Create Thruster Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.CoordinateSystem = " << missionbodies[index_body_visited + 1].name << "J2000Eq" << endl;
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
		++index_body_visited;
	}//end method phase::output_GMAT_fueltank_and_thruster

	void phase::output_GMAT_burn_objects(int& j, int& p, std::ofstream& GMATfile)
	{
		//create two finite burn objects for each phase in each journey (forward + backward)
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << ", forward and backward finite burn objects" << endl;
		GMATfile << "Create FiniteBurn FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
		GMATfile << "FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward.Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward};" << endl;
		GMATfile << "Create FiniteBurn FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward;" << endl;
		GMATfile << "FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward.Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward};" << endl;
		GMATfile << endl;
	}//end method phase::output_GMAT_burn_objects

	void phase::output_GMAT_create_state_and_time_variables(int& j, int& p, std::ofstream& GMATfile)
	{
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << " time variables" << endl;
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
			
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << " forward variables" << endl;
		//created scaled state variables
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_X_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Y_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Z_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_VX_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_VY_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Forward_VZ_Scaled" << endl;
		GMATfile << "Create Variable FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Forward_Scaled" << endl;
		GMATfile << endl;
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << " backward variables" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_X_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Y_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Z_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_VX_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_VY_Scaled" << endl;
		GMATfile << "Create Variable SC_Journey" << j + 1 << "Phase" << p + 1 << "Backward_VZ_Scaled" << endl;
		GMATfile << "Create Variable FuelMass_Journey" << j + 1 << "Phase" << p + 1 << "Backward_Scaled" << endl;
		GMATfile << endl;
	}//end method phase::output_GMAT_create_state_and_time_variables

	void phase::output_GMAT_create_interphase_control_variables(int& j, int& p, missionoptions& options, std::ofstream& GMATfile)
	{
		//create arrays to hold thrust direction vector and unit vectors for each timestep
		GMATfile << "Create Array ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "[" << options.num_timesteps << ", 3];" << endl;
		GMATfile << "Create Array ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "[" << options.num_timesteps << ", 1];" << endl;
	}//end method phase::output_GMAT_create_interphase_control_variables

	void phase::output_GMAT_inter_phase_control_initial_guess(int& j, int& p, missionoptions& options, std::ofstream& GMATfile)
	{
		//initialize thrust vector directions
		//define thrust unit vector bounds to scale to between 0 and 1
		double ThrustUnitVector_lowerbounds = -1;
		double ThrustUnitVector_upperbounds =  1;

		//propagate forward s/c using finite burns (must scale unit vectors to between 0 and 1)
		GMATfile << endl;
		GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Forward s/c scaled thrust vectors" << endl;
		for (int step = 0; step < (options.num_timesteps / 2); ++step)
		{
			GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = " << (control[step][0] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
			GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = " << (control[step][1] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
			GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = " << (control[step][2] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
			GMATfile << endl;
		}
		//propagate backward s/c using finite burns (must scale unit vectors to between 0 and 1)
		GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Backward s/c scaled thrust vectors" << endl;
		for (int step = options.num_timesteps - 1; (step >= options.num_timesteps / 2); --step)
		{
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << ", Time Step #" << step + 1 << ", Backward Propagation" << endl;
			GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = " << (control[step][0] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
			GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = " << (control[step][1] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
			GMATfile << "	ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = " << (control[step][2] - ThrustUnitVector_lowerbounds) / (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << endl;
			GMATfile << endl;
		}
	}//end method phase::output_GMAT_inter_phase_control_initial_guess

	void phase::output_GMAT_propagate_phase_forward_and_back(int& j, int& p, vector<Astrodynamics::body>& missionbodies, int& index_body_visited, vector<string>& SC_created, int& index_SC, missionoptions& options, std::ofstream& GMATfile)
	{
		//define some variables
		double periapse_velocity_magnitude;
		math::Matrix<double> Vinf_out(3, 1), Vinf_in(3, 1);
		double delta_t;
		double boundary_state[6];
		int index_delta_t;

		//define thrust unit vector bounds to scale to between 0 and 1
		double ThrustUnitVector_lowerbounds = -1;
		double ThrustUnitVector_upperbounds =  1;

		//must handle discontinuities in plots
		GMATfile << "	PenUp 'PenUp' " << missionbodies[0].central_body_name << "View";
		for (int b = 0; b < missionbodies.size(); ++b)
		{
			GMATfile << " " << missionbodies[b].name << "View";
		}
		GMATfile << ";" << endl;
		GMATfile << endl;

		//propagate forward s/c using finite burns (must scale unit vectors to between 0 and 1)
		for (int step = 0; step < (options.num_timesteps / 2); ++step)
		{
			GMATfile << "	%Journey #" << j + 1 << ", Phase #" << p + 1 << ", Time Step #" << step + 1 << ", Forward Propagation";
			GMATfile << endl;
			GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	Vary 'VaryThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3)' NLPObject(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) = ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3), {Perturbation = 0.00001,  Lower = " << 0 << ", Upper = " << 1 << ", MaxStep = 0.1})" << endl;
			GMATfile << "	'CalcThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "' ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) = sqrt((ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " << "(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2 + " << "(ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") ^ 2);" << endl;

			//assign thrust unit directions to thruster
			GMATfile << "	NonlinearConstraint 'ConstraintThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)' NLPObject(ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) <= 1)" << endl;
			GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
			GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 2) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
			GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3 = (ThrustVector_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 3) * " << (ThrustUnitVector_upperbounds - ThrustUnitVector_lowerbounds) << " + " << ThrustUnitVector_lowerbounds << ") / ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;
			GMATfile << "	'CalcThruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1' Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1 = ThrusterMaxThrust * ThrustUnitVectorMag_Journey" << j + 1 << "Phase" << p + 1 << "(" << step + 1 << ", 1)" << endl;

			//begin propagation of finite burn

			//locate the position of the bodies at beginning and end of phase
			missionbodies[index_body_visited].locate_body(phase_start_epoch, boundary_state, false, &options);

			//calculate time spent in SOI of body and during which timestep the boundary occurs
			for (int k = 0; k < 3; ++k)
			{
				Vinf_out(k) = state_at_beginning_of_phase[k + 3] - boundary_state[k + 3];
			}
			periapse_velocity_magnitude = sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + flyby_altitude) + Vinf_out.dot(Vinf_out));
			delta_t = (missionbodies[index_body_visited].r_SOI / periapse_velocity_magnitude) / 86400;
			index_delta_t = floor(delta_t / (TOF / options.num_timesteps));

			//check if inside SOI
			//TODO:: remove second part of if-statement after fixing arrival/departures
			if ((step < index_delta_t) && (p != 0))
			{
				GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[index_body_visited].name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int b = 0; b < missionbodies.size(); ++b)
				{
					GMATfile << " " << missionbodies[b].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[index_body_visited].name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << TOF / options.num_timesteps << "};" << endl;
			}

			//check if it leaves SOI during time step  
			else if ((step == index_delta_t) && (p != 0))
			{
				//if so, propagate only until outside estimated SOI of body
				GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[index_body_visited].name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int b = 0; b < missionbodies.size(); ++b)
				{
					GMATfile << " " << missionbodies[b].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[index_body_visited].name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << delta_t - (TOF / options.num_timesteps) * index_delta_t << "};" << endl;
					
				//then propagate the rest of the timestep
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int index_body_visited = 0; index_body_visited < missionbodies.size(); ++index_body_visited)
				{
					GMATfile << " " << missionbodies[index_body_visited].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << (TOF / options.num_timesteps) * (index_delta_t + 1) - delta_t << "};" << endl;
			}

			//otherwise propagate the full timestep
			else
			{
				GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int b = 0; b < missionbodies.size(); ++b)
				{
					GMATfile << " " << missionbodies[b].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << TOF / options.num_timesteps << "};" << endl;
			}
				
			//end finite burn
			GMATfile << "	EndFiniteBurn 'EndFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Forward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Forward(" << SC_created[index_SC] << ");" << endl;
				
			//report things for debugging
			GMATfile << "	Report 'Report_SpacecraftControl' Report_SpacecraftControl " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection1 " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection2 " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.ThrustDirection3 " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward.C1" << endl;
			GMATfile << "	Report 'Report_SpacecraftState' Report_SpacecraftState " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.X " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.Y " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.Z " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VX " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VY " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VZ " << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward.FuelMass;" << endl;
			GMATfile << "	PenUp 'PenUp' " << missionbodies[0].central_body_name << "View";
			for (int b = 0; b < missionbodies.size(); ++b)
			{
				GMATfile << " " << missionbodies[b].name << "View";
			}
			GMATfile << ";" << endl;
			GMATfile << endl;
		}
		++index_SC;


		//propagate backward s/c using finite burns (must scale unit vectors to between 0 and 1)
		for (int step = options.num_timesteps - 1; (step >= options.num_timesteps / 2); --step)
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
			missionbodies[index_body_visited + 1].locate_body(phase_end_epoch, boundary_state, false, &options);

			//calculate time spent in SOI of body and during which timestep the boundary occurs
			for (int k = 0; k < 3; ++k)
			{
				Vinf_in(k) = state_at_end_of_phase[k + 3] - boundary_state[k + 3];
			}
			periapse_velocity_magnitude = sqrt(2 * missionbodies[index_body_visited + 1].mu / (missionbodies[index_body_visited + 1].radius + flyby_altitude) + Vinf_in.dot(Vinf_in));
			delta_t = (missionbodies[index_body_visited + 1].r_SOI / periapse_velocity_magnitude) / 86400;
			index_delta_t = (options.num_timesteps - 1) - floor(delta_t / (TOF / options.num_timesteps));

			//check if inside SOI
			//TODO:: remove second part of if-statement after fixing arrival/departures
			if ((step > index_delta_t) && (p != options.number_of_phases[j] - 1))
			{
				GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[index_body_visited + 1].name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int b = 0; b < missionbodies.size(); ++b)
				{
					GMATfile << " " << missionbodies[b].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[index_body_visited + 1].name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << -TOF / options.num_timesteps << "};" << endl;
			}

			//check if it leaves SOI during time step  
			else if ((step == index_delta_t) && (p != options.number_of_phases[j] - 1))
			{
				//is so, propagate only until outside estimated SOI of body
				GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[index_body_visited + 1].name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int b = 0; b < missionbodies.size(); ++b)
				{
					GMATfile << " " << missionbodies[b].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[index_body_visited + 1].name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << -(delta_t - (TOF / options.num_timesteps) * ((options.num_timesteps - 1) - index_delta_t)) << "};" << endl;
					
				//then propagate the rest of the timestep
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int b = 0; b < missionbodies.size(); ++b)
				{
					GMATfile << " " << missionbodies[b].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << -((TOF / options.num_timesteps) * (options.num_timesteps - index_delta_t) - delta_t) << "};" << endl;
			}

			//otherwise propagate the full timestep
			else
			{
				GMATfile << "	BeginFiniteBurn 'BeginFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ");" << endl;
				GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
				for (int b = 0; b < missionbodies.size(); ++b)
				{
					GMATfile << " " << missionbodies[b].name << "View";
				}
				GMATfile << ";" << endl;
				GMATfile << "	Propagate 'PropJourney" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' BackProp " << missionbodies[0].central_body_name << "Prop(" << SC_created[index_SC] << ") {" << SC_created[index_SC] << ".ElapsedSecs = " << -(TOF / options.num_timesteps) << "};" << endl;
			}
				
			//end finite burn
			GMATfile << "	EndFiniteBurn 'EndFiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Step" << step + 1 << "Backward' FiniteBurn_Journey" << j + 1 << "Phase" << p + 1 << "Backward(" << SC_created[index_SC] << ");" << endl;
				
			//report things for debugging
			GMATfile << "	Report 'Report_SpacecraftControl' Report_SpacecraftControl " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection1 " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection2 " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.ThrustDirection3 " << SC_created[index_SC] << ".Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Backward.C1" << endl;
			GMATfile << "	Report 'Report_SpacecraftState' Report_SpacecraftState " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.X " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.Y " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.Z " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VX " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VY " << SC_created[index_SC] << "." << missionbodies[0].central_body_name << "J2000Eq.VZ " << SC_created[index_SC] << ".FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Backward.FuelMass;" << endl;
			GMATfile << "	PenUp 'PenUp' " << missionbodies[0].central_body_name << "View";
			for (int b = 0; b < missionbodies.size(); ++b)
			{
				GMATfile << " " << missionbodies[b].name << "View";
			}
			GMATfile << ";" << endl;
			GMATfile << endl;
		}
		++index_SC;
		++index_body_visited;

		GMATfile << "	PenDown 'PenDown' " << missionbodies[0].central_body_name << "View";
		for (int b = 0; b < missionbodies.size(); ++b)
		{
			GMATfile << " " << missionbodies[b].name << "View";
		}
		GMATfile << ";" << endl;
	}//end phase::output_GMAT_propagate_phase_forward_and_back

	//********************************************************************methods specific to MGA_DSM

	void MGA_DSM_phase::output_GMAT_burn_objects(int& j, int& p, std::ofstream& GMATfile)
	{
		//create an impulsive burn object for the phase
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << ", impulsive burn object" << endl;
		GMATfile << "Create ImpulsiveBurn  Impulsive_Journey" << j + 1 << "Phase" << p + 1 << ";" << endl;
		GMATfile << "ImpulsiveBurn_Journey" << j + 1 << "Phase" << p + 1 << ".Thrusters = {Thruster_Journey" << j + 1 << "Phase" << p + 1 << "};" << endl;
		GMATfile << endl;
	}//end method MGA_DSM_phase::output_GMAT_burn_objects

	void MGA_DSM_phase::output_GMAT_fueltank_and_thruster(int& j, int& p, vector<EMTG::Astrodynamics::body>& missionbodies, int& index_body_visited, std::ofstream& GMATfile)
	{
		//create two tanks and two thrusters for each phase in each journey (forward + backward)
		GMATfile << "% Journey #"<< j + 1 << ", Phase #"<< p + 1 << ", s/c tanks and thrusters" << endl;
		GMATfile << "Create FuelTank FuelTank_Journey" << j + 1 << "Phase" << p + 1 << ";" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << ".AllowNegativeFuelMass = false;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << ".Volume = 10;" << endl;
		GMATfile << "FuelTank_Journey" << j + 1 << "Phase" << p + 1 << ".FuelMass = " << state_at_beginning_of_phase[6] << endl;
		GMATfile << endl;

		GMATfile << "Create Thruster Thruster_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".CoordinateSystem = " << missionbodies[index_body_visited].name << "J2000Eq" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".ThrustDirection1 = 1;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".ThrustDirection2 = 0;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".ThrustDirection3 = 0;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".DutyCycle = 1;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".Tank = FuelTank_Journey" << j + 1 << "Phase" << p + 1 << "Forward;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".ThrustScaleFactor = 1;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".DecrementMass = true;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".C1 = .1;" << endl;
		GMATfile << "Thruster_Journey" << j + 1 << "Phase" << p + 1 << ".K1 = 3000;" << endl;	
		GMATfile << endl;
		GMATfile << endl;
		++index_body_visited;
	}//end method MGA_DSM_phase::output_GMAT_fueltank_and_thruster

	void MGA_DSM_phase::output_GMAT_create_interphase_control_variables(int& j, int& p, missionoptions& options, std::ofstream& GMATfile)
	{
		//create a variable for the burn index
		GMATfile << "Create Variable BurnIndex_Journey" << j + 1 << "Phase" << p + 1 << ";" << endl;
	}//end method MGA_DSM_phase::output_GMAT_create_interphase_control_variables

	void MGA_DSM_phase::output_GMAT_inter_phase_control_initial_guess(int& j, int& p, missionoptions& options, std::ofstream& GMATfile)
	{
		//for MGADSM missions we must encode the burn index
		GMATfile << endl;
		GMATfile <<"	%Journey #"<< j + 1 <<", Phase #"<< p + 1 << " Burn Index" << endl;
		GMATfile << "	BurnIndex_Journey" << j + 1 << "Phase" << p + 1 << " = " << (burn_index - 0.01) / 0.98 << endl;
	}//end method MGA_DSM_phase::output_GMAT_inter_phase_control_initial_guess

}