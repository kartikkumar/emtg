/*
 * mission.cpp
 *
 *  Created on: Jul 17, 2012
 *      Author: Jacob
 */

#include "journey.h"
#include "mission.h"
#include "missionoptions.h"
#include "file_utilities.h"
#include "universe.h"
#include "EMTG_math.h"
#include "rk8713M.h"
#include "interpolator.h"

#include "SpiceUsr.h"

#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

namespace EMTG {

mission::mission()  :
	number_of_journeys(0)
{
	// default constructor is never used

}

mission::mission(int* Xouter, missionoptions* options_in, boost::ptr_vector<Astrodynamics::universe>& TheUniverse_in, int thread_ID_assigned, int problem_ID_assigned)
{

	//first make a local copy of the options structure
	this->options = *options_in;

	//make a local copy of the Universe vector
	this->TheUniverse = TheUniverse_in;

	//next, parse the outer-loop decision vector
	this->parse_outer_loop(Xouter);

	//next, create the journeys
	this->number_of_journeys = options.number_of_journeys;

	for (int j = 0; j < number_of_journeys; ++j) 
	{
		this->journeys.push_back(new journey(&options, j, this->TheUniverse[j]));
	}

	//calculate the upper and lower bounds on the decision variables and the constraints
	this->calcbounds();

	//size the local "G storage" vector
	this->G.resize(this->Gdescriptions.size());

	//store the indices of G in the options structure, for later use
	this->options.iGfun = this->iGfun;
	this->options.jGvar = this->jGvar;

	//store the scale ranges for all of the decision variables
	for (size_t entry = 0; entry < this->Xupperbounds.size(); ++entry)
		this->options.X_scale_ranges.push_back(this->Xupperbounds[entry] - this->Xlowerbounds[entry]);

	//find the time scale factor
	this->max_TU = 0.0;
	for (int j = 0; j < this->options.number_of_journeys; ++j)
		this->max_TU = max(this->max_TU, this->TheUniverse[j].TU);
}

mission::mission(missionoptions* options_in, boost::ptr_vector<Astrodynamics::universe>& TheUniverse_in)
{
	//first make a local copy of the options structure
	this->options = *options_in;

	//make a local copy of the Universe vector
	this->TheUniverse = TheUniverse_in;

	//next, create the journeys
	this->number_of_journeys = options.number_of_journeys;

	for (int j = 0; j < number_of_journeys; ++j) 
	{
		this->journeys.push_back(new journey(&options, j, this->TheUniverse[j]));
	}

	//calculate the upper and lower bounds on the decision variables and the constraints
	this->calcbounds();

	//size the local "G storage" vector
	this->G.resize(this->Gdescriptions.size());

	//store the indices of G in the options structure, for later use
	this->options.iGfun = this->iGfun;
	this->options.jGvar = this->jGvar;

	//store the scale ranges for all of the decision variables
	for (size_t entry = 0; entry < this->Xupperbounds.size(); ++entry)
		this->options.X_scale_ranges.push_back(this->Xupperbounds[entry] - this->Xlowerbounds[entry]);

	//find the time scale factor
	this->max_TU = 0.0;
	for (int j = 0; j < this->options.number_of_journeys; ++j)
		this->max_TU = max(this->max_TU, this->TheUniverse[j].TU);
}

mission::~mission()
{	
}

//outer-loop parse function
//return 0 for success, 1 for failure
int mission::parse_outer_loop(int* Xouter)
{
	int errcode = 0;


	//parse the outer loop vector
	options.total_number_of_phases = 0;


	//if we have specified the mission type (all MGA, MGA-DSM, or MGA-LT) then it only takes one decision variable to encode a phase
	//if we have NOT specified the mission type, it takes two decision variables to encode a phase
	int phase_encode_length = (options.mission_type > 5 ? 2 : 1);

	//loop through the journeys and figure out how many phases and what type they are
	for (int j = 0; j < options.number_of_journeys; ++j)
	{
		//first, encode in the description the name of the current journey's central body
		if (j > 0) //if not the first journey, insert an underscore
			options.description.append("_");
		options.description.append(TheUniverse[j].central_body_name + "(");
		options.number_of_phases[j] = 1;
		vector<int> temp_sequence;
		vector<int> temp_phase_type;

		temp_sequence.push_back(options.destination_list[j][0]);

		//options.description.push_back(planetcodes[options.destination_list[j]-1]);
		switch (temp_sequence[0])
		{
			case -3: //begin at point on orbit
				{
					options.description.append("o");
					break;
				}
			case -2: //begin at fixed point
				{
					options.description.append("p");
					break;
				}
			case -1: //begin at SOI
				{
					options.description.append("s");
					break;
				}
			default:
				options.description.append(TheUniverse[j].bodies[temp_sequence[0]-1].short_name);
		}
		

		for (int p = 0; p < options.max_phases_per_journey; ++p)
		{
			if (phase_encode_length == 2 && Xouter[j*(2*options.max_phases_per_journey + 1) + p] > 0 && Xouter[j*(2*options.max_phases_per_journey + 1) + p] < (TheUniverse[j].size_of_flyby_menu/2) + 1) //this is a legitimate flyby
			{
				//update the sequence array with a code for the next body
				temp_sequence.push_back(TheUniverse[j].bodies[TheUniverse[j].flyby_menu[Xouter[j * (2*options.max_phases_per_journey + 1) + p] - 1]].body_code);

				//if the outer-loop is choosing phase type, then extract the phase type
				temp_phase_type.push_back(Xouter[(2*j + 1) * options.max_phases_per_journey + j + p] - 1);

				//update the mission description
				options.description.append(TheUniverse[j].bodies[TheUniverse[j].flyby_menu[temp_sequence[temp_sequence.size() - 1]] - 1].short_name);

				//keep track of the number of phases
				++options.number_of_phases[j];
			}
			else if (phase_encode_length == 1 && Xouter[j * options.max_phases_per_journey + p] > 0 && Xouter[j * options.max_phases_per_journey + p] < (TheUniverse[j].size_of_flyby_menu/2) + 1) //this is a legitimate flyby
			{
				//update the sequence array with a code for the next body
				temp_sequence.push_back(TheUniverse[j].bodies[TheUniverse[j].flyby_menu[Xouter[j * options.max_phases_per_journey + p] - 1]].body_code);

				//otherwise use the global phase type specified by the user
				temp_phase_type.push_back(options.mission_type);

				//update the mission description
				options.description.append(TheUniverse[j].bodies[temp_sequence[temp_sequence.size() - 1] - 1].short_name);

				//keep track of the number of phases
				++options.number_of_phases[j];
			}
		}

		//encode the last phase of the journey
		temp_sequence.push_back(options.destination_list[j][1]);

		//if the outer-loop is choosing phase type, then extract the phase type for the last phase
		//otherwise use the global phase type specified by the user
		if (phase_encode_length == 2)
			temp_phase_type.push_back(Xouter[(2*j + 2) * options.max_phases_per_journey + j]);
		else
			temp_phase_type.push_back(options.mission_type);

		//update the mission description
		switch (temp_sequence[temp_sequence.size() - 1])
		{
			case -3: //begin at point on orbit
				{
					options.description.append("o");
					break;
				}
			case -2: //begin at fixed point
				{
					options.description.append("p");
					break;
				}
			case -1: //begin at SOI
				{
					options.description.append("s");
					break;
				}
			default:
				options.description.append(TheUniverse[j].bodies[temp_sequence[temp_sequence.size() - 1]-1].short_name);
		}
		options.description.append(")");

		//store this information in the options structure
		options.sequence.push_back(temp_sequence);
		options.phase_type.push_back(temp_phase_type);

		//track the total number of phases
		options.total_number_of_phases += options.number_of_phases[j];
	}

	return errcode;
}

//bounds calculation function
//return 0 for success, 1 for failure
int mission::calcbounds()
{
	//bounds on the objective function
	Flowerbounds.push_back(-math::LARGE);
	Fupperbounds.push_back(math::LARGE);
	Fdescriptions.push_back("objective function");

	//call the calcbounds() function for each journey
	for (int j = 0; j < number_of_journeys; ++j)
		journeys[j].calcbounds(&Xupperbounds, &Xlowerbounds, &Fupperbounds, &Flowerbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &A, &Adescriptions, &Gdescriptions, &synodic_periods, j, TheUniverse[j], &options);

	//one final constraint, if applicable, for total time bounds
	if (options.global_timebounded)
	{
		Flowerbounds.push_back(options.total_flight_time_bounds[0] / options.total_flight_time_bounds[1] - 1);
		Fupperbounds.push_back(0.0);
		Fdescriptions.push_back("Mission flight time bounds");

		//Generate the Jacobian entries for the mission flight time constraint
		//note:
		//1. ALL time variables present in the decision vector affect the mission flight time constraint
		//2. The mission time constraint is linear with respect to the time variables
		for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
		{
			if (Xdescriptions[entry].find("time") < 1024)
			{
				iGfun.push_back(Fdescriptions.size() - 1);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of mission flight time constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				//A.push_back( ((Xupperbounds)[entry] - (Xlowerbounds)[entry]) / options.total_flight_time_bounds[1]);
				timeconstraints_X_scale_ranges.push_back((Xupperbounds)[entry] - (Xlowerbounds)[entry]);
				timeconstraints_G_indices.push_back(iGfun.size() - 1);		
			}
		}

		//locate dependencies on spirals anywhere in the mission
		bool spirals_exist = false;
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			if (options.journey_departure_type[j] == 5 || options.journey_arrival_type[j] == 7)
			{
				spirals_exist = true;
				break;
			}
		}
		if (spirals_exist)
		{
			//call the usual spiral dependency functions
			this->find_dependencies_due_to_escape_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, Fdescriptions.size() - 1);
			this->find_dependencies_due_to_capture_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, Fdescriptions.size() - 1);
			
			//also, if there is a spiral anywhere in the mission then the flight time constraint has a dependency on launch date
			//this is because time drives spiral propellant use, which in turn drives spiral time
			iGfun.push_back(Fdescriptions.size() - 1);
			jGvar.push_back(0);
			stringstream EntryNameStream;
			EntryNameStream << "Derivative of flight time constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << 0 << "]: " << Xdescriptions[0];
			Gdescriptions.push_back(EntryNameStream.str());
			derivative_of_flight_time_with_respect_to_launch_date_G_index = iGfun.size() - 1;
		}
	}

	//Generate the Jacobian entries for the objective function
	//For MGA-LT, FBLT, FBLT-S, MGA-NDSM the objective function is linear and dependent only on the arrival mass of the last phase
	//For MGA and MGA-DSM, the objective function is noninear and dependent on ALL of the variables
	if (options.mission_type >= 2 && options.objective_type == 2) //MGA-LT, FBLT, FBLT-S, MGA-NDSM maximum final mass optimization
	{
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("arrival mass") < 1024)
			{
				iAfun.push_back(0);
				jAvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Adescriptions.push_back(EntryNameStream.str());
				A.push_back(-1.0 / (options.maximum_mass + journeys[options.number_of_journeys - 1].phases[options.number_of_phases[options.number_of_journeys - 1] - 1].current_mass_increment) * (Xupperbounds[entry] - Xlowerbounds[entry]) + Xlowerbounds[entry]);
				break;
			}
		}

		//dependencies because of a capture spiral in the final journey
		this->find_dependencies_due_to_capture_spiral_in_final_journey(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);

		if (options.journey_arrival_type[options.number_of_journeys - 1] == 1) //chemical rendezvous
		{
			int count = 0;
			for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
			{
				if (Xdescriptions[entry].find("Terminal velocity") < 1024)
				{
					iGfun.push_back(0);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					++count;
				}

				if (count >= 2)
					break;
			}
		}
	}
	else if (options.objective_type == 1) //minimum flight time
	{
        double TU = TheUniverse[options.number_of_journeys - 1].TU;
        for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("flight time") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
                Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
			}
		}

		//locate dependencies on spirals anywhere in the mission
		this->find_dependencies_due_to_escape_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);
		this->find_dependencies_due_to_capture_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);
	}
	else if (options.mission_type >= 2 && options.objective_type == 3) //MGA-LT, FBLT, FBLT-S, MGA-NDSM GTOC1 optimization
	{
		//derivative with respect to arrival mass
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("arrival mass") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				break;
			}
		}

		//derivative with respect to all epochs and times
		for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
		{
			if (Xdescriptions[entry].find("time") < 1024 || Xdescriptions[entry].find("epoch") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
			}
		}

		//derivative with respect to terminal velocity
		int foundcount = 0;
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("terminal velocity") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
				++foundcount;
			}
			if (foundcount >= 3)
				break;
		}
	}
	else if (options.objective_type == 4) //launch as a late as possible
	{
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("launch epoch") < 1024)
			{
				iAfun.push_back(0);
				jAvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
                Adescriptions.push_back(EntryNameStream.str());
				A.push_back(-(Xupperbounds[entry] - Xlowerbounds[entry]) / options.launch_window_open_date);
				break;
			}
		}
	}
    else if (options.objective_type == 5) //launch as a early as possible
	{
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("launch epoch") < 1024)
			{
				iAfun.push_back(0);
				jAvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
                Adescriptions.push_back(EntryNameStream.str());
				A.push_back((Xupperbounds[entry] - Xlowerbounds[entry]) / options.launch_window_open_date);
				break;
			}
		}
	}
    else if (options.objective_type == 6) //maximum orbit energy
	{
        stringstream last_phase_prefixstream;
        last_phase_prefixstream << "j" << options.number_of_journeys - 1 << "p" << options.number_of_phases[options.number_of_journeys - 1];
        string last_phase_prefix = last_phase_prefixstream.str();

		//Jacobian entry for zero energy condition
		//this is a nonlinear constraint dependent on the terminal velocity vector
		//note that it is NOT dependent on position because these phases always end at the SOI and SOI distance is constant
		for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
		{
			if (Xdescriptions[entry].find(last_phase_prefix.c_str()) < 1024 && (Xdescriptions[entry].find("terminal velocity") < 1024 || Xdescriptions[entry].find("right boundary") < 1024))
			{
				iGfun.push_back(Fdescriptions.size() - 1);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
			}
		}
	}
    else if (options.objective_type == 8) //arrive as early as possible
	{
        double TU = TheUniverse[options.number_of_journeys - 1].TU;
        for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("flight time") < 1024 || Xdescriptions[entry].find("launch epoch") < 1024)
			{
				iAfun.push_back(0);
				jAvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
                Adescriptions.push_back(EntryNameStream.str());
				A.push_back(-(Xupperbounds[entry] - Xlowerbounds[entry]) / TU);
			}
		}

		//locate dependencies on spirals anywhere in the mission
		this->find_dependencies_due_to_escape_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);
		this->find_dependencies_due_to_capture_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);
	}
    else if (options.objective_type == 9) //arrive as late as possible
	{
        double TU = TheUniverse[options.number_of_journeys - 1].TU;
        for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("flight time") < 1024 || Xdescriptions[entry].find("launch epoch") < 1024)
			{
				iAfun.push_back(0);
				jAvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
                Adescriptions.push_back(EntryNameStream.str());
				A.push_back((Xupperbounds[entry] - Xlowerbounds[entry]) / TU);
			}
		}

		//locate dependencies on spirals anywhere in the mission
		this->find_dependencies_due_to_escape_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);
		this->find_dependencies_due_to_capture_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);
	}
	else if (options.mission_type >= 2 && options.objective_type == 12) //MGA-LT, FBLT, FBLT-S, MGA-NDSM maximum kinetic-energy
	{
		//derivative with respect to arrival mass
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("arrival mass") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
				break;
			}
		}

		//derivative with respect to terminal velocity
		int foundcount = 0;
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("terminal velocity") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
				++foundcount;
			}
			if (foundcount >= 3)
				break;
		}
	}
	else if (options.mission_type < 2)//MGA, MGA-DSM
	{
		//the objective function of an MGA-DSM mission has derivatives with respect to everything
		for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
		{
			iGfun.push_back(0);
			jGvar.push_back(entry);
			stringstream EntryNameStream;
			EntryNameStream << "Derivative of objectivefunction F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
			Gdescriptions.push_back(EntryNameStream.str());
		}
	}
	else if (options.objective_type == 7 && options.mission_type > 1) //minimize launch mass
	{
        //derivative with respect to initial v-infinity
        //derivative with respect to all epochs and times
		for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
		{
			if (Xdescriptions[entry].find("magnitude of outgoing velocity asymptote") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
				break;
			}
		}

        //if applicable, derivative with respect to initial mass increment for the first journey
        if (options.journey_starting_mass_increment[0] > 0.0 && options.journey_variable_mass_increment[0])
		{
			for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
			{
				if (Xdescriptions[entry].find("journey initial mass scale factor") < 1024)
				{
					iGfun.push_back(0);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					objectivefunction_G_indices.push_back(iGfun.size() - 1);
					objectivefunction_X_indices.push_back(entry);
					objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
                    break;
				}
			} 
		}

        //if we allow initial mass to vary, derivative with respect to initial mass multiplier
        if (options.allow_initial_mass_to_vary)
		{
			for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
			{
				if (Xdescriptions[entry].find("initial mass multiplier (0-1)") < 1024)
				{
					iGfun.push_back(0);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					objectivefunction_G_indices.push_back(iGfun.size() - 1);
					objectivefunction_X_indices.push_back(entry);
					objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
                    break;
				}
			} 
		}         
	}
	else if (options.mission_type > 1 && (options.objective_type == 10 || options.objective_type == 11)) //minimum propellant or maximum dry/wet ratio
	{
		//the minimum propellant and maximum dry/wet ratio objective functions have the same sparsity pattern
		//these objective functions has derivatives with respect to the initial mass, the initial mass multiplier, and the final mass
		//it therefore has the same derivative structure as the combination of the maximum final mass objective function and the
		//minimum launch mass objective function
		//this objective function is nonlinear
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("arrival mass") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				G.push_back(-1.0 / (options.maximum_mass + journeys[options.number_of_journeys - 1].phases[options.number_of_phases[options.number_of_journeys - 1] - 1].current_mass_increment) * (Xupperbounds[entry] - Xlowerbounds[entry]) + Xlowerbounds[entry]);
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
				break;
			}
		}

		//derivative with respect to initial v-infinity
		for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
		{
			if (Xdescriptions[entry].find("magnitude of outgoing velocity asymptote") < 1024)
			{
				iGfun.push_back(0);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				objectivefunction_G_indices.push_back(iGfun.size() - 1);
				objectivefunction_X_indices.push_back(entry);
				objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
				break;
			}
		}

        //if applicable, derivative with respect to initial mass increment for the first journey
        if (options.journey_starting_mass_increment[0] > 0.0 && options.journey_variable_mass_increment[0])
		{
			for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
			{
				if (Xdescriptions[entry].find("journey initial mass scale factor") < 1024)
				{
					iGfun.push_back(0);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					objectivefunction_G_indices.push_back(iGfun.size() - 1);
					objectivefunction_X_indices.push_back(entry);
					objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
                    break;
				}
			} 
		}

        //if we allow initial mass to vary, derivative with respect to initial mass multiplier
        if (options.allow_initial_mass_to_vary)
		{
			for (size_t entry = 0; entry < Xdescriptions.size(); ++entry)
			{
				if (Xdescriptions[entry].find("initial mass multiplier (0-1)") < 1024)
				{
					iGfun.push_back(0);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					objectivefunction_G_indices.push_back(iGfun.size() - 1);
					objectivefunction_X_indices.push_back(entry);
					objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
                    break;
				}
			} 
		}

		if (options.journey_arrival_type[options.number_of_journeys - 1] == 1) //chemical rendezvous
		{
			int count = 0;
			for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
			{
				if (Xdescriptions[entry].find("Terminal velocity") < 1024)
				{
					iGfun.push_back(0);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					objectivefunction_G_indices.push_back(iGfun.size() - 1);
					objectivefunction_X_indices.push_back(entry);
					objectivefunction_X_scale_ranges.push_back(Xupperbounds[entry] - Xlowerbounds[entry]);
					++count;
				}

				if (count >= 2)
					break;
			}
		}

		//dependencies because of a capture spiral in the final journey
		this->find_dependencies_due_to_capture_spiral_in_final_journey(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, 0);
	}
	else if (this->options.objective_type == 13)
	{
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("engine input power (kW)") < 1024)
			{
				iAfun.push_back(0);
				jAvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of objective function F[0] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Adescriptions.push_back(EntryNameStream.str());
				A.push_back(1.0);
				break;
			}
		}
		
	}

	//mission dry mass constraint
	if (options.minimum_dry_mass > 0)
	{
		Flowerbounds.push_back(-math::LARGE);
		Fupperbounds.push_back(0);
		Fdescriptions.push_back("mission dry mass constraint");

		//derivative with respect to final mass
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("arrival mass") < 1024)
			{
				iGfun.push_back(Fdescriptions.size() - 1);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of dry mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				dry_mass_constraint_G_indices.push_back(iGfun.size() - 1);
				dry_mass_constraint_X_indices.push_back(entry);
				break;
			}
		}

		//derivative with respect to v-infinity
		if (!(options.LV_type == 0))
		{
			for (int entry = 0; entry < Xdescriptions.size() - 1; ++entry)
			{
				if (Xdescriptions[entry].find("magnitude of outgoing velocity asymptote") < 1024)
				{
					iGfun.push_back(Fdescriptions.size() - 1);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of dry mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					dry_mass_constraint_G_indices.push_back(iGfun.size() - 1);
					dry_mass_constraint_X_indices.push_back(entry);
					break;
				}
			}
		}

		//derivative with respect to initial mass multiplier
		if (options.allow_initial_mass_to_vary)
		{
			for (int entry = 0; entry < Xdescriptions.size() - 1; ++entry)
			{
				if (Xdescriptions[entry].find("initial mass multiplier") < 1024)
				{
					iGfun.push_back(Fdescriptions.size() - 1);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of dry mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					dry_mass_constraint_G_indices.push_back(iGfun.size() - 1);
					dry_mass_constraint_X_indices.push_back(entry);
					break;
				}
			}
		}

		//derivative entry with respect to final journey mass increment ratio
		if (options.journey_variable_mass_increment[options.number_of_journeys-1])
		{
			for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
			{
				if (Xdescriptions[entry].find("journey initial mass scale factor") < 1024)
				{
					iGfun.push_back(Fdescriptions.size() - 1);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of dry mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					dry_mass_constraint_G_indices.push_back(iGfun.size() - 1);
					dry_mass_constraint_X_indices.push_back(entry);
					break;
				}
			}
		}


		//dependencies due to spirals
		this->find_dependencies_due_to_escape_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, Fdescriptions.size() - 1);
		this->find_dependencies_due_to_capture_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, Fdescriptions.size() - 1);
	}

	if (options.enable_maximum_propellant_mass_constraint)
	{
		Flowerbounds.push_back(0.0);
		Fupperbounds.push_back(1.0);
		Fdescriptions.push_back("mission propellant mass constraint");

		//derivative with respect to final mass
		for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
		{
			if (Xdescriptions[entry].find("arrival mass") < 1024)
			{
				iGfun.push_back(Fdescriptions.size() - 1);
				jGvar.push_back(entry);
				stringstream EntryNameStream;
				EntryNameStream << "Derivative of propellant mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
				Gdescriptions.push_back(EntryNameStream.str());
				propellant_mass_constraint_G_indices.push_back(iGfun.size() - 1);
				propellant_mass_constraint_X_indices.push_back(entry);
				break;
			}
		}

		//derivative with respect to v-infinity
		if (!(options.LV_type == 0))
		{
			for (int entry = 0; entry < Xdescriptions.size() - 1; ++entry)
			{
				if (Xdescriptions[entry].find("magnitude of outgoing velocity asymptote") < 1024)
				{
					iGfun.push_back(Fdescriptions.size() - 1);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of propellant mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					propellant_mass_constraint_G_indices.push_back(iGfun.size() - 1);
					propellant_mass_constraint_X_indices.push_back(entry);
					break;
				}
			}
		}

		//derivative with respect to initial mass multiplier
		if (options.allow_initial_mass_to_vary)
		{
			for (int entry = 0; entry < Xdescriptions.size() - 1; ++entry)
			{
				if (Xdescriptions[entry].find("initial mass multiplier") < 1024)
				{
					iGfun.push_back(Fdescriptions.size() - 1);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of propellant mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					propellant_mass_constraint_G_indices.push_back(iGfun.size() - 1);
					propellant_mass_constraint_X_indices.push_back(entry);
					break;
				}
			}
		}

		//derivative entry with respect to final journey mass increment ratio
		if (options.journey_variable_mass_increment[options.number_of_journeys - 1])
		{
			for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
			{
				if (Xdescriptions[entry].find("journey initial mass scale factor") < 1024)
				{
					iGfun.push_back(Fdescriptions.size() - 1);
					jGvar.push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of propellant mass constraint F[" << Fdescriptions.size() - 1 << "] with respect to X[" << entry << "]: " << Xdescriptions[entry];
					Gdescriptions.push_back(EntryNameStream.str());
					propellant_mass_constraint_G_indices.push_back(iGfun.size() - 1);
					propellant_mass_constraint_X_indices.push_back(entry);
					break;
				}
			}
		}


		//dependencies due to spirals
		this->find_dependencies_due_to_escape_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, Fdescriptions.size() - 1);
		this->find_dependencies_due_to_capture_spiral(&Xupperbounds, &Xlowerbounds, &Flowerbounds, &Fupperbounds, &Xdescriptions, &Fdescriptions, &iAfun, &jAvar, &iGfun, &jGvar, &Adescriptions, &Gdescriptions, &options, Fdescriptions.size() - 1);
	}

	total_number_of_NLP_parameters = Xupperbounds.size();
	total_number_of_constraints = Fupperbounds.size();

	X.resize(total_number_of_NLP_parameters);
	X0.resize(total_number_of_NLP_parameters);
	F.resize(total_number_of_constraints);
	return 0;
}


//evaluate function
//return 0 if successful, 1 if failure
int mission::evaluate(double* X, double* F, double* G, int needG, const vector<int>& iGfun, const vector<int>& jGvar)
{
	int Xindex = 0;
	int Findex = 1; //F[0] is reserved for the objective function
	int Gindex = 0;
	this->current_deltaV = 0.0;
	int errcode = 0;

	//process all of the journeys
	for (int j = 0; j < number_of_journeys; ++j)
	{
		if (!(errcode))
		{
			errcode = journeys[j].evaluate(X, &Xindex, F, &Findex, G, &Gindex, needG, j, &current_epoch, current_state, &current_deltaV, TheUniverse[j], &options);
		}
		else
		{
			cout << "error " << errcode << endl;
			break;
		}
	}
	
	// evaluate global time bounds
	if (options.global_timebounded)
	{
		F[Findex] = (current_epoch - X[0]) / options.total_flight_time_bounds[1] - 1.0;
		++Findex;

		if (options.derivative_type > 0 && needG)
		{
			for (size_t entry = 0; entry < timeconstraints_G_indices.size(); ++entry)
			{
				G[timeconstraints_G_indices[entry]] = timeconstraints_X_scale_ranges[entry] / options.total_flight_time_bounds[1];
			}
		}
	}

	//create a pointer to the first phase
	EMTG::journey* FirstJourney = &journeys[0];
	EMTG::phase* FirstPhase = &FirstJourney->phases[0];

	//create a pointer to the final phase
	EMTG::journey* FinalJourney = &journeys[options.number_of_journeys - 1];
	EMTG::phase* FinalPhase = &FinalJourney->phases[options.number_of_phases[options.number_of_journeys - 1] - 1];

	//evaluate, if applicable, minimum dry mass bound
	if (options.minimum_dry_mass > 0)
	{
		//compute the system and spacecraft masses at the end of the modeled mission
		double final_system_mass = current_state[6];
		double final_spacecraft_mass = final_system_mass - FinalJourney->phases[0].journey_initial_mass_increment_scale_factor * FinalPhase->current_mass_increment;

		//apply the post-mission delta-v to determine the remaining mass of the spacecraft
		double expfun = exp(-1000 * options.post_mission_delta_v / (options.g0 * options.post_mission_Isp));
		double system_mass_after_post_mission_delta_v = final_system_mass * expfun;
		double spacecraft_mass_after_post_mission_delta_v = final_spacecraft_mass - (final_system_mass - system_mass_after_post_mission_delta_v);

		//apply propellant margin
		double initial_spacecraft_mass;
		if (options.journey_departure_type[0] == 5) //if the first journey started with a spiral, get the state before the spiral
			initial_spacecraft_mass = FirstJourney->phases[0].spiral_escape_state_before_spiral[6];
		else
			initial_spacecraft_mass = FirstPhase->state_at_beginning_of_phase[6];
		double propellant_margin_kg = options.propellant_margin * (initial_spacecraft_mass - spacecraft_mass_after_post_mission_delta_v);

		dry_mass = spacecraft_mass_after_post_mission_delta_v - propellant_margin_kg;


		F[Findex] = -dry_mass / options.minimum_dry_mass + 1.0;

		if (options.derivative_type > 0 && needG)
		{
			int whichderiv = 0;
			//derivative with respect to arrival mass
			G[dry_mass_constraint_G_indices[whichderiv]] = -(options.maximum_mass + FinalPhase->current_mass_increment) * expfun * (options.propellant_margin + 1) / options.minimum_dry_mass;
			++whichderiv;

			//derivative with respect to v-infinity
			if (!(options.LV_type == 0))
			{
				G[dry_mass_constraint_G_indices[whichderiv]] = (options.journey_initial_impulse_bounds[0][1] - options.journey_initial_impulse_bounds[0][0]) * FirstPhase->dmdvinf * (FirstPhase->mission_initial_mass_multiplier * options.propellant_margin) / options.minimum_dry_mass;
				++whichderiv;
			}


			//derivative with respect to initial mass scale factor
			if (options.allow_initial_mass_to_vary)
			{
				//the 0.8 is because initial mass scale factor varies in [0.2, 1.0]
				G[dry_mass_constraint_G_indices[whichderiv]] = (0.8 * FirstPhase->unscaled_phase_initial_mass * options.propellant_margin) / options.minimum_dry_mass;
				++whichderiv;
			}

			//derivative with respect to final journey mass increment ratio
			if (options.journey_variable_mass_increment[options.number_of_journeys - 1])
				G[dry_mass_constraint_G_indices[whichderiv]] = FinalPhase->current_mass_increment * (options.propellant_margin + 1) / options.minimum_dry_mass;
		}
	}
	else dry_mass = current_state[6];

	//evaluate, if applicable, the propellant mass constraint
	if (options.enable_maximum_propellant_mass_constraint)
	{
		//compute the system and spacecraft masses at the end of the modeled mission
		double final_system_mass = current_state[6];
		double final_spacecraft_mass = final_system_mass - FinalJourney->phases[0].journey_initial_mass_increment_scale_factor * FinalPhase->current_mass_increment;

		//apply the post-mission delta-v to determine the remaining mass of the spacecraft
		double expfun = exp(-1000 * options.post_mission_delta_v / (options.g0 * options.post_mission_Isp));
		double system_mass_after_post_mission_delta_v = final_system_mass * expfun;
		double spacecraft_mass_after_post_mission_delta_v = final_spacecraft_mass - (final_system_mass - system_mass_after_post_mission_delta_v);

		//apply propellant margin
		double initial_spacecraft_mass;
		if (options.journey_departure_type[0] == 5) //if the first journey started with a spiral, get the state before the spiral
			initial_spacecraft_mass = FirstJourney->phases[0].spiral_escape_state_before_spiral[6];
		else
			initial_spacecraft_mass = FirstPhase->state_at_beginning_of_phase[6];
		double propellant_mass_kg = initial_spacecraft_mass - spacecraft_mass_after_post_mission_delta_v;
		double propellant_margin_kg = options.propellant_margin * propellant_mass_kg;

		total_propellant_mass = propellant_mass_kg + propellant_margin_kg;


		F[Findex] = -total_propellant_mass / options.maximum_propellant_mass + 1.0;
		/*
		if (options.derivative_type > 0 && needG)
		{
			int whichderiv = 0;
			//derivative with respect to arrival mass
			G[propellant_mass_constraint_G_indices[whichderiv]] = -(options.maximum_mass + FinalPhase->current_mass_increment) * expfun * (options.propellant_margin + 1) / options.minimum_dry_mass;
			++whichderiv;

			//derivative with respect to v-infinity
			if (!(options.LV_type == 0))
			{
				G[propellant_mass_constraint_G_indices[whichderiv]] = (options.journey_initial_impulse_bounds[0][1] - options.journey_initial_impulse_bounds[0][0]) * FirstPhase->dmdvinf * (FirstPhase->mission_initial_mass_multiplier * options.propellant_margin) / options.minimum_dry_mass;
				++whichderiv;
			}


			//derivative with respect to initial mass scale factor
			if (options.allow_initial_mass_to_vary)
			{
				//the 0.8 is because initial mass scale factor varies in [0.2, 1.0]
				G[propellant_mass_constraint_G_indices[whichderiv]] = (0.8 * FirstPhase->unscaled_phase_initial_mass * options.propellant_margin) / options.minimum_dry_mass;
				++whichderiv;
			}

			//derivative with respect to final journey mass increment ratio
			if (options.journey_variable_mass_increment[options.number_of_journeys - 1])
				G[propellant_mass_constraint_G_indices[whichderiv]] = FinalPhase->current_mass_increment * (options.propellant_margin + 1) / options.minimum_dry_mass;
		}*/
	}
	else total_propellant_mass = FirstJourney->phases[0].state_at_beginning_of_phase[0] - current_state[6];

	switch (options.objective_type)
	{
	case 0: //minimum deltaV
		{
			F[0] = current_deltaV;
			break;
		}
	case 1:
		{ //minimum flight time
			double TU = TheUniverse[options.number_of_journeys - 1].TU;
			F[0] = (current_epoch - X[0]) / TU;

			if (this->options.derivative_type > 0)
			{
				for (int whichderiv = 0; whichderiv < this->objectivefunction_G_indices.size() - 1; ++ whichderiv)
				{
					G[this->objectivefunction_G_indices[whichderiv]] = objectivefunction_X_scale_ranges[whichderiv] / TU;
				}
			}
			break;
		}
	case 2:
		{ //maximum final mass
			if (options.mission_type > 1)
				F[0] = -current_state[6] / (options.maximum_mass + FinalPhase->current_mass_increment);
			else
				F[0] = -current_state[6];
			break;
		}
	case 3:
		{ //GTOC 1 asteroid deflection function
			if (FinalPhase->boundary2_location_code < 0)
			{
				cout << "GTOC 1 deflection function does not work unless the last boundary point in the mission is a body!" << endl;
				throw 1711;
			}
			//locate the final body
			double final_boundary_state[9];
			FinalPhase->locate_boundary_point(FinalPhase->boundary2_location_code, options.journey_arrival_type[options.number_of_journeys - 1], needG, &TheUniverse[options.number_of_journeys - 1], final_boundary_state, current_state+3, current_epoch, X, &Xindex, F, &Findex, G, &Gindex, true, options.number_of_journeys - 1, options.number_of_phases[options.number_of_journeys - 1] - 1, &options);
			math::Matrix<double> Vsc_final(3, 1, FinalPhase->state_at_end_of_phase + 3);
			math::Matrix<double> Vbody_final(3, 1, final_boundary_state + 3);

			Vsc_final -= Vbody_final;

			double LU = TheUniverse[options.number_of_journeys - 1].LU;
			double TU = TheUniverse[options.number_of_journeys - 1].TU;
			double Vsc_dot_Vbody = Vsc_final.dot(Vbody_final) * (TU * TU / LU / LU);
			
			F[0] = -current_state[6] / (options.maximum_mass + FinalPhase->current_mass_increment) * Vsc_dot_Vbody;

			//derivatives of the objective function
			if (options.derivative_type && needG)
			{
				//derivative with respect to mass
				G[objectivefunction_G_indices[0]] = -1.0 / (options.maximum_mass + FinalPhase->current_mass_increment) * Vsc_dot_Vbody * (Xupperbounds[objectivefunction_X_indices[0]] - Xlowerbounds[objectivefunction_X_indices[0]]);

				//derivatives with respect to time - require finite differencing of the ephemeris, so omitted

				//derivatives with respect to terminal velocity
				G[objectivefunction_G_indices[1]] = -current_state[6] / (options.maximum_mass + FinalPhase->current_mass_increment) * Vbody_final(2) * (TU * TU / LU / LU) * (Xupperbounds[objectivefunction_X_indices[3]] - Xlowerbounds[objectivefunction_X_indices[3]]);
				G[objectivefunction_G_indices[2]] = -current_state[6] / (options.maximum_mass + FinalPhase->current_mass_increment) * Vbody_final(1) * (TU * TU / LU / LU) * (Xupperbounds[objectivefunction_X_indices[2]] - Xlowerbounds[objectivefunction_X_indices[2]]);
				G[objectivefunction_G_indices[3]] = -current_state[6] / (options.maximum_mass + FinalPhase->current_mass_increment) * Vbody_final(0) * (TU * TU / LU / LU) * (Xupperbounds[objectivefunction_X_indices[1]] - Xlowerbounds[objectivefunction_X_indices[1]]);
			}

            break;
		}
	case 4:
		{ //latest possible launch date
			F[0] = -X[0] / options.launch_window_open_date;
            break;
		}
	case 5:
		{
            //earliest possible launch date
            F[0] = X[0] / options.launch_window_open_date;
            break;
		}
	case 6:
		{
            //maximum energy
            //first determine r and v
            double LU = TheUniverse[options.number_of_journeys - 1].LU;
			double TU = TheUniverse[options.number_of_journeys - 1].TU;
            
			double v = math::norm(FinalPhase->state_at_end_of_phase + 3, 3) * TU / LU;
			double r = math::norm(FinalPhase->state_at_end_of_phase, 3) / LU;

			//Energy equation - we want to minimize the negative of energy
			F[0] = -(v*v / 2.0 - 1.0 / r);

            break;
		}
	case 7:
		{
            //minimize launch mass
            F[0] = journeys[0].phases[0].state_at_beginning_of_phase[0] / (options.maximum_mass + FinalPhase->current_mass_increment);

            //**************************************************derivatives
            //at this time analytical derivatives for the minimum launch mass problem are not specified
			//derivatives of the objective function
			if (options.derivative_type && needG)
			{
				//derivative with respect to initial v-infinity
				//G[objectivefunction_G_indices[0]] = objectivefunction_X_scale_ranges[0] * journeys[0].phases[0].dmdvinf / (options.maximum_mass + FinalPhase->current_mass_increment);
            
				//if applicable, derivative with respect to initial mass increment for the first journey
            
				//if we allow initial mass to vary, derivative with respect to initial mass multiplier
				if (options.allow_initial_mass_to_vary)
				{
					G[objectivefunction_G_indices[2]] = 0.8 * FirstPhase->unscaled_phase_initial_mass / (options.maximum_mass + FinalPhase->current_mass_increment);
				}
			}

            break;
		}
    case 8:
	    {
            //earliest possible arrival
            double TU = TheUniverse[options.number_of_journeys - 1].TU;
            F[0] = current_epoch / TU;
            break;
	    }
    case 9:
	    {
            //latest possible arrival
            double TU = TheUniverse[options.number_of_journeys - 1].TU;
            F[0] = -current_epoch / TU;
            break;
	    }
	case 10:
		{
			//minimum propellant
			F[0] = (FirstPhase->state_at_beginning_of_phase[6] - current_state[6] + FinalPhase->current_mass_increment) / (options.maximum_mass + FinalPhase->current_mass_increment);

			//derivatives
			if (options.derivative_type > 0 && options.mission_type >= 2 && needG)
			{
				int current_derivative = 0;
				double m_added = FinalPhase->current_mass_increment;
				double mf = current_state[6];
				double m0 = FirstPhase->state_at_beginning_of_phase[6];
				double gamma_m0 = options.allow_initial_mass_to_vary ? FirstPhase->mission_initial_mass_multiplier : 1.0;
				double gamma_maddj0 = options.journey_starting_mass_increment[0] > 0.0 && options.journey_variable_mass_increment[0] ? FirstPhase->journey_initial_mass_increment_scale_factor : 1.0;

				//derivative with respect to final mass
				G[objectivefunction_G_indices[current_derivative]] = -1.0 / (options.maximum_mass + m_added) * objectivefunction_X_scale_ranges[current_derivative];
				++current_derivative;

				//derivative with respect to initial v-infinity
				G[objectivefunction_G_indices[current_derivative]] = ( (1.0 - options.LV_margin) * gamma_m0 * gamma_maddj0 * FirstPhase->dmdvinf ) / (options.maximum_mass + m_added) * objectivefunction_X_scale_ranges[current_derivative];
				++current_derivative;

				//if applicable, derivative with respect to initial mass increment for the first journey
				if (options.journey_starting_mass_increment[0] > 0.0 && options.journey_variable_mass_increment[0])
				{
					G[objectivefunction_G_indices[current_derivative]] = ( m0 / gamma_maddj0 ) / (options.maximum_mass + m_added) * objectivefunction_X_scale_ranges[current_derivative];
					++current_derivative;
				}

				//if we allow initial mass to vary, derivative with respect to initial mass multiplier
				if (options.allow_initial_mass_to_vary)
				{
					G[objectivefunction_G_indices[current_derivative]] = ( m0 / gamma_m0 ) / (options.maximum_mass + m_added) * objectivefunction_X_scale_ranges[current_derivative];
					++current_derivative;
				}

				//chemical rendezvous (derivatives not currently specified)
				//if (options.journey_arrival_type[options.number_of_journeys - 1] == 1)
				//{
				//}
			}
			break;
		}
	case 11:
		{
			//maximum dry/wet fraction
			F[0] = (-current_state[6] + FinalPhase->current_mass_increment) / FirstPhase->state_at_beginning_of_phase[6];

			//derivatives
			if (options.derivative_type > 0 && options.mission_type >= 2 && needG)
			{
				int current_derivative = 0;
				double m_added = FinalPhase->current_mass_increment;
				double mf = current_state[6];
				double m0 = FirstPhase->state_at_beginning_of_phase[6];
				double gamma_m0 = options.allow_initial_mass_to_vary ? FirstPhase->mission_initial_mass_multiplier : 1.0;
				double gamma_maddj0 = options.journey_starting_mass_increment[0] > 0.0 && options.journey_variable_mass_increment[0] ? FirstPhase->journey_initial_mass_increment_scale_factor : 1.0;

				//derivative with respect to final mass
				G[objectivefunction_G_indices[current_derivative]] = -1.0 / m0 * objectivefunction_X_scale_ranges[current_derivative];
				++current_derivative;

				//derivative with respect to initial v-infinity
				G[objectivefunction_G_indices[current_derivative]] = (mf - m_added) / (m0*m0) * ( (1.0 - options.LV_margin) * gamma_m0 * gamma_maddj0 * FirstPhase->dmdvinf ) * objectivefunction_X_scale_ranges[current_derivative];
				++current_derivative;

				//if applicable, derivative with respect to initial mass increment for the first journey
				if (options.journey_starting_mass_increment[0] > 0.0 && options.journey_variable_mass_increment[0])
				{
					G[objectivefunction_G_indices[current_derivative]] = ( 1.0 / gamma_maddj0 ) * objectivefunction_X_scale_ranges[current_derivative];
					++current_derivative;
				}

				//if we allow initial mass to vary, derivative with respect to initial mass multiplier
				if (options.allow_initial_mass_to_vary)
				{
					G[objectivefunction_G_indices[current_derivative]] = (mf - m_added) / (m0 * gamma_m0) * objectivefunction_X_scale_ranges[current_derivative];
					++current_derivative;
				}

				//chemical rendezvous (derivatives not currently specified)
				//if (options.journey_arrival_type[options.number_of_journeys - 1] == 1)
				//{
				//}
			}
			break;
		}
	case 12: //maximum kinetic energy impact
		{
			double LU = TheUniverse[options.number_of_journeys - 1].LU;
			double TU = TheUniverse[options.number_of_journeys - 1].TU;
			F[0] = -0.5 * current_state[6] * FinalPhase->C3_arrival / (options.maximum_mass + FinalPhase->current_mass_increment) / (LU * LU / TU / TU);

			if (options.derivative_type > 0 && options.mission_type >= 2 && needG)
			{
				//derivative with respect to arrival mass
				G[objectivefunction_G_indices[0]] = -0.5 * FinalPhase->C3_arrival * objectivefunction_X_scale_ranges[0] / (options.maximum_mass + FinalPhase->current_mass_increment)  / (LU * LU / TU / TU);

				//derivatives with respect to terminal velocity
				G[objectivefunction_G_indices[1]] = -current_state[6] * X[objectivefunction_X_indices[1]] * objectivefunction_X_scale_ranges[1] / (options.maximum_mass + FinalPhase->current_mass_increment) / (LU * LU / TU / TU);
				G[objectivefunction_G_indices[2]] = -current_state[6] * X[objectivefunction_X_indices[2]] * objectivefunction_X_scale_ranges[2] / (options.maximum_mass + FinalPhase->current_mass_increment) / (LU * LU / TU / TU);
				G[objectivefunction_G_indices[3]] = -current_state[6] * X[objectivefunction_X_indices[3]] * objectivefunction_X_scale_ranges[3] / (options.maximum_mass + FinalPhase->current_mass_increment) / (LU * LU / TU / TU);
			}

			break;
		}
	case 13: //minimum power
		F[0] = this->options.power_at_1_AU;
	}
	
	//test for errors
	if (failed_c())//test for SPICE errors
	{
		F[0] = 1.0e+20;
		reset_c();
	}
	else if (!(F[0] == F[0])) //NaN trap
	{
		F[0] = 1.0e+20;
	}

 	return errcode;
}

//output function
//return 0 if successful, 1 if failure
int mission::output()
{
	ofstream outputfile(options.outputfile.c_str(), ios::out | ios::trunc);
	outputfile << "Mission " << options.mission_name << endl;
	outputfile << "Written by EMTG_v8 core program compiled " << __DATE__<< " " << __TIME__ << endl;
	outputfile.close();

	//next, output summary lines describing each event in the mission
	int errcode = 0;
	int eventcount = 1;
	for (int j = 0; j < number_of_journeys; ++j) {
		errcode = journeys[j].output(&options, journeys[0].phases[0].phase_start_epoch, j, TheUniverse[j], &eventcount);
		if (!(errcode == 0))
			return errcode;
	}
	

	outputfile.open(options.outputfile.c_str(), ios::out |ios::app);

	//skip 3 lines
	for (int k = 0; k < 2; ++k)
		outputfile << endl;

	//output total deltaV
	outputfile << "Total deltaV: " << current_deltaV << endl;

	//output dry mass
	outputfile << "Spacecraft dry mass: " << dry_mass << endl;

	//output propellant mass
	outputfile << "Total propellant mass including margin: " << total_propellant_mass << endl;
	
	//output flight time in years
	outputfile << "Flight time (y): " << (current_epoch - journeys[0].phases[0].phase_start_epoch) / 86400.0 /365.25 << endl;

	//objective-function specific information
	if (options.objective_type == 12)
	{
		phase* FinalPhase = &journeys[options.number_of_journeys-1].phases[options.number_of_phases[options.number_of_journeys-1]-1];
		double LU = TheUniverse[options.number_of_journeys - 1].LU;
		double TU = TheUniverse[options.number_of_journeys - 1].TU;
		double KE = 0.5 * current_state[6] * FinalPhase->C3_arrival;
		outputfile << "Kinetic energy at impact: " << KE << " kg*km^2/s^2" << endl;
	}
	else if (options.objective_type == 13)
		outputfile << "BOL power at 1 AU: " << options.power_at_1_AU << " kW" << endl;

	//skip 3 lines
	for (int k = 0; k < 3; ++k)
		outputfile << endl;
	

	//finally, output the decision vector
	//upper bounds
	outputfile << "Xupperbounds";
	for (int k = 0; k < total_number_of_NLP_parameters; ++k)
	{
		if (this->Xdescriptions[k].find("epoch") < 1024 || this->Xdescriptions[k].find("time") < 1024)
		{
			outputfile << " " << this->Xupperbounds[k] / 86400.0;
		}
		else
			outputfile << " " << this->Xupperbounds[k];
	}
	outputfile << endl;
	//decision vector
	outputfile << "Decision Vector:";
	outputfile.unsetf(ios::fixed);
	outputfile.precision(20);
	for (int k = 0; k < total_number_of_NLP_parameters; ++k)
	{
		if (this->Xdescriptions[k].find("epoch") < 1024 || this->Xdescriptions[k].find("time") < 1024)
		{
			outputfile << " " << this->Xopt[k] / 86400.0;
		}
		else
			outputfile << " " << this->Xopt[k];
	}
	outputfile << endl;
	//lower bounds
	outputfile << "Xlowerbounds";
	for (int k = 0; k < total_number_of_NLP_parameters; ++k)
	{
		if (this->Xdescriptions[k].find("epoch") < 1024 || this->Xdescriptions[k].find("time") < 1024)
		{
			outputfile << " " << this->Xlowerbounds[k] / 86400.0;
		}
		else
			outputfile << " " << this->Xlowerbounds[k];
	}
	outputfile << endl;
	//descriptions
	outputfile << "Xdescriptions";
	for (int k = 0; k < total_number_of_NLP_parameters; ++k)
		outputfile << "," << Xdescriptions[k];
	outputfile << endl;

	outputfile << endl;
	outputfile << endl;

	//output the constraints
	//upper bounds
	outputfile << "Fupperbounds";
	for (int k = 0; k < total_number_of_constraints; ++k)
		outputfile << " " << Fupperbounds[k];
	outputfile << endl;
	//constraint values
	outputfile << "Constraint_Vector";
	for (int k = 0; k < total_number_of_constraints; ++k)
		outputfile << " " << F[k];
	outputfile << endl;
	//lower bounds
	outputfile << "Flowerbounds";
	for (int k = 0; k < total_number_of_constraints; ++k)
		outputfile << " " << Flowerbounds[k];
	outputfile << endl;
	//descriptions
	outputfile << "Fdescriptions";
	for (int k = 0; k < total_number_of_constraints; ++k)
		outputfile << "," << Fdescriptions[k];
	outputfile << endl;

	outputfile.close();

	return 0;
}

//function to output the mission tree
//return 0 if successful, 1 if failure
int mission::output_mission_tree(string filename)
{
	vector<string> phase_type_codes;
	phase_type_codes.push_back("MGA");
	phase_type_codes.push_back("MGA-DSM");
	phase_type_codes.push_back("MGA-LT");
	phase_type_codes.push_back("FBLT");
	phase_type_codes.push_back("FBLT-S");
	phase_type_codes.push_back("MGA-NDSM");
    phase_type_codes.push_back("NSLT");
	std::ofstream outputfile(filename.c_str(), ios::trunc);

	if (outputfile.is_open())
	{
		outputfile << "Mission name: " << options.mission_name << endl;
		outputfile << endl;

		for (int j = 0; j < number_of_journeys; ++j)
		{
			outputfile << "  Journey #" << j << " in " << journeys[j].central_body_name << " frame:" << endl;
			for (int p = 0; p < journeys[j].number_of_phases; ++p)
			{
				string boundary1_name;
				string boundary2_name;

				switch (options.sequence[j][p])
				{
					case -2: //begin at fixed point
						{
							boundary1_name = "Periapse of inbound hyperbola";
							break;
						}
					case -1: //begin at SOI
						{
							boundary1_name = "Free or fixed boundary orbit";
							break;
						}
					default:
						boundary1_name = (TheUniverse[j].bodies[options.sequence[j][p]-1].name);
				}

				switch (options.sequence[j][p+1])
				{
					case -1: //begin at SOI
						{
							boundary2_name = "Free or fixed boundary orbit";
							break;
						}
					default:
						boundary2_name = (TheUniverse[j].bodies[options.sequence[j][p+1] -1].name);
				}

				outputfile << "    Phase #" << p << ": Type " << phase_type_codes[options.phase_type[j][p]] << ", From " << boundary1_name << " to " << boundary2_name << endl;
			}
		}

		std::cout << "Mission tree printed to '" << filename << "'" << endl;
		return 0;
	}

	std::cout << "Failure to print mission tree" << endl;
	return 1;
}

//function to create an FBLT-S initial guess from an FBLT input
vector<double> mission::create_initial_guess(vector<double> XFBLT, vector<string>& NewXDescriptions)
{
	vector<double> guessX;

	vector<double> transfer_angles;
	integration::rk8713M integrator(8);

	for (int j = 0; j < options.number_of_journeys; ++j)
	{
		for (int p = 0; p < options.number_of_phases[j]; ++p)
		{
			double angle = acos(math::dot(journeys[j].phases[p].state_at_beginning_of_phase, journeys[j].phases[p].spacecraft_state[0].data(), 3) / (math::norm(journeys[j].phases[p].state_at_beginning_of_phase, 3) * math::norm(journeys[j].phases[p].spacecraft_state[0].data(), 3)));

			for (int step = 1; step < options.num_timesteps; ++step)
			{
				angle += acos(math::dot(journeys[j].phases[p].spacecraft_state[step].data(), journeys[j].phases[p].spacecraft_state[step - 1].data(), 3) / (math::norm(journeys[j].phases[p].spacecraft_state[step].data(), 3) * math::norm(journeys[j].phases[p].spacecraft_state[step - 1].data(), 3)));
			}

			angle += acos(math::dot(journeys[j].phases[p].spacecraft_state[options.num_timesteps - 1].data(), journeys[j].phases[p].state_at_end_of_phase, 3) / (math::norm(journeys[j].phases[p].spacecraft_state[options.num_timesteps - 1].data(), 3) * math::norm(journeys[j].phases[p].state_at_end_of_phase, 3)));

			transfer_angles.push_back(angle);
		}
	}
				
	int skipcount = 0;
	for (size_t entry = 0; entry < Xdescriptions.size() + skipcount; ++entry)
	{
		if (NewXDescriptions[entry].find("final value of Sundman variable s_f") < 1024)
		{
			guessX.push_back(transfer_angles[skipcount]);
			++skipcount;
		}
		else
			guessX.push_back(XFBLT[entry - skipcount]);
	}

	return guessX;
}

//function to interpolate an initial guess
void mission::interpolate(int* Xouter, const vector<double>& initialguess)
{
	//the purpose of this function is to interpolate MGALT, FBLT, and FBLTS decision vectors for use as initial guesses for SNOPT
	//linear interpolation will be used
	//the only values that should be interpolated are the unit control vectors and Isps (in the case of VSI propulsion)

	//Step 1: create a new options file that mimics the original, pre-interpolation problem
	missionoptions* OldOptions = new missionoptions(options);
	OldOptions->num_timesteps = options.initial_guess_num_timesteps;

	//Step 2: create a temporary mission object that mimics the original, pre-interpolation problem
	mission* OldMission = new mission(Xouter, OldOptions, TheUniverse, 0, 0);

	//Step 3: construct the interpolated decision vector
	int XOldEntry = 0;
	vector<double> Xnew;
	vector<double> Xold = options.current_trialX;
	double most_recent_time_increment;
	string most_recent_prefix;

	while (XOldEntry < Xold.size())
	{
		if (!(OldMission->Xdescriptions[XOldEntry].find("step") < 1024))
		{
			//for all non-control parameters, just copy them to the new decision vector unaltered
			Xnew.push_back(Xold[XOldEntry]);

			//we need to keep track of the current time, so that we can use it to index our control tables
			if (strcmp(OldMission->Xdescriptions[XOldEntry].c_str(), "time") < 1024 || strcmp(OldMission->Xdescriptions[XOldEntry].c_str(), "epoch") < 1024)
			{
				most_recent_time_increment = Xold[XOldEntry];
			}

			//increment the entry counter
			++XOldEntry;
		}
		else
		{
			//first get the prefix for the current entry, which tells us which phase and journey we are in
			stringstream prefixstream(OldMission->Xdescriptions[XOldEntry]);
			vector<string> parsed_description;
			boost::split(parsed_description, OldMission->Xdescriptions[XOldEntry], boost::is_any_of(":"));
			most_recent_prefix = parsed_description[0];

			//construct a data table of all control values
			//we will do this by looping through the old vector and extracting the control values
			vector< pair<double, double> > OldUx;
			vector< pair<double, double> > OldUy;
			vector< pair<double, double> > OldUz;
			vector< pair<double, double> > OldIsp;

			//first encode the left hand side of the phase
			OldUx.push_back(std::make_pair(0.0, Xold[XOldEntry]));
			OldUy.push_back(std::make_pair(0.0, Xold[XOldEntry+1]));
			OldUz.push_back(std::make_pair(0.0, Xold[XOldEntry+2]));
			if (OldOptions->engine_type == 4 || OldOptions->engine_type == 12 || OldOptions->engine_type == 13)
			{
				OldIsp.push_back(std::make_pair(0.0, Xold[XOldEntry+3]));
			}

			//then each step
			for (int step = 0; step < OldOptions->num_timesteps; ++step)
			{
				OldUx.push_back(std::make_pair(most_recent_time_increment * (step + 0.5) / OldOptions->num_timesteps, Xold[XOldEntry]));
				OldUy.push_back(std::make_pair(most_recent_time_increment * (step + 0.5) / OldOptions->num_timesteps, Xold[XOldEntry+1]));
				OldUz.push_back(std::make_pair(most_recent_time_increment * (step + 0.5) / OldOptions->num_timesteps, Xold[XOldEntry+2]));
				XOldEntry += 3;

				//for variable specific impulse engine types, we must also encode an Isp
				if (OldOptions->engine_type == 4 || OldOptions->engine_type == 12 || OldOptions->engine_type == 13)
				{
					OldIsp.push_back(std::make_pair(most_recent_time_increment * (step + 0.5) / OldOptions->num_timesteps, Xold[XOldEntry]));
					++XOldEntry;
				}
			}

			//and finally the right-hand side
			OldUx.push_back(std::make_pair(most_recent_time_increment, Xold[XOldEntry-4]));
			OldUy.push_back(std::make_pair(most_recent_time_increment, Xold[XOldEntry-3]));
			OldUz.push_back(std::make_pair(most_recent_time_increment, Xold[XOldEntry-2]));
			if (OldOptions->engine_type == 4 || OldOptions->engine_type == 12 || OldOptions->engine_type == 13)
			{
				OldIsp.push_back(std::make_pair(most_recent_time_increment, Xold[XOldEntry-1]));
			}

			//interpolate the data tables to fill in the new control vector
			math::interpolator ux_interpolator(OldUx);
			math::interpolator uy_interpolator(OldUy);
			math::interpolator uz_interpolator(OldUz);
			math::interpolator Isp_interpolator(OldIsp);

			for (int step = 0; step < options.num_timesteps; ++step)
			{
				double current_time = most_recent_time_increment * (step + 0.5) / options.num_timesteps;
				stringstream stepstream;
				stepstream << step;

				//interpolate and encode the control vector
				Xnew.push_back(ux_interpolator.interpolate(current_time));

				Xnew.push_back(uy_interpolator.interpolate(current_time));

				Xnew.push_back(uz_interpolator.interpolate(current_time));

				//for variable specific impulse engine types, we must also encode an Isp
				if (options.engine_type == 4 || options.engine_type == 13)
				{
					Xnew.push_back(Isp_interpolator.interpolate(current_time));
				}
				else if (options.engine_type == 12)
				{
					Xnew.push_back(Isp_interpolator.interpolate(current_time));
				}
			}
		}
	} //end while loop to create new decision vector

	options.current_trialX = Xnew;
	
	//delete the temporary objects
	delete[] OldOptions;
	delete[] OldMission;

	return;
}

	//functions to convert betwen cartesian and polar initial guesses
	void mission::convert_cartesian_solution_to_polar(const vector<double>& initialguess)
	{
		if (this->options.mission_type == 2 || this->options.mission_type == 3)
		{
			int Xindex = 0;
			vector<double> new_initial_guess;
			while (Xindex < initialguess.size())
			{
				if (!(this->Xdescriptions[Xindex].find("step") < 1024)) //non-control parameters
				{
					new_initial_guess.push_back(initialguess[Xindex]);
					++Xindex;
				}
				else
				{
					//convert from cartesian to polar
					double u_x = initialguess[Xindex];
					double u_y = initialguess[Xindex + 1];
					double u_z = initialguess[Xindex + 2];

					Xindex += 3;

					double throttle = sqrt(u_x*u_x + u_y*u_y + u_z*u_z);
					double theta = atan2(u_y, u_x);
					double phi = asin(u_z / throttle);

					new_initial_guess.push_back(throttle);
					new_initial_guess.push_back(theta / math::TwoPI);
					new_initial_guess.push_back((cos(phi) + 1) / 2.0);

					//if this is a VSI mission we leave the Isp for this time-step unchanged
					if (this->options.engine_type == 4 || this->options.engine_type == 12 || this->options.engine_type == 13)
					{
						new_initial_guess.push_back(initialguess[Xindex]);
						++Xindex;
					}
				}
			}

			this->options.current_trialX = new_initial_guess;
		}
		else
		{
			std::cout << "Converting between cartesian and polar coordinates only works for MGALT and FBLT mission types" << std::endl;
			throw 98;
		}

		return;
	}

	void mission::convert_polar_solution_to_cartesian(const vector<double>& initialguess)
	{
		if (this->options.mission_type == 2 || this->options.mission_type == 3)
		{
			int Xindex = 0;
			vector<double> new_initial_guess;
			while (Xindex < initialguess.size())
			{
				if (!(this->Xdescriptions[Xindex].find("step") < 1024)) //non-control parameters
				{
					new_initial_guess.push_back(initialguess[Xindex]);
					++Xindex;
				}
				else
				{
					//convert from cartesian to polar
					double throttle = initialguess[Xindex];
					double theta = initialguess[Xindex + 1] * math::TwoPI;
					double cosphi = 2 * initialguess[Xindex + 2] - 1;
					double sinphi = sqrt(1.0 - cosphi*cosphi);

					Xindex += 3;

					double u_x = throttle * cos(theta) * cosphi;
					double u_y = throttle * sin(theta) * cosphi;
					double u_z = throttle * sinphi;

					new_initial_guess.push_back(u_x);
					new_initial_guess.push_back(u_y);
					new_initial_guess.push_back(u_z);

					//if this is a VSI mission we leave the Isp for this time-step unchanged
					if (this->options.engine_type == 4 || this->options.engine_type == 12 || this->options.engine_type == 13)
					{
						new_initial_guess.push_back(initialguess[Xindex]);
						++Xindex;
					}
				}
			}

			this->options.current_trialX = new_initial_guess;
		}
		else
		{
			std::cout << "Converting between cartesian and polar coordinates only works for MGALT and FBLT mission types" << std::endl;
			throw 98;
		}

		return;
	}

	//function to find constraint/objective function dependencies due to an escape spiral anywhere in the mission
	void mission::find_dependencies_due_to_escape_spiral(vector<double>* Xupperbounds,
														vector<double>* Xlowerbounds,
														vector<double>* Fupperbounds,
														vector<double>* Flowerbounds,
														vector<string>* Xdescriptions,
														vector<string>* Fdescriptions,
														vector<int>* iAfun,
														vector<int>* jAvar,
														vector<int>* iGfun,
														vector<int>* jGvar,
														vector<string>* Adescriptions,
														vector<string>* Gdescriptions,
														missionoptions* options,
														const int& Findex)
	{
		//loop over journeys
		for (int jj = 0; jj < options->number_of_journeys; ++jj)
		{
			//we don't have any dependencies if there is no escape spiral in that journey
			if (options->journey_departure_type[jj] == 5)
			{
				//loop over all variables in the decision vector and determine the first entry in the journey of interest
				int first_entry_in_jj;
				stringstream pjprefix_stream;
				pjprefix_stream << "j" << jj << "p";
				string pjprefix = pjprefix_stream.str();
				for (int Xentry = 0; Xentry < Xdescriptions->size() - 1; ++Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find(pjprefix) < 1024)
					{
						first_entry_in_jj = Xentry;
						break;
					}
				}//end loop over all variables in the decision vector

				//this constraint has a derivative with respect to the mass at the beginning of any journey that has an escape spiral
				//therefore we have derivatives with respect to:
				//1. if the first journey, the initial mass scale factor if enabled
				if (jj == 0 && options->allow_initial_mass_to_vary)
				{
					for (int Xentry = 0; Xentry < Xdescriptions->size() - 1; ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("initial mass multiplier (0-1)") < 1024)
						{
							//first check for duplicates
							bool duplicateflag = false;
							stringstream entry_tag_stream;
							entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
							for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
							{
								if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Findex);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}
				//2. if NOT the first journey, the arrival mass at the end of the previous journey
				if (jj > 0)
				{
					for (int Xentry = first_entry_in_jj; Xentry > 0; --Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("arrival mass") < 1024)
						{
							//first check for duplicates
							bool duplicateflag = false;
							stringstream entry_tag_stream;
							entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
							for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
							{
								if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Findex);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}
				//3. if present, the journey initial mass increment scale factor
				//we must first check for duplicates associated with the current constraint, because this can occur
				
				for (int Xentry = first_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find("journey initial mass scale factor") < 1024)
					{
						bool duplicateflag = false;
						stringstream entry_tag_stream;
						entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
						for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
						{
							if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
							{
								duplicateflag = true;
								break;
							}
						}
						if (!duplicateflag)
						{
							iGfun->push_back(Findex);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
							break;
						}
					}

				}

				//this constraint has a derivative with respect to the Isp at the beginning of any journey that has an escape spiral
				if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
				{
					for (int Xentry = first_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("Escape spiral Isp") < 1024)
						{
							//first check for duplicates
							bool duplicateflag = false;
							stringstream entry_tag_stream;
							entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
							for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
							{
								if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Findex);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}

				//all spirals have a dependency on the BOL power if it is a variable
				if (options->objective_type == 13)
				{
					for (int Xentry = first_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("engine input power (kW)") < 1024 )
						{
							bool duplicateflag = false;
							for (int XXentry = Xdescriptions->size() - 1; XXentry > 0; --XXentry)
							{
								stringstream tempstream;
								tempstream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
								if ( (*Xdescriptions)[XXentry].find("engine input power (kW)") < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}
			}
		}//end loop over journeys
	}

	//function to find constraint/objective function dependencies due to a capture spiral anywhere in the mission
	void mission::find_dependencies_due_to_capture_spiral(vector<double>* Xupperbounds,
														vector<double>* Xlowerbounds,
														vector<double>* Fupperbounds,
														vector<double>* Flowerbounds,
														vector<string>* Xdescriptions,
														vector<string>* Fdescriptions,
														vector<int>* iAfun,
														vector<int>* jAvar,
														vector<int>* iGfun,
														vector<int>* jGvar,
														vector<string>* Adescriptions,
														vector<string>* Gdescriptions,
														missionoptions* options,
														const int& Findex)
	{
		//loop over journeys
		for (int jj = 0; jj < options->number_of_journeys; ++jj)
		{
			if (options->journey_arrival_type[jj] == 7)
			{
				//we have a dependency on the arrival mass and the capture spiral Isp (if applicable) for each preceding journey that has a capture spiral
				//loop over all variables in the decision vector and determine the first entry in the journey of interest
				int last_entry_in_jj;
				stringstream pjprefix_stream;
				pjprefix_stream << "j" << jj << "p";
				string pjprefix = pjprefix_stream.str();
				for (int Xentry = Xdescriptions->size() - 1; Xentry > 0 ; --Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find(pjprefix) < 1024)
					{
						last_entry_in_jj = Xentry;
						break;
					}
				}//end loop over all variables in the decision vector

				//this constraint has a derivative with respect to the arrival mass at the end of any journey that has a capture spiral
				for (int Xentry = last_entry_in_jj; Xentry > 0; --Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find("arrival mass") < 1024)
					{
						//first check for duplicates
						bool duplicateflag = false;
						stringstream entry_tag_stream;
						entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
						for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
						{
							if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
							{
								duplicateflag = true;
								break;
							}
						}
						if (!duplicateflag)
						{
							iGfun->push_back(Findex);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
						}
						break;
					}
				}
				
				//this constraint has a derivative with respect to the Isp at the beginning of any journey that has a capture spiral
				if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
				{
					for (int Xentry = last_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("Capture spiral Isp") < 1024)
						{
							//first check for duplicates
							bool duplicateflag = false;
							stringstream entry_tag_stream;
							entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
							for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
							{
								if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Findex);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " constraint F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}

				//all spirals have a dependency on the BOL power if it is a variable
				if (options->objective_type == 13)
				{
					for (int Xentry = last_entry_in_jj; Xentry < Xdescriptions->size() - 1; --Xentry)
					{
						if ( (*Xdescriptions)[Xentry].find("engine input power (kW)") < 1024 )
						{
							bool duplicateflag = false;
							for (int XXentry = Xdescriptions->size() - 1; XXentry > 0; --XXentry)
							{
								stringstream tempstream;
								tempstream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
								if ( (*Xdescriptions)[XXentry].find("engine input power (kW)") < 1024)
								{
									duplicateflag = true;
									break;
								}
							}
							if (!duplicateflag)
							{
								iGfun->push_back(Fdescriptions->size() - 1);
								jGvar->push_back(Xentry);
								stringstream EntryNameStream;
								EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
								Gdescriptions->push_back(EntryNameStream.str());
								break;
							}
						}
					}
				}
			}
		}//end loop over journeys
	}

	//function to find constraint/objective function dependencies due to a spiral in the last journey
	void mission::find_dependencies_due_to_capture_spiral_in_final_journey(vector<double>* Xupperbounds,
																		vector<double>* Xlowerbounds,
																		vector<double>* Fupperbounds,
																		vector<double>* Flowerbounds,
																		vector<string>* Xdescriptions,
																		vector<string>* Fdescriptions,
																		vector<int>* iAfun,
																		vector<int>* jAvar,
																		vector<int>* iGfun,
																		vector<int>* jGvar,
																		vector<string>* Adescriptions,
																		vector<string>* Gdescriptions,
																		missionoptions* options,
																		const int& Findex)
	{
		//loop over journeys
		int jj = options->number_of_journeys - 1;

		if (options->journey_arrival_type[jj] == 7)
		{
			//we have a dependency on the arrival mass and the capture spiral Isp (if applicable) for each preceding journey that has a capture spiral
			//loop over all variables in the decision vector and determine the first entry in the journey of interest
			int last_entry_in_jj;
			stringstream pjprefix_stream;
			pjprefix_stream << "j" << jj << "p";
			string pjprefix = pjprefix_stream.str();
			for (int Xentry = Xdescriptions->size() - 1; Xentry > 0 ; --Xentry)
			{
				if ( (*Xdescriptions)[Xentry].find(pjprefix) < 1024)
				{
					last_entry_in_jj = Xentry;
					break;
				}
			}//end loop over all variables in the decision vector

			//this constraint has a derivative with respect to the arrival mass at the end of any journey that has a capture spiral
			for (int Xentry = last_entry_in_jj; Xentry > 0; --Xentry)
			{
				if ( (*Xdescriptions)[Xentry].find("arrival mass") < 1024)
				{
					//first check for duplicates
					bool duplicateflag = false;
					stringstream entry_tag_stream;
					entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
					for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
					{
						if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
						{
							duplicateflag = true;
							break;
						}
					}
					if (!duplicateflag)
					{
						iGfun->push_back(Findex);
						jGvar->push_back(Xentry);
						stringstream EntryNameStream;
						EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
						Gdescriptions->push_back(EntryNameStream.str());
					}
					break;
				}
			}
				
			//this constraint has a derivative with respect to the Isp at the beginning of any journey that has a capture spiral
			if (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13)
			{
				for (int Xentry = last_entry_in_jj; Xentry < Xdescriptions->size() - 1; ++Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find("Capture spiral Isp") < 1024)
					{
						//first check for duplicates
						bool duplicateflag = false;
						stringstream entry_tag_stream;
						entry_tag_stream << "F[" << Findex << "] with respect to X[" << Xentry << "]";
						for (int Gentry = Gdescriptions->size()-1; Gentry >=0; --Gentry)
						{
							if ( (*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
							{
								duplicateflag = true;
								break;
							}
						}
						if (!duplicateflag)
						{
							iGfun->push_back(Findex);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Findex] << " constraint F[" << Findex << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
							break;
						}
					}
				}
			}

			//all spirals have a dependency on the BOL power if it is a variable
			if (options->objective_type == 13)
			{
				for (int Xentry = last_entry_in_jj; Xentry < Xdescriptions->size() - 1; --Xentry)
				{
					if ( (*Xdescriptions)[Xentry].find("engine input power (kW)") < 1024 )
					{
						bool duplicateflag = false;
						for (int XXentry = Xdescriptions->size() - 1; XXentry > 0; --XXentry)
						{
							stringstream tempstream;
							tempstream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
							if ( (*Xdescriptions)[XXentry].find("engine input power (kW)") < 1024)
							{
								duplicateflag = true;
								break;
							}
						}
						if (!duplicateflag)
						{
							iGfun->push_back(Fdescriptions->size() - 1);
							jGvar->push_back(Xentry);
							stringstream EntryNameStream;
							EntryNameStream << "Derivative of " << (*Fdescriptions)[Fdescriptions->size()-1] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
							Gdescriptions->push_back(EntryNameStream.str());
							break;
						}
					}
				}
			}
		}
	}

	void mission::extract_objective_function_values(std::vector<double>& objective_functions)
	{
		for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
		{
			switch (options.outerloop_objective_function_choices[objective])
			{
			case 0: //BOL power at 1 AU
				objective_functions[objective] = this->options.power_at_1_AU;
				break;
			case 1: //Launch epoch (MJD)
				objective_functions[objective] = this->Xopt[0] / 86400.0;
				break;
			case 2: //Flight time (days)
				objective_functions[objective] = (this->current_epoch - this->journeys[0].phases[0].phase_start_epoch) / 86400.0;
				break;
			case 3: //number of thrusters
				objective_functions[objective] = this->options.number_of_engines;
				break;
			case 4: //Thruster
				objective_functions[objective] = this->options.engine_type;
				break;
			case 5: //Launch vehicle
				objective_functions[objective] = this->options.LV_type;
				break;
			case 6: //Final mass
				objective_functions[objective] = -(options.minimum_dry_mass > 0 ? this->dry_mass : this->journeys.back().phases.back().state_at_end_of_phase[6]);
				break;
			case 7: //Final journey mass increment (for maximizing sample return)
				objective_functions[objective] = -(this->journeys.back().phases.back().current_mass_increment * this->journeys.back().phases.back().journey_initial_mass_increment_scale_factor);
				break;
			case 8: //first journey departure C3
				objective_functions[objective] = this->journeys[0].phases[0].C3_departure;
				break;
			case 9: //last journey arrival C3
				objective_functions[objective] = this->journeys.back().phases.back().C3_arrival;
			case 10: //delta-V
				objective_functions[objective] = this->current_deltaV;
			}
		}
	}
	
} /* namespace EMTG */
