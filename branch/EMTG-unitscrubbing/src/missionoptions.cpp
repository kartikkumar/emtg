/*
 * missionoptions.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: Jacob
 */

#include "missionoptions.h"
#include "EMTG_math.h"

#include "SpiceUsr.h"

#include "boost/algorithm/string.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

namespace EMTG {

missionoptions::missionoptions() {
	this->outerloop_warmstart = 0;

	this->outerloop_vary_power = false;
	this->outerloop_vary_launch_epoch = false;
	this->outerloop_vary_flight_time_upper_bound = false;
	this->outerloop_vary_thruster_type = false;
	this->outerloop_vary_number_of_thrusters = false;
	this->outerloop_vary_launch_vehicle = false;
	this->outerloop_vary_departure_C3 = false;
	this->outerloop_vary_arrival_C3 = false;
	this->outerloop_restrict_flight_time_lower_bound = false;
	this->outerloop_reevaluate_full_population = false;
	this->outerloop_warm_population = "none";
	this->outerloop_warm_archive = "none";
	this->quiet_outerloop = 1;
	this->quiet_basinhopping = false;
	this->MBH_two_step = false;
	this->FD_stepsize = 1.5e-8;
	this->FD_stepsize_coarse = 1.5e-3;
	this->control_coordinate_system = 0;
	this->initial_guess_control_coordinate_system = 0; 
	this->enable_maximum_propellant_mass_constraint = false;
	this->maximum_propellant_mass = 1000.0;

	this->spiral_model_type = 1;
	this->problem_type = 0;
	this->IspLT_minimum = 3000;
	this->LV_adapter_mass = 0.0;
	this->NLP_solver_type = 0;
	this->NLP_solver_mode = true;
	this->quiet_NLP = false;
	this->ACE_feasible_point_finder = false;
	this->MBH_hop_distribution = 1;
	this->MBH_Pareto_alpha = 3.0;
	this->MBH_time_hop_probability = 0.2;
	this->interpolate_initial_guess = false;
	this->seed_MBH = false;
	this->MBH_zero_control_initial_guess = 1;
	this->AU = 1.49597870691e+8;
	this->snopt_max_run_time = 3600;
	this->power_decay_rate = 0.0;
	this->throttle_logic_mode = 0;
	this->throttle_sharpness = 100.0;
	this->post_mission_delta_v = 0.0;
	this->post_mission_Isp = 3000.0;
	this->propellant_margin = 0.0;
	this->create_GMAT_script = 0;
	this->forced_flyby_coast = 0.0;
	this->forced_post_launch_coast = 0.0;
	this->power_margin = 0.0;
	this->number_of_journeys = 1;

	this->file_status = parse_options_file("options.emtgopt");

	this->construct_thruster_launch_vehicle_name_arrays();
}

missionoptions::missionoptions(string optionsfile) {
	this->outerloop_warmstart = 0;

	this->outerloop_vary_power = false;
	this->outerloop_vary_launch_epoch = false;
	this->outerloop_vary_flight_time_upper_bound = false;
	this->outerloop_vary_thruster_type = false;
	this->outerloop_vary_number_of_thrusters = false;
	this->outerloop_vary_launch_vehicle = false;
	this->outerloop_vary_departure_C3 = false;
	this->outerloop_vary_arrival_C3 = false;
	this->outerloop_restrict_flight_time_lower_bound = false;
	this->outerloop_reevaluate_full_population = false;
	this->outerloop_warm_population = "none";
	this->outerloop_warm_archive = "none";
	this->quiet_outerloop = true;
	this->quiet_basinhopping = false;
	this->MBH_two_step = false;
	this->FD_stepsize = 1.5e-8;
	this->FD_stepsize_coarse = 1.5e-3;
	this->control_coordinate_system = 0;
	this->initial_guess_control_coordinate_system = 0;
	this->enable_maximum_propellant_mass_constraint = false;
	this->maximum_propellant_mass = 1000.0;

	this->spiral_model_type = 1;
	this->problem_type = 0;
	this->IspLT_minimum = 3000;
	this->LV_adapter_mass = 0.0;
	this->NLP_solver_type = 0;
	this->NLP_solver_mode = true;
	this->quiet_NLP = false;
	this->ACE_feasible_point_finder = false;
	this->MBH_hop_distribution = 1;
	this->MBH_Pareto_alpha = 3.0;
	this->MBH_time_hop_probability = 0.1;
	this->interpolate_initial_guess = false;
	this->seed_MBH = false;
	this->MBH_zero_control_initial_guess = 1;
	this->AU = 1.49597870691e+8;
	this->snopt_max_run_time = 3600;
	this->power_decay_rate = 0.0;
	this->throttle_logic_mode = 0;
	this->throttle_sharpness = 100.0;
	this->post_mission_delta_v = 0.0;
	this->post_mission_Isp = 3000.0;
	this->propellant_margin = 0.0;
	this->create_GMAT_script = 0;
	this->forced_flyby_coast = 0.0;
	this->forced_post_launch_coast = 0.0;
	this->power_margin = 0.0;



	this->file_status = parse_options_file(optionsfile);

	this->construct_thruster_launch_vehicle_name_arrays();
}

missionoptions::~missionoptions() {

}

//function to parse an options file
int missionoptions::parse_options_file(string optionsfile) {
	ifstream inputfile(optionsfile.c_str());
	int linenumber = 0;
	string choice;
	string peek;
	double value;
	char dump_buffer[1024];
	
	if (!inputfile.is_open())
	{
		this->error_message = "Cannot find options file: " + optionsfile;
		return 2;
	}

	while (!inputfile.eof()) {
		peek = inputfile.peek();
		if (peek == "#" || peek == "\r" || peek == "\n") {
			//comment or blank line, do not parse
			inputfile.getline(dump_buffer, 1024);
			++linenumber;
		}
		else 
		{
			int returncode = parse_options_line(inputfile, choice, value, dump_buffer, linenumber);

			if (returncode)
			{
				stringstream errorstream;
				errorstream << "Failure reading '" << optionsfile << "' at line " << linenumber + 1 << ", '" << choice <<"'";
				this->error_message = errorstream.str();

				return 1;
			}
		}

		choice = "reset for next line";
	}

	//if we made it this far, this was a successful read
	this->error_message = "Options file '" + optionsfile + "' read successfully";

	//make sure we aren't choosing any engine parameters for impulsive mission types!
	if (this->mission_type == 0 || this->mission_type == 1 || this->mission_type == 5) //impulsive mission types
		this->engine_type = 0;

	return 0;
}

int missionoptions::parse_options_line(ifstream& inputfile, string& choice, double& value, char* dump_buffer, int& linenumber)
{

	inputfile >> choice;
	//handle options that have been removed
	if (choice == "minimum_mass_for_third_body_perturbation" || choice == "adapter_mass" || choice == "NeuroSpiral_neurons_per_layer" || choice == "NeuroSpiral_number_of_layers" || choice == "journey_capture_spiral_starting_radius")
	{
		//do not parse
		inputfile.getline(dump_buffer, 1024);
		return 0;
	}
	//we have to check for string parameters first
	if (choice == "universe_folder") {
		inputfile >> this->universe_folder;
		return 0;
	}	
	if (choice == "mission_name")
	{
		inputfile >> this->mission_name;
		return 0;
	}
	if (choice == "SPICE_leap_seconds_kernel")
	{
		inputfile >> this->SPICE_leap_seconds_kernel;
		return 0;
	}
	if (choice == "SPICE_reference_frame_kernel")
	{
		inputfile >> this->SPICE_reference_frame_kernel;
		return 0;
	}
	if (choice == "journey_names")
	{
		string temp;
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			temp.clear();
			inputfile >> temp;
			if (temp == "\n")
			{
				cout << "Need to assign a name to each journey, line " << linenumber << endl;
				return 1;
			}
			this->journey_names.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_central_body")
	{
		string temp;
		inputfile >> temp;
		this->journey_central_body.push_back(temp);
		for (int j = 1; j < number_of_journeys; ++j)
		{
			temp.clear();
			inputfile >> temp;
			if (temp == "\n")
			{
				cout << "Need to assign a central body to each journey, line " << linenumber << endl;
				return 1;
			}
			this->journey_central_body.push_back(temp);
		}
		return 0;
	}

	if (choice == "outerloop_warm_archive")
	{
		inputfile >> this->outerloop_warm_archive;
		return 0;
	}

	if (choice == "outerloop_warm_population")
	{
		inputfile >> this->outerloop_warm_population;
		return 0;
	}

	inputfile >> value;

	//problem type
	if (choice == "problem_type")
	{
		this->problem_type = (int) value;
		return 0;
	}

	//physical constants
	if (choice ==  "G") 
	{
		this->G = value;
		return 0;
	}
	if (choice ==  "g0") 
	{
		this->g0 = value;
		return 0;
	}

	//ephemeris settings
	if (choice == "ephemeris_source")
	{
		this->ephemeris_source = (int) value;
		return 0;
	}

	//outer loop solver settings
	if (choice ==  "run_outerloop") 
	{
		this->run_outerloop = (bool) value;
		return 0;
	}
	if (choice == "outerloop_popsize")
	{
		this->outerloop_popsize = (int) value;
		return 0;
	}
	if (choice == "outerloop_genmax") 
	{
		this->outerloop_genmax = (int) value;
		return 0;
	}
	if (choice == "outerloop_tournamentsize") 
	{
		this->outerloop_tournamentsize = (int) value;
		return 0;
	}
	if (choice == "outerloop_CR") 
	{
		this->outerloop_CR = value;
		return 0;
	}
	if (choice == "outerloop_mu") 
	{
		this->outerloop_mu = value;
		return 0;
	}
	if (choice == "outerloop_stallmax") 
	{
		this->outerloop_stallmax = (int) value;
		return 0;
	}
	if (choice == "outerloop_tolfit") 
	{
		this->outerloop_tolfit = value;
		return 0;
	}
	if (choice == "outerloop_ntrials") 
	{
		this->outerloop_ntrials = (int) value;
		return 0;
	}
	if (choice == "outerloop_elitecount") 
	{
		this->outerloop_elitecount = (int) value;
		return 0;
	}
	if (choice == "outerloop_useparallel") 
	{
		this->outerloop_useparallel = (bool) value;
		return 0;
	}
	if (choice == "outerloop_warmstart") 
	{
		this->outerloop_warmstart = (int) value;
		return 0;
	}
	if (choice == "outerloop_reevaluate_full_population")
	{
		this->outerloop_reevaluate_full_population = (bool) value;
		return 0;
	}
	if (choice == "quiet_outerloop")
	{
		this->quiet_outerloop = (bool) value;
		return 0;
	}

	//inner loop solver settings
	if (choice == "NLP_solver_type") {
		this->NLP_solver_type = (int) value;
		return 0;
	}
	if (choice == "NLP_solver_mode") {
		this->NLP_solver_mode = (bool) value;
		return 0;
	}
	if (choice == "ACE_feasible_point_finder") {
		this->ACE_feasible_point_finder = (bool) value;
		return 0;
	}
	if (choice == "quiet_NLP") {
		this->quiet_NLP = (bool) value;
		return 0;
	}
	if (choice == "quiet_basinhopping") {
		this->quiet_basinhopping = (bool) value;
		return 0;
	}
	if (choice == "MBH_max_not_improve") {
		this->MBH_max_not_improve = (int) value;
		return 0;
	}
	if (choice == "MBH_max_trials") {
		this->MBH_max_trials = (int) value;
		return 0;
	}
	if (choice == "MBH_max_run_time") {
		this->MBH_max_run_time = (int) value;
		return 0;
	}
	if (choice == "MBH_max_step_size") {
		this->MBH_max_step_size = value;
		return 0;
	}
	if (choice == "MBH_hop_distribution") {
		this->MBH_hop_distribution = value;
		return 0;
	}
	if (choice == "MBH_Pareto_alpha") {
		this->MBH_Pareto_alpha = value;
		return 0;
	}
	if (choice == "MBH_time_hop_probability")
	{
		this->MBH_time_hop_probability = value;
		return 0;
	}
	if (choice == "snopt_feasibility_tolerance") {
		this->snopt_feasibility_tolerance = value;
		return 0;
	}
	if (choice == "snopt_major_iterations") {
		this->snopt_major_iterations = (int) value;
		return 0;
	}
	if (choice == "snopt_max_run_time") {
		this->snopt_max_run_time = (time_t) value;
		return 0;
	}
	if (choice == "derivative_type") {
		this->derivative_type = (int) value;
		return 0;
	}
	if (choice == "seed_MBH") {
		this->seed_MBH = (bool) value;
		return 0;
	}
	if (choice == "interpolate_initial_guess") {
		this->interpolate_initial_guess = (bool) value;
		return 0;
	}
	if (choice == "initial_guess_num_timesteps") {
		this->initial_guess_num_timesteps = (int) value;
		return 0;
	}
	if (choice == "initial_guess_step_size_distribution") {
		this->initial_guess_step_size_distribution = (int) value;
		return 0;
	}
	if (choice == "initial_guess_step_size_stdv_or_scale") {
		this->initial_guess_step_size_stdv_or_scale = value;
		return 0;
	}
	if (choice == "MBH_zero_control_initial_guess") {
		this->MBH_zero_control_initial_guess = (int) value;
		return 0;
	}
	if (choice == "MBH_two_step") {
		this->MBH_two_step = (bool) value;
		return 0;
	}
	if (choice == "FD_stepsize") {
		this->FD_stepsize = value;
		return 0;
	}
	if (choice == "FD_stepsize_coarse") {
		this->FD_stepsize_coarse = value;
		return 0;
	}
				
	//low thrust solver parameters
	if (choice == "num_timesteps") 
	{
		this->num_timesteps = (int) value;
		return 0;
	}
	if (choice == "control_coordinate_system")
	{
		this->control_coordinate_system = (int) value;
		return 0;
	}
	if (choice == "initial_guess_control_coordinate_system")
	{
		this->initial_guess_control_coordinate_system = (int) value;
		return 0;
	}
	if (choice == "step_size_distribution")
	{
		this->step_size_distribution = (int) value;
		return 0;
	}
	if (choice == "step_size_stdv_or_scale") 
	{
		this->step_size_stdv_or_scale = value;
		return 0;
	}
	if (choice == "spiral_model_type")
	{
		this->spiral_model_type = (int) value;
		return 0;
	}

	//vehicle parameters
	if (choice == "maximum_mass") {
		this->maximum_mass = value;
		return 0;
	}
	if (choice == "allow_initial_mass_to_vary") {
		this->allow_initial_mass_to_vary = (bool) value;
		return 0;
	}
	if (choice == "LV_margin") {
		this->LV_margin = value;
		return 0;
	}
	if (choice == "LV_adapter_mass") {
		this->LV_adapter_mass = value;
		return 0;
	}
	if (choice == "custom_LV_coefficients") {
		this->custom_LV_coefficients[0] = value;

		for (int k = 1; k < 6; ++k)
		{
			inputfile >> value;
			this->custom_LV_coefficients[k] = value;
		}
		return 0;
	}
	if (choice == "custom_LV_C3_bounds") {
		this->custom_LV_C3_bounds[0] = value;
		inputfile >> value;
		this->custom_LV_C3_bounds[1] = value;
		return 0;
	}
	if (choice == "parking_orbit_altitude") {
		this->parking_orbit_altitude = (int) value;
		return 0;
	}
	if (choice == "parking_orbit_inclination") {
		this->parking_orbit_inclination = (int) value;
		return 0;
	}
	if (choice == "IspLT") {
		this->IspLT = value;
		return 0;
	}
	if (choice == "IspLT_minimum") {
		this->IspLT_minimum = value;
		return 0;
	}
	if (choice == "IspChem") {
		this->IspChem = value;
		return 0;
	}
	if (choice == "IspDS") {
		this->IspDS = value;
		return 0;
	}
	if (choice == "Thrust") {
		this->Thrust = value;
		return 0;
	}
	if (choice == "LV_type") {
		this->LV_type = (int) value;
		return 0;
	}
	if (choice == "engine_type") {
		this->engine_type = (int) value;
		return 0;
	}
	if (choice == "number_of_engines") {
		this->number_of_engines = (int) value;
		return 0;
	}
	if (choice == "throttle_logic_mode") {
		this->throttle_logic_mode = (int) value;
		return 0;
	}
	if (choice == "throttle_sharpness") {
		this->throttle_sharpness = value;
		return 0;
	}
	if (choice == "power_source_type") {
		this->power_source_type = (int) value;
		return 0;
	}
	if (choice == "power_at_1_AU") {
		this->power_at_1_AU = value;
		return 0;
	}
	if (choice == "solar_power_gamma") {
		this->solar_power_gamma[0] = value;
		for (int k = 1; k < 5; ++k)
		{
			inputfile >> value;
			this->solar_power_gamma[k] = value;
		}
		return 0;
	}
	if (choice == "power_margin") {
		this->power_margin = value;
		return 0;
	}
	if (choice == "spacecraft_power_coefficients") {
		this->spacecraft_power_coefficients[0] = value;
		for (int k = 1; k < 3; ++k)
		{
			inputfile >> value;
			this->spacecraft_power_coefficients[k] = value;
		}
		return 0;
	}
	if (choice == "engine_input_thrust_coefficients") {
		this->engine_input_thrust_coefficients[0] = value;

		string peek;
		peek = inputfile.peek();
		int count = 1;
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->engine_input_thrust_coefficients[count] = value;

			peek = inputfile.peek();
			++count;
		}

		for (int k = count; k < 7; ++k)
			this->engine_input_thrust_coefficients[k] = 0.0;

		return 0;
	}
	if (choice == "engine_input_mass_flow_rate_coefficients") {
		this->engine_input_mass_flow_rate_coefficients[0] = value;

		string peek;
		peek = inputfile.peek();
		int count = 1;
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->engine_input_mass_flow_rate_coefficients[count] = value;

			peek = inputfile.peek();
			++count;
		}

		for (int k = count; k < 7; ++k)
			this->engine_input_mass_flow_rate_coefficients[k] = 0.0;

		return 0;
	}
	if (choice == "engine_input_power_bounds") {
		this->engine_input_power_bounds[0] = value;
		inputfile >> value;
		this->engine_input_power_bounds[1] = value;
		return 0;
	}
	if (choice == "user_defined_engine_efficiency") {
		this->user_defined_engine_efficiency = value;
		return 0;
	}
	if (choice == "spacecraft_power_model_type") {
		this->spacecraft_power_model_type = (int) value;
		return 0;
	}
	if (choice == "EP_dry_mass") {
		this->EP_dry_mass = value;
		return 0;
	}
	if (choice == "engine_duty_cycle") {
		this->engine_duty_cycle = value;
		return 0;
	}
	if (choice == "power_decay_rate") {
		this->power_decay_rate = value;
		return 0;
	}

	//minimum dry mass constraint and related parameters
	if (choice == "minimum_dry_mass") {
		this->minimum_dry_mass = value;
		return 0;
	}
	if (choice == "enable_maximum_propellant_mass_constraint")
	{
		this->enable_maximum_propellant_mass_constraint = (bool)value;
		return 0;
	}
	if (choice == "maximum_propellant_mass")
	{
		this->maximum_propellant_mass = value;
		return 0;
	}
	if (choice == "post_mission_delta_v") {
		this->post_mission_delta_v = value;
		return 0;
	}
	if (choice == "post_mission_Isp") {
		this->post_mission_Isp = value;
		return 0;
	}
	if (choice == "propellant_margin") {
		this->propellant_margin = value;
		return 0;
	}
				
				
	//perturbation settings
	if (choice == "perturb_SRP") {
		this->perturb_SRP = (bool)value;
		return 0;
	}
	if (choice == "perturb_thirdbody") {
		this->perturb_thirdbody = (bool)value;
		return 0;
	}
	if (choice == "journey_perturbation_bodies") {
		//first get the number of perturbation bodies for each journey
		this->journey_number_of_perturbation_bodies.push_back((int) value);

		for (int j = 1; j < this->number_of_journeys; ++j)
		{
			inputfile >> value;
			this->journey_number_of_perturbation_bodies.push_back((int) value);
		}

		//next read the bodies line for each journey
		vector<int> temp_perturbation_list;
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			temp_perturbation_list.clear();
			for (int b = 0; b < this->journey_number_of_perturbation_bodies[j]; ++b)
			{
				inputfile >> value;
				temp_perturbation_list.push_back((int) value);
			}
			this->journey_perturbation_bodies.push_back(temp_perturbation_list);
			++linenumber;
		}

		return 0;
	}
	if (choice == "spacecraft_area") {
		this->spacecraft_area = value;
		return 0;
	}
	if (choice == "coefficient_of_reflectivity") {
		this->coefficient_of_reflectivity = value;
		return 0;
	}

	//outer loop selectable options settings
	if (choice == "outerloop_vary_power") {
		this->outerloop_vary_power = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_launch_epoch") {
		this->outerloop_vary_launch_epoch = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_flight_time_upper_bound") {
		this->outerloop_vary_flight_time_upper_bound = (bool) value;
		return 0;
	}
	if (choice == "outerloop_restrict_flight_time_lower_bound") {
		this->outerloop_restrict_flight_time_lower_bound = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_thruster_type") {
		this->outerloop_vary_thruster_type = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_number_of_thrusters") {
		this->outerloop_vary_number_of_thrusters = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_launch_vehicle") {
		this->outerloop_vary_launch_vehicle = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_departure_C3") {
		this->outerloop_vary_departure_C3 = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_arrival_C3") {
		this->outerloop_vary_arrival_C3 = (bool) value;
		return 0;
	}
	if (choice == "outerloop_vary_journey_destination") {
		this->outerloop_vary_journey_destination.push_back((bool) value);

		for (int j = 1; j < this->number_of_journeys; ++j)
		{
			inputfile >> value;
			this->outerloop_vary_journey_destination.push_back((bool) value);
		}

		return 0;
	}
	if (choice == "outerloop_vary_journey_flyby_sequence") {
		this->outerloop_vary_journey_flyby_sequence.push_back((bool) value);

		for (int j = 1; j < this->number_of_journeys; ++j)
		{
			inputfile >> value;
			this->outerloop_vary_journey_flyby_sequence.push_back((bool) value);
		}

		return 0;
	}
	if (choice == "outerloop_power_choices") {
		this->outerloop_power_choices.push_back(value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_power_choices.push_back(value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_launch_epoch_choices") {
		this->outerloop_launch_epoch_choices.push_back(value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_launch_epoch_choices.push_back(value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_flight_time_upper_bound_choices") {
		this->outerloop_flight_time_upper_bound_choices.push_back(value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_flight_time_upper_bound_choices.push_back(value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_thruster_type_choices") {
		this->outerloop_thruster_type_choices.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_thruster_type_choices.push_back((int) value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_number_of_thrusters_choices") {
		this->outerloop_number_of_thrusters_choices.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_number_of_thrusters_choices.push_back((int) value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_launch_vehicle_choices") {
		this->outerloop_launch_vehicle_choices.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_launch_vehicle_choices.push_back((int) value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_departure_C3_choices") {
		this->outerloop_departure_C3_choices.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_departure_C3_choices.push_back((int) value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_arrival_C3_choices") {
		this->outerloop_arrival_C3_choices.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_arrival_C3_choices.push_back((int) value);

			peek = inputfile.peek();
		}

		return 0;
	}
	if (choice == "outerloop_journey_destination_choices") 
	{
		vector<int> temp;
		temp.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			temp.push_back((int) value);

			peek = inputfile.peek();
		}
		this->outerloop_journey_destination_choices.push_back(temp);
		++linenumber;

		for (int j = 1; j < this->number_of_journeys; ++j)
		{
			temp.clear();
			inputfile >> value;
			temp.push_back((int) value);
			peek = inputfile.peek();
			while (!(peek == "\n" || peek == "#" || peek == "\r")) 
			{
				inputfile >> value;
				temp.push_back((int) value);

				peek = inputfile.peek();
			}
			this->outerloop_journey_destination_choices.push_back(temp);
			++linenumber;
		}
		
		return 0;
	}
	if (choice == "outerloop_journey_flyby_sequence_choices") 
	{
		vector<int> temp;
		temp.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			temp.push_back((int) value);

			peek = inputfile.peek();
		}
		this->outerloop_journey_flyby_sequence_choices.push_back(temp);
		++linenumber;

		for (int j = 1; j < this->number_of_journeys; ++j)
		{
			temp.clear();
			inputfile >> value;
			temp.push_back((int) value);
			peek = inputfile.peek();
			while (!(peek == "\n" || peek == "#" || peek == "\r")) 
			{
				inputfile >> value;
				temp.push_back((int) value);

				peek = inputfile.peek();
			}
			this->outerloop_journey_flyby_sequence_choices.push_back(temp);
			++linenumber;
		}
		
		return 0;
	}
	if (choice == "outerloop_journey_maximum_number_of_flybys") {
		this->outerloop_journey_maximum_number_of_flybys.push_back((int) value);

		for (int j = 1; j < this->number_of_journeys; ++j)
		{
			inputfile >> value;
			this->outerloop_journey_maximum_number_of_flybys.push_back((int) value);
		}

		return 0;
	}

	//outerloop objective settings
	if (choice == "outerloop_objective_function_choices") {
		this->outerloop_objective_function_choices.push_back((int) value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			this->outerloop_objective_function_choices.push_back((int) value);

			peek = inputfile.peek();
		}

		return 0;
	}

	//global problem settings
	if (choice == "number_of_journeys") {
		this->number_of_journeys = (int) value;
		this->number_of_phases.resize(this->number_of_journeys);
		return 0;
	}
	if (choice == "max_phases_per_journey") {
		this->max_phases_per_journey = (int) value;
		return 0;
	}
	if (choice == "destination_list") {
		//encode the first journey
		vector<int> journey_destination_list;
		journey_destination_list.push_back((int) value);

		inputfile >> value;

		journey_destination_list.push_back((int) value);

		this->destination_list.push_back(journey_destination_list);

		//encode successive journeys

		for (int j = 1; j < number_of_journeys; ++j)
		{
			journey_destination_list.clear();
			inputfile >> value;
			journey_destination_list.push_back((int) value);
			inputfile >> value;
			journey_destination_list.push_back((int) value);

			this->destination_list.push_back(journey_destination_list);
		}
		return 0;
	}
	if (choice == "include_initial_impulse_in_cost") {
		this->include_initial_impulse_in_cost = (bool) value;
		return 0;
	}
	if (choice == "global_timebounded") {
		this->global_timebounded = (bool) value;
		return 0;
	}
	if (choice == "launch_window_open_date") {
		this->launch_window_open_date = value * 86400.0;
		return 0;
	}
	if (choice == "total_flight_time_bounds") {
		this->total_flight_time_bounds[0] = value * 86400.0;
		inputfile >> value;
		this->total_flight_time_bounds[1] = value * 86400.0;
		return 0;
	}
	if (choice == "DLA_bounds") {
		this->DLA_bounds.push_back(value);
		inputfile >> value;
		this->DLA_bounds.push_back(value);
		return 0;
	}
	if (choice == "objective_type") {
		this->objective_type = (int) value;
		return 0;
	}
	if (choice == "mission_type") {
		this->mission_type = (int) value;
		return 0;
	}
	if (choice == "initial_V_infinity")
	{
		this->initial_V_infinity.push_back(value);
		inputfile >> value;
		this->initial_V_infinity.push_back(value);
		inputfile >> value;
		this->initial_V_infinity.push_back(value);
		return 0;
	}
	if (choice == "forced_post_launch_coast")
	{
		this->forced_post_launch_coast = value * 86400.0;
		return 0;
	}
	if (choice == "forced_flyby_coast")
	{
		this->forced_flyby_coast = value * 86400.0;
		return 0;
	}

	//settings for each journey
	if (choice == "journey_timebounded") {
		this->journey_timebounded.push_back((int) value);
		for (int k = 1; k < this->number_of_journeys; ++k) {
			inputfile >> value;
			this->journey_timebounded.push_back((int) value);
		}
		return 0;
	}
	if (choice == "journey_wait_time_bounds") {
		vector<double> temp(2);
		temp[0] = value;
		inputfile >> temp[1];
		temp[0] *= 86400.0;
		temp[1] *= 86400.0;
		this->journey_wait_time_bounds.push_back(temp);

		for (int k = 1; k < this->number_of_journeys; ++k) {
			inputfile >> temp[0];
			inputfile >> temp[1];
			temp[0] *= 86400.0;
			temp[1] *= 86400.0;
			this->journey_wait_time_bounds.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_flight_time_bounds") {
		vector<double> temp(2);
		temp[0] = value;
		inputfile >> temp[1];
		temp[0] *= 86400.0;
		temp[1] *= 86400.0;
		this->journey_flight_time_bounds.push_back(temp);

		for (int k=1; k < this->number_of_journeys; ++k) {
			inputfile >> temp[0];
			inputfile >> temp[1];
			temp[0] *= 86400.0;
			temp[1] *= 86400.0;
			this->journey_flight_time_bounds.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_arrival_date_bounds") {
		vector<double> temp(2);
		temp[0] = value;
		inputfile >> temp[1];
		temp[0] *= 86400.0;
		temp[1] *= 86400.0;
		this->journey_arrival_date_bounds.push_back(temp);

		for (int k = 1; k < this->number_of_journeys; ++k) {
			inputfile >> temp[0];
			inputfile >> temp[1];
			temp[0] *= 86400.0;
			temp[1] *= 86400.0;
			this->journey_arrival_date_bounds.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_initial_impulse_bounds") {
		vector<double> temp;
		temp.push_back(value > 1.0e-8 ? value : 1.0e-8);
		inputfile >> value;
		temp.push_back(value > 1.0e-6 ? value : 1.0e-8);
		this->journey_initial_impulse_bounds.push_back(temp);

		for (int k = 1; k < this->number_of_journeys; ++k) {
			temp.clear();
			inputfile >> value;
			temp.push_back(value > 1.0e-8 ? value : 1.0e-8);
			inputfile >> value;
			temp.push_back(value > 1.0e-8 ? value : 1.0e-8);
			this->journey_initial_impulse_bounds.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_arrival_type") {
		stringstream errorstream;
					
		this->journey_arrival_type.push_back((int) value);

		if (value > 7 || value < 0)
		{
			errorstream << "Invalid journey arrival type, journey 1" << endl;
			this->error_message = errorstream.str();
		}
					
		for (int k = 1; k < this->number_of_journeys; ++k) {
			inputfile >> value;
			this->journey_arrival_type.push_back((int) value);

			if (value > 7 || value < 0)
			{
				errorstream << "Invalid journey arrival type, journey " << k+1 << endl;
				this->error_message = errorstream.str();
			}
		}

		if (error_message.size())
			return 1;
		
		return 0;
	}
	if (choice == "journey_departure_type") {
		stringstream errorstream;
					
		this->journey_departure_type.push_back((int) value);

		if (value > 6 || value < 0)
		{
			errorstream << "Invalid journey departure type, journey 1" << endl;
			this->error_message = errorstream.str();
		}
					
		for (int k = 1; k < this->number_of_journeys; ++k) {
			inputfile >> value;
			this->journey_departure_type.push_back((int) value);

			if (value > 6 || value < 0)
			{
				errorstream << "Invalid journey departure type, journey " << k+1 << endl;
				this->error_message = errorstream.str();
			}			
		}

		if (error_message.size())
			return 1;

		
		return 0;
	}
				
	if (choice == "journey_departure_elements_vary_flag") {
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			vector <bool> temp;

			if (j > 0) //if j is 0, the correct value is already loaded
				inputfile >> value;
			temp.push_back((bool) value);

			for (int k = 1; k < 6; ++k)
			{
				inputfile >> value;
				temp.push_back((bool) value);
			}

			this->journey_departure_elements_vary_flag.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_departure_elements_type") {
		this->journey_departure_elements_type.push_back((int) value);
		for (int k = 1; k < this->number_of_journeys; ++k) {
			inputfile >> value;
			this->journey_departure_elements_type.push_back((int) value);
		}
		return 0;
	}
	if (choice == "journey_departure_elements_bounds") {
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			vector < vector <double> > journeyset;
						

			for (int k = 0; k < 6; ++k)
			{
				vector<double> elementset;

				if (j > 0 || k > 0) //if j == 0 and k == 0, the correct value is already loaded
					inputfile >> value;

				elementset.push_back(value);

				if (this->journey_departure_elements_type[j] && k > 1) //element specified in COE as angle, must convert from degrees to radians
				{
					elementset[0] *= math::PI / 180.0;
				}

				inputfile >> value;
				elementset.push_back(value);

				if (this->journey_departure_elements_type[j] && k > 1) //element specified in COE as angle, must convert from degrees to radians
				{
					elementset[1] *= math::PI / 180.0;
				}

				journeyset.push_back(elementset);
			}

			this->journey_departure_elements_bounds.push_back(journeyset);
		}
		return 0;
	}
	if (choice == "journey_departure_elements") {
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			vector <double> temp;
			//semi-major axis or x
			if (j > 0) //if j is 0, the correct value is already loaded
				inputfile >> value;
			temp.push_back(value);

			//eccentricity or y
			inputfile >> value;
			temp.push_back(value);

			if (this->journey_departure_elements_type[j])
			{
				//angles, convert from degrees to radians
				for (int k=2; k<6; ++k) {
					inputfile >> value;
					temp.push_back(value * EMTG::math::PI / 180.0);
				}
			}
			else //inertial elements
			{
							
				for (int k=2; k<6; ++k)
				{
					inputfile >> value;
					temp.push_back(value);
				}
			}

			this->journey_departure_elements.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_arrival_elements_type") {
		this->journey_arrival_elements_type.push_back((int) value);
		for (int k = 1; k < this->number_of_journeys; ++k) {
			inputfile >> value;
			this->journey_arrival_elements_type.push_back((int) value);
		}
		return 0;
	}
	if (choice == "journey_arrival_elements") {
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			vector <double> temp;
			//semi-major axis or x
			if (j > 0) //if j is 0, the correct value is already loaded
				inputfile >> value;
			temp.push_back(value);

			//eccentricity or y
			inputfile >> value;
			temp.push_back(value);

			if (this->journey_arrival_elements_type[j])
			{
				//angles, convert from degrees to radians
				for (int k=2; k<6; ++k) {
					inputfile >> value;
					temp.push_back(value * EMTG::math::PI / 180.0);
				}
			}
			else //inertial elements
			{
				for (int k=2; k<6; ++k)
				{
					inputfile >> value;
					temp.push_back(value);
				}
			}

			this->journey_arrival_elements.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_arrival_elements_vary_flag") {
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			vector <bool> temp;

			if (j > 0) //if j is 0, the correct value is already loaded
				inputfile >> value;
			temp.push_back((bool) value);

			for (int k = 1; k < 6; ++k)
			{
				inputfile >> value;
				temp.push_back((bool) value);
			}

			this->journey_arrival_elements_vary_flag.push_back(temp);
		}
		return 0;
	}
	if (choice == "journey_arrival_elements_bounds") {
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			vector < vector <double> > journeyset;

			for (int k = 0; k < 6; ++k)
			{
				vector<double> elementset;

				if (j > 0 || k > 0) //if j == 0 and k == 0, the correct value is already loaded
					inputfile >> value;

				elementset.push_back(value);

				if (journey_arrival_elements_type[j] && k > 1) //element specified in COE as angle, must convert from degrees to radians
				{
					elementset[0] *= math::PI / 180.0;
				}

				inputfile >> value;
				elementset.push_back(value);

				if (journey_arrival_elements_type[j] && k > 1) //element specified in COE as angle, must convert from degrees to radians
				{
					elementset[1] *= math::PI / 180.0;
				}

				journeyset.push_back(elementset);
			}

			this->journey_arrival_elements_bounds.push_back(journeyset);
		}
		return 0;
	}
	if (choice == "journey_final_velocity") {
		stringstream errorstream;
		vector<double> dummy;
		dummy.push_back(value);
		inputfile >> value;
		dummy.push_back(value);
		inputfile >> value;
		dummy.push_back(value);

		this->journey_final_velocity.push_back(dummy);

		if (this->journey_arrival_type[0] == 2 && dummy[1] < dummy[0])
		{
			errorstream << "Journey final velocity magnitude upper bound cannot be less than lower bound" << endl;
			this->error_message = errorstream.str();		
		}

		for (int j=1; j<number_of_journeys; ++j) 
		{
			dummy.clear();
			inputfile >> value;
			dummy.push_back(value);
			inputfile >> value;
			dummy.push_back(value);
			inputfile >> value;
			dummy.push_back(value);
			journey_final_velocity.push_back(dummy);

			if (this->journey_arrival_type[0] == 2 && dummy[1] < dummy[0])
			{
				errorstream << "Journey final velocity magnitude upper bound cannot be less than lower bound" << endl;
				this->error_message = errorstream.str();		
			}
		}

		if (this->error_message.size())
			return 1;

		
		return 0;
	}
	if (choice == "journey_initial_velocity") {
		stringstream errorstream;
		vector<double> dummy;
		dummy.push_back(value);
		inputfile >> value;
		dummy.push_back(value);
		inputfile >> value;
		dummy.push_back(value);

		this->journey_initial_velocity.push_back(dummy);

		if (this->journey_arrival_type[0] == 2 && dummy[1] < dummy[0])
		{
			errorstream << "Journey initial velocity magnitude upper bound cannot be less than lower bound" << endl;
			this->error_message = errorstream.str();		
		}

		for (int j=1; j<number_of_journeys; ++j) 
		{
			dummy.clear();
			inputfile >> value;
			dummy.push_back(value);
			inputfile >> value;
			dummy.push_back(value);
			inputfile >> value;
			dummy.push_back(value);
			journey_initial_velocity.push_back(dummy);

			if (this->journey_arrival_type[0] == 2 && dummy[1] < dummy[0])
			{
				errorstream << "Journey initial velocity magnitude upper bound cannot be less than lower bound" << endl;
				this->error_message = errorstream.str();		
			}
		}

		if (this->error_message.size())
			return 1;
		
		return 0;
	}
	if (choice == "journey_starting_mass_increment")
	{
		this->journey_starting_mass_increment.push_back(value);

		for (int j = 1; j < this->number_of_journeys; ++j) 
		{
			inputfile >> value;
			this->journey_starting_mass_increment.push_back(value);
		}
		return 0;
	}
	if (choice == "journey_variable_mass_increment")
	{
		this->journey_variable_mass_increment.push_back(value);
					
		for (int j = 1; j < this->number_of_journeys; ++j) 
		{
			inputfile >> value;
			this->journey_variable_mass_increment.push_back(value);
		}
		return 0;
	}
	if (choice == "journey_arrival_declination_constraint_flag")
	{
		this->journey_arrival_declination_constraint_flag.push_back( (bool) value );
		for (int j = 1; j < this->number_of_journeys; ++j) 
		{
			inputfile >> value;
			this->journey_arrival_declination_constraint_flag.push_back( (bool) value );
		}
		return 0;
	}
	if (choice == "journey_arrival_declination_bounds")
	{
		vector<double> temp_bounds;

		temp_bounds.push_back(value * math::PI/180.0);
		inputfile >> value;
		temp_bounds.push_back(value * math::PI/180.0);
		this->journey_arrival_declination_bounds.push_back( temp_bounds );

		for (int j = 1; j < this->number_of_journeys; ++j) 
		{
			vector<double> temp_bounds;
			inputfile >> value;
			temp_bounds.push_back(value * math::PI/180.0);
			inputfile >> value;
			temp_bounds.push_back(value * math::PI/180.0);
			this->journey_arrival_declination_bounds.push_back( temp_bounds );
		}
		return 0;
	}
	if (choice == "journey_escape_spiral_starting_radius")
	{
		this->journey_escape_spiral_starting_radius.push_back(value);
					
		for (int j = 1; j < this->number_of_journeys; ++j) 
		{
			inputfile >> value;
			this->journey_escape_spiral_starting_radius.push_back(value);
		}
		return 0;
	}
	if (choice == "journey_capture_spiral_final_radius")
	{
		this->journey_capture_spiral_final_radius.push_back(value);
					
		for (int j = 1; j < this->number_of_journeys; ++j) 
		{
			inputfile >> value;
			this->journey_capture_spiral_final_radius.push_back(value);
		}
		return 0;
	}
				
	//output format settings
	if (choice == "output_units") {
		this->output_units = (int) value;
		return 0;
	}
	if (choice == "create_GMAT_script") {
		this->create_GMAT_script = (int) value;
		return 0;
	}

	//debug code
	if (choice == "check_derivatives") {
		this->check_derivatives = (bool) value;
		return 0;
	}

	if (choice == "run_inner_loop") {
		this->run_inner_loop = (int) value;

		if (this->run_inner_loop > 4 || this->run_inner_loop < 0)
		{
			stringstream errorstream;
			errorstream << "Invalid inner-loop solver selection: " << run_inner_loop << ", acceptable values are (0-4)";
			this->error_message = errorstream.str();
			return 1;
		}
		return 0;
	}
	if (choice == "sequence") 
	{
		this->number_of_trial_sequences = (int) value;
		for (int seq = 0; seq < this->number_of_trial_sequences; ++seq)
		{
			++linenumber;

			//read in the mission sequence
			vector<int> temp_journey;
			vector<int> temp_number_of_phases;
			vector< vector<int> > temp_mission;

			for (int j= 0; j < this->number_of_journeys; ++j) 
			{
				int nphases = 1;
				temp_journey.clear();
				for (int p = 0; p < this->max_phases_per_journey; ++p) {
					inputfile >> value;
					if ((int) value)
					{
						temp_journey.push_back((int) value);
						++nphases;
					}
					else
						temp_journey.push_back(0);
				}

				temp_mission.push_back(temp_journey);
				temp_number_of_phases.push_back(nphases);
			}

			//append the mission sequence to the list of journeys
			this->sequence_input.push_back(temp_mission);
			this->number_of_phases_input.push_back(temp_number_of_phases);
		}
		return 0;
	}
	if (choice == "phase_type") {

		if (this->run_outerloop == 0 && this->mission_type > 6) //if we need to specify phase types
			{
			vector<int> temp;
			temp.push_back((int) value);

			for (int p = 1; p < this->number_of_phases[0]; ++p) {
				inputfile >> value;
				temp.push_back((int) value);
			}

			this->phase_type_input.push_back(temp);

			for (int j = 1; j < this->number_of_journeys; ++j) {
				temp.clear();
				for (int p =  0; p < this->number_of_phases[j]; ++p) {
					inputfile >> value;

					temp.push_back((int) value);
				}

				this->phase_type_input.push_back(temp);
			}
		}
		else //we don't need to specify phase types, so dump the line and continue
			inputfile.getline(dump_buffer, 1024);

		
		return 0;
	}
	if (choice == "trialX") 
	{
		vector<double> temp;
		temp.push_back(value);

		string peek;
		peek = inputfile.peek();
		while (!(peek == "\n" || peek == "#" || peek == "\r")) 
		{
			inputfile >> value;
			temp.push_back(value);

			peek = inputfile.peek();
		}
		this->trialX.push_back(temp);

		for (int entry = 1; entry < this->number_of_trial_sequences; ++entry)
		{
			temp.clear();
			inputfile >> value;
			temp.push_back(value);
			peek = inputfile.peek();
			while (!(peek == "\n" || peek == "#" || peek == "\r")) 
			{
				inputfile >> value;
				temp.push_back(value);

				peek = inputfile.peek();
			}
			this->trialX.push_back(temp);
		}
		
		return 0;

	}

	if (choice == "reset for next line") {return 0;}

	//or, if we have a string that did not match any option (i.e. we got this far)
	return 1;

}

int missionoptions::print_options_file(string filename) {
	//options print function should print out the entire options file including the automatically generated parameters
	//this allows us to check to make sure the options file was read in and processed correctly

	ofstream outputfile(filename.c_str(), ios::trunc);

	if (outputfile.is_open())
	{

		outputfile << "##Options file for EMTG_v8" << endl;
		outputfile << "##Written by EMTG_v8 core program compiled " << __DATE__ << " " << __TIME__ << endl;
		outputfile << endl;
		outputfile << "##problem type" << endl;
		outputfile << "#0: standard EMTG mission" << endl;
		outputfile << "problem_type " << this->problem_type << endl;
		outputfile << endl;
		outputfile << "##physical constants" << endl;
		outputfile << "#G in km^3/kg/s^2" << endl;
		outputfile << "G " << this->G << endl;
		outputfile << "#gravity at sea level on Earth in m/s^2" << endl;
		outputfile << "g0 " << this->g0 << endl;
		outputfile << endl;

		outputfile << "##outer-loop solver settings" << endl;
		outputfile << "#whether or not to run the outer-loop" << endl;
		outputfile << "run_outerloop " << this->run_outerloop << endl;
		outputfile << "#outer-loop population size" << endl;
		outputfile << "outerloop_popsize " << this->outerloop_popsize << endl;
		outputfile << "#maximum number of outer-loop generations" << endl;
		outputfile << "outerloop_genmax " << this->outerloop_genmax << endl;
		outputfile << "#tournament size for selection" << endl;
		outputfile << "outerloop_tournamentsize " << this->outerloop_tournamentsize << endl;
		outputfile << "#crossover ratio" << endl;
		outputfile << "outerloop_CR " << this->outerloop_CR << endl;
		outputfile << "#mutation rate" << endl;
		outputfile << "outerloop_mu " << this->outerloop_mu << endl;
		outputfile << "#maximum number of stall generations" << endl;
		outputfile << "outerloop_stallmax " << this->outerloop_stallmax << endl;
		outputfile << "#fitness tolerance for the outer-loop" << endl;
		outputfile << "outerloop_tolfit " << this->outerloop_tolfit << endl;
		outputfile << "#how many elite individuals to retain" << endl;
		outputfile << "outerloop_elitecount " << this->outerloop_elitecount << endl;
		outputfile << "#how many times to run the outer-loop" << endl;
		outputfile << "outerloop_ntrials " << this->outerloop_ntrials << endl;
		outputfile << "#whether or not to use the parallel outer-loop" << endl;
		outputfile << "outerloop_useparallel " << this->outerloop_useparallel << endl;
		outputfile << "#whether or not to perform an outer loop warm start" << endl;
		outputfile << "outerloop_warmstart " << this->outerloop_warmstart << endl;
		outputfile << "#Population file for outerloop warm start (set to none if not warm starting)" << endl;
		outputfile << "outerloop_warm_population " << this->outerloop_warm_population << endl;
		outputfile << "#Archive file for outerloop warm start (set to none if not warm starting)" << endl;
		outputfile << "outerloop_warm_archive " << this->outerloop_warm_archive << endl;
		outputfile << "#Re-evaluate the entire outerloop each generation? Otherwise read from the archive." << endl;
		outputfile << "outerloop_reevaluate_full_population " << this->outerloop_reevaluate_full_population << endl;
		outputfile << "#Quiet outer-loop?" << endl;
        outputfile << "quiet_outerloop " << this->quiet_outerloop << endl;
		outputfile << endl;

		outputfile << "##inner-loop solver settings" << endl;
		outputfile << "#NLP solver type" << endl;
		outputfile << "#0: SNOPT" << endl;
		outputfile << "#1: WORHP" << endl;
		outputfile << "NLP_solver_type " << this->NLP_solver_type << endl;
		outputfile << "#NLP solver mode" << endl;
		outputfile << "#0: find feasible point only" << endl;
		outputfile << "#1: find optimal solution" << endl;
		outputfile << "NLP_solver_mode " << this->NLP_solver_mode << endl;
		outputfile << "#Quiet NLP solver?" << endl;
		outputfile << "quiet_NLP " << this->quiet_NLP << endl;
		outputfile << "#Quiet MBH?" << endl;
		outputfile << "quiet_basinhopping " << this->quiet_basinhopping << endl;
		outputfile << "#Enable ACE feasible point finder?" << endl;
		outputfile << "ACE_feasible_point_finder " << this->ACE_feasible_point_finder << endl;
		outputfile << "#quantity 'Max_not_improve' for MBH" << endl;
		outputfile << "MBH_max_not_improve " << this->MBH_max_not_improve << endl;
		outputfile << "#maximum number of trials for MBH" << endl;
		outputfile << "MBH_max_trials " << this->MBH_max_trials << endl;
		outputfile << "#maximum run time for MBH, in seconds" << endl;
		outputfile << "MBH_max_run_time " << this->MBH_max_run_time << endl;
		outputfile << "#maximum step size for uniform MBH, or scaling factor for Cauchy MBH" << endl;
		outputfile << "MBH_max_step_size " << this->MBH_max_step_size << endl;
		outputfile << "#MBH hop probabilty distribution" << endl;
		outputfile << "#0: uniform" << endl;
		outputfile << "#1: Cauchy" << endl;
		outputfile << "#2: Pareto" << endl;
		outputfile << "#3: Gaussian" << endl;
		outputfile << "MBH_hop_distribution " << this->MBH_hop_distribution << endl;
		outputfile << "#Pareto distribution alpha" << endl;
		outputfile << "MBH_Pareto_alpha " << this->MBH_Pareto_alpha << endl;
		outputfile << "#probability of MBH time hop operation" << endl;
		outputfile << "MBH_time_hop_probability " << this->MBH_time_hop_probability << endl;
		outputfile << "#feasibility tolerance" << endl;
		outputfile << "snopt_feasibility_tolerance " << this->snopt_feasibility_tolerance << endl;
		outputfile << "#maximum number of major iterations for SNOPT" << endl;
		outputfile << "snopt_major_iterations " << this->snopt_major_iterations << endl;
		outputfile << "#Maximum run time, in seconds, for a single call to SNOPT" << endl;
		outputfile << "snopt_max_run_time " << this->snopt_max_run_time << endl;
		outputfile << "#method of specifying derivatives" << endl;
		outputfile << "#0: finite difference" << endl;
		outputfile << "#1: analytical flybys and objective function but finite difference the patch points" << endl;
		outputfile << "#2: all but time derivatives" << endl;
		outputfile << "#3: all but current phase flight time derivatives" << endl;
		outputfile << "#4: fully analytical (experimental)" << endl;
		outputfile << "derivative_type " << this->derivative_type << endl;
		outputfile << "#Will MBH be seeded with an initial point? Otherwise MBH starts from a completely random point." << endl;
		outputfile << "seed_MBH " << this->seed_MBH << endl;
		outputfile << "#Will the initial guess be interpolated?" << endl;
		outputfile << "#(i.e. are we solving a problem with a different number of time steps than the initial guess?)" << endl;
		outputfile << "interpolate_initial_guess " << this->interpolate_initial_guess << endl;
		outputfile << "#How many time steps were used to create the initial guess?" << endl;
		outputfile << "initial_guess_num_timesteps " << this->initial_guess_num_timesteps << endl;
		outputfile << "#Distribution from which the initial guess step sizes were drawn" << endl;
        outputfile << "#0: uniform" << endl;
        outputfile << "#1: Gaussian" << endl;
        outputfile << "#2: Cauchy" << endl;
        outputfile << "initial_guess_step_size_distribution " << this->initial_guess_step_size_distribution << endl;
        outputfile << "#What scale width (Cauchy) or standard deviation (Gaussian) was used to create the step sizes in the initial guess" << endl;
        outputfile << "initial_guess_step_size_stdv_or_scale " << this->initial_guess_step_size_stdv_or_scale << endl;
		outputfile << "#Apply zero-control initial guess in MBH?" << endl;
		outputfile << "#0: do not use" << endl;
		outputfile << "#1: zero-control for resets, random perturbations for hops" << endl;
		outputfile << "#2: always use zero-control guess except when seeded" << endl;
		outputfile << "MBH_zero_control_initial_guess " << this->MBH_zero_control_initial_guess << endl;
		outputfile << "#Enable two-step MBH?" << endl;
		outputfile << "MBH_two_step " << this->MBH_two_step << endl;
		outputfile << "#'Fine' finite differencing step size" << endl;
		outputfile << "FD_stepsize " << this->FD_stepsize << endl;
		outputfile << "#'Coarse' finite differencing step size" << endl;
		outputfile << "FD_stepsize_coarse " << this->FD_stepsize_coarse << endl;
        outputfile << endl;

		outputfile << "##low-thrust solver parameters" << endl;
		outputfile << "#number of time steps per phase" << endl;
		outputfile << "num_timesteps " << this->num_timesteps << endl;
		outputfile << "#Control coordinate system" << endl;
		outputfile << "#0: Cartesian" << endl;
		outputfile << "#1: Polar" << endl;
		outputfile << "control_coordinate_system " << this->control_coordinate_system << endl;
		outputfile << "#Initial guess control coordinate system" << endl;
		outputfile << "#0: Cartesian" << endl;
		outputfile << "#1: Polar" << endl;
		outputfile << "initial_guess_control_coordinate_system " << this->initial_guess_control_coordinate_system << endl;
		outputfile << "#Distribution from which to draw the step sizes for each phase" << endl;
        outputfile << "#0: uniform" << endl;
        outputfile << "#1: Gaussian" << endl;
        outputfile << "#2: Cauchy" << endl;
        outputfile << "step_size_distribution " << this->step_size_distribution << endl;
        outputfile << "#What scale width (Cauchy) or standard deviation (Gaussian) is used to create the step sizes" << endl;
        outputfile << "step_size_stdv_or_scale " << this->step_size_stdv_or_scale << endl;
		outputfile << "#Spiral model type" << endl;
		outputfile << "#0: Battin" << endl;
		outputfile << "#1: Edelbaum" << endl;
		outputfile << "spiral_model_type " << this->spiral_model_type << endl;
		outputfile << endl;

		outputfile << "##ephemeris data" << endl;
		outputfile << "#ephemeris source" << endl;
		outputfile << "#0: static" << endl;
		outputfile << "#1: SPICE (default to static if no SPICE file supplied for a body)" << endl;
		outputfile << "ephemeris_source " << this->ephemeris_source << endl;
		outputfile << "#Universe folder" << endl;
		outputfile << "universe_folder " << this->universe_folder << endl;
		outputfile << "#SPICE leap seconds kernel - required for SPICE to work" << endl;
		outputfile << "SPICE_leap_seconds_kernel " << this->SPICE_leap_seconds_kernel << endl;
		outputfile << "#SPICE_reference_frame_kernel" << endl;
		outputfile << "SPICE_reference_frame_kernel " << this->SPICE_reference_frame_kernel << endl;
		outputfile << endl;

		outputfile << "##vehicle parameters" << endl;
		outputfile << "#the maximum possible mass in kg of the spacecraft (negative number means use LV max)" << endl;
		outputfile << "maximum_mass " << this->maximum_mass << endl;
		outputfile << "#Launch vehicle type" << endl;
		outputfile << "#-2: custom launch vehicle" << endl;
		outputfile << "#-1: burn with departure stage engine" << endl;
		outputfile << "#0: fixed initial mass" << endl;
		outputfile << "#1: Atlas V (401)	 NLSII" << endl;
		outputfile << "#2: Atlas V (411)	NLSII" << endl;
		outputfile << "#3: Atlas V (421)	NLSII" << endl;
		outputfile << "#4: Atlas V (431)	NLSII" << endl;
		outputfile << "#5: Atlas V (501)	NLSII" << endl;
		outputfile << "#6: Atlas V (511)	NLSII" << endl;
		outputfile << "#7: Atlas V (521)	NLSII" << endl;
		outputfile << "#8: Atlas V (531)	NLSII" << endl;
		outputfile << "#9: Atlas V (541)	NLSII" << endl;
		outputfile << "#10: Atlas V (551)	NLSII" << endl;
		outputfile << "#11: Falcon 9 (v1.0)	NLSII" << endl;
		outputfile << "#12: Falcon 9 (v1.1)	NLSII" << endl;
		outputfile << "#13: Atlas V (551) w/Star 48	NLSI" << endl;
		outputfile << "#14: Falcon 9 Heavy" << endl;
		outputfile << "#15: Delta IV Heavy	NLSI" << endl;
		outputfile << "#16: SLS Block 1" << endl;
		outputfile << "LV_type " << this->LV_type << endl;
		outputfile << "#Launch vehicle margin (0.0 - 1.0)" << endl;
		outputfile << "LV_margin " << this->LV_margin << endl;
		outputfile << "#Launch vehicle adapter mass (kg)" << endl;
		outputfile << "LV_adapter_mass " << this->LV_adapter_mass << endl;
		outputfile << "#Custom launch vehicle coefficients (must enter 6 coefficients)" << endl;
		outputfile << "#as in a1*C3^5 + a2*C3^4 + a3*C3^3 + a4*C3^2 + a5*C3 + a6" << endl;
		outputfile << "custom_LV_coefficients";
		for (int k = 0; k < 6; ++k)
			outputfile << " " << this->custom_LV_coefficients[k];
		outputfile << endl;
		outputfile << "#Custom launch vehicle C3 bounds (two values)" << endl;
		outputfile << "custom_LV_C3_bounds " << this->custom_LV_C3_bounds[0] << " " << this->custom_LV_C3_bounds[1] << endl;
		outputfile << "#Parking orbit inclination (for use with ''depart from parking orbit'' launch vehicle option or for outputing GMAT scenarios)" << endl;
        outputfile << "parking_orbit_inclination " << this->parking_orbit_inclination << endl;
        outputfile << "#Parking orbit altitude (for use with ''depart from parking orbit'' launch vehicle option or for outputing GMAT scenarios)" << endl;
        outputfile << "parking_orbit_altitude " << this->parking_orbit_altitude << endl;
		outputfile << endl;

		outputfile << "##parameters that are only relevant for missions that use chemical propulsion" << endl;
		outputfile << "##dummy values should be used if the mission does not use chemical propulsion but are not strictly necessary" << endl;
		outputfile << "#specific impulse in seconds of the engine used for impulsive maneuvers" << endl;
		outputfile << "IspChem " << this->IspChem << endl;
		outputfile << endl;

		outputfile << "##parameters that are only relevant for missions that use a chemical EDS" << endl;
		outputfile << "##dummy values should be used if the mission does not use a chemical EDS but are not strictly necessary" << endl;
		outputfile << "#specific impulse in seconds for the earth departure stage, if applicable" << endl;
		outputfile << "IspDS " << this->IspDS << endl;
		outputfile << endl;

		outputfile << "##parameters that are only relevant for missions that use low-thrust" << endl;
		outputfile << "##dummy values should be used if the mission does not use low-thrust but are not strictly necessary" << endl;
		outputfile << "#specific impulse in seconds of the engine used for low-thrust maneuvers." << endl;
		outputfile << "#for VSI systems, this represents maximum Isp" << endl;
		outputfile << "IspLT " << this->IspLT << endl;
		outputfile << "#minimum Isp for VSI systems" << endl;
		outputfile << "IspLT_minimum " << this->IspLT_minimum << endl;
		outputfile << "#thrust of the spacecraft's low-thrust motor, in Newtons" << endl;
		outputfile << "Thrust " << this->Thrust << endl;
		outputfile << "#low-thrust engine type" << endl;
		outputfile << "#0: fixed thrust/Isp" << endl;
		outputfile << "#1: constant Isp, efficiency, EMTG computes input power" << endl;
		outputfile << "#2: choice of power model, constant efficiency, EMTG chooses Isp" << endl;
		outputfile << "#3: choice of power model, constant efficiency and Isp" << endl;
		outputfile << "#4: continuously-varying specific impulse (constant efficiency)" << endl;
		outputfile << "#5: custom thrust and mass flow rate polynomial" << endl;
		outputfile << "#6: NSTAR" << endl;
		outputfile << "#7: XIPS-25" << endl;
		outputfile << "#8: BPT-4000 High-Isp" << endl;
		outputfile << "#9: BPT-4000 High-Thrust" << endl;
		outputfile << "#10: BPT-4000 Ex-High-Isp" << endl;
		outputfile << "#11: NEXT high-Isp Phase 1" << endl;
		outputfile << "#12: VASIMR (argon, using analytical model)" << endl;
		outputfile << "#13: Hall Thruster (Xenon, using analytical model)" << endl;
		outputfile << "#14: NEXT high-ISP v10" << endl;
        outputfile << "#15: NEXT high-thrust v10" << endl;
        outputfile << "#16: BPT-4000 MALTO" << endl;
		outputfile << "#17: NEXIS Cardiff 8-15-2013" << endl;
		outputfile << "#18: H6MS Cardiff 8-15-2013" << endl;
		outputfile << "#19: BHT20K Cardiff 8-16-2013" << endl;
		outputfile << "#20: Aerojet HiVHAC EM" << endl;
		outputfile << "engine_type " << this->engine_type << endl;
		outputfile << "#Custom engine thrust coefficients (T = A + BP + C*P^2 + D*P^3 + E*P^4 + G*P^5 + H*P^6)" << endl;
		outputfile << "engine_input_thrust_coefficients";
		for (int k = 0; k < 7; ++k)
			outputfile << " " << this->engine_input_thrust_coefficients[k];
		outputfile << endl;
		outputfile << "#Custom engine mass flow rate coefficients (mdot = A + BP + C*P^2 + D*P^3 + E*P^4 + G*P^5 + H*P^6)" << endl;
		outputfile << "engine_input_mass_flow_rate_coefficients";
		for (int k = 0; k < 7; ++k)
			outputfile << " " << this->engine_input_mass_flow_rate_coefficients[k];
		outputfile << endl;
		outputfile << "#Custom engine lower and upper bounds on input power (per engine, in kW)" << endl;
		outputfile << "engine_input_power_bounds " << this->engine_input_power_bounds[0] << " " << this->engine_input_power_bounds[1] << endl;
		outputfile << "#Custom engine input efficiency" << endl;
		outputfile << "user_defined_engine_efficiency " << this->user_defined_engine_efficiency << endl;
		outputfile << "#number of low-thrust engines" << endl;
		outputfile << "number_of_engines " << this->number_of_engines << endl;
		outputfile << "#Throttle logic mode" << endl;
		outputfile << "#0: maximum power use" << endl;
		outputfile << "#1: maximum thrust" << endl;
		outputfile << "#2: maximum Isp" << endl;
		outputfile << "#3: maximum efficiency" << endl;
		outputfile << "throttle_logic_mode " << this->throttle_logic_mode << endl;
		outputfile << "#Throttle sharpness (higher means more precise, lower means smoother)" << endl;
		outputfile << "throttle_sharpness " << throttle_sharpness << endl;
		outputfile << "#engine duty cycle [0,1]" << endl;
		outputfile << "engine_duty_cycle " << this->engine_duty_cycle << endl;
		outputfile << "#electrical power available at 1 AU (kW)" << endl;
		outputfile << "power_at_1_AU " << this->power_at_1_AU << endl;
		outputfile << "#power source type, 0: solar, 1: radioisotope (or other fixed power)" << endl;
		outputfile << "power_source_type " << this->power_source_type << endl;
		outputfile << "#solar power coefficients gamma_1 through gamma_5" << endl;
		outputfile << "#if all gamma = [1 0 0 0 0], then solar power is a simple 1/r^2" << endl;
		outputfile << "solar_power_gamma";
		for (int k = 0; k < 5; ++k)
			outputfile << " " << this->solar_power_gamma[k];
		outputfile << endl;
		outputfile << "#Power margin (for thrusters, as a fraction)" << endl;
		outputfile << "power_margin " << this->power_margin << endl;
		outputfile << "#Power system decay rate (per year)" << endl;
		outputfile << "power_decay_rate " << this->power_decay_rate << endl;
		outputfile << "#spacecraft power coefficients A, B, and C" << endl;
		outputfile << "#represent the power requirements of the spacecraft at a distance r from the sun" << endl;
		outputfile << "#i.e. heaters, communications, etc" << endl;
		outputfile << "spacecraft_power_coefficients";
		for (int k = 0; k < 3; ++k)
			outputfile << " " << this->spacecraft_power_coefficients[k];
		outputfile << endl;
		outputfile << "#spacecraft power model type" << endl;
		outputfile << "#0: P_sc = A + B/r + C/r^2" << endl;
		outputfile << "#1: P_sc = A if P > A, A + B(C - P) otherwise" << endl;
		outputfile << "spacecraft_power_model_type " << this->spacecraft_power_model_type << endl;
		outputfile << "#low-thrust propulsion stage dry mass in kg, will be subtracted before chemical arrival or mid-flight switchover to chemical propulsion" << endl;
		outputfile << "EP_dry_mass " << this->EP_dry_mass << endl;
		outputfile << "#Allow initial mass to vary, up to maximum possible mass? (only relevant for MGALT and FBLT)" << endl;
		outputfile << "allow_initial_mass_to_vary " << allow_initial_mass_to_vary << endl;
		outputfile << "#Minimum dry mass" << endl;
		outputfile << "minimum_dry_mass " << this->minimum_dry_mass << endl;
		outputfile << "#Enable maximum propellant mass constraint?" << endl;
		outputfile << "enable_maximum_propellant_mass_constraint " << (int) this->enable_maximum_propellant_mass_constraint << endl;
		outputfile << "#Maximum propellant mass (kg)" << endl;
		outputfile << "maximum_propellant_mass " << this->maximum_propellant_mass << endl;
		outputfile << "#Post-mission delta-v, in km/s (alternatively defined as delta-v margin)" << endl;
		outputfile << "post_mission_delta_v " << this->post_mission_delta_v << endl;
		outputfile << "#Isp used to compute propellant for post-mission delta-v, in seconds" << endl;
		outputfile << "post_mission_Isp " << this->post_mission_Isp << endl;
		outputfile << "#Propellant margin, as a fraction of nominal propellant load" << endl;
		outputfile << "propellant_margin " << this->propellant_margin << endl;
		outputfile << endl;

		outputfile << "##Global problem settings" << endl;
		outputfile << "#mission name" << endl;
		outputfile << "mission_name " << this->mission_name << endl;
		outputfile << "#mission type - you can specify MGA, MGA-DSM, MGA-LT, or allow the outer-loop to choose" << endl;
		outputfile << "#0: MGA" << endl;
		outputfile << "#1: MGA-DSM" << endl;
		outputfile << "#2: MGA-LT" << endl;
		outputfile << "#3: FBLT" << endl;
		outputfile << "#4: MGA-NDSM" << endl;
		outputfile << "#5: DTLT" << endl;
		outputfile << "#6: solver chooses (MGA, MGA-DSM)" << endl;
		outputfile << "#7: solver chooses (MGA, MGA-LT)" << endl;
		outputfile << "#8: solver chooses (MGA-DSM, MGA-LT)" << endl;
		outputfile << "#9: solver chooses (MGA, MGA-DSM, MGA-LT)" << endl;
		outputfile << "mission_type " << this->mission_type << endl;
		outputfile << "#number of journeys (user-defined endpoints)" << endl;
		outputfile << "#Each journey has a central body and two boundary points" << endl;
		outputfile << "#Each central body has a menu of destinations which is used to choose the boundary points. Every menu is structured:" << endl;
		outputfile << "#-1: Boundary at a point in space, either fixed or free" << endl;
		outputfile << "#0: Nothing happens. This code is only used to signify 'no flyby' and should NEVER be coded as a destination." << endl;
		outputfile << "#1: Body 1 (i.e. Mercury, Io, etc)" << endl;
		outputfile << "#2: Body 2 (i.e. Venus, Europa, etc)" << endl;
		outputfile << "#..." << endl;
		outputfile << "#N: Body N " << endl;
		outputfile << "number_of_journeys " << this->number_of_journeys << endl;
		outputfile << "#maximum number of phases allowed per journey" << endl;
		outputfile << "max_phases_per_journey " << this->max_phases_per_journey << endl;
		outputfile << "#destination list (number of journeys + 1)" << endl;
		outputfile << "destination_list";
			for (int j=0; j < this->number_of_journeys; ++j)
				outputfile << " " << this->destination_list[j][0] << " " << this->destination_list[j][1];
		outputfile << endl;
		outputfile << "#the following option is relevant only if optimizing over total deltaV, should the initial impulse be included in the cost?" << endl;
		outputfile << "include_initial_impulse_in_cost " << this->include_initial_impulse_in_cost << endl;
		outputfile << "#global time bounds" << endl;
		outputfile << "#0: unbounded" << endl;
		outputfile << "#1: bounded total time (note that the global arrival date bound is by definition the same as the last journey's arrival date bound and is not duplicated" << endl;
		outputfile << "global_timebounded " << this->global_timebounded << endl;
		outputfile << "#MJD of the opening of the launch window" << endl;
		outputfile << "launch_window_open_date " << this->launch_window_open_date / 86400.0 << endl;
		outputfile << "#total flight time bounds, in days" << endl;
		outputfile << "total_flight_time_bounds " << this->total_flight_time_bounds[0] / 86400.0 << " " << this->total_flight_time_bounds[1] / 86400.0  << endl;
		outputfile << "#objective function type" << endl;
		outputfile << "#0: minimum deltaV" << endl;
		outputfile << "#1: minimum flight time" << endl;
		outputfile << "#2: maximum final mass" << endl;
		outputfile << "#3: GTOC 1 asteroid deflection function" << endl;
		outputfile << "#4: launch as late as possible in the window" << endl;
		outputfile << "#5: launch as early as possible in the window" << endl;
		outputfile << "#6: maximize orbit energy" << endl;
		outputfile << "#7: minimize launch mass" << endl;
        outputfile << "#8: arrive as early as possible" << endl;
        outputfile << "#9: arrive as late as possible" << endl;
		outputfile << "#10: minimum propellant (not the same as #2)" << endl;
		outputfile << "#11: maximum dry/wet ratio" << endl;
		outputfile << "#12: maximum arrival kinetic energy" << endl;
		outputfile << "#13: minimum BOL power" << endl;
		outputfile << "objective_type " << this->objective_type << endl;
		outputfile << "#bounds on the DLA, in degrees (typically set to declination of your launch site)" << endl;
		outputfile << "DLA_bounds " << this->DLA_bounds[0] << " " << this->DLA_bounds[1] << endl;
		outputfile << "#Initial V-Infinity vector (set to zeros unless starting the mission from periapse of a hyperbolic arrival)" << endl;
		outputfile << "initial_V_infinity " << this->initial_V_infinity[0] << " " << this->initial_V_infinity[1] << " " << this->initial_V_infinity[2] << endl;
		outputfile << "#Forced post-launch cost (in days, to be enforced after launch)" << endl;
		outputfile << "forced_post_launch_coast " << this->forced_post_launch_coast/86400.0 << endl;
		outputfile << "#Forced post flyby/intercept coast (in days, to be enforced before/after each flyby/intercept)" << endl;
		outputfile << "forced_flyby_coast " << this->forced_flyby_coast/86400.0 << endl;
		outputfile << endl;

		outputfile << "##Settings for each journey" << endl;
		outputfile << "##dummy values should be used - they should not be necessary but testing was not exhaustive so please use them" << endl;
		outputfile << "#names for each journey" << endl;
		outputfile << "journey_names";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_names[j];
		outputfile << endl;
		outputfile << "#How much mass to add to the spacecraft at the beginning of the journey (a negative number indicates a mass drop)" << endl;
		outputfile << "journey_starting_mass_increment";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_starting_mass_increment[j];
		outputfile << endl;
		outputfile << "#Is the mass increment variable (i.e. can the optimizer choose how much mass to add)" << endl;
		outputfile << "#This option is ignored for journeys with zero or negative mass increment" << endl;
		outputfile << "journey_variable_mass_increment";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_variable_mass_increment[j];
		outputfile << endl;
		outputfile << "#is each journey time bounded (one value per journey)" << endl;
		outputfile << "#0: unbounded" << endl;
		outputfile << "#1: bounded flight time" << endl;
		outputfile << "#2: bounded arrival date" << endl;
		outputfile << "journey_timebounded";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_timebounded[j];
		outputfile << endl;
		outputfile << "#what are the wait time lower and upper bounds, in days, for each journey (two numbers per journey)" << endl;
		outputfile << "journey_wait_time_bounds";
			for (int j=0; j < this->number_of_journeys; ++j)
				outputfile << " " << this->journey_wait_time_bounds[j][0] / 86400.0  << " " << this->journey_wait_time_bounds[j][1] / 86400.0 ;
		outputfile << endl;
		outputfile << "#what are the flight time bounds for each journey (two numbers per journey, use dummy values if no flight time bounds)" << endl;
		outputfile << "journey_flight_time_bounds";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_flight_time_bounds[j][0] / 86400.0  << " " << this->journey_flight_time_bounds[j][1] / 86400.0 ;
		outputfile << endl;
		outputfile << "#what are the arrival date bounds for each journey (two numbers per journey, use dummy values if no flight time bounds)" << endl;
		outputfile << "journey_arrival_date_bounds";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_arrival_date_bounds[j][0] / 86400.0  << " " << this->journey_arrival_date_bounds[j][1] / 86400.0 ;
		outputfile << endl;
		outputfile << "#what are the bounds on the initial impulse for each journey in km/s (two numbers per journey)" << endl;
		outputfile << "#you can set a very high upper bound if you are using a launchy vehicle model - the optimizer will find the correct value" << endl;
		outputfile << "journey_initial_impulse_bounds";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_initial_impulse_bounds[j][0] << " " << this->journey_initial_impulse_bounds[j][1];
		outputfile << endl;
		outputfile << "#journey departure type (one value per journey)" << endl;
		outputfile << "#0: launch or direct insertion" << endl;
		outputfile << "#1: depart from parking orbit (you can use this one in place of a launch vehicle model, and the departure burn will be done with the EDS motor)" << endl;
		outputfile << "#2: 'free' direct departure, i.e. do not burn to get the departure v_infinity (used for when operations about a small body are not modeled but the departure velocity is known)" << endl;
		outputfile << "#3: flyby (only valid for successive journeys)" << endl;
		outputfile << "#4: flyby with fixed v-infinity-out (only valid for successive journeys)" << endl;
		outputfile << "#5: spiral-out from circular orbit (low-thrust missions only)" << endl;
		outputfile << "#6: zero-turn flyby (for small bodies)" << endl;
		outputfile << "journey_departure_type";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_departure_type[j];
		outputfile << endl;
		outputfile << "#journey arrival type (one value per journey)" << endl;
		outputfile << "#0: insertion into parking orbit (use chemical Isp)" << endl;
		outputfile << "#1: rendezvous (use chemical Isp)" << endl;
		outputfile << "#2: intercept with bounded V_infinity" << endl;
		outputfile << "#3: low-thrust rendezvous (does not work if terminal phase is not low-thrust)" << endl;
		outputfile << "#4: match final v-infinity vector" << endl;
		outputfile << "#5: match final v-infinity vector (low-thrust)" << endl;
		outputfile << "#6: escape (E = 0)" << endl;
		outputfile << "#7: capture spiral" << endl;
		outputfile << "journey_arrival_type";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_arrival_type[j];
		outputfile << endl;
		outputfile << "#type of orbit elements specified at beginning of journey(0: inertial, 1: COE)" << endl;
		outputfile << "journey_departure_elements_type";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_departure_elements_type[j];
		outputfile << endl;
		outputfile << "#orbit elements at beginning of journey (a(km), e, i, RAAN, AOP, MA) supply angles in degrees" << endl;
		outputfile << "journey_departure_elements";
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			if (this->journey_departure_elements_type[j])
			{
				outputfile << " " << this->journey_departure_elements[j][0];
				outputfile << " " << this->journey_departure_elements[j][1];
				for (int k = 2; k < 6; ++k)
					outputfile << " " << this->journey_departure_elements[j][k] * 180.0 / math::PI;
			}
			else
			{
				for (int k = 0; k < 6; ++k)
					outputfile << " " << this->journey_departure_elements[j][k];
			}
		}
		outputfile << endl;
		outputfile << "#Vary journey departure elements? (one entry per element per journey: 0 means no, 1 means yes)" << endl;
        outputfile << "journey_departure_elements_vary_flag";
        for (int j = 0; j < this->number_of_journeys; ++j)
		{
			for (int k = 0; k < 6; ++k)
				outputfile << " " << this->journey_departure_elements_vary_flag[j][k];
		}
        outputfile << endl;
		outputfile << "#Lower and upper bounds on journey departure elements (two per element per journey, ignored if vary flag is off for that element)" << endl;
        outputfile << "journey_departure_elements_bounds";
        for (int j = 0; j < this->number_of_journeys; ++j)
		{
            if (this->journey_departure_elements_type[j])
			{
				outputfile << " " << this->journey_departure_elements_bounds[j][0][0] << " " << this->journey_departure_elements_bounds[j][0][1];
				outputfile << " " << this->journey_departure_elements_bounds[j][1][0] << " " << this->journey_departure_elements_bounds[j][1][1];
				for (int k = 2; k < 6; ++k)
					outputfile << " " << this->journey_departure_elements_bounds[j][k][0] * 180.0 / math::PI << " " << this->journey_departure_elements_bounds[j][k][1] * 180.0 / math::PI ;
			}
			else
			{
				for (int k = 0; k < 6; ++k)
					outputfile << " " << this->journey_departure_elements_bounds[j][k][0] << " " << this->journey_departure_elements_bounds[j][k][1];
			}
		}
		outputfile << endl;
		outputfile << "#type of orbit elements specified at end of journey(0: inertial, 1: COE)" << endl;
		outputfile << "journey_arrival_elements_type";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_arrival_elements_type[j];
		outputfile << endl;
		outputfile << "#orbit elements at end of journey (a(km), e, i, RAAN, AOP, MA) supply angles in degrees" << endl;
		outputfile << "journey_arrival_elements";
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			if (this->journey_arrival_elements_type[j])
			{
				outputfile << " " << this->journey_arrival_elements[j][0];
				outputfile << " " << this->journey_arrival_elements[j][1];
				for (int k = 2; k < 6; ++k)
					outputfile << " " << this->journey_arrival_elements[j][k] * 180.0 / math::PI;
			}
			else
			{
				for (int k = 0; k < 6; ++k)
					outputfile << " " << this->journey_arrival_elements[j][k];
			}
		}
		outputfile << endl;
        outputfile << "#Vary journey arrival elements? (one entry per element per journey: 0 means no, 1 means yes)" << endl;
        outputfile << "journey_arrival_elements_vary_flag";
        for (int j = 0; j < this->number_of_journeys; ++j)
		{
			for (int k = 0; k < 6; ++k)
				outputfile << " " << this->journey_arrival_elements_vary_flag[j][k];
		}
        outputfile << endl;
		outputfile << "#Lower and upper bounds on journey arrival elements (two per element per journey, ignored if vary flag is off for that element)" << endl;
        outputfile << "journey_arrival_elements_bounds";
        for (int j = 0; j < this->number_of_journeys; ++j)
		{
			if (this->journey_arrival_elements_type[j])
			{
				outputfile << " " << this->journey_arrival_elements_bounds[j][0][0] << " " << this->journey_arrival_elements_bounds[j][0][1];
				outputfile << " " << this->journey_arrival_elements_bounds[j][1][0] << " " << this->journey_arrival_elements_bounds[j][1][1];
				for (int k = 2; k < 6; ++k)
					outputfile << " " << this->journey_arrival_elements_bounds[j][k][0] * 180.0 / math::PI << " " << this->journey_arrival_elements_bounds[j][k][1] * 180.0 / math::PI ;
			}
			else
			{
				for (int k = 0; k < 6; ++k)
					outputfile << " " << this->journey_arrival_elements_bounds[j][k][0] << " " << this->journey_arrival_elements_bounds[j][k][1];
			}
		}
		outputfile << endl;
		outputfile << "#journey central body" << endl;
		outputfile << "#use SPICE names, as per http://www-int.stsci.edu/~sontag/spicedocs/req/naif_ids.html" << endl;
		outputfile << "journey_central_body";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_central_body[j];
		outputfile << endl;
		outputfile << "#final VHP for journeys that end in intercepts, in km/s (three numbers per journey)" << endl;
		outputfile << "journey_final_velocity";
		for (int j=0; j < this->number_of_journeys; ++j)
		{
			for (int k = 0; k < 3; ++k)
				outputfile << " " << this->journey_final_velocity[j][k];
		}
		outputfile << endl;
		outputfile << "#Impose arrival declination constraint on each journey?" << endl;
		outputfile << "#0: no" << endl;
		outputfile << "#1: yes" << endl;
		outputfile << "journey_arrival_declination_constraint_flag";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_arrival_declination_constraint_flag[j];
		outputfile << endl;
		outputfile << "#Arrival declination bounds for each journey" << endl;
		outputfile << "#Two numbers per journey, in degrees" << endl;
		outputfile << "journey_arrival_declination_bounds";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_arrival_declination_bounds[j][0] * 180.0/math::PI << " " << this->journey_arrival_declination_bounds[j][1] * 180.0/math::PI;
		outputfile << endl;
		outputfile << "#Starting orbital radius for an escape spiral at the beginning of the journey" << endl;
		outputfile << "journey_escape_spiral_starting_radius";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_escape_spiral_starting_radius[j];
		outputfile << endl;
		outputfile << "#Final orbital radius for a capture spiral at the end of the journey" << endl;
		outputfile << "journey_capture_spiral_final_radius";
		for (int j=0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_capture_spiral_final_radius[j];
		outputfile << endl;
		outputfile << endl;

		outputfile << "##Perturbation settings" << endl;
		outputfile << "#Enable solar radiation pressure?" << endl;
		outputfile << "perturb_SRP " << this->perturb_SRP << endl;
		outputfile << "#Enable third-body perturbations?" << endl;
		outputfile << "perturb_thirdbody " << this->perturb_thirdbody << endl;
		outputfile << "#Journey perturbation bodies. One line per journey. The numbers in the line correspond to" << endl;
		outputfile << "#bodies in the journey's Universe file. If perturbations are off, each line should just have a zero" << endl;
		outputfile << "#the numbers in the first line are the number of perturbation bodies for each journey" << endl;
		outputfile << "journey_perturbation_bodies";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->journey_number_of_perturbation_bodies[j];
		outputfile << endl;
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			outputfile << this->journey_perturbation_bodies[j][0];
			for (int b = 1; b < this->journey_number_of_perturbation_bodies[j]; ++b)
				outputfile << " " << this->journey_perturbation_bodies[j][b];
			outputfile << endl;
		}
		outputfile << "#end_journey_perturbation_bodies" << endl;
		outputfile << "#Spacecraft area (in m^2)" << endl;
		outputfile << "spacecraft_area " << this->spacecraft_area << endl;
		outputfile << "#Coefficient of reflectivity" << endl;
		outputfile << "#0.0: perfectly translucent" << endl;
		outputfile << "#1.0: perfectly absorbing" << endl;
		outputfile << "#2.0: perfectly reflecting" << endl;
		outputfile << "coefficient_of_reflectivity " << this->coefficient_of_reflectivity << endl;
		outputfile << endl;
		outputfile << endl;

		outputfile << "##Outer-loop selectable options settings" << endl;
		outputfile << "#Allow outer-loop to vary power level?" << endl;
		outputfile << "outerloop_vary_power " << this->outerloop_vary_power << endl;
		outputfile << "#Allow outer-loop to vary launch epoch?" << endl;
		outputfile << "outerloop_vary_launch_epoch " << this->outerloop_vary_launch_epoch << endl;
		outputfile << "#Allow outer-loop to vary flight time upper bound?" << endl;
		outputfile << "outerloop_vary_flight_time_upper_bound " << this->outerloop_vary_flight_time_upper_bound << endl;
		outputfile << "#Restrict flight-time lower bound when running outer-loop?" << endl;
		outputfile << "outerloop_restrict_flight_time_lower_bound " << this->outerloop_restrict_flight_time_lower_bound << endl;
		outputfile << "#Allow outer-loop to vary thruster type?" << endl;
		outputfile << "outerloop_vary_thruster_type " << this->outerloop_vary_thruster_type << endl;
		outputfile << "#Allow outer-loop to vary number of thrusters?" << endl;
		outputfile << "outerloop_vary_number_of_thrusters " << this->outerloop_vary_number_of_thrusters << endl;
		outputfile << "#Allow outer-loop to vary launch vehicle?" << endl;
		outputfile << "outerloop_vary_launch_vehicle " << this->outerloop_vary_launch_vehicle << endl;
		outputfile << "#Allow outer-loop to vary first journey departure C3?" << endl;
		outputfile << "outerloop_vary_departure_C3 " << this->outerloop_vary_departure_C3 << endl;
		outputfile << "#Allow outer-loop to vary last journey arrival C3?" << endl;
		outputfile << "outerloop_vary_arrival_C3 " << this->outerloop_vary_arrival_C3 << endl;
		outputfile << "#Allow outer-loop to vary journey destination? (one value per journey)" << endl;
		outputfile << "outerloop_vary_journey_destination";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->outerloop_vary_journey_destination[j];
		outputfile << endl;
		outputfile << "#Allow outer-loop to vary journey flyby sequence? (one value per journey)" << endl;
		outputfile << "outerloop_vary_journey_flyby_sequence";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->outerloop_vary_journey_flyby_sequence[j];
		outputfile << endl;
		outputfile << "#Outer-loop power at 1 AU choices (in kW)" << endl;
		outputfile << "outerloop_power_choices";
		for (int entry = 0; entry < this->outerloop_power_choices.size(); ++entry)
			outputfile << " " << this->outerloop_power_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop launch window open epoch choices (in MJD)" << endl;
		outputfile << "outerloop_launch_epoch_choices";
		for (int entry = 0; entry < this->outerloop_launch_epoch_choices.size(); ++entry)
			outputfile << " " << this->outerloop_launch_epoch_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop flight time upper bound choices (in days)" << endl;
		outputfile << "outerloop_flight_time_upper_bound_choices";
		for (int entry = 0; entry < this->outerloop_flight_time_upper_bound_choices.size(); ++entry)
			outputfile << " " << this->outerloop_flight_time_upper_bound_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop thruster type choices (in order of most to least preferable)" << endl;
		outputfile << "outerloop_thruster_type_choices";
		for (int entry = 0; entry < this->outerloop_thruster_type_choices.size(); ++entry)
			outputfile << " " << this->outerloop_thruster_type_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop number of thruster choices" << endl;
		outputfile << "outerloop_number_of_thrusters_choices";
		for (int entry = 0; entry < this->outerloop_number_of_thrusters_choices.size(); ++entry)
			outputfile << " " << this->outerloop_number_of_thrusters_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop launch vehicle choices (in order of most to least preferable)" << endl;
		outputfile << "outerloop_launch_vehicle_choices";
		for (int entry = 0; entry < this->outerloop_launch_vehicle_choices.size(); ++entry)
			outputfile << " " << this->outerloop_launch_vehicle_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop first journey departure C3 choices" << endl;
		outputfile << "outerloop_departure_C3_choices";
		for (int entry = 0; entry < this->outerloop_departure_C3_choices.size(); ++entry)
			outputfile << " " << this->outerloop_departure_C3_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop last arrival departure C3 choices" << endl;
		outputfile << "outerloop_arrival_C3_choices";
		for (int entry = 0; entry < this->outerloop_arrival_C3_choices.size(); ++entry)
			outputfile << " " << this->outerloop_arrival_C3_choices[entry];
		outputfile << endl;
		outputfile << "#Outer-loop maximum number of flybys (one value for each journey)" << endl;
		outputfile << "outerloop_journey_maximum_number_of_flybys";
		for (int j = 0; j < this->number_of_journeys; ++j)
			outputfile << " " << this->outerloop_journey_maximum_number_of_flybys[j];
		outputfile << endl;
		outputfile << "#Outer-loop journey destination choices (one line for each journey)" << endl;
		outputfile << "outerloop_journey_destination_choices" << endl;
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			for (int entry = 0; entry < this->outerloop_journey_destination_choices[j].size(); ++entry)
				outputfile << " " << this->outerloop_journey_destination_choices[j][entry];
			outputfile << endl;
		}
		outputfile << "#Outer-loop flyby sequence choices (one line for each journey)" << endl;
		outputfile << "outerloop_journey_flyby_sequence_choices" << endl;
		for (int j = 0; j < this->number_of_journeys; ++j)
		{
			for (int entry = 0; entry < this->outerloop_journey_flyby_sequence_choices[j].size(); ++entry)
				outputfile << " " << this->outerloop_journey_flyby_sequence_choices[j][entry];
			outputfile << endl;
		}
		outputfile << endl;

		outputfile << "##Outer-loop objective function settings" << endl;
		outputfile << "#Pick as many as you want. The Pareto surface will be generated in these dimensions" << endl;
		outputfile << "#0: BOL power at 1 AU (kW)" << endl;
		outputfile << "#1: Launch epoch (MJD)" << endl;
		outputfile << "#2: Flight time (days)" << endl;
		outputfile << "#3: Thruster preference" << endl;
		outputfile << "#4: Number of thrusters" << endl;
		outputfile << "#5: Launch vehicle preference" << endl;
		outputfile << "#6: Delivered mass to final target (kg)" << endl;
		outputfile << "#7: Final journey mass increment (for maximizing sample return)" << endl;
		outputfile << "#8: First journey departure C3 (km^2/s^2)" << endl;
		outputfile << "#9: Final journey arrival C3 (km^2/s^2)" << endl;
		outputfile << "#10: Total delta-v (km/s)" << endl;
		outputfile << "#11: Inner-loop objective (whatever it was)" << endl;
		outputfile << "outerloop_objective_function_choices";
		for (int entry = 0; entry < this->outerloop_objective_function_choices.size(); ++entry)
			outputfile << " " << this->outerloop_objective_function_choices[entry];
		outputfile << endl;
		outputfile << endl;

		outputfile << "##output format settings" << endl;
		outputfile << "#output units, 0: km and km/s, 1: LU and LU/day" << endl;
		outputfile << "output_units " << this->output_units << endl;
		outputfile << "#Output a GMAT script (not compatible with non-body boundary conditions or thruster/power models)" << endl;
		outputfile << "create_GMAT_script " << this->create_GMAT_script << endl;
		outputfile << endl;

		outputfile << "##debug code" << endl;
		outputfile << "##the purpose of this code is so that you can turn the inner-loop solver on and off, force a sequence of planets and/or phase types" << endl;
		outputfile << "#sequence, must have (max_phases_per_journey) entries for each journey. Use 0 to encode 'no flyby'" << endl;
		outputfile << "#integer codes represent planets" << endl;
		outputfile << "#this option is NOT used if the outer-loop is turned on" << endl;
		outputfile << "#first number is the number of sequences listed, followed by the sequences" << endl;
		outputfile << "sequence " << this->number_of_trial_sequences << endl;
		for (size_t trial = 0; trial < this->sequence_input.size(); ++trial)
		{
			for (int j=0; j < this->number_of_journeys; ++j)
			{
				for (int p = 0; p < max_phases_per_journey; ++p)
					if (p < this->sequence_input[trial][j].size())
					{
						if (p > 0 || j > 0)
							outputfile << " ";
						outputfile << (this->sequence_input[trial][j][p]);
					}
					else
						outputfile << " 0";
			}
			outputfile << endl;
		}
		outputfile << "#phase type, must have one entry for each phase in the mission" << endl;
		outputfile << "#this option allows you to have different phases use different propulsion systems" << endl;
		outputfile << "#0: MGA, 1: MGA-DSM, 2: MGA-LT, 3: FBLT, 4: FBLT-S, 5: MGANDSM, 6: DTLT" << endl;
		outputfile << "#if mission_type is set to 0, 1, 2, 3, 4, 5, 6 then the following option is ignored" << endl;
		outputfile << "#if mission_type > 6 and the outer-loop is ON, then the following option is ignored" << endl;
		outputfile << "#the following option is only used if the outer-loop is OFF and mission_type > 6" << endl;
		if (this->run_outerloop == 0 && this->mission_type > 6) //if we need to specify phase types
		{
			outputfile << "phase_type";
			for (int j=0; j < this->number_of_journeys; ++j)
			{
				for (int p = 0; p < this->number_of_phases[j]; ++p)
					outputfile << " " << this->phase_type_input[j][p];
			}
			outputfile << endl;
		}
		else
		{
			outputfile << "#not specified because either the outer-loop is off or mission_type > 3" << endl;
			//outputfile << "phase_type " << mission_type << endl;
		}

		outputfile << "#Check derivatives against finite differencing?" << endl;
		outputfile << "check_derivatives " << this->check_derivatives << endl;
		outputfile << "#which inner loop solver to run?" << endl;
		outputfile << "#0: none, evaluate trialX, 1: evaluate a batch of decision vectors, 2: run MBH, 3: run constrained DE, 4: run SNOPT using trialX as initial guess" << endl;
		outputfile << "run_inner_loop " << this->run_inner_loop << endl;
		outputfile << "#trial decision vector" << endl;
		outputfile << "#trialX" << endl;
		if (trialX.size() > 0)
		{
			outputfile << "trialX" << endl;

			for (int entry = 0; entry < this->number_of_trial_sequences; ++ entry)
			{
				outputfile << this->trialX[entry][0];
				for (size_t k=1; k < this->trialX[entry].size(); ++k)
					outputfile << " " << this->trialX[entry][k];

				outputfile << endl;
			}
		}
		outputfile << endl;
		outputfile << endl;

		outputfile << "#end options file";
		outputfile.close();

		cout << "Options file written to '" << filename << "'" << endl;
		return 0;
	}

	else cout << "Unable to write to output file " << filename << endl;
	return 1;
}

void missionoptions::construct_thruster_launch_vehicle_name_arrays()
{
	this->LV_names.push_back("fm");
	this->LV_names.push_back("AV401");
	this->LV_names.push_back("AV411");
	this->LV_names.push_back("AV421");
	this->LV_names.push_back("AV431");
	this->LV_names.push_back("AV501");
	this->LV_names.push_back("AV511");
	this->LV_names.push_back("AV521");
	this->LV_names.push_back("AV531");
	this->LV_names.push_back("AV541");
	this->LV_names.push_back("AV551");
	this->LV_names.push_back("F910");
	this->LV_names.push_back("F911");
	this->LV_names.push_back("AV551s48");
	this->LV_names.push_back("F9H");
	this->LV_names.push_back("D4H");
	this->LV_names.push_back("SLSb1");

	this->thruster_names.push_back("fixedthruster");
	this->thruster_names.push_back("constIspEffPickPower");
	this->thruster_names.push_back("ChooseIsp");
	this->thruster_names.push_back("fixedIsp");
	this->thruster_names.push_back("VSI");
	this->thruster_names.push_back("customthruster");
	this->thruster_names.push_back("NSTAR");
	this->thruster_names.push_back("XIPS25");
	this->thruster_names.push_back("BPT4000HIsp");
	this->thruster_names.push_back("BPT4000Hthrust");
	this->thruster_names.push_back("BPT4000XHIsp");
	this->thruster_names.push_back("NEXTHIspv9");
	this->thruster_names.push_back("VASIMRargon");
	this->thruster_names.push_back("VSIxenonhall");
	this->thruster_names.push_back("NEXTHIspv10");
	this->thruster_names.push_back("NEXTHthrustv10");
	this->thruster_names.push_back("BPT4000MALTO");
	this->thruster_names.push_back("NEXIS");
	this->thruster_names.push_back("H6MS");
	this->thruster_names.push_back("BHT20K");
	this->thruster_names.push_back("HiVHAc");
}

} /* namespace EMTG */
