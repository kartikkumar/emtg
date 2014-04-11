/*
 * missionoptions.h
 *
 *  Created on: Jul 13, 2012
 *      Author: Jacob
 */

#ifndef MISSIONOPTIONS_H_
#define MISSIONOPTIONS_H_

#include "boost/ptr_container/ptr_vector.hpp"

#include <vector>
#include <string>

using namespace std;

namespace EMTG {

class missionoptions
{
public:
	//constructor
	missionoptions();
	missionoptions(string optionsfile);

	//destructor
	virtual ~missionoptions();

	//methods
	int parse_options_file(string filename);
	int parse_options_line(ifstream& inputfile, string& choice, double& value, char* dump_buffer, int& linenumber);
	int print_options_file(string filename);
	void construct_thruster_launch_vehicle_name_arrays();

	//check for success
	int file_status; //0: successful read, 1: error reading file, 2: no file found with this file name
	string error_message;


	//fields

	//****************************************************************************************
	//
	// PLEASE NOTE
	//
	// If you add new fields to the missionoptions class, make sure that you also add:
	// 1. code to READ values into the field in missionoptions::parse_options_file(string filename)
	// 2. code to WRITE the field out to the disk in missionoptions::print_options_file(string filename)
	//
	//****************************************************************************************

	//problem type
	//0: standard EMTG mission
	int problem_type;
	
	//physical constants
	double G;
	double TU;
	double g0;
	double AU;

	//outer loop solver settings
	bool run_outerloop; //whether or not to run the outer loop; if false
	int outerloop_popsize; //population size
	int outerloop_genmax; //number of generations
	int outerloop_tournamentsize; //tournament size for selection
	double outerloop_CR; //math::crossover ratio
	double outerloop_mu; //mutation rate
	int outerloop_stallmax; //maximum number of stall generations
	double outerloop_tolfit; //fitness tolerance
	int outerloop_ntrials; //how many times to run the outer loop
	int outerloop_elitecount; //how many elite individuals to retain
	bool outerloop_useparallel; //whether or not to use the parallel outer-loop
	int outerloop_warmstart; //if true, read "population.txt" and "solutions.txt"
	bool outerloop_reevaluate_full_population;//if true, re-evaluate the entire population each generation, otherwise read from the archive
	bool quiet_outerloop; //if true, suppress all text outputs except error catches

	//outer loop selectable options settings
	bool outerloop_vary_power;
	bool outerloop_vary_launch_epoch;
	bool outerloop_vary_flight_time_upper_bound;
	bool outerloop_vary_thruster_type;
	bool outerloop_vary_number_of_thrusters;
	bool outerloop_vary_launch_vehicle;
	bool outerloop_vary_departure_C3;
	bool outerloop_vary_arrival_C3;
	vector<bool> outerloop_vary_journey_destination;
	vector<bool> outerloop_vary_journey_flyby_sequence;
	vector<double> outerloop_power_choices;
	vector<double> outerloop_launch_epoch_choices;
	vector<double> outerloop_flight_time_upper_bound_choices;
	bool outerloop_restrict_flight_time_lower_bound;
	vector<int> outerloop_thruster_type_choices;
	vector<int> outerloop_number_of_thrusters_choices;
	vector<int> outerloop_launch_vehicle_choices;
	vector<double> outerloop_departure_C3_choices;
	vector<double> outerloop_arrival_C3_choices;
	vector< vector<int> > outerloop_journey_destination_choices;
	vector< vector<int> > outerloop_journey_flyby_sequence_choices;
	vector<int> outerloop_journey_maximum_number_of_flybys;

	//outerloop objective settings
	vector<int> outerloop_objective_function_choices;

	//inner loop solver settings
	int NLP_solver_type;
	bool NLP_solver_mode; //false: feasible point only, true: optimize
	bool quiet_NLP;
	bool ACE_feasible_point_finder;
	int MBH_max_not_improve;
	int MBH_max_trials;
	int MBH_max_run_time;
	double MBH_max_step_size;
	int MBH_hop_distribution;
	double MBH_Pareto_alpha;
	double MBH_time_hop_probability;
	double snopt_feasibility_tolerance;
	int snopt_major_iterations;
	time_t snopt_max_run_time;
	int derivative_type;
	bool seed_MBH;
	bool interpolate_initial_guess;
	int initial_guess_num_timesteps;
	int initial_guess_step_size_distribution; //0: uniform, 1: Gaussian, 2: Cauchy
	double initial_guess_step_size_stdv_or_scale;
	int MBH_zero_control_initial_guess; //0: do not use, 1: zero-control for resets, random perturbations for hops, 2: always use zero-control guess except when seeded
	bool MBH_two_step; //whether or not to use the 2-step MBH (coarse then fine derivatives)
	double FD_stepsize; //"fine" finite differencing step size
	double FD_stepsize_coarse; //"coarse" finite differencing step
	


	//problem settings set by the user

	//ephemeris data
	string SPICE_leap_seconds_kernel;
	string SPICE_reference_frame_kernel;
	string universe_folder;
	int ephemeris_source; //0: static, 1: SPICE (default to static if no SPICE file supplied for a body)

	//low thrust solver parameters
	int num_timesteps; //number of timesteps per phase
	int step_size_distribution; //0: uniform, 1: Gaussian, 2: Cauchy
	int spiral_model_type; //0: Battin, 1: Edelbaum
	double step_size_stdv_or_scale;

	//vehicle parameters
	bool allow_initial_mass_to_vary;
	double maximum_mass; //the maximum possible mass of the spacecraft (negative number means use LV max)
	double IspLT; //specific impulse of the engine used for low-thrust maneuvers
	double IspLT_minimum; //minimum specific impulse for variable specific impulse problems. Set to a default value of 3000 for now.
	double IspChem; //specific impulse of the engine used for impulsive maneuvers
	double IspDS; //specific impulse for the earth departure stage, if applicable
	double Thrust; //thrust of the spacecraft, in Newtons
	int LV_type;
	//-1: burn with departure stage engine
	//0: fixed initial mass
	//1: Falcon 9 (v1.0)			NLSII
	//2: Atlas V (501)				NLSII
	//3: Falcon 9 (v1.1)			NLSII
	//4: Atlas V (401)				NLSII
	//5: Atlas V (511)				NLSII
	//6: Atlas V (411)				NLSII
	//7: Atlas V (521)				NLSII
	//8: Atlas V (421)				NLSII
	//9: Atlas V (531)				NLSII
	//10: Atlas V (431)				NLSII
	//11: Atlas V (541)				NLSII
	//12: Atlas V (551)				NLSII
	//13: Atlas V (551) w/Star 48	NLSI
	double LV_margin;
	double LV_adapter_mass;
	double custom_LV_coefficients[6];
	double custom_LV_C3_bounds[2];
	double parking_orbit_inclination;
	double parking_orbit_altitude;
	int engine_type;
	//0: fixed thrust/Isp
	//1: constant Isp, efficiency, EMTG computes input power
	//2: constant power, EMTG chooses exhaust velocity
	//3: constant Isp, efficiency, compute thrust based on available power
	//4: continuously-varying specific impulse
	//5: custom thrust and mass flow rate polynomial
	//6: NSTAR
	//7: XIPS-25
	//8: BPT-4000 High-Isp
	//9: BPT-4000 High-Thrust
	//10: BPT-4000 Ex-High-Isp
	//11: NEXT high-Isp	
	int number_of_engines; //only relevant when modeling an engine
	double engine_duty_cycle; //percentage of time that engine can operate
	double power_at_1_AU; //in kW
	int power_source_type; //0: solar, 1: radioisotope (or other fixed power)
	double solar_power_gamma[5]; //coefficients for solar panels
	double power_margin;
	double power_decay_rate; //in percent per year
	int throttle_logic_mode;
	double throttle_sharpness;
	//0: maximum power use
	//1: maximum thrust
	//2: minimum mass-flow rate
	double spacecraft_power_coefficients[3];
	double engine_input_thrust_coefficients[7];
	double engine_input_mass_flow_rate_coefficients[7];
	double engine_input_power_bounds[2];
	double user_defined_engine_efficiency;
	int spacecraft_power_model_type;
	double EP_dry_mass; //in kg

	//terminal constraints
	double minimum_dry_mass; //in kg
	double post_mission_delta_v; //in km/s
	double post_mission_Isp; //in s, Isp for post_mission_delta_v
	double propellant_margin; //percent of propellant to be held as margin over the nominal propellant load

	//perturbation settings
	bool perturb_SRP; //solar radiation pressure?
	bool perturb_thirdbody; //third body perturbations?
	vector<int> journey_number_of_perturbation_bodies;
	vector< vector<int> > journey_perturbation_bodies;
	double spacecraft_area; //in m^2
	double coefficient_of_reflectivity;

	//global problem settings
	int number_of_journeys;
	int max_phases_per_journey;
	vector< vector <int> > destination_list;
	bool include_initial_impulse_in_cost;
	bool global_timebounded;//0: unbounded, 1: bounded total time (note that the global arrival date bound is by definition the same as the last journey's arrival date bound and is not duplicated
	double launch_window_open_date;//MJD
	double total_flight_time_bounds[2];//days
	int objective_type; //0: minimum deltaV, 1: minimum time, //2: maximum final mass
	vector<double> DLA_bounds; //DLA in degrees
	string mission_name;
	int mission_type; //0: MGA, 1: MGA-DSM; 2: MGA-LT, 3: FBLT, 4: FBLT-S, 5: MGA-NDSM, 6: NSLT, 7: solver chooses (MGA, MGA-DSM), 8: solver chooses (MGA, MGA-LT), 9: solver chooses (MGA-DSM, MGA-LT), 10: solver chooses (MGA, MGA-DSM, MGA-LT)
	vector<double> initial_V_infinity; //initial V-infinity for missions which start at periapse of a hyperbolic arrival
	double forced_post_launch_coast; //in days, to be enforced after launch
	double forced_flyby_coast; //in days, to be enforced before/after each flyby/intercept

	//settings for each journey
	vector<string> journey_names;
	vector<int> journey_timebounded;//0: unbounded, 1: bounded flight time, 2: bounded arrival date
	vector< vector<double> > journey_wait_time_bounds;//days
	vector< vector<double> >  journey_flight_time_bounds;
	vector< vector<double> >  journey_arrival_date_bounds;
	vector< vector<double> >  journey_initial_impulse_bounds; //in km/s
	vector<int> journey_arrival_type; //0: orbit insertion (use chemical Isp), 1: rendezvous (use chemical Isp), 2: flyby with bounded VHP, 3: low-thrust rendezvous (does not work if terminal phase is not low-thrust), 4: match final v-infinity vector
	vector<int> journey_departure_type; //0: launch or direct insertion, 1: depart from parking orbit (you can use this one in place of a launch vehicle model, and the departure burn will be done with the EDS motor), 2: 'free' direct departure, i.e. do not burn to get the departure v_infinity, 3: Start from Sphere of Influence (use SOI angles chosen by previous journey's endpoint, i.e. after a spiral-out or fully modeled departure from parking orbit)
	vector<int> journey_arrival_elements_type; //0: cartesian, 1: COE
	vector< vector <double> > journey_arrival_elements; //a(km), e, i, RAAN, AOP, TA
    vector< vector < vector <double> > > journey_arrival_elements_bounds;
    vector< vector < bool > > journey_arrival_elements_vary_flag;
	vector<int> journey_departure_elements_type; //0: cartesian, 1: COE
	vector< vector <double> > journey_departure_elements; //a(km), e, i, RAAN, AOP, TA
    vector< vector < vector <double> > > journey_departure_elements_bounds;
    vector< vector < bool > > journey_departure_elements_vary_flag;
	vector<string> journey_central_body; //spice names
	vector< vector<double> > journey_initial_velocity; //in km/s
	vector< vector<double> > journey_final_velocity; //in km/s
	vector< double > journey_starting_mass_increment; //how much mass to add to the spacecraft at the beginning of the journey (a negative number indicates a mass drop)
	vector< bool > journey_variable_mass_increment;
	vector< bool > journey_arrival_declination_constraint_flag;
	vector< vector <double> > journey_arrival_declination_bounds; //in degrees
	vector<double> journey_escape_spiral_starting_radius; //in km
	vector<double> journey_capture_spiral_final_radius; //in km

	//output format settings
	int output_units; //0: km and km/s, 1: LU and LU/day
	int create_GMAT_script; //0: no, 1: yes

	//debug code
	int run_inner_loop;
	vector< vector<double> > trialX;
	vector<double> current_trialX;
	vector< vector<int> > phase_type_input;
	int number_of_trial_sequences;
	vector< vector< vector<int> > > sequence_input;
	bool check_derivatives;

	//problem settings which will be automatically generated
	vector< vector<int> > sequence;
	vector< vector<int> > phase_type;
	vector< vector<int> > number_of_phases_input;
	vector<int> number_of_phases;
	int total_number_of_phases;
	bool quiet_basinhopping;
	string outputfile;
	string GMAT_outputfile;
	string description;
	string convergence_file;
	string population_file;
	string solutions_file;
	string working_directory;
	vector<int> iGfun;
	vector<int> jGvar;
	vector<double> X_scale_ranges;
	vector<string> LV_names;
	vector<string> thruster_names;
};

} /* namespace EMTG */

#endif /* MISSIONOPTIONS_H_ */