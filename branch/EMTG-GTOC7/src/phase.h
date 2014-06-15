/*
 * phase.h
 *
 *  Created on: Jul 15, 2012
 *      Author: Jacob
 */



#include "missionoptions.h"
#include "universe.h"
#include "bplane.h"
#include "STM.h"


#include <string>
#include <vector>
#include <fstream>

#ifndef PHASE_H_
#define PHASE_H_

namespace EMTG {
	class phase
	{
	public:
		//constructor
		phase();

		//destructor
		virtual ~phase();

		//methods
		int locate_boundary_point(int location,
								  int boundary_type,
								  bool left_boundary, 
								  EMTG::Astrodynamics::universe* Universe, 
								  double* state, 
								  double* X_infinity, 
								  double epoch, 
								  double* X, 
								  int* Xindex, 
								  double* F, 
								  int* Findex, 
								  double* G, 
								  int* Gindex, 
								  const int& needG, 
								  int j, 
								  int p,
								  missionoptions* options);

		double process_arrival(double* incoming_velocity,
							   double* boundary_state, 
							   double* X_infinity, 
							   double* current_epoch,
							   double mu, 
							   double r_SOI, 
							   double* F, 
							   int* Findex, 
							   int j, 
							   missionoptions* options, 
							   EMTG::Astrodynamics::universe* Universe);

		virtual math::Matrix<double> calculate_flyby_periapse_state(math::Matrix<double>& Vinf_in,
																	math::Matrix<double>& Vinf_out, 
																	const double& flyby_altitude, 
																	EMTG::Astrodynamics::body& TheBody);

		math::Matrix<double> calculate_periapse_state_from_asymptote_and_parking_orbit(math::Matrix<double>& V_infinity,
																					   double parking_orbit_incination,
																					   double parking_orbit_altitude,
																					   double &epoch,
																					   EMTG::Astrodynamics::universe* Universe,
																					   EMTG::Astrodynamics::body &TheBody);

		double compute_timestep_width_from_distribution(double step,
														missionoptions* options,
														double& scale_or_stdv);

		void write_summary_line(missionoptions* options,
								EMTG::Astrodynamics::universe* Universe,
								int* eventcount,
								double current_epoch_MJD,
								string event_type,
								string boundary_name,
								double timestep_size,
								double flyby_altitude,
								double b_plane_insertion_angle,
								double flyby_turn_angle,
								double angle1,
								double angle2,
								double C3,
								double* state,
								double* dV,
								double* ThrustVector,
								double dVmag,
								double Thrust,
								double Isp,
								double AvailPower,
								double mdot,
								int number_of_active_engines,
								double active_power);

		void process_left_boundary_condition(double* X,
											 int* Xindex,
											 double* F,
											 int* Findex,
											 double* G,
											 int* Gindex,
											 const int& needG,
											 double* current_epoch,
											 double* current_state,
											 double* current_deltaV,
											 double* boundary1_state,
											 double* boundary2_state,
											 int j,
											 int p,
											 EMTG::Astrodynamics::universe* Universe,
											 missionoptions* options);

		void process_right_boundary_condition(double* X,
												int* Xindex,
												double* F,
												int* Findex,
												double* G,
												int* Gindex,
												const int& needG,
												double* current_epoch,
												double* current_state,
												double* current_deltaV,
												double* boundary1_state,
												double* boundary2_state,
												int j,
												int p,
												EMTG::Astrodynamics::universe* Universe,
												missionoptions* options);

		void calcbounds_flight_time(const string& prefix,
									 int first_X_entry_in_phase,
									 vector<double>* Xupperbounds,
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
									 vector<double>* synodic_periods,
									 int j,
									 int p,
									 EMTG::Astrodynamics::universe* Universe,
									 missionoptions* options);

		void calcbounds_left_boundary(const string& prefix,
										 int first_X_entry_in_phase,
										 vector<double>* Xupperbounds,
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
										 int j,
										 int p,
										 EMTG::Astrodynamics::universe* Universe,
										 missionoptions* options);

		void calcbounds_right_boundary(const string& prefix,
										 int first_X_entry_in_phase,
										 vector<double>* Xupperbounds,
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
										 int j,
										 int p,
										 EMTG::Astrodynamics::universe* Universe,
										 missionoptions* options);

		void calcbounds_phase_thruster_parameters(const string& prefix,
													 int first_X_entry_in_phase,
													 vector<double>* Xupperbounds,
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
													 int j,
													 int p,
													 EMTG::Astrodynamics::universe* Universe,
													 missionoptions* options);

		void calcbounds_LT_match_points(const string& prefix,
										int first_X_entry_in_phase,
										vector<double>* Xupperbounds,
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
										int j,
										int p,
										EMTG::Astrodynamics::universe* Universe,
										missionoptions* options);

		void find_dependencies_due_to_escape_spiral(vector<double>* Xupperbounds,
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
													int j,
													int p,
													missionoptions* options);

		void find_dependencies_due_to_capture_spiral(vector<double>* Xupperbounds,
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
													int j,
													int p,
													missionoptions* options);

		//GMAT output methods
		void output_GMAT_spacecraft(int& j, 
									int& p,
									vector<string>& SC_created,
									int& index_SC, 
									vector<EMTG::Astrodynamics::body>& missionbodies, 
									int& index_body_visited, 
									std::ofstream& GMATfile);

		void output_GMAT_fueltank_and_thruster(	int& j,
												int& p, 
												vector<EMTG::Astrodynamics::body>& missionbodies,
												int& index_body_visited, 
												std::ofstream& GMATfile);

		void output_GMAT_burn_objects(	int& j, 
										int& p,
										std::ofstream& GMATfile);

		void output_GMAT_create_state_and_time_variables(	int& j,
															int& p,
															std::ofstream& GMATfile);

		void output_GMAT_create_interphase_control_variables(	int& j,
																int& p, 
																missionoptions& options,
																std::ofstream& GMATfile);

		void output_GMAT_inter_phase_control_initial_guess(	int& j,
															int& p,
															missionoptions& options,
															std::ofstream& GMATfile);

		void output_GMAT_propagate_phase_forward_and_back(	int& j, 
															int& p, 
															vector<Astrodynamics::body>& missionbodies,
															int& index_body_visited,
															vector<string>& SC_created,
															int& index_SC,
															missionoptions& options,
															std::ofstream& GMATfile);

		//virtual method templates
		virtual int evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options) = 0;
		virtual int output(missionoptions* options, const double& launchdate, int j, int p,  EMTG::Astrodynamics::universe* Universe, int* eventcount) = 0;
		virtual int calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p,  EMTG::Astrodynamics::universe* Universe, missionoptions* options) = 0;
		virtual void output_GTOC7_format(missionoptions* options, EMTG::Astrodynamics::universe* Universe, const std::string& GTOC_output_file, int j, int p){};
		virtual void output_GTOC7_format_b(missionoptions* options, EMTG::Astrodynamics::universe* Universe, const std::string& GTOC_output_file, int j, int p){};
	
		//b-plane object
		Astrodynamics::bplane Bplane;

		//state information
		vector< vector<double> > spacecraft_state;
		vector< vector<double> > control;
		double state_at_beginning_of_phase[7]; //state vector of the spacecraft at the beginning of the phase, AFTER initial flyby/departure
		double state_at_end_of_phase[7]; //state vector of the spacecraft at the end of the phase, BEFORE any terminal maneuvers
		double state_at_initial_coast_midpoint[7], state_at_terminal_coast_midpoint[7];
		double dVdeparture[3]; //departure deltaV vector
		double dVarrival[3]; //arrival deltaV vector
		double current_mass_increment; //count up the mass increments preceding this phase
		double mission_initial_mass_multiplier;

		//fields containing information about the endpoints
		int boundary1_location_code; //integer code representing the first boundary point in the phase
		int boundary2_location_code; //integer code representing the second body in the phase
		Astrodynamics::body *Body1, *Body2;

		//time
		double phase_start_epoch;
		double phase_end_epoch;
		double phase_time_elapsed_forward;
		double phase_time_elapsed_backward;
		double TOF;
		double time_step_distribution_scale_or_stdv;
		double initial_coast_duration;
		double terminal_coast_duration;
		double total_available_thrust_time; //time after intial and terminal coasts have been removed
		vector<double> time_step_sizes;

		//calculation objects
		double pseudoa1, pseudoa2; //used for calculating bounds on phase flight time and, where applicable, the Sundman variable

		//flyby parameters
		double flyby_turn_angle;
		double flyby_altitude;
		double flyby_outgoing_v_infinity;
		double BdotR;
		double BdotT;
		double Btheta;
		double Bradius;
		double rp_check;

		//parameters for launches, departures, and arrivals
		double RA_departure;
		double RA_arrival;
		double DEC_departure;
		double DEC_arrival;
		double C3_departure;
		double C3_arrival;
		double theta_SOI_departure;
		double theta_SOI_arrival;
		double phi_SOI_departure;
		double phi_SOI_arrival;
		double unscaled_phase_initial_mass;
		double dmdvinf; //derivative of initial mass with respect to v-infinity
		double journey_initial_mass_increment_scale_factor;
		double dV_departure_magnitude; //dV magnitude for impulsive departure burns
		double dV_arrival_magnitude; //dV magnitude for impulsive arrival burns

		//b-plane calculation information
		math::Matrix<double> BoundaryR;
		math::Matrix<double> BoundaryV;
		math::Matrix<double> V_infinity_in;
		math::Matrix<double> V_infinity_out;

		//propulsion information
		double EMTG_chosen_power;
		double EMTG_chosen_Isp;
		vector<double> available_power;
		vector<double> available_thrust;
		vector<double> available_mass_flow_rate;
		vector<double> available_Isp;
		vector<double> active_power;
		vector<int> number_of_active_engines;

		//spiral things
		double spiral_escape_state_before_spiral[7];
		double spiral_escape_mass_before;
		double spiral_escape_mass_after;
		double spiral_escape_time;
		double spiral_escape_dv;
		double spiral_escape_Isp;
		double spiral_escape_thrust;
		double spiral_escape_mdot;
		int spiral_escape_number_of_engines;
		double spiral_escape_power;
		double spiral_escape_active_power;
		double spiral_capture_state_before_spiral[7];
		double spiral_capture_state_after_spiral[7];
		double spiral_capture_mass_before;
		double spiral_capture_mass_after;
		double spiral_capture_time;
		double spiral_capture_dv;
		double spiral_capture_Isp;
		double spiral_capture_thrust;
		double spiral_capture_mdot;
		int spiral_capture_number_of_engines;
		double spiral_capture_power;
		double spiral_capture_active_power;

		//************************derivative information

		//spiral derivatives
		double spiral_escape_dm_after_dm_before;
		double spiral_escape_dt_spiral_dm_before;
		double spiral_escape_dm_after_dIsp;
		double spiral_escape_dt_spiral_dIsp;
		double spiral_capture_dm_after_dm_before;
		double spiral_capture_dt_spiral_dm_before;
		double spiral_capture_dm_after_dIsp;
		double spiral_capture_dt_spiral_dIsp;

		//flyby constraints
		vector<int> flyby_constraints_X_indices; //current initial vx, vy, vz then previous terminal vz, vy, vx, 
		vector<double> flyby_constraints_X_scale_ranges;
		vector<int> flyby_velocity_magnitude_constraint_G_indices;
		vector<int> flyby_altitude_constraint_G_indices;

		//control vector constraints
		vector< vector<int> > control_vector_G_indices;

		//match point constraints
		vector< vector< vector <int> > > match_point_constraint_G_indices;
		vector<double> match_point_constraint_X_scale_ranges;
		vector< Kepler::STM > Forward_STM; //vector of state transition matrices
		vector< Kepler::STM > Backward_STM; //vector of state transition matrices
		Kepler::STM Current_STM;
		vector<double> Kepler_F_Forward, Kepler_Fdot_Forward, Kepler_G_Forward, Kepler_Gdot_Forward, Kepler_Fdotdot_Forward, Kepler_Gdotdot_Forward;
		vector<double> Kepler_F_Backward, Kepler_Fdot_Backward, Kepler_G_Backward, Kepler_Gdot_Backward, Kepler_Fdotdot_Backward, Kepler_Gdotdot_Backward;
		double Kepler_F_Current, Kepler_Fdot_Current, Kepler_G_Current, Kepler_Gdot_Current, Kepler_Fdotdot_Current, Kepler_Gdotdot_Current;
		vector<double> Propagation_Step_Time_Fraction_Forward, Propagation_Step_Time_Fraction_Backward;
		vector<double> Propagation_Step_Time_Fraction_Derivative_Forward, Propagation_Step_Time_Fraction_Derivative_Backward;
		vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass;
		vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass;
		vector< vector<int> > G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables;
		vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier;
		vector<int> G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier;
		vector<int> G_index_of_derivative_of_match_point_with_respect_to_BOL_power;
		double power_range;

		//derivatives of force model
		vector<double> dTdP;
		vector<double> dmdotdP;
		vector<double> dTdIsp;
		vector<double> dmdotdIsp;
		vector<double> dPdr;
		vector<double> dPdt;
		vector<double> dFSRPdr;
		vector< vector<double> > dagravdRvec;
		vector< vector<double> > dagravdtvec;

		//terminal velocity constraint
		vector<int> terminal_velocity_constraint_G_indices;
		vector<int> terminal_velocity_constraint_X_indices;
		vector<double> terminal_velocity_constraint_X_scale_ranges;

		//arrival declination constraint
		vector<int> arrival_declination_constraint_G_indices;
		vector<int> arrival_declination_constraint_X_indices;
		vector<double> arrival_declination_constraint_X_scale_ranges;

		//escape constraint
		vector<int> escape_constraint_G_indices;
		vector<int> escape_constraint_X_indices;
		vector<double> escape_constraint_X_scale_ranges;

		//left boundary central body exclusion radius constraint
		vector<int> left_boundary_central_body_exclusion_radius_constraint_G_indices;
		vector<int> left_boundary_central_body_exclusion_radius_constraint_X_indices;
		vector<double> left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges;


		//right boundary central body exclusion radius constraint
		vector<int> right_boundary_central_body_exclusion_radius_constraint_G_indices;
		vector<int> right_boundary_central_body_exclusion_radius_constraint_X_indices;
		vector<double> right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges;

		//time derivatives
		double left_boundary_state_derivative[6];
		double right_boundary_state_derivative[6];

		//variable boundary orbits
		double left_boundary_orbit_elements[6];
		vector<double> left_boundary_SMA_G_indices;
		vector<double> left_boundary_ECC_G_indices;
		vector<double> left_boundary_INC_G_indices;
		vector<double> left_boundary_RAAN_G_indices;
		vector<double> left_boundary_AOP_G_indices;
		vector<double> left_boundary_TA_G_indices;
		double right_boundary_orbit_elements[6];
		vector<double> right_boundary_SMA_G_indices;
		vector<double> right_boundary_ECC_G_indices;
		vector<double> right_boundary_INC_G_indices;
		vector<double> right_boundary_RAAN_G_indices;
		vector<double> right_boundary_AOP_G_indices;
		vector<double> right_boundary_TA_G_indices;
		double left_boundary_local_frame_state[6];
		vector<double> left_boundary_X_G_indices;
		vector<double> left_boundary_Y_G_indices;
		vector<double> left_boundary_VZ_G_indices;
		vector<double> left_boundary_XDOT_G_indices;
		vector<double> left_boundary_YDOT_G_indices;
		vector<double> left_boundary_ZDOT_G_indices;
		double right_boundary_local_frame_state[6];
		vector<double> right_boundary_X_G_indices;
		vector<double> right_boundary_Y_G_indices;
		vector<double> right_boundary_VZ_G_indices;
		vector<double> right_boundary_XDOT_G_indices;
		vector<double> right_boundary_YDOT_G_indices;
		vector<double> right_boundary_ZDOT_G_indices;
	};

} /* namespace EMTG */
#endif /* PHASE_H_ */
