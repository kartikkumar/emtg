
/*
 * MGAphase.cpp
 *
 *  Created on: Jul 15, 2012
 *      Author: Jacob
 */

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "MGAphase.h"
#include "missionoptions.h"
#include "Astrodynamics.h"

#ifdef _EMTG_proprietary
#include "Lambert.h"
#endif

#include "Lambert_AroraRussell.h"
#include "Kepler_Lagrange_Laguerre_Conway_Der.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "universe.h"

#include "SpiceUsr.h"



using namespace std;

namespace EMTG 
{

    MGA_phase::MGA_phase()
    {
	    //default constructor doesn't do anything
    }

    MGA_phase::MGA_phase(const int& j, const int& p, const missionoptions& options)
    {
        //call phase initialize method
        this->initialize(j, p, options);

	    //if this is a terminal phase, there are two burns
	    //otherwise there is only one
	    if (p < (options.number_of_phases[j] - 1))
		    dVmag.resize(1,0);
	    else
		    dVmag.resize(2,0);
    }

    MGA_phase::~MGA_phase()
    {
	    //destructor does not need to do anything special
    }

    //evaluate function
    //return 0 if successful, 1 if failure
    int MGA_phase::evaluate(const double* X, 
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
                            const int& j,
                            const int& p, 
                            EMTG::Astrodynamics::universe* Universe, 
                            missionoptions* options)
    {
	    //declare some local variables
	    int errcode = 0;
	    double lambert_v1[3];
	    double lambert_v2[3];
	    double hz; //z component of angular momentum

	    //set pointers to the bodies
	    if (boundary1_location_code > 0)
		    Body1 = &Universe->bodies[boundary1_location_code - 1];
	    if (boundary2_location_code > 0)
		    Body2 = &Universe->bodies[boundary2_location_code - 1];

	    //we need to know if we are the first phase of the journey
	    if (p == 0)
	    {
		    //Step 1: extract the journey start epoch
		    //if this is the first journey, the first decision variable is the starting epoch
		    //otherwise, the first decision variable is the journey wait time
		    if (j == 0)
		    {
			    *current_epoch = X[*Xindex];
			    this->phase_wait_time = X[*Xindex] - options->launch_window_open_date;
			    ++(*Xindex);
		    }
		    else
		    {
			    *current_epoch += X[*Xindex];
			    this->phase_wait_time = X[*Xindex];
			    ++(*Xindex);
		    }

		    //Step 2: locate the first body
		    this->locate_boundary_point(boundary1_location_code,
									    options->journey_departure_type[j],
									    true,
									    Universe, 
									    boundary1_state,
									    current_state+3,
									    *current_epoch, 
									    X, 
									    Xindex,
									    F,
									    Findex,
									    G, 
									    Gindex,
									    false,
									    j,
									    p,
									    options);
	    }
	    //there is no alternate step 2
	
	    //Step 3: find the time of flight
	    TOF = X[*Xindex];
	    phase_start_epoch = *current_epoch;
	    phase_end_epoch = *current_epoch + TOF;
	    ++(*Xindex);

	    //Step 4: locate the second body
	    this->locate_boundary_point(boundary2_location_code, 
								    options->journey_arrival_type[j],
								    false,
								    Universe,
								    boundary2_state,
								    current_state+3, 
								    *current_epoch + TOF,
								    X,
								    Xindex, 
								    F, 
								    Findex,
								    G,
								    Gindex, 
								    false, 
								    j,
								    p, 
								    options);

	    //Step 5: solve Lambert's problem between the planets
	    int Nrev;
	    bool ShortPeriod;
	    if (options->maximum_number_of_lambert_revolutions)
	    {
		    int LambertArcType = (int)X[*Xindex];
		    ++(*Xindex);
		    Nrev = LambertArcType / 2;
		    ShortPeriod = LambertArcType % 2;
	    }
	    else
	    {
		    Nrev = 0;
		    ShortPeriod = 0;
	    }

	    double Lam_error;
	    int Lam_iterations;

    #ifdef _EMTG_proprietary
	    if (options->LambertSolver == 1) //Izzo Lambert solver
	    {
		    EMTG::Astrodynamics::Lambert::Lambert_Izzo(boundary1_state,
			                                            boundary2_state,
			                                            this->TOF,
			                                            Universe->mu,
			                                            1,
			                                            Nrev,
			                                            lambert_v1,
			                                            lambert_v2);

		    //check the angular momentum vector to avoid going retrograde
		    hz = boundary1_state[0] * lambert_v1[1] - boundary1_state[1] * lambert_v1[0];

		    if (hz < 0)
		    {
			    EMTG::Astrodynamics::Lambert::Lambert_Izzo(boundary1_state,
				                                            boundary2_state,
				                                            this->TOF,
				                                            Universe->mu,
				                                            0,
				                                            Nrev,
				                                            lambert_v1,
				                                            lambert_v2);
		    }
	    }
	    else
	    {
    #endif
		    //note: if the Lambert solver fails it is almost always because there is no solution for that number of revolutions
		    //therefore wrap the solver in a try-catch block and if it fails, try the zero-rev case
		    try
		    {
			    EMTG::Astrodynamics::Lambert::Lambert_AroraRussell(boundary1_state,
				                                                    boundary2_state,
				                                                    TOF,
				                                                    Universe->mu,
				                                                    Nrev,
				                                                    1,
				                                                    ShortPeriod,
				                                                    1e-13,
				                                                    30,
				                                                    lambert_v1,
				                                                    lambert_v2,
				                                                    Lam_error,
				                                                    Lam_iterations);
		    }
		    catch (int& errorcode)
		    {
			    if (errorcode == 200000) //multi-rev error
				    EMTG::Astrodynamics::Lambert::Lambert_AroraRussell(boundary1_state,
				                                                        boundary2_state,
				                                                        TOF,
				                                                        Universe->mu,
				                                                        0,
				                                                        1,
				                                                        ShortPeriod,
				                                                        1e-13,
				                                                        30,
				                                                        lambert_v1,
				                                                        lambert_v2,
				                                                        Lam_error,
				                                                        Lam_iterations);
		    }


		    //check the angular momentum vector
		    hz = boundary1_state[0] * lambert_v1[1] - boundary1_state[1] * lambert_v1[0];

		    if (hz < 0)
		    {
			    try
			    {
				    EMTG::Astrodynamics::Lambert::Lambert_AroraRussell(boundary1_state,
					                                                    boundary2_state,
					                                                    TOF,
					                                                    Universe->mu,
					                                                    Nrev,
					                                                    1,
					                                                    ShortPeriod,
					                                                    1e-13,
					                                                    30,
					                                                    lambert_v1,
					                                                    lambert_v2,
					                                                    Lam_error,
					                                                    Lam_iterations);
			    }
			    catch (int& errorcode)
			    {
				    if (errorcode == 200000)
					    EMTG::Astrodynamics::Lambert::Lambert_AroraRussell(boundary1_state,
					                                                        boundary2_state,
					                                                        TOF,
					                                                        Universe->mu,
					                                                        0,
					                                                        1,
					                                                        ShortPeriod,
					                                                        1e-13,
					                                                        30,
					                                                        lambert_v1,
					                                                        lambert_v2,
					                                                        Lam_error,
					                                                        Lam_iterations);
			    }
		    }
    #ifdef _EMTG_proprietary
	    }
    #endif

	    //Step 6: compute all parameters of the first event of the phase
	    //if this is the first phase in the journey, compute RA and DEC. Otherwise process the flyby at the beginning of the phase
	    //either way, compute outgoing C3
	    if (p == 0)
	    {
		    for (int k = 0; k < 3; ++k)
			    dVdeparture[k] = lambert_v1[k] - boundary1_state[k+3];

		    //compute C3
		    C3_departure = EMTG::math::dot(dVdeparture, dVdeparture, 3);

		    //compute declination
		    DEC_departure = asin(dVdeparture[2] / sqrt(C3_departure));

		    //there is a DLA constraint for the first journey
		    if (j == 0)
		    {
			    F[*Findex] = DEC_departure * 180.0 / EMTG::math::PI;
			    ++(*Findex);
		    }

		    //compute right ascension
		    RA_departure = atan2(dVdeparture[1], dVdeparture[0]);

		    //compute the burn necessary to enter into the departure orbit
		    if (options->journey_departure_type[j] == 0 || options->journey_departure_type[j] == 2) //this is a direct insertion or launch
		    {
			    dVmag[0] = sqrt(C3_departure);
			
			    if (j == 0) //this is a launch
			    {

				    //figure out how much mass we can launch
				    if (options->LV_type > 0 || options->LV_type == -2)
				    {
					    double launch_mass;
					    EMTG::Astrodynamics::find_mass_to_orbit(C3_departure, DEC_departure * 180.0 / math::PI, options->LV_type, &launch_mass, &dmdvinf, options);
					    state_at_beginning_of_phase[6] = launch_mass > options->maximum_mass ? options->maximum_mass : launch_mass;
					
					    //apply the mass margin
					    state_at_beginning_of_phase[6] *= 1.0 - options->LV_margin;
				    }
				    else if (options->LV_type == 0)
				    {
					    state_at_beginning_of_phase[6] = options->maximum_mass;

					    //add the starting mass increment
					    state_at_beginning_of_phase[6] += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					    dmdvinf = 0.0;
				    }
				    else //chemical burn using departure stage engine
				    {
					    double expfun = exp(-dVmag[0] * 1000 / (options->IspDS * options->g0));

					    double initialmass = options->maximum_mass;

					    //add the starting mass increment
					    initialmass += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					    state_at_beginning_of_phase[6] = initialmass * expfun;
					    dmdvinf = -initialmass * 1000 / (options->IspDS * options->g0) * expfun;
				    }
			    }
			    else if (options->journey_departure_type[j] == 0)
			    {
				    double expfun = exp(-dVmag[0] * 1000 / (options->IspChem * options->g0));
				    if (j > 0)
				    {
					    double initialmass = current_state[6];

					    //add the starting mass increment
					    initialmass += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					    state_at_beginning_of_phase[6] = initialmass * expfun;
					    dmdvinf = -initialmass * 1000 / (options->IspChem * options->g0) * expfun;
				    }
				    else
				    {
					    double initialmass = options->maximum_mass;

					    //add the starting mass increment
					    initialmass += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

					    state_at_beginning_of_phase[6] = initialmass * expfun;
					    dmdvinf = -options->maximum_mass * 1000 / (options->IspChem * options->g0) * expfun;
				    }
			    }
			    else //"free" departure
			    {
				    state_at_beginning_of_phase[6] = (j == 0 ? options->maximum_mass : current_state[6]);

				    //add the starting mass increment
				    state_at_beginning_of_phase[6] += journey_initial_mass_increment_scale_factor * options->journey_starting_mass_increment[j];

				    dmdvinf = 0.0;
			    }
		    }
		    else if (options->journey_departure_type[j] == 1) //depart from a parking orbit
		    {
			    double vinf;
			    dVmag[0] = EMTG::Astrodynamics::insertion_burn(state_at_beginning_of_phase+3, boundary1_state+3, Universe->bodies[boundary1_location_code-1].mu, Universe->bodies[boundary1_location_code-1].r_SOI, options->journey_departure_elements[j][0], options->journey_departure_elements[j][1], &vinf);
			    *current_deltaV += dVmag[0];
			    C3_departure = vinf*vinf;

			    if (boundary1_location_code < 0) //cannot do a simplified parking orbit departure unless the departure point is a body!
			    {
				    cout << "Cannot do a simplified parking orbit departure unless the departure point is a body. For other departure point types, use a 'direct' departure" << endl;
				    throw 20;
			    }

			    if (j == 0) //if this is the first journey, no mass is assigned yet
				    state_at_beginning_of_phase[6] = options->maximum_mass * exp(-dVmag[0] * 1000/ (options->IspDS * options->g0));
			    else
				    state_at_beginning_of_phase[6] = current_state[6] * exp(-dVmag[0] * 1000/ (options->IspDS * options->g0));
		    }
	    }
	    else
	    {
		    //phase starts with a flyby
		    double flyby_orbit_energy;
		    EMTG::Astrodynamics::powered_flyby (current_state+3, lambert_v1, boundary1_state+3, Body1->mu, Body1->radius, Body1->r_SOI,
				    &dVmag[0], &flyby_altitude, &flyby_turn_angle, &C3_departure, &flyby_orbit_energy);

		    //calculate the b-plane parameters, check the periapse altitude
		    for (int k = 0; k < 3; ++k)
		    {
			    V_infinity_in(k) = current_state[k+3] - boundary1_state[k+3];
			    V_infinity_out(k) = lambert_v1[k] - boundary1_state[k+3];
		    }
		    BoundaryR.assign_all(boundary1_state);
		    BoundaryV.assign_all(boundary1_state+3);

		    //apply flyby altitude constraint
		    F[*Findex] = (Body1->minimum_safe_flyby_altitude - flyby_altitude) / (Body1->minimum_safe_flyby_altitude + Body1->radius);
		    ++(*Findex);

		    //apply flyby orbit energy constraint
		    F[*Findex] = flyby_orbit_energy;
		    ++(*Findex);

		    //compute the outgoing mass
		    state_at_beginning_of_phase[6] = current_state[6] * exp(-dVmag[0] * 1000 / (options->IspChem * options->g0));
	    }
	    //update the current total deltaV
	    *current_deltaV += dVmag[0];

	    //compute the outgoing state vector from the first body
	    for (int k = 0; k < 3; ++k)
	    {
		    state_at_beginning_of_phase[k] = boundary1_state[k];
		    state_at_beginning_of_phase[k+3] = lambert_v1[k];
	    }
	
	    //Step 7: Process the end of the phase and compute the final state vector. Is this an arrival?
	    if (p == options->number_of_phases[j] - 1)
	    {
		    if (boundary2_location_code > 0) //ending at body
			    dVmag[1] = this->process_arrival(	lambert_v2, 
												    boundary2_state,
												    current_state+3, 
												    current_epoch,
												    Body2->mu, 
												    Body2->r_SOI,
												    F,
												    Findex,
												    j,
												    options, 
												    Universe);
		    else //ending at point on central body SOI, fixed point, or fixed orbit
			    dVmag[1] = this->process_arrival(	lambert_v2, 
												    boundary2_state,
												    current_state+3, 
												    current_epoch,
												    Universe->mu,
												    Universe->r_SOI, 
												    F,
												    Findex,
												    j,
												    options,
												    Universe);
		
		    state_at_end_of_phase[6] = state_at_beginning_of_phase[6] * exp(-dVmag[1] * 1000/ (options->IspChem * options->g0));
		    *current_deltaV += dVmag[1];
	    }
	    else
	    {
		    state_at_end_of_phase[6] = state_at_beginning_of_phase[6];
	    }
	    for (int k = 0; k < 3; ++k)
	    {
		    state_at_end_of_phase[k] = boundary2_state[k];
		    state_at_end_of_phase[k+3] = lambert_v2[k];
	    }

	    //Step 8: advance the current epoch
	    *current_epoch += TOF;
	
	    //Step 9: update the current state
	    for (int k = 0; k < 7; ++k)
		    current_state[k] = state_at_end_of_phase[k];

	    return errcode;
    }

    //output function
    //return 0 if successful, 1 if failure
    void MGA_phase::output(missionoptions* options,
                            const double& launchdate,
                            const int& j,
                            const int& p,
                            EMTG::Astrodynamics::universe* Universe,
                            int* eventcount) 
    {
	    //Step 1: store data that will be used for the printing
	    double empty_vector[] = {0,0,0};
	    string event_type;

	    if (p > 0 || (options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4))
	    {
		    event_type = "upwr_flyby";
		    math::Matrix<double> periapse_state = calculate_flyby_periapse_state(V_infinity_in, V_infinity_out, flyby_altitude, *Body1);
		    math::Matrix<double> periapse_R(3, 1);
		    periapse_R(0) = periapse_state(0);
		    periapse_R(1) = periapse_state(1);
		    periapse_R(2) = periapse_state(2);
		    Bplane.define_bplane(V_infinity_in, BoundaryR, BoundaryV);
		    Bplane.compute_BdotR_BdotT_from_periapse_position(Body1->mu, V_infinity_in, periapse_R, &BdotR, &BdotT);

		    //compute RA and DEC in the frame of the target body
		    this->Body1->J2000_body_equatorial_frame.construct_rotation_matrices(this->phase_start_epoch / 86400.0 + 2400000.5);
		    math::Matrix<double> rot_out_vec = this->Body1->J2000_body_equatorial_frame.R_from_ICRF_to_local * V_infinity_in;

		    this->RA_departure = atan2(rot_out_vec(1), rot_out_vec(0));

		    this->DEC_departure = asin(rot_out_vec(2) / V_infinity_in.norm());
	    }
	    else if (j == 0 && boundary1_location_code > 0 && (options->LV_type >= 0 || options->LV_type == -2))
		    event_type = "launch";
	    else if (options->journey_departure_type[j] == 6)
	    {
		    event_type = "zeroflyby";
		    //compute RA and DEC in the frame of the target body
		    this->Body1->J2000_body_equatorial_frame.construct_rotation_matrices(this->phase_start_epoch / 86400.0 + 2400000.5);
		    math::Matrix<double> rot_out_vec = this->Body1->J2000_body_equatorial_frame.R_from_ICRF_to_local * V_infinity_in;

		    this->RA_departure = atan2(rot_out_vec(1), rot_out_vec(0));

		    this->DEC_departure = asin(rot_out_vec(2) / V_infinity_in.norm());
	    }
	    else
		    event_type = "departure";

	    string boundary1_name;
	    string boundary2_name;

		    switch (boundary1_location_code)
	    {
		    case -1:
			    {
				    boundary2_name = "free point";
				    break;
			    }
		    case -2: //begin at SOI
			    {
				    boundary1_name = "Hyp-arrival";
				    break;
			    }
		    default:
			    boundary1_name = (Universe->bodies[boundary1_location_code - 1].name);
	    }

	    switch (boundary2_location_code)
	    {
		    case -1:
			    {
				    boundary2_name = "free point";
				    break;
			    }
		    default:
			    boundary2_name = (Universe->bodies[boundary2_location_code - 1].name);
	    }
	    double initial_Isp;
	    if (j == 0 && p == 0)
	    {
		    if (options->journey_departure_type[j] == 1 || options->LV_type == -1)
		    {
			    initial_Isp =  options->IspDS;
		    }
		    else
		    {
			    initial_Isp =  -1;
		    }
	    }
	    else
	    {
		    initial_Isp = options->IspChem;
	    }
	
	    //*****************************************************************************
	    //first let's print the departure/flyby
	
	    write_summary_line(options,
						    Universe,
						    eventcount,
						    phase_start_epoch / 86400.0,
						    event_type,
						    boundary1_name,
						    0,
                            ( (p > 0 || options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4) ? flyby_altitude : Bradius),
						    (Btheta),
                            ( (p > 0 || options->journey_departure_type[j] == 3 || options->journey_departure_type[j] == 4) ? flyby_turn_angle : -1),
						    RA_departure,
						    DEC_departure,
						    C3_departure,
						    state_at_beginning_of_phase,
						    dVdeparture,
						    empty_vector,
						    dVmag[0],
						    -1,
						    initial_Isp,
						    -1,
						    0,
						    0,
						    0);

	    //*****************************************************************************
	    //coast until the end of the phase
	    double timestep = TOF / options->num_timesteps;
	    double output_state[7];
	    double nodV[] = {0,0,0};
	    output_state[6] = state_at_beginning_of_phase[6];
	
	    for (int step = 0; step < options->num_timesteps; ++step)
	    {
		    //compute the current epoch
		    double epoch = phase_start_epoch + timestep * (step + 0.5);

		    //propagate the spacecraft
		    Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_beginning_of_phase,
													    output_state,
													    Universe->mu,
													    Universe->LU,
													    (epoch - phase_start_epoch));

		    //write the summary line
		    write_summary_line(options,
						    Universe,
						    eventcount,
						    epoch / 86400.0,
						    "coast",
						    "deep-space",
						    timestep / 86400.0,
						    -1,
						    -1,
						    -1,
						    0,
						    0,
						    0,
						    output_state,
						    nodV,
						    empty_vector,
						    0,
						    -1,
						    -1,
						    -1,
						    0,
						    0,
						    0);
	    }

	    //*****************************************************************************
	    //finally, terminal phases have an arrival maneuver
	
	

	    if (p == options->number_of_phases[j] - 1)
	    {
		    if (options->journey_arrival_type[j] == 0)
			    event_type = "insertion";
		    else if (options->journey_arrival_type[j] == 1 || options->journey_arrival_type[j] == 3)
			    event_type = "rendezvous";
		    else if (options->journey_arrival_type[j] == 2)
			    event_type = "intercept";
		    else if (options->journey_arrival_type[j] == 4)
			    event_type = "match-vinf";

		    //compute RA and DEC in the frame of the target body
		    if (options->destination_list[j][1] > 0)
		    {
			    this->Body2->J2000_body_equatorial_frame.construct_rotation_matrices((this->phase_start_epoch + this->TOF) / 86400.0 + 2400000.5);
			    math::Matrix<double> rot_in_vec(3, 1, this->dVarrival);
			    math::Matrix<double> rot_out_vec = this->Body2->J2000_body_equatorial_frame.R_from_ICRF_to_local * rot_in_vec;

			    this->RA_arrival = atan2(rot_out_vec(1), rot_out_vec(0));

			    this->DEC_arrival = asin(rot_out_vec(2) / sqrt(this->C3_arrival));
		    }
		    else
		    {
			    this->RA_arrival = 0.0;
			    this->DEC_arrival = 0.0;
		    }
	
		    double dV_arrival_mag;
		    if (options->journey_arrival_type[j] == 2)
		    {
			    dV_arrival_mag = sqrt(C3_arrival);
		    }
		    else if (options->journey_arrival_type[j] == 4)
		    {
			    dV_arrival_mag = 0;
			    dVarrival[0] = 0;
			    dVarrival[1] = 0;
			    dVarrival[2] = 0;
		    }
		    else
		    {
			    dV_arrival_mag = dVmag[1];
		    }

		    write_summary_line(options,
						    Universe,
						    eventcount,
						    (phase_start_epoch + TOF) / 86400.0,
						    event_type,
						    boundary2_name,
						    0,
						    -1,
						    -1,
						    -1,
						    RA_arrival,
						    DEC_arrival,
						    C3_arrival,
						    state_at_end_of_phase,
						    dVarrival,
						    empty_vector,
						    dV_arrival_mag,
						    -1,
						    options->IspChem,
						    -1,
						    0,
						    0,
						    0);
	    }
    }

    //bounds calculation function
    //return 0 if successful, 1 if failure
    void MGA_phase::calcbounds(vector<double>* Xupperbounds,
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
                                const int& j,
                                const int& p, 
                                EMTG::Astrodynamics::universe* Universe,
                                missionoptions* options)
    {
	    //this function calculates the upper and lower bounds for the decision and constraint vectors for MGA

	    //create a prefix string with journey and phase information
	    stringstream prefixstream;
	    prefixstream << "j" << j << "p" << p << ": ";
	    string prefix = prefixstream.str();
	    int first_X_entry_in_phase = Xupperbounds->size();

	    //first, we need to know if we are the first phase in the journey
	    if (p == 0)
	    {
		    //if we are the first phase, we also need to know if we are the first journey
		    if (j == 0)
		    {
			    //if so, then the first decision variable is the launch epoch in MJD
			    Xlowerbounds->push_back(options->launch_window_open_date + options->journey_wait_time_bounds[j][0]);
			    Xupperbounds->push_back(options->launch_window_open_date + options->journey_wait_time_bounds[j][1]);
			    Xdescriptions->push_back(prefix + "launch epoch (MJD)");

			    //and we have a DLA constraint
			    if (j == 0 && boundary1_location_code > 0) //if this is the first journey and we are leaving from a planet, i.e. if this is a launch
			    {
				    Flowerbounds->push_back(options->DLA_bounds[0] * math::PI / 180.0);
				    Fupperbounds->push_back(options->DLA_bounds[1] * math::PI / 180.0);
				    Fdescriptions->push_back(prefix + "DLA constraint (degrees)");
			    }
			    else
			    {
				    Flowerbounds->push_back(-math::PI / 2.0);
				    Fupperbounds->push_back(math::PI / 2.0);
				    Fdescriptions->push_back(prefix + "DLA constraint (degrees)");
			    }
		    }
		    else
		    {
			    //if we are not the first journey, are we starting from a hyperbolic arrival?
			    //if so, then there is no wait time. If not, then the first decision variable is the stay time at the first body in the journey (i.e. at the asteroid for sample return)
			    if (!(options->sequence[j-1][p+1] == -2 || boundary1_location_code == -2))
			    {
				    Xlowerbounds->push_back(options->journey_wait_time_bounds[j][0]);
				    Xupperbounds->push_back(options->journey_wait_time_bounds[j][1]);
				    Xdescriptions->push_back(prefix + "stay time (days)");
			    }
		    }

		    if  (boundary1_location_code == -1) //if this boundary point is at a free point in space, with the various elements either fixed or free
		    {
			    vector<string> CartesianElementNames;
			    CartesianElementNames.push_back("x (km)");
			    CartesianElementNames.push_back("y (km)");
			    CartesianElementNames.push_back("z (km)");
			    CartesianElementNames.push_back("xdot (km/s)");
			    CartesianElementNames.push_back("ydot (km/s)");
			    CartesianElementNames.push_back("zdot (km/s)");

			    vector<string> ClassicalOrbitElementNames;
			    ClassicalOrbitElementNames.push_back("SMA (km)");
			    ClassicalOrbitElementNames.push_back("ECC (km)");
			    ClassicalOrbitElementNames.push_back("INC (rad)");
			    ClassicalOrbitElementNames.push_back("RAAN (rad)");
			    ClassicalOrbitElementNames.push_back("AOP (rad)");
			    ClassicalOrbitElementNames.push_back("TA (rad)");

			    for (int k = 0; k < 6; ++k)
			    {
				    if (options->journey_arrival_elements_vary_flag[j][k])
				    {
					    Xlowerbounds->push_back(options->journey_arrival_elements_bounds[j][k][0]);
					    Xupperbounds->push_back(options->journey_arrival_elements_bounds[j][k][1]);
					    if (options->journey_arrival_elements_type[j])
						    Xdescriptions->push_back(prefix + " left boundary point " + ClassicalOrbitElementNames[k]);
					    else
						    Xdescriptions->push_back(prefix + " left boundary point " + CartesianElementNames[k]);
				    }
			    }

			    //if it is possible for the optimizer to select a point inside the exclusion zone of the central body
				    //then there must be a nonlinear constraint to prevent this
				    if (options->journey_departure_elements_type[j]) //classical orbit elements
				    {
					    //this constraint is applied if we are varying SMA or ECC
					    if (options->journey_departure_elements_vary_flag[j][0] || options->journey_departure_elements_vary_flag[j][1])
					    {
						    Flowerbounds->push_back(0.0);
						    Fupperbounds->push_back(math::LARGE);
						    Fdescriptions->push_back(prefix + " left boundary central body exclusion radius constraint");

						    //this constraint has derivatives with respect to SMA, ECC, and TA
						    //only create a derivative entry with respect to an orbit element if that element is being varied
						    if (options->journey_departure_elements_vary_flag[j][0]) //SMA
						    {
							    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							    {
								    if ((*Xdescriptions)[entry].find("left boundary point SMA") < 1024)
								    {
									    iGfun->push_back(Fdescriptions->size() - 1);
									    jGvar->push_back(entry);
									    stringstream EntryNameStream;
									    EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									    Gdescriptions->push_back(EntryNameStream.str());
									    left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									    left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									    left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								    }
							    }
						    }
						    if (options->journey_departure_elements_vary_flag[j][1]) //ECC
						    {
							    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							    {
								    if ((*Xdescriptions)[entry].find("left boundary point ECC") < 1024)
								    {
									    iGfun->push_back(Fdescriptions->size() - 1);
									    jGvar->push_back(entry);
									    stringstream EntryNameStream;
									    EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									    Gdescriptions->push_back(EntryNameStream.str());
									    left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									    left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									    left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								    }
							    }
						    }
						    if (options->journey_departure_elements_vary_flag[j][5]) //TA
						    {
							    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
							    {
								    if ((*Xdescriptions)[entry].find("left boundary point TA") < 1024)
								    {
									    iGfun->push_back(Fdescriptions->size() - 1);
									    jGvar->push_back(entry);
									    stringstream EntryNameStream;
									    EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
									    Gdescriptions->push_back(EntryNameStream.str());
									    left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
									    left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
									    left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
								    }
							    }
						    }
					    }
				    }
				    else //cartesian orbit elements
				    {
					    //this constraint is applied if we are varying x, y, or z
					    if (options->journey_departure_elements_vary_flag[j][0] || options->journey_departure_elements_vary_flag[j][1] || options->journey_departure_elements_vary_flag[j][2])
					    {
						    Flowerbounds->push_back(0.0);
						    Fupperbounds->push_back(math::LARGE);
						    Fdescriptions->push_back(prefix + " left boundary central body exclusion radius constraint");
					    }

					    //this constraint has derivatives with respect to x, y, and z
					    //only create a derivative entry with respect to an orbit element if that element is being varied
					    if (options->journey_departure_elements_vary_flag[j][0]) //x
					    {
						    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						    {
							    if ((*Xdescriptions)[entry].find("left boundary point x (km)") < 1024)
							    {
								    iGfun->push_back(Fdescriptions->size() - 1);
								    jGvar->push_back(entry);
								    stringstream EntryNameStream;
								    EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								    Gdescriptions->push_back(EntryNameStream.str());
								    left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								    left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								    left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							    }
						    }
					    }
					    if (options->journey_departure_elements_vary_flag[j][1]) //y
					    {
						    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						    {
							    if ((*Xdescriptions)[entry].find("left boundary point y (km)") < 1024)
							    {
								    iGfun->push_back(Fdescriptions->size() - 1);
								    jGvar->push_back(entry);
								    stringstream EntryNameStream;
								    EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								    Gdescriptions->push_back(EntryNameStream.str());
								    left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								    left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								    left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							    }
						    }
					    }
					    if (options->journey_departure_elements_vary_flag[j][2]) //z
					    {
						    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						    {
							    if ((*Xdescriptions)[entry].find("left boundary point z (km)") < 1024)
							    {
								    iGfun->push_back(Fdescriptions->size() - 1);
								    jGvar->push_back(entry);
								    stringstream EntryNameStream;
								    EntryNameStream << "Derivative of " << prefix << " left boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								    Gdescriptions->push_back(EntryNameStream.str());
								    left_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								    left_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								    left_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							    }
						    }
					    }
				    }
		    }
		    //if we are starting at periapse of an arrival hyperbola, we must choose Bradius and Btheta
		    else if (boundary1_location_code == -2)
		    {
			    Xlowerbounds->push_back(Universe->minimum_safe_distance);
			    Xupperbounds->push_back(10 * Universe->minimum_safe_distance);
			    Xdescriptions->push_back(prefix + "left boundary B-plane 'Bradius'");

			    Xlowerbounds->push_back(-2 * EMTG::math::PI);
			    Xupperbounds->push_back(2 * EMTG::math::PI);
			    Xdescriptions->push_back(prefix + "left boundary B-plane 'Btheta'");

			    //we must introduce a constraint that the periapse point must not be too close to the planet
			    Flowerbounds->push_back(-EMTG::math::LARGE);
			    Fupperbounds->push_back(0.0);
			    Fdescriptions->push_back(prefix + "left boundary periapse point must be a safe distance from the planet");
		    }
	    }
	    else
	    {
		    //if we are NOT the first phase then we don't need to worry about additional decision variables
		    //BUT we do need to implement three constraints: the flyby altitude constraint and the flyby hyperbolic constraint


		    //we also need to encode a no-collision constraint
		    if (Universe->bodies[boundary1_location_code-1].mass < 1.0e+25)
			    Flowerbounds->push_back(-10.0);
		    else
			    Flowerbounds->push_back(-300.0);
		    Fupperbounds->push_back(0.0);
		    Fdescriptions->push_back(prefix + "flyby altitude constraint (above minimum altitude but below [10x, 300x] altitude for [rocky, gas] planets");

		    //flyby hyperbolic constraint
		    Flowerbounds->push_back(0.0);
		    Fupperbounds->push_back(math::LARGE);
		    Fdescriptions->push_back(prefix + "flyby orbit energy plus 10% fudge factor must be positive");
	    }

	    //**************************************************************************
	    //next, we need to encode the phase flight time
	    calcbounds_flight_time(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, synodic_periods, j, p, Universe, options);

	    //**************************************************************************
	    //if the allowed number of Lambert revolutions is greater than zero then we must encode the Lambert type variable
	    if (options->maximum_number_of_lambert_revolutions)
	    {
		    Xlowerbounds->push_back(0.0);
		    Xupperbounds->push_back(ceil((double)2 * options->maximum_number_of_lambert_revolutions + 1));
		    Xdescriptions->push_back(prefix + "Lambert arc type");
	    }

	    //******************
	    //if we are the last journey, then encode any variables necessary for the right hand boundary condition
	    if (p == (options->number_of_phases[j] - 1))
	    {
		    //if this boundary point is at a free point in space, with the various elements either fixed or free
		    //this is only relevant for the first journey - succcessive journeys will start from the right hand boundary of the previous journey
		    if (boundary1_location_code == -1 && j == 0)
		    {
			    vector<string> CartesianElementNames;
			    CartesianElementNames.push_back("x (km)");
			    CartesianElementNames.push_back("y (km)");
			    CartesianElementNames.push_back("z (km)");
			    CartesianElementNames.push_back("xdot (km/s)");
			    CartesianElementNames.push_back("ydot (km/s)");
			    CartesianElementNames.push_back("zdot (km/s)");

			    vector<string> ClassicalOrbitElementNames;
			    ClassicalOrbitElementNames.push_back("SMA (km)");
			    ClassicalOrbitElementNames.push_back("ECC (km)");
			    ClassicalOrbitElementNames.push_back("INC (rad)");
			    ClassicalOrbitElementNames.push_back("RAAN (rad)");
			    ClassicalOrbitElementNames.push_back("AOP (rad)");
			    ClassicalOrbitElementNames.push_back("TA (rad)");

			    for (int k = 0; k < 6; ++k)
			    {
				    if (options->journey_arrival_elements_vary_flag[j][k])
				    {
					    Xlowerbounds->push_back(options->journey_arrival_elements_bounds[j][k][0]);
					    Xupperbounds->push_back(options->journey_arrival_elements_bounds[j][k][1]);
					    if (options->journey_arrival_elements_type[j])
						    Xdescriptions->push_back(prefix + " right boundary point " + ClassicalOrbitElementNames[k]);
					    else
						    Xdescriptions->push_back(prefix + " right boundary point " + CartesianElementNames[k]);
				    }
			    }


			    //if it is possible for the optimizer to select a point inside the exclusion zone of the central body
			    //then there must be a nonlinear constraint to prevent this
			    if (options->journey_arrival_elements_type[j]) //classical orbit elements
			    {
				    //this constraint is applied if we are varying SMA or ECC
				    if (options->journey_arrival_elements_vary_flag[j][0] || options->journey_arrival_elements_vary_flag[j][1])
				    {
					    Flowerbounds->push_back(0.0);
					    Fupperbounds->push_back(math::LARGE);
					    Fdescriptions->push_back(prefix + " right boundary central body exclusion radius constraint");

					    //this constraint has derivatives with respect to SMA, ECC, and TA
					    //only create a derivative entry with respect to an orbit element if that element is being varied
					    if (options->journey_arrival_elements_vary_flag[j][0]) //SMA
					    {
						    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						    {
							    if ((*Xdescriptions)[entry].find("right boundary point SMA") < 1024)
							    {
								    iGfun->push_back(Fdescriptions->size() - 1);
								    jGvar->push_back(entry);
								    stringstream EntryNameStream;
								    EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								    Gdescriptions->push_back(EntryNameStream.str());
								    right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								    right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								    right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							    }
						    }
					    }
					    if (options->journey_arrival_elements_vary_flag[j][1]) //ECC
					    {
						    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						    {
							    if ((*Xdescriptions)[entry].find("right boundary point ECC") < 1024)
							    {
								    iGfun->push_back(Fdescriptions->size() - 1);
								    jGvar->push_back(entry);
								    stringstream EntryNameStream;
								    EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								    Gdescriptions->push_back(EntryNameStream.str());
								    right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								    right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								    right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							    }
						    }
					    }
					    if (options->journey_arrival_elements_vary_flag[j][5]) //TA
					    {
						    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
						    {
							    if ((*Xdescriptions)[entry].find("right boundary point TA") < 1024)
							    {
								    iGfun->push_back(Fdescriptions->size() - 1);
								    jGvar->push_back(entry);
								    stringstream EntryNameStream;
								    EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
								    Gdescriptions->push_back(EntryNameStream.str());
								    right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
								    right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
								    right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
							    }
						    }
					    }
				    }
			    }
			    else //cartesian orbit elements
			    {
				    //this constraint is applied if we are varying x, y, or z
				    if (options->journey_arrival_elements_vary_flag[j][0] || options->journey_arrival_elements_vary_flag[j][1] || options->journey_arrival_elements_vary_flag[j][2])
				    {
					    Flowerbounds->push_back(0.0);
					    Fupperbounds->push_back(math::LARGE);
					    Fdescriptions->push_back(prefix + " right boundary central body exclusion radius constraint");
				    }

				    //this constraint has derivatives with respect to x, y, and z
				    //only create a derivative entry with respect to an orbit element if that element is being varied
				    if (options->journey_arrival_elements_vary_flag[j][0]) //x
				    {
					    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
					    {
						    if ((*Xdescriptions)[entry].find("right boundary point x (km)") < 1024)
						    {
							    iGfun->push_back(Fdescriptions->size() - 1);
							    jGvar->push_back(entry);
							    stringstream EntryNameStream;
							    EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
							    Gdescriptions->push_back(EntryNameStream.str());
							    right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
							    right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
							    right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
						    }
					    }
				    }
				    if (options->journey_arrival_elements_vary_flag[j][1]) //y
				    {
					    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
					    {
						    if ((*Xdescriptions)[entry].find("right boundary point y (km)") < 1024)
						    {
							    iGfun->push_back(Fdescriptions->size() - 1);
							    jGvar->push_back(entry);
							    stringstream EntryNameStream;
							    EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
							    Gdescriptions->push_back(EntryNameStream.str());
							    right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
							    right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
							    right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
						    }
					    }
				    }
				    if (options->journey_arrival_elements_vary_flag[j][2]) //z
				    {
					    for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
					    {
						    if ((*Xdescriptions)[entry].find("right boundary point z (km)") < 1024)
						    {
							    iGfun->push_back(Fdescriptions->size() - 1);
							    jGvar->push_back(entry);
							    stringstream EntryNameStream;
							    EntryNameStream << "Derivative of " << prefix << " right boundary central body exclusion radius constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
							    Gdescriptions->push_back(EntryNameStream.str());
							    right_boundary_central_body_exclusion_radius_constraint_X_indices.push_back(entry);
							    right_boundary_central_body_exclusion_radius_constraint_G_indices.push_back(iGfun->size() - 1);
							    right_boundary_central_body_exclusion_radius_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
						    }
					    }
				    }
			    }
		    }

		    //variables and constraints necessary to process the arrival type
		    if (options->journey_arrival_type[j] == 2) //flyby with bounded v-infinity
		    {
			    Flowerbounds->push_back((options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1])*(options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1]) - 1);
			    Fupperbounds->push_back(0.0);
			    Fdescriptions->push_back(prefix + "arrival C3 constraint");

			    //Jacobian entry for a bounded v-infinity intercept
			    //this is a nonlinear constraint dependent on all previous variables
			    for (size_t entry = 0; entry < Xdescriptions->size(); ++entry)
			    {
				    iGfun->push_back(Fdescriptions->size() - 1);
				    jGvar->push_back(entry);
				    stringstream EntryNameStream;
				    EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
				    Gdescriptions->push_back(EntryNameStream.str());
			    }
		    }
	    }
    }

} /* namespace EMTG */

