/*
 * MGANDSMphase.cpp
 *
 *  Created on: April 22, 2013
 *      Author: Jacob
 */

#include <sstream>
#include <fstream>

#include "MGANDSMphase.h"
#include "Astrodynamics.h"
#include "Kepler_Lagrange_Laguerre_Conway_Der.h"
#include "missionoptions.h"
#include "mjd_to_mdyhms.h"
#include "EMTG_math.h"
#include "universe.h"
#include "EMTG_Matrix.h"

#include "SpiceUsr.h"



namespace EMTG {

MGA_NDSM_phase::MGA_NDSM_phase() {
//default constructor does nothing
}

MGA_NDSM_phase::MGA_NDSM_phase(int j, int p, missionoptions* options) {

	//must resize all data vectors to the correct length
	dV.resize(3, 1);

	match_point_state.resize(7);

	//size the vectors that will be used to calculate the b-plane
	V_infinity_in.resize(3, 1);
	V_infinity_out.resize(3, 1);
	BoundaryR.resize(3, 1);
	BoundaryV.resize(3, 1);

	//set the bodies
	boundary1_location_code = options->sequence[j][p];
	boundary2_location_code = options->sequence[j][p+1];

	current_mass_increment = 0.0;
	journey_initial_mass_increment_scale_factor = 1.0;

	//set up the forward and backward STMs
	Forward_STM.resize(1);
	Backward_STM.resize(1);
	Kepler_F_Forward.resize(1);
	Kepler_Fdot_Forward.resize(1);
	Kepler_G_Forward.resize(1);
	Kepler_Gdot_Forward.resize(1);
	Kepler_F_Backward.resize(1);
	Kepler_Fdot_Backward.resize(1);
	Kepler_G_Backward.resize(1);
	Kepler_Gdot_Backward.resize(1);
	Kepler_Fdotdot_Forward.resize(1);
	Kepler_Gdotdot_Forward.resize(1);
	Kepler_Fdotdot_Backward.resize(1);
	Kepler_Gdotdot_Backward.resize(1);
}

MGA_NDSM_phase::~MGA_NDSM_phase() {
	//destructor doesn't have to do anything
}

//evaluate function
//return 0 if successful, 1 if failure
int MGA_NDSM_phase::evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options) 
{
	//declare some local variables
	int errcode = 0;

	//******************************************************************
	//Steps 1-4: Process the left boundary condition
	this->process_left_boundary_condition(	X,
											Xindex,
											F, 
											Findex,
											G, 
											Gindex,
											needG,
											current_epoch, 
											current_state, 
											current_deltaV, 
											boundary1_state,
											boundary2_state, 
											j, 
											p,
											Universe, 
											options);
	
	//******************************************************************
	//Step 5: we need to know the state of the spacecraft at the right hand side (end) of the phase in order to propagate backward
	this->process_right_boundary_condition(	X,
											Xindex,
											F,
											Findex,
											G,
											Gindex,
											needG, 
											current_epoch, 
											current_state,
											current_deltaV, 
											boundary1_state, 
											boundary2_state, 
											j,
											p, 
											Universe,
											options);

	//******************************************************************
	//Step 6: propagate forward and back

	//Step 6.1: extract the burn index
	eta = X[*Xindex];
	++(*Xindex);

	//Step 6.2: propagate forward to the DSM
	double spacecraft_state_forward[7];
	spacecraft_state_forward[6] = state_at_beginning_of_phase[6];

	if (options->derivative_type > 1 && needG)
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_beginning_of_phase,
													spacecraft_state_forward,
													Universe->mu,
													Universe->LU,
													this->eta * this->TOF,
													this->Kepler_F_Forward[0],
													this->Kepler_Fdot_Forward[0],
													this->Kepler_G_Forward[0], 
													this->Kepler_Gdot_Forward[0],
													this->Kepler_Fdotdot_Forward[0],
													this->Kepler_Gdotdot_Forward[0],
													this->Forward_STM[0], 
													true);
	else
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_beginning_of_phase,
													spacecraft_state_forward,
													Universe->mu,
													Universe->LU,
													this->eta * this->TOF,
													this->Kepler_F_Forward[0],
													this->Kepler_Fdot_Forward[0],
													this->Kepler_G_Forward[0], 
													this->Kepler_Gdot_Forward[0],
													this->Kepler_Fdotdot_Forward[0],
													this->Kepler_Gdotdot_Forward[0],
													this->Forward_STM[0], 
													false);

	//Step 6.3: propagate backward
	double spacecraft_state_backward[7];
	spacecraft_state_backward[6] = state_at_end_of_phase[6];

	if (options->derivative_type > 1 && needG)
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_end_of_phase,
													spacecraft_state_backward,
													Universe->mu,
													Universe->LU,
													-((1 - this->eta) * this->TOF),
													this->Kepler_F_Backward[0],
													this->Kepler_Fdot_Backward[0],
													this->Kepler_G_Backward[0], 
													this->Kepler_Gdot_Backward[0],
													this->Kepler_Fdotdot_Backward[0],
													this->Kepler_Gdotdot_Backward[0],
													this->Backward_STM[0], 
													true);
	else
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(this->state_at_end_of_phase,
													spacecraft_state_backward,
													Universe->mu,
													Universe->LU,
													-((1 - this->eta) * this->TOF),
													this->Kepler_F_Backward[0],
													this->Kepler_Fdot_Backward[0],
													this->Kepler_G_Backward[0], 
													this->Kepler_Gdot_Backward[0],
													this->Kepler_Fdotdot_Backward[0],
													this->Kepler_Gdotdot_Backward[0],
													this->Backward_STM[0], 
													false);

	//Step 6.4: enforce match point constraint

	//Step 6.4.1 position constraints and dV calculation
	for (size_t k=0; k<3; ++k)
	{
		//position
		F[*Findex+k] = (spacecraft_state_backward[k] - spacecraft_state_forward[k]) / Universe->LU;

		//compute velocity
		dV(k) = spacecraft_state_backward[k+3] - spacecraft_state_forward[k+3];

		//match point state
		match_point_state[k] = spacecraft_state_backward[k];
		match_point_state[k+3] = spacecraft_state_backward[k+3];
	}
	(*Findex) += 3;

	//Step 6.4.2 compute dV magnitude and post-burn mass
	DSM_magnitude = dV.norm();
	*current_deltaV += DSM_magnitude;
	match_point_state[6] = spacecraft_state_forward[6] * exp(-DSM_magnitude * 1000 / (options->IspChem * options->g0));

	//mass
	F[*Findex] = (spacecraft_state_backward[6] - match_point_state[6])/(options->maximum_mass + journey_initial_mass_increment_scale_factor * current_mass_increment);
	++(*Findex);

	//Step 6.5: If required, compute the match point derivatives
#ifdef _MGANDSM_STM
	if (options->derivative_type > 0 && needG)
		calculate_match_point_derivatives(G, Gindex, j, p, options, Universe, spacecraft_state_forward, spacecraft_state_backward);
#endif

	//******************************************************************
	//Step 7: process the arrival, if applicable
	if (p == options->number_of_phases[j] - 1)
	{
		if (options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 4)
			dV_arrival_magnitude = 0.0;
		
		//note that "3" for journey_arrival_type indicates a "low-thrust rendezvous," which means we are there already and we don't need to do anything
		else
		{
			//compute the arrival deltaV
			if (boundary2_location_code > 0) //ending at body
				dV_arrival_magnitude = process_arrival(	state_at_end_of_phase+3,
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
				dV_arrival_magnitude = process_arrival(	state_at_end_of_phase+3,
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

			//apply the arrival burn
			state_at_end_of_phase[6] *= exp(-dV_arrival_magnitude * 1000/ ((options->IspChem > 0 ? options->IspChem : 1e-6)* options->g0));
			*current_deltaV += dV_arrival_magnitude;

			//if this is a terminal intercept, we need to compute derivatives
			if (needG && options->journey_arrival_type[j] == 2)
			{
				double C3_desired = options->journey_final_velocity[j][1] * options->journey_final_velocity[j][1];
				//derivative with respect to x component of terminal velocity
				G[terminal_velocity_constraint_G_indices[0]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[0] * X[terminal_velocity_constraint_X_indices[0]] / C3_desired;

				//derivative with respect to y component of terminal velocity
				G[terminal_velocity_constraint_G_indices[1]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[1] * X[terminal_velocity_constraint_X_indices[1]] / C3_desired;

				//derivative with respect to z component of terminal velocity
				G[terminal_velocity_constraint_G_indices[2]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[2] * X[terminal_velocity_constraint_X_indices[2]] / C3_desired;
			}
			/*if ( needG && options->journey_arrival_declination_constraint_flag[j] && (options->journey_arrival_type[j] == 0 || options->journey_arrival_type[j] == 2) ) //intercept with bounded v-infinity or orbit insertion
			{
				double Vx = X[arrival_declination_constraint_X_indices[0]];
				double Vy = X[arrival_declination_constraint_X_indices[1]];
				double Vz = X[arrival_declination_constraint_X_indices[2]];
				double A = sqrt(Vx*Vx + Vy*Vy);
				double B = Vx*Vx + Vy*Vy + Vz*Vz;

				//derivative with respect to x component of terminal velocity
				G[arrival_declination_constraint_G_indices[0]] = arrival_declination_constraint_X_scale_ranges[0] * Vx*Vz / (A*B) / options->journey_arrival_declination_bounds[j][1];

				//derivative with respect to y component of terminal velocity
				G[arrival_declination_constraint_G_indices[1]] = arrival_declination_constraint_X_scale_ranges[0] * Vy*Vz / (A*B) / options->journey_arrival_declination_bounds[j][1];

				//derivative with respect to z component of terminal velocity
				G[arrival_declination_constraint_G_indices[2]] = arrival_declination_constraint_X_scale_ranges[0] * -A/B / options->journey_arrival_declination_bounds[j][1];
			}*/
            else if (needG && options->journey_arrival_type[j] == 6)
			{
				double r = math::norm(boundary2_state, 3) / Universe->LU;
				double denominator = r * r * r * Universe->LU * Universe->LU;
				double LUTU2 = (Universe->TU / Universe->LU) * (Universe->TU / Universe->LU);
				double sqLU = sqrt(Universe->LU);

				//derivative with respect to x component of position
				G[escape_constraint_G_indices[0]] = escape_constraint_X_scale_ranges[0] * X[escape_constraint_X_indices[0]] / denominator;

				//derivative with respect to x component of position
				G[escape_constraint_G_indices[1]] = escape_constraint_X_scale_ranges[1] * X[escape_constraint_X_indices[1]] / denominator;

				//derivative with respect to x component of position
				G[escape_constraint_G_indices[2]] = escape_constraint_X_scale_ranges[2] * X[escape_constraint_X_indices[2]] / denominator;

				//derivative with respect to x component of terminal velocity
				G[escape_constraint_G_indices[3]] = escape_constraint_X_scale_ranges[3] * X[escape_constraint_X_indices[3]] * LUTU2;

				//derivative with respect to y component of terminal velocity
				G[escape_constraint_G_indices[4]] = escape_constraint_X_scale_ranges[4] * X[escape_constraint_X_indices[4]] * LUTU2;

				//derivative with respect to z component of terminal velocity
				G[escape_constraint_G_indices[5]] = escape_constraint_X_scale_ranges[5] * X[escape_constraint_X_indices[5]] * LUTU2;
			}
		}
			
	}
	
	//******************************************************************
	//Step 8: update the current epoch
	*current_epoch += TOF;

	//******************************************************************
	//Step 9: update the current state
	for (int k = 0; k < 7; ++k)
		current_state[k] = state_at_end_of_phase[k];

	return 0;
}


//bounds calculation function
//return 0 if successful, 1 if failure
int MGA_NDSM_phase::calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options)
{
	//this function calculates the upper and lower bounds for the decision and constraint vectors for MGA-NDSM
	//create a prefix string with journey and phase information
	stringstream prefixstream;
	prefixstream << "j" << j << "p" << p << ": ";
	string prefix = prefixstream.str();
	int first_X_entry_in_phase = Xupperbounds->size();

	//**************************************************************************
	//calculate bounds on variables and constraints governing the left boundary
	calcbounds_left_boundary(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	//**************************************************************************
	//if EMTG is choosing an input power or Isp for the phase (for REP/NEP models), then this information must be encoded
	calcbounds_phase_thruster_parameters(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);

	//**************************************************************************
	//next, we need to encode the phase flight time
	calcbounds_flight_time(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, synodic_periods, j, p, Universe, options);

	//**************************************************************************
	//calculate bounds on variables and constraints governing the right boundary
	calcbounds_right_boundary(prefix, first_X_entry_in_phase, Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, Universe, options);


	//**************************************************************************
	//next, we need to include the burn index
	Xlowerbounds->push_back(0.05);
	Xupperbounds->push_back(0.95);
	Xdescriptions->push_back(prefix + "burn index");

	//**************************************************************************
	//finally, we encode the match point continuity constraints and their Jacobian entries,
	//noting that every patch point constraint in the phase has a derivative with respect to every variable in the phase
	//in addition, the patch point constraints have a derivative with respect to the previous phase's arrival mass
	//and the patch point constraints have a derivatives with respect to all previous time variables, including the launch date

	G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass.resize(7);
	G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass.resize(7);

	vector<string> statename;
	statename.push_back("x");
	statename.push_back("y");
	statename.push_back("z");
	statename.push_back("xdot");
	statename.push_back("ydot");
	statename.push_back("zdot");
	statename.push_back("m");

	Flowerbounds->push_back(-math::SMALL);
	Fupperbounds->push_back(math::SMALL);
	Fdescriptions->push_back(prefix + "match point x");
	for (size_t entry =  first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
	{
		iGfun->push_back(Fdescriptions->size() - 1);
		jGvar->push_back(entry);
		stringstream EntryNameStream;
		EntryNameStream << "Derivative of " << prefix << " patch point x constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
		Gdescriptions->push_back(EntryNameStream.str());
	}
	for (int entry =  first_X_entry_in_phase; entry >= 0; --entry)
	{
		if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
		{
			iGfun->push_back(Fdescriptions->size() - 1);
			jGvar->push_back(entry);
			stringstream EntryNameStream;
			EntryNameStream << "Derivative of " << prefix << " patch point x constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
			Gdescriptions->push_back(EntryNameStream.str());
			break;
		}
	}
	//derivative with respect to times and epochs in previous journeys/phases
	for (int pj = 0; pj <= j; ++pj)
	{
		for (int pp = 0; pp < (pj == j ? p : options->number_of_phases[pj]); ++pp)
		{
			stringstream pprefix_stream;
			pprefix_stream << "j" << pj << "p" << pp;
			string pprefix = pprefix_stream.str();

			for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
			{
				if ( (*Xdescriptions)[entry].find(pprefix) < 1024 && ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024) )
				{
					iGfun->push_back(Fdescriptions->size() - 1);
					jGvar->push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of " << prefix << " patch point x constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					Gdescriptions->push_back(EntryNameStream.str());
				}
			}
		}
	}

	Flowerbounds->push_back(-math::SMALL);
	Fupperbounds->push_back(math::SMALL);
	Fdescriptions->push_back(prefix + "match point y");
	for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
	{
		iGfun->push_back(Fdescriptions->size() - 1);
		jGvar->push_back(entry);
		stringstream EntryNameStream;
		EntryNameStream << "Derivative of " << prefix << " patch point y constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
		Gdescriptions->push_back(EntryNameStream.str());
	}
	for (int entry =  first_X_entry_in_phase; entry >= 0; --entry)
	{
		if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
		{
			iGfun->push_back(Fdescriptions->size() - 1);
			jGvar->push_back(entry);
			stringstream EntryNameStream;
			EntryNameStream << "Derivative of " << prefix << " patch point y constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
			Gdescriptions->push_back(EntryNameStream.str());
			break;
		}
	}
	//derivative with respect to times and epochs in previous journeys/phases
	for (int pj = 0; pj <= j; ++pj)
	{
		for (int pp = 0; pp < (pj == j ? p : options->number_of_phases[pj]); ++pp)
		{
			stringstream pprefix_stream;
			pprefix_stream << "j" << pj << "p" << pp;
			string pprefix = pprefix_stream.str();

			for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
			{
				if ( (*Xdescriptions)[entry].find(pprefix) < 1024 && ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024) )
				{
					iGfun->push_back(Fdescriptions->size() - 1);
					jGvar->push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of " << prefix << " patch point y constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					Gdescriptions->push_back(EntryNameStream.str());
				}
			}
		}
	}

	Flowerbounds->push_back(-math::SMALL);
	Fupperbounds->push_back(math::SMALL);
	Fdescriptions->push_back(prefix + "match point z");
	for (size_t entry =  first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
	{
		iGfun->push_back(Fdescriptions->size() - 1);
		jGvar->push_back(entry);
		stringstream EntryNameStream;
		EntryNameStream << "Derivative of " << prefix << " patch point z constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
		Gdescriptions->push_back(EntryNameStream.str());
	}
	for (int entry =  first_X_entry_in_phase; entry >= 0; --entry)
	{
		if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
		{
			iGfun->push_back(Fdescriptions->size() - 1);
			jGvar->push_back(entry);
			stringstream EntryNameStream;
			EntryNameStream << "Derivative of " << prefix << " patch point z constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
			Gdescriptions->push_back(EntryNameStream.str());
			break;
		}
	}
	//derivative with respect to times and epochs in previous journeys/phases
	for (int pj = 0; pj <= j; ++pj)
	{
		for (int pp = 0; pp < (pj == j ? p : options->number_of_phases[pj]); ++pp)
		{
			stringstream pprefix_stream;
			pprefix_stream << "j" << pj << "p" << pp;
			string pprefix = pprefix_stream.str();

			for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
			{
				if ( (*Xdescriptions)[entry].find(pprefix) < 1024 && ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024) )
				{
					iGfun->push_back(Fdescriptions->size() - 1);
					jGvar->push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of " << prefix << " patch point z constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					Gdescriptions->push_back(EntryNameStream.str());
				}
			}
		}
	}

	Flowerbounds->push_back(-math::SMALL);
	Fupperbounds->push_back(math::SMALL);
	Fdescriptions->push_back(prefix + "match point m");
	for (size_t entry =  first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
	{
		iGfun->push_back(Fdescriptions->size() - 1);
		jGvar->push_back(entry);
		stringstream EntryNameStream;
		EntryNameStream << "Derivative of " << prefix << " patch point m constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
		Gdescriptions->push_back(EntryNameStream.str());

		if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
		{
			G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass[6] = Gdescriptions->size() - 1;
		}
	}
	for (int entry =  first_X_entry_in_phase; entry >= 0; --entry)
	{
		if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
		{
			iGfun->push_back(Fdescriptions->size() - 1);
			jGvar->push_back(entry);
			stringstream EntryNameStream;
			EntryNameStream << "Derivative of " << prefix << " patch point m constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
			Gdescriptions->push_back(EntryNameStream.str());

			G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass[6] = Gdescriptions->size() - 1;
			break;
		}
	}
	//derivative with respect to times and epochs in previous journeys/phases
	for (int pj = 0; pj <= j; ++pj)
	{
		for (int pp = 0; pp < (pj == j ? p : options->number_of_phases[pj]); ++pp)
		{
			stringstream pprefix_stream;
			pprefix_stream << "j" << pj << "p" << pp;
			string pprefix = pprefix_stream.str();

			for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
			{
				if ( (*Xdescriptions)[entry].find(pprefix) < 1024 && ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024) )
				{
					iGfun->push_back(Fdescriptions->size() - 1);
					jGvar->push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of " << prefix << " patch point m constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					Gdescriptions->push_back(EntryNameStream.str());
				}
			}
		}
	}
	
	//***************************************************************************
	//Construct a helper array "rectangular prism" of G indices for match point derivatives with respect to state variables within the phase
	vector < vector <int> > timeslice;
	vector <int> scanline;
	vector<string> constraint_type;
	constraint_type.push_back(" patch point x constraint");
	constraint_type.push_back(" patch point y constraint");
	constraint_type.push_back(" patch point z constraint");
	constraint_type.push_back(" patch point m constraint");

	for (int step = 0; step < 2; ++step)
	{
		timeslice.clear();
		
		//loop over constraint type dimension
		for (int c = 0; c < 4; ++c)
		{
			scanline.clear();

			string constraintname = prefix + constraint_type[c];

			for (int entry = 0; entry < Gdescriptions->size(); ++entry)
			{
				if ( (*Gdescriptions)[entry].find(constraintname) < 1024)
				{
					if (step == 0) //derivatives with respect to initial velocity
					{
						//the first step of the first phase of a journey is abnormal because instead of encoding an XYZ vector, we have encoded a magnitude and two angles
						if (p == 0)
						{
							if ( (*Gdescriptions)[entry].find("magnitude of outgoing velocity asymptote") < 1024)
							{
								scanline.push_back(entry);     //derivative with respect to magnitude
								scanline.push_back(entry + 1); //derivative with respect to RA
								scanline.push_back(entry + 2); //derivative with respect to DEC
							}
						}
						//otherwise look for an XYZ initial velocity increment
						else
						{
							if ( (*Gdescriptions)[entry].find("initial velocity increment x") < 1024)
							{
								scanline.push_back(entry);     //derivative with respect to initial velocity increment x
								scanline.push_back(entry + 1); //derivative with respect to initial velocity increment y
								scanline.push_back(entry + 2); //derivative with respect to initial velocity increment z
							}
						}
					}
					//derivatives with respect to final velocity - these do not exist for a terminal rendezvous phase
					else if (step == 1 && (!(p == options->number_of_phases[j] - 1 && options->journey_arrival_type[j] == 4)))
					{
						if ( (*Gdescriptions)[entry].find("terminal velocity increment x") < 1024)
						{
							scanline.push_back(entry);	   //derivative with respect to terminal velocity increment x
							scanline.push_back(entry + 1); //derivative with respect to terminal velocity increment y
							scanline.push_back(entry + 2); //derivative with respect to terminal velocity increment z
						}
					}
				} //end if ( (*Gdescriptions)[entry].find(constraintname) )
			} //end loop over Gdescriptions
			
			timeslice.push_back(scanline);
			
		} //end loop over constraint entries
		match_point_constraint_G_indices.push_back(timeslice);
	
	} //end loop over time steps

	//***************************************************************************
	//if this is the last phase, encode any constraints for the arrival processing
	if (p == options->number_of_phases[j] - 1)
	{
		if (options->journey_arrival_type[j] == 2) //intercept with bounded v-infinity
		{
			Flowerbounds->push_back((options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1])*(options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1]) - 1);
			Fupperbounds->push_back(0.0);
			Fdescriptions->push_back(prefix + "arrival C3 constraint");

			//Jacobian entry for a bounded v-infinity intercept
			//this is a nonlinear constraint dependent only on the terminal velocity vector for this phase
			for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
			{
				if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
				{
					iGfun->push_back(Fdescriptions->size() - 1);
					jGvar->push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					Gdescriptions->push_back(EntryNameStream.str());
					terminal_velocity_constraint_G_indices.push_back(iGfun->size() - 1);
					terminal_velocity_constraint_X_indices.push_back(entry);
					terminal_velocity_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
				}
			}
		}

		if ( options->journey_arrival_declination_constraint_flag[j] && (options->journey_arrival_type[j] == 0 || options->journey_arrival_type[j] == 2) ) //intercept with bounded v-infinity or orbit insertion
		{
			Flowerbounds->push_back(options->journey_arrival_declination_bounds[j][0] / options->journey_arrival_declination_bounds[j][1] - 1.0);
			Fupperbounds->push_back(0.0);
			Fdescriptions->push_back(prefix + "arrival declination constraint");

			//Jacobian entry for arrival declination constraint
			//this is a nonlinear constraint dependent only on the terminal velocity vector for this phase
			for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
			{
				if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
				{
					iGfun->push_back(Fdescriptions->size() - 1);
					jGvar->push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					Gdescriptions->push_back(EntryNameStream.str());
					arrival_declination_constraint_G_indices.push_back(iGfun->size() - 1);
					arrival_declination_constraint_X_indices.push_back(entry);
					arrival_declination_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
				}
			}
		}

        if (options->journey_arrival_type[j] == 6) //enforce escape constraint, E = 0
		{
			Flowerbounds->push_back(-math::SMALL);
			Fupperbounds->push_back(math::LARGE);
			Fdescriptions->push_back(prefix + "arrival escape condition (zero energy constraint)");

			//Jacobian entry for zero energy condition
			//this is a nonlinear constraint dependent on the terminal velocity vector
			//note that it is NOT dependent on position because these phases always end at the SOI and SOI distance is constant
			for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
			{
				if ((*Xdescriptions)[entry].find("terminal velocity") < 1024 || (*Xdescriptions)[entry].find("right boundary point") < 1024)
				{
					iGfun->push_back(Fdescriptions->size() - 1);
					jGvar->push_back(entry);
					stringstream EntryNameStream;
					EntryNameStream << "Derivative of " << prefix << " escape constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
					Gdescriptions->push_back(EntryNameStream.str());
					escape_constraint_G_indices.push_back(iGfun->size() - 1);
					escape_constraint_X_indices.push_back(entry);
					escape_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
				}
			}
		}
	}
		
	return 0;
}

//output function
//return 0 if successful, 1 if failure
int MGA_NDSM_phase::output(missionoptions* options, const double& launchdate, int j, int p, EMTG::Astrodynamics::universe* Universe, int* eventcount)
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
	else if (j == 0 && boundary1_location_code > 0 && options->LV_type >= 0)
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

	//Step 2: all phases have at least two events: departure/flyby and burn
	//*****************************************************************************
	//first let's print the departure/flyby
	
	for (int k = 0; k < 3; ++k)
		dVdeparture[k] = V_infinity_out(k);

	this->write_summary_line(options,
							Universe,
							eventcount,
							phase_start_epoch / 86400.0,
							event_type,
							boundary1_name,
							0,
							(p > 0 ? flyby_altitude : Bradius),
							(Btheta),
							(p > 0 ? flyby_turn_angle : -1),
							RA_departure,
							DEC_departure,
							C3_departure,
							state_at_beginning_of_phase,
							dVdeparture,
							empty_vector,
							(p == 0 ? dV_departure_magnitude : flyby_outgoing_v_infinity),
							-1,
							initial_Isp,
							-1,
							0,
							0,
							0);



	//*****************************************************************************
	//Next, all phases have a deep-space maneuver (DSM). Coast until the DSM
	double timestep = eta * TOF / options->num_timesteps;
	double output_state[7];
	output_state[6] = state_at_beginning_of_phase[6];
	
	for (int step = 0; step < options->num_timesteps; ++step)
	{
		//compute the current epoch
		double epoch = phase_start_epoch + timestep * (step + 0.5);

		//propagate the spacecraft
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(state_at_beginning_of_phase,
													output_state,
													Universe->mu,
													Universe->LU,
													epoch - phase_start_epoch, 
													this->Kepler_F_Current, 
													this->Kepler_Fdot_Current,
													this->Kepler_G_Current,
													this->Kepler_Gdot_Current, 
													this->Kepler_Fdotdot_Current,
													this->Kepler_Gdotdot_Current,
													this->Current_STM,
													false);

		//write the summary line
		this->write_summary_line(options,
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
								empty_vector,
								empty_vector,
								0,
								-1,
								-1,
								-1,
								0,
								0,
								0);
	}

	//then print the DSM
	double dVout[3];
	dVout[0] = dV(0);
	dVout[1] = dV(1);
	dVout[2] = dV(2);

	this->write_summary_line(options,
							Universe,
							eventcount,
							(phase_start_epoch + eta * TOF) / 86400.0,
							"chem_burn",
							"deep-space",
							0,
							-1,
							-1,
							-1,
							atan2(dV(0), dV(1)),
							asin(dV(2) / DSM_magnitude),
							0,
							match_point_state.data(),
							dVout,
							empty_vector,
							DSM_magnitude,
							-1,
							options->IspChem,
							-1,
							0,
							0,
							0);

	//now print the coast after the DSM
	timestep = (1 - eta) * TOF / options->num_timesteps;
	output_state[6] = match_point_state[6];

	for (int step = 0; step < options->num_timesteps; ++step)
	{
		//compute the current epoch
		double epoch = (phase_start_epoch + eta * TOF + timestep * (step + 0.5));

		//propagate the spacecraft
		Kepler::Kepler_Lagrange_Laguerre_Conway_Der(match_point_state.data(),
													output_state,
													Universe->mu,
													Universe->LU,
													timestep * (step + 0.5), 
													this->Kepler_F_Current,
													this->Kepler_Fdot_Current,
													this->Kepler_G_Current, 
													this->Kepler_Gdot_Current, 
													this->Kepler_Fdotdot_Current,
													this->Kepler_Gdotdot_Current, 
													this->Current_STM,
													false);

		//write the summary line
		this->write_summary_line(options,
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
								empty_vector,
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
		else if (options->journey_arrival_type[j] == 1)
			event_type = "rendezvous";
		else if (options->journey_arrival_type[j] == 2)
			event_type = "intercept";
		else if (options->journey_arrival_type[j] == 3)
			event_type = "LT_rndzvs";
		else if (options->journey_arrival_type[j] == 5 || options->journey_arrival_type[j] == 4)
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
		else if (options->journey_arrival_type[j] == 4 || options->journey_arrival_type[j] == 3)
		{
			dV_arrival_mag = 0;
			dVarrival[0] = 0;
			dVarrival[1] = 0;
			dVarrival[2] = 0;
		}
		else
		{
			dV_arrival_mag = dV_arrival_magnitude;
		}

		this->write_summary_line(options,
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

	return 0;
}

	//function to calculate the patch point derivatives
	int MGA_NDSM_phase::calculate_match_point_derivatives(double* G, int* Gindex, int j, int p, missionoptions* options, EMTG::Astrodynamics::universe* Universe, double* match_point_state_forward, double* match_point_state_backward)
	{
		//forward derivative with respect to the initial v-infinity
		if (p == 0)
		{
			//first phase, match point forward derivative
			double cosRA = cos(RA_departure);
			double sinRA = sin(RA_departure);
			double cosDEC = cos(DEC_departure);
			double sinDEC = cos(DEC_departure);
			double vinf = sqrt(C3_departure);
				
			for (int c = 0; c < 4; ++c)
			{
				double dFdvx = Forward_STM[0](c, 3);
				double dFdvy = Forward_STM[0](c, 4);
				double dFdvz = Forward_STM[0](c, 5);

				//derivatives with respect to v-infinity
				G[match_point_constraint_G_indices[0][c][0]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][c][0]]] 
																* (cosRA*cosDEC*dFdvx + sinRA*cosDEC*dFdvy + sinDEC*dFdvz) / Universe->LU;

				//derivatives with respect to RA_departure
				G[match_point_constraint_G_indices[0][c][1]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][c][1]]] 
																* vinf * (-sinRA*cosDEC*dFdvx + cosRA*cosDEC*dFdvy) / Universe->LU;

				//derivatives with respect to DEC_departure
				G[match_point_constraint_G_indices[0][c][2]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][c][2]]] 
																* vinf * (-cosRA*sinDEC*dFdvx - sinRA*sinDEC*dFdvy + cosDEC*dFdvz) / Universe->LU;
			}
		}
		//if this is not the first phase, then we have derivatives with respect to the initial velocity increment
		else
		{
			for (int c = 0; c < 3; ++c)
			{
				G[match_point_constraint_G_indices[0][c][0]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][c][0]]] * Forward_STM[0](c, 3) / Universe->LU;
				G[match_point_constraint_G_indices[0][c][1]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][c][1]]] * Forward_STM[0](c, 4) / Universe->LU;
				G[match_point_constraint_G_indices[0][c][2]] = -options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[0][c][2]]] * Forward_STM[0](c, 5) / Universe->LU;
			}
		}

		//derivatives with respect to the terminal velocity increment
		if (!(p == options->number_of_phases[j] - 1 && (options->journey_arrival_type[j] == 4)))
		{
			for (int c = 0; c < 3; ++c)
			{
				G[match_point_constraint_G_indices[1][c][0]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[1][c][0]]] * Backward_STM[0](c, 3) / Universe->LU;
				G[match_point_constraint_G_indices[1][c][1]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[1][c][1]]] * Backward_STM[0](c, 4) / Universe->LU;
				G[match_point_constraint_G_indices[1][c][2]] = options->X_scale_ranges[options->jGvar[match_point_constraint_G_indices[1][c][2]]] * Backward_STM[0](c, 5) / Universe->LU;
			}
		}

		return 0;
	}
} /* namespace EMTG */
