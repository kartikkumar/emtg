/*
 * MGANDSMphase.h
 *
 *  Created on: April 22, 2013
 *      Author: Jacob
 */

#ifndef MGANDSMPHASE_H_
#define MGANDSMPHASE_H_

#include <vector>

#include "phase.h"
#include "journey.h"
#include "missionoptions.h"
#include "universe.h"



namespace EMTG {

class MGA_NDSM_phase: public EMTG::phase {
public:
	//constructor
	MGA_NDSM_phase();
	MGA_NDSM_phase(int j, int p, missionoptions* options);

	//destructor
	virtual ~MGA_NDSM_phase();

	//evaluate function
	//return 0 if successful, 1 if failure
	int evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options);

	//output function
	//return 0 if successful, 1 if failure
	int output(missionoptions* options, const double& launchdate, int j, int p, EMTG::Astrodynamics::universe* Universe, int* eventcount);

	//bounds calculation function
	//return 0 if successful, 1 if failure
	int calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options);

	//function to calculate the match point derivatives
	int calculate_match_point_derivatives(double* G, int* Gindex, int j, int p, missionoptions* options, EMTG::Astrodynamics::universe* Universe, double* match_point_state_forward, double* match_point_state_backward);

	//time information
	double eta;

	//state information
	vector<double> match_point_state;
	math::Matrix<double> dV;
	double DSM_magnitude;
};

} /* namespace EMTG */
#endif /* MGANDSMPHASE_H_ */
