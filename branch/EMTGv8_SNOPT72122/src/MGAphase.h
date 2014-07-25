/*
 * MGAphase.h
 *
 *  Created on: Jul 15, 2012
 *      Author: Jacob
 */

#ifndef MGAPHASE_H_
#define MGAPHASE_H_

#include <vector>

#include "phase.h"
#include "journey.h"
#include "missionoptions.h"
#include "universe.h"



namespace EMTG {

class MGA_phase: public EMTG::phase {
public:
	//constructor
	MGA_phase();
	MGA_phase(int j, int p, missionoptions* options);

	//destructor
	virtual ~MGA_phase();

	//evaluate function
	//return 0 if successful, 1 if failure
	int evaluate(double* X, int* Xindex, double* F, int* Findex, double* G, int* Gindex, int needG, double* current_epoch, double* current_state, double* current_deltaV, double* boundary1_state, double* boundary2_state, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options);

	//output function
	//return 0 if successful, 1 if failure
	int output(missionoptions* options, const double& launchdate, int j, int p,  EMTG::Astrodynamics::universe* Universe, int* eventcount);

	//bounds calculation function
	//return 0 if successful, 1 if failure
	int calcbounds(vector<double>* Xupperbounds, vector<double>* Xlowerbounds, vector<double>* Fupperbounds, vector<double>* Flowerbounds, vector<string>* Xdescriptions, vector<string>* Fdescriptions, vector<int>* iAfun, vector<int>* jAvar, vector<int>* iGfun, vector<int>* jGvar, vector<string>* Adescriptions, vector<string>* Gdescriptions, vector<double>* synodic_periods, int j, int p, EMTG::Astrodynamics::universe* Universe, missionoptions* options);

	//burn information
	//for a phase that begins with a departure or a flyby, there is one burn
	//for a phase that ends with an arrival, there is one additional burn
	vector<double> dVmag;
};

} /* namespace EMTG */
#endif /* MGAPHASE_H_ */
