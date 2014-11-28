/*
 * MGADSMphase.h
 *
 *  Created on: Jul 15, 2012
 *      Author: Jacob
 */

#ifndef MGADSMPHASE_H_
#define MGADSMPHASE_H_

#include <vector>

#include "phase.h"
#include "journey.h"
#include "missionoptions.h"
#include "universe.h"

namespace EMTG {

class MGA_DSM_phase: public EMTG::phase {
public:
	//constructor
	MGA_DSM_phase();
    MGA_DSM_phase(const int& j, const int& p, const missionoptions& options);

	//destructor
	virtual ~MGA_DSM_phase();

	//evaluate function
	//return 0 if successful, 1 if failure
    int evaluate(const double* X,
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
        missionoptions* options);

	//output function
	//return 0 if successful, 1 if failure
    int output(missionoptions* options,
        const double& launchdate,
        const int& j,
        const int& p,
        EMTG::Astrodynamics::universe* Universe,
        int* eventcount);

	//bounds calculation function
	//return 0 if successful, 1 if failure
    void calcbounds(vector<double>* Xupperbounds,
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
                    missionoptions* options);

	//function to create an initial guess for another mission type
	virtual void create_initial_guess(	const int& desired_mission_type,
										const bool& VSI,
										double& current_epoch,
										const int& j,
										const int& p,
										vector<double>& NewX,
										int& NewXIndex,
										const vector<string>& NewXDescriptions,
										const missionoptions& options,
                                        const Astrodynamics::universe& Universe);

	//GMAT output methods
	void output_GMAT_fueltank_and_thruster(int& j, int& p, vector<EMTG::Astrodynamics::body>& missionbodies, int& index_body_visited, std::ofstream& GMATfile);
	void output_GMAT_burn_objects(int& j, int& p, std::ofstream& GMATfile);
	void output_GMAT_create_interphase_control_variables(int& j, int& p, missionoptions& options, std::ofstream& GMATfile);
	void output_GMAT_inter_phase_control_initial_guess(int& j, int& p, missionoptions& options, std::ofstream& GMATfile);

	//burn information
	//all phases have at least one burn
	//phases beginning wih a departure have an additional burn
	//phases ending in an arrival have an additional burn
	vector<double> dVmag;

	//burn index
	double burn_index;
	double midcourse_dV[3];

	//time information
	double time_before_burn;
	double time_after_burn;

	//flyby information
	double b_plane_insertion_angle;

	//additional state information
	double state_before_burn[7];
	double state_after_burn[7];
	double lambert_v1[3];
	double lambert_v2[3];

	//step information for outputs
	int steps_before_burn;
	double timestep_length_before_burn;
	double timestep_length_after_burn;
};

} /* namespace EMTG */
#endif /* MGADSMPHASE_H_ */
