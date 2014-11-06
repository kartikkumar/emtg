/*
 * FBLTphase.h
 *
 *  Created on: September 17, 2012
 *      Author: Jacob
 */

#ifndef FBLTPHASE_H_
#define FBLTPHASE_H_

#include <vector>

#include "TwoPointShootingPhase.h"
#include "journey.h"
#include "missionoptions.h"
#include "rk8713M.h"
#include "equations_of_motion.h"
#include "universe.h"

namespace EMTG {

    class FBLT_phase : public EMTG::TwoPointShootingPhase {
public:
	//constructor
	FBLT_phase();
	FBLT_phase(int j, int p, missionoptions* options);

	//destructor
	virtual ~FBLT_phase();

	//evaluate function
	//return 0 if successful, 1 if failure
	int evaluate(double* X,
                int* Xindex,
                double* F,
                int* Findex, 
                double* G, 
                int* Gindex,
                int needG, 
                double* current_epoch, 
                double* current_state,
                double* current_deltaV,
                double* boundary1_state,
                double* boundary2_state, 
                int j, 
                int p, 
                EMTG::Astrodynamics::universe* Universe, 
                missionoptions* options);

	//output function
	//return 0 if successful, 1 if failure
	int output(missionoptions* options,
                const double& launchdate,
                int j,
                int p,
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
                    int j,
                    int p,
                    EMTG::Astrodynamics::universe* Universe,
                    missionoptions* options);

    //methods to output a .e ephemeris file
    void write_ephemeris_file(const missionoptions& options,
                                const EMTG::Astrodynamics::universe& Universe,
                                const double& launch_epoch,
                                const double& journey_starting_epoch,
                                vector< vector<string> >& output_line_array,
                                const int& j,
                                const int& p);

    void propagate_forward_ephemeris(const missionoptions& options,
                                    const EMTG::Astrodynamics::universe& Universe,
                                    const double& launch_epoch,
                                    double& current_epoch,
                                    double* current_state,
                                    const int& control_step,
                                    vector< vector<string> >& output_line_array,
                                    const int& j,
                                    const int& p,
                                    const double& propagation_time,
                                    const double& journey_starting_epoch);

	void calculate_match_point_derivatives(double* G,
		int * Gindex,
		const int & j,
		const int & p,
		std::vector < EMTG::math::Matrix< double > > & STM_archive_forward,
		std::vector < EMTG::math::Matrix< double > > & STM_archive_backward,
		missionoptions * options,
		EMTG::Astrodynamics::universe * Universe);

	//FBLT STMs
	std::vector < EMTG::math::Matrix< double > > STM_archive_forward;
	EMTG::math::Matrix< double > initial_coast_STM;
	std::vector < EMTG::math::Matrix< double > > STM_archive_backward;
	EMTG::math::Matrix< double > terminal_coast_STM;

	//FBLT STM dimension information
	int STMrows;
	int STMcolumns;
	int num_states;
	
	//time information
	vector <double> event_epochs;

	//state information
	vector<double> match_point_state;

	//integrator
	EMTG::integration::rk8713M *integrator;

    //dummy controller pointer
    void* DummyControllerPointer;
};

} /* namespace EMTG */
#endif /* FBLTPHASE_H_ */
