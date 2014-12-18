/*
 * MGALTphase.h
 *
 *  Created on: Jul 15, 2012
 *      Author: Jacob
 */

#ifndef MGALTPHASE_H_
#define MGALTPHASE_H_

#include <vector>

#include "TwoPointShootingPhase.h"
#include "journey.h"
#include "missionoptions.h"
#include "universe.h"
#include "STM.h"

namespace EMTG {

    class MGA_LT_phase : public EMTG::TwoPointShootingPhase {
	public:
		//constructor
		MGA_LT_phase();
        MGA_LT_phase(const int& j, const int& p, const missionoptions& options);

		//destructor
		virtual ~MGA_LT_phase();

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
        void output(missionoptions* options,
            const double& launchdate,
            const int& j,
            const int& p,
            EMTG::Astrodynamics::universe* Universe,
            int* eventcount);

		//bounds calculation function
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

        //top-level function to calculate the derivatives of the distance constraints
        void calculate_distance_constraint_derivatives( double* G,
                                                        const int& j,
                                                        const int& p,
                                                        missionoptions* options,
                                                        EMTG::Astrodynamics::universe* Universe);

		//top-level function to calculate the match point derivatives
		void calculate_match_point_derivatives(	double* G,
												const int& j, 
												const int& p,
												missionoptions* options, 
												EMTG::Astrodynamics::universe* Universe);

		//function to calculate the derivative of a match point constraint with respect to a decision variable in the forward propagation
		void calculate_match_point_forward_propagation_derivatives(	double* G,
																	const int& j, 
																	const int& p,
																	missionoptions* options, 
																	EMTG::Astrodynamics::universe* Universe,
																	const int& step,
																	const int& stepnext,
																	double& dxdu,
																	double& dydu,
																	double& dzdu,
																	double& dxdotdu,
																	double& dydotdu,
																	double& dzdotdu,
																	double& dmdu,
																	double& dtdu,
																	double& dtotal_available_thrust_time_du,
																	double& dPdu);

		//function to calculate the derivative of a match point constraint with respect to a decision variable in the backward propagation
		void calculate_match_point_backward_propagation_derivatives(double* G,
																	const int& j, 
																	const int& p,
																	missionoptions* options, 
																	EMTG::Astrodynamics::universe* Universe,
																	const int& backstep,
																	const int& stepnext,
																	double& dxdu,
																	double& dydu,
																	double& dzdu,
																	double& dxdotdu,
																	double& dydotdu,
																	double& dzdotdu,
																	double& dmdu,
																	double& dtdu,
																	double& dtotal_available_thrust_time_du,
																	double& dPdu);
	};

} /* namespace EMTG */
#endif /* MGALTPHASE_H_ */
