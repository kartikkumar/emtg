//header file for EMTG frame class
//Jacob Englander 1/3/2013

#pragma once

#include "EMTG_Matrix.h"

namespace EMTG { namespace Astrodynamics {

	class frame
	{
	public:
		//constructor
		frame(void);

		frame(const double& alpha0_in, const double& alphadot_in, const double& delta0_in, const double& deltadot_in, const double& W_in, const double& Wdot_in);

		//destructor
		virtual ~frame(void);

		//methods
		void initialize();
		void initialize(const double& alpha0_in, const double& alphadot_in, const double& delta0_in, const double& deltadot_in, const double& W_in, const double& Wdot_in);
		double convert_JED_to_TDB(const double& JEDepoch);
		void construct_rotation_matrices(const double& JEDepoch);

		//fields
		double alpha0, alphadot, delta0, deltadot, W, Wdot;
		math::Matrix<double> R_from_local_to_ICRF;
		math::Matrix<double> R_from_ICRF_to_local;	
	};

}} //close namespace