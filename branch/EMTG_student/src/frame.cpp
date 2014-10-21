//EMTG frame class
//Jacob Englander 1/3/2013

#include <cmath>

#include "frame.h"
#include "EMTG_math.h"
#include "EMTG_time_utilities.h"


namespace EMTG { namespace Astrodynamics {
	//default constructor
	frame::frame(void)
	{
	}

	//constructor with data
	frame::frame(const double& alpha0_in, const double& alphadot_in, const double& delta0_in, const double& deltadot_in, const double& W_in, const double& Wdot_in)
	{
		initialize(alpha0_in, alphadot_in, delta0_in, deltadot_in, W_in, Wdot_in);
	}

	//destructor
	frame::~frame(void)
	{
	}

	//*************************************************************
	//methods

	//initialization method
	void frame::initialize(const double& alpha0_in, const double& alphadot_in, const double& delta0_in, const double& deltadot_in, const double& W_in, const double& Wdot_in)
	{
		alpha0 = alpha0_in;
		alphadot = alphadot_in;
		delta0 = delta0_in;
		deltadot = deltadot_in;
		W = W_in;
		Wdot = Wdot_in;
	}

	//initialization method for Earth J2000 frame, i.e. for when rotation matrices should be the identity matrix
	void frame::initialize()
	{
		alpha0 = -math::PIover2;
		alphadot = 0.0;
		delta0 = math::PIover2;
		deltadot = 0.0;
		W = 3.3187;
		Wdot = 6.3004;
	}

	//construct the rotation matrices between ICRF and the local frame
	void frame::construct_rotation_matrices(const double& ETepoch)
	{
		//first, get the current date in TBB
		double TDBepoch = time_utilities::convert_ET_to_TDB(ETepoch);

		double days_since_reference_epoch = TDBepoch - 2451545.0;
		double centuries_since_reference_epoch = days_since_reference_epoch / 36525;

		//compute the current values of the angles alpha and delta
		double alpha = alpha0 + alphadot * centuries_since_reference_epoch;
		double delta = delta0 + deltadot * centuries_since_reference_epoch;

		//next, construct the rotation matrix
		math::Matrix<double> Rx (3, (math::PIover2 - delta), math::Rxhat);
		math::Matrix<double> Rz (3, (math::PIover2 + alpha), math::Rzhat);
		R_from_local_to_ICRF = Rz*Rx;
		R_from_ICRF_to_local = R_from_local_to_ICRF.transpose();
	}
	
}} //close namespace