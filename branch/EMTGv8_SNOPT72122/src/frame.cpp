//EMTG frame class
//Jacob Englander 1/3/2013

#include <cmath>

#include "frame.h"
#include "EMTG_math.h"


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

	//convert from JED to TDB
	//from http://sourcecodebrowser.com/astronomical-almanac/5.6/tdb_8c.html
	double frame::convert_ET_to_TDB(const double& ETepoch)
	{
		double M, T;

		/* Find time T in Julian centuries from J2000.  */
		T = (ETepoch - 2451545.0)/36525.0;

		/* Mean anomaly of sun = l' (J. Laskar) */
		M = 129596581.038354 * T +  1287104.76154;

		/* Reduce arc seconds mod 360 degrees.  */
		M = M - 1296000.0 * floor( M/1296000.0 );

		M += ((((((((
			  1.62e-20 * T
			- 1.0390e-17 ) * T
			- 3.83508e-15 ) * T
			+ 4.237343e-13 ) * T
			+ 8.8555011e-11 ) * T
			- 4.77258489e-8 ) * T
			- 1.1297037031e-5 ) * T
			+ 1.4732069041e-4 ) * T
			- 0.552891801772 ) * T * T;

		M *= 4.8481368110953599359e-6;
		/* TDB - TDT, in seconds.  */
		T = 0.001658 * sin(M) + 0.000014 * sin(M+M);

		T = ETepoch + T / 86400.0;

		return T;
	}

	//construct the rotation matrices between ICRF and the local frame
	void frame::construct_rotation_matrices(const double& ETepoch)
	{
		//first, get the current date in TBB
		double TDBepoch = convert_ET_to_TDB(ETepoch);

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