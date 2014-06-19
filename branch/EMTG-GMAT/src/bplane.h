//header file for b-plane class
//for use with EMTGv8
//Jacob Englander 1-10-2013


#ifndef _EMTG_BPLANECLASS
#define _EMTG_BPLANECLASS

#include "EMTG_Matrix.h"

namespace EMTG { namespace Astrodynamics
{
	class bplane
	{
	public:
		//constructor
		bplane();

		//destructor
		virtual ~bplane();

		//methods

		//function to define the b-plane coordinate system
		void define_bplane(math::Matrix<double>& V_infinity_in, math::Matrix<double>& R_body, math::Matrix<double>& V_body);

		//function to compute the b-plane Bradius, Btheta, BdotR, BdotT, and periapse distance from an outbound V-infinity vector, for flybys
		void compute_bplane_coordinates_from_Vinfinity_out(const double mu, math::Matrix<double>& V_infinity_out, double* Bradius, double* Btheta, double* BdotR, double* BdotT, double* rp);

		//function to compute the outbound V-infinity vector and periapse distance from Bradius and Btheta
		void compute_Vinfinity_out_from_Bradius_Btheta(const double mu, const double Bradius, const double Btheta, math::Matrix<double>& V_infinity_out, double* rp);

		//function to compute the outbound V-infinity vector and periapse distance from BdotR and BdotT
		void compute_Vinfinity_out_from_BdotR_BdotT(const double mu, const double BdotR, const double BdotT, math::Matrix<double>& V_infinity_out, double* rp);

		//function to compute Bradius and Btheta from periapse position vector
		void compute_Bradius_Btheta_from_periapse_position(const double mu, math::Matrix<double>& Vinfinity_in, math::Matrix<double>& Rp, double* Bradius, double* Btheta);

		//function to compute BdotR and BdotT from periapse position vector
		void compute_BdotR_BdotT_from_periapse_position(const double mu, math::Matrix<double>& Vinfinity_in, math::Matrix<double>& Rp, double* BdotR, double* BdotT);

		//function to convert from Bradius, Btheta to BdotR, BdotT
		void convert_polar_to_cartesian(const double Bradius, const double Btheta, double* BdotR, double* BdotT);

		//function to convert from BdotR, BdotT to Bradius, Btheta
		void convert_cartesian_to_polar(const double BdotR, const double BdotT, double* Bradius, double* Btheta);


	private:
		//fields
		math::Matrix<double> S_hat, h_hat, k_hat, B_hat, T_hat, R_hat, B;
		double Bmag;
		double V_infinity_mag;
		double mu;
		double rbody;

	};

}}//close namespace

#endif //_EMTG_BPLANECLASS