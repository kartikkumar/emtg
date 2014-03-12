//source file for b-plane class
//for use with EMTGv8
//Jacob Englander 1-10-2013

#include "EMTG_Matrix.h"
#include "EMTG_math.h"
#include "Astrodynamics.h"
#include "bplane.h"

namespace EMTG { namespace Astrodynamics {
	//default constructor
	bplane::bplane() {}

	//destructor doesn't currently do anything
	bplane::~bplane() {}


	//function to define the b-plane coordinate system. B-plane is defined by the vectors T_hat and R_hat
	//S_hat is the incoming velocity asymptote unit vector
	//from http://ccar.colorado.edu/asen5519/imd/documents/BPlaneHandout.pdf
	void bplane::define_bplane(math::Matrix<double>& V_infinity_in, math::Matrix<double>& R_body, math::Matrix<double>& V_body)
	{
		//compute S_hat
		V_infinity_mag = V_infinity_in.norm();
		S_hat = V_infinity_in.unitize();

		//compute the body orbit normal vector - this is the STK/GMAT "reference vector" from which the b-plane is constructed
		//note: reference frame is ICRF
		k_hat = R_body.unitcross(V_body);

		//T_hat and R_hat
		T_hat = S_hat.unitcross(k_hat);
		R_hat = S_hat.unitcross(T_hat);
	}

	//function to convert from Bradius, Btheta to BdotR, BdotT
	void bplane::convert_polar_to_cartesian(const double Bradius, const double Btheta, double* BdotR, double* BdotT)
	{
		*BdotT = Bradius * cos(Btheta);
		*BdotR = Bradius * sin(Btheta);
	}

	//function to convert from BdotR, BdotT to Bradius, Btheta
	void bplane::convert_cartesian_to_polar(const double BdotR, const double BdotT, double* Bradius, double* Btheta)
	{
		*Bradius = sqrt(BdotR*BdotR + BdotT*BdotT);
		*Btheta = atan2(BdotR, BdotT);
	}

	//function to compute the b-plane Bradius, Btheta, BdotR, BdotT, and periapse distance from an outbound V-infinity vector, for flybys
	//from http://ccar.colorado.edu/asen5519/imd/documents/BPlaneHandout.pdf
	void bplane::compute_bplane_coordinates_from_Vinfinity_out(const double mu, math::Matrix<double>& V_infinity_out, double* Bradius, double* Btheta, double* BdotR, double* BdotT, double* rp)
	{
		//compute the square of the V_infinity magnitude
		double C3 = V_infinity_mag*V_infinity_mag;
		
		//compute h_hat, the unit vector in the direction of spacecraft orbit normal
		//note that both the incoming and outgoing V_infinity will be in the plane of the orbit, so the angular momentum unit vector is the normalized 
		//cross product of the two V_infinity vectors
		h_hat = S_hat.unitcross(V_infinity_out);
		
		//compute unit vector pointing toward b-plane crossing
		B_hat = S_hat.unitcross(h_hat);

		//compute phi
		double phi = acos(S_hat.dot(V_infinity_out) / V_infinity_mag);

		//compute rp
		*rp = mu / C3 * ( 1 / cos( (math::PI - phi) / 2) - 1);

		//compute Bmag and B
		double A = 1 + C3 * (*rp) / mu;

		*Bradius = mu / C3 * sqrt(A*A - 1);
		B = B_hat * *Bradius;

		//compute Btheta
		if (B_hat.dot(R_hat) >= 0)
			*Btheta = acos(T_hat.dot(B_hat));
		else
			*Btheta = 2 * math::PI - acos(T_hat.dot(B_hat));

		//compute BdotR, BdotT
		*BdotR = B.dot(R_hat);
		*BdotT = B.dot(T_hat);
	}

	//function to compute the outbound V-infinity vector and periapse distance from Bradius and Btheta
	void bplane::compute_Vinfinity_out_from_Bradius_Btheta(const double mu, const double Bradius, const double Btheta, math::Matrix<double>& V_infinity_out, double* rp)
	{
		//first convert from (Bradius, Btheta) to (BdotR, BdotT)
		double BdotR, BdotT;
		convert_polar_to_cartesian(Bradius, Btheta, &BdotR, &BdotT);

		compute_Vinfinity_out_from_BdotR_BdotT(mu, BdotR, BdotT, V_infinity_out, rp);
	}

	//function to compute the outbound V-infinity vector and periapse distance from BdotR and BdotT
	void bplane::compute_Vinfinity_out_from_BdotR_BdotT(const double mu, const double BdotR, const double BdotT, math::Matrix<double>& V_infinity_out, double* rp)
	{
		cout << "bplane::compute_Vinfinity_out_from_BdotR_BdotT is not yet implemented!" << endl;
		cin.ignore();

		//basically, compute B_hat and negate the component of V_infinity_in in the direction of Bhat
		//i.e. (V_infinity_in dot B_hat) = (V_infinity_out dot B_hat)
	}

	//function to compute Bradius and Btheta from periapse position vector
	void bplane::compute_Bradius_Btheta_from_periapse_position(const double mu, math::Matrix<double>& Vinfinity_in, math::Matrix<double>& Rp, double* Bradius, double* Btheta)
	{
		//first compute periapse distance
		double rp = Rp.norm();

		//compute C3
		double C3 = V_infinity_mag*V_infinity_mag;
		
		//compute velocity magnitude at periapse
		double vp = sqrt(C3 + 2*mu/rp);

		//compute Bradius
		*Bradius = vp * rp / V_infinity_mag;

		//compute h_hat
		h_hat = Rp.unitcross(Vinfinity_in);

		//compute B
		B_hat = S_hat.unitcross(h_hat);
		B = B_hat * *Bradius;

		//compute Btheta
		if (B_hat.dot(R_hat) >= 0)
			*Btheta = acos(T_hat.dot(B_hat));
		else
			*Btheta = 2 * math::PI - acos(T_hat.dot(B_hat));
	}

	//function to compute BdotR and BdotT from periapse position vector
	void bplane::compute_BdotR_BdotT_from_periapse_position(const double mu, math::Matrix<double>& Vinfinity_in, math::Matrix<double>& Rp, double* BdotR, double* BdotT)
	{
		double Bradius, Btheta;
		compute_Bradius_Btheta_from_periapse_position(mu, Vinfinity_in, Rp, &Bradius, &Btheta);

		bplane::convert_polar_to_cartesian(Bradius, Btheta, BdotR, BdotT);
	}
}}//close namespaces