//header file for EMTG body class

#include <string>
#include <vector>

#include "missionoptions.h"
#include "frame.h"


using namespace std;

#ifndef _EMTG_BODY
#define _EMTG_BODY

namespace EMTG {namespace Astrodynamics {
	class body
	{
	public:
		//constructor
		body(); //default constructor
		body(const int& ibody_code, 
            const string& iname, 
            const string& ishortname,
            const int& ispice_ID,
            const double& imininum_altitude,
            const double& imass, 
            const double& iradius,
            const double& iepoch,
            vector<double>& ireference_angles, 
            vector<double>& iclassical_orbit_elements,
            const double& iuniverse_mu, 
            const int& icentral_body_SPICE_ID,
            const string& icentral_body_name, 
            const double& icentral_body_radius,
            missionoptions* options);

		//destructor
		virtual ~body();

		//**************************************
		//methods
		
		//function to load new data into the body
		void load_body_data(const int& ibody_code,
            const string& iname, 
            const string& ishortname, 
            const int& ispice_ID, 
            const double& imininum_altitude,
            const double& imass, 
            const double& iradius, 
            const double& iepoch, 
            vector<double>& ireference_angles,
            vector<double>& iclassical_orbit_elements,
            const double& iuniverse_mu, 
            const int& icentral_body_SPICE_ID,
            const string& icentral_body_name, 
            const double& icentral_body_radius,
            missionoptions* options);

		//function to find the body state vector at epoch
		int locate_body(const double& epoch, 
            double* state,
            const bool& need_deriv,
            missionoptions* options) const;

		//function to locate a point on the sphere of influence in cartesian coordinates (Earth Equatorial J2000, measured from central body of current universe)
		int locate_point_on_SOI(const double& theta,
            const double& phi,
            double* point_relative_to_body) const;

		//function to print body to screen, for debug purposes
        void print_body_to_screen(string filename) const;

		//comparisons
        virtual bool operator== (const body& OtherBody) const;
        virtual bool operator!= (const body& OtherBody) const;

		//**************************************
		//fields

		//the following fields are read in
		int spice_ID;
		int body_code;
		int central_body_spice_ID;
		double central_body_radius;
		string central_body_name;
		string name;
		string short_name;
		double mass; //in kg
		double radius; //in km
		double reference_epoch; //in MJD
		double SMA; //in km
		double ECC;
		double INC; //in radians, convert this from input degrees
		double RAAN; //in radians, convert this from input degrees
		double AOP; //in radians, convert this from input degrees
		double MA; //in radians, convert this from input degrees
		double universe_mu; //gravitational constant for the local universe
		double minimum_safe_flyby_altitude; //in km, if this number is <= 0, then the object will not be a member of the flyby menu

		frame J2000_body_equatorial_frame; //local reference frame, for use in calculating rotation matrices

		//the following fields are computed internal to the class
		string SPICE_frame;
		double mu; //gravitational constant
		double r_SOI; //radius of the sphere of influence/hill sphere
	private:
		int body_ephemeris_source; //0: DE-405, 1: SPICE, 2: static
		double ephemeris_start_date;
		double ephemeris_end_date;
		double state_at_beginning_of_ephemeris[6];
		double state_at_end_of_ephemeris[6];
	};//end class body

}}//close namespace

#endif //_EMTG_BODY