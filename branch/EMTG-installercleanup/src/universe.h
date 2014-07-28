//header file for EMTG universe class

#include <string>
#include <vector>

#include "missionoptions.h"
#include "body.h"
#include "frame.h"

#include "boost/ptr_container/ptr_vector.hpp"

using namespace std;

#ifndef _EMTG_UNIVERSE
#define _EMTG_UNIVERSE

namespace EMTG {namespace Astrodynamics {
	class universe
	{
	public:
		//constructor
		universe(); //default constructor
		universe(const int& j, string universefile, missionoptions* options);

		//destructor
		virtual ~universe();

		//**************************************
		//methods
		
		//function to load new data into the universe
		int load_universe_data(const int& j, string universefile, missionoptions* options);

		//function to find the central body state vector relative to the sun at epoch
		int locate_central_body(const double& epoch, double* state, missionoptions* options);

		//function to create the flyby menu - creates a list of bodies, by SPICE ID, which are flyby capable
		void create_flyby_and_perturbation_menus(const int& j, missionoptions* options);

		//function to print the flyby menu
		void print_flyby_and_perturbation_menus(string filename, missionoptions* options);

		//function to print the universe
		void print_universe(string filename, missionoptions* options);

		//**************************************
		//fields

		//the following fields are read in
		string central_body_name;
		int central_body_SPICE_ID;
		double central_body_radius;
		double mu;
		double LU;
		double r_SOI; //radius of the central body's sphere of influence
		double minimum_safe_distance; //minimum safe distance from the central body, in km
		double reference_epoch; //reference epoch for rotational state


		//the following fields are computed internal to the class
		double TU;
		boost::ptr_vector<body> bodies;
		vector<int> flyby_menu; //vector containing the list of bodies, in the order that they appear in the Universe list, which can be used for flybys
		int size_of_flyby_menu; //this integer will be twice the number of flyby-capable bodies. If it is zero, then no flybys are possible in this universe.
		vector<int> perturbation_menu; //vector containing the list of bodies large enough to be considered for third-body perturbations
		frame LocalFrame; //local reference frame, for use in calculating rotation matrices

	};//end class body

}}//close namespace

#endif //_EMTG_UNIVERSE