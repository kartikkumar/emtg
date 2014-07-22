#include <vector>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include "FBLT_EOM_T.h"
#include "file_utilities.h"
#include "rk8713M_T.h"
#include "SpiceUsr.h"
#include "EMTG_math.h"
#include "EMTG_Matrix.h"

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"


#include "GSAD.h"

void FBLT_chaintest()
{

	//Canonical Unit Normalization?
	double DU;
	double TU;
	double mu_cb;
	bool canonical = true;
	double mu_cb_mks = 132712440018.0;
	double mu_3B;

	if (canonical)
	{
		DU = 149597870.691;
		TU = sqrt(DU*DU*DU / mu_cb_mks);
		mu_cb = 1.0; //normalize the central body's mu
		mu_3B = 126686534.0 / 132712440018.0; //126686534.0
		//mu_3B = 0.0;
	}
	else
	{
		DU = 1.0;
		TU = 1.0;
		mu_cb = mu_cb_mks;
		mu_3B = 126686534.0;
	}

	//STM size
	int rows = 11;
	int columns = 11;

	//masses
	double m_sc = 3000.0; //spacecraft
	double m_3B = 1.8986e+27; //Jupiter
	double m_cb = 1.98855e+30; //Sun

	//create an integrator object
	int ns = 11 + 11 * 11; //spacecraft state vector length + number of entries in the STM 
	int statecount = 11;

	//load all SPICE ephemeris data
	std::vector<fs::path> SPICE_files;
	//F:/research/EMTG_v8/Universe/ephemeris_files/
	//C:/research/Universe/ephemeris_files/
	EMTG::filesystem::get_all_files_with_extension(fs::path("F:/research/EMTG_v8/Universe/ephemeris_files/"), ".bsp", SPICE_files);
	EMTG::filesystem::get_all_files_with_extension(fs::path("F:/research/EMTG_v8/Universe/ephemeris_files/"), ".cmt", SPICE_files);

	std::string filestring;
	for (size_t k = 0; k < SPICE_files.size(); ++k)
	{
		filestring = "F:/research/EMTG_v8/Universe/ephemeris_files/" + SPICE_files[k].string();
		furnsh_c(filestring.c_str());
	}

	//SPICE reference frame kernel
	std::string leapsecondstring = "F:/research/EMTG_v8/Universe/ephemeris_files/naif0009.tls";
	std::string referenceframestring = "F:/research/EMTG_v8/Universe/ephemeris_files/pck00010.tpc";
	furnsh_c(leapsecondstring.c_str());
	furnsh_c(referenceframestring.c_str());

	//disable SPICE errors. This is because we can, and will often, go off the edge of an ephemeris file.
	errprt_c("SET", 100, "NONE");
	erract_c("SET", 100, "RETURN");


	//declare adouble containers for the beginning and end spacecraft states

#ifdef _AD_VERIFICATION
	std::vector <GSAD::adouble> state0(ns,0.0); //initial state
	std::vector <GSAD::adouble> state1(ns, 0.0); //state at the end of the half phase
	std::vector <GSAD::adouble> dstate1dTOF(ns, 0.0); //final state partials w.r.t. TOF
	std::vector <GSAD::adouble> state_current(ns, 0.0);
	std::vector <GSAD::adouble> dstate_currentdTOF(ns, 0.0); //partials of the left hand state vector w.r.t. TOF
	std::vector <GSAD::adouble> state_next(ns, 0.0);
	std::vector <GSAD::adouble> dstate_nextdTOF(ns, 0.0); //partials of the right hand state vector w.r.t. TOF
	
	double launch_epoch = (56774.0 * 86400.0) / TU;
	GSAD::adouble t_left = (56774.0 * 86400.0) / TU; //April 27th 2014
	GSAD::adouble dt_leftdTOF = 0.0;
#else
	std::vector <double> state0(ns,0.0);
	std::vector <double> state1(ns,0.0);
	std::vector <double> dstate1dTOF(ns, 0.0);
	std::vector <double> state_current(ns,0.0);
	std::vector <double> dstate_currentdTOF(ns, 0.0);
	std::vector <double> state_next(ns,0.0);
	std::vector <double> dstate_nextdTOF(ns, 0.0);

	double launch_epoch = (56774.0 * 86400.0) / TU;
	double t_left = (56774.0 * 86400.0) / TU;
	double dt_leftdTOF = 0.0;
#endif

	//Initial spacecraft state: same as Jupiter but 95% sma and 1 degree lag in true anomaly
	state0[0] = -2.956394096440888e+08 / DU; //x
	state0[1] = 6.242131720747155e+08 / DU;  //y
	state0[2] = 2.747527803730505e+08 / DU;  //z
	state0[3] = -12.471875356709877 / DU * TU; //vx
	state0[4] = -4.432681177438947 / DU * TU;  //vy
	state0[5] = -1.596329404751935 / DU * TU;  //vz
	state0[6] = m_sc; //mass
	state0[7] = 0.1;  //ux
	state0[8] = 0.1;  //uy
	state0[9] = 0.1;  //uz
	state0[10] = (200.0 * 86400.0) / TU; // time of flight (TOF) -- also total propagation time
	
	//STM is initialized to the identity matrix
	for (size_t i = statecount; i < ns; ++i)
		state0[i] = 0.0;

	for (size_t i = statecount; i < ns; i = i + statecount + 1)
		state0[i] = 1.0;

	//set derivative of each state variable with respect to itself equal to 1 (to seed the auto-diff routine)	
#ifdef _AD_VERIFICATION
	for (size_t i = 0; i < statecount; ++i)
		state0[i].setDerivative(i, 1.0);
#endif

	//**************************
	//
	//Numerically integrated STM
	//
	//**************************

	//instantiate an integrator
#ifdef _AD_VERIFICATION
	EMTG::integration::rk8713M<GSAD::adouble> *integrator;
	integrator = new EMTG::integration::rk8713M<GSAD::adouble>(ns);
#else
	EMTG::integration::rk8713M<double> *integrator;
	integrator = new EMTG::integration::rk8713M<double>(ns);
#endif

	//Propagate spacecraft for all time steps
	int numsteps = 10;

	//some integrator settings
#ifdef _AD_VERIFICATION
	GSAD::adouble resumeH = state0[10] / numsteps; //take a guess at a starting integration step size, try to do the whole integration in one step
#else
	double resumeH = state0[10] / numsteps;
#endif

	double resumeError = 1.0e-13;
	double PRECISION_TARGET = 1.0e-11; //integration error tolerance

	//create a vector of STM's.....one STM per time step
#ifdef _AD_VERIFICATION
	std::vector < EMTG::math::Matrix<GSAD::adouble> > Phi_archive(numsteps, EMTG::math::Matrix<GSAD::adouble> (rows, columns, 0.0));
#else
	std::vector < EMTG::math::Matrix<double> > Phi_archive(numsteps, EMTG::math::Matrix< double >(rows, columns, 0.0));
#endif

	std::cout << "INTEGRATING TRAJECTORY, BUILDING STM: " << std::endl << std::endl;
	clock_t tStart = clock();

	state_current = state0;

	//setup some initial TOF derivatives
#ifdef _AD_VERIFICATION
	GSAD::adouble steptime = state0[10] / numsteps;
	GSAD::adouble dsteptimedTOF = 1.0 / numsteps;
#else
	double steptime = state0[10] / numsteps;
	double dsteptimedTOF = 1.0 / numsteps;
#endif


	//FBLT propagation loop

	for (size_t step = 0; step < numsteps; ++step)
	{
		state_next = integrator->adaptive_step_int(state_current, dstate_currentdTOF, dstate_nextdTOF, t_left, dt_leftdTOF, launch_epoch, mu_3B, mu_cb, steptime, dsteptimedTOF, &resumeH, &resumeError, PRECISION_TARGET, DU, TU, EMTG::Astrodynamics::EOM::FBLT_EOM);
		//Build and store numerically integrated STM for this FBLT time step
		statecount = 11;
		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < columns; ++j)
			{
				Phi_archive[step](i, j) = state_next[statecount];
				++statecount;
			}
		}

#ifdef _AD_VERIFICATION
		//AD calculated derivatives
		//double dvxdy0 = state_next[1].getDerivative(2);
		//std::cout << dvxdy0 << std::endl;

		std::cout << "STM for time step " << step << std::endl;
		std::cout << setprecision(16);
		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < columns; ++j)
				std::cout << Phi_archive[step](i, j) << "     ";

			std::cout << std::endl << std::endl;
		}

		std::cout << std::endl << std::endl;
#endif	

		//state at the end of the time step becomes the new current state
		state_current = state_next;
		
		//state TOF partial derivatives at the right of the current step become the partials for the left of the next 
		dstate_currentdTOF = dstate_nextdTOF;

		//STM for the next step is initialized to the identity matrix
		for (size_t i = 11; i < ns; ++i)
			state_current[i] = 0.0;

		for (size_t i = 11; i < ns; i = i + 12)
			state_current[i] = 1.0;

		//change control ?
		//state_current[7] += 0.0005;
		//state_current[8] -= 0.0005;
		//state_current[9] -= 0.0005 ;

		t_left += steptime; //advance the mission clock
		dt_leftdTOF += dsteptimedTOF;

	} //end FBLT propagation loop


	double elapsed_time = double((clock() - tStart)) / CLOCKS_PER_SEC;
	std::cout.precision(16);
	std::cout << "Trajectory propagated and STM constructed in: " << std::to_string(elapsed_time) << " s" << std::endl << std::endl;


	state1 = state_next;
	dstate1dTOF = dstate_nextdTOF;

	//Build match-point derivatives
#ifdef _AD_VERIFICATION
	EMTG::math::Matrix <GSAD::adouble> derivatives(rows, columns, 0.0);
#else
	EMTG::math::Matrix <double> derivatives(rows, columns, 0.0);
#endif
	//initialize the derivatives matrix to the identity
	for (size_t i = 0; i < rows; ++i)
		derivatives(i, i) = 1.0;
	//build the derivatives matrix through successive STM multiplication
	for (int step = numsteps-1; step >= 0; --step)
		derivatives *= Phi_archive[step];
	


	//Insert the TOF partials into the final numerically integrated STM so we have everything in one place
	for (size_t i = 0; i < rows; ++i)
		derivatives(i, 10) = dstate1dTOF[i];

	std::cout << "Numerically integrated resultant STM:" << std::endl;
	std::cout << setprecision(16);
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
			std::cout << derivatives(i,j) << "     ";

		std::cout << std::endl << std::endl;
	}

	std::cout << std::endl << std::endl;

#ifdef _AD_VERIFICATION
	//AD calculated derivatives
	EMTG::math::Matrix<GSAD::adouble> PhiAD(rows, columns, 0.0);
	EMTG::math::Matrix<GSAD::adouble> PhiERROR(rows, columns, 0.0);

	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < columns; ++j)
			PhiAD(i, j) = state1[i].getDerivative(j);

	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < columns; ++j)
			PhiERROR(i, j) = PhiAD(i, j) - derivatives(i, j);

	std::cout << "AD calculated resultant STM:" << std::endl;
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
		{
			std::cout << PhiAD(i, j) << "     ";
		}

		std::cout << std::endl << std::endl;
	}
	
	std::cout << std::endl;

	std::cout << "Error between AD STM and integrated STM:" << std::endl;
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
		{
			std::cout << PhiERROR(i, j) << "     ";
		}

		std::cout << std::endl << std::endl;
	}

	std::cout << std::endl;
#endif
	
	getchar();

}