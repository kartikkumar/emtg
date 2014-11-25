// Dormand-Prince (DOPRI) 8th(7th) order 13 step algorithm
// DOPRI constants are from Numerical Recipies
// Jacob Englander 9/10/2012

#include "missionoptions.h"
#include "universe.h"

#ifndef _RK8713M
#define _RK8713M

namespace EMTG { namespace integration
{

class rk8713M
{
	//fields
	int ns;
	
	std::vector <double> f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, y, x_left, x_right;
		                 
    EMTG::math::Matrix <double> df1dTOF, df2dTOF, df3dTOF, df4dTOF, df5dTOF, df6dTOF, df7dTOF, df8dTOF, df9dTOF, df10dTOF, df11dTOF, df12dTOF, df13dTOF, dydTOF, dx_leftdTOF, dx_rightdTOF;


	//constructor
public:
	rk8713M();
	rk8713M(int ns_in, const int & number_of_phases);

	//destructor
	~rk8713M();

	//methods
	void rk8713M_step(const int & phase_num, 
					  const std::vector <double> & u, 
		              const double & t_left, 
		              std::vector <double> & dt_left_stepdTOF, 
		              const double & t_0, 
		              const double & h, 
		              const double & dhdTOF, 
		              double * error,

					  void(*EOM)(std::vector <double> & x,
					  EMTG::math::Matrix <double> & dx_dTOF,
					  //std::vector <double> & dx_dTOF,
					  const double & t_left_step,
					  std::vector <double> & dt_left_stepdTOF,
					  const double & c,
					  const double & h,
					  const double & dhdTOF,
					  const double & t0,
					  const std::vector <double> & u,
					  std::vector <double> & f,
					  EMTG::math::Matrix <double> & dfdTOF,
					  const int & phase_num,
					  double * thrust,
					  double * mdot,
					  double * Isp,
					  double * power,
					  double * active_power,
					  int * number_of_active_engines,
					  int & STMrows,
					  int & STMcolumns,
					  void * optionsvoidpointer,
					  void * Universepointer,
					  void * ControllerPointer),
					  
					  double* thrust,
					  double* mdot,
					  double* Isp,
					  double* power,
					  double* active_power,
					  int* number_of_active_engines,
					  int &STMrows,
					  int &STMcolumns,
					  void* optionspointer, void* Universepointer, void* ControllerPointer);

	void adaptive_step_int(std::vector <double> & x_left_in, 
						   EMTG::math::Matrix <double> & dx_left_indTOF, 
						   std::vector<double> & x_right_out, 
						   EMTG::math::Matrix <double> & dx_right_outdTOF, 
						   const int & phase_num, 
						   const std::vector <double> & uleft, 
						   const double & t_left_in, 
						   std::vector <double> & dt_left_indTOF, 
						   const double & t_0, 
						   const double & local_step, 
						   double & dlocal_stepdTOF, 
						   double * resumeH, 
						   double * resumeError, 
						   const double & PRECISION_TARGET,
		                                        
						   void(*EOM)(std::vector <double> & x,
						   EMTG::math::Matrix <double> & dx_dTOF,
						   const double & t_left_step,
						   std::vector <double> & dt_left_stepdTOF,
						   const double & c,
						   const double & h,
						   const double & dhdTOF,
						   const double & t0,
						   const std::vector <double> & u,
						   std::vector <double> & f,
						   EMTG::math::Matrix <double> & dfdTOF,
						   const int & phase_num,
						   double * thrust,
						   double * mdot,
						   double * Isp,
						   double * power,
						   double * active_power,
						   int * number_of_active_engines,
						   int & STMrows,
						   int & STMcolumns,
						   void * optionsvoidpointer,
						   void * Universepointer,
						   void * ControllerPointer),
						   
						   double * thrust,
						   double * mdot,
						   double * Isp,
						   double * power,
						   double * active_power,
						   int * number_of_active_engines,
						   int & STMrows,
						   int & STMcolumns,
						   void * optionspointer, void * Universepointer, void * ControllerPointer);
};
	
}} //namespace EMTG::integration

#endif //_RK8713M