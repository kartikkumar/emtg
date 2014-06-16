#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <vector>

#include <string>
#include <sstream>
#include <iostream>
#include <ostream>
#include <istream>
#include <fstream>
#include <iomanip>

//asteroid structure
struct Body{
	double a;
	double ecc;
	double inc;
	double LAN;
	double omega;
	double M;
	double tru;
	double E;
	double h;
	double r;
	double v;
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double time_stamp;
	std::string name;
	int GTOC7_num;
};

//spacecraft structure
struct Spacecraft{
	std::vector <std::vector <double>> a;
	std::vector <std::vector <double>> ecc;
	std::vector <std::vector <double>> inc;
	std::vector <std::vector <double>> LAN;
	std::vector <std::vector <double>> omega;
	std::vector <std::vector <double>> M;
	std::vector <std::vector <double>> tru;
	std::vector <std::vector <double>> E;
	std::vector <std::vector <double>> h;
	std::vector <std::vector <double>> r;
	std::vector <std::vector <double>> v;
	std::vector <std::vector <double>> x;
	std::vector <std::vector <double>> y;
	std::vector <std::vector <double>> z;
	std::vector <std::vector <double>> vx;
	std::vector <std::vector <double>> vy;
	std::vector <std::vector <double>> vz;
	std::vector <std::vector <double>> mass;
	std::vector <std::vector <double>> Tx;
	std::vector <std::vector <double>> Ty;
	std::vector <std::vector <double>> Tz;
	std::vector <std::vector <double>> time_stamp;
	std::string name;
};


const double PI = 3.14159265358979323;

//function prototypes
double mag(std::vector <double> & vector);
double dot(std::vector <double> & vectorA, std::vector <double> & vectorB);
std::vector <double> cross(std::vector <double> & A, std::vector <double> & B);

std::vector <Body> read_ephemeris(std::string & ephemeris_file, const double & mu_sun);
Spacecraft read_probe_summary_file(std::string & probe_summary_file_name, const double & mu_sun);

double laguerreConway(double &, double &);
void coe2cartesian(Body * body, const double & mu_sun);
void cartesian2coe(Body * body, const double & mu_sun);
void cartesian2coe_probe(Spacecraft * probe, const int & phase, const int & timestep, const double & mu_sun);
void coe2cartesian(Spacecraft * probe, const int & phase, const int & timestep, const double & mu_sun);
void orbitprop(Body * body, const double & delta_t, const double & mu_sun);
void orbitprop(Spacecraft * probe, const double & delta_t, const double & phase, const double & timestep, const double & mu_sun);

std::vector <double> rk8713M(std::vector <double> x_left, std::vector <double> Tvec, std::vector <double> f1, double & h, const int & ns, double & error, const double & DU, const double & TU, const double & mu_sun);
std::vector <double> GTOC7EOM(std::vector <double> & X, std::vector <double> Tvec, const double & DU, const double & TU, const double & mu_sun);
std::vector <double> adaptive_step_int(std::vector <double> x_left, std::vector <double> Tvec, double & h, const int & ns, const double & precisionTarget, const double & DU, const double & TU, const double & mu_sun);