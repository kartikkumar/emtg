//Lambert solver by Arora and Russell
//A FAST AND ROBUST MULTIPLE REVOLUTION LAMBERT ALGORITHM USING A COSINE TRANSFORMATION
//AAS Hilton Head 2013
//by Matthew Vavrina, matthew.vavrina@ai-solutions.com and Jacob Englander (6-11-2014)

#include <cmath>
#include <vector>
#include <cmath> 
#include <complex>

#include <string>
#include <iostream>
#include <fstream>

#include "Lambert_AroraRussell.h"
#include "EMTG_math.h"



using namespace std;


namespace EMTG {
    namespace Astrodynamics {
        namespace Lambert {
            //assume input and output units are consistent (i.e. designed for km and s)
            void Lambert_AroraRussell(const double* R1,
                                    const double* R2,
                                    const double& TOFin,
                                    const double& mu,
                                    const int& Nrev,
                                    const bool& LongWay,
                                    const bool& ShortPeriod,
                                    const double& tolerance,
                                    const int& max_iterations,
                                    double* V1,
                                    double* V2,
                                    double& error,
                                    int& iterations)
            {
                //Step 0: declare some variables
                double k; //iteration variable
                double deltak = 1.0e+10; //error in current solution for k
                double r1 = sqrt(R1[0] * R1[0] + R1[1] * R1[1] + R1[2] * R1[2]); //magnitude of initial position
                double r2 = sqrt(R2[0] * R2[0] + R2[1] * R2[1] + R2[2] * R2[2]); //magnitude of initial position
                double LU = r1; // length unit based on magnitude of r1
                double TU = sqrt((1 / mu)*LU*LU*LU); // time unit set so that mu = 1;

                // normalize r1, r2, and mu
                double R1_n[3] = { R1[0] / LU, R1[1] / LU, R1[2] / LU };
                double R2_n[3] = { R2[0] / LU, R2[1] / LU, R2[2] / LU };
                double r1_n = 1.0;
                double r2_n = r2 / r1;
                double TOF = TOFin / TU;
                double mu_n = mu*(TU*TU) / (LU*LU*LU);

                // define transfer angle based on LongWay flag
                double ctheta = math::dot(R1_n, R2_n, 3) / (r1_n * r2_n); //cosine of the transfer angle
                double theta = acos(ctheta); //transfer angle
                if (LongWay == 1)
                {
                    theta = 2 * math::PI - theta;
                }
                else
                {
                    theta = theta;
                }

                // calculate S and tau
                double S = sqrt((r1_n + r2_n) * (r1_n + r2_n) * (r1_n + r2_n) / mu_n); //semi-perimeter
                iterations = 0; //iterations counter
                double sq2 = sqrt(2.0);
                double eps = 2.0e-2; //to prevent singularity at k = sqrt(2)
                double d = theta <= math::PI ? 1.0 : -1.0;
                double tau = d * sqrt(r1_n*r2_n*(1 + ctheta)) / (r1_n + r2_n); //lambert geometry parameter
                double r_buff = 0.2; //user-defined parameter to determine when to skip k_bi root solve

                //Step 1: generate appropriate initial guess
                //declare some variables that will be used in the initial guess

                //Step 1.1 compare the desired time of flight to the parabolic time of flight
                double T_parabolic = S * sqrt(1 - sq2*tau) * (tau + sq2) / 3.0;

                if (TOF <= T_parabolic) //a hyperbolic trajectory is desired
                {
                    double k_n, k_m, k_i, Z, alpha, F_0, F_1, F_i, F_star;
                    double TOF20 = S * sqrt(1.0 - 20.0 * tau) * (tau + 0.04940968903 * (1.0 - 20.0 * tau));
                    double TOF100 = S * sqrt(1.0 - 100.0 * tau) * (tau + 0.00999209404 * (1.0 - 100.0 * tau));
                    if (d > 0)
                    {
                        k_n = sq2;
                        k_m = 1.0 / tau;
                        k_i = (k_n + k_m) / 2.0;
                        Z = 1.0 / sq2;
                        alpha = 0.5;
                        F_0 = T_parabolic;
                        F_1 = 0.0;
                        double m_i = 2 - k_i * k_i;
                        double W = compute_W(k_i, m_i, Nrev);
                        F_i = S * sqrt(1 - k_i*tau) * (tau + (1 - k_i*tau) * W);
                        F_star = TOF;

                        double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                        k = k_n + (k_m - k_n) * x_star;
                        //cout << "H, d>0; x_star, k: " << x_star << ", " << k << "\n"; //remove
                    }
                    else if (TOF > TOF20) //Arora's "H1" region
                    {
                        //cout << "H1 \n"; //remove
                        k_n = sq2;
                        k_m = 20.0;
                        k_i = (2 * k_n + k_m) / 3.0;
                        Z = 1.0 / 3.0;
                        alpha = 1.0;
                        F_0 = T_parabolic;
                        F_1 = TOF20;
                        double m_i = 2 - k_i * k_i;
                        double W = compute_W(k_i, m_i, Nrev);
                        F_i = S * sqrt(1 - k_i*tau) * (tau + (1 - k_i*tau) * W);
                        F_star = TOF;

                        double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                        k = k_n + (k_m - k_n) * x_star;
                    }
                    else //Arora's "H2" region
                    {
                        //cout << "H2 \n"; //remove
                        double T0 = TOF20;
                        double T1 = TOF100;
                        double k_check = pow((TOF100*(TOF20 - TOF) * 10 - TOF20*sqrt(20)*(TOF100 - TOF)) / TOF*(TOF20 - TOF100), 2); //remove
                        k = ((T1*(T0 - TOF)*10.0 - T0*sqrt(20)*(T1 - TOF)) / (TOF*(T0 - T1))) * ((T1*(T0 - TOF)*10.0 - T0*sqrt(20)*(T1 - TOF)) / (TOF*(T0 - T1)));
                    }
                }
                else if (Nrev >= 1.0) //a multi-revolution elliptical orbit is desired
                {
                    // precomputed values for for delE_bi0 for first 20 revs (deltaE_b point where tau crosses zero from Arora eqn 55)
                    double delE_bi0[20] = { 2.848574, 2.969742, 3.019580, 3.046927, 3.064234, 3.076182, 3.084929, 3.091610, 3.096880, 3.101145,
                        3.104666, 3.107623, 3.110142, 3.112312, 3.114203, 3.115864, 3.117335, 3.118646, 3.119824, 3.120886 };

                    // calculate estimate for k_bi (k_biGuess)
                    double sgn_tau = tau >= 0 ? 1.0 : -1.0;
                    double var2 = delE_bi0[Nrev - 1]; // dummy variable 2
                    double var1 = 8.0*fabs(tau) / (var2*(sq2 - 2.0*fabs(tau))); // dummy variable 1
                    double delE_biGuess = var2*(1 - sgn_tau) + var2*sgn_tau*pow((1.0 / (1.0 + var1)), 0.25);
                    double var3 = EMTG::math::PI - delE_biGuess;
                    double sgn_var3 = var3 >= 0 ? 1.0 : -1.0;
                    double k_biGuess = sgn_var3*sqrt(cos(delE_biGuess) + 1.0);

                    // calculate T_biGuess
                    double m_kbiGuess = 2.0 - k_biGuess*k_biGuess; // m based on k_biGuess
                    double W_kbiGuess = compute_W(k_biGuess, m_kbiGuess, Nrev);
                    double T_biGuess = S*sqrt(1 - k_biGuess*tau)*(tau + (1.0 - k_biGuess*tau)*W_kbiGuess);

                    // root solve to find k_bi and T_bi
                    double k_bi;
                    double T_bi;
                    if ((fabs(TOF - T_biGuess)) > (r_buff*TOF) && (TOF > T_biGuess)) //do not need to root solve
                    {
                        k_bi = k_biGuess;
                        T_bi = T_biGuess;
                    }
                    else // find k_bi using Newton Raphson
                    {
                        k_bi = compute_kb(k_biGuess, tau, S, Nrev, tolerance, max_iterations, sq2, eps);  // solve via Newton Raphson
                        double m_bi = 2.0 - k_bi*k_bi;
                        double W_bi = compute_W(k_bi, m_bi, Nrev);
                        T_bi = S*sqrt(1 - k_bi*tau)*(tau + (1.0 - k_bi*tau)*W_bi);
                    }

                    if (TOF < T_bi)
                    {
                        //return - no solution for this Nrev
                        //cout << "No solution for this Nrev \n"; //remove
                        throw 200000;
                        return;
                    }

                    double W_km1 = 5.71238898 + 2 * EMTG::math::PI*Nrev;
                    double W_kmp5 = 1.95494660 + 2.71408094*Nrev;
                    double W_k0 = sq2 / 4.0*(EMTG::math::PI + 2 * EMTG::math::PI*Nrev);
                    double W_kp5 = 0.75913433 + 2.71408094*Nrev;
                    double W_k1 = 0.57079632 + 2 * EMTG::math::PI*Nrev;

                    double TOF_km1 = compute_TOF(-1.0, S, tau, W_km1);  // (k,S,tau,W)
                    double TOF_kmp5 = compute_TOF(-0.5, S, tau, W_kmp5);  // (k,S,tau,W)
                    double TOF_k0 = compute_TOF(0.0, S, tau, W_k0);  // (k,S,tau,W)
                    double TOF_kp5 = compute_TOF(0.5, S, tau, W_kp5);  // (k,S,tau,W)
                    double TOF_k1 = compute_TOF(1.0, S, tau, W_k1);  // (k,S,tau,W)

                    // generate initial guess for k (k_star) for long period and short period solutions for a given Nrev
                    if (k_bi >= 0.0)
                    {
                        if (ShortPeriod == 0) // long period solution
                        {
                            if ((TOF >= TOF_k1) && (k_bi >= 1.0))// use first row of table 5
                            {
                                // compute intial guess for k
                                double k_n = k_bi;
                                double k_m = sq2;
                                double k_i = (k_bi + sq2)*0.5;
                                double Z = 0.25;
                                double alpha = 2.0;
                                double F_0 = 1.0 / T_bi;
                                double F_1 = 0.0;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = 1.0 / TOF_ki;
                                double F_star = 1.0 / TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 1, table 5 \n";//remove
                            }
                            else if ((TOF >= TOF_k1) && (k_bi <= 1.0)) // use second row of table 5
                            {
                                // compute intial guess for k
                                double k_n = 1.0;
                                double k_m = sq2;
                                double k_i = (1.0 + 2.0*sq2) / 3.0;
                                double Z = 4.0 / 9.0;
                                double alpha = 2.0;
                                double F_0 = 1.0 / TOF_k1;
                                double F_1 = 0.0;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = 1.0 / TOF_ki;
                                double F_star = 1.0 / TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 2, table 5 \n";//remove
                            }
                            else //((TOF < TOF_k1) && (k_bi <= 1.0)) // use third row of table 5
                            {
                                // compute intial guess for k
                                double k_n = k_bi;
                                double k_m = 1.0;
                                double k_i = (1.0 + k_bi)*0.5;
                                double Z = 0.25;
                                double alpha = 2.0;
                                double F_0 = T_bi;
                                double F_1 = TOF_k1;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = TOF_ki;
                                double F_star = TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 3, table 5 \n";//remove
                            }
                        }
                        else // long period solution
                        {
                            if (TOF < TOF_k0) // use fourth row of table 5
                            {
                                // compute intial guess for k
                                double k_n = 0.0;
                                double k_m = k_bi;
                                double k_i = k_bi*0.5;
                                double Z = pow(0.5, 6.0 / 5.0);
                                double alpha = 6.0 / 5.0;
                                double F_0 = TOF_k0;
                                double F_1 = T_bi;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = TOF_ki;
                                double F_star = TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 4, table 5 \n";//remove
                            }
                            else if ((TOF > TOF_k0) && (TOF < TOF_km1)) // use fifth row of table 5
                            {
                                // compute intial guess for k
                                double k_n = -1.0;
                                double k_m = 0.0;
                                double k_i = -0.5;
                                double Z = 0.5;
                                double alpha = 1.0;
                                double F_0 = TOF_km1;
                                double F_1 = TOF_k0;
                                double F_i = TOF_kmp5;
                                double F_star = TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 5, table 5 \n";//remove
                            }
                            else //(TOF > TOF_km1) // use sixth row of table 5
                            {
                                // compute intial guess for k
                                double k_n = -1.0;
                                double k_m = -1.0*sq2;
                                double k_i = (-1.0 - 2 * sq2) / 3.0;
                                double Z = 4.0 / 9.0;
                                double alpha = 2.0;
                                double F_0 = 1.0 / TOF_km1;
                                double F_1 = 0.0;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = 1.0 / TOF_ki;
                                double F_star = 1.0 / TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 6, table 5 \n";//remove
                            }
                        }
                    }

                    else // k_bi < 0 
                    {
                        if (ShortPeriod == 1) // short period solution
                        {
                            if ((TOF >= TOF_km1) && (k_bi <= -1.0))// use first row of table 6
                            {
                                // compute intial guess for k
                                double k_n = k_bi;
                                double k_m = -1.0*sq2;
                                double k_i = (k_bi - sq2)*0.5;
                                double Z = 0.25;
                                double alpha = 2.0;
                                double F_0 = 1.0 / T_bi;
                                double F_1 = 0.0;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = 1.0 / TOF_ki;
                                double F_star = 1.0 / TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 1, table 6 \n"; //remove
                            }
                            else if ((TOF >= TOF_km1) && (k_bi >= -1.0))// use second row of table 6
                            {
                                // compute intial guess for k
                                double k_n = -1.0;
                                double k_m = -sq2;
                                double k_i = (-1.0 - 2.0*sq2) / 3.0;
                                double Z = 4.0 / 9.0;
                                double alpha = 2.0;
                                double F_0 = 1.0 / TOF_km1;
                                double F_1 = 0.0;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = 1.0 / TOF_ki;
                                double F_star = 1.0 / TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 2, table 6 \n"; //remove
                            }
                            else // ((TOF <= TOF_km1) && (k_bi >= -1.0)) // use third row of table 6
                            {
                                // compute intial guess for k
                                double k_n = k_bi;
                                double k_m = -1.0;
                                double k_i = (-1.0 + k_bi)*0.5;
                                double Z = 0.25;
                                double alpha = 2.0;
                                double F_0 = T_bi;
                                double F_1 = TOF_km1;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = TOF_ki;
                                double F_star = TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 3, table 6 \n"; //remove
                            }
                        }
                        else // long period solution
                        {
                            if (TOF <= TOF_k0) // use fourth row of table 6
                            {
                                // compute intial guess for k
                                double k_n = k_bi;
                                double k_m = 0.0;
                                double k_i = k_bi*0.5;
                                double Z = 0.25;
                                double alpha = 2.0;
                                double F_0 = T_bi;
                                double F_1 = TOF_k0;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = TOF_ki;
                                double F_star = TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 4, table 6 \n"; //remove
                            }
                            else if ((TOF > TOF_k0) && (TOF < TOF_k1)) // use fifth row of table 6
                            {
                                // compute intial guess for k
                                double k_n = 0.0;
                                double k_m = 1.0;
                                double k_i = 0.5;
                                double Z = pow(0.5, 6.0 / 5.0);
                                double alpha = 6.0 / 5.0;
                                double F_0 = TOF_k0;
                                double F_1 = TOF_k1;
                                double F_i = TOF_kp5;
                                double F_star = TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 5, table 6 \n"; //remove
                            }
                            else //(TOF > TOF_k1) // use sixth row of table 6
                            {
                                // compute intial guess for k
                                double k_n = 1.0;
                                double k_m = sq2;
                                double k_i = (1.0 + 2.0*sq2) / 3.0;
                                double Z = 4.0 / 9.0;
                                double alpha = 2.0;
                                double F_0 = 1.0 / TOF_k1;
                                double F_1 = 0.0;
                                double m_ki = 2.0 - k_i*k_i;
                                double W_ki = compute_W(k_i, m_ki, Nrev);
                                double TOF_ki = compute_TOF(k_i, S, tau, W_ki);
                                double F_i = 1.0 / TOF_ki;
                                double F_star = 1.0 / TOF;
                                double x_star = pow((Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star)), 1.0 / alpha);
                                k = k_n + (k_m - k_n) * x_star;
                                //cout << "row 6, table 6 \n"; //remove
                            }
                        }
                    }
                }
                else if (Nrev < 1.0) // a single-revolution elliptical orbit is desired
                {
                    //determine elliptical region by comparing actual TOF (i.e., Tstar) with TOF(k) function value
                    double best_k_test = 9e9;
                    int k_ind = 0;
                    double W_km1 = 5.712388981; // W calculated at k=-1
                    double W_km1p41 = 4839.684497246; // W calculated at k=-1.41
                    double W_km1p38 = 212.087279879; // W calculated at k=-1.38
                    double m_km1p38 = 2.0 - (1.38*1.38); //remove
                    double W_kmp5 = 1.954946607; // W calculated at k=-0.5
                    double W_k1dsq2 = 0.6686397730; // W calculate at k=a/sqrt(2)
                    double TOF_km1p41 = S*sqrt(1.0 + 1.41*tau)*(tau + (1 + 1.41*tau)*W_km1p41);
                    double TOF_km1p38 = S*sqrt(1.0 + 1.38*tau)*(tau + (1 + 1.38*tau)*W_km1p38);
                    double TOF_km1 = S*sqrt(1.0 + tau)*(tau + (1 + tau)*W_km1);
                    double TOF_k0 = S*(sq2 / 4 * EMTG::math::PI + tau);
                    double TOF_k1dsq2 = S*sqrt(1.0 - 1.0 / sq2*tau)*(tau + (1.0 - 1.0 / sq2*tau)*W_k1dsq2);
                    double TOF_kmp5 = S*sqrt(1.0 + 0.5*tau)*(tau + (1 + 0.5*tau)*W_kmp5);
                    double k_n, k_m, k_i, alpha, F_1, F_i, F_star;

                    if (TOF >= TOF_km1p38) // use region E4 for guess
                    {
                        //cout << "E4 \n"; //remove
                        k_n = -1.38;
                        k_m = -sq2;
                        k_i = -1.41;
                        double c_1 = 49267 / 27059;
                        double c_2 = 67286 / 17897;
                        double c_3 = 2813 / 287443;
                        double c_4 = 4439 / 3156;
                        alpha = 243;
                        double F_n = 1 / TOF_km1p38;
                        F_i = 1 / TOF_km1p41;
                        F_star = 1 / TOF;
                        double gamma_1 = F_i*(F_star - F_n);
                        double gamma_2 = F_star*(F_n - F_i);
                        double gamma_3 = F_n*(F_star - F_i);
                        k = -c_4*pow(((gamma_1*c_1 - c_3*gamma_3)*c_2 + c_3*c_1*gamma_2) / (gamma_3*c_1 - gamma_1*c_3 - gamma_2*c_2), 1 / alpha);
                    }
                    else if ((TOF_km1p38 >= TOF) && (TOF >= TOF_km1)) // use region E3 for guess
                    {
                        //cout << "E3 \n"; //remove
                        k_n = -1.0;
                        k_m = -1.0*sq2;
                        k_i = -1.38;
                        double c_1 = 540649.0 / 3125.0;
                        double c_2 = 256.0;
                        double c_3 = 1.0;
                        double c_4 = 1.0;
                        double alpha = 16.0;
                        double F_n = 1 / TOF_km1;
                        F_i = 1 / TOF_km1p38;
                        F_star = 1 / TOF;
                        double gamma_1 = F_i*(F_star - F_n);
                        double gamma_2 = F_star*(F_n - F_i);
                        double gamma_3 = F_n*(F_star - F_i);
                        k = -c_4*pow(((gamma_1*c_1 - c_3*gamma_3)*c_2 + c_3*c_1*gamma_2) / (gamma_3*c_1 - gamma_1*c_3 - gamma_2*c_2), (1.0 / alpha));
                    }
                    else if ((TOF_km1 >= TOF) && (TOF >= TOF_k0)) // use region E2 for guess
                    {
                        //cout << "E2 \n"; //remove
                        k_n = 0.0;
                        k_m = -1.0;
                        k_i = -0.5;
                        double Z = 0.5;
                        alpha = 1.0;
                        double F_0 = TOF_k0;
                        double F_1 = TOF_km1;
                        double F_i = TOF_kmp5;
                        double F_star = TOF;
                        double x_star = (Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star));
                        k = k_n + (k_m - k_n) * x_star;
                    }
                    else // use region E1 for guess
                    {
                        //cout << "E1 \n"; //remove
                        k_n = 0.0;
                        k_m = sq2;
                        k_i = 1 / sq2;
                        double Z = 0.5;
                        alpha = 1.0;
                        double F_0 = TOF_k0;
                        double F_1 = T_parabolic;
                        double F_i = TOF_k1dsq2;
                        double F_star = TOF;
                        double x_star = (Z * (F_0 - F_star)*(F_1 - F_i)) / ((F_i - F_star)*(F_1 - F_0)*Z + (F_0 - F_i)*(F_1 - F_star));
                        k = k_n + (k_m - k_n) * x_star;
                    }
                }

                //Step 2: iterate to find k
                while (fabs(deltak) > tolerance && iterations < max_iterations)
                {
                    //Step 2.1 increment the iterations counter
                    ++iterations;

                    //Step 2.2 compute W, dW, ddW
                    double m = 2 - k*k;
                    double sgnk = k >= 0 ? 1.0 : -1.0;
                    double W, dW, ddW;
                    W = compute_W(k, m, Nrev);

                    //dW and ddW
                    if ((k < sq2 - eps) || (k < sq2 && Nrev > 0))
                    {
                        dW = (-2.0 + 3.0 * W * k) / m;
                        ddW = (5.0 * dW * k + 3.0 * W) / m;
                    }
                    else if (k > sq2 + eps)
                    {
                        dW = (-2.0 + 3.0 * W * k) / m;
                        ddW = (5.0 * dW * k + 3.0 * W) / m;
                    }
                    else
                    {
                        double v = k - sq2;
                        double v2 = v*v;
                        double v3 = v*v2;
                        double v4 = v3*v;
                        double v5 = v4*v;
                        double v6 = v5*v;
                        double v7 = v6*v;
                        dW = -1 / 5.0
                            + sq2 * v * (4.0 / 35.0)
                            - v2 * (6.0 / 63.0)
                            + sq2 * v3 * (8.0 / 231.0)
                            - v4 * (10.0 / 429.0)
                            + sq2 * v5 * (48.0 / 6435.0)
                            - v6 * (56.0 / 12155.0)
                            + sq2 * v7 * (64.0 / 46189.0);
                        ddW = sq2 * (4.0 / 35.0)
                            - v * (12.0 / 63.0)
                            + sq2 * v2 * (24.0 / 231.0)
                            - v3 * (40.0 / 429.0)
                            + sq2 * v4 * (240.0 / 6435.0)
                            - v5 * (336.0 / 12155.0)
                            + sq2 * v6 * (448.0 / 46189.0);
                    }


                    //Step 2.3 compute TOFc, dTOFc, ddTOFc
                    double TOFc = S * sqrt(1 - k*tau) * (tau + (1 - k*tau) * W);
                    if ((fabs(TOF - TOFc) < tolerance))
                    {
                        break;
                    }

                    double c = (1 - k * tau) / tau;
                    double sqrtctau = sqrt(1 - k*tau);
                    double dTOFc = -TOFc / (2.0 * c) + S * tau * sqrtctau * (dW * c - W);
                    double ddTOFc = -TOFc / (4.0 * c*c) + S * tau * sqrtctau * (W / c + c*ddW - 3.0*dW);

                    //Step 2.4 compute deltak
                    deltak = -(TOFc - TOF) / (dTOFc - (TOFc - TOF) * ddTOFc / (2.0 * dTOFc));

                    //Step 2.5 update k from deltak
                    k += deltak;

                    // Step 2.6 bound k 
                    // ensure k is not less than -sqrt(2)
                    if (k < -sq2)
                    {
                        k = -sq2 + 1e-14;
                    }
                    // ensure k doesn't bleed into to elliptic area when hyperbolic
                    if ((TOF < T_parabolic) && (k < sq2))
                    {
                        k = sq2 + 1e-14;
                    }
                    // ensure k doesn't bleed into to hyperbolic area when elliptic
                    if ((TOF > T_parabolic) && (k > sq2))
                    {
                        k = sq2 - 1e-14;
                    }
                    // ensure TOF doesn't become indeterminate when d=1
                    if ((TOF < T_parabolic) && (d > 0) && ((1.0 - tau*k) < 0.0))
                    {
                        k = 1.0 / tau - 1e-14;
                    }
                }
                double m_k = 2.0 - k*k;
                double W_k = compute_W(k, m_k, Nrev);
                error = TOF - compute_TOF(k, S, tau, W_k);


                // Step 3: calculate f & g (typos in Arora-Russell 2013 AAS paper)
                double sgnk = k >= 0 ? 1.0 : -1.0;
                std::complex <double> q_comp, x_comp, p_comp, a_comp;
                if (k <= sq2)
                {
                    q_comp = (1 - sgnk)*EMTG::math::PI + sgnk*acos(k*k - 1.0) + 2.0*EMTG::math::PI*Nrev; //ellipse only
                }
                else
                {
                    q_comp = (0, -acoshAR(k*k - 1));
                }
                //double q = real(q_comp);
                double q = real(q_comp);
                std::complex <double> r1_n_c(r1_n, 0.0);
                std::complex <double> r2_n_c(r2_n, 0.0);
                std::complex <double> tau_c(tau, 0.0);
                std::complex <double> k_c(k, 0.0);
                x_comp = q_comp*sqrt((r1_n_c + r2_n_c)*(1.0 - k_c*tau_c) / (2.0 - k_c*k_c));
                p_comp = r1_n*r2_n*q_comp*q_comp*(1.0 - cos(theta)) / (x_comp*x_comp*(2.0 - k*k));
                a_comp = x_comp*x_comp / (q_comp*q_comp);
                double a = real(a_comp);
                double p = real(p_comp);

                double f = 1.0 - r2_n / p*(1 - cos(theta));
                double g = r1_n*r2_n*sin(theta) / sqrt(mu_n*p);
                double gdot = 1 - r1_n / p*(1 - cos(theta));

                // Step 4: calculate V1 and V2
                double V1_n[3];
                double V2_n[3];
                V1_n[0] = (R2_n[0] - f * R1_n[0]) / g;
                V1_n[1] = (R2_n[1] - f * R1_n[1]) / g;
                V1_n[2] = (R2_n[2] - f * R1_n[2]) / g;

                V2_n[0] = (gdot*R2_n[0] - R1_n[0]) / g;
                V2_n[1] = (gdot*R2_n[1] - R1_n[1]) / g;
                V2_n[2] = (gdot*R2_n[2] - R1_n[2]) / g;

                //Step 5: transform to input length and time units (final output)
                V1[0] = V1_n[0] * LU / TU;
                V1[1] = V1_n[1] * LU / TU;
                V1[2] = V1_n[2] * LU / TU;

                V2[0] = V2_n[0] * LU / TU;
                V2[1] = V2_n[1] * LU / TU;
                V2[2] = V2_n[2] * LU / TU;

                double ecc = sqrt(1.0 - p / a); //remove
                //cout << "error: " << error << "\n"; //remove
                //cout << "iterations:                                            " << iterations << "\n"; //remove
                //cout << "eccentricity: " << ecc << "\n"; //remove
                //cout << "\n"; //remove
            }

            //fast computation of acosh
            double acoshAR(const double& b)
            {
                return log(b + sqrt(b*b - 1.0));
            }


            // Newtons method to find k_bi given tau, S, Nrev, and initial guess for k_bi 
            double compute_kb(const double &k_bGuess,
                const double &tau,
                const double &S,
                const double &Nrev,
                const double &tolerance,
                const double &max_iterations,
                const double &sq2,
                const double &eps)
            {
                // initialize
                int iterations = 0;
                double deltak = 1.0e+10;
                double k = k_bGuess;

                // perform iteration loop
                while (iterations < max_iterations)
                {
                    //Step 1.1 increment the iterations counter
                    ++iterations;

                    // Step 1.2 bound k 
                    // ensure k is not less than -sqrt(2)
                    if (k < -sq2)
                    {
                        k = -sq2 + 0.00001;
                    }
                    if (k > sq2)
                    {
                        k = sq2 - 0.00001;
                    }

                    //Step 1.3 compute W, dW, ddW
                    double m = 2 - k*k;
                    double sgnk = k >= 0 ? 1.0 : -1.0;
                    double W, dW, ddW;
                    W = compute_W(k, m, Nrev);
                    dW = (-2.0 + 3.0 * W * k) / m;
                    ddW = (5.0 * dW * k + 3.0 * W) / m;

                    //Step 2.3 compute TOFc, dTOFc, ddTOFc
                    double TOFc = S * sqrt(1 - k*tau) * (tau + (1 - k*tau) * W);
                    double c = (1 - k * tau) / tau;
                    double sqrtctau = sqrt(1 - k*tau);
                    double dTOFc = -TOFc / (2.0 * c) + S * tau * sqrtctau * (dW * c - W);

                    // check for convergence
                    if (fabs(dTOFc) < tolerance)
                    {
                        break;
                    }

                    double ddTOFc = -TOFc / (4.0 * c*c) + S * tau * sqrtctau * (W / c + c*ddW - 3.0*dW);

                    //Step 2.4 compute deltak
                    deltak = -1.0*dTOFc / ddTOFc;

                    //Step 2.5 update k from deltak
                    k += deltak;

                }
                return k;
            }


            // function calculate TOF function for a given k, tau,
            double compute_TOF(const double &k,
                const double &S,
                const double &tau,
                const double &W)
            {
                double TOFc = S*sqrt(1.0 - k*tau)*(tau + (1.0 - k*tau)*W);
                return TOFc;
            }



            //fast, rough computation of acos
            double acosAR(const double& x)
            {
                double coeff;
                double fx = fabs(x);
                double sgnx = x >= 0.0 ? 1.0 : -1.0;

                if (fx <= 0.6)
                    coeff = (0.000014773722 + (1.1782782 - 0.52020038 * fx) * fx) / (1.1793469 + (-0.53277664 - 0.14454764 * fx) * fx);
                else if (fx <= 0.97)
                    coeff = (0.011101554 + (8.9810074 + (-14.816468 + 5.9249913 * fx) * fx) * fx) / (9.2299851 + (-16.001036 + 6.8381053 * fx) * fx);
                else if (fx <= 0.99)
                    coeff = (-35.750586 + (107.24325 - 70.780244 * fx) * fx) / (27.105764 - 26.638535 * fx);

                else
                    coeff = asin(fx);

                return math::PIover2 - sgnx * coeff;
            }



            //function to compute the parameter W
            double compute_W(const double& k,
                const double& m,
                const int& Nrev)
            {
                const double sq2 = sqrt(2.0);
                const double eps = 2.0e-2;
                int sgnk = k < 0.0 ? -1 : 1;
                //cout << "k, Wcompute: " << k << "\n"; //remove
                if (-sq2 <= k && k < sq2 - eps) //elliptical orbit case
                    return ((1 - sgnk) * math::PI + sgnk*acos(1 - m) + 2 * math::PI*Nrev) / sqrt(m*m*m) - k / m;
                else if (k < sq2 && Nrev > 0)
                    return ((1 - sgnk) * math::PI + sgnk*acos(1 - m) + 2 * math::PI*Nrev) / sqrt(m*m*m) - k / m;
                else if (k > sq2 + eps) //hyperbolic orbits
                    return -1.0*acoshAR(1 - m) / sqrt(-m*m*m) - k / m;
                else if (sq2 - eps <= k && k <= sq2 + eps) //Nrev = 0 case
                {
                    double v = k - sq2;
                    double v2 = v*v;
                    double v3 = v*v2;
                    double v4 = v3*v;
                    double v5 = v4*v;
                    double v6 = v5*v;
                    double v7 = v6*v;
                    double v8 = v7*v;
                    return sq2 / 3.0
                        - v / 5.0
                        + sq2 * v2 * (2.0 / 35.0)
                        - v3 * (2.0 / 63.0)
                        + sq2 * v4 * (2.0 / 231.0)
                        - v5 * (2.0 / 429.0)
                        + sq2 * v6 * (8.0 / 6435.0)
                        - v7 * (8.0 / 12155.0)
                        + sq2 * v8 * (8.0 / 46189.0);
                }
                else
                {
                    //cout << "Error on W compute *************************" << k << "\n"; //remove
                    throw 200000;
                    return math::LARGE;
                }
            }

        }
    }
} //close namespace