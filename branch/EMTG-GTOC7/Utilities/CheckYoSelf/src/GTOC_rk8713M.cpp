// adaptive_step_int by Alexander Ghosh
// Modified by Donald Ellison
// March 5th 2014

#include "GTOC7_solution_check.h"

std::vector <double> rk8713M(std::vector <double> x_left, std::vector <double> Tvec, std::vector <double> f1, double & h, int & ns, double & error, double & DU, double & TU, double & mu_sun)
{
	std::vector <double> x_right (ns, 0.0);

	//double * f1 = new double[ns]; this is passed in

	std::vector <double> f2(ns, 0.0);
	std::vector <double> f3(ns, 0.0);
	std::vector <double> f4(ns, 0.0);
	std::vector <double> f5(ns, 0.0);
	std::vector <double> f6(ns, 0.0);
	std::vector <double> f7(ns, 0.0);
	std::vector <double> f8(ns, 0.0);
	std::vector <double> f9(ns, 0.0);
	std::vector <double> f10(ns, 0.0);
	std::vector <double> f11(ns, 0.0);
	std::vector <double> f12(ns, 0.0);
	std::vector <double> f13(ns, 0.0);
	std::vector <double> y(ns, 0.0);
	/*double * f2 = new double[ns];
	double * f3 = new double[ns];
	double * f4 = new double[ns];
	double * f5 = new double[ns];
	double * f6 = new double[ns];
	double * f7 = new double[ns];
	double * f8 = new double[ns];
	double * f9 = new double[ns];
	double * f10 = new double[ns];
	double * f11 = new double[ns];
	double * f12 = new double[ns];
	double * f13 = new double[ns];*/
	//double * y = new double[ns];

	error = 0.0; //sufficiently high

	//prep all local variables, note we're forward using variables that we won't need until later to minimize memory usage


	const static double a21 = 1.0 / 18.0;
	const static double c2 = 1.0 / 18.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + (h*a21)*f1[i];
	}

	f2 = GTOC7EOM(y, Tvec, DU, TU, mu_sun); //This is the call to the EOM function, replace this call with the call to your EOM function of choice.

	const static double a31 = 1.0 / 48.0;
	const static double a32 = 1.0 / 16.0;
	const static double c3 = 1.0 / 12.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a31*f1[i] + h*a32*f2[i];
	}

	f3 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a41 = 1.0 / 32.0;
	//const static double a42 = 0;
	const static double a43 = 3.0 / 32.0;
	const static double c4 = 1.0 / 8.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + (h*a41)*f1[i] + h*a43*f3[i];
	}

	f4 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a51 = 5.0 / 16.0;
	//const static double a52 = 0;
	const static double a53 = -75.0 / 64.0;
	const static double a54 = 75.0 / 64.0;
	const static double c5 = 5.0 / 16.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a51*f1[i] + h*a53*f3[i] + h*a54*f4[i];
	}

	f5 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a61 = 3.0 / 80.0;
	//const static double a62 = 0;
	//const static double a63 = 0;
	const static double a64 = 3.0 / 16.0;
	const static double a65 = 3.0 / 20.0;
	const static double c6 = 3.0 / 8.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a61*f1[i] + h*a64*f4[i] + h*a65*f5[i];
	}

	f6 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a71 = 29443841.0 / 614563906.0;
	//const static double a72 = 0;
	//const static double a73 = 0;
	const static double a74 = 77736538.0 / 692538347.0;
	const static double a75 = -28693883.0 / 1125000000.0;
	const static double a76 = 23124283.0 / 1800000000.0;
	const static double c7 = 59.0 / 400.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a71*f1[i] + h*a74*f4[i] + h*a75*f5[i] + h*a76*f6[i];
	}

	f7 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a81 = 16016141.0 / 946692911.0;
	//const static double a82 = 0;
	//const static double a83 = 0;
	const static double a84 = 61564180.0 / 158732637.0;
	const static double a85 = 22789713.0 / 633445777.0;
	const static double a86 = 545815736.0 / 2771057229.0;
	const static double a87 = -180193667.0 / 1043307555.0;
	const static double c8 = 93.0 / 200.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a81*f1[i] + h*a84*f4[i] + h*a85*f5[i] + h*a86*f6[i] + h*a87*f7[i];
	}

	f8 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a91 = 39632708.0 / 573591083.0;
	//const static double a92 = 0;
	//const static double a93 = 0;
	const static double a94 = -433636366.0 / 683701615.0;
	const static double a95 = -421739975.0 / 2616292301.0;
	const static double a96 = 100302831.0 / 723423059.0;
	const static double a97 = 790204164.0 / 839813087.0;
	const static double a98 = 800635310.0 / 3783071287.0;
	const static double c9 = 5490023248.0 / 9719169821.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a91*f1[i] + h*a94*f4[i] + h*a95*f5[i] + h*a96*f6[i] + h*a97*f7[i] + h*a98*f8[i];
	}

	f9 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a10_1 = 246121993.0 / 1340847787.0;
	//const static double a10_2 = 0;
	//const static double a10_3 = 0;
	const static double a10_4 = -37695042795.0 / 15268766246.0;
	const static double a10_5 = -309121744.0 / 1061227803.0;
	const static double a10_6 = -12992083.0 / 490766935.0;
	const static double a10_7 = 6005943493.0 / 2108947869.0;
	const static double a10_8 = 393006217.0 / 1396673457.0;
	const static double a10_9 = 123872331.0 / 1001029789.0;
	const static double c10 = 13.0 / 20.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a10_1*f1[i] + h*a10_4*f4[i] + h*a10_5*f5[i] + h*a10_6*f6[i] + h*a10_7*f7[i] + h*a10_8*f8[i] + h*a10_9*f9[i];
	}

	f10 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a11_1 = -1028468189.0 / 846180014.0;
	//const static double a11_2 = 0;
	//const static double a11_3 = 0;
	const static double a11_4 = 8478235783.0 / 508512852.0;
	const static double a11_5 = 1311729495.0 / 1432422823.0;
	const static double a11_6 = -10304129995.0 / 1701304382.0;
	const static double a11_7 = -48777925059.0 / 3047939560.0;
	const static double a11_8 = 15336726248.0 / 1032824649.0;
	const static double a11_9 = -45442868181.0 / 3398467696.0;
	const static double a11_10 = 3065993473.0 / 597172653.0;
	const static double c11 = 1201146811.0 / 1299019798.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a11_1*f1[i] + h*a11_4*f4[i] + h*a11_5*f5[i] + h*a11_6*f6[i] + h*a11_7*f7[i] + h*a11_8*f8[i] + h*a11_9*f9[i] + h*a11_10*f10[i];
	}

	f11 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a12_1 = 185892177.0 / 718116043.0;
	//const static double a12_2 = 0;
	//const static double a12_3 = 0;
	const static double a12_4 = -3185094517.0 / 667107341.0;
	const static double a12_5 = -477755414.0 / 1098053517.0;
	const static double a12_6 = -703635378.0 / 230739211.0;
	const static double a12_7 = 5731566787.0 / 1027545527.0;
	const static double a12_8 = 5232866602.0 / 850066563.0;
	const static double a12_9 = -4093664535.0 / 808688257.0;
	const static double a12_10 = 3962137247.0 / 1805957418.0;
	const static double a12_11 = 65686358.0 / 487910083.0;
	const static double c12 = 1.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a12_1*f1[i] + h*a12_4*f4[i] + h*a12_5*f5[i] + h*a12_6*f6[i] + h*a12_7*f7[i] + h*a12_8*f8[i] + h*a12_9*f9[i] + h*a12_10*f10[i] + h*a12_11*f11[i];
	}

	f12 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);

	const static double a13_1 = 403863854.0 / 491063109.0;
	//const static double a13_2 = 0;
	//const static double a13_3 = 0;
	const static double a13_4 = -5068492393.0 / 434740067.0;
	const static double a13_5 = -411421997.0 / 543043805.0;
	const static double a13_6 = 652783627.0 / 914296604.0;
	const static double a13_7 = 11173962825.0 / 925320556.0;
	const static double a13_8 = -13158990841.0 / 6184727034.0;
	const static double a13_9 = 3936647629.0 / 1978049680.0;
	const static double a13_10 = -160528059.0 / 685178525.0;
	const static double a13_11 = 248638103.0 / 1413531060.0;
	const static double c13 = 1.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*a13_1*f1[i] + h*a13_4*f4[i] + h*a13_5*f5[i] + h*a13_6*f6[i] + h*a13_7*f7[i] + h*a13_8*f8[i] + h*a13_9*f9[i] + h*a13_10*f10[i] + h*a13_11*f11[i];
	}


	f13 = GTOC7EOM(y, Tvec, DU, TU, mu_sun);


	//8th order solution -- do this first because we need x_left unchanged; when I did 5th order first I was then accumulating on the value of x_left twice for y
	const static double b1upper = 14005451.0 / 335480064.0;
	//const static double b2upper = 0;
	//const static double b3upper = 0;
	//const static double b4upper = 0;
	//const static double b5upper = 0;
	const static double b6upper = -59238493.0 / 1068277825.0;
	const static double b7upper = 181606767.0 / 758867731.0;
	const static double b8upper = 561292985.0 / 797845732.0;
	const static double b9upper = -1041891430.0 / 1371343529.0;
	const static double b10upper = 760417239.0 / 1151165299.0;
	const static double b11upper = 118820643.0 / 751138087.0;
	const static double b12upper = -528747749.0 / 2220607170.0;
	const static double b13upper = 1.0 / 4.0;

	for (size_t i = 0; i < ns; ++i)
	{
		y[i] = x_left[i] + h*b1upper*f1[i] + h*b6upper*f6[i] + h*b7upper*f7[i] + h*b8upper*f8[i] + h*b9upper*f9[i] + h*b10upper*f10[i] + h*b11upper*f11[i] + h*b12upper*f12[i] + h*b13upper*f13[i];
	}



	//7th order solution (compared against the 8th order solution to determine the error due to one RK step)
	const static double b1lower = 13451932.0 / 455176623.0;
	//const static double b2lower = 0; 
	//const static double b3lower = 0; 
	//const static double b4lower = 0; 
	//const static double b5lower = 0; 
	const static double b6lower = -808719846.0 / 976000145.0;
	const static double b7lower = 1757004468.0 / 5645159321.0;
	const static double b8lower = 656045339.0 / 265891186.0;
	const static double b9lower = -3867574721.0 / 1518517206.0;
	const static double b10lower = 465885868.0 / 322736535.0;
	const static double b11lower = 53011238.0 / 667516719.0;
	const static double b12lower = 2.0 / 45.0;
	//const static double b13lower = 0;


	for (size_t i = 0; i < ns; ++i)
	{
		//Determine right hand values of states, take 8th order as truth and compare it with 7th order to quantify error
		x_right[i] = x_left[i] + h*b1lower*f1[i] + h*b6lower*f6[i] + h*b7lower*f7[i] + h*b8lower*f8[i] + h*b9lower*f9[i] + h*b10lower*f10[i] + h*b11lower*f11[i] + h*b12lower*f12[i];
		//error = std::max(error, fabs(x_right[i]-y[i]));
		error = error > fabs(x_right[i] - y[i]) ? error : fabs(x_right[i] - y[i]);
	}



	//clean everything up
	/*delete[] y;
	delete[] f13;
	delete[] f12;
	delete[] f11;
	delete[] f10;
	delete[] f9;
	delete[] f8;
	delete[] f7;
	delete[] f6;
	delete[] f5;
	delete[] f4;
	delete[] f3;
	delete[] f2;*/

	return x_right;
}


std::vector <double> GTOC7EOM(std::vector <double> & X, std::vector <double> Tvec, double & DU, double & TU, double & mu_sun)
{
	double r = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	double rcubed = r*r*r;
	double one_over_mass = 1.0 / X[6];
	double T = sqrt(Tvec[0] * Tvec[0] + Tvec[1] * Tvec[1] + Tvec[2] * Tvec[2]);
	double Isp = 3000.0/TU;
	double g0 = 9.80665 / 1000.0 * TU * TU / DU;
	std::vector <double> dX(7, 0.0);

	dX[0] = X[3];
	dX[1] = X[4];
	dX[2] = X[5];
	dX[3] = -mu_sun*X[0] / rcubed + Tvec[0] * one_over_mass;
	dX[4] = -mu_sun*X[1] / rcubed + Tvec[1] * one_over_mass;
	dX[5] = -mu_sun*X[2] / rcubed + Tvec[2] * one_over_mass;
	dX[6] = -T / (g0*Isp);

	return dX;
}




std::vector <double> adaptive_step_int(std::vector <double> x_left, std::vector <double> Tvec, double & h, int & ns, double & precisionTarget, double & DU, double & TU, double & mu_sun)
{
	double nseg = 1; //let the integrator try to make it in one step -- won't happen but Alex magic will compensate
	double  resumeH = h; //set to the first segment's step size
	double accumulatedH = 0.0, effectiveH = h;
	double resumeError = 1.0e+20;
	double precisionError = 1.0e+20;
	double percentagecomplete = 0.0, lastpercentcomplete = 0.0;
	bool last_substep;
	bool print_progress = false;

	double t = 0.0;

	std::vector <double> f1 (ns, 0.0);
	std::vector <double> x_right(ns, 0.0);


	//adaptive step-size integration routine
	for (int segment = 0; segment < nseg; ++segment)
	{
		if (fabs(resumeH) >= 1e-14)
		{   //this deals with an issue where the resume is 0 from the previous step (not sure how that got set) and it kills this new step
			effectiveH = resumeH;
			precisionError = resumeError;
		}

		else if (resumeH == 0.0 && effectiveH == 0.0)
		{   //this situation can come up if the first segment is zeroed out and we have initialized both to 0.
			effectiveH = h; //set it to be h
		}

		if (fabs(effectiveH) > fabs(h)) { //historical step suggests we jump too far
			effectiveH = h;
		};

		accumulatedH = 0.0; //at the beginning of a segment we have not moved along the step yet
		last_substep = false; //at the beginning of a segment we are on the FIRST substep

		//Forward Integration of the states
		do
		{ //loop until we get all the way through a full h (a full segment)

			//-> INSERT EOM 1st call and STORE f1 values. Replace this with your EOM call of choice.	
			f1 = GTOC7EOM(x_left, Tvec, DU, TU, mu_sun);


			//take a trial substep
			do
			{ //cycle until trial substep is done with sufficient accuracy
				//copy in for the current step our left point

				if (!last_substep)
				{
					//no error!  give it a real value so we don't divide by zero.
					if (precisionError == 0.0)
					{
						precisionError = 1e-15; //Almost zero!
					}

					if (precisionError >= precisionTarget)
						effectiveH = 0.98*effectiveH*pow(precisionTarget / precisionError, 0.17);

					else
					{
						effectiveH = 1.01*effectiveH*pow(precisionTarget / precisionError, 0.18);

						//if our increased step kicks us too long, make it shorter anyway and just run it
						if (fabs(h - accumulatedH) < fabs(effectiveH) && !last_substep)
						{
							effectiveH = h - accumulatedH;
						}

					}

					if (fabs(effectiveH) < 1e-10)
					{//H is too small.....Alexing
						std::cout << "H Got too Small.  Aborting (Alexing)." << std::endl;
						throw 13;
					}

				}

				else if (precisionError > precisionTarget && last_substep)
				{
					//we got here because we thought it was the last substep, and upon calculation it 
					//was too big and not precise enough so we have to shrink it
					last_substep = false; //not last step after all
					effectiveH = 0.98*effectiveH*pow(precisionTarget / precisionError, 0.17);
				}


				//rk takes in x1 as the left-hand side, but returns it as the answer (right hand side)
				//rk8713M (x_left, x_right, uleft, f1, effectiveH, precision_error, segment, ns, nc, nseg, sat_num, *design);
				x_right = rk8713M(x_left, Tvec, f1, effectiveH, ns, precisionError, DU, TU, mu_sun);

				//check for NaN's, assign huge error if they happen
				for (size_t j = 0; j < ns; ++j)
				{
					if (x_left[j] != x_left[j] || x_right[j] != x_right[j] || f1[j] != f1[j])
					{
						precisionError = 1.0e-7;
						break;
					}
				}

			} while (precisionError > precisionTarget);


			//trial substep was accurate enough; it becomes the new left
			for (int statenum = 0; statenum < ns; ++statenum)
			{
				x_left[statenum] = x_right[statenum];
			}


			accumulatedH += effectiveH;
			
			if (print_progress)
			{
				percentagecomplete = accumulatedH / h*100.0;
				if (percentagecomplete - lastpercentcomplete >= 1.0) {
					std::cout << accumulatedH / h*100.0 << "% ";
					lastpercentcomplete = percentagecomplete;
				}
			}

			//if our next step will push us over, reduce it down to be as small as necessary to hit target exactly
			if (fabs(h - accumulatedH) < fabs(effectiveH) && fabs(h - accumulatedH) > 0.0 && !last_substep)
			{
				resumeH = effectiveH;
				effectiveH = h - accumulatedH;
				resumeError = precisionError;
				last_substep = true; //assume that the next substep will be the last substep now
			}

			std::cout << effectiveH << ' ' << accumulatedH << std::endl;

		} while (fabs(accumulatedH) < fabs(h));


		t = t + h; //update the independent variable
		//segment = probinfo->nseg;
		/*
		//write all states and the time to file THIS IS VERY COSTLY in terms of execution time and should only be done on a final "plotting" run
		//where all intermediate values of the states are desired
		std::ofstream outputfile(probinfo->datafile, std::ios::app);
		outputfile.precision(10);
		outputfile.width(14); outputfile << std::left << t;
		outputfile.width(14); outputfile << std::left << x_left[0];
		outputfile.width(14); outputfile << std::left << x_left[1];
		outputfile.width(14); outputfile << std::left << x_left[2];
		outputfile.width(14); outputfile << std::left << x_left[3];
		outputfile.width(14); outputfile << std::left << x_left[4];
		outputfile.width(14); outputfile << std::left << x_left[5];
		outputfile << std::endl;
		outputfile.close();
		*/

	}// end for loop 

	return x_right;

	//end of integration routine
}