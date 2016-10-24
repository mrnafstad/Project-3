#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <time.h>
#include "planet.h"
#include "solver.h"
using namespace std;


int main(int argc, char * argv[])
{

	int dim = 3, N = 1e5;
	double final_time = 2.0;
	//double final_time = 20.0;

	double conv = 365.25;

	bool energy = false;
	bool stationary = false;
	bool relativity = false;
	bool MercPeri = false;

	//PLANETS WITH FIXED INITIAL VALUES
	//-----------------------------------
	//planet Sun(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//planet Earth(0.000003, 1.0, 0.0, 0.0, 0.0, 2*M_PI, 0.0);
	//planet Mercury(1.652*pow(10, -7), 0.3075, 0.0, 0.0, 0.0, 12.44, 0.0);
	//planet Jupiter(1.0, 0.0, -5.2, 0.0, 2*5.2*M_PI/11.8, 0.0, 0.0);	
	
	//PLANETS WITH ACTUAL DATA FROM NASA (23. oct 2016 00:00:00)
	//----------------------------------------------------------
	planet Sun(1.0, 3.5493e-03, 3.4641e-03, -1.5946e-04, 0.0, 0.0, 0.0);
	planet Mercury(1.652e-7, -0.3882, -0.132, 0.0247, 3.368e-3*conv, -2.5365e-2*conv, -2.382e-3*conv);
	planet Venus(2.447e-6, 0.3494, -0.6365, -0.0289, 1.765e-2*conv, 9.554e-3*conv, -8.88e-4*conv);
	planet Earth(0.000003, 0.8670, 0.4976, -0.00018, -0.00822*conv, 0.0149*conv, -9.495e-07*conv);
	planet Mars(3.226e-7, 1.223, -0.6442, -0.0436, 7.094e-3*conv, 1.356e-2*conv, 1.1e-4*conv);
	planet Jupiter(0.00095, -5.424, -0.5108, 0.1234, 0.00062*conv, -0.007155*conv, 1.5882e-5*conv);
 	planet Saturn(0.00028, -2.226, -9.785, 0.2587, 5.134e-3*conv, -1.2544e-3*conv, -1.828e-4*conv);
 	planet Uranus(4.365e-5, 18.451, 7.589, -0.2108, -1.5249e-3*conv, 3.454e-3*conv, 3.247e-5*conv);
 	planet Neptune(5.149e-5, 28.27, -9.898, -0.4476, 1.017e-3*conv, 2.981e-3*conv, -8.516e-5*conv);
	
	/*
	solver earthsun_Euler;
	earthsun_Euler.add(Earth);
	earthsun_Euler.add(Sun);

	earthsun_Euler.ForwardEuler(dim, N, final_time);
	*/
	
	solver system_VV;
	system_VV.add(Sun);
	system_VV.add(Mercury);
	system_VV.add(Venus);
	system_VV.add(Earth);
	system_VV.add(Mars);
	system_VV.add(Jupiter);
	system_VV.add(Saturn);
	system_VV.add(Uranus);
	system_VV.add(Neptune);

	system_VV.velVerlet( dim, N, final_time, energy, stationary, relativity, MercPeri);
	

	return 0;
}