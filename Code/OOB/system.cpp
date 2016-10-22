#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
//#include <random>
//#include <chrono>
#include <time.h>
#include "planet.h"
#include "solver.h"
using namespace std;


int main(int argc, char * argv[])
{

	int dim = 3, N = 100000;
	double final_time = 1.0;
	double b = atof(argv[1]);
	double conv = 365.25;
	bool energy = true;

	planet Earth(0.000003, 1.0, 0.0, 0.0, 0.0, 2*M_PI, 0.0);
	//planet Earth(0.000003, 0.8757, 0.4827, -0.00018, -0.00856*conv, 0.015*conv, -0.000000846*conv);
	planet Sun(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	planet Jupiter(0.001, 5.2, 0.0, 0.0, 0.0, 2*0.45*M_PI, 0.0);	
	//planet Jupiter(0.001, -5.42, -0.509, 0.1234, 0.00061*conv, -0.007157*conv, 0.0000161*conv);
 	

	/*
	solver earthsun_Euler;
	earthsun_Euler.add(Earth);
	earthsun_Euler.add(Sun);

	earthsun_Euler.ForwardEuler(dim, N, final_time);
	*/
	
	solver earthsun_VV;
	earthsun_VV.add(Sun);
	earthsun_VV.add(Earth);
	earthsun_VV.add(Jupiter);

	earthsun_VV.velVerlet( dim, N, final_time, 1, energy);
	

	return 0;
}