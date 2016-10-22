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

	int dim = 3, N = 10000;
	double final_time = 1.0;
	double b = atof(argv[1]);

	planet Earth(0.000003, 1.0, 0.0, 0.0, 0.0, 2*M_PI, 0.0);
	planet Sun(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);	
	planet Jupiter(0.001, 5.2, 0.0, 0.0, 0.0, 2*0.45*M_PI, 0.0);
 	

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

	earthsun_VV.velVerlet( dim, N, final_time, 1);
	

	return 0;
}