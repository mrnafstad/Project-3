
#include "solver.h"
#include "planet.h"
#include <iostream>
#include <cmath>
#include "time.h"

solver::solver()
{
	radius = 1.0;
	G = 4*M_PI*M_PI;
	Kinetic = 0.0;
	Potential = 0.0;
}

solver::solver( double radi )
{
	radius = radi;
	G = 4*M_PI*M_PI;
	Kinetic = 0.0;
	Potential = 0.0;
}

void solver::Gravitationalconstant()
{
	G = 4*M_PI*M_PI/32 * radius * radius * radius / mass;
}

void Eulerf( int dim, int N, double final_time )
{
	time = 0.0
	double h = final_time/N;
	while ( time < final_time) {
		
		
	}
}