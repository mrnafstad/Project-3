#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double acceleration( double p, double r )
{
	double pi = M_PI;
	double rcube = r*r*r;
	double acc = -4 * pi * p / rcube;
	return acc;
}


void forwardeuler ( double **Position, double **velocity, int N )
{
	//Not even close to being something real
	double h = 1/(double) N;
	double r;
	double ax, ay;
	for ( int i = 0; i < N - 1; i++) {
		r = sqrt( Position [0][i] * Position [0][i] + Position [1][i] * Position [1][i] );
		ax = acceleration( Position [0][i], r );
		ay = acceleration( Position [1][i], r );
		velocity [0][i+1] = velocity[0][i] + h * ax;
		velocity [1][i+1] = velocity[1][i] + h * ay;
		Position [0][i+1] = Position [0][i] + h * velocity [0][i];
		Position [1][i+1] = Position [1][i] + h * velocity [1][i];
	}
	return;	
}

int main( int argc, char * argv[] )
{
	//initialize position matrix and velocity matrix, to go 3d, change first double and add two elements to forloop
	int N = atoi(argv[1]);
	int dim = 2;
	double ** Position, **velocity;
	Position = new double*[dim];
	velocity = new double*[dim];
	for ( int i = 0; i < dim; i++ ) {
		Position[i] = new double [N];
		velocity[i] = new double [N];
	}
	Position[0][0] = 1.0;
	Position[1][0] = 0.0;
	velocity[0][0] = 0.0;
	//velocity[1][0] = //velocity in y-dir, but what should it be?
	forwardeuler ( Position, velocity, N);
	return 0;
}