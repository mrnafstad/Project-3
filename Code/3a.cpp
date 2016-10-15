#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

void forwardeuler ( double **Position, double **velocity, int N )
{
	//Not even close to being something real
	double h = 1/(double) N
	double a;
	for ( int i = 0; i < N - 1; i++) {
		a = acceleration( velocity[1][i])
		velocity [1][i+1] = velocity[1][i] + h * a
		Position [1][i+1] = Position [1][i] + h * velocity [1][i]
	}
	return;	
}