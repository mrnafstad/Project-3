
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

void solver::add(planet newplanet)
{
    total_planets += 1;
    total_mass += newplanet.mass;
    all_planets.push_back(newplanet);
}

void solver::addM(planet newplanet)
{
    total_planets +=1;
    all_planets.push_back(newplanet);
}

void Eulerf( int dim, int N, double final_time, planet )
{
	time = 0.0
	double h = final_time/N;
	for (int i = 0; i < N; i++) {
		
	}
}

void velVerlet( int dim, int N, double final_time, )
{
	
}

double ** solver::setup_matrix(int height,int width)
{   // Function to set up a 2D array

    // Set up matrix
    double **matrix;
    matrix = new double*[height];

    // Allocate memory
    for(int i=0;i<height;i++)
        matrix[i] = new double[width];

    // Set values to zero
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}