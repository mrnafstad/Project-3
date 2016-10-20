
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

void solver::delete_matrix(double **matrix)
{   // Function to deallocate memory of a 2D array

    for (int i=0; i<total_planets; i++)
        delete [] matrix[i];
    delete [] matrix;
}

void solver::Gravitationalconstant()
{
	G = 4*M_PI*M_PI/32 * radius * radius * radius / mass;
}

void solver::ForwardEuler( int dim, int N, double final_time )
{
	double **Acc = setup_matrix(N, dim);
	double **Vel = setup_matrix(N, dim);
	double **Pos = setup_matrix(N, dim);

	double h = final_time/(N -1);

	//MÅ FÅ INN INITIALVERDIEN PÅ EN ELLER ANNEN MÅTE

	for(int i = 0; i < N; i = i + h){
		for(int j = 0; i < dim, j++){
			Vel[i+1][j] = Vel[i][j] + h*Acc[i][j];
			Pos[i+1][j] = Pos[i][j] + h*Vel[i][j];
		}
		
		
	}
}