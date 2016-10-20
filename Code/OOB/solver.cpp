
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

void velVerlet( int dim, int N, double final_time )
{
	//Will only work for binary system atm
	double h = final_time/(double)N;
	double time = 0.0;
	while ( time < final_time) {
		for (int nr1 = 0; nr1 < total_planets; nr1++) {
			planet &current = all_planets[nr1];
			Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0;
			for ( int nr2 = nr1 + 1; nr2 < total_planets; nr2++) {
				planet &other = all_planets[nr2];
				GravitationalForce(other, current, Fx, Fy, Fz);
			}
			//next define acceleration, then the real algo

		}
	}	
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

void solver::GravitationalForce(planet &current,planet &other,double &Fx,double &Fy,double &Fz)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);
    //double smoothing = epsilon*epsilon*epsilon;

    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r))// + smoothing);
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r))// + smoothing);
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r))// + smoothing);
}