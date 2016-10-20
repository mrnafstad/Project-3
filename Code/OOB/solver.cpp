
#include "solver.h"
#include "planet.h"
#include <iostream>
#include <cmath>
#include "time.h"

solver::solver()
{
    total_planets = 0;
    radius = 100.0;
    total_mass = 0;
    G = 4*M_PI*M_PI;

}

solver::solver( double radi )
{
    total_planets = 0;
    radius = radi;
    total_mass = 0;
    G = 4*M_PI*M_PI;
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

void solver::print_position(std::ofstream &output, int dimension, double time,int number)
{   // Writes mass, position and velocity to a file "output"
    if(dimension > 3 || dimension <= 0) dimension = 3;
    else{
        for(int i=0;i<number;i++){
            planet &current = all_planets[i];
            output << time;
            for(int j=0;j<dimension;j++) output << "\t" << current.position[j];
            //for(int j=0;j<dimension;j++) output << "\t" << current.velocity[j];
            output << std::endl;
        }
    }
}

void solver::velVerlet( int dim, int N, double final_time, int print_number )
{
	//Will only work for binary system atm
	double h = final_time/(double)N;
	double time = 0.0;
	double** acceleration = setup_matrix(total_planets, 3);
	double** acceleration_new = setup_matrix(total_planets, 3);

	double Fx, Fy, Fz, Fxnew, Fynew, Fznew;

	// Create files for data storage
    char *filename = new char[1000];

    sprintf( filename, "clusterVV_%d_%.3f.txt", total_planets, h );
    std::ofstream output_file(filename);

	print_position( output_file, dim, time, print_number );

	while ( time < final_time) {
		for (int nr1 = 0; nr1 < total_planets; nr1++) {
			planet &current = all_planets[nr1];
			Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0;
			for ( int nr2 = nr1 + 1; nr2 < total_planets; nr2++) {
				planet &other = all_planets[nr2];
				GravitationalForce(other, current, Fx, Fy, Fz);
			}
			//next define acceleration, then the real algo
			acceleration[nr1][0] = -Fx/current.mass;
			acceleration[nr1][1] = -Fy/current.mass;
			acceleration[nr1][2] = -Fz/current.mass;
			for ( int j = 0; j < dim; j++ ) {
				current.position[j] += h*current.velocity[j] + 0.5*h*h*acceleration[nr1][j];
			}

			for ( int nr2 = nr1 + 1; nr2 < total_planets; nr2++) {
				planet &other = all_planets[nr2];
				GravitationalForce(other, current, Fxnew, Fynew, Fznew);
			}
			//next define acceleration, then the real algo
			acceleration_new[nr1][0] = -Fxnew/current.mass;
			acceleration_new[nr1][1] = -Fynew/current.mass;
			acceleration_new[nr1][2] = -Fznew/current.mass;
			for ( int j = 0; j < dim; j++ ) {
				current.velocity[j] += 0.5*h*h*( acceleration[nr1][j] + acceleration_new[nr1][j]);
			}
		}
		time += h;
		print_position(output_file, dim, time, print_number );
		

	}
	output_file.close();
	delete_matrix(acceleration);
	delete_matrix(acceleration_new);	
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

/*
void solver::Gravitationalconstant()
{
	G = 4*M_PI*M_PI/32 * radius * radius * radius / mass;
}*/


void solver::ForwardEuler( int dim, int N, double final_time )
{
	double **Acc = setup_matrix(N, dim);
	double **Vel = setup_matrix(N, dim);
	double **Pos = setup_matrix(N, dim);

	double h = final_time/(N -1);

	//MÅ FÅ INN INITIALVERDIEN PÅ EN ELLER ANNEN MÅTE

	for(int i = 0; i < N; i = i + h){
		for(int j = 0; i < dim; j++){
			Vel[i+1][j] = Vel[i][j] + h*Acc[i][j];
			Pos[i+1][j] = Pos[i][j] + h*Vel[i][j];
		}
	}
		
}

void solver::GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz){   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);
    //double smoothing = epsilon*epsilon*epsilon;

    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r));// + smoothing);
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r));// + smoothing);
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r));// + smoothing);
}