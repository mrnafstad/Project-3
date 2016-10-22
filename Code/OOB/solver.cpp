
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
{/*
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
	
	FILE *fp;
	fp = fopen("VerletTest.txt", "w+");
	while ( time < final_time) {
		for (int nr1 = 0; nr1 < total_planets; nr1++) {
			planet &current = all_planets[nr1];
			Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0;
			for ( int nr2 = nr1 + 1; nr2 < total_planets; nr2++) {
				planet &other = all_planets[nr2];
				GravitationalForce(current, other, Fx, Fy, Fz);
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
				GravitationalForce(current, other, Fxnew, Fynew, Fznew);
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
		//print_position(output_file, dim, time, print_number );
		fprintf(fp, "%f %f %f %f\n", time, current.position[0], current.position[1], earth.position[2]);		

	}
	fclose(fp);

	//output_file.close();
	delete_matrix(acceleration);
	delete_matrix(acceleration_new);*/

	double time = 0.0;
	double h = final_time/(double)N;

	planet &earth = all_planets[1];
	planet &sun = all_planets[0];

	double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;
	double acc[3];
	double acc_new[3];

	FILE *fp;
	fp = fopen("VerletTest.txt", "w+");
	

	while(time < final_time){
		for ( int j = 1; j < total_planets; j++ ) {
			planet &thisplanet = all_planets[j];
			Fx = 0; Fy = 0; Fz = 0;
			GravitationalForce(thisplanet, sun, Fx, Fy, Fz);
			for ( int k = 1; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx, Fy, Fz );
				}

			}

			acc[0] = Fx/thisplanet.mass; acc[1] = Fy/thisplanet.mass; acc[2] = Fz/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.position[i] += h*thisplanet.velocity[i] + 0.5*acc[i]*h*h;
			}

			Fx_new = 0; Fy_new = 0; Fz_new = 0;
			GravitationalForce(thisplanet, sun, Fx_new, Fy_new, Fz_new);
			for ( int k = 1; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx_new, Fy_new, Fz_new );
				}

			}
			acc_new[0] = Fx_new/thisplanet.mass; acc_new[1] = Fy_new/thisplanet.mass; acc_new[2] = Fz_new/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.velocity[i] += 0.5*(acc[i] + acc_new[i])*h;
			}

			fprintf(fp, "%f %f %f %f\n", time, thisplanet.position[0], thisplanet.position[1], thisplanet.position[2]);
		}
		time += h;

	}

	fclose(fp);	
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
	double time = 0.0;
	double h = final_time/(double)N;

	planet &earth = all_planets[0];
	planet &sun = all_planets[1];

	double Fx, Fy, Fz;
	double acc[3];

	FILE *fp;
	fp = fopen("EulerTest.txt", "w+");
	fprintf(fp, "%f %f %f %f\n", time, earth.position[0], earth.position[1], earth.position[2]);

	while(time < final_time){

		Fx = 0, Fy = 0, Fz = 0;
		GravitationalForce(earth, sun, Fx, Fy, Fz);

		acc[0] = Fx/earth.mass; acc[1] = Fy/earth.mass; acc[2] = Fz/earth.mass;

		for(int i = 0; i < dim; i++){

			earth.position[i] += h*earth.velocity[i];
			earth.velocity[i] += h*acc[i];
		}
		fprintf(fp, "%f %f %f %f\n", time, earth.position[0], earth.position[1], earth.position[2]);
			
		time += h;
	}

	fclose(fp);
		
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