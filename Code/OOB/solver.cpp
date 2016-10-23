
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
    totalKinetic = 0;
    totalPotential = 0;
    totalAngularMomentum = 0;

}

solver::solver( double radi )
{
    total_planets = 0;
    radius = radi;
    total_mass = 0;
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
    totalAngularMomentum = 0;
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

void solver::velVerlet( int dim, int N, double final_time, int print_number, bool energy, bool stationary, bool relativity)
{
	double time = 0.0;
	double h = final_time/(double)N;

	//planet &earth = all_planets[1];
	planet &sun = all_planets[0];

	double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;
	double acc[3];
	double acc_new[3];

	FILE *fp;
	fp = fopen("VerletTest.txt", "w+");
	
	int counter = 0;

	if(energy) printf("Time       Total Kinetic Energy  Total Potential Energy  Total Angular Momentum\n");

	int j, k;
	if(stationary) j = 1;
	else {
		j = 0;
		for ( int i = 1; i < total_planets; i++ ) {
			for (int k = 0; k < dim; k++ ) {
				planet &thisother = all_planets[i];
				sun.velocity[k] += thisother.velocity[k]*thisother.mass/sun.mass;
			}
		}
	}

	while(time < final_time){

		fprintf(fp, "%f ", time);
		if ( stationary ) j = 1;
		else j = 0;
		
		for ( j; j < total_planets; j++ ) {
			planet &thisplanet = all_planets[j];
			Fx = 0; Fy = 0; Fz = 0;
			

			if (stationary) {
				GravitationalForce(thisplanet, sun, Fx, Fy, Fz, relativity);
				k = 1;
			}
			else k = 0;

			for ( k; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx, Fy, Fz, relativity );
				}

			}

			acc[0] = Fx/thisplanet.mass; acc[1] = Fy/thisplanet.mass; acc[2] = Fz/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.position[i] += h*thisplanet.velocity[i] + 0.5*acc[i]*h*h;
			}
			
			Fx_new = 0; Fy_new = 0; Fz_new = 0;
			if (stationary) {
				GravitationalForce(thisplanet, sun, Fx_new, Fy_new, Fz_new, relativity);
				k = 1;
			}
			else k = 0;			
			
			for ( k; k < total_planets; k++ ) {
				if ( k != j ) {
					planet other_planet = all_planets[k];
					GravitationalForce( thisplanet, other_planet, Fx_new, Fy_new, Fz_new, relativity );
				}

			}
			acc_new[0] = Fx_new/thisplanet.mass; acc_new[1] = Fy_new/thisplanet.mass; acc_new[2] = Fz_new/thisplanet.mass;

			for(int i = 0; i < dim; i++){
				thisplanet.velocity[i] += 0.5*(acc[i] + acc_new[i])*h;
			}

			fprintf(fp, "%f %f %f ", thisplanet.position[0], thisplanet.position[1], thisplanet.position[2]);
		}
		fprintf(fp, "\n");


		if(energy){
			if(counter == 0){
				KineticEnergySystem();
				PotentialEnergySystem(0.0);
				AngularMomentumSystem();
				printf("%f      %e        %e            %e\n", time, totalKinetic, totalPotential, totalAngularMomentum);
			}
		}

		counter += 1;
		if(N/counter == 10) counter = 0;
		time += h;

	}

	fclose(fp);	
}


/*
void solver::Gravitationalconstant()
{
	G = 4*M_PI*M_PI/32 * radius * radius * radius / mass;
}*/


void solver::ForwardEuler( int dim, int N, double final_time, bool relativity)
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
		GravitationalForce(earth, sun, Fx, Fy, Fz, relativity);

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

void solver::GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz, bool relativity){   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);

    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r));
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r));
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r));

    if(relativity){
    	double angmom2 = current.AngularMomentum()*current.AngularMomentum()/(current.mass*current.mass);
    	double c = 63197.8;  //Speed of light AU/yr
    	Fx *= (1 + 3*angmom2/(r*r*c*c));
    	Fy *= (1 + 3*angmom2/(r*r*c*c));
    	Fz *= (1 + 3*angmom2/(r*r*c*c));
    }
}

void solver::KineticEnergySystem()
{
    totalKinetic = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.kinetic = Current.KineticEnergy();
        totalKinetic += Current.kinetic;
    }
}

void solver::PotentialEnergySystem(double epsilon)
{
    totalPotential = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.potential = 0;
    }
    for(int nr1=0;nr1<total_planets;nr1++){
        planet &Current = all_planets[nr1];
        for(int nr2=nr1+1;nr2<total_planets;nr2++){
            planet &Other = all_planets[nr2];
            Current.potential += Current.PotentialEnergy(Other,G,epsilon);
            Other.potential += Other.PotentialEnergy(Current,G,epsilon);
            totalPotential += Current.potential;
        }
    }
}

void solver::AngularMomentumSystem(){
	totalAngularMomentum = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.ang_mom = Current.AngularMomentum();
        totalAngularMomentum += Current.ang_mom;
    }
}