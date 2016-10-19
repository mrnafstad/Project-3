#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

//Just for the sun earth system

class solver
{
public:
	friend class planet;

	//properties
	double radius, G;


	//initializer
	solver();
	solver( double radii );

	//functions;
	void Gravitiationalconstant();
    void Eulerf(int dim, int N, double final_time);
    void VelocityVerlet(int dim, int N, double final_time);
};




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

class planet
{
public:
	friend class solver;
	
	//Properties
	double mass;
	double position[3];
	double velocity[3];
	double potential;
	double kinetic;

	//initializers
	planet();
    planet( double M, double x, double y, double z, double vx, double vy, double vz );

	//functions
    double distance(planet otherPlanet);
    double GravitationalForce(planet otherPlanet, double Gconst);
    double Acceleration(planet otherPlanet, double Gconst);
    double KineticEnergy();
    double PotentialEnergy(planet &otherPlanet, double Gconst, double epsilon);
		
};

planet::planet()
{
    mass = 1.;
    position[0] = 1.;
    position[1] = 0.;
    position[2] = 0.;
    velocity[0] = 0.;
    velocity[1] = 0.;
    velocity[2] = 0.;
    potential = 0.;
    kinetic = 0.;
}

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    potential = 0.;
    kinetic = 0.;
}


double planet::distance(planet Sun)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->position[0];
    y1 = this->position[1];
    z1 = this->position[2];

    x2 = Sun.position[0];
    y2 = Sun.position[1];
    z2 = Sun.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
 }