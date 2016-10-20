#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <vector>
#include <fstream>

using std::vector;

class solver
{
public:
	friend class planet;

	//properties
	double radius, G;
	vector<planet> all_planets;
	int total_planets;
	double total_mass;


	//initializer
	solver();
	solver( double radii );

	//functions;
	void Gravitiationalconstant();
	void add(planet newPlanet);
	void addM(planet newPlanet);
	void print_position(std::ofstream &output, int dim, double time, int number);
    void ForwardEuler(int dim, int N, double final_time);
    void velVerlet(int dim, int N, double final_time, int print_number);
    double **setup_matrix(int height, int width);
    void delete_matrix(double **matrix);
    void GravitationalForce(planet current, planet other, double Fx, double Fy, double Fz);
};

#endif //SOLVER_H