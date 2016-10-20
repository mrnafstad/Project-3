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


	//initializer
	solver();
	solver( double radii );

	//functions;
	void Gravitiationalconstant();
	void add(planet newPlanet);
	void addM(planet newPlanet);
    void Eulerf(int dim, int N, double final_time);
    void VelocityVerlet(int dim, int N, double final_time);
    double setup_matrix(int height, int width)
};

#endif //SOLVER_H