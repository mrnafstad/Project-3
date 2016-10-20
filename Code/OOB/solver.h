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
    void ForwardEuler(int dim, int N, double final_time);
    void VelocityVerlet(int dim, int N, double final_time);
<<<<<<< HEAD
    double setup_matrix(int height, int width);
    void delete_matrix(double **matrix);
=======
    void GravitationalForce(planet current, planet other, double Fx, double Fy, double Fz);
    double setup_matrix(int height, int width)
>>>>>>> a8671c564768fee1bd2132e05f2b7a33a225918f
};

#endif //SOLVER_H