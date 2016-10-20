#include "solver.h"
#include "planet.h"

int main(){

	int dim = 2, N = 10000;
	double final_time = 1.0;

	double Pos_init[3], Vel_init[3];

	planet Earth(0.000003, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	planet Sun(1, 0, 0, 0, 0, 0, 0);

	/*for(int i = 0; i < dim; i++){
		Pos_init[0, i] = Earth.position[i];
		Vel_init[0, i] = Earth.velocity[i];
	}*/

	solver earthsun_Euler;
	earthsun_Euler.add(Earth);
	earthsun_Euler.add(Sun)



}