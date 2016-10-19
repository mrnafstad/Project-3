#include <stdlib.h>
#include <math.h>
#include <stdio.h>


double *Acceleration(double *pos){

	double M_earth = 1.0; //Earth Masses
	double GM_sun = 4*M_PI*M_PI;  // AU^3/yr^2

	double r3 = pow(sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]), 3);
	double *Acc = new double[3];

	Acc[0] = -GM_sun*(pos[0])/fabs(r3);
	Acc[1] = -GM_sun*(pos[1])/fabs(r3);
	Acc[2] = -GM_sun*(pos[2])/fabs(r3);

	return Acc;
}

void ForwardEuler(double **Pos, double **Vel, double **Acc, double h, int n){

	for(int i = 0; i < n-1; i++){
		for(int j = 0; j < 3; j++){
			Vel[i+1][j] = Vel[i][j] + h*Acc[i][j];
			Pos[i+1][j] = Pos[i][j] + h*Vel[i][j];
		}
		Acc[i+1] = Acceleration(Pos[i+1]);
	}

	return;
}

void VelocityVerlet(double **Pos, double **Vel, double **Acc, double h, double n){

	for(int i = 0; i < n-1; i++){
		for(int j = 0; j < 3; j++){
			Pos[i+1][j] = Pos[i][j] + h*Vel[i][j] + 0.5*Acc[i][j]*h*h;
		}
		Acc[i+1] = Acceleration(Pos[i+1]);

		for(int j = 0; j < 3; j++){
			Vel[i+1][j] = Vel[i][j] + 0.5*(Acc[i][j] + Acc[i+1][j])*h;
		}
	}

	printf("Verlet\n");

	return;

}


int main(int argc, char* argv[]){

	double **Pos, **Vel, **Acc;
	int n = 10000;

	Pos = new double*[n]; Vel = new double*[n]; Acc = new double*[n];
	for(int i = 0; i < n; i++){
		Pos[i] = new double[3];
		Vel[i] = new double[3];
		Acc[i] = new double[3];
	}

	for(int i = 0; i < n; i++){
		for(int j = 0; j < 3; j++){
			Pos[i][j] = 0.0;
			Vel[i][j] = 0.0;
			Acc[i][j] = 0.0;
		}
	}

	Pos[0][0] = 1.0 ; Vel[0][1] = 2*M_PI ; 
	Acc[0] = Acceleration(Pos[0]);

	for(int j = 0; j < 3; j++){
		printf("%f\n", Acc[0][j]);
	}


	double h = 1.0/(n-1);


	//ForwardEuler(Pos, Vel, Acc, h, n);

	VelocityVerlet(Pos, Vel, Acc, h, n);

	printf("%f\n", Pos[100][0]);
	printf("%f\n", Vel[100][1]);

	/*FILE *fp;

	fp = fopen("jordsol.txt", "w+");


	for (int i = 0; i < n; i++){
		fprintf(fp, "%f %f\n", Pos[i][0], Pos[i][1]);
	}

	fclose(fp);*/



	for(int i = 0; i < n; i++){
		delete[] Pos[i];
		delete[] Vel[i];
		delete[] Acc[i];
	}
	delete[] Pos; delete[] Vel; delete[] Acc;
	


	return 0;
}