
#include "planet.h"

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