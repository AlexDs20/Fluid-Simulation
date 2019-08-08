#pragma once

#include <string>

#include "matrix.h"

using namespace std;

#define real float

class Parameters{
  private:
    real rhoInf, L, uInf, pInf;
  public:
    real xlength, ylength, dx, dy, dt, t_end, tau, eps, omega;
    real Re, gx, gy, UI,VI,PI,rho,vis,inflow;
    int wW, wE, wN, wS;
    int imax, jmax, itermax;

    Parameters(string file);
    void setScale(real uMax, real vMax, real pMax);
    void toDimensionless(matrix<real>* U, matrix<real>* V, matrix<real>* P);
    void toDimensional(matrix<real>* U, matrix<real>* V, matrix<real>* P);
};
