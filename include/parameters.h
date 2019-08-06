#pragma once

#include <string>

#include "matrix.h"

using namespace std;

#define Real float

class Parameters{
  private:
    Real rhoInf, L, uInf, pInf;
  public:
    Real xlength, ylength, dx, dy, dt, t_end, tau, eps, omega;
    Real Re, gx, gy, UI,VI,PI,rho,vis;
    int wW, wE, wN, wS;
    int imax, jmax, itermax;

    Parameters(string file);
    void setScale(Real uMax, Real vMax, Real pMax);
    void toDimensionless(Matrix* U, Matrix* V, Matrix* P);
    void toDimensional(Matrix* U, Matrix* V, Matrix* P);
};
