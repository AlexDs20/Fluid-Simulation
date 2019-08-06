#pragma once
#include <string>

#include "definitions.h"
#include "matrix.h"

using namespace std;

class Parameters{
  public:
    Real xlength, ylength, dx, dy, dt, t_end, tau, eps, omega;
    Real Re, gx, gy, UI,VI,PI,rho,vis;
    Real L, uInf, vInf, pInf, rhoInf;
    int wW, wE, wN, wS, imax, jmax, itermax;

    Parameters(string file);
    void setScale(Real umax, Real vmax, Real pmax);
    void toDimensionless(Matrix* U, Matrix* V, Matrix* P);
    void toDimensional(Matrix* U, Matrix* V, Matrix* P);
};
