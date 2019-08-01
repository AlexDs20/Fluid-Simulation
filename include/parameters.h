#pragma once

#include <string>

#include "matrix.h"

using namespace std;

#define Real float


class Parameters{
  public:
    Real xlength, ylength, dx, dy, dt, t_end, tau, eps, omega;
    Real gx, gy, UI,VI,PI, rho, vis, wW, wE, wN, wS;
    Real Re, L, uInf, pInf, rhoInf;
    int imax, jmax, itermax;

    Parameters(string file);
    void setScale(Real uMax, Real vMax, Real pMax);
    void toDimensionless(Matrix* U, Matrix* V, Matrix* P);
    void toDimensional(Matrix* U, Matrix* V, Matrix* P);
};
