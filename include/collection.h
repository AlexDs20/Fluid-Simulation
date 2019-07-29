#pragma once
#include <vector>
using namespace std;

#include "parameters.h"
#include "matrix.h"

#define Real float

// Compute the time step
Real calcDT(Parameters p, Real umax, Real vmax);
// write vector data to a file
void writeData(string file, std::vector<Real> v);
// Compute the gamma factor needed in F and G
Real computeGamma(Matrix* U, Matrix* V, Real dx, Real dy, Real dt);
// Compute F and G which depends on u,v, their derivatives
void computeF(Parameters p, Real dt, Matrix* U, Matrix* V, Matrix* F);
void computeG(Parameters p, Real dt, Matrix* U, Matrix* V, Matrix* G);
// Compute RHS of the system to obtain the pressure (eq. 3.38)
void computeRHS(Parameters p, Real dt, Matrix* F, Matrix* G,Matrix* RHS);
// Compute the value of the pressure for the following time step
// eq 3.44 without eps and with 3.48 boundary conditions
void computePt1(Parameters p, Matrix* rhs, Matrix* pt1);
