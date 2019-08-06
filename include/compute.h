#pragma once
#include <vector>
using namespace std;

#include "definitions.h"
#include "parameters.h"
#include "matrix.h"

// Compute the time step
real computeDT(Parameters p, real umax, real vmax);
// write vector data to a file
void writeData(string file, std::vector<real> v);
void writeOutput(string file, matrix* U, matrix* V, matrix* P);
// Compute the gamma factor needed in F and G
real computeGamma(matrix* U, matrix* V, real dx, real dy, real dt);
// Compute F and G which depends on u,v, their derivatives
void computeF(Parameters p, real dt, matrix* U, matrix* V, matrix* F);
void computeG(Parameters p, real dt, matrix* U, matrix* V, matrix* G);
// Compute RHS of the system to obtain the pressure (eq. 3.38)
void computeRHS(Parameters p, real dt, matrix* F, matrix* G,matrix* RHS);
// Compute the value of the pressure for the following time step
// eq 3.44 without eps and with 3.48 boundary conditions
void computeP(Parameters p, matrix* rhs, matrix* pt1);
// Compute the velocities at the following time step
void computeNewVel(Parameters p, real dt, matrix* F, matrix* G, matrix* P, matrix* U, matrix* V);
