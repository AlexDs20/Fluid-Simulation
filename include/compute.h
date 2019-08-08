#pragma once
#include <vector>
using namespace std;

#include "definitions.h"
#include "parameters.h"
#include "matrix.h"
#include "obstacle.h"

// Compute the time step
real computeDT(Parameters p, real umax, real vmax);
// write vector data to a file
void writeData(string file, std::vector<real> v);
void writeOutput(string file, matrix<real>* U, matrix<real>* V, matrix<real>* P);
// Compute the gamma factor needed in F and G
real computeGamma(matrix<real>* U, matrix<real>* V, real dx, real dy, real dt);
// Compute F and G which depends on u,v, their derivatives
void computeF(Parameters p, real dt, obstacle* Obs, matrix<real>* U, matrix<real>* V, matrix<real>* F);
void computeG(Parameters p, real dt, obstacle* Obs, matrix<real>* U, matrix<real>* V, matrix<real>* G);
// Compute RHS of the system to obtain the pressure (eq. 3.38)
void computeRHS(Parameters p, real dt, obstacle* Obs, matrix<real>* F, matrix<real>* G,matrix<real>* RHS);
// Compute the value of the pressure for the following time step
// eq 3.44 without eps and with 3.48 boundary conditions
void computeP(Parameters p, obstacle* Obs, matrix<real>* rhs, matrix<real>* pt1);
// Compute the velocities at the following time step
void computeNewVel(Parameters p, real dt, matrix<real>* F, matrix<real>* G,
                         matrix<real>* P, matrix<real>* U, matrix<real>* V);
