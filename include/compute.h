#pragma once
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include "definitions.h"
#include "parameters.h"
#include "matrix.h"
#include "obstacle.h"

// Compute the time step
real computeDT(const Parameters& p, const real umax, const real vmax);
// write vector data to a file
void writeData(std::ofstream& file, std::vector<real> v);
// Compute the gamma factor needed in F and G
real computeGamma(matrix<real>* U, matrix<real>* V, real dx, real dy, real dt);
// Compute F and G which depends on u,v, their derivatives
void computeF(const Parameters& p, const real dt, obstacle* const Obs, matrix<real>* const U, matrix<real>* const V, matrix<real>* const F);
void computeG(const Parameters& p, const real dt, obstacle* const Obs, matrix<real>* const U, matrix<real>* const V, matrix<real>* const G);
void computeFG(const Parameters& p, const real dt, obstacle* const Obs, matrix<real>* const U, matrix<real>* const V, matrix<real>* const F, matrix<real>* const G);
// Compute RHS of the system to obtain the pressure (eq. 3.38)
void computeRHS(const Parameters& p, real dt, obstacle* const Obs, matrix<real>* const F, matrix<real>* const G,matrix<real>* const RHS);
// Compute the value of the pressure for the following time step
// eq 3.44 without eps and with 3.48 boundary conditions
void computeP(const Parameters& p, obstacle* Obs, matrix<real>* const rhs, matrix<real>* const pt1);
// Compute the velocities at the following time step
void computeNewVel(const Parameters& p, real dt, matrix<real>* const F, matrix<real>* const G,
                         matrix<real>* const P, matrix<real>* const U, matrix<real>* const V);
