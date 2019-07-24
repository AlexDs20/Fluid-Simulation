#pragma once
#include <vector>
using namespace std;

#include "parameters.h"
#include "matrix.h"

#define Real float

Real calcDT(Parameters p, Real umax, Real vmax);
void writeData(string file, std::vector<Real> v);
//Real computeGamma(std::vector<Real> U, std::vector<Real> V, Real dx, Real dy, Real dt);
std::vector<Real> computeF(Parameters p, Real dt, Matrix U, Matrix V);
std::vector<Real> computeG(Parameters p, Real dt, Matrix U, Matrix V);
