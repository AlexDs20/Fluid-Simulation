#pragma once
#include <vector>
using namespace std;

#include "parameters.h"

#define Real float

Real calcDT(Parameters p, Real umax, Real vmax);
void writeData(string file, std::vector<Real> v);
