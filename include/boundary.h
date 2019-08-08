#pragma once
#include <iostream>

using namespace std;

#include "parameters.h"
#include "matrix.h"

void setBoundaries(Parameters p, matrix<real>* U, matrix<real>* V);
void setSpecificBoundaries(Parameters p, matrix<real>* U, matrix<real>* V);
