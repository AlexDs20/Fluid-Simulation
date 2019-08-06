#pragma once
#include <iostream>

using namespace std;

#include "parameters.h"
#include "matrix.h"

void setBoundaries(Parameters p, matrix* U, matrix* V);
void setSpecificBoundaries(Parameters p, matrix* U, matrix*);
