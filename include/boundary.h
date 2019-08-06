#pragma once
#include <iostream>

using namespace std;

#include "parameters.h"
#include "matrix.h"

void setBoundaries(Parameters p, Matrix* U, Matrix* V);
void setSpecificBoundaries(Parameters p, Matrix* U, Matrix*);
