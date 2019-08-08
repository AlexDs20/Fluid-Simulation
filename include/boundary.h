#pragma once
#include <iostream>

using namespace std;

#include "matrix.h"
#include "obstacle.h"
#include "parameters.h"

void setBoundaries(Parameters p, matrix<real>* U, matrix<real>* V);
void setSpecificBoundaries(Parameters p, matrix<real>* U, matrix<real>* V);
void setObstaclesBoundaries(matrix<real>* U, matrix<real>* V, obstacle* Obs);
