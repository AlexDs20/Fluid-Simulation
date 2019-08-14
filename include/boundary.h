#pragma once
#include <iostream>

using namespace std;

#include "matrix.h"
#include "obstacle.h"
#include "parameters.h"

void setBoundaries(Parameters p, matrix<real>* U, matrix<real>* V);
void setSpecificBoundaries(Parameters p, obstacle* Obs, matrix<real>* U, matrix<real>* V);
void setObsVelBoundaries(obstacle* Obs, matrix<real>* U, matrix<real>* V);
void setObsFBoundaries(int i, int j, int bType, matrix<real>* U, matrix<real>* F);
void setObsGBoundaries(int i, int j, int bType, matrix<real>* V, matrix<real>* G);
void setObsFGBoundaries(int i, int j, int bType, matrix<real>* U, matrix<real>* V, matrix<real>* F, matrix<real>* G);
void setObsPBoundaries(obstacle* Obs, matrix<real>* P);
