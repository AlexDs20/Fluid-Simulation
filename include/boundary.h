#pragma once
#include <iostream>

using namespace std;

#include "matrix.h"
#include "obstacle.h"
#include "parameters.h"

void setBoundaries(const Parameters& p, matrix<real>* const U, matrix<real>* const V);
void setSpecificBoundaries(const Parameters& p, obstacle* const Obs, matrix<real>* const U, matrix<real>* const V);
void setObsVelBoundaries(obstacle* const Obs, matrix<real>* const U, matrix<real>* const V);
void setObsFBoundaries(const int& i, const int& j, const int& bType, matrix<real>* const U, matrix<real>* const F);
void setObsGBoundaries(const int& i, const int& j, const int& bType, matrix<real>* const V, matrix<real>* const G);
void setObsFGBoundaries(const int& i, const int& j, const int& bType, matrix<real>* F, matrix<real>* const G,
                       matrix<real>* const U, matrix<real>* const V);
void setObsPBoundaries(obstacle* const Obs, matrix<real>* const P);
