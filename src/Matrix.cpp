#include <iostream>
#include "Matrix.h"
using namespace std;

#define Real float

Matrix::Matrix(int w, int h){
  width = w;
  height = h;
  val = new Real[w*h];
}

Matrix::~Matrix(void){
  delete [] val;
}

// Set all the values to the initial values
// These will be changed depending on the boundary conditions
void Matrix::setVal(Real in){
  inVal = in;
  for (int i = 0; i<width; ++i){
    for (int j = 0; j<height; ++j){
      val[idx(i,j)] = inVal;
    };
  };
}

// Return the linear index of the pair i,j
int Matrix::idx(int i, int j){
  return i+j*width;
}

// Return the element i,j of the Matrix
Real Matrix::getVal(int i, int j){
  return val[idx(i,j)];
}
