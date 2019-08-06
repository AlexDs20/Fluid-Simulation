#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

#include "matrix.h"

//--------------------------------------------------
//  Constructors
Matrix::Matrix(int w, int h, Real inVal){
  width = w;
  height = h;
  val.resize(w*h,inVal);
}
Matrix::Matrix(int w, int h){
  width = w;
  height = h;
  val.resize(w*h);
}

//--------------------------------------------------
//  Linear index of the (i,j)
int Matrix::idx(int i, int j){
  if (i>width-1 || j>height-1){
    std::cout << "requested index out of bound!" << std::endl;
  }
  return i+j*width;
}

//--------------------------------------------------
//  Element i,j
Real Matrix::get(int i, int j){
  return val[idx(i,j)];
}

//--------------------------------------------------
//  Return all values
vector<Real> Matrix::getAll(){
  return val;
}

//--------------------------------------------------
//  Set Specific value
void Matrix::set(int i, int j, Real inVal){
  val[idx(i,j)] = inVal;
}

//--------------------------------------------------
//  Maximum
Real Matrix::max(){
  return *std::max_element(std::begin(val),std::end(val));
}
//  Minimum
Real Matrix::min(){
  return *std::min_element(std::begin(val),std::end(val));
}

//  Maximum in absolute value
Real Matrix::amax(){
  return std::max(std::abs(Matrix::max()),std::abs(Matrix::min()));
}
