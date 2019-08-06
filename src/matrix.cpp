#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

#include "matrix.h"

//--------------------------------------------------
//  Constructors
matrix::matrix(int w, int h, real inVal){
  width = w;
  height = h;
  val.resize(w*h,inVal);
}
matrix::matrix(int w, int h){
  width = w;
  height = h;
  val.resize(w*h);
}

//--------------------------------------------------
//  Linear index of the (i,j)
int matrix::idx(int i, int j){
  if (i>width-1 || j>height-1){
    std::cout << "requested index out of bound!" << std::endl;
  }
  return i+j*width;
}

//--------------------------------------------------
//  Element i,j
real matrix::get(int i, int j){
  return val[idx(i,j)];
}

//--------------------------------------------------
//  Return all values
vector<real> matrix::getAll(){
  return val;
}

//--------------------------------------------------
//  Set Specific value
void matrix::set(int i, int j, real inVal){
  val[idx(i,j)] = inVal;
}

//--------------------------------------------------
//  Maximum
real matrix::max(){
  return *std::max_element(std::begin(val),std::end(val));
}
//  Minimum
real matrix::min(){
  return *std::min_element(std::begin(val),std::end(val));
}

//  Maximum in absolute value
real matrix::amax(){
  return std::max(std::abs(matrix::max()),std::abs(matrix::min()));
}
