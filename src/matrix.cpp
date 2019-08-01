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

//--------------------------------------------------
//  Set Boundary Conditions for the U field
//  TBD: Later on, I probably want the boundary to be a grid/vector from which I check each value
//--------------------------------------------------
void Matrix::setUBoundCond(int wW, int wE, int wS, int wN){
  Real inflow=10;
  if (wW > 5 || wE > 5 || wS > 5 || wN > 5){
    std::cout << "Invalid boundary request!" << std::endl;
  }
  for (int j = 1; j<height-1; ++j){
    if (wW == 1){             // No slip
      val[idx(0,j)] = 0;
    }else if (wW == 2){       // Free-slip
      val[idx(0,j)] = 0;
    }else if (wW == 3){       // Outflow
      val[idx(0,j)] = val[idx(1,j)];
    }else if (wW == 4){       // Inflow
      val[idx(0,j)] = inflow;
    }else if (wW == 5 && wE == 5){       // Periodic
      val[idx(0,j)] = val[idx(width-3,j)];
      val[idx(width-2,j)] = val[idx(1,j)];
    }
    if (wE == 1){
      val[idx(width-2,j)] = 0;
    }else if (wE == 2){
      val[idx(width-2,j)] = 0;
    }else if (wE == 3){
      val[idx(width-2,j)] = val[idx(width-3,j)];
    }else if (wE == 4){
      val[idx(width-2,j)] = -abs(inflow);
    }
  }
  for (int i = 1; i<width-1; ++i){
    if (wS == 1){
      val[idx(i,0)] = -val[idx(i,1)];
    }else if (wS == 2){
      val[idx(i,0)] = val[idx(i,1)];
    }else if (wS == 3){
      val[idx(i,0)] = val[idx(i,1)];
    }else if (wS == 4){
      val[idx(i,0)] = -val[idx(i,1)];
    }else if (wS == 5 && wN == 5){
      val[idx(i,0)] = val[idx(i,height-3)];
      val[idx(i,1)] = val[idx(i,height-2)];
      val[idx(i,height-1)] = val[idx(i,2)];
    }
    if (wN == 1){
      val[idx(i,height-1)] = -val[idx(i,height-2)];
    }else if (wN == 2){
      val[idx(i,height-1)] = val[idx(i,height-2)];
    }else if (wN == 3){
      val[idx(i,height-1)] = val[idx(i,height-2)];
    }else if (wN == 4){
      val[idx(i,height-1)] = -val[idx(i,height-2)];
    }
  }
}

void Matrix::setVBoundCond(int wW, int wE, int wS, int wN){
   Real inflow=10;
  if (wW > 5 || wE > 5 || wS > 5 || wN > 5){
    std::cout << "Invalid boundary request!" << std::endl;
  }
  for (int i = 1; i<width-1; ++i){
    if (wS == 1){
      val[idx(i,0)] = 0;
    }else if (wS == 2){
      val[idx(i,0)] = 0;
    }else if (wS == 3){
      val[idx(i,0)] = val[idx(i,1)];
    }else if (wS == 4){
      val[idx(i,0)] = abs(inflow);
    }else if (wS == 5 && wN == 5){
      val[idx(i,0)] = val[idx(i,height-3)];
      val[idx(i,height-2)] = val[idx(i,1)];
    }
    if (wN == 1){
      val[idx(i,height-2)] = 0;
    }else if (wN == 2){
      val[idx(i,height-2)] = 0;
    }else if (wN == 3){
      val[idx(i,height-2)] = val[idx(i,height-3)];
    }else if (wN == 4){
      val[idx(i,height-2)] = -abs(inflow);
    }
  }
  for (int j = 1; j<height-1; ++j){
    if (wW == 1){             // No slip
      val[idx(0,j)] = -val[idx(1,j)];
    }else if (wW == 2){       // Free-slip
      val[idx(0,j)] = val[idx(1,j)];
    }else if (wW == 3){       // Outflow
      val[idx(0,j)] = val[idx(1,j)];
    }else if (wW == 4){       // Inflow
      val[idx(0,j)] = -val[idx(1,j)];
    }else if (wW == 5 && wE == 5){       // Inflow
      val[idx(0,j)] = val[idx(width-3,j)];
      val[idx(1,j)] = val[idx(width-2,j)];
      val[idx(width-1,j)] = val[idx(2,j)];
    }
    if (wE == 1){
      val[idx(width-1,j)] = -val[idx(width-2,j)];
    }else if (wE == 2){
      val[idx(width-1,j)] = val[idx(width-2,j)];
    }else if (wE == 3){
      val[idx(width-1,j)] = val[idx(width-2,j)];
    }else if (wE == 4){
      val[idx(width-1,j)] = -val[idx(width-2,j)];
    }
  }
}
