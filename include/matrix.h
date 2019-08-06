#pragma once
#include <vector>
#include "definitions.h"

class Matrix{
  private:
    int width, height;
    Real inVal;
    std::vector<Real> val;
    int idx(int i, int j);

  public:
    Matrix(int w, int h, Real inVal);
    Matrix(int w, int h);

    Real get(int i, int j);
    std::vector<Real> getAll();
    void set(int i, int j, Real inVal);   // Set specific index val
    Real max();
    Real min();
    Real amax();
};
