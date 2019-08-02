#pragma once
#include <vector>
#define Real float

class Matrix{
  private:
    int width, height;
    Real inVal;
    std::vector<Real> val;

  public:
    Matrix(int w, int h, Real inVal);
    Matrix(int w, int h);

    int idx(int i, int j);
    Real get(int i, int j);
    std::vector<Real> getAll();
    void set(int i, int j, Real inVal);   // Set specific index val
    Real max();
    Real min();
};
