#pragma once
#include <vector>
#define Real float

class Matrix{
  private:
    int width, height;
    Real inVal;
    std::vector<Real> val;
    //Real* val;

  public:
    Matrix(int w, int h, Real inVal);
    int idx(int i, int j);
    Real get(int i, int j);
    vector<Real> getAll();
    void set(int i, int j, Real inVal);   // Set specific index val
    Real max();
    void setUBoundCond(int wW, int wE, int wS, int wN);
    void setVBoundCond(int wW, int wE, int wS, int wN);
};
