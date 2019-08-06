#pragma once
#include <vector>
#include "definitions.h"

class matrix{
  private:
    int width, height;
    real inVal;
    std::vector<real> val;
    int idx(int i, int j);

  public:
    matrix(int w, int h, real inVal);
    matrix(int w, int h);

    real get(int i, int j);
    std::vector<real> getAll();
    void set(int i, int j, real inVal);   // Set specific index val
    real max();
    real min();
    real amax();
};
