#pragma once

#include <string>

using namespace std;

#define Real float

class Parameters{
  public:
    Real xlength, ylength, dx, dy, dt, t_end, tau, eps, omega;
    Real Re, gx, gy, UI,VI,PI, wW, wE, wN, wS;
    int imax, jmax, itermax;

    Parameters(string file);
};
