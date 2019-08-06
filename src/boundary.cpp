#include <iostream>

using namespace std;

#include "boundary.h"

void setBoundaries(Parameters p, matrix* U, matrix* V){
  real inflow = 10.0;

// Problem in the corners:
// if a corner has 1 bonudary with no slip and one with free slip,
// the value of the velocity for the boundary bin will not necessarily be correct,
// it depends on the order in which the boundaries are treated ...
// An easy fix, though calculating a lot for nothing, is to run this twice.
// This is what the "it loop" is for.
  for (int it=1;it<3;it++){
    //-------------------------
    //  Left and Right
    for (int j=1;j<p.jmax+1;j++){
      switch (p.wW){
        case 1:                 // No slip
          U->set(0,j,         0.0);
          V->set(0,j,        -V->get(1,j));
          break;
        case 2:                 // Free slip
          U->set(0,j,         0.0);
          V->set(0,j,         V->get(1,j));
          break;
        case 3:                 // Outflow
          U->set(0,j,         U->get(1,j));
          V->set(0,j,         V->get(1,j));
          break;
        case 4:                 // Inflow
          U->set(0,j,         std::abs(inflow));
          V->set(0,j,        -V->get(1,j));
          break;
        case 5:                 // Periodic
          U->set(0,j,         U->get(p.imax-1,j));
          U->set(p.imax,j,    U->get(1,j));
          V->set(0,j,         V->get(p.imax-1,j));
          V->set(1,j,         V->get(p.imax,j));
          V->set(p.imax+1,j,  V->get(2,j));
          break;
      }
      switch (p.wE){
        case 1:
          U->set(p.imax,j,    0.0);
          V->set(p.imax+1,j, -V->get(p.imax,j));
          break;
        case 2:
          U->set(p.imax,j,    0.0);
          V->set(p.imax+1,j,  V->get(p.imax,j));
          break;
        case 3:
          U->set(p.imax,j,    U->get(p.imax-1,j));
          V->set(p.imax+1,j,  V->get(p.imax,j));
          break;
        case 4:
          U->set(p.imax,j,   -std::abs(inflow));
          V->set(p.imax+1,j, -V->get(p.imax,j));
          break;
      }
    }
    //-------------------------
    //  Top And Bottom
    for (int i=1;i<p.imax+1;i++){
      switch (p.wN){
        case 1:
          U->set(i,p.jmax+1, -U->get(i,p.jmax));
          V->set(i,p.jmax,    0.0);
          break;
        case 2:
          U->set(i,p.jmax+1,  U->get(i,p.jmax));
          V->set(i,p.jmax,    0.0);
          break;
        case 3:
          U->set(i,p.jmax+1,  U->get(i,p.jmax));
          V->set(i,p.jmax,    V->get(i,p.jmax-1));
          break;
        case 4:
          U->set(i,p.jmax+1, -U->get(i,p.jmax));
          V->set(i,p.jmax,   -std::abs(inflow));
          break;
        case 5:
          U->set(i,0,         U->get(i,p.jmax-1));
          U->set(i,1,         U->get(i,p.jmax));
          U->set(i,p.jmax+1,  U->get(i,2));
          V->set(i,0,         V->get(i,p.jmax-1));
          V->set(i,p.jmax,    V->get(i,1));
          break;
      }
      switch (p.wS){
        case 1:
          U->set(i,0,        -U->get(i,1));
          V->set(i,0,         0.0);
          break;
        case 2:
          U->set(i,0,         U->get(i,1));
          V->set(i,0,         0.0);
          break;
        case 3:
          U->set(i,0,         U->get(i,1));
          V->set(i,0,         V->get(i,1));
          break;
        case 4:
          U->set(i,0,        -U->get(i,1));
          V->set(i,0,         std::abs(inflow));
          break;
      }
    }
  }
}

void setSpecificBoundaries(Parameters p, matrix* U, matrix* V){
  if (true){
    real inflow = 10;
    for (int j=p.jmax/2-2;j<p.jmax/2+2;j++){
      U->set(0,j,         std::abs(inflow));
      V->set(0,j,        -V->get(1,j));
    }
  }
  if (false){
    for (int i=1;i<=p.imax;i++){
      V->set(i,0, 0);
      V->set(i,p.jmax,0);
    }
    for (int j=1;j<=p.jmax;j++){
      U->set(0,j, 0);
      U->set(p.imax,j, 0);
    }
    for (int i=1;i<=p.imax;i++){
      U->set(i,0, -U->get(i,1));
      U->set(i,p.jmax+1, 2-U->get(i,p.jmax));
    }
    for (int j=1;j<=p.jmax;j++){
      V->set(0,j, -V->get(1,j));
      V->set(p.imax+1,j, -V->get(p.imax,j));
    }
  }
}
