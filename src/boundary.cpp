#include <iostream>

using namespace std;

#include "boundary.h"

void setBoundaries(Parameters p, matrix<real>* U, matrix<real>* V){
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

void setSpecificBoundaries(Parameters p, matrix<real>* U, matrix<real>* V){
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

void setObsVelBoundaries(obstacle* Obs, matrix<real>* U, matrix<real>* V){
  int i,j;
  for (int it=1;it<3;it++){
    //-------------------------------
    //  FREE SLIP CONDITIONS
    //-------------------------------
    // North (i.e. values should be similar as to wS)
    for (int k=0; k<=(int)Obs->FS_N.size();k++){
      Obs->getIJ(B_FS_N,k,&i,&j);

      U->set(i,j,         U->get(i,j+1));
      U->set(i-1,j,       U->get(i-1,j+1));
      V->set(i,j,         0.0);
    }
    // South
    for (int k=0; k<=(int)Obs->FS_S.size();k++){
      Obs->getIJ(B_FS_S,k,&i,&j);

      U->set(i,j,         U->get(i,j-1));
      U->set(i-1,j,       U->get(i-1,j-1));
      V->set(i,j-1,       0.0);
    }
    // West
    for (int k=0; k<=(int)Obs->FS_W.size();k++){
      Obs->getIJ(B_FS_W,k,&i,&j);

      U->set(i-1,j,       0.0);
      V->set(i,j,         V->get(i-1,j));
      V->set(i,j-1,       V->get(i-1,j-1));
    }
    // East
    for (int k=0; k<=(int)Obs->FS_E.size();k++){
      Obs->getIJ(B_FS_E,k,&i,&j);

      U->set(i,j,         0.0);
      V->set(i,j,         V->get(i+1,j));
      V->set(i,j-1,       V->get(i+1,j-1));
    }
    // NW
    for (int k=0; k<=(int)Obs->FS_NW.size();k++){
      Obs->getIJ(B_FS_NW,k,&i,&j);

      U->set(i-1,j,       0.0);
      U->set(i,j,         U->get(i,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       V->get(i-1,j-1));
    }
    // NE
    for (int k=0; k<=(int)Obs->FS_NE.size();k++){
      Obs->getIJ(B_FS_NE,k,&i,&j);

      U->set(i,j,         0.0);
      U->set(i-1,j,       U->get(i-1,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       V->get(i+1,j-1));
    }
    // SW
    for (int k=0; k<=(int)Obs->FS_SW.size();k++){
      Obs->getIJ(B_FS_SW,k,&i,&j);

      U->set(i,j,         U->get(i,j-1));
      U->set(i-1,j,       0.0);
      V->set(i,j-1,       0.0);
      V->set(i,j,         V->get(i-1,j));
    }
    // SE
    for (int k=0; k<=(int)Obs->FS_SE.size();k++){
      Obs->getIJ(B_FS_SE,k,&i,&j);

      U->set(i,j,         0.0);
      U->set(i-1,j,       U->get(i-1,j-1));
      V->set(i,j,         V->get(i,j+1));
      V->set(i,j-1,       0.0);
    }
    //-------------------------------
    //  NO SLIP CONDITIONS
    //-------------------------------
    // North (i.e. values should be similar as to wS)
    for (int k=0; k<=(int)Obs->NS_N.size();k++){
      Obs->getIJ(B_NS_N,k,&i,&j);

      U->set(i,j,         -U->get(i,j+1));
      U->set(i-1,j,       -U->get(i-1,j+1));
      V->set(i,j,         0.0);
    }
    // South
    for (int k=0; k<=(int)Obs->NS_S.size();k++){
      Obs->getIJ(B_NS_S,k,&i,&j);

      U->set(i,j,         -U->get(i,j-1));
      U->set(i-1,j,       -U->get(i-1,j-1));
      V->set(i,j-1,       0.0);
    }
    // West
    for (int k=0; k<=(int)Obs->NS_W.size();k++){
      Obs->getIJ(B_NS_W,k,&i,&j);

      U->set(i-1,j,       0.0);
      V->set(i,j,         -V->get(i-1,j));
      V->set(i,j-1,       -V->get(i-1,j-1));
    }
    // East
    for (int k=0; k<=(int)Obs->NS_E.size();k++){
      Obs->getIJ(B_NS_E,k,&i,&j);

      U->set(i,j,         0.0);
      V->set(i,j,         -V->get(i+1,j));
      V->set(i,j-1,       -V->get(i+1,j-1));
    }
    // NW
    for (int k=0; k<=(int)Obs->NS_NW.size();k++){
      Obs->getIJ(B_NS_NW,k,&i,&j);

      U->set(i-1,j,       0.0);
      U->set(i,j,         -U->get(i,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       -V->get(i-1,j-1));
    }
    // NE
    for (int k=0; k<=(int)Obs->NS_NE.size();k++){
      Obs->getIJ(B_NS_NE,k,&i,&j);

      U->set(i,j,         0.0);
      U->set(i-1,j,       -U->get(i-1,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       -V->get(i+1,j-1));
    }
    // SW
    for (int k=0; k<=(int)Obs->NS_SW.size();k++){
      Obs->getIJ(B_NS_SW,k,&i,&j);

      U->set(i,j,         -U->get(i,j-1));
      U->set(i-1,j,       0.0);
      V->set(i,j-1,       0.0);
      V->set(i,j,         -V->get(i-1,j));
    }
    // SE
    for (int k=0; k<=(int)Obs->NS_SE.size();k++){
      Obs->getIJ(B_NS_SE,k,&i,&j);

      U->set(i,j,         0.0);
      U->set(i-1,j,       -U->get(i-1,j-1));
      V->set(i,j,         -V->get(i,j+1));
      V->set(i,j-1,       0.0);
    }

  }
}

void setObsFBoundaries(int i, int j, int bType, matrix<real>* F, matrix<real>* U){
  switch(bType){
    case B_FS_W:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_FS_E:
      F->set(i,j,     U->get(i,j));
      break;
    case B_FS_NW:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_FS_NE:
      F->set(i,j,     U->get(i,j));
      break;
    case B_FS_SW:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_FS_SE:
      F->set(i,j,     U->get(i,j));
      break;
    case B_NS_W:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_NS_E:
      F->set(i,j,     U->get(i,j));
      break;
    case B_NS_NW:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_NS_NE:
      F->set(i,j,     U->get(i,j));
      break;
    case B_NS_SW:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_NS_SE:
      F->set(i,j,     U->get(i,j));
      break;
  }
}

void setObsGBoundaries(int i, int j, int bType, matrix<real>* G, matrix<real>* V){
  switch(bType){
    case B_FS_N:
      G->set(i,j,     V->get(i,j));
      break;
    case B_FS_S:
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_FS_NW:
      G->set(i,j,     V->get(i,j));
      break;
    case B_FS_NE:
      G->set(i,j,     V->get(i,j));
      break;
    case B_FS_SW:
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_FS_SE:
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_NS_N:
      G->set(i,j,     V->get(i,j));
      break;
    case B_NS_S:
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_NS_NW:
      G->set(i,j,     V->get(i,j));
      break;
    case B_NS_NE:
      G->set(i,j,     V->get(i,j));
      break;
    case B_NS_SW:
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_NS_SE:
      G->set(i,j-1,   V->get(i,j-1));
      break;
  }
}
