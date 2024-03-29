#include <cmath>
#include <iostream>

using namespace std;

#include "boundary.h"

void setBoundaries(const Parameters& p, matrix<real>* const U, matrix<real>* const V){

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
          U->set(0,j,         std::abs(p.inflow));
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
          U->set(p.imax,j,   -std::abs(p.inflow));
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
          V->set(i,p.jmax,   -std::abs(p.inflow));
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
          V->set(i,0,         std::abs(p.inflow));
          break;
      }
    }
  }
}

void setSpecificBoundaries(const Parameters& p, obstacle* const Obs, matrix<real>* const U, matrix<real>* const V){
  //  Lid driven flow
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
  // Backward facing step
  if (true){
    for (int j = 1;j<=p.jmax;j++){
      if (Obs->get(1,j)==FC){
          U->set(0,j,         std::abs(p.inflow));
          V->set(0,j,        -V->get(1,j));
      }
    }
  }
  // Narrow inflow on the left boundary
  if (false){
    for (int j=floor((p.jmax+2)/2)-1;j<ceil((p.jmax+2)/2)+1;j++){
      U->set(0,j,         std::abs(p.inflow));
      V->set(0,j,        -V->get(1,j));
    }
  }
}

void setObsVelBoundaries(obstacle* const Obs, matrix<real>* const U, matrix<real>* const V){
  int i,j;
  for (int it=1;it<3;it++){
    //-------------------------------
    //  FREE SLIP CONDITIONS
    //-------------------------------
    // North (i.e. values should be similar as to wS)
    for (int k=0; k<(int)Obs->FS_N.size();k++){
      Obs->getIJ(B_FS_N,k,&i,&j);

      U->set(i,j,         U->get(i,j+1));
      U->set(i-1,j,       U->get(i-1,j+1));
      V->set(i,j,         0.0);
    }
    // South
    for (int k=0; k<(int)Obs->FS_S.size();k++){
      Obs->getIJ(B_FS_S,k,&i,&j);

      U->set(i,j,         U->get(i,j-1));
      U->set(i-1,j,       U->get(i-1,j-1));
      V->set(i,j-1,       0.0);
    }
    // West
    for (int k=0; k<(int)Obs->FS_W.size();k++){
      Obs->getIJ(B_FS_W,k,&i,&j);

      U->set(i-1,j,       0.0);
      V->set(i,j,         V->get(i-1,j));
      V->set(i,j-1,       V->get(i-1,j-1));
    }
    // East
    for (int k=0; k<(int)Obs->FS_E.size();k++){
      Obs->getIJ(B_FS_E,k,&i,&j);

      U->set(i,j,         0.0);
      V->set(i,j,         V->get(i+1,j));
      V->set(i,j-1,       V->get(i+1,j-1));
    }
    // NW
    for (int k=0; k<(int)Obs->FS_NW.size();k++){
      Obs->getIJ(B_FS_NW,k,&i,&j);

      U->set(i-1,j,       0.0);
      U->set(i,j,         U->get(i,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       V->get(i-1,j-1));
    }
    // NE
    for (int k=0; k<(int)Obs->FS_NE.size();k++){
      Obs->getIJ(B_FS_NE,k,&i,&j);

      U->set(i,j,         0.0);
      U->set(i-1,j,       U->get(i-1,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       V->get(i+1,j-1));
    }
    // SW
    for (int k=0; k<(int)Obs->FS_SW.size();k++){
      Obs->getIJ(B_FS_SW,k,&i,&j);

      U->set(i,j,         U->get(i,j-1));
      U->set(i-1,j,       0.0);
      V->set(i,j-1,       0.0);
      V->set(i,j,         V->get(i-1,j));
    }
    // SE
    for (int k=0; k<(int)Obs->FS_SE.size();k++){
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
    for (int k=0; k<(int)Obs->NS_N.size();k++){
      Obs->getIJ(B_NS_N,k,&i,&j);

      U->set(i,j,         -U->get(i,j+1));
      U->set(i-1,j,       -U->get(i-1,j+1));
      V->set(i,j,         0.0);
    }
    // South
    for (int k=0; k<(int)Obs->NS_S.size();k++){
      Obs->getIJ(B_NS_S,k,&i,&j);

      U->set(i,j,         -U->get(i,j-1));
      U->set(i-1,j,       -U->get(i-1,j-1));
      V->set(i,j-1,       0.0);
    }
    // West
    for (int k=0; k<(int)Obs->NS_W.size();k++){
      Obs->getIJ(B_NS_W,k,&i,&j);

      U->set(i-1,j,       0.0);
      V->set(i,j,         -V->get(i-1,j));
      V->set(i,j-1,       -V->get(i-1,j-1));
    }
    // East
    for (int k=0; k<(int)Obs->NS_E.size();k++){
      Obs->getIJ(B_NS_E,k,&i,&j);

      U->set(i,j,         0.0);
      V->set(i,j,         -V->get(i+1,j));
      V->set(i,j-1,       -V->get(i+1,j-1));
    }
    // NW
    for (int k=0; k<(int)Obs->NS_NW.size();k++){
      Obs->getIJ(B_NS_NW,k,&i,&j);

      U->set(i-1,j,       0.0);
      U->set(i,j,         -U->get(i,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       -V->get(i-1,j-1));
    }
    // NE
    for (int k=0; k<(int)Obs->NS_NE.size();k++){
      Obs->getIJ(B_NS_NE,k,&i,&j);

      U->set(i,j,         0.0);
      U->set(i-1,j,       -U->get(i-1,j+1));
      V->set(i,j,         0.0);
      V->set(i,j-1,       -V->get(i+1,j-1));
    }
    // SW
    for (int k=0; k<(int)Obs->NS_SW.size();k++){
      Obs->getIJ(B_NS_SW,k,&i,&j);

      U->set(i,j,         -U->get(i,j-1));
      U->set(i-1,j,       0.0);
      V->set(i,j-1,       0.0);
      V->set(i,j,         -V->get(i-1,j));
    }
    // SE
    for (int k=0; k<(int)Obs->NS_SE.size();k++){
      Obs->getIJ(B_NS_SE,k,&i,&j);

      U->set(i,j,         0.0);
      U->set(i-1,j,       -U->get(i-1,j-1));
      V->set(i,j,         -V->get(i,j+1));
      V->set(i,j-1,       0.0);
    }

  }
}

void setObsFBoundaries(const int& i, const int& j, const int& bType, matrix<real>* const F, matrix<real>* const U){
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

void setObsGBoundaries(const int& i, const int& j, const int& bType, matrix<real>* const G, matrix<real>* const V){
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

// FG inner boundaries
void setObsFGBoundaries(const int& i, const int& j, const int& bType, matrix<real>* const F, matrix<real>* const G,
                        matrix<real>* const U, matrix<real>* const V){
  switch(bType){
    // Free-slip
    case B_FS_N:
      G->set(i,j,     V->get(i,j));
      break;
    case B_FS_S:
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_FS_W:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_FS_E:
      F->set(i,j,     U->get(i,j));
      break;
    case B_FS_NW:
      F->set(i-1,j,   U->get(i-1,j));
      G->set(i,j,     V->get(i,j));
      break;
    case B_FS_NE:
      F->set(i,j,     U->get(i,j));
      G->set(i,j,     V->get(i,j));
      break;
    case B_FS_SW:
      F->set(i-1,j,   U->get(i-1,j));
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_FS_SE:
      F->set(i,j,     U->get(i,j));
      G->set(i,j-1,   V->get(i,j-1));
      break;
    // No-slip
    case B_NS_N:
      G->set(i,j,     V->get(i,j));
      break;
    case B_NS_S:
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_NS_W:
      F->set(i-1,j,   U->get(i-1,j));
      break;
    case B_NS_E:
      F->set(i,j,     U->get(i,j));
      break;
    case B_NS_NW:
      F->set(i-1,j,   U->get(i-1,j));
      G->set(i,j,     V->get(i,j));
      break;
    case B_NS_NE:
      F->set(i,j,     U->get(i,j));
      G->set(i,j,     V->get(i,j));
      break;
    case B_NS_SW:
      F->set(i-1,j,   U->get(i-1,j));
      G->set(i,j-1,   V->get(i,j-1));
      break;
    case B_NS_SE:
      F->set(i,j,     U->get(i,j));
      G->set(i,j-1,   V->get(i,j-1));
      break;
  }
}

void setObsPBoundaries(obstacle* const Obs, matrix<real>* const P){
  int i,j;
  for (int it=1;it<3;it++){
    //-------------------------------
    //  FREE SLIP CONDITIONS
    //-------------------------------
    // North (i.e. values should be similar as to wS)
    for (int k=0; k<(int)Obs->FS_N.size();k++){
      Obs->getIJ(B_FS_N,k,&i,&j);
      P->set(i,j,         P->get(i,j+1));
    }
    // South
    for (int k=0; k<(int)Obs->FS_S.size();k++){
      Obs->getIJ(B_FS_S,k,&i,&j);
      P->set(i,j,         P->get(i,j-1));
    }
    // West
    for (int k=0; k<(int)Obs->FS_W.size();k++){
      Obs->getIJ(B_FS_W,k,&i,&j);
      P->set(i,j,         P->get(i-1,j));
    }
    // East
    for (int k=0; k<(int)Obs->FS_E.size();k++){
      Obs->getIJ(B_FS_E,k,&i,&j);
      P->set(i,j,         P->get(i+1,j));
    }
    // NW
    for (int k=0; k<(int)Obs->FS_NW.size();k++){
      Obs->getIJ(B_FS_NW,k,&i,&j);
      P->set(i,j,         (P->get(i-1,j)+P->get(i,j+1))/2.0 );
    }
    // NE
    for (int k=0; k<(int)Obs->FS_NE.size();k++){
      Obs->getIJ(B_FS_NE,k,&i,&j);
      P->set(i,j,         (P->get(i+1,j)+P->get(i,j+1))/2.0 );
    }
    // SW
    for (int k=0; k<(int)Obs->FS_SW.size();k++){
      Obs->getIJ(B_FS_SW,k,&i,&j);
      P->set(i,j,         (P->get(i-1,j)+P->get(i,j-1))/2.0 );
    }
    // SE
    for (int k=0; k<(int)Obs->FS_SE.size();k++){
      Obs->getIJ(B_FS_SE,k,&i,&j);
      P->set(i,j,         (P->get(i+1,j)+P->get(i,j-1))/2.0 );
    }
    //-------------------------------
    //  NO SLIP CONDITIONS
    //-------------------------------
    // North (i.e. values should be similar as to wS)
    for (int k=0; k<(int)Obs->NS_N.size();k++){
      Obs->getIJ(B_NS_N,k,&i,&j);
      P->set(i,j,         P->get(i,j+1));
    }
    // South
    for (int k=0; k<(int)Obs->NS_S.size();k++){
      Obs->getIJ(B_NS_S,k,&i,&j);
      P->set(i,j,         P->get(i,j-1));
    }
    // West
    for (int k=0; k<(int)Obs->NS_W.size();k++){
      Obs->getIJ(B_NS_W,k,&i,&j);
      P->set(i,j,         P->get(i-1,j));
    }
    // East
    for (int k=0; k<(int)Obs->NS_E.size();k++){
      Obs->getIJ(B_NS_E,k,&i,&j);
      P->set(i,j,         P->get(i+1,j));
    }
    // NW
    for (int k=0; k<(int)Obs->NS_NW.size();k++){
      Obs->getIJ(B_NS_NW,k,&i,&j);
      P->set(i,j,         (P->get(i-1,j)+P->get(i,j+1))/2.0 );
    }
    // NE
    for (int k=0; k<(int)Obs->NS_NE.size();k++){
      Obs->getIJ(B_NS_NE,k,&i,&j);
      P->set(i,j,         (P->get(i+1,j)+P->get(i,j+1))/2.0 );
    }
    // SW
    for (int k=0; k<(int)Obs->NS_SW.size();k++){
      Obs->getIJ(B_NS_SW,k,&i,&j);
      P->set(i,j,         (P->get(i-1,j)+P->get(i,j-1))/2.0 );
    }
    // SE
    for (int k=0; k<(int)Obs->NS_SE.size();k++){
      Obs->getIJ(B_NS_SE,k,&i,&j);
      P->set(i,j,         (P->get(i+1,j)+P->get(i,j-1))/2.0 );
    }

  }
}
