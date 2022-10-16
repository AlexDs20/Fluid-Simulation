#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "definitions.h"
#include "obstacle.h"
#include "matrix.h"

using namespace std;

obstacle::obstacle(const int w, const int h, const string file) : matrix<int>::matrix(w, h, file){
  for (int k=0; k<this->width*this->height; k++) {
    switch (this->val[k]){
      case B_FS_N:
        FS_N.push_back(k);
        break;
      case B_FS_S:
        FS_S.push_back(k);
        break;
      case B_FS_W:
        FS_W.push_back(k);
        break;
      case B_FS_E:
        FS_E.push_back(k);
        break;
      case B_FS_NW:
        FS_NW.push_back(k);
        break;
      case B_FS_NE:
        FS_NE.push_back(k);
        break;
      case B_FS_SW:
        FS_SW.push_back(k);
        break;
      case B_FS_SE:
        FS_SE.push_back(k);
        break;
      case B_NS_N:
        NS_N.push_back(k);
        break;
      case B_NS_S:
        NS_S.push_back(k);
        break;
      case B_NS_W:
        NS_W.push_back(k);
        break;
      case B_NS_E:
        NS_E.push_back(k);
        break;
      case B_NS_NW:
        NS_NW.push_back(k);
        break;
      case B_NS_NE:
        NS_NE.push_back(k);
        break;
      case B_NS_SW:
        NS_SW.push_back(k);
        break;
      case B_NS_SE:
        NS_SE.push_back(k);
        break;
    }
  }
}

//--------------------------------------------------
//  Convert the linear index k of the vector
//  into the pair i,j
void obstacle::toIJ(const int k, int* const i, int* const j){
  *j = floor(k/this->width);
  *i = k - *j*this->width;
}


//--------------------------------------------------
//  Get the (i,j) pair for the element idxB in B_FLAG
void obstacle::getIJ(const int B_FLAG, const int idxB, int* const i ,int* const j){
    switch (B_FLAG){
      case B_FS_N:
        toIJ(FS_N[idxB],i,j);
        break;
      case B_FS_S:
        toIJ(FS_S[idxB],i,j);
        break;
      case B_FS_W:
        toIJ(FS_W[idxB],i,j);
        break;
      case B_FS_E:
        toIJ(FS_E[idxB],i,j);
        break;
      case B_FS_NW:
        toIJ(FS_NW[idxB],i,j);
        break;
      case B_FS_NE:
        toIJ(FS_NE[idxB],i,j);
        break;
      case B_FS_SW:
        toIJ(FS_SW[idxB],i,j);
        break;
      case B_FS_SE:
        toIJ(FS_SE[idxB],i,j);
        break;
      case B_NS_N:
        toIJ(NS_N[idxB],i,j);
        break;
      case B_NS_S:
        toIJ(NS_S[idxB],i,j);
        break;
      case B_NS_W:
        toIJ(NS_W[idxB],i,j);
        break;
      case B_NS_E:
        toIJ(NS_E[idxB],i,j);
        break;
      case B_NS_NW:
        toIJ(NS_NW[idxB],i,j);
        break;
      case B_NS_NE:
        toIJ(NS_NE[idxB],i,j);
        break;
      case B_NS_SW:
        toIJ(NS_SW[idxB],i,j);
        break;
      case B_NS_SE:
        toIJ(NS_SE[idxB],i,j);
        break;
    }
}
