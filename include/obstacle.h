#pragma once

#include <string>
#include <vector>

#include "matrix.h"

using namespace std;

class obstacle : public matrix<int>{
  private:
    void toIJ(int k, int* i, int* j);

  public:
    // Vectors that will contain the indices of the boundary
    vector<int> FS_N;
    vector<int> FS_S;
    vector<int> FS_W;
    vector<int> FS_E;
    vector<int> FS_NW;
    vector<int> FS_NE;
    vector<int> FS_SW;
    vector<int> FS_SE;

    vector<int> NS_N;
    vector<int> NS_S;
    vector<int> NS_W;
    vector<int> NS_E;
    vector<int> NS_NW;
    vector<int> NS_NE;
    vector<int> NS_SW;
    vector<int> NS_SE;


    obstacle(int w, int h, string file);

    void getIJ(int B_FLAG, int idxInB, int*i, int* j);
};
