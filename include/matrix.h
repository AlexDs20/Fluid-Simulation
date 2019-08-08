#pragma once
#include <string>
#include <vector>

#include "definitions.h"

using namespace std;

template <typename T>
class matrix{
  private:
    void assignValues(string file);

  protected:
    int width, height;
    int idx(int i, int j);
    vector<T> val;

  public:
    matrix(int w, int h);   // simple mem alloc
    matrix(int w, int h, T inVal); // mem alloc with all values set to inVal
    matrix(int w, int h, string file);// mem alloc with values given in the file

    T get(int i, int j);
    vector<T> getAll();
    void set(int i, int j, T inVal);   // Set specific index val
    T max();
    T min();
    T amax();
};
