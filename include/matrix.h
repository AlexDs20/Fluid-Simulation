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
    // Var
    vector<T> val;
    int width, height;
    // Methods
    int idx(const int i, const int j);

  public:
    // Methods
    // Constructors
    matrix(const int w, const int h);   // simple mem alloc
    matrix(const int w, const int h, const T inVal); // mem alloc with all values set to inVal
    matrix(const int w, const int h, const string file);// mem alloc with values given in the file

    T get(const int i, const int j);
    vector<T> getAll();
    void set(const int i, const int j, const T inVal);   // Set specific index val
    T max();
    T min();
    T amax();
};
