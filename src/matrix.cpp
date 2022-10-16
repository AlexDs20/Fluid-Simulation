#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "matrix.h"

using namespace std;

//--------------------------------------------------
//  Template implementation
template matrix<int>::matrix (const int,const int);
template matrix<real>::matrix(const int,const int);
template matrix<int>::matrix (const int,const int,const int);
template matrix<real>::matrix(const int,const int,const real);
template matrix<int>::matrix (const int,const int,const string);
template matrix<real>::matrix(const int,const int,const string);

template int matrix<int>::idx(const int, const int);
template int matrix<real>::idx(const int, const int);

template int matrix<int>::get  (const int, const int);
template real matrix<real>::get(const int, const int);

template vector<int> matrix<int>::getAll ();
template vector<real> matrix<real>::getAll();

template void matrix<int>::set (const int, const int, const int);
template void matrix<real>::set(const int, const int, const real);

template int matrix<int>::max  ();
template real matrix<real>::max();

template int matrix<int>::min  ();
template real matrix<real>::min();

template int matrix<int>::amax  ();
template real matrix<real>::amax();
//--------------------------------------------------


//--------------------------------------------------
//  Constructors
template <typename T>
matrix<T>::matrix(const int w, const int h){
  this->width = w;
  this->height = h;
  this->val.resize(w*h);
}
template <typename T>
matrix<T>::matrix(const int w, const int h, const T inVal){
  this->width = w;
  this->height = h;
  this->val.resize(w*h,inVal);
}
template <typename T>
matrix<T>::matrix(const int w, const int h, string file){
  this->width = w;
  this->height = h;
  this->val.resize(w*h);
  assignValues(file);
}

//--------------------------------------------------
//  Index pair (i,j) to linear
template <typename T>
int matrix<T>::idx(const int i, const int j){
  if (i>width-1 || j>height-1){
    std::cout << "requested index out of bound!" << std::endl;
  }
  return i+j*width;
}

//--------------------------------------------------
//  Element i,j
template <typename T>
T matrix<T>::get(const int i, const int j){
  return val[idx(i,j)];
}

//--------------------------------------------------
//  Return all values
template <typename T>
vector<T> matrix<T>::getAll(){
  return val;
}

//--------------------------------------------------
//  Set Specific value
template <typename T>
void matrix<T>::set(const int i, const int j, T inVal) {
  val[idx(i,j)] = inVal;
}

//--------------------------------------------------
// Assign values everywhere within the domain from the file
template <typename T>
void matrix<T>::assignValues(string file){
  ifstream in;
  int readVal;
  int i = 1;
  int j = 1;

  in.open(file);
  if (in.is_open()){
    while (in >> readVal){
      // Assign the value
      this->val[idx(i,j)] = readVal;

      // Update indices:
      if (i==width-2){
        i = 1;
        j++;
      }else{
        i++;
      }
    }
  in.close();
  }else{
    cout << "Did not manage to open the file to initialise the matrix!" << endl;
  }
}

//--------------------------------------------------
//  Maximum
template <typename T>
T matrix<T>::max(){
  return *std::max_element(std::begin(val),std::end(val));
}
//  Minimum
template <typename T>
T matrix<T>::min(){
  return *std::min_element(std::begin(val),std::end(val));
}

//  Maximum in absolute value
template <typename T>
T matrix<T>::amax(){
  return std::max(std::abs(matrix<T>::max()),std::abs(matrix<T>::min()));
}
