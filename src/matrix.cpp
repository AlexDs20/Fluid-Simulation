#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "matrix.h"

using namespace std;

//--------------------------------------------------
//  Template implementation
template matrix<int>::matrix (int,int);
template matrix<real>::matrix(int,int);
template matrix<int>::matrix (int,int,int);
template matrix<real>::matrix(int,int,real);
template matrix<int>::matrix (int,int,string);
template matrix<real>::matrix(int,int,string);

template int matrix<int>::get(int,int);
template real matrix<real>::get(int,int);

template vector<int> matrix<int>::getAll();
template vector<real> matrix<real>::getAll();

template void matrix<int>::set(int,int,int);
template void matrix<real>::set(int,int,real);

template int matrix<int>::max();
template real matrix<real>::max();

template int matrix<int>::min();
template real matrix<real>::min();

template int matrix<int>::amax();
template real matrix<real>::amax();


//--------------------------------------------------
//  Constructors
template <typename T>
matrix<T>::matrix(int w, int h){
  width = w;
  height = h;
  val.resize(w*h);
}
template <typename T>
matrix<T>::matrix(int w, int h, T inVal){
  width = w;
  height = h;
  val.resize(w*h,inVal);
}
template <typename T>
matrix<T>::matrix(int w, int h, string file){
  width = w;
  height = h;
  val.resize(w*h);
  assignValues(file);
}

//--------------------------------------------------
//  Linear index of the (i,j)
template <typename T>
int matrix<T>::idx(int i, int j){
  if (i>width-1 || j>height-1){
    std::cout << "requested index out of bound!" << std::endl;
  }
  return i+j*width;
}

//--------------------------------------------------
//  Element i,j
template <typename T>
T matrix<T>::get(int i, int j){
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
void matrix<T>::set(int i, int j, T inVal){
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
      // Assign the value to the
      val[idx(i,j)] = readVal;

      // The indices: when reaching the end of the x line (at imax i.e. w-2),
      // jump to the next y line
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
