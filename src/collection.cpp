#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "collection.h"

using namespace std;

//-------------------------
// Calculates time step
Real calcDT(Parameters p, Real umax, Real vmax){
  Real A = (p.Re/2)*1/( pow(p.dx,-2)+pow(p.dy,-2) );
  Real B = p.dx/abs(umax);
  Real C = p.dy/abs(vmax);
  return p.tau*min({A,B,C});
}

//-------------------------
// Write data to file
void writeData(string file, std::vector<Real> v){
  ofstream myfile;
  myfile.open(file, ios::out | ios::app);

  if (myfile.is_open()){
    for (auto i = v.begin(); i!=v.end(); ++i){
      myfile << *i << " ";
    }
    myfile << "\n";
  }else{
    cout << "Did not manage to open the output file!" << endl;
  }
}
