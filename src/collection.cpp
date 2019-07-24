#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "collection.h"
#include "matrix.h"

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

//-------------------------
//  Compute gamma
// Real computeGamma(std::vector<Real> U, std::vector<Real> V, Real dx, Real dy, Real dt){
//   cout *std::max_elements(std::begin(U),std::end(U)) << endl;
//   return 1;
//   // Real safetyFac = 2;
//
//   // Real maxU = *std::max_element(std::begin(std::abs(U)),std::end(std::abs(U)));
//   // maxU = maxU*dt/dx;
//
//   // Real maxV = *std::max_element(std::begin(std::abs(V)),std::end(std::abs(V)));
//   // maxV = maxV*dt/dx;
//
//   // return safetyFac* *std::max(maxU,maxV);
// }


//-------------------------
//  Compute F
std::vector<Real> computeF(Parameters p, Real dt, Matrix U, Matrix V){
  //-------------------------
  //  Compute gamma factor
  Real gamma = 1;
  // Real gamma = computeGamma(U,V,p.dx,p.dy,dt);

  // Initialize output vector
  std::vector<Real> F(U.getAll().size());
  Real dx2u, dy2u, dxu2, dyuv;

    for (int j = 1; j<p.jmax+1; ++j){
      for (int i = 1; i<p.imax; ++i){
        // Computes dx2u
        dx2u = (U.get(i+1,j)-2*U.get(i,j)+U.get(i-1,j))/pow(p.dx,2);
        // Computes dy2u
        dy2u = (U.get(i,j+1)-2*U.get(i,j)+U.get(i,j-1))/pow(p.dy,2);
        // Computes dxu2
        dxu2 = (1/p.dx)*
          ( pow((U.get(i,j)+U.get(i+1,j))/2,2) - pow((U.get(i-1,j)+U.get(i,j))/2,2) )
          + (gamma/p.dx)*
          ( abs(U.get(i,j)+U.get(i+1,j))*(U.get(i,j)-U.get(i+1,j))/4 -
            abs(U.get(i-1,j)+U.get(i,j))*(U.get(i-1,j)-U.get(i,j))/4
          );
        // Computes dyuv
        dyuv = (1/p.dy)*
          ( (V.get(i,j)+V.get(i+1,j)) * (U.get(i,j)+U.get(i,j+1)) / 4 -
            (V.get(i,j-1)+V.get(i+1,j-1)) * (U.get(i,j-1)+U.get(i,j)) / 4
          )
          + (gamma/p.dy)*
          ( abs( (V.get(i,j)+V.get(i+1,j)) ) * (U.get(i,j)-U.get(i,j+1)) / 4 -
            abs( (V.get(i,j-1)+V.get(i+1,j-1)) ) * (U.get(i,j-1)-U.get(i,j)) / 4
          );

        // Computes F
        F[U.idx(i,j)] = U.get(i,j) +
                  dt*( (1/p.Re)* ( dx2u+dy2u ) - dxu2 - dyuv + p.gx );
      }
    }
  return F;
}

//-------------------------
//  Compute G
std::vector<Real> computeG(Parameters p, Real dt, Matrix U, Matrix V){
  //-------------------------
  //  Compute gamma factor
  Real gamma = 1;
  // Real gamma = computeGamma(U,V,p.dx,p.dy,dt);

  // Initialize output vector
  std::vector<Real> G(U.getAll().size());
  Real dx2v, dy2v, dyv2, dxuv;

    for (int j = 1; j<p.jmax; ++j){
      for (int i = 1; i<p.imax+1; ++i){
        // Computes dx2u
        dx2v = (V.get(i+1,j)-2*V.get(i,j)+V.get(i-1,j))/pow(p.dx,2);
        // Computes dy2u
        dy2v = (V.get(i,j+1)-2*V.get(i,j)+V.get(i,j-1))/pow(p.dy,2);
        // Computes dxu2
        dyv2 = (1/p.dy)*
          ( pow((V.get(i,j)+V.get(i,j+1))/2,2) - pow((V.get(i,j-1)+V.get(i,j))/2,2) )
          + (gamma/p.dy)*
          ( abs(V.get(i,j)+V.get(i,j+1))*(V.get(i,j)-V.get(i,j+1))/4 -
            abs(V.get(i,j-1)+V.get(i,j))*(V.get(i,j-1)-V.get(i,j))/4
          );
        // Computes dyuv
        dxuv = (1/p.dx)*
          ( (U.get(i,j)+U.get(i,j+1)) * (V.get(i,j)+V.get(i+1,j)) / 4 -
            (U.get(i-1,j)+U.get(i-1,j+1)) * (V.get(i-1,j)+V.get(i,j)) / 4
          )
          + (gamma/p.dx)*
          ( abs( (U.get(i,j)+U.get(i,j+1)) ) * (V.get(i,j)-V.get(i+1,j)) / 4 -
            abs( (U.get(i-1,j)+U.get(i-1,j+1)) ) * (V.get(i-1,j)-V.get(i,j)) / 4
          );

        // Computes G
        G[U.idx(i,j)] = V.get(i,j) +
                  dt*( (1/p.Re)* ( dx2v+dy2v ) - dyv2 - dxuv + p.gy );
      }
    }
  return G;
}
