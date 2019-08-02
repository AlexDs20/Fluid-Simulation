#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "collection.h"
#include "matrix.h"

using namespace std;

//--------------------------------------------------
// Calculates time step
//--------------------------------------------------
Real calcDT(Parameters p, Real umax, Real vmax){
  Real A = (p.Re/2.0) * 1.0/( 1.0/(p.dx*p.dx) + 1.0/(p.dy*p.dy) );
  Real B = p.dx/abs(umax);
  Real C = p.dy/abs(vmax);
  Real D = p.dt;
  Real out = D;
  if (p.tau > 0){
    out =  p.tau*min({A,B,C,D});
  }
  return out;
}
//--------------------------------------------------

//--------------------------------------------------
// Write data to file
//--------------------------------------------------
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
//--------------------------------------------------

//--------------------------------------------------
//  Compute gamma
//--------------------------------------------------
Real computeGamma(Matrix* U, Matrix * V, Real dx, Real dy, Real dt){
  Real safetyFac = 2;

  Real minU = std::abs(U->min());
  Real maxU = std::abs(U->max());
  Real maxAbsU = (minU < maxU)? maxU*dt/dx:minU*dt/dx;

  Real minV = std::abs(V->min());
  Real maxV = std::abs(V->max());
  Real maxAbsV = (minV < maxV)? maxV*dt/dy:minV*dt/dy;

  return safetyFac* std::fmax(maxAbsU,maxAbsV);
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute F
//--------------------------------------------------
void computeF(Parameters p, Real dt, Matrix* U, Matrix* V, Matrix* F){
  Real gamma = computeGamma(U,V,p.dx,p.dy,dt);
  Real dx2u, dy2u, dxu2, dyuv;

    for (int j = 1; j<p.jmax+1; ++j){
      for (int i = 1; i<p.imax; ++i){
        // Computes dx2u
        dx2u = (U->get(i+1,j)-2*U->get(i,j)+U->get(i-1,j))/pow(p.dx,2);
        // Computes dy2u
        dy2u = (U->get(i,j+1)-2*U->get(i,j)+U->get(i,j-1))/pow(p.dy,2);
        // Computes dxu2
        dxu2 = (1/p.dx)*
          ( pow((U->get(i,j)+U->get(i+1,j))/2,2) - pow((U->get(i-1,j)+U->get(i,j))/2,2) )
          + (gamma/p.dx)*
          ( abs(U->get(i,j)+U->get(i+1,j))*(U->get(i,j)-U->get(i+1,j))/4 -
            abs(U->get(i-1,j)+U->get(i,j))*(U->get(i-1,j)-U->get(i,j))/4
          );
        // Computes dyuv
        dyuv = (1/p.dy)*
          ( (V->get(i,j)   + V->get(i+1,j)  ) * (U->get(i,j)   + U->get(i,j+1)) -
            (V->get(i,j-1) + V->get(i+1,j-1)) * (U->get(i,j-1) + U->get(i,j)  )
          ) / 4
          + (gamma/p.dy)*
          ( abs( (V->get(i,j)   + V->get(i+1,j)  ) ) * (U->get(i,j)   - U->get(i,j+1)) -
            abs( (V->get(i,j-1) + V->get(i+1,j-1)) ) * (U->get(i,j-1) - U->get(i,j)  )
          ) / 4;

        // Computes F
        F->set(i,j, U->get(i,j)+dt*( (1/p.Re)*( dx2u+dy2u ) -dxu2-dyuv+p.gx ));
      }
      F->set(0,j,U->get(0,j));
      F->set(p.imax,j,U->get(p.imax,j));
    }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute G
//--------------------------------------------------
void computeG(Parameters p, Real dt, Matrix* U, Matrix* V, Matrix* G){
  Real gamma = computeGamma(U,V,p.dx,p.dy,dt);
  Real dx2v, dy2v, dyv2, dxuv;

    for (int i = 1; i<p.imax+1; ++i){
      for (int j = 1; j<p.jmax; ++j){
        // Computes dx2u = (d^2/dx^2)u
        dx2v = (V->get(i+1,j)-2*V->get(i,j)+V->get(i-1,j))/pow(p.dx,2);
        // Computes dy2u
        dy2v = (V->get(i,j+1)-2*V->get(i,j)+V->get(i,j-1))/pow(p.dy,2);
        // Computes dxu2
        dyv2 = (1/p.dy)*
          ( pow((V->get(i,j)+V->get(i,j+1))/2,2) - pow((V->get(i,j-1)+V->get(i,j))/2,2) )
          + (gamma/p.dy)*
          ( abs(V->get(i,j)+V->get(i,j+1))*(V->get(i,j)-V->get(i,j+1))/4 -
            abs(V->get(i,j-1)+V->get(i,j))*(V->get(i,j-1)-V->get(i,j))/4
          );
        // Computes dyuv
        dxuv = (1/p.dx)*
          ( (U->get(i,j)+U->get(i,j+1)) * (V->get(i,j)+V->get(i+1,j)) / 4 -
            (U->get(i-1,j)+U->get(i-1,j+1)) * (V->get(i-1,j)+V->get(i,j)) / 4
          )
          + (gamma/p.dx)*
          ( abs( (U->get(i,j)+U->get(i,j+1)) ) * (V->get(i,j)-V->get(i+1,j)) / 4 -
            abs( (U->get(i-1,j)+U->get(i-1,j+1)) ) * (V->get(i-1,j)-V->get(i,j)) / 4
          );

        // Computes G
        G->set(i,j,V->get(i,j)+dt*( (1/p.Re)*(dx2v+dy2v) -dyv2-dxuv+p.gy ) );
      }
      G->set(i,0,V->get(i,0));
      G->set(i,p.jmax,V->get(i,p.jmax));
    }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute RHS
//--------------------------------------------------
void computeRHS(Parameters p, Real dt, Matrix* F, Matrix* G, Matrix* RHS){

  // Initialize output vector
    for (int i = 1; i<p.imax+1; ++i){
      for (int j = 1; j<p.jmax+1; ++j){
        RHS->set(i,j, 1/dt*( (F->get(i,j)-F->get(i-1,j))/p.dx +
              (G->get(i,j) - G->get(i,j-1))/p.dy ) );
      }
    }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute pt1 = pressure time step +1
// Use pt1 as input as the time step t and
// as output for time step t+1
//--------------------------------------------------
void computePt1(Parameters p, Matrix* rhs, Matrix* pt1){
  Matrix pt(p.imax+2,p.jmax+2);
  pt = *pt1;
  int it = 1;
  Real res=0;

  do{
    // Set boundary conditions before the iterations (eq. 3.48)
    for (int j = 1; j<= p.jmax; ++j){
      pt1->set(0,j,pt.get(1,j));
      pt1->set(p.imax+1,j,pt.get(p.imax,j));
    }
    for (int i = 1; i<= p.imax; ++i){
      pt1->set(i,0,pt.get(i,1));
      pt1->set(i,p.jmax+1,pt.get(i,p.jmax));
    }
    for (int i = 1; i <= p.imax; ++i){
      for (int j = 1; j <= p.jmax; ++j){
        pt1->set(i,j, (1-p.omega)*pt.get(i,j) +
                p.omega/( (1+1)/pow(p.dx,2) + (1+1)/pow(p.dy,2) ) *
                ( (1*pt.get(i+1,j)+1*pt1->get(i-1,j))/pow(p.dx,2) +
                  (1*pt.get(i,j+1)+1*pt1->get(i,j-1))/pow(p.dy,2) -
                  rhs->get(i,j)
                )
            );
        // Residual:
        // If L2 norm
        res = res +pow( ( pt.get(i+1,j) - 2*pt.get(i,j) + pt.get(i-1,j) )/(p.dx*p.dx)
                       +( pt.get(i,j+1) - 2*pt.get(i,j) + pt.get(i,j-1) )/(p.dy*p.dy)
                       -rhs->get(i,j)
                      ,2);
      }
    }
    res = sqrt( res/(p.imax*p.jmax) );

    // reassign the names: pt as to become the pt1 that was just calculated
    pt = *pt1;
    it++;
  }while ( it<=p.itermax && res > p.eps);
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute u, v at time step +1
//--------------------------------------------------
void computeNewVel(Parameters p, Matrix* F, Matrix* G, Matrix* P, Matrix* U, Matrix* V){
  for (int i = 1 ; i<=p.imax-1; ++i){
    for (int j = 1; j<=p.jmax; ++j){
      U->set(i,j, F->get(i,j) - (p.dt/p.dx)*(P->get(i+1,j)-P->get(i,j)));
    }
  }
  for (int i = 1 ; i<=p.imax; ++i){
    for (int j = 1; j<=p.jmax-1; ++j){
      V->set(i,j, G->get(i,j) - (p.dt/p.dy)*(P->get(i,j+1)-P->get(i,j)));
    }
  }
}
