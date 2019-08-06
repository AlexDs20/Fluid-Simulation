#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "compute.h"
#include "matrix.h"

using namespace std;

//--------------------------------------------------
// Write all output data to files
//--------------------------------------------------
void writeOutput(string file, Matrix* U, Matrix* V, Matrix* P){
  string file1 = file+"U";
  string file2 = file+"V";
  string file3 = file+"P";

  writeData(file1,U->getAll());
  writeData(file2,V->getAll());
  writeData(file3,P->getAll());
}
//--------------------------------------------------

//--------------------------------------------------
// Write data to file
//--------------------------------------------------
void writeData(string file, std::vector<Real> v){
  ofstream myfile;
  myfile.open(file, ios::out | ios::app);

  if (myfile.is_open()){
    for (auto i = v.begin(); i!=v.end(); i++){
      myfile << *i << " ";
    }
    myfile << "\n";
  }else{
    cout << "Did not manage to open the output file!" << endl;
  }
}
//--------------------------------------------------

//--------------------------------------------------
// Calculates time step (eq.3.50, page 39)
//--------------------------------------------------
Real computeDT(Parameters p, Real umax, Real vmax){
  Real A = (p.Re/2.0) * 1.0/( 1.0/(p.dx*p.dx) + 1.0/(p.dy*p.dy) );
  Real B = p.dx/umax;
  Real C = p.dy/vmax;
  Real D = p.dt;
  Real out = D;
  if (p.tau > 0){
    out =  p.tau*min({A,B,C,D});
  }
  return out;
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute gamma (eq 3.20 page 30)
//--------------------------------------------------
Real computeGamma(Matrix* U, Matrix * V, Real dx, Real dy, Real dt){
  Real minU = std::abs(U->min());
  Real maxU = U->max();
  Real maxAbsU = (minU < maxU)? maxU*dt/dx:minU*dt/dx;

  Real minV = std::abs(V->min());
  Real maxV = V->max();
  Real maxAbsV = (minV < maxV)? maxV*dt/dy:minV*dt/dy;

  return (1-std::max(maxAbsU,maxAbsV))/2;    // So that the value is right between the max and 1
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute F
//--------------------------------------------------
void computeF(Parameters p, Real dt, Matrix* U, Matrix* V, Matrix* F){
  Real gamma = computeGamma(U,V,p.dx,p.dy,dt);
  Real dx2u, dy2u, dxu2, dyuv;
  Real idx = 1.0/p.dx;
  Real idy = 1.0/p.dy;
  Real idx2 = 1.0/(p.dx*p.dx);
  Real idy2 = 1.0/(p.dy*p.dy);

  for (int j=1; j<=p.jmax; j++){
    for (int i=1; i<=p.imax-1; i++){
      // Computes d^2u/dx^2
      dx2u = ( U->get(i+1,j) - 2.0*U->get(i,j) + U->get(i-1,j) ) * idx2;
      // Computes d^2u/dy^2
      dy2u = ( U->get(i,j+1) - 2.0*U->get(i,j) + U->get(i,j-1) ) * idy2;
      // Computes du^2/dx
      dxu2 = (idx/4.0)*
        ( pow( U->get(i,j)+U->get(i+1,j) ,2) - pow( U->get(i-1,j)+U->get(i,j) ,2) )
        + (gamma*idx/4.0)*
        ( abs(U->get(i,j)   + U->get(i+1,j)) * (U->get(i,j)   - U->get(i+1,j)) -
          abs(U->get(i-1,j) + U->get(i,j)  ) * (U->get(i-1,j) - U->get(i,j)  )
        );
      // Computes duv/dy
      dyuv = (idy/4.0)*
        ( (V->get(i,j)   + V->get(i+1,j)  ) * (U->get(i,j)   + U->get(i,j+1)) -
          (V->get(i,j-1) + V->get(i+1,j-1)) * (U->get(i,j-1) + U->get(i,j)  )
        )
        + (gamma*idy/4.0)*
        ( abs(V->get(i,j)   + V->get(i+1,j)  ) * (U->get(i,j)   - U->get(i,j+1)) -
          abs(V->get(i,j-1) + V->get(i+1,j-1)) * (U->get(i,j-1) - U->get(i,j)  )
        );

      // Computes F
      F->set(i,j, U->get(i,j)+dt*( (1.0/p.Re)*( dx2u+dy2u ) -dxu2-dyuv+p.gx ));
    }
    F->set(0,j,      U->get(0,j));
    F->set(p.imax,j, U->get(p.imax,j));
  }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute G
//--------------------------------------------------
void computeG(Parameters p, Real dt, Matrix* U, Matrix* V, Matrix* G){
  Real gamma = computeGamma(U,V,p.dx,p.dy,dt);
  Real dx2v, dy2v, dyv2, dxuv;
  Real idx = 1.0/p.dx;
  Real idy = 1.0/p.dy;
  Real idx2 = 1.0/(p.dx*p.dx);
  Real idy2 = 1.0/(p.dy*p.dy);

  for (int i=1; i<=p.imax; i++){
    for (int j=1; j<=p.jmax-1; j++){
      // Computes d^2v/dx^2
      dx2v = ( V->get(i+1,j) - 2.0*V->get(i,j) + V->get(i-1,j) ) * idx2;
      // Computes d^2v/dy^2
      dy2v = ( V->get(i,j+1) - 2.0*V->get(i,j) + V->get(i,j-1) ) * idy2;
      // Computes dv^2/dy
      dyv2 = (idy/4.0)*
        ( pow( V->get(i,j)+V->get(i,j+1) ,2) - pow( V->get(i,j-1)+V->get(i,j) ,2) )
        + (gamma*idy/4.0)*
        ( abs( V->get(i,j)   + V->get(i,j+1)) * (V->get(i,j)   - V->get(i,j+1)) -
          abs( V->get(i,j-1) + V->get(i,j)  ) * (V->get(i,j-1) - V->get(i,j)  )
        );
      // Computes duv/dx
      dxuv = (idx/4.0)*
        ( (U->get(i,j)   + U->get(i,j+1)  ) * (V->get(i,j)   + V->get(i+1,j)) -
          (U->get(i-1,j) + U->get(i-1,j+1)) * (V->get(i-1,j) + V->get(i,j)  )
        )
        + (gamma*idx/4.0)*
        ( abs( U->get(i,j)   + U->get(i,j+1)   ) * (V->get(i,j)   - V->get(i+1,j)) -
          abs( U->get(i-1,j) + U->get(i-1,j+1) ) * (V->get(i-1,j) - V->get(i,j)  )
        );

      // Computes G
      G->set(i,j, V->get(i,j)+dt* ((1.0/p.Re)*(dx2v+dy2v) -dyv2-dxuv+p.gy) );
    }
    G->set(i,0,      V->get(i,0));
    G->set(i,p.jmax, V->get(i,p.jmax));
  }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute RHS: 1/dt * ((Fij - Fi-1,j)/dx + (Gij-Gij-1)/dy)
//--------------------------------------------------
void computeRHS(Parameters p, Real dt, Matrix* F, Matrix* G, Matrix* RHS){
  Real idx = 1.0/p.dx;
  Real idy = 1.0/p.dy;
  Real idt = 1.0/dt;

  for (int i=1; i<=p.imax; i++){
    for (int j=1; j<=p.jmax; j++){
      RHS->set(i,j, idt*( (F->get(i,j)-F->get(i-1,j)) *idx +
                          (G->get(i,j)-G->get(i,j-1)) *idy
                        )
              );
    }
  }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute pressure at time step +1
// Use pt1 as input as the time step t and
// as output for time step t+1
//--------------------------------------------------
void computeP(Parameters p, Matrix* rhs, Matrix* P){
  Real idx2 = 1.0/(p.dx*p.dx);
  Real idy2 = 1.0/(p.dy*p.dy);
  Real coeff = p.omega/( 2.0*(idx2+idy2) );

  int it = 1;
  Real res;

  do{
    res = 0.0;

    // Boundary conditions to be computed before the inner cells
    // eq. 3.48 page 38
    // Bottom & Top
    for (int i=1; i<=p.imax; i++){
      P->set(i,0,        P->get(i,1));
      P->set(i,p.jmax+1, P->get(i,p.jmax));
    }
    // Left & Right
    for (int j=1; j<=p.jmax; j++){
      P->set(0,j,        P->get(1,j));
      P->set(p.imax+1,j, P->get(p.imax,j));
    }

    // Compute inside the domain
    for (int i=1; i<=p.imax; i++){
      for (int j=1; j<=p.jmax; j++){
        P->set(i,j,
                (1.0-p.omega)*P->get(i,j) +
                coeff*( (P->get(i+1,j) + P->get(i-1,j)) * idx2 +
                        (P->get(i,j+1) + P->get(i,j-1)) * idy2 -
                        rhs->get(i,j) )
              );
      }
    }
    // Residual using L2 norm:
    for (int i=1; i<=p.imax; i++){
      for (int j=1; j<=p.jmax; j++){
        res = res +pow( ( P->get(i+1,j) - 2.0*P->get(i,j) + P->get(i-1,j) )*idx2
                       +( P->get(i,j+1) - 2.0*P->get(i,j) + P->get(i,j-1) )*idy2
                       -rhs->get(i,j)
                      ,2);
      }
    }
    res = sqrt( res/(p.imax*p.jmax) );

    it++;
  }while ( it<=p.itermax && res > p.eps);
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute u, v at time step +1 (eq 3.34, 3.35, page 34)
//--------------------------------------------------
void computeNewVel(Parameters p, Real delt, Matrix* F, Matrix* G, Matrix* P, Matrix* U, Matrix* V){
  Real dtdx = delt/p.dx;
  Real dtdy = delt/p.dy;

  for (int i=1 ; i<=p.imax-1; i++){
    for (int j=1; j<=p.jmax; j++){
      U->set(i,j, F->get(i,j) - dtdx*(P->get(i+1,j)-P->get(i,j)));
    }
  }
  for (int i=1 ; i<=p.imax; i++){
    for (int j=1; j<=p.jmax-1; j++){
      V->set(i,j, G->get(i,j) - dtdy*(P->get(i,j+1)-P->get(i,j)));
    }
  }
}
