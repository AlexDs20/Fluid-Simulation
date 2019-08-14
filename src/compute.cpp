#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "compute.h"
#include "boundary.h"
#include "matrix.h"
#include "obstacle.h"

using namespace std;

//--------------------------------------------------
// Write data to file
//--------------------------------------------------
void writeData(std::ofstream& file, std::vector<real> v){
  if (file.is_open()){
    for (auto i = v.begin(); i!=v.end(); i++){
      file << *i << " ";
    }
    file << "\n";
  }else{
    cout << "The file is not open!" << endl;
  }
}
//--------------------------------------------------

//--------------------------------------------------
// Calculates time step (eq.3.50, page 39)
//--------------------------------------------------
real computeDT(Parameters p, real umax, real vmax){
  real A = (p.Re/2.0) * 1.0/( 1.0/(p.dx*p.dx) + 1.0/(p.dy*p.dy) );
  real B = p.dx/umax;
  real C = p.dy/vmax;
  real D = p.dt;
  real out = D;
  if (p.tau > 0){
    out =  p.tau*min({A,B,C,D});
  }
  return out;
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute gamma (eq 3.20 page 30)
//--------------------------------------------------
real computeGamma(matrix<real>* U, matrix<real>* V, real dx, real dy, real dt){
  real minU = std::abs(U->min());
  real maxU = U->max();
  real maxAbsU = (minU < maxU)? maxU*dt/dx:minU*dt/dx;

  real minV = std::abs(V->min());
  real maxV = V->max();
  real maxAbsV = (minV < maxV)? maxV*dt/dy:minV*dt/dy;

  return (1-std::max(maxAbsU,maxAbsV))/2;    // So that the value is right between the max and 1
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute F
//--------------------------------------------------
void computeF(Parameters p, real dt, obstacle* Obs, matrix<real>* U, matrix<real>* V, matrix<real>* F){
  real gamma = computeGamma(U,V,p.dx,p.dy,dt);
  real dx2u, dy2u, dxu2, dyuv;
  real idx = 1.0/p.dx;
  real idy = 1.0/p.dy;
  real idx2 = 1.0/(p.dx*p.dx);
  real idy2 = 1.0/(p.dy*p.dy);

  for (int j=1; j<=p.jmax; j++){
    for (int i=1; i<=p.imax-1; i++){

      // Check if fluid cell
      if (Obs->get(i,j)==FC){

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

      // Set boundary values for F
      }else{
        setObsFBoundaries(i,j,Obs->get(i,j),F,U);
      }
    }
    F->set(0,j,      U->get(0,j));
    F->set(p.imax,j, U->get(p.imax,j));
  }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute G
//--------------------------------------------------
void computeG(Parameters p, real dt, obstacle* Obs, matrix<real>* U, matrix<real>* V, matrix<real>* G){
  real gamma = computeGamma(U,V,p.dx,p.dy,dt);
  real dx2v, dy2v, dyv2, dxuv;
  real idx = 1.0/p.dx;
  real idy = 1.0/p.dy;
  real idx2 = 1.0/(p.dx*p.dx);
  real idy2 = 1.0/(p.dy*p.dy);

  for (int i=1; i<=p.imax; i++){
    for (int j=1; j<=p.jmax-1; j++){

      // Check if fluid fluid cell
      if (Obs->get(i,j)==FC){

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

      // Set boundary values for G
      }else{
        setObsGBoundaries(i,j,Obs->get(i,j),G,V);
      }

    }
    G->set(i,0,      V->get(i,0));
    G->set(i,p.jmax, V->get(i,p.jmax));
  }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute F and G
//--------------------------------------------------
void computeFG(Parameters p, real dt, obstacle* Obs, matrix<real>* U, matrix<real>* V,
    matrix<real>* F, matrix<real>* G){

  real gamma = computeGamma(U,V,p.dx,p.dy,dt);
  real dx2u, dy2u, dxu2, dyuv;  // F
  real dx2v, dy2v, dyv2, dxuv;  // G
  real idx = 1.0/p.dx;
  real idy = 1.0/p.dy;
  real idx2 = 1.0/(p.dx*p.dx);
  real idy2 = 1.0/(p.dy*p.dy);

  for (int j=1; j<=p.jmax; j++){
    for (int i=1; i<=p.imax; i++){

      // Check if fluid cell
      if (Obs->get(i,j)==FC){
        //------------------------------
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
        //------------------------------

        //------------------------------
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
        //------------------------------

      // Set inner boundary values for F and G
      }else{
        setObsFGBoundaries(i,j,Obs->get(i,j),F,G,U,V);
      }
    }
    // Set Edge of domain values
    F->set(0,j,      U->get(0,j));
    F->set(p.imax,j, U->get(p.imax,j));
  }
  for (int i=1; i<=p.imax; i++){
    // Set Edge of domain values
    G->set(i,0,      V->get(i,0));
    G->set(i,p.jmax, V->get(i,p.jmax));
  }
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute RHS: 1/dt * ((Fij - Fi-1,j)/dx + (Gij-Gij-1)/dy)
//--------------------------------------------------
void computeRHS(Parameters p, real dt, obstacle* Obs, matrix<real>* F, matrix<real>* G, matrix<real>* RHS){
  real idx = 1.0/p.dx;
  real idy = 1.0/p.dy;
  real idt = 1.0/dt;

  for (int i=1; i<=p.imax; i++){
    for (int j=1; j<=p.jmax; j++){

      // Calculate only if fluid cell
      if (Obs->get(i,j)==FC){
        RHS->set(i,j, idt*( (F->get(i,j)-F->get(i-1,j)) *idx +
                            (G->get(i,j)-G->get(i,j-1)) *idy
                          )
                );
      }

    }
  }

}
//--------------------------------------------------

//--------------------------------------------------
//  Compute pressure at time step +1
//  Use P as input at the time step t and
//  as output for time step t+1
//--------------------------------------------------
void computeP(Parameters p, obstacle* Obs, matrix<real>* rhs, matrix<real>* P){
  real idx2 = 1.0/(p.dx*p.dx);
  real idy2 = 1.0/(p.dy*p.dy);
  real coeff = p.omega/( 2.0*(idx2+idy2) );

  int it = 1;
  real res;

  int numcell=0;

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
    setObsPBoundaries(Obs,P);

    // Compute inside the domain
    for (int i=1; i<=p.imax; i++){
      for (int j=1; j<=p.jmax; j++){
        if (Obs->get(i,j)==FC){
          P->set(i,j,
                  (1.0-p.omega)*P->get(i,j) +
                  coeff*( (P->get(i+1,j) + P->get(i-1,j)) * idx2 +
                          (P->get(i,j+1) + P->get(i,j-1)) * idy2 -
                          rhs->get(i,j) )
                );
          numcell++;
        }
      }
    }
    // Residual using L2 norm:
    for (int i=1; i<=p.imax; i++){
      for (int j=1; j<=p.jmax; j++){
        if (Obs->get(i,j)==FC){
          res = res +pow( ( P->get(i+1,j) - 2.0*P->get(i,j) + P->get(i-1,j) )*idx2
                         +( P->get(i,j+1) - 2.0*P->get(i,j) + P->get(i,j-1) )*idy2
                         -rhs->get(i,j)
                        ,2);
        }
      }
    }
    res = sqrt( res/numcell );

    it++;
  }while ( it<=p.itermax && res > p.eps);
}
//--------------------------------------------------

//--------------------------------------------------
//  Compute u, v at time step +1 (eq 3.34, 3.35, page 34)
//--------------------------------------------------
void computeNewVel(Parameters p, real delt, matrix<real>* F, matrix<real>* G, matrix<real>* P, matrix<real>* U, matrix<real>* V){
  real dtdx = delt/p.dx;
  real dtdy = delt/p.dy;

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
