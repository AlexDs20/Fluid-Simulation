#include <algorithm>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using namespace std;

#include "boundary.h"
#include "collection.h"
#include "matrix.h"
#include "parameters.h"

#define Real float

int main(int , char *argv[]){
//-------------------------
// MPI STUFF
MPI_Init(NULL, NULL);
int nproc;
MPI_Comm_size(MPI_COMM_WORLD, &nproc);
int pid;
MPI_Comm_rank(MPI_COMM_WORLD, &pid);

//-------------------------
// Internal variables
string inputfile  = argv[1];
string outputfile = argv[2];

//-------------------------
// Parameters
Parameters p(inputfile);
Real delt;
Real t = 0;

//-------------------------
// Set initial UVP arrays
Matrix U(p.imax+2,p.jmax+2,p.UI);
Matrix V(p.imax+2,p.jmax+2,p.VI);
Matrix P(p.imax+2,p.jmax+2,p.PI);

Matrix F((p.imax+2),(p.jmax+2));
Matrix G((p.imax+2),(p.jmax+2));
Matrix RHS((p.imax+2),(p.jmax+2));

//-------------------------
//  Set scale parameters
//  and Reynolds number
p.setScale(std::max({std::abs(U.max()),std::abs(U.min())}),
           std::max({std::abs(V.max()),std::abs(V.min())}),
           P.max());
//-------------------------
//  Go to dimensionless
//  variables
p.toDimensionless(&U, &V, &P);

//-------------------------
// Create a boundary Matrix
/*
To be done: Create a vector that would contain values indicating whether it is a boundary or not.
The boundary values for the field should then be changed accordingly to that vector
*/

while (t<p.t_end){
  //-------------------------
  // Calculate time step
  std::cout << t/p.t_end << "%" << std::endl;
  delt = calcDT(p,U.max(),V.max());

  //-------------------------
  // Boundary Conditions
  setBoundaries(p,&U,&V);
//  U.setUBoundCond(p.wW,p.wE,p.wS,p.wN);
//  V.setVBoundCond(p.wW,p.wE,p.wS,p.wN);

  //-------------------------
  //  Compute F
  computeF(p,delt,&U,&V,&F);
  computeG(p,delt,&U,&V,&G);

  //-------------------------
  //  Compute RHS
  computeRHS(p,delt,&F,&G,&RHS);

  //-------------------------
  //  Compute Pressure at
  //  the next time step
  //
  //  &P serves as input
  //  and output
  //-------------------------
  computePt1(p,&RHS,&P);

  //-------------------------
  //  Compute new velocities
  computeNewVel(p,&F,&G,&P,&U,&V);

  //-------------------------
  // Write output
  writeData(outputfile,U.getAll());

  t += delt;
}

// Split work between processors
if (pid == 1){

}

// Finalize the MPI environment.
MPI_Finalize();
return 0;
}
