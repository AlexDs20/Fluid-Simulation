#include <algorithm>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using namespace std;

#include "definitions.h"
#include "boundary.h"
#include "compute.h"
#include "matrix.h"
#include "parameters.h"

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
// p.setScale(std::max(U.max(),std::abs(U.min())),
//            std::max(V.max(),std::abs(V.min())),
//            P.max());
//-------------------------
//  Go to dimensionless
//  variables
// p.toDimensionless(&U, &V, &P);

//-------------------------
// Create a boundary Matrix
/*
To be done: Create a matrix that would contain values indicating whether it is a boundary or not.
The boundary values for the field should then be changed accordingly to the boudaries
*/

while (t<p.t_end){
  //-------------------------
  // Calculate time step
  std::cout << t*100.0/p.t_end << "%" << std::endl;
  delt = computeDT(p,U.amax(),V.amax());

  //-------------------------
  // Boundary Conditions
  setBoundaries(p,&U,&V);
  setSpecificBoundaries(p,&U,&V);


  //-------------------------
  //  Compute F and G
  computeF(p,delt,&U,&V,&F);
  computeG(p,delt,&U,&V,&G);

  //-------------------------
  //  Compute RHS
  computeRHS(p,delt,&F,&G,&RHS);

  //-------------------------
  //  Push Pressure to
  //  the next time step
  computeP(p,&RHS,&P);

  //-------------------------
  //  Compute new velocities
  computeNewVel(p,delt,&F,&G,&P,&U,&V);

  //-------------------------
  // Write output
  writeOutput(outputfile,&U,&V,&P);

  t += delt;
}

// Split work between processors
if (pid == 1){

}

// Finalize the MPI environment.
MPI_Finalize();
return 0;
}
