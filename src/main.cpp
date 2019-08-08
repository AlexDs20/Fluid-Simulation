#include <algorithm>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <stdio.h>  // to remove files
#include <vector>

#include "boundary.h"
#include "compute.h"
#include "definitions.h"
#include "matrix.h"
#include "obstacle.h"
#include "parameters.h"

using namespace std;


int main(int , char *argv[]){

//-------------------------
// MPI STUFF
MPI_Init(NULL, NULL);
int nproc;
MPI_Comm_size(MPI_COMM_WORLD, &nproc);
int pid;
MPI_Comm_rank(MPI_COMM_WORLD, &pid);

//-------------------------
// Input variables
string inputfile  = argv[1];
string outputfile = argv[2];
string obsfile    = argv[3];
// Parameters
Parameters p(inputfile);
real delt;
real t = 0;
real t_print = 0;
// Set initial UVP arrays
matrix<real> U(p.imax+2,p.jmax+2,p.UI);
matrix<real> V(p.imax+2,p.jmax+2,p.VI);
matrix<real> P(p.imax+2,p.jmax+2,p.PI);
// Set F,G,RHS
matrix<real> F((p.imax+2),(p.jmax+2));
matrix<real> G((p.imax+2),(p.jmax+2));
matrix<real> RHS((p.imax+2),(p.jmax+2));
// Create the obstacles
obstacle Obs(p.imax+2,p.jmax+2,obsfile);

// Delete the existing output files
remove((outputfile+'U').c_str());
remove((outputfile+'V').c_str());
remove((outputfile+'P').c_str());

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

while (t<p.t_end){
  // Show runtime
  std::cout.flush();
  std::cout << "   " << t*100.0/p.t_end << "  %  \r";
  std::cout.flush();

  //-------------------------
  // Calculate time step
  delt = computeDT(p,U.amax(),V.amax());

  //-------------------------
  // Boundary Conditions
  setBoundaries(p,&U,&V);
  setSpecificBoundaries(p,&U,&V);
  setObsVelBoundaries(&Obs,&U,&V);

  //-------------------------
  //  Compute F and G
  computeF(p,delt,&Obs,&U,&V,&F);
  computeG(p,delt,&Obs,&U,&V,&G);

  //-------------------------
  //  Compute RHS
  computeRHS(p,delt,&Obs,&F,&G,&RHS);

  //-------------------------
  //  Push Pressure to
  //  the next time step
  computeP(p,&Obs,&RHS,&P);

  //-------------------------
  //  Compute new velocities
  computeNewVel(p,delt,&F,&G,&P,&U,&V);

  //-------------------------
  // Write output
  t_print += delt;
  if (t_print>= p.dt_out){
    writeOutput(outputfile,&U,&V,&P);
    t_print = 0;
  }

  t += delt;
}

// Split work between processors
if (pid == 1){

}

// Finalize the MPI environment.
MPI_Finalize();
return 0;
}
