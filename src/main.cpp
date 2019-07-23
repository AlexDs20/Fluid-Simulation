#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using namespace std;

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

//-------------------------
// Set initial UVP arrays
Matrix U(p.imax+2,p.jmax+2,p.UI);
Matrix V(p.imax+2,p.jmax+2,p.VI);
Matrix P(p.imax+2,p.jmax+2,p.PI);

//-------------------------
// Create a boundary vector
/*
To be done: Create a vector that would contain values indicating whether it is a boundary or not.
The boundary values for the field should then be changed accordingly to that vector
*/

//-------------------------
// Calculate time step
Real delt=0;
delt = calcDT(p,U.max(),V.max());

//-------------------------
// Boundary Conditions
U.setUBoundCond(p.wW,p.wE,p.wS,p.wN);
V.setVBoundCond(p.wW,p.wE,p.wS,p.wN);

//-------------------------
// Write output
std::vector<Real> dataOut;
dataOut = U.getAll();

writeData(outputfile,dataOut);


// Split work between processors
if (pid == 1){

}

// Finalize the MPI environment.
MPI_Finalize();
return 0;
}
