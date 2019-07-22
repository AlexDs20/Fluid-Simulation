#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
using namespace std;

#include "Parameters.h"
#include "Matrix.h"
#define Real float

int main(int argc, char *argv[]){
//-------------------------
// MPI STUFF
MPI_Init(NULL, NULL);
int nproc;
MPI_Comm_size(MPI_COMM_WORLD, &nproc);
int pid;
MPI_Comm_rank(MPI_COMM_WORLD, &pid);

//-------------------------
// Parameter file
//string inputfile="input";
string inputfile=argv[argc-1];
// Read input parameters in the structure P
Parameters p(inputfile);
// Allocate the memory for UVP arrays
Matrix U(p.imax+2,p.jmax+2);
U.setVal(p.UI);

// Split work between processors
if (pid == 1){

}

// Finalize the MPI environment.
MPI_Finalize();
return 0;
}

