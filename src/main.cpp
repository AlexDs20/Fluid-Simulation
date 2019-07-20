#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
using namespace std;

#include "Parameters.h"
#define Real float

int main(){
//-------------------------
// MPI STUFF
MPI_Init(NULL, NULL);
int nproc;
MPI_Comm_size(MPI_COMM_WORLD, &nproc);
int pid;
MPI_Comm_rank(MPI_COMM_WORLD, &pid);

//-------------------------
// Parameter file
string inputfile="input";


if (pid == 1){
// Read input parameters in the structure P
Parameters p(inputfile);
cout << p.Re << endl;

}

// Finalize the MPI environment.
MPI_Finalize();
return 0;
}

