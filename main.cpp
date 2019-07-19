#include <iostream>
#include <fstream>
#include <mpi.h>

#define REAL float

using namespace std;

int main(){
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Get the rank of the process
  int pid;
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  if (pid == 1){
  // Print off a hello world message
  printf("Hello world, this is rank %d out of %d processors\n",
  pid, nproc);
  }

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
