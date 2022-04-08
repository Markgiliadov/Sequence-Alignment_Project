#include "/usr/include/x86_64-linux-gnu/mpich/mpi.h"
#include "mySeqSol.h"

char** readFile(char* NS1, double* weights, int* number_of_sequences);
void writeToFile(FinalResult* results, int number_of_sequences);
MPI_Datatype newType();