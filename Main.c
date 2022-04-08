#include <stdio.h>
#include <stdlib.h>
#include "Main.h"
#include "/usr/include/x86_64-linux-gnu/mpich/mpi.h"
#include <time.h>
#include "myFunctions.h"
#include "mySeqSol.h"

#define MASTER 0
#define NUM_WEIGHTS 4

int main(int argc, char* argv[]) {
	int my_rank, number_of_processes;
	clock_t begin_parallel = clock(); // Get time for comparison

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Datatype MPI_RESULT = newType();
	MPI_Status status;
	// Get values from input file
	double* weights = (double*)malloc(NUM_WEIGHTS * sizeof(double));
	if (weights == NULL) {
		error_handler("Failed to allocate memory -> weights");
	}
	char* NS1 = (char*)malloc(3000 * sizeof(char));
	if (NS1 == NULL) {
		error_handler("Failed to allocate memory -> NS1");
	}
	int number_of_sequences = 0;
	char** NS2 = readFile(NS1, weights, &number_of_sequences);
	// Initialize parallel result array for all sequences results
	FinalResult* parallel_results = (FinalResult*)malloc(number_of_sequences * sizeof(FinalResult));
	if (parallel_results == NULL) {
		error_handler("Failed to allocate memory -> parallel results");
	}
	// For every SN2
	for (int seq = 0; seq < number_of_sequences; seq++) {
		FinalResult final_result = getBestScore(NS2[seq], NS1, weights, my_rank, number_of_processes);
		if (my_rank != MASTER) {
			// Slaves send final result - highest result for all mutants checked
			MPI_Send(&final_result, 1, MPI_RESULT, MASTER, 0, MPI_COMM_WORLD);
		}

		if (my_rank == MASTER) {
			// Master receive from all slaves their highest result
			for (int i = 1; i < number_of_processes; i++) {
				FinalResult result;
				MPI_Recv(&result, 1, MPI_RESULT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				// Update to final result if higher score
				final_result = compareFinalResults(result, final_result);
			}
			// Master gets final result for all calculations to current SN2 to parallel results
			parallel_results[seq] = final_result;
		}
		// All processes work on each SN2, so need to wait before continuing to next SN2
		MPI_Barrier(MPI_COMM_WORLD);
	}
	// After getting all SN2 results
	printf("\nafter SN2 calculations finished\n");
	if (my_rank == MASTER) {
		writeToFile(parallel_results, number_of_sequences);  // Master prints results to output file
		
		clock_t end_parallel = clock();
		double time_parallel = (double)(end_parallel - begin_parallel) / CLOCKS_PER_SEC;

		testAndCompareTime(NS2, NS1, weights, number_of_sequences, parallel_results);
		// Print parallel calculation execution time		
		printf("Total time for parallel calculation: %.4f minutes\n", time_parallel / 60);
	}

	free(NS1);
	free(NS2);
	free(weights);
	free(parallel_results);

	MPI_Finalize();

	return 0;
}

char** readFile(char* NS1, double* weights, int* number_of_sequences) {
	FILE* inputFile = fopen("input.txt", "r");

	if (inputFile == NULL) {
		error_handler("Failed to open the file");
	}
	fscanf(inputFile, "%lf %lf %lf %lf", &weights[0], &weights[1], &weights[2],
		&weights[3]);
	fscanf(inputFile, "%s", NS1);
	fscanf(inputFile, "%d", number_of_sequences);

	char** NS2 = (char**)malloc(sizeof(char*) * (*number_of_sequences));
	if (NS2 == NULL) {
		error_handler("Failed to allocate memory for NS2\n");
	}
	for (int i = 0; i < *number_of_sequences; i++) {
		NS2[i] = (char*)malloc(sizeof(char) * 2000);
		if (NS2[i] == NULL) {
			error_handler("Failed to allocate memory for NS2");
		}
		fscanf(inputFile, "%s", NS2[i]);
	}

	fclose(inputFile);
	return NS2;
}

void writeToFile(FinalResult* results, int number_of_sequences)
{
	FILE* outputFile;
	outputFile = fopen("output.txt", "w");

	if (outputFile == NULL) {
		error_handler("Failed to open the file");
	}

	for (int seq = 0; seq < number_of_sequences; seq++) {
		fprintf(outputFile, "n = %d\t k = %d \n", results[seq].offset, results[seq].mutant);
	}
	fclose(outputFile);
}

MPI_Datatype newType() {
	// Create new mpi data type for Result to transfer between processes
	MPI_Datatype MPI_RESULT;
	int lengths[3] = { 1, 1, 1 }; // none is array - each contains 1 value
	// Where each Result variable starts 
	// First one is score -  starts at 0
	// Second is offset - starts after score which is double
	// Last is mutant - starts after offset and score which are int and double
	const MPI_Aint displacements[3] = { 0, sizeof(double), sizeof(double) + sizeof(int) };
	// Result variables types	
	MPI_Datatype types[3] = { MPI_DOUBLE, MPI_INT, MPI_INT };
	MPI_Type_create_struct(3, lengths, displacements, types, &MPI_RESULT);
	MPI_Type_commit(&MPI_RESULT);

	return MPI_RESULT;
}
