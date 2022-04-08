#pragma once
#include "myFunctions.h"

FinalResult* sequential(char* NS1, char** NS2, double* weights, int num_seqs);
void testResults(FinalResult* parallel_res, FinalResult* sequential_res, int num_seqs);
void testAndCompareTime(char** NS2, char* NS1, double* weights, int num_seqs, FinalResult* parallel_results);
