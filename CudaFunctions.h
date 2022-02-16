__global__ void calcScoresKernel(double* dev_results, double* dev_scores, int num_rows, int num_cols);
__host__ void checkErrors(cudaError_t err, const char* error_msg);