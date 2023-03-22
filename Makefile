build:
	mpicxx -fopenmp -c Main.c -o Main.o
	mpicxx -fopenmp -c myFunctions.c -o myFunctions.o
	mpicxx -fopenmp -c mySeqSol.c -o mySeqSol.o
	nvcc -gencode arch=compute_61,code=sm_61 -c CudaFunctions.cu -o CudaFunctions.o
	mpicxx -fopenmp -o SeqAlignment Main.o myFunctions.o mySeqSol.o CudaFunctions.o  -L/usr/local/cuda/lib -L/usr/local/cuda/lib64 -lcudart
run:
	mpiexec -np 4 ./SeqAlignment

run2Machines:
	mpiexec -np 15 --machinefile hosts.txt --map-by node ./SeqAlignment

clean:
	rm -f *.o SeqAlignment
	
