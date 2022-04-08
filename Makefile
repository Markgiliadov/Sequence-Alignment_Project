build:
	mpicxx -fopenmp -c Main.c -o Main.o
	mpicxx -fopenmp -c myFunctions.c -o myFunctions.o
	mpicxx -fopenmp -c mySeqSol.c -o mySeqSol.o
	nvcc -I./inc -c CudaFunctions.cu -o CudaFunctions.o
	mpicxx -fopenmp -o SeqAlignment Main.o myFunctions.o mySeqSol.o CudaFunctions.o  /usr/local/cuda-11.0/lib64/libcudart_static.a -ldl -lrt

run:
	mpiexec -np 6 ./SeqAlignment

run2Machines:
	mpiexec -np 15 --machinefile hosts.txt --map-by node ./SeqAlignment

clean:
	rm -f *.o SeqAlignment
	