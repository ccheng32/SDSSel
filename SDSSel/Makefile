DEBUG=n
STDFLAG=--std=c++11
CUDAPATH=/usr/local/cuda
CUDALIB=-L$(CUDAPATH)/lib64/ -lcuda -lcudart
NVCC=$(CUDAPATH)/bin/nvcc
CUDACC=sm_61

ifeq ($(DEBUG),n)
	OPTFLAG= -O3
else
	OPTFLAG= -O0 -g -DDEBUG
	CUDADEBUGFLAG = -G
endif

sdssel: SDSSel.cpp sdssel_gpu.h all_pairs_bfs.h macros.h sdssel_gpu.o all_pairs_bfs.o
	g++ -fopenmp $(STDFLAG) $(OPTFLAG) -Wall sdssel_gpu.o all_pairs_bfs.o SDSSel.cpp $(CUDALIB) -o sdssel

sdssel_gpu.o: sdssel_gpu.cu sdssel_gpu.h
	$(NVCC) --default-stream per-thread $(CUDADEBUGFLAG) -Xcompiler -fopenmp $(OPTFLAG) $(STDFLAG) --gpu-architecture=$(CUDACC) -c sdssel_gpu.cu

all_pairs_bfs.o: all_pairs_bfs.cu all_pairs_bfs.h cuda_queue/cuda_queue.cuh
	$(NVCC) --default-stream per-thread $(CUDADEBUGFLAG) -Xcompiler -fopenmp $(OPTFLAG) $(STDFLAG) --gpu-architecture=$(CUDACC) -c all_pairs_bfs.cu -o all_pairs_bfs.o


clean:
	rm -f sdssel *.o
