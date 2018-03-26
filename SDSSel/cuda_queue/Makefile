DEBUG=y
STDFLAG=--std=c++11
CUDALIB=-L/usr/local/cuda/lib64/ -lcuda -lcudart

ifeq ($(DEBUG),n)
	OPTFLAG= -O3
else
	OPTFLAG= -O0 -g -DDEBUG
endif

test_cuda_queue: test_cuda_queue.cu cuda_queue.cuh
	nvcc --default-stream per-thread -Xcompiler $(OPTFLAG) $(STDFLAG) --gpu-architecture=sm_61 test_cuda_queue.cu -o test_cuda_queue

clean:
	rm -f test_cuda_queue *.o
