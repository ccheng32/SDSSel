#ifndef CUDA_QUEUE_CUH
#define CUDA_QUEUE_CUH

typedef struct {
  int* node_array;
  int* hop_array;
  int head_index;
  int tail_index;
  int total_length;
} cuda_queue;

__device__ __host__ void init_cuda_queue(cuda_queue* q, int size) {
  q->node_array = new int[size];
  q->hop_array = new int[size];
  q->tail_index = 1;
  q->head_index = 0;
  q->total_length = size;
}

__device__ __host__ void free_cuda_queue(cuda_queue* q) {
  delete[] q->node_array;
  delete[] q->hop_array;
}

__device__ __host__ bool isEmpty_cuda_queue(cuda_queue* q) {
  return q->tail_index == (q->head_index + 1) % q->total_length;
}

__device__ __host__ int peekNode_cuda_queue(cuda_queue* q) {
  return q->node_array[q->head_index];
}

__device__ __host__ int peekHop_cuda_queue(cuda_queue* q) {
  return q->hop_array[q->head_index];
}

__device__ __host__ void push_cuda_queue(cuda_queue* q, int node_id, int hops) {
  
  if (q->tail_index == q->head_index) {
    int length = q->total_length;
    int ti = q->tail_index;
    int hi = q->head_index;
    int empty_index = ((ti - 1) % length + length) % length;
    // array is full, need to double its size
    int* node_array_2x = new int[length << 1];
    int* hop_array_2x = new int[length << 1];
    for (int i = hi, j = 0; i != empty_index; i = (i + 1) % length, j++) {
      node_array_2x[j] = q->node_array[i];
      hop_array_2x[j] = q->hop_array[i];
    }
    
    delete[] q->node_array;
    delete[] q->hop_array;

    // update variables
    q->node_array = node_array_2x; 
    q->hop_array =  hop_array_2x;
    q->head_index = 0;
    q->tail_index = length;
    q->total_length = length << 1;
  }
  
  int length = q->total_length;
  int empty_index = ((q->tail_index - 1) % length + length) % length;
  q->node_array[empty_index] = node_id;
  q->hop_array[empty_index] = hops;
  q->tail_index = (q->tail_index + 1) % length;
}

__device__ __host__ void pop_cuda_queue(cuda_queue* q) {
  
  if (!isEmpty_cuda_queue(q)) {
    q->head_index = (q->head_index + 1) % q->total_length;
  }
}

#endif
