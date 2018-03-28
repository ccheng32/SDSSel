#include "cuda_queue.cuh"
#include <cstdio>

int main() {

  cuda_queue q;
  init_cuda_queue(&q, 5);

  for (int i  = 0; i < 100; i++) {
    push_cuda_queue(&q, i, i+1);
  }

  if (peekNode_cuda_queue(&q) != 0 || peekHop_cuda_queue(&q) != 1) {
    printf("fail %d\n", __LINE__);
    return 1;
  }


  free_cuda_queue(&q);

  return 0;
}
