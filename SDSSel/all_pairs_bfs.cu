#include <cstdint>
#include "all_pairs_bfs.h"
#include "cuda_queue/cuda_queue.cuh"
#include "macros.h"
#define BLOCKDIM 128

__global__ void all_pairs_bfs_cuda(int* d_an, uint64_t* d_sp, uint64_t* d_l,
                                   int* d_ed, int h, size_t num_an) {
  size_t id = threadIdx.x + blockIdx.x * blockDim.x;
  if (id < num_an) {
    int node_id = d_an[id];
    // start BFS for each node here

    cuda_queue q;
    init_cuda_queue(&q, 16);
    push_cuda_queue(&q, node_id, 0);
    free_cuda_queue(&q);
  }
}

void all_pairs_bfs(std::vector<int>& active_nodes,
                   std::vector<std::vector<std::pair<int, double>>>& edges,
                   int v, int c, int h) {
  // copy active nodes onto gpu
  int* d_active_nodes;
  cudaMalloc(&d_active_nodes, sizeof(int) * active_nodes.size());
  cudaMemcpy(d_active_nodes, active_nodes.data(),
             sizeof(int) * active_nodes.size(), cudaMemcpyHostToDevice);

  // preprocess the edges
  // the original one-indexing of nodes has been converted to zero indexing
  std::vector<uint64_t> start_positions(v, 0);
  std::vector<uint64_t> lengths(v, 0);
  std::vector<int> edge_data;
  uint64_t start_position = 0;
  uint64_t length = 0;
  for (size_t i = 1; i <= v; i++) {
    start_positions[i - 1] = start_position;
    length = 0;
    for (auto neighbor_pair : edges[i]) {
      if (neighbor_pair.F <= v) {
        start_position++;
        length++;
        edge_data.PB(neighbor_pair.F - 1);
      }
    }
    lengths[i - 1] = length;
  }

  // copy edges onto gpu
  int* d_edge_data;
  cudaMalloc(&d_edge_data, sizeof(int) * edge_data.size());
  cudaMemcpy(d_edge_data, edge_data.data(), sizeof(int) * edge_data.size(),
             cudaMemcpyHostToDevice);
  uint64_t* d_lengths;
  cudaMalloc(&d_lengths, sizeof(uint64_t) * lengths.size());
  cudaMemcpy(d_lengths, lengths.data(), sizeof(uint64_t) * lengths.size(),
             cudaMemcpyHostToDevice);
  uint64_t* d_start_positions;
  cudaMalloc(&d_start_positions, sizeof(uint64_t) * start_positions.size());
  cudaMemcpy(d_start_positions, start_positions.data(),
             sizeof(uint64_t) * start_positions.size(), cudaMemcpyHostToDevice);

  // launch kernel
  int blocknum = (active_nodes.size() - 1) / BLOCKDIM + 1;
  all_pairs_bfs_cuda<<<blocknum, BLOCKDIM>>>(d_active_nodes, d_start_positions,
                                             d_lengths, d_edge_data, h,
                                             active_nodes.size());

  cudaFree(d_active_nodes);
  cudaFree(d_edge_data);
  cudaFree(d_lengths);
  cudaFree(d_start_positions);
}
