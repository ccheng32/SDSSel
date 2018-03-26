#include <omp.h>
#include <cfloat>
#include <vector>
#include "macros.h"
#include "sdssel_gpu.h"
#define BLOCKDIM 128

__device__ int dst_in_h2c(int dst, int* h2c, int length) {
  for (int i = 0; i < length; i++) {
    if (dst == h2c[i]) {
      return i;
    }
  }
  return -1;
}

__global__ void sdssel_gpu(int* an, int num_an, int* h2c_data, int* h2c_pos,
                           int* h2c_length, int* c_dst_data, double* c_wt_data,
                           int* c_pos, int* c_length, int c_offset,
                           // output
                           int* best_node_ids, double* max_omegas) {
#define EMPTY 0xE321
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int node_id = -1;
  double max_omega = 0.0;

  __shared__ double max_omega_array[BLOCKDIM];
  __shared__ int best_user_array[BLOCKDIM];

  if (tid < num_an) {
    node_id = an[tid];
    int* h2c = h2c_data + h2c_pos[tid];
    int num_h2c = h2c_length[tid];
    double total_edge_sum = 0.0;
    double* h2c_edge_sum = new double[num_h2c];
    // calculate edge sum for all channels in h2c
    for (int i = 0; i < num_h2c; i++) {
      h2c_edge_sum[i] = 0.0;
      int chan_id = h2c[i] - c_offset;
      int* dsts = c_dst_data + c_pos[chan_id];
      double* wts = c_wt_data + c_pos[chan_id];

      // add the weight if its destination is in
      // h2c
      for (int j = 0; j < c_length[chan_id]; j++) {
        if (dst_in_h2c(dsts[j], h2c, num_h2c) != -1) {
          h2c_edge_sum[i] += wts[j];
          // add to total sum, same edge is counted twice
          total_edge_sum += wts[j];
        }
      }
    }

    // halve it because each edge is counted twice
    total_edge_sum /= 2;

    // pinch off the channel with smallest edge sum
    // in every iteration
    max_omega = total_edge_sum / num_h2c;
    for (int i = 0; i < num_h2c; i++) {
      double min_edge_sum = DBL_MAX;
      int min_edge_index = -1;

      // find channel with minimum edge sum
      for (int j = 0; j < num_h2c; j++) {
        if (h2c[j] != EMPTY) {
          if (h2c_edge_sum[j] < min_edge_sum) {
            min_edge_sum = h2c_edge_sum[j];
            min_edge_index = j;
          }
        }
      }

      // remove that channel
      int min_edge_chan_index = h2c[min_edge_index] - c_offset;
      int* min_chan_dsts = c_dst_data + c_pos[min_edge_chan_index];
      double* min_chan_wts = c_wt_data + c_pos[min_edge_chan_index];
      for (int j = 0; j < num_h2c; j++) {
        if (h2c[j] != EMPTY) {
          int this_idx;
          if ((this_idx = dst_in_h2c(h2c[j], min_chan_dsts,
                                     c_length[min_edge_chan_index])) != -1) {
            h2c_edge_sum[j] -= min_chan_wts[this_idx];
            total_edge_sum -= min_chan_wts[this_idx];
          }
        }
      }
      h2c[min_edge_index] = EMPTY;
      max_omega = max(max_omega, total_edge_sum / max(num_h2c - i - 1, 1));
    }
    delete[] h2c_edge_sum;
  }

  int tidx = threadIdx.x;
  max_omega_array[tidx] = max_omega;
  best_user_array[tidx] = node_id;
  __syncthreads();

  // start reduction to find the best center user
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tidx < s) {
      if (max_omega_array[tidx] < max_omega_array[tidx + s]) {
        max_omega_array[tidx] = max_omega_array[tidx + s];
        best_user_array[tidx] = best_user_array[tidx + s];
      }
    }
    __syncthreads();
  }

  if (tidx == 0) {
    max_omegas[blockIdx.x] = max_omega_array[0];
    best_node_ids[blockIdx.x] = best_user_array[0];
  }
}

/**
 * returns the best node to use as center
 **/
int sdssel_gpu_wrapper(int gpu_id, std::vector<int>& active_nodes,
                       std::vector<std::vector<int>>& h2v_group,
                       std::vector<std::vector<int>>& h2c_group,
                       std::vector<std::vector<int>>& channel_edge_dst,
                       std::vector<std::vector<double>>& channel_edge_wt,
                       int channel_offset) {
  // choose device
  cudaSetDevice(gpu_id);

  // copy active nodes onto gpu
  int* d_an;
  cudaMalloc(&d_an, sizeof(int) * active_nodes.size());
  cudaMemcpy(d_an, active_nodes.data(), sizeof(int) * active_nodes.size(),
             cudaMemcpyHostToDevice);

  // preprocess 2h dense set data
  std::vector<int> h2c_data;
  std::vector<int> h2c_start_poss;
  std::vector<int> h2c_lengths;
  int h2c_start_pos = 0, h2c_length = 0;

  for (int i = 0; i < h2v_group.size(); i++) {
    h2c_start_poss.PB(h2c_start_pos);
    h2c_length = 0;
    for (int h2c_mem : h2c_group[i]) {
      h2c_start_pos++;
      h2c_length++;
      h2c_data.PB(h2c_mem);
    }
    h2c_lengths.PB(h2c_length);
  }

  // copy 2h dense set data onto gpu
  int* d_h2c_data;
  cudaMalloc(&d_h2c_data, sizeof(int) * h2c_data.size());
  cudaMemcpy(d_h2c_data, h2c_data.data(), sizeof(int) * h2c_data.size(),
             cudaMemcpyHostToDevice);
  int* d_h2c_lengths;
  cudaMalloc(&d_h2c_lengths, sizeof(int) * h2c_lengths.size());
  cudaMemcpy(d_h2c_lengths, h2c_lengths.data(),
             sizeof(int) * h2c_lengths.size(), cudaMemcpyHostToDevice);
  int* d_h2c_start_poss;
  cudaMalloc(&d_h2c_start_poss, sizeof(int) * h2c_start_poss.size());
  cudaMemcpy(d_h2c_start_poss, h2c_start_poss.data(),
             sizeof(int) * h2c_start_poss.size(), cudaMemcpyHostToDevice);

  // pre process channel edges
  std::vector<int> channel_edge_dst_data;
  std::vector<double> channel_edge_wt_data;
  std::vector<int> channel_start_poss;
  std::vector<int> channel_lengths;

  int channel_start_pos = 0;
  int channel_length = 0;
  for (int i = 0; i < channel_edge_dst.size(); i++) {
    channel_start_poss.PB(channel_start_pos);
    channel_length = 0;
    for (int j = 0; j < channel_edge_dst[i].size(); j++) {
      channel_start_pos++;
      channel_length++;
      channel_edge_dst_data.PB(channel_edge_dst[i][j]);
      channel_edge_wt_data.PB(channel_edge_wt[i][j]);
    }
    channel_lengths.PB(channel_length);
  }

  // copy channel edges onto gpu
  int* d_channel_edge_dst_data;
  cudaMalloc(&d_channel_edge_dst_data,
             sizeof(int) * channel_edge_dst_data.size());
  cudaMemcpy(d_channel_edge_dst_data, channel_edge_dst_data.data(),
             sizeof(int) * channel_edge_dst_data.size(),
             cudaMemcpyHostToDevice);

  double* d_channel_edge_wt_data;
  cudaMalloc(&d_channel_edge_wt_data,
             sizeof(double) * channel_edge_wt_data.size());
  cudaMemcpy(d_channel_edge_wt_data, channel_edge_wt_data.data(),
             sizeof(double) * channel_edge_wt_data.size(),
             cudaMemcpyHostToDevice);

  int* d_channel_start_poss;
  cudaMalloc(&d_channel_start_poss, sizeof(int) * channel_start_poss.size());
  cudaMemcpy(d_channel_start_poss, channel_start_poss.data(),
             sizeof(int) * channel_start_poss.size(), cudaMemcpyHostToDevice);

  int* d_channel_lengths;
  cudaMalloc(&d_channel_lengths, sizeof(int) * channel_lengths.size());
  cudaMemcpy(d_channel_lengths, channel_lengths.data(),
             sizeof(int) * channel_lengths.size(), cudaMemcpyHostToDevice);

  int block_num = (active_nodes.size() - 1) / BLOCKDIM + 1;
  int* d_best_node_ids;
  double* d_max_omegas;
  cudaMalloc(&d_best_node_ids, sizeof(int) * block_num);
  cudaMalloc(&d_max_omegas, sizeof(double) * block_num);

  // launch kernel
  sdssel_gpu<<<block_num, BLOCKDIM>>>(
      d_an, active_nodes.size(), d_h2c_data, d_h2c_start_poss, d_h2c_lengths,
      d_channel_edge_dst_data, d_channel_edge_wt_data, d_channel_start_poss,
      d_channel_lengths, channel_offset, d_best_node_ids, d_max_omegas);

  // get output results
  int* best_node_ids = new int[block_num];
  double* max_omegas = new double[block_num];
  cudaMemcpy(best_node_ids, d_best_node_ids, sizeof(int) * block_num,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(max_omegas, d_max_omegas, sizeof(double) * block_num,
             cudaMemcpyDeviceToHost);

  double max_omega = 0.0;
  int best_node_id = -1;
  for (int i = 0; i < block_num; i++) {
    if (max_omega < max_omegas[i]) {
      max_omega = max_omegas[i];
      best_node_id = best_node_ids[i];
    }
  }

  delete[] best_node_ids;
  delete[] max_omegas;
  return best_node_id;
}
