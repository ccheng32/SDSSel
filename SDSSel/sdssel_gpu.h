#ifndef SDSSEL_GPU_H
#define SDSSEL_GPU_H
int sdssel_gpu_wrapper(int gpu_id, std::vector<int>& active_nodes,
                       std::vector<std::vector<int>>& h2v_group,
                       std::vector<std::vector<int>>& h2c_group,
                       std::vector<std::vector<int>>& channel_edge_dst,
                       std::vector<std::vector<double>>& channel_edge_wt,
                       int channel_offset);
#endif
