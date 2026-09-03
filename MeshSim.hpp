
#pragma once
/*
  MeshSim.hpp
  Discrete-Event Simulation for Mesh Conveyor System (C++ translation)

  This is a translation of the Python implementation in MeshSim.py. fileciteturn0file1

  Notes:
  - The simulation uses a ring-buffer (head index) representation for both feeder belts
    (circular roll) and primary belts (shift with zero injection) to avoid per-step memmoves.
  - Buffer stats collection can be disabled (collect_buffer_stats=false) for performance.
*/

#include <cstdint>
#include <vector>
#include <string>
#include <optional>

namespace meshsim {

struct BufferBlocks {
  // Dimensions: [n_blocks][m][n][dirs]
  int m = 0;
  int n = 0;
  int dirs = 0;     // 1 for one-way, 2 for bidirectional
  int n_blocks = 0;

  // Flattened arrays in row-major order:
  // idx = (((b*m + f)*n + p)*dirs + d)
  std::vector<double> sum_blocks;   // per-block sum of buffer_levels over time-steps
  std::vector<int>    max_blocks;   // per-block max of buffer_levels
  std::vector<int>    hits_blocks;  // per-block count of time-steps when buffer achieved its (global) max within the block (same semantics as Python)

  bool empty() const { return n_blocks == 0 || m == 0 || n == 0 || dirs == 0; }
};

struct SimulationResult {
  int total_entered = 0;                        // total items that entered feeders (across all feeders)
  std::vector<std::vector<int>> feeder_primary_blocks; // per-block arrivals per feeder: [block][feeder]
  std::vector<int> block_counts;                // per-block total arrivals across all feeders
  std::vector<int> block_sizes;                 // per-block length in time steps

  std::optional<BufferBlocks> buffer_blocks;    // present if collect_buffer_stats=true
};

// buffer_capacity is interpreted like the Python code:
// buffer_capacity[0][f] = capacity at forward intersection of feeder f
// buffer_capacity[1][f] = capacity at backward intersection of feeder f (if bidirectional)
SimulationResult run_simulation(
    int NumberOfFeeders,
    int NumberOfPrimary,
    bool BiDirectional,
    int T,
    const std::vector<std::vector<int>>& buffer_capacity,
    int seed,
    bool Verbal,
    int block_size,
    bool collect_buffer_stats,
    std::string arrival_sequence_csv = "",
    int delta_f = 4,
    int delta_p = 2,
    int delta_L = 0);

// Helpers (translated from Python)
std::vector<double> autocorr(const std::vector<int>& blocks, int max_lag = 10);
int detect_warmup_from_blocks(const std::vector<int>& block_counts,
                             const std::vector<int>& block_sizes,
                             bool strict = true);

// Post-warmup aggregation (compatible with Python function aggregate_post_warmup)
struct PostWarmupAgg {
  int total_steps = 0;
  int total_entered = 0;
  double throughput = 0.0;   // items per time step
  double utilization = 0.0;  // throughput / n
  // avg_buffer/max_buffer/max_hits only meaningful if buffer stats were collected.
  // Flattened: [m][n][dirs]
  std::vector<double> avg_buffer;
  std::vector<int> max_buffer;
  std::vector<int> max_hits;
  std::vector<int> post_blocks; // block_counts[k:]
};

PostWarmupAgg aggregate_post_warmup(
    int k_warm,
    const std::vector<int>& block_sizes,
    const std::vector<int>& block_counts,
    const BufferBlocks& buffer_blocks,
    int NumberOfPrimary);

// Approximation model (translated from Python MeshSim.py)
struct ApproxResult {
  double utilization = 0.0;
  double throughput = 0.0;
  std::string case_name;
};

ApproxResult approximate_model(int NumberOfFeeders,
                              int NumberOfPrimary,
                              bool BiDirectional,
                              const std::vector<std::vector<int>>& buffer_capacity);

std::vector<double> approximate_model_per_feeder(int NumberOfFeeders,
                                                 int NumberOfPrimary,
                                                 bool BiDirectional,
                                                 const std::vector<std::vector<int>>& buffer_capacity);

} // namespace meshsim
