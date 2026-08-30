#include "MeshSim.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <fstream>
#include <sstream>

namespace meshsim {

namespace {

// Safe access to buffer_capacity[dir][f]; returns 0 if missing.
inline int cap_at(const std::vector<std::vector<int>>& cap, int dir, int f) {
  if (dir < 0) return 0;
  if (dir >= static_cast<int>(cap.size())) return 0;
  if (f < 0) return 0;
  if (f >= static_cast<int>(cap[dir].size())) return 0;
  return std::max(0, cap[dir][f]);
}

struct Junction {
  int f;
  int fp;
  int p;
  int pp;
  int dir;
};

inline int mod_pos(int x, int L) {
  // L must be > 0
  x %= L;
  if (x < 0) x += L;
  return x;
}

} // namespace

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
    std::string arrival_sequence_csv) {

  if (NumberOfFeeders <= 0 || NumberOfPrimary <= 0) {
    throw std::invalid_argument("NumberOfFeeders and NumberOfPrimary must be positive.");
  }
  if (T < 0) throw std::invalid_argument("T must be >= 0");
  if (block_size <= 0) throw std::invalid_argument("block_size must be positive.");

  const int m = NumberOfFeeders;
  const int n = NumberOfPrimary;
  const int PrimaryLength = 4 * m;
  const int FeederLength  = 4 * n + 8;
  const int dirs = BiDirectional ? 2 : 1;

  std::vector<int> empirical_sequence;
  int seq_idx = 0;
  if (!arrival_sequence_csv.empty()) {
    std::ifstream fs(arrival_sequence_csv);
    if (!fs.is_open()) {
        throw std::runtime_error("Cannot open arrival sequence file: " + arrival_sequence_csv);
    }
    std::string line;
    std::getline(fs, line); // header
    while (std::getline(fs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        int col = 0;
        int belt_val = -1;
        while (std::getline(ss, token, ',')) {
            if (col == 8) { // belt is 9th column
                belt_val = std::stoi(token);
                break;
            }
            col++;
        }
        if (belt_val != -1) {
            empirical_sequence.push_back(belt_val);
        }
    }
    if (empirical_sequence.empty()) {
        throw std::runtime_error("No valid arrivals found in " + arrival_sequence_csv);
    }
  }

  std::mt19937 rng(static_cast<uint32_t>(seed));
  std::uniform_int_distribution<int> dest_dist(1, n);

  // Belts represented as ring buffers via "head" offsets.
  std::vector<std::vector<uint8_t>> primary_belts(n, std::vector<uint8_t>(PrimaryLength, 0));
  std::vector<int> primary_head(n, 0);

  std::vector<std::vector<int>> feeder_belts(m, std::vector<int>(FeederLength, 0));
  std::vector<int> feeder_head(m, 0);

  // buffer_levels: [m][n][dirs] flattened
  const int levels_dim = m * n * dirs;
  std::vector<int> buffer_levels(levels_dim, 0);

  auto idx_level = [&](int f, int p, int d) -> int {
    return ((f * n + p) * dirs + d);
  };

  // Block stats
  int total_entered = 0;
  std::vector<int> block_counts;
  std::vector<int> block_sizes;
  block_counts.reserve((T / block_size) + 2);
  block_sizes.reserve((T / block_size) + 2);

  std::vector<std::vector<int>> feeder_primary_blocks;
  feeder_primary_blocks.reserve((T / block_size) + 2);
  std::vector<int> feeder_primary_curblock(m, 0);

  std::vector<double> cur_sum;
  std::vector<int> cur_max;
  std::vector<int> cur_hits;
  std::vector<double> sum_blocks;
  std::vector<int> max_blocks;
  std::vector<int> hits_blocks;

  if (collect_buffer_stats) {
    cur_sum.assign(levels_dim, 0.0);
    cur_max.assign(levels_dim, 0);
    cur_hits.assign(levels_dim, 0);
    sum_blocks.reserve(static_cast<size_t>(levels_dim) * ((T / block_size) + 2));
    max_blocks.reserve(static_cast<size_t>(levels_dim) * ((T / block_size) + 2));
    hits_blocks.reserve(static_cast<size_t>(levels_dim) * ((T / block_size) + 2));
  }

  // Pre-build junction list
  std::vector<Junction> junction;
  junction.reserve(static_cast<size_t>(m) * static_cast<size_t>(n) * (BiDirectional ? 2 : 1));

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      junction.push_back(Junction{i, 1 + 2 * j, j, 4 * i, 0});
      if (BiDirectional) {
        junction.push_back(Junction{i, FeederLength - (3 + 2 * j), j, 2 + 4 * i, 1});
      }
    }
  }

  auto feeder_at = [&](int f, int logical_pos) -> int& {
    const int phys = (feeder_head[f] + logical_pos) % FeederLength;
    return feeder_belts[f][phys];
  };
  auto primary_at = [&](int p, int logical_pos) -> uint8_t& {
    const int phys = (primary_head[p] + logical_pos) % PrimaryLength;
    return primary_belts[p][phys];
  };

  int cur_size = 0;
  int inblock = 0;

  for (int t = 0; t < T; ++t) {

    // Arrivals at feeder entrances
    for (int f = 0; f < m; ++f) {
      int& cell0 = feeder_at(f, 0);
      if (cell0 == 0) {
        if (!empirical_sequence.empty()) {
          cell0 = empirical_sequence[seq_idx];
          seq_idx = (seq_idx + 1) % empirical_sequence.size();
        } else {
          cell0 = dest_dist(rng);
        }
        total_entered += 1;
        inblock += 1;
        feeder_primary_curblock[f] += 1;
      } else if (Verbal) {
        std::cout << "Feeder " << f << " blocked at entrance, time " << t << "\n";
      }
    }

    // Transfers at junctions
    for (const auto& j : junction) {
      const int f = j.f;
      const int fp = j.fp;
      const int p = j.p;
      const int pp = j.pp;
      const int dir = j.dir;

      uint8_t& pcell = primary_at(p, pp);
      const int id = idx_level(f, p, dir);

      // Buffer to primary
      if (pcell == 0 && buffer_levels[id] > 0) {
        pcell = 1;
        buffer_levels[id] -= 1;
      }

      // From feeder to primary or buffer
      int& fcell = feeder_at(f, fp);
      if (fcell == (p + 1)) {
        if (pcell == 0) {
          fcell = 0;
          pcell = 1;
        } else {
          const int cap = cap_at(buffer_capacity, dir, f);
          if (cap > 0 && buffer_levels[id] < cap) {
            buffer_levels[id] += 1;
            fcell = 0;
          } else if (Verbal) {
            std::cout << "Feeder " << f << "->Primary " << p << " (dir=" << dir
                      << ") blocked at t=" << t << " (buffer=" << buffer_levels[id] << ")\n";
          }
        }
      }
    }

    // Advance primary belts: shift with zero injection
    for (int p = 0; p < n; ++p) {
      primary_head[p] = mod_pos(primary_head[p] - 1, PrimaryLength);
      primary_belts[p][primary_head[p]] = 0;
    }

    // Advance feeder belts: circular roll
    for (int f = 0; f < m; ++f) {
      feeder_head[f] = mod_pos(feeder_head[f] - 1, FeederLength);
    }

    // Block collection
    cur_size += 1;
    if (collect_buffer_stats) {
      for (int id = 0; id < levels_dim; ++id) {
        const int curr = buffer_levels[id];
        cur_sum[id] += static_cast<double>(curr);

        if (curr > cur_max[id]) {
          cur_max[id] = curr;
          cur_hits[id] = (curr > 0) ? 1 : 0;
        } else if (curr == cur_max[id] && cur_max[id] > 0) {
          cur_hits[id] += 1;
        }
      }
    }

    if (cur_size == block_size) {
      block_counts.push_back(inblock);
      block_sizes.push_back(cur_size);
      feeder_primary_blocks.push_back(feeder_primary_curblock);
      std::fill(feeder_primary_curblock.begin(), feeder_primary_curblock.end(), 0);

      if (collect_buffer_stats) {
        sum_blocks.insert(sum_blocks.end(), cur_sum.begin(), cur_sum.end());
        max_blocks.insert(max_blocks.end(), cur_max.begin(), cur_max.end());
        hits_blocks.insert(hits_blocks.end(), cur_hits.begin(), cur_hits.end());
        std::fill(cur_sum.begin(), cur_sum.end(), 0.0);
        std::fill(cur_max.begin(), cur_max.end(), 0);
        std::fill(cur_hits.begin(), cur_hits.end(), 0);
      }

      inblock = 0;
      cur_size = 0;
    }
  }

  // Tail block (if T not multiple of block_size)
  if (cur_size > 0) {
    block_counts.push_back(inblock);
    block_sizes.push_back(cur_size);
    feeder_primary_blocks.push_back(feeder_primary_curblock);
    // no need to reset

    if (collect_buffer_stats) {
      sum_blocks.insert(sum_blocks.end(), cur_sum.begin(), cur_sum.end());
      max_blocks.insert(max_blocks.end(), cur_max.begin(), cur_max.end());
      hits_blocks.insert(hits_blocks.end(), cur_hits.begin(), cur_hits.end());
    }
  }

  SimulationResult res;
  res.total_entered = total_entered;
  res.feeder_primary_blocks = std::move(feeder_primary_blocks);
  res.block_counts = std::move(block_counts);
  res.block_sizes = std::move(block_sizes);

  if (collect_buffer_stats) {
    BufferBlocks bb;
    bb.m = m;
    bb.n = n;
    bb.dirs = dirs;
    bb.n_blocks = static_cast<int>(res.block_counts.size());
    bb.sum_blocks = std::move(sum_blocks);
    bb.max_blocks = std::move(max_blocks);
    bb.hits_blocks = std::move(hits_blocks);
    res.buffer_blocks = std::move(bb);
  }

  return res;
}

std::vector<double> autocorr(const std::vector<int>& blocks, int max_lag) {
  std::vector<double> out;
  const int N = static_cast<int>(blocks.size());
  if (N < 2) return out;

  std::vector<double> x(N);
  double mean = 0.0;
  for (int i = 0; i < N; ++i) mean += static_cast<double>(blocks[i]);
  mean /= static_cast<double>(N);
  for (int i = 0; i < N; ++i) x[i] = static_cast<double>(blocks[i]) - mean;

  double denom = 0.0;
  for (double v : x) denom += v * v;
  if (denom == 0.0) {
    out.assign(std::min(max_lag, N - 1), 0.0);
    return out;
  }

  const int L = std::min(max_lag, N - 1);
  out.reserve(L);
  for (int k = 1; k <= L; ++k) {
    double num = 0.0;
    for (int i = 0; i < N - k; ++i) num += x[i] * x[i + k];
    out.push_back(num / denom);
  }
  return out;
}

int detect_warmup_from_blocks(const std::vector<int>& block_counts,
                             const std::vector<int>& block_sizes,
                             bool strict) {
  const int n = static_cast<int>(block_counts.size());
  if (n < 2) return 0;

  std::vector<double> rates(n, 0.0);
  for (int i = 0; i < n; ++i) {
    const int s = (i < static_cast<int>(block_sizes.size())) ? block_sizes[i] : 0;
    const int c = block_counts[i];
    rates[i] = (s > 0) ? (static_cast<double>(c) / static_cast<double>(s)) : 0.0;
  }

  std::vector<double> suf_sum(n + 1, 0.0);
  for (int i = n - 1; i >= 0; --i) {
    suf_sum[i] = suf_sum[i + 1] + rates[i];
  }

  for (int k = 0; k < n - 1; ++k) {
    const int after_len = n - (k + 1);
    const double after_avg = suf_sum[k + 1] / static_cast<double>(after_len);
    if (strict) {
      if (rates[k] < after_avg) return k;
    } else {
      if (rates[k] <= after_avg) return k;
    }
  }
  return 0;
}

PostWarmupAgg aggregate_post_warmup(
    int k_warm,
    const std::vector<int>& block_sizes,
    const std::vector<int>& block_counts,
    const BufferBlocks& buffer_blocks,
    int NumberOfPrimary) {

  PostWarmupAgg out;
  const int n_blocks = static_cast<int>(block_counts.size());
  if (k_warm < 0) k_warm = 0;
  if (k_warm > n_blocks) k_warm = n_blocks;

  const int m = buffer_blocks.m;
  const int n = buffer_blocks.n;
  const int dirs = buffer_blocks.dirs;
  const int dims = m * n * dirs;

  int total_steps = 0;
  int total_entered = 0;
  for (int b = k_warm; b < n_blocks; ++b) {
    const int bs = (b < static_cast<int>(block_sizes.size())) ? block_sizes[b] : 0;
    total_steps += bs;
    total_entered += block_counts[b];
  }

  out.total_steps = total_steps;
  out.total_entered = total_entered;
  out.throughput = (total_steps > 0) ? (static_cast<double>(total_entered) / static_cast<double>(total_steps)) : 0.0;
  out.utilization = (total_steps > 0 && NumberOfPrimary > 0)
      ? (static_cast<double>(total_entered) / (static_cast<double>(total_steps) * static_cast<double>(NumberOfPrimary)))
      : 0.0;

  out.post_blocks.assign(block_counts.begin() + k_warm, block_counts.end());

  // If buffer data isn't present or dims mismatch, return without buffer stats.
  if (buffer_blocks.sum_blocks.size() < static_cast<size_t>(dims) * static_cast<size_t>(n_blocks) ||
      buffer_blocks.max_blocks.size() < static_cast<size_t>(dims) * static_cast<size_t>(n_blocks) ||
      buffer_blocks.hits_blocks.size() < static_cast<size_t>(dims) * static_cast<size_t>(n_blocks)) {
    return out;
  }

  out.avg_buffer.assign(dims, 0.0);
  out.max_buffer.assign(dims, 0);
  out.max_hits.assign(dims, 0);

  auto base = [&](int b) -> int { return b * dims; };

  // avg_buffer: sum of per-block sums, divided by total_steps
  if (total_steps > 0) {
    for (int b = k_warm; b < n_blocks; ++b) {
      const int off = base(b);
      for (int i = 0; i < dims; ++i) {
        out.avg_buffer[i] += buffer_blocks.sum_blocks[off + i];
      }
    }
    for (int i = 0; i < dims; ++i) out.avg_buffer[i] /= static_cast<double>(total_steps);
  }

  // max_buffer: max across blocks
  for (int b = k_warm; b < n_blocks; ++b) {
    const int off = base(b);
    for (int i = 0; i < dims; ++i) {
      out.max_buffer[i] = std::max(out.max_buffer[i], buffer_blocks.max_blocks[off + i]);
    }
  }

  // max_hits: sum hits over blocks where block_max equals global max (>0)
  for (int b = k_warm; b < n_blocks; ++b) {
    const int off = base(b);
    for (int i = 0; i < dims; ++i) {
      const int gmax = out.max_buffer[i];
      if (gmax > 0 && buffer_blocks.max_blocks[off + i] == gmax) {
        out.max_hits[i] += buffer_blocks.hits_blocks[off + i];
      }
    }
  }

  return out;
}

// ---------------- Approximation model (pure math) ----------------

static void extract_buffer_caps(int m,
                                const std::vector<std::vector<int>>& buffer_capacity,
                                std::vector<int>& c_f,
                                std::vector<int>& c_b) {
  c_f.assign(m, 0);
  c_b.assign(m, 0);
  if (m <= 0) return;

  if (!buffer_capacity.empty()) {
    const auto& forward = buffer_capacity[0];
    for (int j = 0; j < m && j < static_cast<int>(forward.size()); ++j) {
      c_f[j] = std::max(0, forward[j]);
    }
  }
  if (buffer_capacity.size() >= 2) {
    const auto& backward = buffer_capacity[1];
    for (int j = 0; j < m && j < static_cast<int>(backward.size()); ++j) {
      c_b[j] = std::max(0, backward[j]);
    }
  }

  // theory: first feeder is never blocked -> no buffers there
  c_f[0] = 0;
  c_b[0] = 0;
}

static std::pair<double,double> geo_geo_pi0_pi_c(double rho, int c) {
  if (c <= 0 || rho <= 0.0) return {1.0, 0.0};
  if (std::abs(rho - 1.0) < 1e-12) {
    const double pi0 = 1.0 / (static_cast<double>(c) + 1.0);
    return {pi0, pi0};
  }
  double rho_pow_c1 = 0.0;
  try {
    rho_pow_c1 = std::pow(rho, static_cast<double>(c + 1));
  } catch (...) {
    return {0.0, 1.0};
  }
  const double denom = 1.0 - rho_pow_c1;
  if (std::abs(denom) < 1e-18) {
    if (rho > 1.0) return {0.0, 1.0};
    return {1.0, 0.0};
  }
  double pi0 = (1.0 - rho) / denom;
  double pic = (1.0 - rho) * std::pow(rho, static_cast<double>(c)) / denom;
  pi0 = std::max(0.0, std::min(1.0, pi0));
  pic = std::max(0.0, std::min(1.0, pic));
  return {pi0, pic};
}

static double approx_one_way_no_buffers(int n, int m) {
  if (n <= 0 || m <= 0) return 0.0;
  return 1.0 - std::pow(1.0 - 1.0 / static_cast<double>(n), static_cast<double>(m));
}

static double approx_one_way_forward_buffers(int n, int m, const std::vector<int>& c_forward) {
  if (n <= 0 || m <= 0) return 0.0;
  double u_prev = 0.0;
  for (int j = 1; j <= m; ++j) {
    const int c_j = (j - 1 < static_cast<int>(c_forward.size())) ? c_forward[j - 1] : 0;
    double pi0 = 1.0;
    if (c_j > 0) {
      const double denom = (static_cast<double>(n) - 1.0) * (1.0 - u_prev);
      const double rho = (denom <= 0.0) ? 1e9 : (u_prev / denom);
      pi0 = geo_geo_pi0_pi_c(rho, c_j).first;
    }
    const double u = u_prev + (1.0 - u_prev) * ((1.0 - pi0) + pi0 / static_cast<double>(n));
    u_prev = u;
  }
  return u_prev;
}

static std::vector<double> approx_one_way_no_buffers_trace(int n, int m) {
  std::vector<double> out;
  if (n <= 0 || m <= 0) return out;
  out.reserve(m);
  double u = 0.0;
  for (int j = 1; j <= m; ++j) {
    u = u + (1.0 - u) * (1.0 / static_cast<double>(n));
    out.push_back(u);
  }
  return out;
}

static std::vector<double> approx_one_way_forward_buffers_trace(int n, int m, const std::vector<int>& c_forward) {
  std::vector<double> out;
  if (n <= 0 || m <= 0) return out;
  out.reserve(m);
  double u_prev = 0.0;
  for (int j = 1; j <= m; ++j) {
    const int c_j = (j - 1 < static_cast<int>(c_forward.size())) ? c_forward[j - 1] : 0;
    double pi0 = 1.0;
    if (c_j > 0) {
      const double denom = (static_cast<double>(n) - 1.0) * (1.0 - u_prev);
      const double rho = (denom <= 0.0) ? 1e9 : (u_prev / denom);
      pi0 = geo_geo_pi0_pi_c(rho, c_j).first;
    }
    const double u = u_prev + (1.0 - u_prev) * ((1.0 - pi0) + pi0 / static_cast<double>(n));
    u_prev = u;
    out.push_back(u);
  }
  return out;
}

static std::pair<std::vector<double>, std::vector<double>> approx_two_way_per_feeder_buffers_trace(
    int n, int m, const std::vector<int>& c_forward, const std::vector<int>& c_backward) {

  std::vector<double> u_f;
  std::vector<double> u_b;
  if (n <= 0 || m <= 0) return {u_f, u_b};

  u_f.assign(m + 1, 0.0);
  u_b.assign(m + 1, 0.0);

  u_f[1] = 1.0 / static_cast<double>(n);
  u_b[1] = 1.0 / static_cast<double>(n);

  for (int j = 2; j <= m; ++j) {
    const int c_fj = (j - 1 < static_cast<int>(c_forward.size())) ? c_forward[j - 1] : 0;
    const int c_bj = (j - 1 < static_cast<int>(c_backward.size())) ? c_backward[j - 1] : 0;

    const double u_prev_b = u_b[j - 1];

    // Forward
    double pi0_f = 1.0;
    double pi_c_f = 0.0;
    double P_jf = 1.0 / static_cast<double>(n);

    if (c_fj > 0) {
      const double denom = (static_cast<double>(n) - 1.0) * (1.0 - u_prev_b);
      const double rho_f = (denom <= 0.0) ? 1e9 : (u_prev_b / denom);
      auto [pi0, pic] = geo_geo_pi0_pi_c(rho_f, c_fj);
      pi0_f = pi0;
      pi_c_f = pic;
      P_jf = (1.0 - pi0_f) + (pi0_f / static_cast<double>(n));
    } else {
      P_jf = 1.0 / static_cast<double>(n);
    }

    u_f[j] = u_prev_b + (1.0 - u_prev_b) * P_jf;

    // Backward
    double P_jb = 0.0;
    if (c_fj <= 0 && c_bj <= 0) {
      P_jb = (u_prev_b / static_cast<double>(n));
    } else if (c_fj > 0 && c_bj <= 0) {
      P_jb = (u_prev_b * pi_c_f) / static_cast<double>(n);
    } else {
      if (c_bj <= 0) {
        P_jb = (u_prev_b / static_cast<double>(n));
      } else {
        double b_jb = 0.0;
        double d_jb = 0.0;
        if (c_fj > 0) {
          b_jb = (pi_c_f * u_prev_b * u_f[j]) / static_cast<double>(n);
          d_jb = (1.0 - u_f[j]) * (1.0 - (u_prev_b * pi_c_f) / static_cast<double>(n));
        } else {
          b_jb = (u_f[j] * u_prev_b) / static_cast<double>(n);
          d_jb = (1.0 - u_f[j]) * ((1.0 - 1.0 / static_cast<double>(n)) + (1.0 / static_cast<double>(n)) * (1.0 - u_prev_b));
        }
        const double rho_b = (d_jb <= 0.0) ? 1e9 : (b_jb / d_jb);
        const double pi0_b = geo_geo_pi0_pi_c(rho_b, c_bj).first;
        P_jb = (1.0 - pi0_b) + pi0_b * (u_prev_b / static_cast<double>(n));
      }
    }

    u_b[j] = u_f[j] + (1.0 - u_f[j]) * P_jb;
  }

  std::vector<double> u_f_list;
  std::vector<double> u_b_list;
  u_f_list.reserve(m);
  u_b_list.reserve(m);
  for (int j = 1; j <= m; ++j) {
    u_f_list.push_back(u_f[j]);
    u_b_list.push_back(u_b[j]);
  }
  return {u_f_list, u_b_list};
}

static double approx_two_way_per_feeder_buffers(int n, int m, const std::vector<int>& c_forward, const std::vector<int>& c_backward) {
  auto [u_f_list, u_b_list] = approx_two_way_per_feeder_buffers_trace(n, m, c_forward, c_backward);
  if (u_b_list.empty()) return 0.0;
  return u_b_list.back();
}

ApproxResult approximate_model(int NumberOfFeeders,
                              int NumberOfPrimary,
                              bool BiDirectional,
                              const std::vector<std::vector<int>>& buffer_capacity) {

  const int m = NumberOfFeeders;
  const int n = NumberOfPrimary;

  std::vector<int> c_f, c_b;
  extract_buffer_caps(m, buffer_capacity, c_f, c_b);

  const bool has_f = std::any_of(c_f.begin(), c_f.end(), [](int x){ return x > 0; });
  const bool has_b = std::any_of(c_b.begin(), c_b.end(), [](int x){ return x > 0; });

  double u = 0.0;
  std::string case_name;

  if (!BiDirectional) {
    if (has_f) {
      u = approx_one_way_forward_buffers(n, m, c_f);
      case_name = "one_way_forward_buffers";
    } else {
      u = approx_one_way_no_buffers(n, m);
      case_name = "one_way_no_buffers";
    }
  } else {
    u = approx_two_way_per_feeder_buffers(n, m, c_f, c_b);
    if (!has_f && !has_b) case_name = "two_way_no_buffers";
    else if (has_f && !has_b) case_name = "two_way_forward_buffers_only";
    else if (!has_f && has_b) case_name = "two_way_backward_buffers_only";
    else case_name = "two_way_forward_and_backward_buffers";
  }

  u = std::max(0.0, std::min(1.0, u));
  ApproxResult out;
  out.utilization = u;
  out.throughput = static_cast<double>(n) * u;
  out.case_name = case_name;
  return out;
}

std::vector<double> approximate_model_per_feeder(int NumberOfFeeders,
                                                 int NumberOfPrimary,
                                                 bool BiDirectional,
                                                 const std::vector<std::vector<int>>& buffer_capacity) {
  const int m = NumberOfFeeders;
  const int n = NumberOfPrimary;

  std::vector<int> c_f, c_b;
  extract_buffer_caps(m, buffer_capacity, c_f, c_b);

  const bool has_f = std::any_of(c_f.begin(), c_f.end(), [](int x){ return x > 0; });

  if (!BiDirectional) {
    return has_f ? approx_one_way_forward_buffers_trace(n, m, c_f)
                 : approx_one_way_no_buffers_trace(n, m);
  }
  auto [_u_f, u_b] = approx_two_way_per_feeder_buffers_trace(n, m, c_f, c_b);
  return u_b;
}

} // namespace meshsim
