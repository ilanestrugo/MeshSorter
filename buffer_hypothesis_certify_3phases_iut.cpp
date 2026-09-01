
/*
  buffer_hypothesis_certify_3phases_iut.cpp
  C++ translation of buffer_hypothesis_certify_3phases_iut.py. fileciteturn0file0

  Build:
    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    cmake --build . -j

  Run (example):
    ./buffer_certify --m 6 --n 8 --B 12 --subset one_way_monotone --Num_of_blocks 150 --block_size 200 --eps 0.001 --alpha 0.01 --log_prefix run

  Phase-1 pilot horizons:
    --Num_of_blocks    blocks per Phase-1 pilot run for the structured set S (also the
                       length of each Phase-2 extension run of x_S).
    --Num_of_blocks_C  blocks per Phase-1 pilot run for the competitors C.
                       Defaults to --Num_of_blocks (i.e. the symmetric design).

  Notes:
  - This translation preserves the algorithmic structure. If the Python is slow because the number of
    configurations explodes combinatorially (especially for --bidirectional true), C++ will help
    but may not fix the core combinatorial growth.
*/

#include "MeshSim.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Config {
  std::vector<int> c_f;
  std::vector<int> c_b;
};

std::string vec_to_string(const std::vector<int>& v) {
  std::ostringstream oss;
  oss << "(";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i) oss << ", ";
    oss << v[i];
  }
  oss << ")";
  return oss.str();
}

std::string fmt_cfg(const Config& cfg) {
  std::ostringstream oss;
  oss << "forward=" << "["; for (size_t i=0;i<cfg.c_f.size();++i){ if(i)oss<<","; oss<<cfg.c_f[i]; } oss << "]";
  oss << " | backward=" << "["; for (size_t i=0;i<cfg.c_b.size();++i){ if(i)oss<<","; oss<<cfg.c_b[i]; } oss << "]";
  return oss.str();
}

// ----- Normal inverse CDF (Acklam approximation) -----
double normal_inv_cdf(double p) {
  // Returns z such that Phi(z) = p, for 0<p<1.
  // For p==0 or p==1, returns +/-inf.
  if (p <= 0.0) return -std::numeric_limits<double>::infinity();
  if (p >= 1.0) return  std::numeric_limits<double>::infinity();

  static const double a1 = -3.969683028665376e+01;
  static const double a2 =  2.209460984245205e+02;
  static const double a3 = -2.759285104469687e+02;
  static const double a4 =  1.383577518672690e+02;
  static const double a5 = -3.066479806614716e+01;
  static const double a6 =  2.506628277459239e+00;

  static const double b1 = -5.447609879822406e+01;
  static const double b2 =  1.615858368580409e+02;
  static const double b3 = -1.556989798598866e+02;
  static const double b4 =  6.680131188771972e+01;
  static const double b5 = -1.328068155288572e+01;

  static const double c1 = -7.784894002430293e-03;
  static const double c2 = -3.223964580411365e-01;
  static const double c3 = -2.400758277161838e+00;
  static const double c4 = -2.549732539343734e+00;
  static const double c5 =  4.374664141464968e+00;
  static const double c6 =  2.938163982698783e+00;

  static const double d1 =  7.784695709041462e-03;
  static const double d2 =  3.224671290700398e-01;
  static const double d3 =  2.445134137142996e+00;
  static const double d4 =  3.754408661907416e+00;

  const double plow  = 0.02425;
  const double phigh = 1.0 - plow;

  double q = 0.0;
  double r = 0.0;
  if (p < plow) {
    q = std::sqrt(-2.0 * std::log(p));
    return (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
           ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
  }
  if (p > phigh) {
    q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
             ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
  }

  q = p - 0.5;
  r = q * q;
  return (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6) * q /
         (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0);
}

// ----- smallest epsilon satisfying the one-sided validation inequality -----
// Smallest eps >= 0 such that
//   mu_S * (1 + eps) - mu_C >= z * sqrt( (1+eps)^2 * var_S/m_S + var_C/m_C ).
// Epsilon appears on both sides, so the root is obtained numerically. The gap
// function is continuous and increasing in eps whenever mu_S > z*sqrt(var_S/m_S),
// so we bracket it by doubling and then bisect. Returns +infinity when no finite
// eps satisfies the inequality.
double smallest_eps_one_sided(double muS, double varS, int mS,
                              double muC, double varC, int mC,
                              double z) {
  const double inf = std::numeric_limits<double>::infinity();
  if (!(muS > 0.0)) return inf;

  const double aS = varS / static_cast<double>(std::max(1, mS));
  const double aC = varC / static_cast<double>(std::max(1, mC));

  auto gap = [&](double e) {
    const double one_plus_e = 1.0 + e;
    return muS * one_plus_e - muC - z * std::sqrt(one_plus_e * one_plus_e * aS + aC);
  };

  if (gap(0.0) >= 0.0) return 0.0;

  double hi = 1e-6;
  for (int i = 0; gap(hi) < 0.0; ++i) {
    hi *= 2.0;
    if (i > 200 || !std::isfinite(hi)) return inf;
  }

  double lo = 0.0;
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (gap(mid) >= 0.0) hi = mid; else lo = mid;
    if (hi - lo <= 1e-15 * std::max(1.0, hi)) break;
  }
  return hi;
}

// ----- compositions / allocation enumerators -----
template <class Fn>
void for_each_composition_rec(int total, int parts, std::vector<int>& cur, Fn&& fn) {
  if (parts <= 0) return;
  if (parts == 1) {
    cur.push_back(total);
    fn(cur);
    cur.pop_back();
    return;
  }
  for (int x = 0; x <= total; ++x) {
    cur.push_back(x);
    for_each_composition_rec(total - x, parts - 1, cur, fn);
    cur.pop_back();
  }
}

template <class Fn>
void for_each_composition(int total, int parts, Fn&& fn) {
  std::vector<int> cur;
  cur.reserve(static_cast<size_t>(std::max(0, parts)));
  for_each_composition_rec(total, parts, cur, fn);
}

std::vector<Config> enumerate_allocations_one_way(int m, int B) {
  std::vector<Config> allocs;
  if (m <= 0) return allocs;
  if (B < 0) throw std::invalid_argument("B must be >= 0");
  if (m == 1) {
    Config cfg;
    cfg.c_f = {0};
    cfg.c_b = {0};
    allocs.push_back(cfg);
    return allocs;
  }

  for_each_composition(B, m - 1, [&](const std::vector<int>& comp){
    Config cfg;
    cfg.c_f.assign(m, 0);
    cfg.c_b.assign(m, 0);
    for (int i = 1; i < m; ++i) cfg.c_f[i] = comp[static_cast<size_t>(i - 1)];
    allocs.push_back(std::move(cfg));
  });
  return allocs;
}

void enumerate_two_way_splits_rec(
    int idx,
    const std::vector<int>& totals, // length m-1
    std::vector<int>& f,
    std::vector<int>& b,
    std::vector<Config>& out) {

  if (idx == static_cast<int>(totals.size())) {
    Config cfg;
    cfg.c_f = f;
    cfg.c_b = b;
    out.push_back(std::move(cfg));
    return;
  }
  const int x = totals[static_cast<size_t>(idx)];
  for (int k = 0; k <= x; ++k) {
    f[static_cast<size_t>(idx + 1)] = k;
    b[static_cast<size_t>(idx + 1)] = x - k;
    enumerate_two_way_splits_rec(idx + 1, totals, f, b, out);
  }
}

std::vector<Config> enumerate_allocations_two_way(int m, int B) {
  std::vector<Config> allocs;
  if (m <= 0) return allocs;
  if (B < 0) throw std::invalid_argument("B must be >= 0");
  if (m == 1) {
    Config cfg;
    cfg.c_f = {0};
    cfg.c_b = {0};
    allocs.push_back(cfg);
    return allocs;
  }

  for_each_composition(B, m - 1, [&](const std::vector<int>& comp_total){
    std::vector<int> f(m, 0);
    std::vector<int> b(m, 0);
    enumerate_two_way_splits_rec(0, comp_total, f, b, allocs);
  });
  return allocs;
}

// ----- subset definitions -----
bool is_monotone_non_decreasing(const std::vector<int>& x) {
  for (size_t i = 0; i + 1 < x.size(); ++i) {
    if (x[i] > x[i + 1]) return false;
  }
  return true;
}

bool subset_one_way_monotone(const Config& cfg) {
  if (cfg.c_f.empty() || cfg.c_b.empty()) return false;
  if (cfg.c_f[0] != 0) return false;
  for (int v : cfg.c_b) if (v != 0) return false;
  return is_monotone_non_decreasing(cfg.c_f);
}

bool subset_two_way_monotone_backward_only(const Config& cfg) {
  if (cfg.c_f.empty() || cfg.c_b.empty()) return false;
  if (cfg.c_f[0] != 0 || cfg.c_b[0] != 0) return false;
  for (int v : cfg.c_f) if (v != 0) return false;
  return is_monotone_non_decreasing(cfg.c_b);
}

bool subset_two_way_monotone(const Config& cfg) {
  if (cfg.c_f.empty() || cfg.c_b.empty()) return false;
  if (cfg.c_f[0] != 0 || cfg.c_b[0] != 0) return false;
  return is_monotone_non_decreasing(cfg.c_f) && is_monotone_non_decreasing(cfg.c_b);
}

// ----- stats helpers -----
struct RunAgg {
  int n = 0;
  double mean = 0.0;
  double m2 = 0.0;
  void update(double x) {
    n += 1;
    const double d = x - mean;
    mean += d / static_cast<double>(n);
    m2 += d * (x - mean);
  }
  double var() const {
    if (n <= 1) return std::numeric_limits<double>::quiet_NaN();
    return m2 / static_cast<double>(n - 1);
  }
};

struct BlockAgg {
  int n = 0;
  double mean = 0.0;
  double m2 = 0.0;
  void update_many(const std::vector<double>& xs) {
    for (double v : xs) {
      n += 1;
      const double d = v - mean;
      mean += d / static_cast<double>(n);
      m2 += d * (v - mean);
    }
  }
  double var() const {
    if (n <= 1) return std::numeric_limits<double>::quiet_NaN();
    return m2 / static_cast<double>(n - 1);
  }
};

double sample_mean(const std::vector<double>& x) {
  if (x.empty()) return 0.0;
  double s = 0.0;
  for (double v : x) s += v;
  return s / static_cast<double>(x.size());
}

double sample_var_ddof1(const std::vector<double>& x) {
  const size_t n = x.size();
  if (n < 2) return std::numeric_limits<double>::quiet_NaN();
  const double mu = sample_mean(x);
  double s2 = 0.0;
  for (double v : x) {
    const double d = v - mu;
    s2 += d * d;
  }
  return s2 / static_cast<double>(n - 1);
}

// Return post-warmup per-block utilization samples u_b.
std::vector<double> simulate_blocks_series(
    int NumberOfFeeders,
    int NumberOfPrimary,
    bool BiDirectional,
    int T_time,
    int seed,
    const Config& cfg,
    int block_size,
    bool fast_mode) {

  std::vector<std::vector<int>> cap;
  cap.reserve(2);
  cap.push_back(cfg.c_f);
  cap.push_back(cfg.c_b);

  const bool collect = !fast_mode;

  auto sim = meshsim::run_simulation(
      NumberOfFeeders, NumberOfPrimary, BiDirectional,
      T_time, cap, seed,
      /*Verbal=*/false,
      block_size,
      collect
  );

  const int k_warm = meshsim::detect_warmup_from_blocks(sim.block_counts, sim.block_sizes, /*strict=*/true);

  const size_t nb = sim.block_counts.size();
  std::vector<double> u;
  u.reserve(nb);

  for (size_t i = 0; i < nb; ++i) {
    const double bc = static_cast<double>(sim.block_counts[i]);
    const double bs = static_cast<double>(sim.block_sizes[i]);
    if (bs <= 0.0) u.push_back(0.0);
    else u.push_back(bc / (bs * static_cast<double>(NumberOfPrimary)));
  }

  if (static_cast<size_t>(k_warm) >= u.size()) return {};
  return std::vector<double>(u.begin() + k_warm, u.end());
}

// ----- core 3-phase procedure -----
struct ThreePhaseResult {
  int K = 0;
  int S_size = 0;
  int C_size = 0;

  int S_star_k = -1;
  Config S_star_cfg;
  int m_S_star_blocks = 0;

  int worst_C_k = -1;
  Config worst_C_cfg;
  int m_worst_C_blocks = 0;

  double epsilon_hat = 0.0;
  bool certified = false;

  std::string phase1_S_csv;
  std::string phase3_detail_csv;

  double muS3 = 0.0;
  double varS3 = 0.0;
  int M_S_eff3 = 1;

  double worst_mean = 0.0;
  double worst_var = 0.0;
  int worst_m_eff = 1;

  int phase2_rounds = 0;
  int phase2_capped_infeasible = 0;
};

ThreePhaseResult three_phase_iut_epsilon(
    const std::vector<Config>& arm_configs,
    const std::function<bool(const Config&)>& is_in_subset,
    int m,
    int n,
    bool bidirectional,
    int Num_of_blocks,
    int block_size,
    double alpha,
    double eps_target,
    int n0,
    int seed0,
    bool fast_mode,
    const std::string& log_csv_prefix,
    double var_fallback = 0.25,
    double m_mult = 1.0,
    int min_blocks_C = 500,
    int Num_of_blocks_C = 0) {

  ThreePhaseResult res;

  const int K = static_cast<int>(arm_configs.size());
  if (K == 0) throw std::invalid_argument("No configurations provided.");

  std::vector<int> S_ids;
  std::vector<int> C_ids;
  S_ids.reserve(K);
  C_ids.reserve(K);

  for (int k = 0; k < K; ++k) {
    if (is_in_subset(arm_configs[k])) S_ids.push_back(k);
  }
  std::set<int> S_set(S_ids.begin(), S_ids.end());
  for (int k = 0; k < K; ++k) {
    if (S_set.find(k) == S_set.end()) C_ids.push_back(k);
  }
  if (S_ids.empty() || C_ids.empty()) {
    throw std::invalid_argument("Both subset S and complement C must be non-empty.");
  }

  const int base_blocks = Num_of_blocks;
  if (base_blocks <= 0) throw std::invalid_argument("Num_of_blocks must be positive.");
  const int T_base = base_blocks * block_size;

  // Phase-1 pilot horizon for the competitors. The reference set S keeps the long
  // horizon (it selects x_S and feeds the sampling plan); C may be piloted much more
  // cheaply, since its pilot only plans m(c) and is superseded by the Phase-3 run.
  const int base_blocks_C = (Num_of_blocks_C > 0) ? Num_of_blocks_C : base_blocks;
  const int T_base_C = base_blocks_C * block_size;
  std::cout << "[Phase 1] pilot horizons: S=" << base_blocks
            << " blocks, C=" << base_blocks_C << " blocks (block_size=" << block_size << ")\n";

  // Phase 1: pilot runs
  int next_seed = seed0;

  std::vector<double> pilot_mean(K, 0.0);
  std::vector<double> pilot_var(K, std::numeric_limits<double>::quiet_NaN());
  std::vector<int> pilot_m_eff(K, 0);

  for (int k = 0; k < K; ++k) {
    RunAgg run_stats;
    double last_block_var = std::numeric_limits<double>::quiet_NaN();
    int last_m_eff = 0;

    const int T_pilot = (S_set.find(k) != S_set.end()) ? T_base : T_base_C;

    for (int r = 0; r < n0; ++r) {
      const int seed = next_seed++;
      const auto u_post = simulate_blocks_series(m, n, bidirectional, T_pilot, seed, arm_configs[k], block_size, fast_mode);
      double mu = 0.0;
      double v = var_fallback;
      int m_eff = 0;

      if (!u_post.empty()) {
        mu = sample_mean(u_post);
        m_eff = static_cast<int>(u_post.size());
        if (u_post.size() >= 2) {
          const double vv = sample_var_ddof1(u_post);
          v = (std::isfinite(vv) ? vv : var_fallback);
        } else {
          v = var_fallback;
        }
        v = std::max(0.0, std::min(var_fallback, v));
      } else {
        mu = 0.0;
        v = var_fallback;
        m_eff = 0;
      }

      run_stats.update(mu);
      last_block_var = v;
      last_m_eff = m_eff;
    }

    pilot_mean[k] = run_stats.mean;
    pilot_var[k]  = std::isfinite(last_block_var) ? last_block_var : var_fallback;
    pilot_m_eff[k] = last_m_eff;
  }

  // Select S*
  int S_star = S_ids[0];
  for (int k : S_ids) {
    if (pilot_mean[k] > pilot_mean[S_star]) S_star = k;
  }
  const Config S_star_cfg = arm_configs[S_star];

  std::cout << "[Phase 1] Selected S*: k=" << S_star
            << " pilot_mean=" << std::fixed << std::setprecision(6) << pilot_mean[S_star]
            << " cfg=" << fmt_cfg(S_star_cfg) << "\n";

  // Phase-1 S CSV
  res.phase1_S_csv = log_csv_prefix + "_phase1_S_detail.csv";
  if (fs::path(res.phase1_S_csv).has_parent_path()) {
    fs::create_directories(fs::path(res.phase1_S_csv).parent_path());
  }
  {
    std::ofstream f(res.phase1_S_csv);
    f << "k,Num_of_blocks,T_time,mean,var,m_eff_blocks,cfg_forward,cfg_backward\n";
    for (int k : S_ids) {
      f << k << "," << base_blocks << "," << T_base << ","
        << std::setprecision(10) << pilot_mean[k] << ","
        << std::setprecision(10) << pilot_var[k] << ","
        << pilot_m_eff[k] << ","
        << "\"" << vec_to_string(arm_configs[k].c_f) << "\"" << ","
        << "\"" << vec_to_string(arm_configs[k].c_b) << "\"" << "\n";
    }
  }

  // Phase 2: extend S* until max m(c) <= m(S*)
  const double z = normal_inv_cdf(1.0 - alpha);
  BlockAgg S_blocks_agg;

  // Start S* with one base run after pilot
  {
    const int seed = next_seed++;
    const auto uS = simulate_blocks_series(m, n, bidirectional, T_base, seed, S_star_cfg, block_size, fast_mode);
    S_blocks_agg.update_many(uS);
  }
  if (S_blocks_agg.n == 0) {
    // fallback: use pilot mean as one sample
    S_blocks_agg.update_many(std::vector<double>{pilot_mean[S_star]});
  }

  std::map<int, int> m_blocks;
  std::vector<int> infeasible;
  int extension_rounds = 0;

  while (true) {
    extension_rounds += 1;

    const double muS = S_blocks_agg.mean;
    double varS = (S_blocks_agg.n > 1) ? S_blocks_agg.var() : var_fallback;
    if (!std::isfinite(varS)) varS = var_fallback;
    varS = std::max(0.0, std::min(var_fallback, varS));

    const int M_S_eff = std::max(1, S_blocks_agg.n);
    const int M_S_blocks = std::max(1, static_cast<int>(std::ceil(static_cast<double>(M_S_eff))));

    infeasible.clear();
    m_blocks.clear();

    int max_m = 1;
    for (int j : C_ids) {
      const double mu_j = pilot_mean[j];
      double var_j = std::isfinite(pilot_var[j]) ? pilot_var[j] : var_fallback;
      var_j = std::max(0.0, std::min(var_fallback, var_j));

      const double Delta = muS * (1.0 + eps_target) - mu_j;
      const double one_plus_eps2 = (1.0 + eps_target) * (1.0 + eps_target);
      const double denom = std::pow(Delta / z, 2.0) - one_plus_eps2 * (varS / static_cast<double>(M_S_eff));

      int m_req = 1;
      const bool infeasible_j = (Delta <= 0.0 || denom <= 0.0);
      if (infeasible_j) {
        infeasible.push_back(j);
        m_req = M_S_blocks;
      } else {
        m_req = static_cast<int>(std::ceil(var_j / denom));
      }

      m_req = std::max(1, m_req);
      if (m_mult != 1.0) {
        m_req = static_cast<int>(std::ceil(m_mult * static_cast<double>(m_req)));
      }

      m_req = std::max(min_blocks_C, m_req);

      // The planned effort for an infeasible competitor is capped at m(x_S). Applying
      // that cap before the m_mult inflation would leave m_req = m_mult * M_S_blocks,
      // which the extension loop below can never satisfy: each round raises M_S_blocks
      // and the requirement rises with it. Re-apply the cap after inflation.
      if (infeasible_j) m_req = std::min(m_req, M_S_blocks);

      m_blocks[j] = m_req;
      if (m_req > max_m) max_m = m_req;
    }

    std::cout << "[Phase 2] round=" << extension_rounds
              << " S*_blocks_eff=" << M_S_eff
              << " max_m(C)=" << max_m
              << " capped_infeasible=" << infeasible.size()
              << "\n";

    if (max_m <= M_S_blocks) break;

    // Extend S* by another base run
    const int seed = next_seed++;
    auto uS_more = simulate_blocks_series(m, n, bidirectional, T_base, seed, S_star_cfg, block_size, fast_mode);
    if (uS_more.empty()) uS_more = std::vector<double>{muS};
    S_blocks_agg.update_many(uS_more);
  }

  const int M_S_blocks_final = std::max(1, static_cast<int>(std::ceil(static_cast<double>(std::max(1, S_blocks_agg.n)))));

  // Phase 3: independent validation
  int next_seed3 = seed0 + 10'000'000;

  // Re-run S*
  std::vector<double> uS3 = simulate_blocks_series(m, n, bidirectional, M_S_blocks_final * block_size, next_seed3++, S_star_cfg, block_size, fast_mode);
  double muS3 = 0.0;
  double varS3 = var_fallback;
  int M_S_eff3 = 1;
  if (!uS3.empty()) {
    muS3 = sample_mean(uS3);
    M_S_eff3 = static_cast<int>(uS3.size());
    if (uS3.size() >= 2) {
      const double vv = sample_var_ddof1(uS3);
      varS3 = std::isfinite(vv) ? vv : var_fallback;
    } else {
      varS3 = var_fallback;
    }
    varS3 = std::max(0.0, std::min(var_fallback, varS3));
  } else {
    muS3 = S_blocks_agg.mean;
    varS3 = var_fallback;
    M_S_eff3 = 1;
  }

  // Detail CSV
  res.phase3_detail_csv = log_csv_prefix + "_phase3_detail.csv";
  if (fs::path(res.phase3_detail_csv).has_parent_path()) {
    fs::create_directories(fs::path(res.phase3_detail_csv).parent_path());
  }

  int worst_k = -1;
  double worst_eps = -1.0;
  double worst_mean = 0.0;
  double worst_var = 0.0;
  int worst_mblocks = 0;
  int worst_m_eff = 1;
  Config worst_cfg;

  {
    std::ofstream f(res.phase3_detail_csv);
    f << "k,m_blocks,T_time,mean,var,m_eff_blocks,eps_i,cfg_forward,cfg_backward\n";
    for (int j : C_ids) {
      const int seed = next_seed3++;
      int m_j = m_blocks[j];
      if (m_j > M_S_blocks_final) m_j = M_S_blocks_final;
      const int T_j = m_j * block_size;

      auto uj = simulate_blocks_series(m, n, bidirectional, T_j, seed, arm_configs[j], block_size, fast_mode);

      double muj = 0.0;
      double varj = var_fallback;
      int m_j_eff = 1;
      if (!uj.empty()) {
        muj = sample_mean(uj);
        m_j_eff = static_cast<int>(uj.size());
        if (uj.size() >= 2) {
          const double vv = sample_var_ddof1(uj);
          varj = std::isfinite(vv) ? vv : var_fallback;
        } else {
          varj = var_fallback;
        }
        varj = std::max(0.0, std::min(var_fallback, varj));
      } else {
        muj = pilot_mean[j];
        varj = var_fallback;
        m_j_eff = 1;
      }

      const double eps_req = smallest_eps_one_sided(muS3, varS3, M_S_eff3,
                                                    muj, varj, m_j_eff, z);

      if (eps_req > worst_eps) {
        worst_eps = eps_req;
        worst_k = j;
        worst_mean = muj;
        worst_var = varj;
        worst_mblocks = m_j;
        worst_m_eff = m_j_eff;
        worst_cfg = arm_configs[j];
      }

      f << j << "," << m_j << "," << T_j << ","
        << std::setprecision(10) << muj << ","
        << std::setprecision(10) << varj << ","
        << m_j_eff << ","
        << std::setprecision(10) << eps_req << ","
        << "\"" << vec_to_string(arm_configs[j].c_f) << "\"" << ","
        << "\"" << vec_to_string(arm_configs[j].c_b) << "\"" << "\n";
    }
  }

  const double epsilon_hat = worst_eps;
  const bool certified = (epsilon_hat <= eps_target);

  res.K = K;
  res.S_size = static_cast<int>(S_ids.size());
  res.C_size = static_cast<int>(C_ids.size());
  res.S_star_k = S_star;
  res.S_star_cfg = S_star_cfg;
  res.m_S_star_blocks = M_S_blocks_final;
  res.worst_C_k = worst_k;
  res.worst_C_cfg = worst_cfg;
  res.m_worst_C_blocks = worst_mblocks;
  res.epsilon_hat = epsilon_hat;
  res.certified = certified;
  res.muS3 = muS3;
  res.varS3 = varS3;
  res.M_S_eff3 = M_S_eff3;
  res.worst_mean = worst_mean;
  res.worst_var = worst_var;
  res.worst_m_eff = worst_m_eff;
  res.phase2_rounds = extension_rounds;
  res.phase2_capped_infeasible = static_cast<int>(infeasible.size());

  return res;
}

// ----- basic CLI parsing -----
std::unordered_map<std::string, std::string> parse_args(int argc, char** argv) {
  std::unordered_map<std::string, std::string> out;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--", 0) == 0) {
      std::string key = a.substr(2);
      std::string val = "true";
      if (i + 1 < argc) {
        std::string nxt = argv[i + 1];
        if (nxt.rfind("--", 0) != 0) {
          val = nxt;
          i += 1;
        }
      }
      out[key] = val;
    }
  }
  return out;
}

bool as_bool(const std::string& s) {
  std::string t;
  t.reserve(s.size());
  for (char c : s) t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  return (t == "1" || t == "true" || t == "t" || t == "yes" || t == "y");
}

int as_int(const std::unordered_map<std::string,std::string>& a, const std::string& k, std::optional<int> def = std::nullopt) {
  auto it = a.find(k);
  if (it == a.end()) {
    if (def.has_value()) return *def;
    throw std::runtime_error("Missing required --" + k);
  }
  return std::stoi(it->second);
}

double as_double(const std::unordered_map<std::string,std::string>& a, const std::string& k, std::optional<double> def = std::nullopt) {
  auto it = a.find(k);
  if (it == a.end()) {
    if (def.has_value()) return *def;
    throw std::runtime_error("Missing required --" + k);
  }
  return std::stod(it->second);
}

std::string as_string(const std::unordered_map<std::string,std::string>& a, const std::string& k, std::optional<std::string> def = std::nullopt) {
  auto it = a.find(k);
  if (it == a.end()) {
    if (def.has_value()) return *def;
    throw std::runtime_error("Missing required --" + k);
  }
  return it->second;
}

std::string now_iso_seconds() {
  auto tp = std::chrono::system_clock::now();
  std::time_t tt = std::chrono::system_clock::to_time_t(tp);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &tt);
#else
  localtime_r(&tt, &tm);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S");
  return oss.str();
}

} // namespace

int main(int argc, char** argv) {
  try {
    auto args = parse_args(argc, argv);

    const int m = as_int(args, "m");
    const int n = as_int(args, "n");
    const int B = as_int(args, "B");

    const int Num_of_blocks = as_int(args, "Num_of_blocks", 150);
    const int Num_of_blocks_C = as_int(args, "Num_of_blocks_C", Num_of_blocks);
    const int block_size = as_int(args, "block_size", 200);
    const bool bidirectional = as_bool(as_string(args, "bidirectional", std::string("false")));

    const std::string subset_name = as_string(args, "subset");
    const double alpha = as_double(args, "alpha", 0.01);
    const double eps = as_double(args, "eps", 0.001);
    const int n0 = as_int(args, "n0", 1);

    const double m_mult = as_double(args, "m_mult", 1.0);
    const int min_blocks_C = as_int(args, "min_blocks_C", 500);
    const int seed0 = as_int(args, "seed0", 123);
    const bool fast = (args.find("fast") != args.end()) ? as_bool(args["fast"]) : false;

    const std::string log_prefix = as_string(args, "log_prefix", std::string("run"));

    const std::optional<std::string> best_summary_csv =
        (args.find("best_summary_csv") != args.end()) ? std::optional<std::string>(args["best_summary_csv"])
                                                      : std::nullopt;

    // Enumerate allocations
    std::vector<Config> arm_configs = bidirectional
      ? enumerate_allocations_two_way(m, B)
      : enumerate_allocations_one_way(m, B);

    std::function<bool(const Config&)> subset_fn;
    std::string subset_lower = subset_name;
    std::transform(subset_lower.begin(), subset_lower.end(), subset_lower.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

    if (subset_lower == "one_way_monotone") subset_fn = subset_one_way_monotone;
    else if (subset_lower == "two_way_monotone_backward_only") subset_fn = subset_two_way_monotone_backward_only;
    else if (subset_lower == "two_way_monotone") subset_fn = subset_two_way_monotone;
    else throw std::runtime_error("Unknown --subset. Use one_way_monotone | two_way_monotone_backward_only | two_way_monotone");

    auto t0 = std::chrono::steady_clock::now();

    auto res = three_phase_iut_epsilon(
        arm_configs,
        subset_fn,
        m, n, bidirectional,
        Num_of_blocks,
        block_size,
        alpha,
        eps,
        n0,
        seed0,
        fast,
        log_prefix,
        /*var_fallback=*/0.25,
        m_mult,
        min_blocks_C,
        Num_of_blocks_C
    );

    auto t1 = std::chrono::steady_clock::now();
    const double elapsed_sec = std::chrono::duration<double>(t1 - t0).count();

    std::cout << "\n=== RESULT ===\n";
    std::cout << "Certified at eps=" << eps << ": " << (res.certified ? "true" : "false") << "\n";
    std::cout << "epsilon_hat=" << std::fixed << std::setprecision(6) << res.epsilon_hat << "\n";
    std::cout << "K=" << res.K << " |S|=" << res.S_size << " |C|=" << res.C_size << "\n";
    std::cout << "S*: k=" << res.S_star_k << " blocks=" << res.m_S_star_blocks << " cfg=" << fmt_cfg(res.S_star_cfg) << "\n";
    std::cout << "Worst in C: k=" << res.worst_C_k << " blocks=" << res.m_worst_C_blocks << "\n";
    std::cout << "Phase1 S-detail CSV: " << res.phase1_S_csv << "\n";
    std::cout << "Phase3 detail CSV: " << res.phase3_detail_csv << "\n";
    std::cout << "Elapsed: " << std::fixed << std::setprecision(3) << elapsed_sec << " sec\n";

    // Reporting CI half-widths (Normal)
    const double z_two = normal_inv_cdf(1.0 - alpha / 2.0);
    const double z_one = normal_inv_cdf(1.0 - alpha);
    const int ms_eff = std::max(1, res.M_S_eff3);
    const int wc_eff = std::max(1, res.worst_m_eff);

    const double varS = std::max(0.0, res.varS3);
    const double varC = std::max(0.0, res.worst_var);

    const double s_hw_two = z_two * std::sqrt(varS / static_cast<double>(ms_eff));
    const double s_hw_one = z_one * std::sqrt(varS / static_cast<double>(ms_eff));
    const double c_hw_two = z_two * std::sqrt(varC / static_cast<double>(wc_eff));
    const double c_hw_one = z_one * std::sqrt(varC / static_cast<double>(wc_eff));

    if (best_summary_csv.has_value()) {
      fs::path p(*best_summary_csv);
      if (p.has_parent_path()) fs::create_directories(p.parent_path());

      const bool file_exists = fs::exists(p);
      const bool write_header = (!file_exists) || (fs::file_size(p) == 0);

      std::ofstream f(*best_summary_csv, std::ios::app);
      if (write_header) {
        f << "timestamp,elapsed_sec,m,n,B,bidirectional,K,"
          << "s_size,s_mean,s_CI_halfwidth_onesided,s_CI_halfwidth_twosided,s_buffers,"
          << "c_size,c_mean,c_CI_halfwidth_onesided,c_CI_halfwidth_twosided,c_buffers,"
          << "m_S_star_blocks,m_worst_C_blocks,epsilon_hat,certified\n";
      }

      f << now_iso_seconds() << ","
        << std::setprecision(6) << elapsed_sec << ","
        << m << "," << n << "," << B << "," << (bidirectional ? 1 : 0) << ","
        << res.K << ","
        << res.S_size << ","
        << std::setprecision(10) << res.muS3 << ","
        << std::setprecision(10) << s_hw_one << ","
        << std::setprecision(10) << s_hw_two << ","
        << "\"" << fmt_cfg(res.S_star_cfg) << "\"" << ","
        << res.C_size << ","
        << std::setprecision(10) << res.worst_mean << ","
        << std::setprecision(10) << c_hw_one << ","
        << std::setprecision(10) << c_hw_two << ","
        << "\"" << fmt_cfg(res.worst_C_cfg) << "\"" << ","
        << res.m_S_star_blocks << ","
        << res.m_worst_C_blocks << ","
        << std::setprecision(10) << res.epsilon_hat << ","
        << (res.certified ? 1 : 0)
        << "\n";

      std::cout << "Run-summary appended to: " << *best_summary_csv << "\n";
    }

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    std::cerr << "Usage example:\n"
              << "  ./buffer_certify --m 6 --n 8 --B 12 --subset one_way_monotone\n";
    return 1;
  }
}
