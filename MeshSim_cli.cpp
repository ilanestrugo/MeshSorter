/*
  MeshSim_cli.cpp
  Standalone CLI for the C++ MeshSim library.

  This mirrors the CLI behavior in MeshSim.py (run_mode: sim/model/both),
  but uses the C++ implementation for speed. fileciteturn0file1

  Build (CMake):
    cmake -DBUILD_MESHSIM_CLI=ON ..

  Example:
    ./meshsim_cli --NumberOfFeeders 6 --NumberOfPrimary 8 --BiDirectional false \
      --buffer_capacity "[[0,1,2,3,4,5],[0,0,0,0,0,0]]" --T 30000 --seed 1 --block_size 200 \
      --run_mode both --ci_enable --ci_level 0.95 --output_csv MeshSim_results.csv
*/

#include "MeshSim.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

static std::string to_lower(std::string s) {
  for (char& c : s) {
    if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
  }
  return s;
}

static bool parse_bool(const std::string& s) {
  const std::string v = to_lower(s);
  if (v == "1" || v == "true" || v == "t" || v == "yes" || v == "y" || v == "on") return true;
  if (v == "0" || v == "false" || v == "f" || v == "no" || v == "n" || v == "off") return false;
  throw std::runtime_error("Cannot parse boolean value from: '" + s + "'");
}

static std::vector<int> parse_int_list_loose(const std::string& s) {
  // Extract all integers from a string like "[1, 2, 3]" or "1,2,3".
  std::vector<int> out;
  long long cur = 0;
  bool in_num = false;
  bool neg = false;

  auto flush = [&]() {
    if (in_num) {
      long long v = neg ? -cur : cur;
      out.push_back(static_cast<int>(v));
      cur = 0;
      in_num = false;
      neg = false;
    }
  };

  for (size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (!in_num) {
      if (c == '-') {
        neg = true;
        in_num = true;
        cur = 0;
      } else if (c >= '0' && c <= '9') {
        neg = false;
        in_num = true;
        cur = static_cast<long long>(c - '0');
      }
    } else {
      if (c >= '0' && c <= '9') {
        cur = cur * 10 + static_cast<long long>(c - '0');
      } else {
        flush();
      }
    }
  }
  flush();
  return out;
}

static std::string trim(const std::string& s) {
  size_t a = 0;
  while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\n' || s[a] == '\r')) ++a;
  size_t b = s.size();
  while (b > a && (s[b-1] == ' ' || s[b-1] == '\t' || s[b-1] == '\n' || s[b-1] == '\r')) --b;
  return s.substr(a, b - a);
}

static bool parse_nested_2d_list(const std::string& s_in, std::vector<std::vector<int>>& out) {
  // Very small parser for formats like: [[1,2,3],[4,5,6]]
  // Accepts optional whitespace.
  const std::string s = trim(s_in);
  if (s.size() < 4) return false;
  if (s[0] != '[') return false;

  // If it's a 1D list like [1,2,3], not 2D.
  size_t first_non_ws = 1;
  while (first_non_ws < s.size() && (s[first_non_ws] == ' ' || s[first_non_ws] == '\t')) ++first_non_ws;
  if (first_non_ws >= s.size() || s[first_non_ws] != '[') return false;

  size_t i = 0;
  auto skip_ws = [&]() {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r')) ++i;
  };

  out.clear();
  i = 0;
  skip_ws();
  if (i >= s.size() || s[i] != '[') return false;
  ++i;

  while (true) {
    skip_ws();
    if (i >= s.size()) return false;
    if (s[i] == ']') { ++i; break; }

    if (s[i] != '[') return false;
    ++i;
    skip_ws();

    std::vector<int> inner;
    while (true) {
      skip_ws();
      if (i >= s.size()) return false;
      if (s[i] == ']') { ++i; break; }

      bool neg = false;
      if (s[i] == '-') { neg = true; ++i; }
      if (i >= s.size() || !(s[i] >= '0' && s[i] <= '9')) return false;
      long long v = 0;
      while (i < s.size() && (s[i] >= '0' && s[i] <= '9')) {
        v = v * 10 + static_cast<long long>(s[i] - '0');
        ++i;
      }
      if (neg) v = -v;
      inner.push_back(static_cast<int>(v));

      skip_ws();
      if (i >= s.size()) return false;
      if (s[i] == ',') { ++i; continue; }
      if (s[i] == ']') { continue; }
      return false;
    }

    out.push_back(std::move(inner));

    skip_ws();
    if (i >= s.size()) return false;
    if (s[i] == ',') { ++i; continue; }
    if (s[i] == ']') { ++i; break; }
    return false;
  }

  return !out.empty();
}

static std::vector<std::vector<int>> parse_buffer_capacity(const std::string& user_str, int m) {
  // Returns cap[2][m]. If user provides only one list, it's treated as forward.
  // Supports:
  //   - "[[...],[...]]"  (Python-like)
  //   - "a,b,c; d,e,f"   (two lists separated by ';' or '|')
  //   - "a,b,c"          (forward only)
  std::vector<std::vector<int>> cap(2, std::vector<int>(std::max(0, m), 0));

  const std::string s = trim(user_str);
  if (s.empty()) {
    if (m > 0) { cap[0][0] = 0; cap[1][0] = 0; }
    return cap;
  }

  std::vector<std::vector<int>> parsed2d;
  if (parse_nested_2d_list(s, parsed2d)) {
    if (!parsed2d.empty()) {
      auto f = parsed2d[0];
      auto b = (parsed2d.size() >= 2) ? parsed2d[1] : std::vector<int>();
      for (int i = 0; i < m; ++i) {
        cap[0][i] = (i < static_cast<int>(f.size())) ? std::max(0, f[i]) : 0;
        cap[1][i] = (i < static_cast<int>(b.size())) ? std::max(0, b[i]) : 0;
      }
    }
  } else {
    // Try split on ';' or '|'
    size_t sep = s.find(';');
    if (sep == std::string::npos) sep = s.find('|');

    if (sep != std::string::npos) {
      const std::string a = trim(s.substr(0, sep));
      const std::string b = trim(s.substr(sep + 1));
      const auto f = parse_int_list_loose(a);
      const auto bb = parse_int_list_loose(b);
      for (int i = 0; i < m; ++i) {
        cap[0][i] = (i < static_cast<int>(f.size())) ? std::max(0, f[i]) : 0;
        cap[1][i] = (i < static_cast<int>(bb.size())) ? std::max(0, bb[i]) : 0;
      }
    } else {
      const auto f = parse_int_list_loose(s);
      for (int i = 0; i < m; ++i) cap[0][i] = (i < static_cast<int>(f.size())) ? std::max(0, f[i]) : 0;
      // cap[1] already zeros
    }
  }

  // by theory, first feeder is never blocked
  if (m > 0) {
    cap[0][0] = 0;
    cap[1][0] = 0;
  }

  return cap;
}

static std::string join_ints(const std::vector<int>& v, const std::string& sep) {
  std::ostringstream oss;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i) oss << sep;
    oss << v[i];
  }
  return oss.str();
}

static std::string join_doubles(const std::vector<double>& v, const std::string& sep, int prec = 6) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss << std::setprecision(prec);
  for (size_t i = 0; i < v.size(); ++i) {
    if (i) oss << sep;
    oss << v[i];
  }
  return oss.str();
}

static double tcrit_approx(double level, int n) {
  // Mirrors MeshSim.py _tcrit(). For df<=30 uses table, else uses normal approx.
  // Supported levels: 0.90, 0.95, 0.99.
  if (n <= 1) return std::numeric_limits<double>::infinity();
  const int df = n - 1;

  auto almost = [](double a, double b) { return std::abs(a - b) < 1e-12; };

  if (almost(level, 0.95)) {
    static const double T95[] = {
      0.0,
      12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262,
      2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093,
      2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
    };
    if (df <= 30) return T95[df];
    return 1.96;
  }

  if (almost(level, 0.90)) {
    static const double T90[] = {
      0.0,
      6.314, 2.920, 2.353, 2.132, 2.015, 1.943, 1.895, 1.860, 1.833,
      1.812, 1.796, 1.782, 1.771, 1.761, 1.753, 1.746, 1.740, 1.734, 1.729,
      1.725, 1.721, 1.717, 1.714, 1.711, 1.708, 1.706, 1.703, 1.701, 1.699, 1.697
    };
    if (df <= 30) return T90[df];
    return 1.645;
  }

  if (almost(level, 0.99)) {
    static const double T99[] = {
      0.0,
      63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250,
      3.169, 3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861,
      2.845, 2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750
    };
    if (df <= 30) return T99[df];
    return 2.576;
  }

  // Fallback: default to 0.95 behavior.
  return tcrit_approx(0.95, n);
}

static double tcrit_one_sided_approx(double level, int n) {
  // One-sided critical value for Student-t at confidence 'level'.
  // For example, level=0.95 returns t_{0.95, df=n-1}.
  // We reuse the two-sided tables when possible (one-sided 0.95 == two-sided 0.90),
  // otherwise use a normal quantile with a small df correction (Cornish-Fisher).
  if (n <= 1) return std::numeric_limits<double>::infinity();
  const int df = n - 1;

  auto almost = [](double a, double b) { return std::abs(a - b) < 1e-12; };

  // Exact via existing table: one-sided 0.95 corresponds to two-sided 0.90.
  if (almost(level, 0.95)) return tcrit_approx(0.90, n);

  // Normal quantiles for one-sided levels.
  double z = 1.6448536269514722; // default ~0.95
  if (almost(level, 0.90)) z = 1.2815515655446004;
  else if (almost(level, 0.95)) z = 1.6448536269514722;
  else if (almost(level, 0.99)) z = 2.3263478740408408;

  // Cornish-Fisher correction for Student-t quantile.
  const double d = static_cast<double>(df);
  const double z2 = z * z;
  const double z3 = z2 * z;
  const double z5 = z3 * z2;
  const double t = z + (z3 + z) / (4.0 * d) + (5.0 * z5 + 16.0 * z3 + 3.0 * z) / (96.0 * d * d);
  return t;
}

static void print_usage() {
  std::cout
    << "MeshSim CLI (C++)\n"
    << "\nRequired:\n"
    << "  --NumberOfFeeders <int>\n"
    << "  --NumberOfPrimary <int>\n"
    << "  --T <int>\n"
    << "\nOptional:\n"
    << "  --BiDirectional <true|false>        (default: true)\n"
    << "  --buffer_capacity <str>             e.g. \"[[0,1,2],[0,0,0]]\" or \"0,1,2;0,0,0\"\n"
    << "  --seed <int>                        (default: 0)\n"
    << "  --block_size <int>                  (default: 100)\n"
    << "  --Verbal                            (default: false)\n"
    << "  --run_mode <sim|model|both>         (default: both)\n"
    << "  --no_buffer_stats                   disable buffer stats collection (faster)\n"
    << "  --output_csv <path>                 (default: MeshSim_results.csv)\n"
    << "  --ci_enable                         enable throughput CI on post-warmup blocks\n"
    << "  --ci_level <0.90|0.95|0.99>         (default: 0.95)\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    int NumberOfFeeders = -1;
    int NumberOfPrimary = -1;
    bool BiDirectional = true;
    int T = -1;
    int seed = 0;
    int block_size = 100;
    bool Verbal = false;
    bool collect_buffer_stats = true;
    bool ci_enable = false;
    double ci_level = 0.95;
    bool ci_level_set = false;
    std::string run_mode = "both";
    std::string output_csv = "MeshSim_results.csv";
    std::string buffer_capacity_str;

    // Minimal argument parser
    for (int i = 1; i < argc; ++i) {
      std::string a = argv[i];
      if (a == "-h" || a == "--help") {
        print_usage();
        return 0;
      }

      auto need = [&](const std::string& flag) -> std::string {
        if (i + 1 >= argc) throw std::runtime_error("Missing value after " + flag);
        return std::string(argv[++i]);
      };

      if (a == "--NumberOfFeeders") {
        NumberOfFeeders = std::stoi(need(a));
      } else if (a == "--NumberOfPrimary") {
        NumberOfPrimary = std::stoi(need(a));
      } else if (a == "--BiDirectional") {
        BiDirectional = parse_bool(need(a));
      } else if (a == "--T") {
        T = std::stoi(need(a));
      } else if (a == "--seed") {
        seed = std::stoi(need(a));
      } else if (a == "--block_size") {
        block_size = std::stoi(need(a));
      } else if (a == "--buffer_capacity") {
        buffer_capacity_str = need(a);
      } else if (a == "--run_mode") {
        run_mode = to_lower(need(a));
      } else if (a == "--output_csv") {
        output_csv = need(a);
      } else if (a == "--Verbal") {
        Verbal = true;
      } else if (a == "--ci_enable") {
        ci_enable = true;
      } else if (a == "--ci_level") {
        ci_level = std::stod(need(a));
        ci_level_set = true;
        ci_enable = true; // if user sets level, assume they want CI
      } else if (a == "--no_buffer_stats" || a == "--fast") {
        collect_buffer_stats = false;
      } else {
        throw std::runtime_error("Unknown argument: " + a);
      }
    }

    if (NumberOfFeeders <= 0 || NumberOfPrimary <= 0 || T < 0) {
      print_usage();
      throw std::runtime_error("Missing required args (NumberOfFeeders/NumberOfPrimary/T). Use --help.");
    }

    if (run_mode != "sim" && run_mode != "model" && run_mode != "both") {
      throw std::runtime_error("--run_mode must be one of: sim, model, both");
    }

    const auto buffer_capacity = parse_buffer_capacity(buffer_capacity_str, NumberOfFeeders);

    const bool run_sim = (run_mode == "sim" || run_mode == "both");
    const bool run_model = (run_mode == "model" || run_mode == "both");

    // ---------------- Sim ----------------
    int warmup_blocks = 0;
    int warmup_steps = 0;
    double warmup_rate_curr = 0.0;
    double warmup_rate_avg = 0.0;

    int total_steps = 0;
    int total_entered = 0;
    double throughput = 0.0;
    double utilization = 0.0;

    std::vector<int> post_blocks;
    std::string acf_str;

    std::vector<double> avg_buffer;
    std::vector<int> max_buffer;
    std::vector<int> max_hits;

    std::vector<double> util_primary_per_feeder;
    std::string util_primary_str;

    double ci_half_two_sided = 0.0;
    double ci_half_one_sided = 0.0;
    int n_blocks_ci = 0;

    std::string buffers_str;

    if (run_sim) {
      std::cout << "\nSeed " << seed << ":\n";

      auto sim = meshsim::run_simulation(
        NumberOfFeeders,
        NumberOfPrimary,
        BiDirectional,
        T,
        buffer_capacity,
        seed,
        Verbal,
        block_size,
        collect_buffer_stats
      );

      warmup_blocks = meshsim::detect_warmup_from_blocks(sim.block_counts, sim.block_sizes, /*strict=*/true);
      warmup_steps = 0;
      for (int i = 0; i < warmup_blocks && i < static_cast<int>(sim.block_sizes.size()); ++i) warmup_steps += sim.block_sizes[i];

      // Warmup rates diagnostic
      std::vector<double> rates;
      rates.reserve(sim.block_counts.size());
      for (size_t i = 0; i < sim.block_counts.size(); ++i) {
        const int s = (i < sim.block_sizes.size()) ? sim.block_sizes[i] : 0;
        rates.push_back(s > 0 ? (static_cast<double>(sim.block_counts[i]) / static_cast<double>(s)) : 0.0);
      }
      warmup_rate_curr = (warmup_blocks < static_cast<int>(rates.size())) ? rates[warmup_blocks] : 0.0;
      if (warmup_blocks < static_cast<int>(rates.size()) - 1) {
        double sum_after = 0.0;
        int cnt = 0;
        for (size_t i = static_cast<size_t>(warmup_blocks + 1); i < rates.size(); ++i) { sum_after += rates[i]; ++cnt; }
        warmup_rate_avg = (cnt > 0) ? (sum_after / static_cast<double>(cnt)) : 0.0;
      }

      // Aggregate post-warmup
      if (collect_buffer_stats && sim.buffer_blocks.has_value()) {
        auto agg = meshsim::aggregate_post_warmup(
          warmup_blocks,
          sim.block_sizes,
          sim.block_counts,
          *sim.buffer_blocks,
          NumberOfPrimary
        );
        total_steps = agg.total_steps;
        total_entered = agg.total_entered;
        throughput = agg.throughput;
        utilization = agg.utilization;
        post_blocks = agg.post_blocks;
        avg_buffer = agg.avg_buffer;
        max_buffer = agg.max_buffer;
        max_hits = agg.max_hits;
      } else {
        // Compute only throughput/utilization without buffer stats.
        total_steps = 0;
        total_entered = 0;
        for (int b = warmup_blocks; b < static_cast<int>(sim.block_counts.size()); ++b) {
          const int bs = (b < static_cast<int>(sim.block_sizes.size())) ? sim.block_sizes[b] : 0;
          total_steps += bs;
          total_entered += sim.block_counts[b];
        }
        throughput = (total_steps > 0) ? (static_cast<double>(total_entered) / static_cast<double>(total_steps)) : 0.0;
        utilization = (total_steps > 0 && NumberOfPrimary > 0)
          ? (static_cast<double>(total_entered) / (static_cast<double>(total_steps) * static_cast<double>(NumberOfPrimary)))
          : 0.0;
        post_blocks.assign(sim.block_counts.begin() + std::min<int>(warmup_blocks, (int)sim.block_counts.size()), sim.block_counts.end());
      }
      // ACF string (post-warmup blocks)
      const auto acf = meshsim::autocorr(post_blocks, 2);
      {
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(6);
        for (size_t i = 0; i < acf.size(); ++i) {
          if (i) oss << ";";
          oss << acf[i];
        }
        acf_str = oss.str();
      }

      // Primary utilization after each feeder (correct, matches Python)
      util_primary_per_feeder.assign(NumberOfFeeders, 0.0);
      if (total_steps > 0) {
        // Sum arrivals per feeder over post-warmup blocks
        std::vector<long long> feeder_post(NumberOfFeeders, 0);
        for (int b = warmup_blocks; b < static_cast<int>(sim.feeder_primary_blocks.size()); ++b) {
          const auto& v = sim.feeder_primary_blocks[b];
          for (int f = 0; f < NumberOfFeeders && f < static_cast<int>(v.size()); ++f) {
            feeder_post[f] += static_cast<long long>(v[f]);
          }
        }
        // Cumulative arrivals up to feeder f
        long long cum = 0;
        for (int f = 0; f < NumberOfFeeders; ++f) {
          cum += feeder_post[f];
          util_primary_per_feeder[f] = static_cast<double>(cum) / (static_cast<double>(total_steps) * static_cast<double>(NumberOfPrimary));
        }
      }
      {
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(6);
        for (int f = 0; f < NumberOfFeeders; ++f) {
          if (f) oss << ";";
          oss << "F" << (f + 1) << ":" << util_primary_per_feeder[f];
        }
        util_primary_str = oss.str();
      }

      // CI for throughput over post-warmup blocks
      if (ci_enable) {
        const int n_blocks = static_cast<int>(post_blocks.size());
        n_blocks_ci = n_blocks;
        if (n_blocks >= 2) {
          // rates per block = block_counts / block_sizes
          std::vector<double> post_rates;
          post_rates.reserve(n_blocks);
          for (int i = 0; i < n_blocks; ++i) {
            const int b = warmup_blocks + i;
            const int bs = (b < static_cast<int>(sim.block_sizes.size())) ? sim.block_sizes[b] : 0;
            const int bc = post_blocks[i];
            post_rates.push_back(bs > 0 ? (static_cast<double>(bc) / static_cast<double>(bs)) : 0.0);
          }
          double mean_rate = 0.0;
          for (double r : post_rates) mean_rate += r;
          mean_rate /= static_cast<double>(n_blocks);
          double m2 = 0.0;
          for (double r : post_rates) m2 += (r - mean_rate) * (r - mean_rate);
          const double sd = std::sqrt(m2 / static_cast<double>(n_blocks - 1));
          const double tcrit = tcrit_approx(ci_level, n_blocks);
          ci_half_two_sided = (n_blocks > 1) ? (tcrit * sd / std::sqrt(static_cast<double>(n_blocks))) : 0.0;
          const double tcrit_one = tcrit_one_sided_approx(ci_level, n_blocks);
          ci_half_one_sided = (n_blocks > 1) ? (tcrit_one * sd / std::sqrt(static_cast<double>(n_blocks))) : 0.0;
          std::cout << "  CI(" << std::fixed << std::setprecision(2) << ci_level
                    << ") for mean throughput over post-WU blocks: two_sided=±" << std::setprecision(6) << ci_half_two_sided
                    << ", one_sided=±" << ci_half_one_sided
                    << " (n=" << n_blocks << ")\n";
        }
      }

      std::cout << "  Warm-up detected: blocks=" << warmup_blocks
                << ", steps=" << warmup_steps
                << ", rate_avg=" << std::fixed << std::setprecision(6) << warmup_rate_avg
                << ", rate_curr=" << warmup_rate_curr << "\n";

      std::cout << "  Post-warmup throughput=" << std::fixed << std::setprecision(6) << throughput
                << "  utilization=" << utilization * 100.0 << "%\n";

      // Buffer summary string (Avg:Max:Hits)
      if (collect_buffer_stats && !avg_buffer.empty()) {
        const int dirs = BiDirectional ? 2 : 1;
        auto idx = [&](int f, int p, int d) { return (f * NumberOfPrimary + p) * dirs + d; };

        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        for (int f = 0; f < NumberOfFeeders; ++f) {
          if (f) oss << " | ";
          oss << "F" << (f + 1) << "[";
          for (int p = 0; p < NumberOfPrimary; ++p) {
            if (p) oss << "/";
            for (int d = 0; d < dirs; ++d) {
              if (d) oss << "/";
              const int id = idx(f, p, d);
              const double av = (id < static_cast<int>(avg_buffer.size())) ? avg_buffer[id] : 0.0;
              const int mx = (id < static_cast<int>(max_buffer.size())) ? max_buffer[id] : 0;
              const int ht = (id < static_cast<int>(max_hits.size())) ? max_hits[id] : 0;
              oss << std::setprecision(3) << av << ":" << mx << ":" << ht;
            }
          }
          oss << "]";
        }
        buffers_str = oss.str();
      }

    } else {
      // Model-only mode placeholders (match Python behavior)
      util_primary_per_feeder.assign(NumberOfFeeders, 0.0);
      {
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(6);
        for (int f = 0; f < NumberOfFeeders; ++f) {
          if (f) oss << ";";
          oss << "F" << (f + 1) << ":0.000000";
        }
        util_primary_str = oss.str();
      }
    }

    // ---------------- Model ----------------
    double approx_util = 0.0;
    double approx_throughput = 0.0;
    std::string approx_case;
    std::vector<double> approx_u_per_feeder;

    if (run_model) {
      auto ar = meshsim::approximate_model(
        NumberOfFeeders,
        NumberOfPrimary,
        BiDirectional,
        buffer_capacity
      );
      approx_util = ar.utilization;
      approx_throughput = ar.throughput;
      approx_case = ar.case_name;

      approx_u_per_feeder = meshsim::approximate_model_per_feeder(
        NumberOfFeeders,
        NumberOfPrimary,
        BiDirectional,
        buffer_capacity
      );

      std::cout << "  Approximation model [" << approx_case << "]: "
                << "util=" << std::fixed << std::setprecision(2) << approx_util * 100.0 << "%, "
                << "throughput=" << std::setprecision(4) << approx_throughput << " items/step\n";
    }

    std::string approx_util_per_feeder_str;
    std::string util_gap_per_feeder_str;

    if (run_model) {
      {
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(6);
        for (int f = 0; f < static_cast<int>(approx_u_per_feeder.size()); ++f) {
          if (f) oss << ";";
          oss << "F" << (f + 1) << ":" << approx_u_per_feeder[f];
        }
        approx_util_per_feeder_str = oss.str();
      }

      if (run_sim) {
        const int L = std::min<int>(static_cast<int>(approx_u_per_feeder.size()), NumberOfFeeders);
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(6);
        for (int f = 0; f < L; ++f) {
          if (f) oss << ";";
          const double gap = approx_u_per_feeder[f] - util_primary_per_feeder[f];
          oss << "F" << (f + 1) << ":" << (gap >= 0 ? "+" : "") << gap;
        }
        util_gap_per_feeder_str = oss.str();
      }
    }

    // ---------------- CSV output ----------------
    const std::vector<std::string> header = {
      "seed", "NumberOfFeeders", "NumberOfPrimary", "BiDirectional",
      "T", "buffer_capacity", "TotalEntered", "Throughput", "Utilization",
      "BlockSize", "WarmupBlocks", "WarmupSteps", "WarmupAvgAfter",
      "WarmupRateAtK", "ACF", "CI_Level", "CI_TwoSided_HalfWidth", "CI_OneSided_HalfWidth",
      "PrimaryUtilPerFeeder",
      "ApproxUtil_Model", "ApproxUtilPerFeeder_Model", "GapUtilPerFeeder_ModelMinusSim",
      "ApproxThroughput_Model", "ApproxCase",
      "Buffers(Avg:Max:Hits)"
    };

    const bool file_exists = fs::exists(output_csv);
    bool write_header = true;
    if (file_exists) {
      std::error_code ec;
      auto sz = fs::file_size(output_csv, ec);
      if (!ec && sz > 0) write_header = false;
    }

    {
      std::ofstream f(output_csv, std::ios::app);
      if (!f) throw std::runtime_error("Cannot open output_csv for append: " + output_csv);

      if (write_header) {
        for (size_t i = 0; i < header.size(); ++i) {
          if (i) f << ",";
          f << header[i];
        }
        f << "\n";
      }

      // buffer_capacity as a string
      std::ostringstream cap_ss;
      cap_ss << "[[" << join_ints(buffer_capacity[0], ",") << "],[" << join_ints(buffer_capacity[1], ",") << "]]";

      // Values row (CSV-escaped for a couple fields)
      auto csv_escape = [](const std::string& x) {
        bool need_q = false;
        for (char c : x) {
          if (c == ',' || c == '"' || c == '\n' || c == '\r') { need_q = true; break; }
        }
        if (!need_q) return x;
        std::string y;
        y.reserve(x.size() + 2);
        y.push_back('"');
        for (char c : x) {
          if (c == '"') y += "\"\"";
          else y.push_back(c);
        }
        y.push_back('"');
        return y;
      };

      const int bidir_int = BiDirectional ? 1 : 0;

      // Note: TotalEntered/Throughput/Utilization are post-warmup when run_sim=true.
      f << seed << "," << NumberOfFeeders << "," << NumberOfPrimary << "," << bidir_int << ","
        << T << "," << csv_escape(cap_ss.str()) << ","
        << total_entered << "," << std::fixed << std::setprecision(6) << throughput << "," << utilization << ","
        << block_size << "," << warmup_blocks << "," << warmup_steps << ","
        << warmup_rate_avg << "," << warmup_rate_curr << ","
        << csv_escape(acf_str) << ",";

      if (ci_enable) {
        f << ci_level << "," << ci_half_two_sided << "," << ci_half_one_sided << ",";
      } else {
        f << "NA,NA,NA,"; // CI_Level, CI_TwoSided_HalfWidth, CI_OneSided_HalfWidth
      }

      f << csv_escape(util_primary_str) << ",";

      if (run_model) {
        f << approx_util << "," << csv_escape(approx_util_per_feeder_str) << "," << csv_escape(util_gap_per_feeder_str) << ","
          << approx_throughput << "," << csv_escape(approx_case) << ",";
      } else {
        f << ",,,,,"; // placeholders (ApproxUtil_Model, ApproxUtilPerFeeder_Model, Gap..., ApproxThroughput_Model, ApproxCase)
      }

      f << csv_escape(buffers_str) << "\n";
    }

    std::cout << "  Results appended to: " << output_csv << "\n";

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
