
/*
  recheck_candidates.cpp
  High-precision re-simulation of the strongest Phase-3 competitors.

  Phase 3 of buffer_certify measures most competitors with a short run (often the
  minimum budget M), so the reported epsilon_hat = max_c eps_c is taken over many
  noisy estimates and is inflated by that maximisation. This tool re-simulates the
  leading competitors, and the reference allocation x_S, at a long common horizon,
  so that each candidate's standing can be judged on a precise estimate:

    - is mu_c above mu_{x_S} at all?
    - is mu_c above (1+eps) mu_{x_S}, i.e. a genuine counterexample?
    - what is eps_c when the competitor is measured precisely?

  Candidates are taken from a Phase-3 detail CSV: the top --top rows by eps_i and
  the top --top rows by mean (union, de-duplicated), so both the allocations that
  drove epsilon_hat and the allocations that simply look strong are re-examined.

  Build:
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel
  (or: c++ -O2 -std=c++17 -pthread MeshSim.cpp recheck_candidates.cpp -o recheck_candidates)

  Run (example):
    ./recheck_candidates --phase3_csv out/m5_n4_B10_twoway_phase3_detail.csv \
                         --phase1_S_csv out/m5_n4_B10_twoway_phase1_S_detail.csv \
                         --m 5 --n 4 --bidirectional true \
                         --top 50 --blocks 250000 --block_size 250 \
                         --alpha 0.01 --eps 0.001 --threads 0 \
                         --out out/m5_n4_B10_twoway_recheck.csv

  The runs use a seed range (--seed0, default 20,000,000) disjoint from the ones
  buffer_certify consumes, so this is an independent experiment, not a replay.
*/

#include "MeshSim.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ---------------- small helpers ----------------

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool in_q = false;
  for (char ch : line) {
    if (ch == '"') { in_q = !in_q; continue; }
    if (ch == ',' && !in_q) { out.push_back(cur); cur.clear(); continue; }
    cur.push_back(ch);
  }
  out.push_back(cur);
  return out;
}

// "(0, 0, 1, 2, 4)" -> {0,0,1,2,4}
std::vector<int> parse_vec(const std::string& s) {
  std::vector<int> v;
  std::string cur;
  for (char ch : s) {
    if (std::isdigit(static_cast<unsigned char>(ch)) || ch == '-') cur.push_back(ch);
    else if (!cur.empty()) { v.push_back(std::stoi(cur)); cur.clear(); }
  }
  if (!cur.empty()) v.push_back(std::stoi(cur));
  return v;
}

std::string vec_to_string(const std::vector<int>& v) {
  std::ostringstream oss;
  oss << "(";
  for (size_t i = 0; i < v.size(); ++i) { if (i) oss << ", "; oss << v[i]; }
  oss << ")";
  return oss.str();
}

double normal_inv_cdf(double p) {
  if (p <= 0.0) return -std::numeric_limits<double>::infinity();
  if (p >= 1.0) return  std::numeric_limits<double>::infinity();
  static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                              -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00};
  static const double d[4] = { 7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
                               3.754408661907416e+00};
  const double plow = 0.02425, phigh = 1.0 - plow;
  if (p < plow) {
    const double q = std::sqrt(-2.0 * std::log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
  }
  if (p > phigh) {
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
  }
  const double q = p - 0.5, r = q * q;
  return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q /
         (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
}

// Smallest eps >= 0 with mu_S(1+eps) - mu_C >= z*sqrt((1+eps)^2 var_S/m_S + var_C/m_C).
// Identical to the rule buffer_certify applies in Phase 3.
double smallest_eps_one_sided(double muS, double varS, long long mS,
                              double muC, double varC, long long mC, double z) {
  const double inf = std::numeric_limits<double>::infinity();
  if (!(muS > 0.0)) return inf;
  const double aS = varS / static_cast<double>(std::max<long long>(1, mS));
  const double aC = varC / static_cast<double>(std::max<long long>(1, mC));
  auto gap = [&](double e) {
    const double u = 1.0 + e;
    return muS * u - muC - z * std::sqrt(u * u * aS + aC);
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

// ---------------- one long run ----------------

struct Estimate {
  double mean = 0.0;
  double var  = 0.0;
  long long m_eff = 0;
  int k_warm = 0;
};

Estimate run_once(int m, int n, bool bidirectional, int blocks, int block_size,
                  const std::vector<int>& c_f, const std::vector<int>& c_b, int seed) {
  std::vector<std::vector<int>> cap{c_f, c_b};
  auto sim = meshsim::run_simulation(m, n, bidirectional, blocks * block_size, cap, seed,
                                     /*Verbal=*/false, block_size, /*collect_buffer_stats=*/false);
  const int k_warm = meshsim::detect_warmup_from_blocks(sim.block_counts, sim.block_sizes, /*strict=*/true);

  Estimate e;
  e.k_warm = k_warm;
  const size_t nb = sim.block_counts.size();
  if (static_cast<size_t>(k_warm) >= nb) return e;

  long double s = 0.0L;
  std::vector<double> u;
  u.reserve(nb - k_warm);
  for (size_t i = k_warm; i < nb; ++i) {
    const double bs = static_cast<double>(sim.block_sizes[i]);
    const double val = (bs > 0.0) ? (sim.block_counts[i] / (bs * n)) : 0.0;
    u.push_back(val);
    s += val;
  }
  e.m_eff = static_cast<long long>(u.size());
  if (e.m_eff == 0) return e;
  e.mean = static_cast<double>(s / static_cast<long double>(u.size()));
  if (e.m_eff > 1) {
    long double q = 0.0L;
    for (double v : u) { const long double d = v - e.mean; q += d * d; }
    e.var = static_cast<double>(q / static_cast<long double>(u.size() - 1));
  }
  return e;
}

// ---------------- CLI ----------------

std::unordered_map<std::string, std::string> parse_args(int argc, char** argv) {
  std::unordered_map<std::string, std::string> out;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--", 0) == 0) {
      std::string key = a.substr(2), val = "true";
      if (i + 1 < argc) { std::string nxt = argv[i + 1]; if (nxt.rfind("--", 0) != 0) { val = nxt; ++i; } }
      out[key] = val;
    }
  }
  return out;
}
bool as_bool(const std::string& s) {
  std::string t; for (char c : s) t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  return (t == "1" || t == "true" || t == "t" || t == "yes" || t == "y");
}
std::string need(const std::unordered_map<std::string,std::string>& a, const std::string& k) {
  auto it = a.find(k);
  if (it == a.end()) throw std::runtime_error("Missing required --" + k);
  return it->second;
}
std::string opt(const std::unordered_map<std::string,std::string>& a, const std::string& k, const std::string& d) {
  auto it = a.find(k); return (it == a.end()) ? d : it->second;
}

struct Candidate {
  int k = -1;
  std::vector<int> c_f, c_b;
  double phase3_mean = 0.0, phase3_eps = 0.0;
  long long phase3_m_eff = 0;
  Estimate est;
};

} // namespace

int main(int argc, char** argv) {
  try {
    auto args = parse_args(argc, argv);

    const int m = std::stoi(need(args, "m"));
    const int n = std::stoi(need(args, "n"));
    const bool bidirectional = as_bool(opt(args, "bidirectional", "false"));

    const std::string phase3_csv = need(args, "phase3_csv");
    const int top        = std::stoi(opt(args, "top", "50"));
    const int blocks     = std::stoi(opt(args, "blocks", "250000"));
    const int block_size = std::stoi(opt(args, "block_size", "250"));
    const double alpha   = std::stod(opt(args, "alpha", "0.01"));
    const double eps     = std::stod(opt(args, "eps", "0.001"));
    const int seed0      = std::stoi(opt(args, "seed0", "20000000"));
    int threads          = std::stoi(opt(args, "threads", "0"));
    const std::string out_csv = opt(args, "out", "recheck_results.csv");

    if (threads <= 0) threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    if (blocks <= 1 || block_size <= 0) throw std::runtime_error("--blocks must be > 1 and --block_size > 0");
    if (static_cast<long long>(blocks) * block_size > 2000000000LL)
      throw std::runtime_error("blocks * block_size exceeds the simulator's 32-bit horizon; lower --blocks");

    // ---- reference allocation x_S ----
    std::vector<int> s_f, s_b;
    if (args.count("s_forward")) {
      s_f = parse_vec(args["s_forward"]);
      s_b = args.count("s_backward") ? parse_vec(args["s_backward"]) : std::vector<int>(s_f.size(), 0);
    } else {
      const std::string p1 = need(args, "phase1_S_csv");
      std::ifstream f(p1);
      if (!f) throw std::runtime_error("Cannot open " + p1);
      std::string line; std::getline(f, line);           // header
      double best = -1.0;
      while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto t = split_csv_line(line);
        if (t.size() < 8) continue;
        const double mu = std::stod(t[3]);
        if (mu > best) { best = mu; s_f = parse_vec(t[6]); s_b = parse_vec(t[7]); }
      }
      if (best < 0.0) throw std::runtime_error("No rows found in " + p1);
      std::cout << "x_S taken from " << p1 << " (largest Phase-1 mean " << std::fixed
                << std::setprecision(7) << best << ")\n";
    }

    // ---- candidates from the Phase-3 detail CSV ----
    std::vector<Candidate> all;
    {
      std::ifstream f(phase3_csv);
      if (!f) throw std::runtime_error("Cannot open " + phase3_csv);
      std::string line; std::getline(f, line);           // header
      while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto t = split_csv_line(line);
        if (t.size() < 9) continue;
        Candidate c;
        c.k            = std::stoi(t[0]);
        c.phase3_mean  = std::stod(t[3]);
        c.phase3_m_eff = std::stoll(t[5]);
        c.phase3_eps   = std::stod(t[6]);
        c.c_f          = parse_vec(t[7]);
        c.c_b          = parse_vec(t[8]);
        all.push_back(std::move(c));
      }
    }
    if (all.empty()) throw std::runtime_error("No competitor rows in " + phase3_csv);

    std::vector<int> by_eps(all.size()), by_mean(all.size());
    for (size_t i = 0; i < all.size(); ++i) by_eps[i] = by_mean[i] = static_cast<int>(i);
    const int take = std::min<int>(top, static_cast<int>(all.size()));
    std::partial_sort(by_eps.begin(), by_eps.begin() + take, by_eps.end(),
                      [&](int a, int b){ return all[a].phase3_eps > all[b].phase3_eps; });
    std::partial_sort(by_mean.begin(), by_mean.begin() + take, by_mean.end(),
                      [&](int a, int b){ return all[a].phase3_mean > all[b].phase3_mean; });

    std::unordered_set<int> chosen;
    std::vector<Candidate> cand;
    for (int i = 0; i < take; ++i) {
      for (int idx : {by_eps[i], by_mean[i]}) {
        if (chosen.insert(idx).second) cand.push_back(all[idx]);
      }
    }

    std::cout << "candidates: " << cand.size() << " (union of top " << take
              << " by eps_i and top " << take << " by Phase-3 mean, out of " << all.size() << ")\n"
              << "horizon   : " << blocks << " blocks x " << block_size << " steps = "
              << static_cast<long long>(blocks) * block_size << " time steps per allocation\n"
              << "threads   : " << threads << "\n" << std::flush;

    // ---- simulate: x_S first, then every candidate ----
    Estimate S = run_once(m, n, bidirectional, blocks, block_size, s_f, s_b, seed0);
    std::cout << "x_S = forward" << vec_to_string(s_f) << " backward" << vec_to_string(s_b)
              << "  mean " << std::fixed << std::setprecision(7) << S.mean
              << "  var " << std::setprecision(6) << std::scientific << S.var
              << "  batches " << S.m_eff << "\n" << std::defaultfloat << std::flush;
    if (S.m_eff < 2) throw std::runtime_error("x_S run produced no usable batches.");

    std::atomic<size_t> next{0};
    std::atomic<size_t> done{0};
    std::mutex io;
    auto worker = [&]() {
      for (;;) {
        const size_t i = next++;
        if (i >= cand.size()) return;
        cand[i].est = run_once(m, n, bidirectional, blocks, block_size,
                               cand[i].c_f, cand[i].c_b, seed0 + 1 + static_cast<int>(i));
        const size_t d = ++done;
        std::lock_guard<std::mutex> lk(io);
        std::cerr << "\r  simulated " << d << "/" << cand.size() << std::flush;
      }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    std::cerr << "\n";

    // ---- verdicts ----
    const double z = normal_inv_cdf(1.0 - alpha);
    std::sort(cand.begin(), cand.end(),
              [](const Candidate& a, const Candidate& b){ return a.est.mean > b.est.mean; });

    fs::path op(out_csv);
    if (op.has_parent_path()) fs::create_directories(op.parent_path());
    std::ofstream f(out_csv);
    f << "k,cfg_forward,cfg_backward,phase3_mean,phase3_m_eff,phase3_eps_i,"
      << "recheck_mean,recheck_var,recheck_m_eff,rel_gap_to_xS,z_vs_xS,z_vs_xS_scaled,"
      << "beats_xS,counterexample,eps_c_precise\n";

    int n_beats = 0, n_counter = 0;
    double eps_max = 0.0;
    std::cout << "\n  k        forward        backward        recheck mean   rel. gap    z      eps_c\n";
    for (const auto& c : cand) {
      const double se  = std::sqrt(S.var / static_cast<double>(S.m_eff) +
                                   c.est.var / static_cast<double>(std::max<long long>(1, c.est.m_eff)));
      const double se2 = std::sqrt((1.0 + eps) * (1.0 + eps) * S.var / static_cast<double>(S.m_eff) +
                                   c.est.var / static_cast<double>(std::max<long long>(1, c.est.m_eff)));
      const double gap = (S.mean > 0.0) ? (c.est.mean - S.mean) / S.mean : 0.0;
      const double zb  = (se  > 0.0) ? (c.est.mean - S.mean) / se : 0.0;
      const double zc  = (se2 > 0.0) ? (c.est.mean - (1.0 + eps) * S.mean) / se2 : 0.0;
      const bool beats = zb > z;
      const bool counter = zc > z;
      const double epsc = smallest_eps_one_sided(S.mean, S.var, S.m_eff,
                                                 c.est.mean, c.est.var, c.est.m_eff, z);
      if (beats) ++n_beats;
      if (counter) ++n_counter;
      if (std::isfinite(epsc)) eps_max = std::max(eps_max, epsc);

      f << c.k << ",\"" << vec_to_string(c.c_f) << "\",\"" << vec_to_string(c.c_b) << "\","
        << std::setprecision(10) << c.phase3_mean << "," << c.phase3_m_eff << ","
        << std::setprecision(10) << c.phase3_eps << ","
        << std::setprecision(10) << c.est.mean << "," << std::setprecision(10) << c.est.var << ","
        << c.est.m_eff << "," << std::setprecision(10) << gap << ","
        << std::setprecision(6) << zb << "," << std::setprecision(6) << zc << ","
        << (beats ? 1 : 0) << "," << (counter ? 1 : 0) << ","
        << std::setprecision(10) << epsc << "\n";

      std::cout << std::setw(6) << c.k << "  " << std::setw(14) << vec_to_string(c.c_f)
                << " " << std::setw(14) << vec_to_string(c.c_b)
                << "  " << std::fixed << std::setprecision(7) << c.est.mean
                << "  " << std::setw(9) << std::setprecision(5) << 100.0 * gap << "%"
                << "  " << std::setw(6) << std::setprecision(2) << zb
                << "  " << std::setw(9) << std::setprecision(6) << epsc
                << (counter ? "   <-- counterexample" : (beats ? "   <-- above x_S" : "")) << "\n";
    }

    std::cout << "\n=== SUMMARY ===\n"
              << "x_S mean                      : " << std::setprecision(7) << S.mean
              << "  (" << S.m_eff << " batches)\n"
              << "candidates above x_S          : " << n_beats << " / " << cand.size()
              << "   (one-sided, alpha=" << std::defaultfloat << alpha << std::fixed << ")\n"
              << "candidates above (1+eps) x_S  : " << n_counter << " / " << cand.size()
              << "   <- genuine counterexamples to the conjecture\n"
              << "max eps_c among candidates    : " << std::setprecision(6) << eps_max
              << "   (epsilon_hat restricted to these allocations, measured at "
              << blocks << " batches)\n"
              << "results written to            : " << out_csv << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n"
              << "Usage:\n  ./recheck_candidates --phase3_csv <file> --phase1_S_csv <file> "
                 "--m 5 --n 4 --bidirectional true [--top 50] [--blocks 250000] "
                 "[--block_size 250] [--alpha 0.01] [--eps 0.001] [--threads 0] [--out <file>]\n";
    return 1;
  }
}
