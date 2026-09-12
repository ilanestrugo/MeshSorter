// -----------------------------------------------------------------------------
//  sweep_rep - evaluate a set of buffer allocations by independent replications
//
//  certify_rep answers whether the structured class contains the best allocation,
//  and to do that it must enumerate every feasible allocation.  That is only
//  possible for small systems: fifteen feeders and a budget of twenty admit
//  about 1.1 billion allocations.  This program evaluates a stated set instead,
//  by default the structured class itself, which is what a designer following
//  the rules of Section 5.1 would actually search.
//
//  It writes the same CSV as certify_rep --dump, so approx_eval.py reads the
//  output of either without knowing which produced it.  Each allocation is run
//  twice, on the same two disjoint stream families certify_rep uses, so the two
//  columns of an allocation are independent and their difference measures the
//  simulation noise.
//
//  Run:    ./sweep_rep -n 15 -m 15 -B 15 --dual --dump out.csv
//  Help:   ./sweep_rep --help
// -----------------------------------------------------------------------------

#include "meshsorter_core.hpp"
#include <atomic>
#include <fstream>
#include <iomanip>

//  An allocation is the vector of capacities along one primary belt: m entries
//  for the single-drop system, and 2m for the dual-drop system with the forward
//  capacity of feeder j at index 2j and the backward capacity at 2j+1.  Feeder 1
//  is never blocked and always receives zero.
struct Alloc { std::vector<int> c; bool structured = true; };

//  Every partition of B into at most k parts, each returned in nondecreasing
//  order.  These are exactly the capacity vectors the design rules admit.
static void partitions(int B, int k, int cap, std::vector<int> &cur,
                       std::vector<std::vector<int>> &out) {
    if (k == 0) { if (B == 0) out.push_back(cur); return; }
    for (int v = 0; v <= std::min(B, cap); v++) {
        cur.push_back(v);
        partitions(B - v, k - 1, v, cur, out);   // nonincreasing as built
        cur.pop_back();
    }
}

static std::vector<Alloc> enumerate_structured(int m, int B, bool dual) {
    std::vector<std::vector<int>> parts;
    std::vector<int> cur;
    partitions(B, m - 1, B, cur, parts);
    std::vector<Alloc> out;
    out.reserve(parts.size());
    for (auto v : parts) {
        std::reverse(v.begin(), v.end());        // nondecreasing along the feeders
        Alloc a;
        a.c.assign((size_t)(dual ? 2 : 1) * m, 0);
        for (int j = 1; j < m; j++) {
            if (dual) a.c[(size_t)2 * j + 1] = v[(size_t)(j - 1)];
            else      a.c[(size_t)j]         = v[(size_t)(j - 1)];
        }
        out.push_back(std::move(a));
    }
    return out;
}

struct Est { double mean = 0.0, var = 0.0; int reps = 0; };

static Est summarize(const std::vector<double> &y) {
    Est e; e.reps = (int)y.size();
    for (double v : y) e.mean += v;
    e.mean /= e.reps;
    for (double v : y) e.var += (v - e.mean) * (v - e.mean);
    e.var /= (e.reps - 1);
    return e;
}

struct Runner {
    Config base;
    Layout ly0;

    void apply(Config &cfg, const Alloc &a) const {
        const int n = cfg.n, m = cfg.m;
        cfg.bufF.assign((size_t)m * n, 0);
        cfg.bufB.assign((size_t)m * n, 0);
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++) {
                if (cfg.dual) {
                    cfg.bufF[(size_t)j * n + i] = a.c[(size_t)2 * j];
                    cfg.bufB[(size_t)j * n + i] = a.c[(size_t)2 * j + 1];
                } else {
                    cfg.bufF[(size_t)j * n + i] = a.c[(size_t)j];
                }
            }
    }

    std::vector<double> run(const Alloc &a, uint64_t streamBase, int count) const {
        Config cfg = base;
        apply(cfg, a);
        Layout ly = buildLayout(cfg);
        long long warm = (base.warmup > 0)
                       ? base.warmup
                       : std::max(20000LL, 10LL * (ly.cmax + 1) * ly.Lmax);
        std::vector<double> Y((size_t)count, 0.0);
        for (int k = 0; k < count; k++) {
            uint64_t s = streamBase + 0x1000193ULL * (uint64_t)(k + 1);
            Replication rep(ly, warm, base.steps, splitmix64(s));
            Y[(size_t)k] = rep.run();
        }
        return Y;
    }
};

//  Allocations are spread over the pool rather than the replications of any one
//  of them, for the same reason as in certify_rep: there are many allocations
//  and few replications each.
static void sweep(const Runner &run, const std::vector<Alloc> &A,
                  std::vector<Est> &out, uint64_t stream, int reps,
                  int threads, const char *label) {
    std::atomic<size_t> next{0}, done{0};
    std::vector<std::thread> pool;
    const int W = (int)std::max<size_t>(1, std::min<size_t>((size_t)threads, A.size()));
    pool.reserve((size_t)W);
    for (int w = 0; w < W; w++)
        pool.emplace_back([&]() {
            for (;;) {
                const size_t i = next.fetch_add(1);
                if (i >= A.size()) return;
                out[i] = summarize(run.run(A[i], stream + 1000003ULL * i, reps));
                const size_t d = done.fetch_add(1) + 1;
                if (label && d % 100 == 0)
                    std::fprintf(stderr, "  %s: %zu/%zu\n", label, d, A.size());
            }
        });
    for (auto &t : pool) t.join();
}

static void usage() {
    std::printf(R"(sweep_rep - evaluate buffer allocations by independent replications

  -n, --belts N          primary belts                                (required)
  -m, --feeders M        feeder loops                                 (required)
  -B, --budget B         per-primary-belt buffer budget               (required)
      --dual             dual-drop                                    (default)
      --single           single-drop
      --spacing S        slots between consecutive primary belts along
                         a feeder loop                                (default 4)
      --turn E           slots in each end curve of the loop          (default 4)
  -L, --loop LEN         set the loop length directly
      --df F             spacing of feeders along a primary belt      (default 4)
      --width W          forward-to-backward distance on a belt       (default 2)
  -T, --steps T          measured time steps per replication     (default 1330000)
  -R, --reps R           replications in each of the two passes       (default 10)
  -w, --warmup W         warm-up steps discarded from each replication
                         (default max(20000, 10*(c+1)*L))
      --seed S           base seed; the two passes use disjoint streams
  -t, --threads K        worker threads                    (default hardware - 2)
      --dump FILE        write every allocation to FILE as CSV        (required)
      --json             machine-readable summary
      --verbose          report progress
  -h, --help             this text

The allocations evaluated are the structured class of Section 5.1: capacities
nondecreasing along the feeders, and, in the dual-drop system, nothing at a
forward drop point.  Their number is the partitions of B into at most m-1 parts,
which stays manageable where the full allocation space does not.
)");
}

int main(int argc, char **argv) {
    Config cfg;
    cfg.steps = 1330000;
    cfg.dp = 8;
    cfg.turn = 4;
    int B = -1, R = 10;
    bool json = false, verbose = false;
    std::string dumpPath;
    try {
        for (int a = 1; a < argc; a++) {
            std::string k = argv[a];
            auto need = [&]() -> std::string {
                if (a + 1 >= argc) throw std::runtime_error(k + " needs a value");
                return std::string(argv[++a]);
            };
            if      (k == "-h" || k == "--help") { usage(); return 0; }
            else if (k == "-n" || k == "--belts")    cfg.n = std::stoi(need());
            else if (k == "-m" || k == "--feeders")  cfg.m = std::stoi(need());
            else if (k == "-B" || k == "--budget")   B = std::stoi(need());
            else if (k == "--dual")                  cfg.dual = true;
            else if (k == "--single")                cfg.dual = false;
            else if (k == "--spacing")               cfg.dp = 2 * std::stoi(need());
            else if (k == "--turn")                  cfg.turn = std::stoi(need());
            else if (k == "-L" || k == "--loop")     cfg.Lbase = std::stoll(need());
            else if (k == "--df")                    cfg.df = std::stoi(need());
            else if (k == "--width")                 cfg.width = std::stoi(need());
            else if (k == "-T" || k == "--steps")    cfg.steps = std::stoll(need());
            else if (k == "-R" || k == "--reps")     R = std::stoi(need());
            else if (k == "-w" || k == "--warmup")   cfg.warmup = std::stoll(need());
            else if (k == "--seed")                  cfg.seed = std::stoull(need());
            else if (k == "-t" || k == "--threads")  cfg.threads = std::stoi(need());
            else if (k == "--dump")                  dumpPath = need();
            else if (k == "--json")                  json = true;
            else if (k == "--verbose")               verbose = true;
            else throw std::runtime_error("unknown option " + k);
        }
        if (cfg.n < 1 || cfg.m < 2 || B < 0) { usage(); return 2; }
        if (R < 3) throw std::runtime_error("--reps must be at least 3");
        if (dumpPath.empty()) throw std::runtime_error("--dump is required");
        if (cfg.threads <= 0) {
            const unsigned hw = std::thread::hardware_concurrency();
            cfg.threads = (hw > 3u) ? (int)hw - 2 : 1;
        }

        Runner run;
        run.base = cfg;
        { Config c0 = cfg; c0.bufF.assign((size_t)cfg.m * cfg.n, 0);
          c0.bufB = c0.bufF; run.ly0 = buildLayout(c0); }

        std::vector<Alloc> A = enumerate_structured(cfg.m, B, cfg.dual);
        if (A.empty()) throw std::runtime_error("no allocations to evaluate");

        //  The same two stream families certify_rep uses, so that a cell
        //  produced here and a cell produced there are directly comparable.
        const uint64_t STREAM1 = cfg.seed;
        const uint64_t STREAM3 = cfg.seed ^ 0xD1B54A32D192ED03ULL;
        std::vector<Est> e1(A.size()), e3(A.size());
        sweep(run, A, e1, STREAM1, R, cfg.threads, verbose ? "pass 1" : nullptr);
        sweep(run, A, e3, STREAM3, R, cfg.threads, verbose ? "pass 2" : nullptr);

        size_t best = 0;
        for (size_t i = 0; i < A.size(); i++)
            if (e3[i].mean > e3[best].mean) best = i;

        {
            std::ofstream out(dumpPath);
            if (!out) throw std::runtime_error("cannot write " + dumpPath);
            out << "n,m,B,dual,c_forward,c_backward,structured,"
                   "r1,mean1,sd1,r3,mean3,sd3\n";
            out.setf(std::ios::fixed);
            for (size_t i = 0; i < A.size(); i++) {
                out << cfg.n << ',' << cfg.m << ',' << B << ','
                    << (cfg.dual ? 1 : 0) << ',';
                for (int j = 0; j < cfg.m; j++)
                    out << (j ? "|" : "")
                        << (cfg.dual ? A[i].c[(size_t)2 * j] : A[i].c[(size_t)j]);
                out << ',';
                if (cfg.dual)
                    for (int j = 0; j < cfg.m; j++)
                        out << (j ? "|" : "") << A[i].c[(size_t)2 * j + 1];
                out << ",1," << e1[i].reps << ','
                    << std::setprecision(6) << e1[i].mean << ','
                    << std::setprecision(6) << std::sqrt(e1[i].var) << ','
                    << e3[i].reps << ','
                    << std::setprecision(6) << e3[i].mean << ','
                    << std::setprecision(6) << std::sqrt(e3[i].var) << '\n';
            }
            if (!out) throw std::runtime_error("error while writing " + dumpPath);
        }

        std::string bs;
        for (int j = 0; j < cfg.m; j++) {
            if (j) bs += ",";
            bs += std::to_string(cfg.dual ? A[best].c[(size_t)2 * j + 1]
                                          : A[best].c[(size_t)j]);
        }
        if (json)
            std::printf("{\"n\":%d,\"m\":%d,\"B\":%d,\"dual\":%s,\"loop\":%lld,"
                        "\"allocations\":%zu,\"reps\":%d,\"steps\":%lld,"
                        "\"best\":\"%s\",\"throughput\":%.6f}\n",
                        cfg.n, cfg.m, B, cfg.dual ? "true" : "false",
                        run.ly0.L[0], A.size(), R, cfg.steps, bs.c_str(),
                        e3[best].mean);
        else
            std::printf("sweep_rep: n=%d m=%d B=%d %s, %zu structured allocations,"
                        " %d replications in each of two passes\n"
                        "  best %s, throughput %.6f\n",
                        cfg.n, cfg.m, B, cfg.dual ? "dual-drop" : "single-drop",
                        A.size(), R, bs.c_str(), e3[best].mean);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "sweep_rep: %s\n", e.what());
        return 1;
    }
    return 0;
}
