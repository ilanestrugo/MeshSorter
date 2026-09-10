// -----------------------------------------------------------------------------
//  certify_rep: does a structured class of buffer allocations contain the best
//  allocation?
//
//  This is the replication-based form of the three-phase certification
//  procedure of Supplement S2.  The earlier version measured simulation effort
//  in batches of a single long run and had to argue that successive batches were
//  approximately uncorrelated.  Here the unit of effort is an independent
//  replication, so the observations entering every variance estimate and every
//  comparison are independent by construction and no such argument is needed.
//
//  Build:  c++ -O2 -std=c++17 -pthread -o certify_rep certify_rep.cpp
//  Run:    ./certify_rep -n 4 -m 4 -B 10 --dual
//  Help:   ./certify_rep --help
//
//  THE QUESTION
//  ------------
//  Fix a system and a per-primary-belt budget B.  An allocation assigns the B
//  units among the drop points along one primary belt; the same assignment is
//  replicated on all n belts, so the installed capacity is nB.  Feeder 1 is
//  never blocked, so it receives nothing.  The feasible allocations split into a
//  structured class S, defined by the design rules of the manuscript, and its
//  complement C.  Writing mu_x for the steady-state throughput under x and
//  x_S* for the best allocation in S, we test
//
//      H0: some c in C has mu_c >= (1+eps) mu_{x_S*}
//      H1: every c in C has mu_c <  (1+eps) mu_{x_S*}
//
//  and report the smallest eps for which no competitor is supported by the data.
//
//  THE THREE PHASES
//  ----------------
//  Phase 1 (pilot).  Every allocation, in S and in C alike, is run for R0
//    independent replications.  Each yields a mean and a variance across those
//    replications.  The allocation of S with the largest mean is the reference
//    x_hat_S.
//
//  Phase 2 (plan).  For each competitor c, the number of replications R(c) is
//    the smallest that makes the one-sided comparison against x_hat_S
//    informative at the target indifference zone eps and level alpha:
//
//      t * sqrt( (1+eps)^2 s2_S / R_S  +  s2_c / R(c) )  <=  mu_S (1+eps) - mu_c
//
//    solved for R(c), rounded up, multiplied by a safety factor kappa that
//    guards against a pilot underestimate of s2_c, and floored at R0.  A
//    competitor whose planned effort exceeds that of the reference triggers an
//    extension of the reference by R0 further replications, after which every
//    R(c) is recomputed; the loop ends when no competitor needs more effort than
//    the reference.  Effort here is a count of replications, where the earlier
//    procedure counted batches.
//
//  Phase 3 (validate).  Everything is simulated again on a disjoint set of
//    random number streams, the reference for R_S replications and each
//    competitor for its planned R(c).  For each competitor the smallest eps
//    satisfying the inequality above is found by bisection, and the largest of
//    those is reported.  Running the validation on fresh streams is what keeps
//    the selection of x_hat_S in Phase 1 from biasing the conclusion.
//
//  All comparisons use Welch's t with Satterthwaite degrees of freedom rather
//  than a normal quantile, because the number of replications behind a variance
//  estimate is tens rather than thousands.
// -----------------------------------------------------------------------------

#include "meshsorter_core.hpp"
#include <atomic>
#include <mutex>
#include <numeric>

// ------------------------------------------------------------- allocations ---
//
//  An allocation is the vector of capacities along ONE primary belt.  In the
//  single-drop system it has m entries, one per feeder.  In the dual-drop system
//  it has 2m, the forward capacity of feeder j at index 2j and the backward
//  capacity at 2j+1.  Feeder 1 always receives zero.

struct Alloc {
    std::vector<int> c;
    bool structured = false;
};

//  All ways of writing B as an ordered sum of k non-negative integers.
static void compositions(int B, int k, std::vector<int> &cur,
                         std::vector<std::vector<int>> &out) {
    if (k == 1) { cur.push_back(B); out.push_back(cur); cur.pop_back(); return; }
    for (int v = 0; v <= B; v++) {
        cur.push_back(v);
        compositions(B - v, k - 1, cur, out);
        cur.pop_back();
    }
}

//  Enumerate the feasible allocations and mark the structured ones.
//    single-drop: S is the set of nondecreasing vectors (design rule 1)
//    dual-drop:   S is the set with no forward capacity and nondecreasing
//                 backward capacity (design rule 2)
static std::vector<Alloc> enumerate_allocs(int m, int B, bool dual) {
    const int slots = (dual ? 2 : 1) * (m - 1);       // feeder 1 gets nothing
    std::vector<std::vector<int>> comps;
    std::vector<int> cur;
    compositions(B, slots, cur, comps);
    std::vector<Alloc> out;
    out.reserve(comps.size());
    for (const auto &v : comps) {
        Alloc a;
        a.c.assign((size_t)(dual ? 2 : 1) * m, 0);
        for (int j = 1; j < m; j++) {
            if (dual) {
                a.c[(size_t)2 * j]     = v[(size_t)2 * (j - 1)];
                a.c[(size_t)2 * j + 1] = v[(size_t)2 * (j - 1) + 1];
            } else {
                a.c[(size_t)j] = v[(size_t)(j - 1)];
            }
        }
        bool ok = true;
        if (dual) {
            for (int j = 0; j < m && ok; j++) if (a.c[(size_t)2 * j] != 0) ok = false;
            for (int j = 1; j + 1 < m && ok; j++)
                if (a.c[(size_t)2 * j + 1] > a.c[(size_t)2 * (j + 1) + 1]) ok = false;
        } else {
            for (int j = 1; j + 1 < m && ok; j++)
                if (a.c[(size_t)j] > a.c[(size_t)j + 1]) ok = false;
        }
        a.structured = ok;
        out.push_back(std::move(a));
    }
    return out;
}

//  Backward capacities alone when there is no forward capacity, which is the
//  convention the manuscript uses for the structured class; otherwise the full
//  vector as forward/backward pairs, so an outside allocation is never printed
//  in a form that hides part of its budget.
static std::string alloc_str(const Alloc &a, int m, bool dual) {
    std::string s;
    if (!dual) {
        for (int j = 0; j < m; j++) { if (j) s += ","; s += std::to_string(a.c[(size_t)j]); }
        return s;
    }
    bool anyF = false;
    for (int j = 0; j < m; j++) if (a.c[(size_t)2 * j]) anyF = true;
    for (int j = 0; j < m; j++) {
        if (j) s += ",";
        if (anyF) s += std::to_string(a.c[(size_t)2 * j]) + "/"
                     + std::to_string(a.c[(size_t)2 * j + 1]);
        else      s += std::to_string(a.c[(size_t)2 * j + 1]);
    }
    return s;
}

// ------------------------------------------------------------- statistics ----

//  Two-sided Student t quantile at level 1-alpha for one-sided use, by
//  bisection on the incomplete beta function.  Only alpha and df vary here, so
//  a small table plus interpolation would do, but this keeps the program exact
//  for any alpha the user asks for.
static double betacf(double a, double b, double x) {
    const int MAXIT = 300; const double EPS = 3e-14, FPMIN = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0, c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d; double h = d;
    for (int mI = 1; mI <= MAXIT; mI++) {
        int m2 = 2 * mI;
        double aa = mI * (b - mI) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c; if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; h *= d * c;
        aa = -(a + mI) * (qab + mI) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c; if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; double del = d * c; h *= del;
        if (std::fabs(del - 1.0) < EPS) break;
    }
    return h;
}
static double betai(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double lbeta = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
                 + a * std::log(x) + b * std::log1p(-x);
    return (x < (a + 1.0) / (a + b + 2.0))
         ? std::exp(lbeta) * betacf(a, b, x) / a
         : 1.0 - std::exp(lbeta) * betacf(b, a, 1.0 - x) / b;
}
//  P(T_df > t) for t >= 0
static double t_sf(double t, double df) {
    return 0.5 * betai(0.5 * df, 0.5, df / (df + t * t));
}
//  smallest t with P(T_df > t) <= alpha
static double t_quantile(double alpha, double df) {
    if (df < 1.0) df = 1.0;
    double lo = 0.0, hi = 200.0;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (t_sf(mid, df) > alpha) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}
//  Welch-Satterthwaite degrees of freedom for the difference of two means
static double welch_df(double v1, int r1, double v2, int r2) {
    const double a = v1 / r1, b = v2 / r2;
    const double num = (a + b) * (a + b);
    const double den = a * a / (r1 - 1) + b * b / (r2 - 1);
    return (den > 0.0) ? num / den : 1.0;
}

// ---------------------------------------------------------------- runner -----

struct Est {                       // what a set of replications tells us
    double mean = 0.0, var = 0.0;
    int reps = 0;
};

struct Runner {
    Config base;                   // geometry, run length, warm-up
    Layout ly0;                    // layout without buffers, for Lmax
    int threads = 1;

    //  Fill the per-crossing capacities of `cfg` from an allocation along one
    //  primary belt, replicated across the n belts.
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

    //  Replications `from` .. `from+count-1` of one allocation.  The stream of a
    //  replication depends only on (streamBase, index), never on the thread that
    //  runs it or on how the work was split, so results are reproducible and
    //  Phase 1 and Phase 3 draw from disjoint streams.
    std::vector<double> run(const Alloc &a, uint64_t streamBase,
                            int from, int count) const {
        Config cfg = base;
        apply(cfg, a);
        Layout ly = buildLayout(cfg);
        long long warm = (base.warmup > 0)
                       ? base.warmup
                       : std::max(20000LL, 10LL * (ly.cmax + 1) * ly.Lmax);
        std::vector<double> Y((size_t)count, 0.0);
        std::atomic<int> next{0};
        const int W = std::max(1, std::min(threads, count));
        std::vector<std::thread> pool;
        pool.reserve((size_t)W);
        for (int w = 0; w < W; w++)
            pool.emplace_back([&]() {
                for (;;) {
                    const int k = next.fetch_add(1);
                    if (k >= count) return;
                    uint64_t s = streamBase + 0x1000193ULL * (uint64_t)(from + k + 1);
                    Replication rep(ly, warm, base.steps, splitmix64(s));
                    Y[(size_t)k] = rep.run();
                }
            });
        for (auto &t : pool) t.join();
        return Y;
    }
};

static Est summarize(const std::vector<double> &y) {
    Est e; e.reps = (int)y.size();
    for (double v : y) e.mean += v;
    e.mean /= e.reps;
    for (double v : y) e.var += (v - e.mean) * (v - e.mean);
    e.var /= (e.reps - 1);
    return e;
}

//  The smallest eps >= 0 at which the reference is not beaten by this
//  competitor, at level alpha.  Zero when the data already separate them.
static double implied_gap(const Est &S, const Est &c, double alpha) {
    auto shortfall = [&](double eps) {
        const double df = welch_df((1 + eps) * (1 + eps) * S.var, S.reps, c.var, c.reps);
        const double t = t_quantile(alpha, df);
        const double se = std::sqrt((1 + eps) * (1 + eps) * S.var / S.reps
                                    + c.var / c.reps);
        return S.mean * (1 + eps) - c.mean - t * se;       // >= 0 means separated
    };
    if (shortfall(0.0) >= 0.0) return 0.0;
    double lo = 0.0, hi = 1.0;
    while (shortfall(hi) < 0.0 && hi < 1e6) hi *= 2.0;
    for (int i = 0; i < 100; i++) {
        double mid = 0.5 * (lo + hi);
        if (shortfall(mid) < 0.0) lo = mid; else hi = mid;
    }
    return hi;
}

// ------------------------------------------------------------------ main -----

static void usage() {
    std::printf(R"(certify_rep - replication-based certification of a structured
                buffer-allocation class (Supplement S2)

USAGE
  certify_rep -n BELTS -m FEEDERS -B BUDGET [options]

SYSTEM
  -n, --belts N          primary belts                                (required)
  -m, --feeders M        feeder loops                                 (required)
  -B, --budget B         per-primary-belt buffer budget               (required)
      --dual             dual-drop                                    (default)
      --single           single-drop
      --spacing S        slots between consecutive primary belts along
                         a feeder loop; the loop length is 2*S*N + 2*E (default 4)
      --turn E           slots in each end curve of the loop          (default 4)
  -L, --loop LEN         set the loop length directly, overriding the
                         two lengths above
      --df F             spacing of feeders along a primary belt      (default 4)
      --width W          forward-to-backward distance on a belt       (default 2)

  The loop lengths are common to all feeders.  Section 5.1 of the manuscript
  shows that staggering them is worth almost nothing once buffers are installed,
  so it is not offered here.

DESIGN
  -T, --steps T          measured time steps per replication     (default 1330000)
  -R, --reps R0          replications in the pilot, and the floor on
                         the planned effort of any allocation         (default 10)
      --kappa K          safety factor on the planned replication
                         count, against a pilot underestimate of the
                         competitor variance                          (default 3)
      --eps E            target indifference zone                     (default 0.001)
      --alpha A          one-sided level of every comparison          (default 0.01)
  -w, --warmup W         warm-up steps discarded from each replication
                         (default max(20000, 10*(c+1)*L))
      --seed S           base seed; Phase 1 and Phase 3 use disjoint
                         stream ranges derived from it            (default 20260901)
  -t, --threads K        worker threads                    (default hardware - 2)
      --max-reps R       refuse to plan more than this many replications
                         for one allocation, as a guard                (default 4000)

OUTPUT
      --json             machine-readable summary
      --verbose          report progress as the phases run
  -h, --help             this text
)");
}

int main(int argc, char **argv) {
    Config cfg;
    cfg.steps = 1330000;
    cfg.dp = 8;                    // spacing 4 between consecutive primary belts
    cfg.turn = 4;
    int B = -1, R0 = 10, kappa = 3, maxReps = 4000;
    double eps = 0.001, alpha = 0.01;
    bool json = false, verbose = false;
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
            else if (k == "-R" || k == "--reps")     R0 = std::stoi(need());
            else if (k == "--kappa")                 kappa = std::stoi(need());
            else if (k == "--eps")                   eps = std::stod(need());
            else if (k == "--alpha")                 alpha = std::stod(need());
            else if (k == "-w" || k == "--warmup")   cfg.warmup = std::stoll(need());
            else if (k == "--seed")                  cfg.seed = std::stoull(need());
            else if (k == "-t" || k == "--threads")  cfg.threads = std::stoi(need());
            else if (k == "--max-reps")              maxReps = std::stoi(need());
            else if (k == "--json")                  json = true;
            else if (k == "--verbose")               verbose = true;
            else throw std::runtime_error("unknown option " + k);
        }
        if (cfg.n < 1 || cfg.m < 2 || B < 0) { usage(); return 2; }
        if (R0 < 3) throw std::runtime_error("--reps must be at least 3");
        if (cfg.threads <= 0) {
            const unsigned hw = std::thread::hardware_concurrency();
            cfg.threads = (hw > 3u) ? (int)hw - 2 : 1;
        }

        Runner run;
        run.base = cfg;
        run.threads = cfg.threads;
        { Config c0 = cfg; c0.bufF.assign((size_t)cfg.m * cfg.n, 0);
          c0.bufB = c0.bufF; run.ly0 = buildLayout(c0); }

        std::vector<Alloc> A = enumerate_allocs(cfg.m, B, cfg.dual);
        size_t nS = 0;
        for (const auto &a : A) nS += a.structured;
        const long long warmUsed = (cfg.warmup > 0) ? cfg.warmup
                                 : std::max(20000LL, 10LL * (B + 1) * run.ly0.Lmax);

        if (!json) {
            std::printf("certify_rep: n=%d m=%d B=%d %s, loop length %lld\n",
                        cfg.n, cfg.m, B, cfg.dual ? "dual-drop" : "single-drop",
                        run.ly0.L[0]);
            std::printf("  %zu feasible allocations, %zu structured, %zu competitors\n",
                        A.size(), nS, A.size() - nS);
            std::printf("  pilot %d replications of %lld steps, warm-up at most %lld,"
                        " kappa %d, eps %g, alpha %g, %d threads\n\n",
                        R0, cfg.steps, warmUsed, kappa, eps, alpha, cfg.threads);
        }
        if (nS == 0) throw std::runtime_error("the structured class is empty");

        // ---------------- Phase 1: pilot every allocation ------------------
        const uint64_t STREAM1 = cfg.seed;
        const uint64_t STREAM3 = cfg.seed ^ 0xD1B54A32D192ED03ULL;
        std::vector<Est> e1(A.size());
        for (size_t i = 0; i < A.size(); i++) {
            e1[i] = summarize(run.run(A[i], STREAM1 + 1000003ULL * i, 0, R0));
            if (verbose && (i % 200 == 199))
                std::fprintf(stderr, "  phase 1: %zu/%zu\n", i + 1, A.size());
        }
        size_t ref = 0; double best = -1.0;
        for (size_t i = 0; i < A.size(); i++)
            if (A[i].structured && e1[i].mean > best) { best = e1[i].mean; ref = i; }
        if (!json)
            std::printf("phase 1: reference allocation %s, pilot mean %.6f\n",
                        alloc_str(A[ref], cfg.m, cfg.dual).c_str(), e1[ref].mean);

        // ---------------- Phase 2: plan the effort -------------------------
        int Rs = R0;
        std::vector<double> refY = run.run(A[ref], STREAM1 + 1000003ULL * ref, 0, R0);
        Est eS = summarize(refY);
        std::vector<int> plan(A.size(), 0);
        for (int round = 0; ; round++) {
            int needed = 0;
            for (size_t i = 0; i < A.size(); i++) {
                if (i == ref) { plan[i] = Rs; continue; }
                const double delta = eS.mean * (1 + eps) - e1[i].mean;
                const double df = welch_df((1 + eps) * (1 + eps) * eS.var, Rs,
                                           e1[i].var, R0);
                const double t = t_quantile(alpha, df);
                const double lhs = (delta / t) * (delta / t)
                                 - (1 + eps) * (1 + eps) * eS.var / Rs;
                double r;
                if (delta <= 0.0 || lhs <= 0.0) r = (double)Rs;       // infeasible: cap
                else r = std::min((double)Rs,
                                  (double)kappa * std::ceil(e1[i].var / lhs));
                plan[i] = std::max(R0, (int)std::min((double)maxReps, r));
                needed = std::max(needed, plan[i]);
            }
            if (needed <= Rs || Rs >= maxReps) break;
            const int add = std::min(R0, maxReps - Rs);
            auto extra = run.run(A[ref], STREAM1 + 1000003ULL * ref, Rs, add);
            refY.insert(refY.end(), extra.begin(), extra.end());
            Rs += add;
            eS = summarize(refY);
            if (verbose)
                std::fprintf(stderr, "  phase 2: reference extended to %d replications\n", Rs);
        }
        long long planned = 0;
        for (size_t i = 0; i < A.size(); i++) planned += plan[i];
        if (!json)
            std::printf("phase 2: reference at %d replications, "
                        "validation plans %lld replications in total\n", Rs, planned);

        // ---------------- Phase 3: independent validation ------------------
        std::vector<Est> e3(A.size());
        for (size_t i = 0; i < A.size(); i++) {
            e3[i] = summarize(run.run(A[i], STREAM3 + 1000003ULL * i, 0, plan[i]));
            if (verbose && (i % 200 == 199))
                std::fprintf(stderr, "  phase 3: %zu/%zu\n", i + 1, A.size());
        }
        size_t refv = ref; double bestv = -1.0;
        for (size_t i = 0; i < A.size(); i++)
            if (A[i].structured && e3[i].mean > bestv) { bestv = e3[i].mean; refv = i; }

        double epsHat = 0.0; size_t worst = ref;
        for (size_t i = 0; i < A.size(); i++) {
            if (A[i].structured) continue;
            const double g = implied_gap(e3[ref], e3[i], alpha);
            if (g > epsHat) { epsHat = g; worst = i; }
        }
        const double hw = t_quantile(0.025, Rs - 1)
                        * std::sqrt(e3[ref].var / Rs);

        if (json) {
            std::printf("{\"n\":%d,\"m\":%d,\"B\":%d,\"dual\":%s,\"loop\":%lld,"
                        "\"allocations\":%zu,\"structured\":%zu,\"R0\":%d,"
                        "\"reps_reference\":%d,\"steps\":%lld,\"warmup\":%lld,"
                        "\"eps_target\":%g,\"alpha\":%g,\"kappa\":%d,"
                        "\"reference\":\"%s\",\"throughput\":%.6f,\"halfwidth\":%.6f,"
                        "\"eps_hat\":%.6f,\"worst_competitor\":\"%s\","
                        "\"reference_confirmed\":%s,\"replications_total\":%lld}\n",
                        cfg.n, cfg.m, B, cfg.dual ? "true" : "false", run.ly0.L[0],
                        A.size(), nS, R0, Rs, cfg.steps, warmUsed, eps, alpha, kappa,
                        alloc_str(A[ref], cfg.m, cfg.dual).c_str(),
                        e3[ref].mean, hw, epsHat,
                        alloc_str(A[worst], cfg.m, cfg.dual).c_str(),
                        (refv == ref) ? "true" : "false", planned);
        } else {
            std::printf("phase 3: throughput of the reference %.6f +- %.6f\n",
                        e3[ref].mean, hw);
            std::printf("         smallest supported relative gap  eps_hat = %.6f\n", epsHat);
            if (epsHat > 0.0)
                std::printf("         attained by the outside allocation %s\n",
                            alloc_str(A[worst], cfg.m, cfg.dual).c_str());
            std::printf("         the validation run %s the Phase 1 choice of reference\n",
                        (refv == ref) ? "confirms" : "DOES NOT confirm");
            std::printf("\n%s at eps = %g: %s\n",
                        epsHat <= eps ? "H0 rejected" : "H0 not rejected", eps,
                        epsHat <= eps
                          ? "no allocation outside the structured class is supported"
                          : "at least one outside allocation is supported");
        }
    } catch (const std::exception &e) {
        std::fprintf(stderr, "certify_rep: %s\n", e.what());
        return 1;
    }
    return 0;
}
