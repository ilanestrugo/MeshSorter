// -----------------------------------------------------------------------------
//  meshsorter: a standalone replication-based simulator of the MeshSorter
//
//  Accompanies the manuscript
//      "MeshSorter: A Two-Layer Conveyor Architecture for High-Throughput
//       Sortation"
//
//  The program estimates the steady-state throughput of one MeshSorter
//  configuration by running R statistically independent replications in
//  parallel, one worker thread per replication, discarding a warm-up period
//  from each, and forming an ordinary t confidence interval on the R
//  replication averages.  The output analysis protocol is the one described in
//  Supplement S1 of the manuscript.
//
//  Build:      c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp
//  Run:        ./meshsorter -n 4 -m 4 --dual
//  Help:       ./meshsorter --help
//
//  This file is self-contained: it depends on nothing but the C++17 standard
//  library, and it shares no code with the other programs in this repository.
// -----------------------------------------------------------------------------
//
//  THE MODEL
//  ---------
//  Time is discrete.  Every conveyor advances exactly one slot per time step.
//
//  There are n primary belts, numbered 1..n, and m circular feeder loops,
//  numbered 1..m.  Feeder j has L_j slots.  Slot 0 of a feeder is its loading
//  station.  A feeder crosses every primary belt once on its forward run and,
//  in the dual-drop design, once again on its backward run.
//
//  A crossing is identified by two coordinates:
//
//    q^f_i, q^b_i   the position of the crossing with primary belt i on the
//                   feeder loop, measured in slots from the loading station;
//    s^f_j, s^b_j   the position of feeder j's forward and backward crossings
//                   along a primary belt.
//
//  In this program the loop coordinates q are the same for every feeder, and
//  the belt coordinates s are the same for every primary belt, which is the
//  regular layout analysed in the manuscript.
//
//  One time step consists of three phases, in this order:
//
//    1. Loading.     For every feeder, if slot 0 is empty, a new item is
//                    admitted and given a destination drawn uniformly from the
//                    n primary belts.  A feeder whose slot 0 is occupied
//                    admits nothing; blocked items are not queued outside the
//                    system.
//    2. Transfers.   Every crossing is examined.  See transfer() below.
//    3. Advancement. Every conveyor moves one slot.  An item that reaches the
//                    end of a primary belt leaves the system.
//
//  Throughput is measured as the number of admissions per time step during the
//  measurement window.  In steady state admissions and departures balance, so
//  this is the throughput of the system.
//
//  GEOMETRY
//  --------
//  The default geometry is the one used throughout the manuscript.  It is
//  built from two lengths:
//
//    dp    slots of the feeder loop consumed by each primary belt.  The feeder
//          crosses a belt twice, so consecutive belts sit dp/2 slots apart on
//          the forward run and dp/2 slots apart on the backward run.
//          Default 4, hence a spacing of 2 slots between consecutive belts.
//    turn  slots of the feeder loop consumed by each end of the loop.
//          Default 4.
//
//  giving the base loop length
//
//        L = n * dp + 2 * turn                        ( = 4n + 8 by default )
//
//  and the crossing positions
//
//        q^f_i = 1 + (i-1) * dp/2                     ( = 2i - 1 )
//        q^b_i = L - turn + 1 - (i-1) * dp/2          ( = 4n + 7 - 2i )
//
//  With the defaults the last forward crossing and the first backward crossing
//  of a feeder are dp + turn = 8 slots apart, and the last backward crossing is
//  turn = 4 slots from the loading station, going forward around the loop.
//
//  Along a primary belt, consecutive feeders are df slots apart and the two
//  crossings of one feeder are width slots apart:
//
//        s^f_j = (j-1) * df ,      s^b_j = s^f_j + width
//
//  with df = 4 and width = 2 by default.  A single-drop system has no backward
//  crossings at all.
//
//  STAGGERED LOOP LENGTHS
//  ----------------------
//  With --stagger the loop lengths are
//
//        L_1 = L_2 = L ,   L_j = L + step * (j - 2)   for j >= 3
//
//  The added slots have to go somewhere on the loop, and where they go
//  matters, because it changes the distance a blocked item travels between its
//  forward and its backward crossing with the same belt.  --extra selects the
//  placement:
//
//    turnaround  (default)  the added slots extend the far end of the loop, so
//                           the backward crossings move away from the forward
//                           ones: q^b_i is computed from L_j.  This is the
//                           layout of panel (c) of Table 2.
//    return                 the added slots extend the return run between the
//                           last backward crossing and the loading station, so
//                           q^b_i keeps the value computed from the base
//                           length L.  This is the layout of panel (b).
//
//  In a single-drop system the two placements coincide.
//
//  BUFFERS
//  -------
//  A buffer of capacity c may be placed at any crossing.  It holds items that
//  were ready to leave the feeder but found the primary slot occupied.  See
//  transfer() for the exact discipline.
//
//  OUTPUT ANALYSIS
//  ---------------
//  R replications differ only in their random number streams.  Each discards
//  its first W time steps and then averages the throughput over the next T
//  steps.  With Y_1..Y_R the replication averages, the program reports
//
//        Ybar +- t_{R-1,0.975} * s / sqrt(R) .
//
//  The default warm-up is the rule of Supplement S1,
//
//        W = max( 20000 , 10 * (c_max + 1) * L_max ) ,
//
//  where L_max is the longest feeder loop and c_max the largest buffer
//  capacity anywhere in the configuration.
// -----------------------------------------------------------------------------

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iostream>

// ---------------------------------------------------------------- random ----
//
// Each (replication, feeder) pair gets its own stream.  The stream seed is
// derived from the base seed and the two indices with SplitMix64, and the
// stream itself is a 64-bit Mersenne twister.  Two consequences worth stating:
// the streams do not depend on how the work is distributed over threads, so
// the result is the same for any --threads setting; and rerunning with the
// same --seed reproduces every number exactly, on any platform.

static inline uint64_t splitmix64(uint64_t &x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

struct Rng {                       // xoshiro-free, dependency-free 64-bit LCG-mix
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) { for (int i = 0; i < 8; i++) next(); }
    inline uint64_t next() {       // SplitMix64 as a generator
        return splitmix64(s);
    }
    // Unbiased uniform on {0,...,n-1} by rejection.  Portable and exactly
    // reproducible; std::uniform_int_distribution is not, because its
    // implementation differs between standard libraries.
    inline uint32_t below(uint32_t n) {
        const uint64_t limit = (UINT64_MAX / n) * n;
        uint64_t x;
        do { x = next(); } while (x >= limit);
        return static_cast<uint32_t>(x % n);
    }
};

// ------------------------------------------------------------ parameters ----

struct Config {
    int  n = 0;                    // primary belts
    int  m = 0;                    // feeder loops
    int  dp = 4;                   // loop slots per primary belt (even)
    int  turn = 4;                 // loop slots per end of the loop
    int  df = 4;                   // spacing of feeders along a primary belt
    int  width = 2;                // s^b_j - s^f_j
    bool dual = true;              // dual-drop if true, single-drop if false
    long long Lbase = -1;          // base loop length, -1 = n*dp + 2*turn
    std::vector<long long> loops;  // explicit per-feeder lengths, empty if unused
    bool stagger = false;          // derive L_j from Lbase
    long long stagStep = 1;        // slots added per feeder when staggering
    bool extraTurnaround = true;   // where slots beyond Lbase are inserted

    std::vector<int> bufF, bufB;   // capacity per crossing, indexed [j*n + i]

    long long steps = 4000000;     // measured steps per replication
    long long warmup = -1;         // -1 = the rule of Supplement S1
    int  reps = 10;
    int  threads = 0;              // 0 = one per replication
    uint64_t seed = 20260901ULL;
    bool json = false;
    bool perFeeder = false;
};

// Everything the simulation needs, derived once from Config and shared,
// read-only, by all worker threads.
struct Layout {
    int n, m, P, pmin;
    bool dual;
    std::vector<long long> L;      // L[j]
    std::vector<long long> qf, qb; // [j*n + i]
    std::vector<int> sf, sb;       // [j]
    std::vector<int> cf, cb;       // [j*n + i]
    long long Lmax = 0;
    int cmax = 0;
};

static Layout buildLayout(const Config &cfg) {
    Layout ly;
    ly.n = cfg.n; ly.m = cfg.m; ly.dual = cfg.dual;
    const int n = cfg.n, m = cfg.m;

    const long long Lbase = (cfg.Lbase > 0) ? cfg.Lbase
                                            : (long long)n * cfg.dp + 2LL * cfg.turn;
    ly.L.assign(m, Lbase);
    if (!cfg.loops.empty()) {
        if ((int)cfg.loops.size() != m)
            throw std::runtime_error("--loops needs exactly m values");
        ly.L = cfg.loops;
    } else if (cfg.stagger) {
        for (int j = 0; j < m; j++)
            ly.L[j] = Lbase + cfg.stagStep * std::max(0, j - 1);   // j is 0-based
    }
    for (int j = 0; j < m; j++)
        if (ly.L[j] < Lbase)
            throw std::runtime_error("a feeder loop is shorter than the base length");

    ly.qf.assign((size_t)m * n, 0);
    ly.qb.assign((size_t)m * n, 0);
    const int half = cfg.dp / 2;
    for (int j = 0; j < m; j++) {
        // The backward crossings are anchored either to this feeder's own
        // length (extra slots at the far turnaround) or to the base length
        // (extra slots on the return run).
        const long long anchor = cfg.extraTurnaround ? ly.L[j] : Lbase;
        for (int i = 0; i < n; i++) {
            ly.qf[(size_t)j * n + i] = 1 + (long long)i * half;
            ly.qb[(size_t)j * n + i] = anchor - cfg.turn + 1 - (long long)i * half;
        }
    }

    ly.sf.assign(m, 0); ly.sb.assign(m, 0);
    for (int j = 0; j < m; j++) {
        ly.sf[j] = j * cfg.df;
        ly.sb[j] = ly.sf[j] + (cfg.dual ? cfg.width : 0);
    }
    ly.pmin = ly.sf[0];
    ly.P = ly.sb[m - 1] - ly.pmin + 1;

    ly.cf = cfg.bufF; ly.cb = cfg.bufB;
    if (!cfg.dual) std::fill(ly.cb.begin(), ly.cb.end(), 0);

    for (int j = 0; j < m; j++) ly.Lmax = std::max(ly.Lmax, ly.L[j]);
    for (size_t k = 0; k < ly.cf.size(); k++)
        ly.cmax = std::max(ly.cmax, std::max(ly.cf[k], ly.cb[k]));

    // ---- validation: the layout must be physically meaningful -------------
    for (int j = 0; j < m; j++) {
        if (ly.qf[(size_t)j * n] < 1)
            throw std::runtime_error("the first forward crossing must follow the loading station");
        for (int i = 0; i + 1 < n; i++)
            if (ly.qf[(size_t)j*n+i] >= ly.qf[(size_t)j*n+i+1])
                throw std::runtime_error("forward crossings must increase along the loop");
        if (cfg.dual) {
            if (ly.qf[(size_t)j*n+n-1] >= ly.qb[(size_t)j*n+n-1])
                throw std::runtime_error("the last forward crossing must precede the last backward one");
            for (int i = 0; i + 1 < n; i++)
                if (ly.qb[(size_t)j*n+i] <= ly.qb[(size_t)j*n+i+1])
                    throw std::runtime_error("backward crossings must decrease in the belt index");
            if (ly.qb[(size_t)j*n] >= ly.L[j])
                throw std::runtime_error("the first backward crossing must precede the loading station");
        }
        if (j + 1 < m && ly.sb[j] > ly.sf[j + 1])
            throw std::runtime_error("consecutive feeders overlap on the primary belts");
    }
    if (cfg.dp % 2 != 0) throw std::runtime_error("--dp must be even");
    return ly;
}

// ------------------------------------------------------------ simulation ----
//
//  One replication.  The state is
//
//    feeder[j][.]  L_j cells, 0 for empty and 1..n for an item destined for
//                  that primary belt;
//    prim[i][.]    a circular window of P slots of primary belt i, covering
//                  the positions from the first to the last crossing, 1 if the
//                  slot carries an item;
//    bufF, bufB    the current contents of the buffers, one counter per
//                  crossing.
//
//  Conveyor motion is represented by rotating the index offsets foff and poff
//  rather than by moving the contents.  Advancing a primary belt therefore
//  means incrementing poff and clearing the one slot that has just passed the
//  last crossing: that cell is reused as the fresh empty slot entering at the
//  upstream end, which is exactly the assumption that items reaching the end of
//  a primary belt leave the system at once.

struct Replication {
    const Layout &ly;
    long long warmup, steps;
    uint64_t seed;

    std::vector<std::vector<int>> feeder;
    std::vector<int> prim, bufF, bufB;
    std::vector<long long> foff;
    int poff = 0;
    std::vector<Rng> rng;
    std::vector<long long> admissions;   // per feeder, during the window only

    Replication(const Layout &L, long long w, long long t, uint64_t sd)
        : ly(L), warmup(w), steps(t), seed(sd) {}

    inline int fidx(int j, long long q) const {
        long long a = q - foff[j];
        if (a < 0) a += ly.L[j];
        return (int)a;
    }
    inline int pidx(int i, int p) const {
        int a = p - ly.pmin - poff;
        if (a < 0) a += ly.P;
        return i * ly.P + a;
    }

    //  One crossing of feeder j with primary belt i, at loop position q and
    //  belt position p, with buffer capacity c and buffer content *buf.
    //
    //  The discipline is:
    //    - if the primary slot is empty, the buffer discharges into it first;
    //      the feeder cell may then take the place freed in the buffer;
    //    - if the primary slot is empty and the buffer is empty, the feeder
    //      cell transfers directly;
    //    - if the primary slot is occupied, the feeder cell joins the buffer
    //      when there is room, and otherwise stays on the feeder for another
    //      revolution.
    //  With c = 0 this reduces to a direct transfer when the slot is free.
    inline void transfer(int j, int i, long long q, int p, int c, int *buf) {
        const int fa = fidx(j, q);
        const int pa = pidx(i, p);
        int &cell = feeder[j][fa];
        const bool has = (cell == i + 1);
        int b = buf ? *buf : 0;
        if (!prim[pa]) {
            if (b > 0) {
                prim[pa] = 1; b--;
                if (has && b < c) { b++; cell = 0; }
            } else if (has) {
                prim[pa] = 1; cell = 0;
            }
        } else if (has && b < c) {
            b++; cell = 0;
        }
        if (buf) *buf = b;
    }

    //  Returns the average throughput over the measurement window.
    double run() {
        const int n = ly.n, m = ly.m;
        feeder.resize(m);
        for (int j = 0; j < m; j++) feeder[j].assign((size_t)ly.L[j], 0);
        prim.assign((size_t)n * ly.P, 0);
        bufF.assign((size_t)m * n, 0);
        bufB.assign((size_t)m * n, 0);
        foff.assign(m, 0);
        admissions.assign(m, 0);
        poff = 0;

        rng.clear();
        rng.reserve(m);
        for (int j = 0; j < m; j++) {
            uint64_t s = seed ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(j + 1));
            rng.emplace_back(splitmix64(s));
        }

        const long long total = warmup + steps;
        for (long long t = 0; t < total; t++) {
            const bool counting = (t >= warmup);

            // 1. loading
            for (int j = 0; j < m; j++) {
                const int fa = fidx(j, 0);
                if (feeder[j][fa] == 0) {
                    feeder[j][fa] = (int)rng[j].below((uint32_t)n) + 1;
                    if (counting) admissions[j]++;
                }
            }
            // 2. transfers
            for (int j = 0; j < m; j++) {
                const size_t base = (size_t)j * n;
                for (int i = 0; i < n; i++) {
                    const int c = ly.cf[base + i];
                    transfer(j, i, ly.qf[base + i], ly.sf[j], c, c ? &bufF[base + i] : nullptr);
                }
                if (ly.dual) {
                    for (int i = 0; i < n; i++) {
                        const int c = ly.cb[base + i];
                        transfer(j, i, ly.qb[base + i], ly.sb[j], c, c ? &bufB[base + i] : nullptr);
                    }
                }
            }
            // 3. advancement
            if (++poff == ly.P) poff = 0;
            for (int i = 0; i < n; i++) prim[pidx(i, ly.pmin)] = 0;
            for (int j = 0; j < m; j++)
                if (++foff[j] == ly.L[j]) foff[j] = 0;
        }

        long long sum = 0;
        for (int j = 0; j < m; j++) sum += admissions[j];
        return (double)sum / (double)steps;
    }
};

// ------------------------------------------------------------ statistics ----

// Two-sided 97.5% Student t quantile.  Exact table up to 30 degrees of
// freedom, Cornish-Fisher expansion beyond it (error below 1e-4).
static double tquantile975(int df) {
    static const double tab[31] = {
        0.0,
        12.706, 4.3027, 3.1824, 2.7764, 2.5706, 2.4469, 2.3646, 2.3060,
        2.2622, 2.2281, 2.2010, 2.1788, 2.1604, 2.1448, 2.1314, 2.1199,
        2.1098, 2.1009, 2.0930, 2.0860, 2.0796, 2.0739, 2.0687, 2.0639,
        2.0595, 2.0555, 2.0518, 2.0484, 2.0452, 2.0423 };
    if (df <= 0) return 0.0;
    if (df <= 30) return tab[df];
    const double z = 1.959963985, z2 = z * z, d = df;
    return z + (z2 * z + z) / (4 * d)
             + (5 * z2 * z2 * z + 16 * z2 * z + 3 * z) / (96 * d * d);
}

// --------------------------------------------------------------- parsing ----

static std::vector<long long> parseList(const std::string &s) {
    std::vector<long long> v;
    std::string tok;
    std::istringstream in(s);
    while (std::getline(in, tok, ',')) {
        if (tok.empty()) continue;
        v.push_back(std::stoll(tok));
    }
    return v;
}

static void usage() {
    std::printf(R"(meshsorter - replication-based simulator of the MeshSorter

USAGE
  meshsorter -n BELTS -m FEEDERS [options]

GEOMETRY
  -n, --belts N          number of primary belts                      (required)
  -m, --feeders M        number of feeder loops                       (required)
      --dual             dual-drop: every feeder crosses every belt
                         twice, forward and backward                  (default)
      --single           single-drop: forward crossings only
      --dp D             feeder-loop slots consumed by each primary
                         belt; consecutive belts sit D/2 slots apart
                         on each run.  D must be even                 (default 4)
      --turn E           feeder-loop slots consumed by each end of
                         the loop                                     (default 4)
  -L, --loop LEN         base feeder-loop length          (default N*D + 2*E)
      --loops a,b,...    explicit lengths for the M feeders; overrides
                         -L and --stagger
      --stagger          use L, L, L+s, L+2s, ... for the M feeders
      --stagger-step s   the increment s                              (default 1)
      --extra WHERE      where slots beyond the base length are
                         inserted: turnaround or return               (default turnaround)
      --df F             spacing of consecutive feeders along a
                         primary belt                                 (default 4)
      --width W          distance along a primary belt between the
                         forward and the backward crossing of one
                         feeder                                       (default 2)

BUFFERS
  -b, --buffers SPEC     buffer capacity at the crossings.  SPEC is
                         either one value for every crossing, or M
                         values, one per feeder, or 2*M*N values
                         listed feeder by feeder, forward crossings
                         1..N then backward crossings 1..N            (default 0)

EXPERIMENT
  -T, --steps T          measured time steps per replication          (default 4000000)
  -R, --reps R           independent replications                     (default 10)
  -w, --warmup W         time steps discarded from each replication
                         (default max(20000, 10*(c+1)*L), the rule of
                         Supplement S1, where c is the largest buffer
                         capacity and L the longest feeder loop)
      --seed S           base seed                                    (default 20260901)
  -t, --threads K        worker threads; the result does not depend
                         on this                                      (default R)

OUTPUT
      --json             machine-readable output
      --per-feeder       also report the loader utilization of each feeder
  -h, --help             this text

EXAMPLES
  Table 1(b), the 4 by 4 cell:
      meshsorter -n 4 -m 4 --dual
  Table 2(c), the 4 by 9 cell:
      meshsorter -n 4 -m 9 --dual --stagger --extra turnaround
  A buffered system with three places at every crossing:
      meshsorter -n 4 -m 6 --dual -b 3
)");
}

int main(int argc, char **argv) {
    Config cfg;
    std::string bufSpec;
    try {
        for (int a = 1; a < argc; a++) {
            std::string k = argv[a];
            auto need = [&](void) -> std::string {
                if (a + 1 >= argc) throw std::runtime_error(k + " needs a value");
                return std::string(argv[++a]);
            };
            if      (k == "-h" || k == "--help")  { usage(); return 0; }
            else if (k == "-n" || k == "--belts")        cfg.n = std::stoi(need());
            else if (k == "-m" || k == "--feeders")      cfg.m = std::stoi(need());
            else if (k == "--dual")                      cfg.dual = true;
            else if (k == "--single")                    cfg.dual = false;
            else if (k == "--dp")                        cfg.dp = std::stoi(need());
            else if (k == "--turn")                      cfg.turn = std::stoi(need());
            else if (k == "-L" || k == "--loop")         cfg.Lbase = std::stoll(need());
            else if (k == "--loops")                     cfg.loops = parseList(need());
            else if (k == "--stagger")                   cfg.stagger = true;
            else if (k == "--stagger-step")              cfg.stagStep = std::stoll(need());
            else if (k == "--extra") {
                std::string w = need();
                if (w == "turnaround")   cfg.extraTurnaround = true;
                else if (w == "return")  cfg.extraTurnaround = false;
                else throw std::runtime_error("--extra takes turnaround or return");
            }
            else if (k == "--df")                        cfg.df = std::stoi(need());
            else if (k == "--width")                     cfg.width = std::stoi(need());
            else if (k == "-b" || k == "--buffers")      bufSpec = need();
            else if (k == "-T" || k == "--steps")        cfg.steps = std::stoll(need());
            else if (k == "-R" || k == "--reps")         cfg.reps = std::stoi(need());
            else if (k == "-w" || k == "--warmup")       cfg.warmup = std::stoll(need());
            else if (k == "--seed")                      cfg.seed = std::stoull(need());
            else if (k == "-t" || k == "--threads")      cfg.threads = std::stoi(need());
            else if (k == "--json")                      cfg.json = true;
            else if (k == "--per-feeder")                cfg.perFeeder = true;
            else throw std::runtime_error("unknown option " + k);
        }
        if (cfg.n < 1 || cfg.m < 1) { usage(); return 2; }
        if (cfg.reps < 2) throw std::runtime_error("--reps must be at least 2");
        if (cfg.steps < 1) throw std::runtime_error("--steps must be positive");
        if (cfg.width > cfg.df)
            throw std::runtime_error("--width may not exceed --df, or feeders would overlap");

        // buffer capacities
        const size_t K = (size_t)cfg.m * cfg.n;
        cfg.bufF.assign(K, 0); cfg.bufB.assign(K, 0);
        if (!bufSpec.empty()) {
            std::vector<long long> v = parseList(bufSpec);
            if (v.size() == 1) {
                std::fill(cfg.bufF.begin(), cfg.bufF.end(), (int)v[0]);
                std::fill(cfg.bufB.begin(), cfg.bufB.end(), (int)v[0]);
            } else if ((int)v.size() == cfg.m) {
                for (int j = 0; j < cfg.m; j++)
                    for (int i = 0; i < cfg.n; i++)
                        cfg.bufF[(size_t)j*cfg.n+i] = cfg.bufB[(size_t)j*cfg.n+i] = (int)v[j];
            } else if (v.size() == 2 * K) {
                for (int j = 0; j < cfg.m; j++)
                    for (int i = 0; i < cfg.n; i++) {
                        cfg.bufF[(size_t)j*cfg.n+i] = (int)v[(size_t)2*j*cfg.n + i];
                        cfg.bufB[(size_t)j*cfg.n+i] = (int)v[(size_t)2*j*cfg.n + cfg.n + i];
                    }
            } else {
                throw std::runtime_error("--buffers needs 1, m, or 2*m*n values");
            }
        }

        Layout ly = buildLayout(cfg);

        if (cfg.warmup < 0)
            cfg.warmup = std::max(20000LL, 10LL * (ly.cmax + 1) * ly.Lmax);
        if (cfg.threads <= 0) cfg.threads = cfg.reps;

        // ---- run the replications, one worker thread at a time ------------
        std::vector<double> Y((size_t)cfg.reps, 0.0);
        std::vector<std::vector<long long>> perFeeder((size_t)cfg.reps);
        int done = 0;
        while (done < cfg.reps) {
            const int batch = std::min(cfg.threads, cfg.reps - done);
            std::vector<std::thread> pool;
            pool.reserve(batch);
            for (int b = 0; b < batch; b++) {
                const int r = done + b;
                pool.emplace_back([&, r]() {
                    uint64_t s = cfg.seed + 0x1000193ULL * (uint64_t)(r + 1);
                    Replication rep(ly, cfg.warmup, cfg.steps, splitmix64(s));
                    Y[(size_t)r] = rep.run();
                    perFeeder[(size_t)r] = rep.admissions;
                });
            }
            for (auto &th : pool) th.join();
            done += batch;
        }

        // ---- statistics ----------------------------------------------------
        const int R = cfg.reps;
        double mean = 0.0;
        for (double y : Y) mean += y;
        mean /= R;
        double ss = 0.0;
        for (double y : Y) ss += (y - mean) * (y - mean);
        const double sd = std::sqrt(ss / (R - 1));
        const double se = sd / std::sqrt((double)R);
        const double hw = tquantile975(R - 1) * se;

        std::vector<double> util((size_t)cfg.m, 0.0);
        for (int r = 0; r < R; r++)
            for (int j = 0; j < cfg.m; j++)
                util[(size_t)j] += (double)perFeeder[(size_t)r][(size_t)j] / (double)cfg.steps / R;

        // ---- output --------------------------------------------------------
        if (cfg.json) {
            std::printf("{\"n\":%d,\"m\":%d,\"dual\":%s,\"loops\":[", cfg.n, cfg.m,
                        cfg.dual ? "true" : "false");
            for (int j = 0; j < cfg.m; j++) std::printf("%s%lld", j ? "," : "", ly.L[j]);
            std::printf("],\"buffer_max\":%d,\"reps\":%d,\"steps\":%lld,\"warmup\":%lld,"
                        "\"seed\":%llu,\"throughput\":%.9f,\"halfwidth\":%.9f,\"sd\":%.9f,"
                        "\"replications\":[",
                        ly.cmax, R, cfg.steps, cfg.warmup,
                        (unsigned long long)cfg.seed, mean, hw, sd);
            for (int r = 0; r < R; r++) std::printf("%s%.9f", r ? "," : "", Y[(size_t)r]);
            std::printf("]");
            if (cfg.perFeeder) {
                std::printf(",\"loader_utilization\":[");
                for (int j = 0; j < cfg.m; j++) std::printf("%s%.6f", j ? "," : "", util[(size_t)j]);
                std::printf("]");
            }
            std::printf("}\n");
        } else {
            std::printf("MeshSorter, %s, n = %d primary belts, m = %d feeder loops\n",
                        cfg.dual ? "dual-drop" : "single-drop", cfg.n, cfg.m);
            std::printf("  feeder loop lengths   ");
            for (int j = 0; j < cfg.m; j++) std::printf(" %lld", ly.L[j]);
            std::printf("\n  crossings on a loop    forward");
            for (int i = 0; i < cfg.n; i++) std::printf(" %lld", ly.qf[(size_t)i]);
            if (cfg.dual) {
                std::printf("\n                         backward, feeder 1");
                for (int i = 0; i < cfg.n; i++) std::printf(" %lld", ly.qb[(size_t)i]);
                const size_t last = (size_t)(cfg.m - 1) * cfg.n;
                if (cfg.m > 1 && ly.qb[last] != ly.qb[0]) {
                    std::printf("\n                         backward, feeder %d", cfg.m);
                    for (int i = 0; i < cfg.n; i++) std::printf(" %lld", ly.qb[last + i]);
                }
            }
            std::printf("\n  crossings on a belt    s^f");
            for (int j = 0; j < cfg.m; j++) std::printf(" %d", ly.sf[j]);
            if (cfg.dual) {
                std::printf("   s^b");
                for (int j = 0; j < cfg.m; j++) std::printf(" %d", ly.sb[j]);
            }
            std::printf("\n  buffers                %s (largest capacity %d)\n",
                        ly.cmax ? "yes" : "none", ly.cmax);
            std::printf("  design                 %d replications of %lld steps, "
                        "warm-up %lld each, seed %llu\n",
                        R, cfg.steps, cfg.warmup, (unsigned long long)cfg.seed);
            std::printf("\n  throughput   %.5f   95%% half-width %.5f   [%.5f, %.5f]\n",
                        mean, hw, mean - hw, mean + hw);
            std::printf("  replications ");
            for (int r = 0; r < R; r++) std::printf(" %.5f", Y[(size_t)r]);
            std::printf("\n");
            if (cfg.perFeeder) {
                std::printf("  loader utilization");
                for (int j = 0; j < cfg.m; j++) std::printf(" %.4f", util[(size_t)j]);
                std::printf("\n");
            }
        }
    } catch (const std::exception &e) {
        std::fprintf(stderr, "meshsorter: %s\n", e.what());
        return 1;
    }
    return 0;
}
