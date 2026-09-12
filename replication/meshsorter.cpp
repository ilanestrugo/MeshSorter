// -----------------------------------------------------------------------------
//  meshsorter: command-line front end for the MeshSorter replication design.
//  The model itself is in meshsorter_core.hpp; see README.md.
//
//  Build:  c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp
// -----------------------------------------------------------------------------
#include <fstream>
#include "meshsorter_core.hpp"

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
      --spacing S        the same thing stated as the distance between
                         consecutive primary belts along the loop;
                         --spacing S is exactly --dp 2S               (default 2)
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
  -T, --steps T          measured time steps per replication          (default 1330000)
  -R, --reps R           independent replications                     (default 30)
  -w, --warmup W         time steps discarded from each replication
                         (default max(20000, 10*(c+1)*L), the rule of
                         Supplement S1, where c is the largest buffer
                         capacity and L the longest feeder loop)
      --seed S           base seed                                    (default 20260901)
      --sequence FILE    read destination labels from FILE instead of
                         drawing them uniformly: one integer in 1..N per
                         line, the order-derived sequence of Supplement
                         S4.  The feeders consume it in index order and
                         it repeats cyclically, and each replication
                         starts at its own position in the cycle, so the
                         replications differ in phase rather than in the
                         labels they see
  -t, --threads K        worker threads; the result does not depend on
                         this, only the wall clock does
                         (default: two fewer than the hardware
                          threads this machine reports, capped at R)

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
  One row of Table 6, on the order-derived sequence:
      meshsorter -n 4 -m 4 --dual --spacing 4 -b 0,0,1,2 \
                 --sequence results/order_derived_sequence.txt
)");
}

int main(int argc, char **argv) {
    Config cfg;
    std::string seqFile;
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
            else if (k == "--spacing")                   cfg.dp = 2 * std::stoi(need());
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
            else if (k == "--sequence")                  seqFile = need();
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

        //  Order-derived destination sequence, if one was named.  Labels are
        //  1..n, one per line; blank lines and lines starting with # are
        //  ignored, so the file can carry a provenance header.
        std::vector<int> seq;
        if (!seqFile.empty()) {
            std::ifstream in(seqFile);
            if (!in) throw std::runtime_error("cannot open " + seqFile);
            std::string line;
            while (std::getline(in, line)) {
                size_t a = line.find_first_not_of(" \t\r");
                if (a == std::string::npos || line[a] == '#') continue;
                int v = std::stoi(line.substr(a));
                if (v < 1 || v > cfg.n)
                    throw std::runtime_error("destination label out of range in "
                                             + seqFile + ": " + std::to_string(v));
                seq.push_back(v);
            }
            if (seq.empty()) throw std::runtime_error("no labels in " + seqFile);
        }

        if (cfg.warmup < 0)
            cfg.warmup = std::max(20000LL, 10LL * (ly.cmax + 1) * ly.Lmax);
        if (cfg.threads <= 0) {
            // Leave two hardware threads for the operating system and for
            // whatever else the machine is doing, and never start more workers
            // than there are replications to run.
            const unsigned hw = std::thread::hardware_concurrency();
            const int avail = (hw > 3u) ? (int)hw - 2 : 1;
            cfg.threads = std::min(cfg.reps, avail);
        }

        // ---- run the replications ------------------------------------------
        // The workers pull replications from a shared counter, so a worker that
        // finishes early picks up the next one instead of idling until the rest
        // of its wave is done.  Replication r always uses the same stream, so
        // the result does not depend on how many workers there are.
        std::vector<double> Y((size_t)cfg.reps, 0.0);
        std::vector<std::vector<long long>> perFeeder((size_t)cfg.reps);
        std::atomic<int> nextRep{0};
        {
            std::vector<std::thread> pool;
            pool.reserve((size_t)cfg.threads);
            for (int w = 0; w < cfg.threads; w++) {
                pool.emplace_back([&]() {
                    for (;;) {
                        const int r = nextRep.fetch_add(1);
                        if (r >= cfg.reps) return;
                        uint64_t s = cfg.seed + 0x1000193ULL * (uint64_t)(r + 1);
                        Replication rep(ly, cfg.warmup, cfg.steps, splitmix64(s));
                        if (!seq.empty()) rep.seq = &seq;
                        Y[(size_t)r] = rep.run();
                        perFeeder[(size_t)r] = rep.admissions;
                    }
                });
            }
            for (auto &th : pool) th.join();
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
                        "\"seed\":%llu,\"threads\":%d,\"throughput\":%.9f,"
                        "\"halfwidth\":%.9f,\"sd\":%.9f,\"replications\":[",
                        ly.cmax, R, cfg.steps, cfg.warmup,
                        (unsigned long long)cfg.seed, cfg.threads, mean, hw, sd);
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
                        "warm-up %lld each, seed %llu\n"
                        "                         %d worker thread%s\n",
                        R, cfg.steps, cfg.warmup, (unsigned long long)cfg.seed,
                        cfg.threads, cfg.threads == 1 ? "" : "s");
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
