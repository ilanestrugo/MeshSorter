// transient_pilot.cpp
//
// The transient pilot of Section S2 of the Supplemental Online Material.
//
// Every replication of the simulator starts from an empty system, which is
// optimistic: in the first round every feeder meets an empty slot on every
// primary belt, so the throughput starts at its largest possible value and then
// settles.  This program measures that settling.  It runs R independent
// replications from the empty state, deletes nothing, records the admissions of
// every consecutive block of B steps, and prints the mean and the standard
// deviation over replications of the throughput of each block, together with
// the number of items held in buffers at the end.  transient_pilot.py turns
// that into the figures the supplement quotes.
//
//   transient_pilot -n 4 -m 4 --dual -L 100 -b 3 -R 300 -T 200000 --block 20
//
// The system is built exactly as meshsorter builds it, from the same Config and
// buildLayout, and the replication seeds are derived in the same way, so the
// two programs share a model and a stream family.  -L is the base loop length
// on the crossings of the built-in geometry, and -b the capacity at every drop
// point, forward and backward.

#include "meshsorter_core.hpp"
#include <iostream>

int main(int argc, char **argv) {
    Config cfg;
    cfg.dual = true;
    int reps = 300, cap = 0;
    long long steps = 200000, block = 20;
    try {
        for (int a = 1; a < argc; a++) {
            std::string k = argv[a];
            auto need = [&]() -> std::string {
                if (a + 1 >= argc) throw std::runtime_error(k + " needs a value");
                return argv[++a];
            };
            if (k == "-n") cfg.n = std::stoi(need());
            else if (k == "-m") cfg.m = std::stoi(need());
            else if (k == "-L") cfg.Lbase = std::stoll(need());
            else if (k == "-b") cap = std::stoi(need());
            else if (k == "-R") reps = std::stoi(need());
            else if (k == "-T") steps = std::stoll(need());
            else if (k == "--block") block = std::stoll(need());
            else if (k == "--seed") cfg.seed = std::stoull(need());
            else if (k == "-t") cfg.threads = std::stoi(need());
            else if (k == "--dual") cfg.dual = true;
            else if (k == "--single") cfg.dual = false;
            else throw std::runtime_error("unknown option " + k);
        }
        if (cfg.n < 1 || cfg.m < 1) throw std::runtime_error("-n and -m are required");
        cfg.bufF.assign((size_t)cfg.m * cfg.n, cap);
        cfg.bufB.assign((size_t)cfg.m * cfg.n, cap);
        const Layout ly = buildLayout(cfg);
        if (steps % block != 0) throw std::runtime_error("-T must be a multiple of --block");
        if (cfg.threads <= 0) {
            const unsigned hw = std::thread::hardware_concurrency();
            cfg.threads = std::min(reps, (hw > 3u) ? (int)hw - 2 : 1);
        }

        std::vector<std::vector<int>> trace((size_t)reps);
        std::vector<long long> occ((size_t)reps, 0);
        std::atomic<int> next{0};
        std::vector<std::thread> pool;
        for (int w = 0; w < cfg.threads; w++)
            pool.emplace_back([&]() {
                for (;;) {
                    const int r = next.fetch_add(1);
                    if (r >= reps) return;
                    uint64_t s = cfg.seed + 0x1000193ULL * (uint64_t)(r + 1);
                    Replication rep(ly, 0, steps, splitmix64(s));
                    rep.trace = &trace[(size_t)r];
                    rep.traceBlock = block;
                    rep.run();
                    occ[(size_t)r] = rep.bufferContentAtEnd;
                }
            });
        for (auto &t : pool) t.join();

        const size_t B = (size_t)(steps / block);
        std::vector<double> mean(B, 0.0), sd(B, 0.0);
        for (size_t b = 0; b < B; b++) {
            double s = 0.0;
            for (int r = 0; r < reps; r++) s += (double)trace[(size_t)r][b] / (double)block;
            mean[b] = s / reps;
            double v = 0.0;
            for (int r = 0; r < reps; r++) {
                const double d = (double)trace[(size_t)r][b] / (double)block - mean[b];
                v += d * d;
            }
            sd[b] = (reps > 1) ? std::sqrt(v / (reps - 1)) : 0.0;
        }
        double occMean = 0.0;
        for (int r = 0; r < reps; r++) occMean += (double)occ[(size_t)r];
        occMean /= reps;

        std::printf("{\"n\":%d,\"m\":%d,\"dual\":%s,\"loop\":%lld,\"capacity\":%d,"
                    "\"reps\":%d,\"steps\":%lld,\"block\":%lld,\"seed\":%llu,"
                    "\"occupancy_end\":%.4f,\"mean\":[",
                    cfg.n, cfg.m, cfg.dual ? "true" : "false", ly.Lmax, cap, reps, steps,
                    block, (unsigned long long)cfg.seed, occMean);
        for (size_t b = 0; b < B; b++) std::printf("%s%.7f", b ? "," : "", mean[b]);
        std::printf("],\"sd\":[");
        for (size_t b = 0; b < B; b++) std::printf("%s%.7f", b ? "," : "", sd[b]);
        std::printf("]}\n");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "transient_pilot: %s\n", e.what());
        return 1;
    }
    return 0;
}
