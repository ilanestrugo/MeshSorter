#!/usr/bin/env python3
"""Figure S1: the efficient frontier of throughput against load balancing.

At a per-primary-belt budget of ten, every allocation the design rules leave
open is evaluated and plotted as one point, throughput against the minimum
loader utilization.  The non-dominated subset is the frontier drawn over them.

The space is larger than the structured class that Tables S5 and S6 search.  It
fixes nothing on feeder 1, which is never blocked, and under the dual-drop
mechanism it places capacity only at the backward drop points; but it does not
require the capacities to be nondecreasing, so it holds all 66 compositions of
the budget rather than the 14 partitions.  That is the point of the figure: it
shows what the designer gives up by moving along the frontier, over the whole
space the rules admit.

    python3 frontier.py           both mechanisms, about ten minutes
    python3 frontier.py --quick   a cheap pass, to check the wiring

Every allocation is evaluated under the protocol of Section S1, the same
geometry, run length, warm-up rule and thirty replications used everywhere else,
so the figure and Tables S5 and S6 are directly comparable.

Writes the four data files the manuscript's figure reads, which belong beside
the manuscript in its data/ directory:

    m4_n4_B10_OneWay_efficiency_front_data.dat   every allocation, single-drop
    m4_n4_B10_OneWay_efficiency_front.dat        its non-dominated subset
    m4_n4_B10_TwoWay_efficiency_front_data.dat   every allocation, dual-drop
    m4_n4_B10_TwoWay_efficiency_front.dat        its non-dominated subset
"""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN  = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
OUT  = os.path.join(HERE, "results")
CELLS = os.path.join(OUT, "frontier")

# ---------------------------------------------------------------- parameters
BELTS    = 4
FEEDERS  = 4
BUDGET   = 10
import geometry as _geo
SPACING  = _geo.SPACING              # slots between consecutive primary belts along a loop
TURN     = _geo.TURN
DF     = _geo.DF
WIDTH     = _geo.WIDTH
STEPS    = 1_330_000
REPS     = 30
SEED     = 20260901
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, REPS, BUDGET = 100_000, 6, 3

if not os.path.exists(BIN):
    sys.exit(f"{BIN} not found.  Build it first:\n"
             f"    c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp")
os.makedirs(CELLS, exist_ok=True)


def design():
    return dict(belts=BELTS, feeders=FEEDERS, budget=BUDGET, spacing=SPACING,
                turn=TURN, steps=STEPS, reps=REPS, seed=SEED)


def allocations(B, m):
    """Nothing on feeder 1; the budget spread over the remaining m-1 positions
    in any order.  These are the compositions of B into m-1 parts."""
    res = []
    def rec(k, remaining, cur):
        if k == 1:
            res.append([0] + cur + [remaining]); return
        for v in range(remaining + 1):
            rec(k - 1, remaining - v, cur + [v])
    rec(m - 1, B, [])
    return res


def buffer_spec(c, dual):
    out = []
    for cj in c:
        out += ([0] * BELTS + [cj] * BELTS) if dual else ([cj] * BELTS + [0] * BELTS)
    return ",".join(map(str, out))


def evaluate(c, dual):
    cmd = [BIN, "-n", str(BELTS), "-m", str(FEEDERS),
           "--dual" if dual else "--single",
           "--spacing", str(SPACING), "--turn", str(TURN),
           "--df", str(DF), "--width", str(WIDTH),
           "-b", buffer_spec(c, dual),
           "-T", str(STEPS), "-R", str(REPS), "--seed", str(SEED),
           "--per-feeder", "--json"]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"meshsorter failed on {c}: {p.stderr.strip()}")
    d = json.loads(p.stdout)
    return dict(alloc=",".join(map(str, c)), th=d["throughput"],
                load=min(d["loader_utilization"]))


def pareto(rows):
    """Non-dominated on both objectives, largest throughput first."""
    out = []
    for r in sorted(rows, key=lambda r: (-r["th"], -r["load"])):
        if all(r["load"] > s["load"] for s in out):
            out.append(r)
    return out


def run(mech, dual):
    path = os.path.join(CELLS, f"{mech}_B{BUDGET}.json")
    if os.path.exists(path):
        old = json.load(open(path))
        if old.get("design") == design():
            return old["rows"]
        print(f"  {mech}: cached under a different design, rerunning", file=sys.stderr)
    t0 = time.time()
    rows = []
    for i, c in enumerate(allocations(BUDGET, FEEDERS), 1):
        rows.append(evaluate(c, dual))
        if i % 10 == 0:
            print(f"  {mech}: {i} evaluated", file=sys.stderr, flush=True)
    json.dump(dict(design=design(), rows=rows), open(path, "w"), indent=1)
    print(f"  {mech}: {len(rows)} allocations in {time.time() - t0:.0f}s",
          file=sys.stderr, flush=True)
    return rows


t0 = time.time()
for mech, dual, stem in (("single", False, "OneWay"), ("dual", True, "TwoWay")):
    rows = run(mech, dual)
    front = pareto(rows)
    for name, data in ((f"m4_n4_B{BUDGET}_{stem}_efficiency_front_data.dat", rows),
                       (f"m4_n4_B{BUDGET}_{stem}_efficiency_front.dat",
                        sorted(front, key=lambda r: r["th"]))):
        with open(os.path.join(OUT, name), "w") as f:
            f.write("Throughput\tLoad\n")
            for r in data:
                f.write(f"{r['th']:.6f}\t{r['load']:.6f}\n")
    best_th   = max(rows, key=lambda r: r["th"])
    best_load = max(rows, key=lambda r: r["load"])
    print(f"  {mech}: throughput optimum {best_th['alloc']} "
          f"(TH {best_th['th']:.3f}, load {best_load and best_th['load']:.3f}); "
          f"best load {best_load['alloc']} (load {best_load['load']:.3f}, "
          f"TH {best_load['th']:.3f}); frontier holds {len(front)}", file=sys.stderr)

print(f"\nwrote four .dat files to results/ in {time.time() - t0:.0f}s.  Copy them "
      f"into the manuscript's data/ directory.", file=sys.stderr)
