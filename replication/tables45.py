#!/usr/bin/env python3
"""Tables S4 and S5: throughput against load balancing, over the structured class.

At each per-primary-belt budget the two objectives are optimized separately and
the cost of preferring one to the other is reported:

    c   the allocation with the largest minimum loader utilization
    c'  the allocation with the largest throughput

The search is confined to the structured class of allocations the design rules
admit, which Section 5.1 of the main manuscript certifies as containing a
near-optimal design.  That is a deliberate restriction: the load-balancing
optimum is not what the design rules were derived for, so this reports the best
balance available to a designer who follows them, not the best in the whole
allocation space.

    python3 tables45.py            both mechanisms, about ten minutes
    python3 tables45.py --quick    a cheap pass, to check the wiring
    python3 tables45.py --only dual

Every allocation is evaluated under the protocol of Section S1: the geometry,
run length, warm-up rule and thirty replications used by the certification and
by Section 7, so these tables and Table 4 are directly comparable.  The earlier
version of these tables was run at the shorter loop length and is not.

Results are cached per cell, so an interrupted run resumes and a rerun is free.
"""
import csv, itertools, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN  = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
OUT  = os.path.join(HERE, "results")
CELLS = os.path.join(OUT, "loadbalance")

# ---------------------------------------------------------------- parameters
BELTS    = 4
FEEDERS  = 4
BUDGETS  = list(range(0, 11))
SPACING  = 4              # slots between consecutive primary belts along a loop
TURN     = 4              # slots in each end curve; loop length is 2*S*n + 2*E
STEPS    = 1_330_000
REPS     = 30
SEED     = 20260901
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, REPS, BUDGETS = 100_000, 6, [0, 1, 2, 3]
ONLY = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None

if not os.path.exists(BIN):
    sys.exit(f"{BIN} not found.  Build it first:\n"
             f"    c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp")
os.makedirs(CELLS, exist_ok=True)


def design():
    return dict(belts=BELTS, feeders=FEEDERS, spacing=SPACING, turn=TURN,
                steps=STEPS, reps=REPS, seed=SEED)


def structured(B, m):
    """The structured class: nothing on feeder 1, capacities nondecreasing.

    Both design rules reduce to the same object, the partitions of B into at
    most m-1 parts, laid out in nondecreasing order along the feeders.  Under
    the dual-drop mechanism these are the capacities at the backward drop
    points, the forward ones being zero.
    """
    res = []
    def rec(k, remaining, lo, cur):
        if k == 1:
            if remaining >= lo:
                res.append([0] + cur + [remaining])
            return
        for v in range(lo, remaining // k + 1):
            rec(k - 1, remaining - v, v, cur + [v])
    if B == 0:
        return [[0] * m]
    rec(m - 1, B, 0, [])
    return res


def buffer_spec(c, dual):
    """The 2*m*n form: capacity at the backward drop points under dual-drop,
    at the forward ones under single-drop, replicated across the belts."""
    out = []
    for cj in c:
        out += ([0] * BELTS + [cj] * BELTS) if dual else ([cj] * BELTS + [0] * BELTS)
    return ",".join(map(str, out))


def evaluate(c, dual):
    cmd = [BIN, "-n", str(BELTS), "-m", str(FEEDERS),
           "--dual" if dual else "--single",
           "--spacing", str(SPACING), "--turn", str(TURN),
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
                hw=d["halfwidth"], loads=d["loader_utilization"],
                minload=min(d["loader_utilization"]))


def cell(dual, B):
    tag = f"{'dual' if dual else 'single'}_n{BELTS}_m{FEEDERS}_B{B}"
    path = os.path.join(CELLS, tag + ".json")
    if os.path.exists(path):
        old = json.load(open(path))
        if old.get("design") == design():
            return old
        print(f"  {tag}: cached under a different design, rerunning", file=sys.stderr)
    t0 = time.time()
    rows = [evaluate(c, dual) for c in structured(B, FEEDERS)]
    #  Ties: prefer the allocation that is also better on the other objective,
    #  so a cell where one allocation is best at both reports a gap of zero
    #  rather than two different names for the same design.
    balanced = max(rows, key=lambda r: (r["minload"], r["th"]))
    fastest  = max(rows, key=lambda r: (r["th"], r["minload"]))
    r = dict(n=BELTS, m=FEEDERS, B=B, dual=dual, allocations=len(rows),
             balanced=balanced, fastest=fastest, rows=rows,
             seconds=round(time.time() - t0, 1), design=design())
    json.dump(r, open(path, "w"), indent=1)
    print(f"  {tag}: {len(rows):2d} allocations  balance {balanced['alloc']} "
          f"(L {balanced['minload']:.3f})  throughput {fastest['alloc']} "
          f"(TH {fastest['th']:.3f})  ({r['seconds']:.0f}s)", file=sys.stderr, flush=True)
    return r


def body(cells):
    """One LaTeX row per budget, in the nine-column form of Tables S4 and S5."""
    out = []
    for r in cells:
        b, f = r["balanced"], r["fastest"]
        same = b["alloc"] == f["alloc"]
        thg = "0" if same else f"{f['th'] - b['th']:.3f}"
        lg  = "0" if same else f"{b['minload'] - f['minload']:.3f}"
        out.append(f"{r['B']} & {b['alloc']} & {b['th']:.3f} & {b['minload']:.3f} & "
                   f"{f['alloc']} & {f['th']:.3f} & {f['minload']:.3f} & {thg} & {lg} \\\\")
    return "\n".join(out) + "\n"


t0 = time.time()
made = {}
for mech, dual in (("single", False), ("dual", True)):
    if ONLY and ONLY != mech:
        continue
    print(f"{mech}-drop:", file=sys.stderr)
    made[mech] = [cell(dual, B) for B in BUDGETS]
    open(os.path.join(OUT, f"table45_{mech}.tex"), "w").write(body(made[mech]))

with open(os.path.join(OUT, "table45.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["mechanism", "n", "m", "B", "allocations", "role",
                "alloc", "throughput", "halfwidth", "min_load", "loader_utilizations"])
    for mech, cells in made.items():
        for r in cells:
            for role in ("balanced", "fastest"):
                x = r[role]
                w.writerow([mech, r["n"], r["m"], r["B"], r["allocations"], role,
                            x["alloc"], f"{x['th']:.6f}", f"{x['hw']:.6f}",
                            f"{x['minload']:.6f}",
                            "|".join(f"{v:.6f}" for v in x["loads"])])

print(f"\nwrote {' '.join('results/table45_%s.tex' % m for m in made)} and "
      f"results/table45.csv in {time.time() - t0:.0f}s", file=sys.stderr)
