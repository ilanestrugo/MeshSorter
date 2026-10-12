"""How far is the independent-draws approximation from the unbuffered single-drop system?

Section 6.5 of the manuscript.  The approximation A(n,m) = n(1-(1-1/n)^m) treats
the destinations of the m feeders as independent.  Two statements about it need
numbers:

  1. Against the exact throughput Q(n,m) of the common-loop system
     (exact_single_drop.py, rational arithmetic), A never underestimates:
     checked over every pair with 2 <= n <= 16 and 1 <= m <= 16, 240 pairs.
     The script reports the smallest and largest relative overestimate, the pairs
     where A equals Q, and the pairs where it falls below Q, if any.  The 49
     pairs of Table 2(a), n, m = 3..9, are reported separately; their range is
     the one quoted from Table 4.

  2. Against the simulated throughput of the staggered configuration, over the
     same 49 pairs.  Nothing is simulated here: the inputs are the cells that
     stagger_grid.py writes to results/table2.csv (panel table2a_single_staggered,
     Table S5(a)), so this script reads what that table prints.  The script
     reports the largest relative departure of the simulated value from A.

    python3 approximation_gap.py            about ten minutes for part 1
    python3 approximation_gap.py --table    part 2 only, instant

Output
    results/approximation_gap.json   both parts, cell by cell
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "results", "table2.csv")
PANEL = "table2a_single_staggered"
GRID = range(3, 10)


def approx(n, m):
    return n * (1.0 - (1.0 - 1.0 / n) ** m)


def pct(x):
    return 100.0 * x


out = {}

# ------------------------------------------------- 1. approximation against exact
if "--table" not in sys.argv:
    from exact_single_drop import exact
    gap = {}
    for n in range(2, 17):
        for m in range(1, 17):
            q = exact(n, m)                       # a Fraction
            gap[(n, m)] = pct(approx(n, m) / float(q) - 1.0)
        print("n = %2d done" % n, file=sys.stderr, flush=True)
    below = [k for k, v in gap.items() if v < -1e-9]
    equal = sorted(k for k, v in gap.items() if abs(v) <= 1e-9)
    table = {k: v for k, v in gap.items() if k[0] in GRID and k[1] in GRID}
    print("exact against approximation, %d pairs with 2<=n<=16, 1<=m<=16" % len(gap))
    print("  A below Q in %d pairs%s" % (len(below), (": %s" % below[:5]) if below else ""))
    print("  A equals Q in %d pairs, all with m <= 2: %s"
          % (len(equal), all(m <= 2 for _, m in equal)))
    print("  largest overestimate %.2f%% at (n, m) = %s"
          % (max(gap.values()), max(gap, key=gap.get)))
    print("the %d pairs of Table 2(a), n, m = 3..9" % len(table))
    print("  overestimate %.2f%% to %.2f%%, above Q in %d of %d"
          % (min(table.values()), max(table.values()),
             sum(1 for v in table.values() if v > 1e-9), len(table)))
    out["exact"] = {"pairs": len(gap), "below": len(below),
                    "max_overestimate_percent": max(gap.values()),
                    "table2a_min_percent": min(table.values()),
                    "table2a_max_percent": max(table.values()),
                    "percent": {"%d,%d" % k: v for k, v in sorted(gap.items())}}

# ---------------------------------------- 2. approximation against the staggered runs
rows = [r for r in csv.DictReader(open(SRC)) if r["panel"] == PANEL]
sim = {(int(r["n"]), int(r["m"])): float(r["throughput"]) for r in rows}
missing = [(n, m) for n in GRID for m in GRID if (n, m) not in sim]
if missing:
    raise SystemExit("%s: %s has no cells %s" % (SRC, PANEL, missing[:5]))
dev = {k: pct(sim[k] / approx(*k) - 1.0) for k in sim if k[0] in GRID and k[1] in GRID}
worst = max(dev, key=lambda k: abs(dev[k]))
print("staggered single-drop simulation against the approximation, %d cells (Table S5(a))" % len(dev))
print("  largest departure %.2f%% at (n, m) = %s; simulated above A in %d cells"
      % (abs(dev[worst]), worst, sum(1 for v in dev.values() if v > 0)))
out["staggered"] = {"cells": len(dev), "max_abs_percent": abs(dev[worst]),
                    "percent": {"%d,%d" % k: v for k, v in sorted(dev.items())}}

if "--table" not in sys.argv:
    json.dump(out, open(os.path.join(HERE, "results", "approximation_gap.json"), "w"), indent=1)
    print("wrote results/approximation_gap.json")
