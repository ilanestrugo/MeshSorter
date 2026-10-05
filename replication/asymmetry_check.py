#!/usr/bin/env python3
"""What does the symmetry restriction cost?  (Section 6.1 of the manuscript)

The buffer-allocation experiments restrict attention to symmetric allocations:
the same capacity at the same feeder on every one of the n primary belts.  This
script perturbs the symmetric optimum of the 4 x 4 dual-drop system at a budget
of ten units per belt (forty system-wide), allowing the capacity of each of the
individual drop points to differ across belts, and measures what that costs.

Geometry, run length and warm-up are those of every other experiment: loops of
20n + 24 = 104 slots, 1,330,000 steps per replication, the warm-up rule of
Supplement S2.  Every allocation receives 100 independent replications on its
own seed, so every reported difference carries a confidence interval.  (An
earlier version of this check was a single run of 2 x 10^7 steps on loops of 24
slots, with no uncertainty estimate.)

Backward drop points only, as in the structured class for the dual-drop system.
Capacities are written per feeder, one value per primary belt.

    python3 asymmetry_check.py
    python3 asymmetry_check.py --quick     6 short replications per allocation;
                                           writes results/asymmetry_quick.*
"""
import csv, json, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
OUT = os.path.join(HERE, "results")
import geometry as _geo

BELTS, FEEDERS = 4, 4
STEPS, REPS, SEED0 = 1_330_000, 100, 20261004
SUFFIX = ""
if "--quick" in sys.argv:
    STEPS, REPS, SUFFIX = 100_000, 6, "_quick"

# the symmetric optimum, (0,1,2,7) on every belt: feeder j, then capacity per belt
BASE = [[0] * 4, [1] * 4, [2] * 4, [7] * 4]


def with_(base, edits):
    out = [row[:] for row in base]
    for (j, i), v in edits.items():
        out[j][i] = v
    return out


def spec(alloc):
    flat = []
    for row in alloc:
        flat += [0] * BELTS + list(row)
    return ",".join(map(str, flat))


def simulate(alloc, seed):
    cmd = [BIN, "-n", str(BELTS), "-m", str(FEEDERS), "--dual",
           "--spacing", str(_geo.SPACING), "--turn", str(_geo.TURN),
           "--df", str(_geo.DF), "--width", str(_geo.WIDTH),
           "-b", spec(alloc), "-T", str(STEPS), "-R", str(REPS),
           "--seed", str(seed), "--json"]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("meshsorter failed: " + p.stderr.strip())
    return json.loads(p.stdout)


def total(alloc):
    return sum(sum(r) for r in alloc)


def write_tex(rows):
    """The body of the supplement's table: one row per allocation."""
    with open(os.path.join(OUT, f"asymmetry{SUFFIX}.tex"), "w") as f:
        for r in rows:
            d = ("--" if r["name"] == "base" else
                 f"{float(r['diff_from_base']):+.4f} $\\pm$ {float(r['diff_halfwidth']):.4f}")
            f.write(f"{r['description']} & {r['budget']} & "
                    f"{float(r['throughput']):.4f} $\\pm$ {float(r['halfwidth']):.4f} & {d} \\\\\n")
        f.write("\\hline\n")


if "--tex-only" in sys.argv:           # rebuild the LaTeX body from the saved CSV
    write_tex(list(csv.DictReader(open(os.path.join(OUT, f"asymmetry{SUFFIX}.csv")))))
    sys.exit(0)

cases = [("base", "symmetric (0,1,2,7) on every belt", BASE)]

# 1. the twelve single-unit moves between belts, within the most buffered feeder
for a in range(BELTS):
    for b in range(BELTS):
        if a != b:
            cases.append((f"move_{a + 1}to{b + 1}", f"feeder 4: one unit from belt {a + 1} to belt {b + 1}",
                          with_(BASE, {(3, a): 6, (3, b): 8})))

# 2. large tilts of the whole allocation toward one end, budget held at 40
cases.append(("tilt_to_belt4", "tilt toward belt 4 across feeders 2 to 4",
              [[0] * 4, [0, 1, 1, 2], [1, 2, 2, 3], [4, 6, 8, 10]]))
cases.append(("tilt_to_belt1", "tilt toward belt 1 across feeders 2 to 4",
              [[0] * 4, [2, 1, 1, 0], [3, 2, 2, 1], [10, 8, 6, 4]]))

# 3. a residual unit, for budgets that are not a multiple of n: where does it belong?
cases.append(("plus1_f4", "one extra unit at feeder 4, belt 4", with_(BASE, {(3, 3): 8})))
cases.append(("plus1_f3", "one extra unit at feeder 3, belt 4", with_(BASE, {(2, 3): 3})))
cases.append(("plus1_f2", "one extra unit at feeder 2, belt 4", with_(BASE, {(1, 3): 2})))

t0 = time.time()
results = {}
for k, (name, desc, alloc) in enumerate(cases):
    r = simulate(alloc, SEED0 + 1000 * k)
    results[name] = dict(desc=desc, alloc=alloc, total=total(alloc), seed=SEED0 + 1000 * k,
                         mean=r["throughput"], hw=r["halfwidth"], sd=r["sd"], reps=REPS,
                         y=r["replications"])
    print(f"  {name:14s} total {total(alloc):2d}  {r['throughput']:.6f} +- {r['halfwidth']:.6f}",
          file=sys.stderr, flush=True)

# differences from the symmetric base, on independent streams: Welch interval
from statistics import mean, variance

b = results["base"]


def welch(x, y):
    nx, ny = len(x), len(y)
    vx, vy = variance(x) / nx, variance(y) / ny
    se = math.sqrt(vx + vy)
    df = (vx + vy) ** 2 / (vx * vx / (nx - 1) + vy * vy / (ny - 1))
    return se, df


# 97.5% point of Student t for the degrees of freedom met here (about 190)
def t975(df):
    from math import exp
    # Cornish-Fisher, ample for df > 30
    z = 1.959964
    return z + (z ** 3 + z) / (4 * df) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * df ** 2)


rows = []
for name, r in results.items():
    d = r["mean"] - b["mean"]
    if name == "base":
        se, df, hw = 0.0, 0.0, 0.0
    else:
        se, df = welch(r["y"], b["y"])
        hw = t975(df) * se
    rows.append(dict(name=name, description=r["desc"], budget=r["total"],
                     throughput=r["mean"], halfwidth=r["hw"],
                     diff_from_base=d, diff_halfwidth=hw))

with open(os.path.join(OUT, f"asymmetry{SUFFIX}.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
with open(os.path.join(OUT, f"asymmetry{SUFFIX}.json"), "w") as f:
    json.dump(dict(steps=STEPS, reps=REPS, seed0=SEED0, spacing=_geo.SPACING, loop=_geo.loop_length(4),
                   seconds=round(time.time() - t0, 1),
                   rows=rows,
                   raw={k: {kk: vv for kk, vv in v.items() if kk != "y"} for k, v in results.items()}),
              f, indent=1)

print("\n%-14s %-48s %9s %9s %10s %9s" % ("case", "description", "budget", "thrput", "diff", "+-95%"))
for r in rows:
    print("%-14s %-48s %9d %9.5f %+10.5f %9.5f" % (r["name"], r["description"][:48], r["budget"],
                                                  r["throughput"], r["diff_from_base"], r["diff_halfwidth"]))
print(f"\ndone in {time.time() - t0:.0f}s")
