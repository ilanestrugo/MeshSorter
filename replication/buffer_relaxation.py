#!/usr/bin/env python3
"""How long does the buffer content take to settle?  (Section S2.2)

The transient pilot times the throughput series.  Buffers add a second, slower
transient: the number of items held in them first overshoots its steady-state
level and then returns to it.  This script measures that on the most heavily
loaded configuration of the pilot, four belts and nine feeders with loops of 416
slots and ten places at every drop point, by running 200 replications from the
empty state for a range of lengths and reporting the mean number of items held in
buffers when each run ends.  A run of length T started at step 0 is the first T
steps of a longer run on the same seed, so the figures are the occupancy at
those times.

    python3 buffer_relaxation.py            a few minutes
    python3 buffer_relaxation.py --quick    shorter runs, to check the wiring
"""
import json, os, subprocess, sys

import geometry as geo

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get("MESHSORTER_PILOT_BIN", os.path.join(HERE, "transient_pilot"))
OUT = os.path.join(HERE, "results")

N, M, L, C, REPS, SEED = 4, 9, 4 * geo.loop_length(4), 10, 200, 20260901
TIMES = [2_000, 10_000, 20_000, 100_000, 200_000, 400_000, 1_000_000, 2_000_000, 4_000_000]
if "--quick" in sys.argv:
    REPS, TIMES = 20, TIMES[:5]

rows = []
for T in TIMES:
    cmd = [BIN, "-n", str(N), "-m", str(M), "--dual", "-L", str(L), "-b", str(C),
           "--spacing", str(geo.SPACING), "--turn", str(geo.TURN),
           "--df", str(geo.DF), "--width", str(geo.WIDTH),
           "-R", str(REPS), "-T", str(T), "--block", str(T), "--seed", str(SEED)]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("transient_pilot failed: " + p.stderr.strip())
    occ = json.loads(p.stdout)["occupancy_end"]
    rows.append(dict(steps=T, items_in_buffers=occ))
    print(f"  T = {T:>9,}   items held in buffers {occ:8.1f}", flush=True)

steady = rows[-1]["items_in_buffers"]
peak = max(r["items_in_buffers"] for r in rows)
print(f"\npeak {peak:.1f}, level at the longest run {steady:.1f}: "
      f"overshoot {100 * (peak / steady - 1):.0f} percent")
for r in rows:
    r["relative_to_last"] = r["items_in_buffers"] / steady
if "--quick" not in sys.argv:
    with open(os.path.join(OUT, "buffer_relaxation.json"), "w") as f:
        json.dump(dict(n=N, m=M, loop=L, capacity=C, reps=REPS, seed=SEED, rows=rows,
                       overshoot_percent=100 * (peak / steady - 1)), f, indent=1)
