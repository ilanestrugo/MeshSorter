#!/usr/bin/env python3
"""Does the estimated throughput depend on how much of the start of a run is
deleted?  (Table S3 of the Supplemental Online Material)

Each row is one configuration, simulated with 40 replications of 4 x 10^6 steps
in total.  The columns re-estimate the throughput from those same runs after
deleting the first 0, 10^3, 10^4, 10^5 and 10^6 steps.  A replication with seed s
follows the same trajectory whatever the warm-up, so deleting w steps from a run
of 4 x 10^6 is the same as running with warm-up w over the remaining
4 x 10^6 - w steps on the same seed, and the entries of a row share their random
numbers and differ only in the deletion.

The configurations are the seven of the supplement's pilot range: dual-drop,
loops of one common length, buffers of one capacity at every drop point.  The
loop length is given explicitly, on the crossings of the simulator's built-in
geometry, as in the pilot the table accompanies; the table is about how long the
transient lasts, not about the layout drawn in Figure 3.

    python3 warmup_sensitivity.py            about ten minutes
    python3 warmup_sensitivity.py --quick    a few replications, to check the wiring
"""
import csv, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
OUT = os.path.join(HERE, "results")

TOTAL, REPS, SEED = 4_000_000, 40, 20260901
DELETIONS = [0, 1_000, 10_000, 100_000, 1_000_000]
CASES = [(4, 4, 24, 0), (9, 9, 44, 0), (4, 9, 24, 0), (4, 9, 200, 0),
         (4, 9, 24, 10), (4, 9, 200, 10), (9, 9, 200, 10)]       # n, m, L, c
if "--quick" in sys.argv:
    TOTAL, REPS, DELETIONS, CASES = 400_000, 6, [0, 1_000, 10_000], CASES[:2]

os.makedirs(OUT, exist_ok=True)


def simulate(n, m, L, c, w):
    cmd = [BIN, "-n", str(n), "-m", str(m), "--dual", "-L", str(L),
           "-b", str(c), "-w", str(w), "-T", str(TOTAL - w), "-R", str(REPS),
           "--seed", str(SEED), "--json"]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("meshsorter failed: " + p.stderr.strip())
    return json.loads(p.stdout)


t0 = time.time()
rows = []
for (n, m, L, c) in CASES:
    row = dict(n=n, m=m, L=L, c=c)
    for w in DELETIONS:
        r = simulate(n, m, L, c, w)
        row[f"w{w}"] = r["throughput"]
        if w == 10_000:
            row["halfwidth"] = r["halfwidth"]
    rows.append(row)
    print("  n=%d m=%d L=%3d c=%2d  " % (n, m, L, c)
          + "  ".join("%.5f" % row[f"w{w}"] for w in DELETIONS)
          + "   hw %.5f" % row["halfwidth"], file=sys.stderr, flush=True)

with open(os.path.join(OUT, "warmup_sensitivity.csv"), "w", newline="") as f:
    w_ = csv.DictWriter(f, fieldnames=list(rows[0]))
    w_.writeheader()
    w_.writerows(rows)
with open(os.path.join(OUT, "warmup_sensitivity.tex"), "w") as f:
    for r in rows:
        f.write(f"{r['n']} & {r['m']} & {r['L']:3d} & {r['c']:2d} & "
                + " & ".join(f"{r[f'w{w}']:.5f}" for w in DELETIONS)
                + f" & {r['halfwidth']:.5f} \\\\\n")

# the two quantities the supplement quotes
spread0 = max(abs(r["w0"] - r["w1000000"]) for r in rows) if 1_000_000 in DELETIONS else None
shift0 = max(abs(r["w0"] - r["w10000"]) for r in rows)
shift4 = (max(abs(r["w10000"] - r["w1000000"]) for r in rows)
          if 1_000_000 in DELETIONS else None)
below = all(abs(r["w10000"] - r["w1000000"]) < r["halfwidth"] for r in rows) \
    if 1_000_000 in DELETIONS else None
with open(os.path.join(OUT, "warmup_sensitivity.json"), "w") as f:
    json.dump(dict(total=TOTAL, reps=REPS, seed=SEED, deletions=DELETIONS,
                   seconds=round(time.time() - t0, 1), rows=rows,
                   largest_shift_deleting_nothing_vs_1e4=shift0,
                   largest_shift_1e4_vs_1e6=shift4,
                   shift_1e4_vs_1e6_below_halfwidth_in_every_row=below), f, indent=1)
print(f"\ndeleting nothing instead of 1e4 steps shifts an estimate by at most {shift0:.5f}",
      file=sys.stderr)
if shift4 is not None:
    print(f"deleting 1e4 instead of 1e6 shifts it by at most {shift4:.5f}; "
          f"below the half-width in every row: {below}", file=sys.stderr)
print(f"done in {time.time() - t0:.0f}s", file=sys.stderr)
