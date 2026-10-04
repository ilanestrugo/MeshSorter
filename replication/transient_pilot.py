#!/usr/bin/env python3
"""The transient pilot of Section S2 of the Supplemental Online Material.

For each of 48 configurations -- four system sizes (n, m), three loop lengths L
and four buffer capacities c, dual-drop, loops of one common length, the same
capacity at every drop point -- the C++ driver transient_pilot runs 300
independent replications of 2 x 10^5 steps from the empty state, deletes nothing,
and returns the ensemble mean and standard deviation of the throughput of every
block of 20 consecutive steps.  This script reads those series and reports, for
each configuration

  thr     the long-run level: the mean of the second half of the block means
  occ     the mean number of items held in buffers at the end of a run
  settle  the end of the last block whose ensemble mean differs from thr by more
          than 5.5 standard errors (SD/sqrt(300)) and by more than 3e-4
  A(w)    the excess of the ensemble mean over thr, summed over the steps from w
          onwards, in item-steps: a run that deletes w steps and averages over T
          more carries a bias of A(w)/T.  A is resolved only to about +-100 items.

and the three summary lines the supplement quotes.

    python3 transient_pilot.py            about half an hour on six cores
    python3 transient_pilot.py --quick    two configurations, 60 runs, to check the wiring
    python3 transient_pilot.py --from-cache   re-analyse series saved by an earlier run

The loop length -L is the base length on the crossings of the simulator's built-in
geometry, as in the pilot the supplement describes: the pilot is about how long
the transient lasts, not about the layout drawn in Figure 3.  The empty start is
what makes the transient: in the first round every feeder meets an empty slot on
every belt, so the series starts at its maximum, n, and decays.
"""
import json, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get("MESHSORTER_PILOT_BIN", os.path.join(HERE, "transient_pilot"))
OUT = os.path.join(HERE, "results", "transient_pilot")
CACHE = os.path.join(HERE, "results", "cache", "transient_pilot")

REPS, STEPS, BLOCK, SEED = 300, 200_000, 20, 20260901
Z, FLOOR = 5.5, 3e-4
SIZES = [(4, 4), (9, 9), (4, 9), (9, 4)]              # (n, m)
LOOPS = {4: 24, 9: 44}                                # built-in shortest loop, by n
CAPS = [0, 1, 3, 10]
WINDOWS = [0, 2_000, 20_000]                          # w in A(w)

cases = [(n, m, L, c) for (n, m) in SIZES for L in (LOOPS[n], 100, 200) for c in CAPS]
if "--quick" in sys.argv:
    REPS, cases = 60, [(4, 4, 24, 0), (4, 9, 24, 1)]


def run(n, m, L, c):
    path = os.path.join(CACHE, f"n{n}_m{m}_L{L}_c{c}.json")
    if "--from-cache" in sys.argv or os.path.exists(path) and "--quick" not in sys.argv:
        with open(path) as f:
            return json.load(f)
    cmd = [BIN, "-n", str(n), "-m", str(m), "--dual", "-L", str(L), "-b", str(c),
           "-R", str(REPS), "-T", str(STEPS), "--block", str(BLOCK), "--seed", str(SEED)]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("transient_pilot failed: " + p.stderr.strip())
    r = json.loads(p.stdout)
    if "--quick" not in sys.argv:
        os.makedirs(CACHE, exist_ok=True)
        with open(path, "w") as f:
            json.dump(r, f)
    return r


def analyse(r):
    mean, sd, reps, block = r["mean"], r["sd"], r["reps"], r["block"]
    B = len(mean)
    thr = sum(mean[B // 2:]) / (B - B // 2)
    se = [s / math.sqrt(reps) for s in sd]
    settle = 0
    for b in range(B):
        d = abs(mean[b] - thr)
        if d > Z * se[b] and d > FLOOR:
            settle = (b + 1) * block
    A = {}
    for w in WINDOWS:
        A[w] = sum(mean[b] - thr for b in range(w // block, B)) * block
    return dict(thr=thr, occ=r["occupancy_end"], settle=settle, A=A)


t0 = time.time()
rows = []
print("       n,m,L,c      thr      occ  settle     A(0)    A(2k)   A(20k)")
for (n, m, L, c) in cases:
    a = analyse(run(n, m, L, c))
    a.update(n=n, m=m, L=L, c=c)
    rows.append(a)
    print(f"{n:>7},{m},{L},{c}".rjust(14)
          + f"{a['thr']:9.4f}{a['occ']:9.2f}{a['settle']:8d}"
          + "".join(f"{a['A'][w]:9.1f}" for w in WINDOWS), flush=True)

slow = max(rows, key=lambda a: a["settle"])
within = sum(1 for a in rows if a["settle"] <= 1200)
slowest_within = max((a for a in rows if a["settle"] <= 1200), key=lambda a: a["settle"])
exceptions = [a for a in rows if a["settle"] > 1200]
a0 = max(rows, key=lambda a: abs(a["A"][0]))
a2 = max(rows, key=lambda a: abs(a["A"][2_000]))
summary = dict(
    reps=REPS, steps=STEPS, block=BLOCK, seed=SEED, configurations=len(rows),
    settle_within_1200=within,
    slowest_within_1200=dict(case=[slowest_within[k] for k in "nmLc"],
                             steps=slowest_within["settle"]),
    exceptions=[dict(case=[a[k] for k in "nmLc"], steps=a["settle"]) for a in exceptions],
    largest_cumulative_excess_from_0=dict(case=[a0[k] for k in "nmLc"], items=abs(a0["A"][0])),
    largest_cumulative_excess_from_2000=dict(case=[a2[k] for k in "nmLc"],
                                             items=abs(a2["A"][2_000])),
    seconds=round(time.time() - t0, 1))
print()
print(f"{within} of {len(rows)} settle within 1,200 steps; the slowest of them takes "
      f"{slowest_within['settle']} ({','.join(str(slowest_within[k]) for k in 'nmLc')})")
for a in exceptions:
    print(f"  exception: settle at {a['settle']} steps, {','.join(str(a[k]) for k in 'nmLc')}")
print(f"largest |cumulative excess| from t=0: {abs(a0['A'][0]):.0f} items")
print(f"largest |cumulative excess| from t=2,000 or later: {abs(a2['A'][2_000]):.0f} items")

if "--quick" not in sys.argv:
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "pilot_summary.json"), "w") as f:
        json.dump(dict(summary=summary,
                       rows=[dict(n=a["n"], m=a["m"], L=a["L"], c=a["c"], thr=a["thr"],
                                  occ=a["occ"], settle=a["settle"],
                                  A={str(w): v for w, v in a["A"].items()}) for a in rows]),
                  f, indent=1)
