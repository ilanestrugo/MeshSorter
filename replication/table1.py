#!/usr/bin/env python3
"""Table 1 of the manuscript: steady-state throughput of the unbuffered
MeshSorter with identical feeder-loop lengths.

    panel (a)  single-drop, computed exactly       (exact_single_drop.py)
    panel (b)  dual-drop, simulated                (meshsorter)

Every parameter of the experiment is written out below.  Run

    python3 table1.py

to reproduce the panel bodies in results/table1.tex and the caption figures
printed at the end.  Expect a few hours on two cores; pass --quick for a
low-precision run that finishes in minutes.
"""
import sys, itertools, math
import common
from exact_single_drop import exact

# ---------------------------------------------------------------- parameters
BELTS   = [4, 5, 6, 7, 8, 9]      # n
FEEDERS = [4, 5, 6, 7, 8, 9]      # m
DP      = 4                       # loop slots per primary belt
TURN    = 4                       # loop slots per end of the loop
DF      = 4                       # feeder spacing along a primary belt
WIDTH   = 2                       # forward-to-backward distance on a belt
                                  # loop length is therefore 4n + 8
BUFFERS = None                    # unbuffered
STEPS   = 4_000_000               # measured steps per replication
REPS    = 10                      # independent replications
WARMUP  = None                    # None = the rule of Supplement S1
SEED    = 20260901
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, REPS = 200_000, 4

common.require_binary()

print("panel (a): exact, single-drop, identical loop lengths", file=sys.stderr)
A = {(n, m): float(exact(n, m)) for n in BELTS for m in FEEDERS}

print("panel (b): simulated, dual-drop, identical loop lengths", file=sys.stderr)
B = common.grid(BELTS, FEEDERS, dual=True, stagger=False, buffers=BUFFERS,
                steps=STEPS, reps=REPS, warmup=WARMUP, seed=SEED,
                dp=DP, turn=TURN, df=DF, width=WIDTH)

body_a = common.latex_rows(BELTS, FEEDERS, lambda n, m: A[(n, m)])
body_b = common.latex_rows(BELTS, FEEDERS, lambda n, m: B[(n, m)]["throughput"])
hw = common.max_halfwidth(B)
warm = next(iter(B.values()))["warmup"]

common.save_raw("table1b_dual_identical", B)
common.write("results/table1.tex",
             "%% panel (a), single-drop, exact\n" + body_a +
             "\n%% panel (b), dual-drop, simulated\n" + body_b + "\n")

# ------------------------------------------------------- caption arithmetic
print(f"\nmaximum half-width in panel (b): {hw:.5f}")
print(f"warm-up used: {warm} steps, {REPS} replications of {STEPS} steps")

keys = list(B)
P = len(keys) * (len(keys) - 1) // 2
tcrit = 2.262 if REPS == 10 else None
zbonf = 3.95                       # two-sided normal quantile for 0.05/630
fails_pc = fails_bf = 0
closest = None
for a, b in itertools.combinations(keys, 2):
    sa = B[a]["halfwidth"] / (tcrit or 1) ; sb = B[b]["halfwidth"] / (tcrit or 1)
    se = math.sqrt(sa * sa + sb * sb)
    d = abs(B[a]["throughput"] - B[b]["throughput"])
    if d <= (tcrit or 2) * se: fails_pc += 1
    if d <= zbonf * se:        fails_bf += 1
    if closest is None or d < closest[0]: closest = (d, a, b, zbonf * se)
print(f"pairwise differences in panel (b): {P} pairs, "
      f"{fails_pc} not separated at 5%, {fails_bf} not separated under Bonferroni")
print(f"closest pair: {closest[1]} vs {closest[2]}, gap {closest[0]:.5f}, "
      f"Bonferroni threshold {closest[3]:.5f}")
