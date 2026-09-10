#!/usr/bin/env python3
"""Table 2 of the manuscript: steady-state throughput of the unbuffered
MeshSorter under the staggered configuration, with the gain over identical
loop lengths in parentheses.

    panel (a)  single-drop, staggered
    panel (b)  dual-drop, staggered, added slots on the return run
    panel (c)  dual-drop, staggered, added slots at the turnaround

The gains of panel (a) are measured against the exact single-drop values of
Table 1(a); those of panels (b) and (c) against the simulated dual-drop values
of Table 1(b), which this script recomputes so that the two tables are always
consistent with each other.

    python3 table2.py            full run, several hours on two cores
    python3 table2.py --quick    low precision, minutes
"""
import sys, math
import common
from exact_single_drop import exact

# ---------------------------------------------------------------- parameters
BELTS   = [4, 5, 6, 7, 8, 9]      # n
FEEDERS = [4, 5, 6, 7, 8, 9]      # m
DP      = 4
TURN    = 4
DF      = 4
WIDTH   = 2                       # base loop length is 4n + 8
STAGGER_STEP = 1                  # L_1 = L_2 = L, L_j = L + (j-2) for j >= 3
BUFFERS = None                    # unbuffered
STEPS   = 4_000_000
REPS    = 10
WARMUP  = None                    # None = the rule of Supplement S1
SEED    = 20260901
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, REPS = 200_000, 4

common.require_binary()
geo = dict(steps=STEPS, reps=REPS, warmup=WARMUP, seed=SEED,
           dp=DP, turn=TURN, df=DF, width=WIDTH, buffers=BUFFERS)

EXACT = {(n, m): float(exact(n, m)) for n in BELTS for m in FEEDERS}

print("reference: simulated, dual-drop, identical loop lengths [Table 1(b)]",
      file=sys.stderr)
REF = common.grid(BELTS, FEEDERS, dual=True, stagger=False, **geo)

print("panel (a): single-drop, staggered", file=sys.stderr)
A = common.grid(BELTS, FEEDERS, dual=False, stagger=True, extra="turnaround", **geo)

print("panel (b): dual-drop, staggered, slots on the return run", file=sys.stderr)
B = common.grid(BELTS, FEEDERS, dual=True, stagger=True, extra="return", **geo)

print("panel (c): dual-drop, staggered, slots at the turnaround", file=sys.stderr)
C = common.grid(BELTS, FEEDERS, dual=True, stagger=True, extra="turnaround", **geo)

common.save_raw("table1b_dual_identical", REF)
common.save_raw("table2a_single_staggered", A)
common.save_raw("table2b_dual_staggered_return", B)
common.save_raw("table2c_dual_staggered_turnaround", C)

body = []
body.append("%% panel (a), single-drop, gain over the exact identical-length values\n"
            + common.latex_rows(BELTS, FEEDERS,
                                lambda n, m: A[(n, m)]["throughput"],
                                lambda n, m: EXACT[(n, m)]))
for name, G in (("(b), added slots on the return run", B),
                ("(c), added slots at the turnaround", C)):
    body.append(f"%% panel {name}, gain over Table 1(b)\n"
                + common.latex_rows(BELTS, FEEDERS,
                                    lambda n, m, G=G: G[(n, m)]["throughput"],
                                    lambda n, m: REF[(n, m)]["throughput"]))
common.write("results/table2.tex", "\n".join(body) + "\n")

hw = common.max_halfwidth(A, B, C)
print(f"\nmaximum half-width over panels (a) to (c): {hw:.5f}")
print(f"warm-up used: {next(iter(A.values()))['warmup']} steps")

worse = [(n, m) for n in BELTS for m in FEEDERS
         if C[(n, m)]["throughput"] <= B[(n, m)]["throughput"]]
tcrit = 2.262 if REPS == 10 else 2.0
notsig = [(n, m) for n in BELTS for m in FEEDERS
          if C[(n, m)]["throughput"] - B[(n, m)]["throughput"] <=
             math.sqrt((C[(n,m)]["halfwidth"])**2 + (B[(n,m)]["halfwidth"])**2)]
print(f"cells where the turnaround placement does not beat the return run: {worse}")
print(f"cells where that advantage is not significant: {notsig}")
