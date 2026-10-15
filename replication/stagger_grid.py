#!/usr/bin/env python3
"""Tables S5 and S6 of the Supplemental Online Material: steady-state throughput of the unbuffered
MeshSorter under the staggered configuration.

    panel (a)  single-drop, staggered                              Table S5(a)
    panel (b)  dual-drop, staggered, added slots on the return run   Table S6
    panel (c)  dual-drop, staggered, added slots at the turnaround   Table S5(b)

In Table S5 the parentheses give the gain over identical loop lengths: for
panel (a) against the exact single-drop values of Table 2(a), for panel (c)
against the simulated dual-drop values of Table 2(b), which this script
recomputes so that the tables are always consistent with each other.  In
Table S6 they give the shortfall of the return-run placement from the
turnaround placement of the same cell, in percent with two decimals, since
the smallest is a few hundredths of a percent.  results/table2.csv keeps the
gain over identical loop lengths for all three panels.

    python3 stagger_grid.py            full run, several hours on two cores
    python3 stagger_grid.py --quick    low precision, minutes
"""
import sys, math
import common
from exact_single_drop import exact

# ---------------------------------------------------------------- parameters
BELTS   = [3, 4, 5, 6, 7, 8, 9]      # n
FEEDERS = [3, 4, 5, 6, 7, 8, 9]      # m
import geometry as _geo
DP      = _geo.DP
TURN    = _geo.TURN
DF      = _geo.DF
WIDTH   = _geo.WIDTH                    # base loop length is 20n + 24
STAGGER_STEP = 1                  # L_1 = L_2 = L, L_j = L + (j-2) for j >= 3
BUFFERS = None                    # unbuffered
STEPS   = 1_330_000
REPS    = 30
WARMUP  = None                    # None = the rule of Section S2
SEED    = 20260901
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, REPS = 100_000, 8

common.require_binary()
geo = dict(steps=STEPS, reps=REPS, warmup=WARMUP, seed=SEED,
           dp=DP, turn=TURN, df=DF, width=WIDTH, buffers=BUFFERS)

EXACT = {(n, m): float(exact(n, m)) for n in BELTS for m in FEEDERS}

print("reference: simulated, dual-drop, identical loop lengths [Table 2(b)]",
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

common.write_table_csv("results/table2.csv",
                       [("table2a_single_staggered", A, EXACT),
                        ("table2b_dual_staggered_return", B,
                         {k: REF[k]["throughput"] for k in REF}),
                        ("table2c_dual_staggered_turnaround", C,
                         {k: REF[k]["throughput"] for k in REF}),
                        ("table1b_dual_identical_reference", REF, None)])

body = []
body.append("%% panel (a), single-drop, gain over the exact identical-length values\n"
            + common.latex_rows(BELTS, FEEDERS,
                                lambda n, m: A[(n, m)]["throughput"],
                                lambda n, m: EXACT[(n, m)]))
rows_b = []
for n in BELTS:
    cells = []
    for m in FEEDERS:
        b, c = B[(n, m)]["throughput"], C[(n, m)]["throughput"]
        cells.append(f"{b:.4f} ({100 * (1 - b / c):.2f})")
    rows_b.append(f"{n} & " + " & ".join(cells) + r"\\")
body.append("%% panel (b), added slots on the return run (Table S6), "
            "shortfall from panel (c) in percent\n" + "\n".join(rows_b))
body.append("%% panel (c), added slots at the turnaround, gain over Table 2(b)\n"
            + common.latex_rows(BELTS, FEEDERS,
                                lambda n, m: C[(n, m)]["throughput"],
                                lambda n, m: REF[(n, m)]["throughput"]))
common.write("results/table2.tex", "\n".join(body) + "\n")

hw = common.max_halfwidth(A, B, C)
print(f"\nmaximum half-width over panels (a) to (c): {hw:.5f}")
print(f"warm-up used: {next(iter(A.values()))['warmup']} steps")

worse = [(n, m) for n in BELTS for m in FEEDERS
         if C[(n, m)]["throughput"] <= B[(n, m)]["throughput"]]
tcrit = {10: 2.262, 20: 2.093, 30: 2.045, 40: 2.023}.get(REPS, 2.0)
notsig = [(n, m) for n in BELTS for m in FEEDERS
          if C[(n, m)]["throughput"] - B[(n, m)]["throughput"] <=
             math.sqrt((C[(n,m)]["halfwidth"])**2 + (B[(n,m)]["halfwidth"])**2)]
print(f"cells where the turnaround placement does not beat the return run: {worse}")
print(f"cells where that advantage is not significant: {notsig}")
