#!/usr/bin/env python3
"""Check the simulator against the exact single-drop values.

The unbuffered single-drop system is solvable in closed form, so it gives an
independent reference for the simulator over the whole grid.  This script runs
the single-drop configuration for a set of cells and reports whether the exact
value falls inside the confidence interval, and by how many standard errors it
misses if it does not.

    python3 validate.py             a handful of cells, a few minutes
    python3 validate.py --full      the whole 6 by 6 grid
"""
import sys, math
import common
from exact_single_drop import exact

CELLS = [(4, 4), (4, 9), (6, 6), (9, 4), (9, 9)]
if "--full" in sys.argv:
    CELLS = [(n, m) for n in range(4, 10) for m in range(4, 10)]

STEPS, REPS, SEED = 4_000_000, 10, 20260901
TCRIT = 2.262                      # t_{9, 0.975}

common.require_binary()
print(f"{'n':>3} {'m':>3} {'exact':>12} {'simulated':>12} {'half-width':>11} "
      f"{'error/se':>9}  verdict")
worst = 0.0
for n, m in CELLS:
    e = float(exact(n, m))
    r = common.run(n, m, dual=False, stagger=False,
                   steps=STEPS, reps=REPS, seed=SEED)
    se = r["halfwidth"] / TCRIT
    z = (r["throughput"] - e) / se if se > 0 else 0.0
    worst = max(worst, abs(z))
    print(f"{n:>3} {m:>3} {e:>12.6f} {r['throughput']:>12.6f} "
          f"{r['halfwidth']:>11.6f} {z:>9.2f}  "
          f"{'inside' if abs(z) <= TCRIT else 'OUTSIDE'}")
print(f"\nlargest deviation: {worst:.2f} standard errors over {len(CELLS)} cells")
print("A few deviations beyond 2.26 standard errors are expected by chance;")
print("a systematic pattern of large ones would indicate a modelling error.")
