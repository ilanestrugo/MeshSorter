#!/usr/bin/env python3
"""Loop length and staggering under buffering, 4 x 4 system.

Produces the table that reports the throughput of the 4 x 4 MeshSorter as a
function of the base feeder-loop length, for three per-primary-belt budgets and
both drop mechanisms, with common and with staggered loop lengths.

Two questions.
  1. With buffers in place, does staggering the loop lengths still help?
  2. Does making all the loops longer help, and does that depend on the budget?

The buffered runs use the near-optimal allocations reported for n=4, m=4 in the
buffer-allocation table: c = (0,1,1,3) at B=5 and c = (0,1,2,7) at B=10, along
each primary belt, at the backward drop points in the dual-drop system and at
the single drop points in the single-drop system.

Output
    results/loops_buffered.json   one record per cell
    results/loops_buffered.tex    the body of the table, ready to paste
  and, on standard output, the comparison summary and the caption figures.

Protocol.  Thirty independent replications of 1,330,000 measured time steps,
the protocol of Supplement S1 and of the other tables in the manuscript.  The
earlier version of this script left `common.run` at its own defaults, ten
replications of 4,000,000 steps, which is the same total work but not the
protocol the manuscript states; pass --legacy to reproduce those numbers.

    python3 loops_buffered.py            # 96 cells, about 15 minutes on 8 cores
    python3 loops_buffered.py --quick    # low precision, a couple of minutes
"""
import sys, json, math
from fractions import Fraction
import common

# ---------------------------------------------------------------- parameters
BELTS, FEEDERS = 4, 4
LOOPS   = [24, 32, 40, 56, 72, 100, 140, 200]   # base feeder-loop length
BUDGETS = [0, 5, 10]                            # per-primary-belt budget B
ALLOC   = {0: None,                             # allocation used at each budget
           5: [0, 1, 1, 3],
           10: [0, 1, 2, 7]}
DP      = 4                       # loop slots per primary belt
TURN    = 4                       # loop slots per end of the loop
DF      = 4                       # feeder spacing along a primary belt
WIDTH   = 2                       # forward-to-backward distance on a belt
EXTRA   = "turnaround"            # where the staggered slots are inserted
STEPS   = 1_330_000               # measured steps per replication
REPS    = 30                      # independent replications
WARMUP  = None                    # None = the rule of Supplement S1
SEED    = 20260901
# ---------------------------------------------------------------------------

if "--legacy" in sys.argv:                      # the protocol of the first run
    STEPS, REPS = 4_000_000, 10
if "--quick" in sys.argv:
    STEPS, REPS = 100_000, 8

TQ = {8: 2.365, 10: 2.262, 30: 2.045}.get(REPS, 2.045)   # t_{0.975, R-1}


def buffers(dual, B):
    """The -b argument for budget B: one value per feeder in the single-drop
    system, 2*m*n values with the backward crossings buffered in the dual-drop
    system."""
    if B == 0:
        return None
    c = ALLOC[B]
    assert sum(c) == B, f"allocation {c} does not spend the budget {B}"
    if not dual:
        return c
    v = []
    for cj in c:
        v += [0] * BELTS + [cj] * BELTS
    return v


common.require_binary()
rows = []
for dual in (False, True):
    for B in BUDGETS:
        for stag in (False, True):
            for L in LOOPS:
                r = common.run(BELTS, FEEDERS, dual=dual, stagger=stag,
                               extra=EXTRA, buffers=buffers(dual, B), loop=L,
                               steps=STEPS, reps=REPS, warmup=WARMUP, seed=SEED,
                               dp=DP, turn=TURN, df=DF, width=WIDTH)
                rows.append(dict(dual=dual, B=B, stagger=stag, L=L,
                                 mu=r["throughput"], hw=r["halfwidth"],
                                 sd=r["sd"], warmup=r["warmup"],
                                 loops=r["loops"], reps=r["reps"],
                                 steps=r["steps"]))
                print(f"{'dual' if dual else 'single':6s} B={B:<3d} "
                      f"{'staggered' if stag else 'common   '} "
                      f"L={L:4d}  {r['throughput']:.5f} +-{r['halfwidth']:.5f}",
                      file=sys.stderr, flush=True)
json.dump(rows, open("results/loops_buffered.json", "w"), indent=1)

cell = {(r["dual"], r["B"], r["stagger"], r["L"]): r for r in rows}

# ------------------------------------------------------------ the comparison
print(f"\n{'mechanism':10s} {'budget':7s} {'L':>5s} {'common':>18s} "
      f"{'staggered':>18s} {'gain %':>8s}")
for dual in (False, True):
    for B in BUDGETS:
        for L in LOOPS:
            c = cell[(dual, B, False, L)]
            s = cell[(dual, B, True, L)]
            se = math.sqrt((c["sd"] ** 2 + s["sd"] ** 2) / REPS)
            star = "*" if abs(s["mu"] - c["mu"]) > TQ * se else " "
            print(f"{'dual' if dual else 'single':10s} B={B:<5d} "
                  f"{L:5d} {c['mu']:11.5f}+-{c['hw']:.5f} "
                  f"{s['mu']:11.5f}+-{s['hw']:.5f} "
                  f"{100*(s['mu']/c['mu']-1):7.2f}{star}")
print("\n* the difference between common and staggered exceeds a 95% interval")

# ------------------------------------------------------------- the LaTeX body
# The single-drop unbuffered system with a common loop length is exact and does
# not depend on L, so that column is one multirow entry rather than eight
# estimates.  Everything else is simulated.
EXACT_SD = Fraction(325, 124)
SKIP = (False, 0, False)          # (dual, B, stagger) of the exact column

body = []
for k, L in enumerate(LOOPS):
    cells = []
    for B in BUDGETS:
        for dual in (False, True):
            for stag in (False, True):
                if (dual, B, stag) == SKIP:
                    # one multirow entry, the exact value, rather than eight estimates
                    cells.append(
                        rf"\multirow{{{len(LOOPS)}}}{{*}}{{${float(EXACT_SD):.4f}^{{\dagger}}$}}"
                        if k == 0 else "")
                else:
                    cells.append(f"{cell[(dual, B, stag, L)]['mu']:.4f}")
    body.append(f"{L} & " + " & ".join(cells) + r" \\")
body = "\n".join(body)
common.write("results/loops_buffered.tex", body + "\n")

# The whole table, caption and all, ready to replace the table environment in the
# manuscript.  Thirteen columns do not fit at \footnotesize with the headings
# spelled out; measured against the 408pt text width of interact.cls, "Com." and
# "Stag." with \tabcolsep at 3pt come to 370pt and "Common" to 434pt.
hw = max(r["hw"] for r in rows if (r["dual"], r["B"], r["stagger"]) != SKIP)
parts = [rf"$\mathbf{{c}}=({','.join(str(x) for x in ALLOC[B])})$ at $B={B}$"
         for B in BUDGETS if B]
alloc = (" and ".join(parts) if len(parts) < 3
         else ", ".join(parts[:-1]) + " and " + parts[-1])
span = " & ".join(rf"\multicolumn{{4}}{{c}}{{{'No buffers' if B == 0 else f'$B={B}$'}}}"
                  for B in BUDGETS)
top = "".join(rf"\cline{{{2+4*i}-{5+4*i}}}" for i in range(len(BUDGETS)))
mech = " & ".join(r"\multicolumn{2}{c}{Single-drop} & \multicolumn{2}{c}{Dual-drop}"
                  for _ in BUDGETS)
mid = "".join(rf"\cline{{{2+2*i}-{3+2*i}}}" for i in range(2 * len(BUDGETS)))
head = " & ".join("Com. & Stag." for _ in range(2 * len(BUDGETS)))
table = rf"""\begin{{table}}[!htbp]
\centering
\caption{{\textcolor{{cyan}}{{Steady-state throughput of the ${BELTS}\times{FEEDERS}$ MeshSorter
as a function of the base feeder-loop length $L$, in items per time step, at
{'three' if len(BUDGETS) == 3 else len(BUDGETS)} per-primary-belt budgets. The buffered panels use the allocations
{alloc}, placed at the backward drop points in the dual-drop system. Com.\ denotes
a common loop length for every feeder and Stag.\ the staggered lengths
(\ref{{eq: staggered lengths}}). Reading a pair of columns left to right gives the
value of staggering, reading a column downward the value of a longer loop, and
reading across the panels the value of a budget. The largest half-width of a
$95\%$ confidence interval over the simulated entries is ${hw:.5f}$.
$^{{\dagger}}$This entry is not simulated: by
Proposition~\ref{{prop: exact single drop}} the unbuffered single-drop system with a
common loop length has throughput ${EXACT_SD.numerator}/{EXACT_SD.denominator}$
whatever the value of $L$}}}}
\label{{tab: loop length buffered}}
\footnotesize
\setlength{{\tabcolsep}}{{3pt}}
\begin{{tabular}}{{c {' '.join(['cccc'] * len(BUDGETS))}}}
\hline
 & {span}\\
{top}
 & {mech}\\
{mid}
$L$ & {head}\\
\hline
{body}
\hline
\end{{tabular}}
\end{{table}}
"""
common.write("results/loops_buffered_table.tex", table)

# ---------------------------------------------------------- caption figures
sim = [r for r in rows if (r["dual"], r["B"], r["stagger"]) != SKIP]
print(f"\nlargest half-width over the simulated entries: "
      f"{max(r['hw'] for r in sim):.5f}")
print(f"warm-up used: {sorted({r['warmup'] for r in rows})}, "
      f"{REPS} replications of {STEPS} steps")
dev = max(abs(cell[(False, 0, False, L)]["mu"] - float(EXACT_SD)) for L in LOOPS)
print(f"exact single-drop unbuffered value {EXACT_SD} = {float(EXACT_SD):.6f}; "
      f"simulated at the eight lengths anyway, largest deviation {dev:.5f}")

print("\nstaggering gain, smallest and largest over L (percent):")
for dual in (False, True):
    for B in BUDGETS:
        g = [100 * (cell[(dual, B, True, L)]["mu"] / cell[(dual, B, False, L)]["mu"] - 1)
             for L in LOOPS]
        print(f"  {'dual' if dual else 'single':6s} B={B:<3d} "
              f"{min(g):6.2f} to {max(g):6.2f}   (at L=200: {g[-1]:.2f})")

print("\nvalue of lengthening the loop from 24 to 200 (percent):")
for dual in (False, True):
    for B in BUDGETS:
        for stag in (False, True):
            a = cell[(dual, B, stag, LOOPS[0])]["mu"]
            b = cell[(dual, B, stag, LOOPS[-1])]["mu"]
            print(f"  {'dual' if dual else 'single':6s} B={B:<3d} "
                  f"{'staggered' if stag else 'common':9s} "
                  f"{a:.4f} -> {b:.4f}  {100*(b/a-1):5.2f}")
