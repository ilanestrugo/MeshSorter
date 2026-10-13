#!/usr/bin/env python3
"""What staggering and the second drop point are worth, with and without buffers.

Figure 7 of the manuscript (Section 6.2).  The 4 x 4 MeshSorter at the shortest
buildable loop length (104 slots), at the three per-primary-belt budgets of
Table S12, B = 0, 5 and 10.  Three gains, each in percent of the throughput of
the same system without the device in question:

    second drop point     dual-drop over single-drop, both with common loop lengths
    staggering, single    staggered over common loop lengths, single-drop
    staggering, dual      staggered over common loop lengths, dual-drop

Nothing is simulated here.  The inputs are the cells that loops_buffered.py
writes to results/loops_buffered.json (Table S12), so this script reads what
that table prints and cannot disagree with it.

Output
    results/lever_gains.json   the gains at each budget, and the range of each
                               staggering gain over all eight loop lengths
    results/lever_gains.tex    the pgfplots coordinates of Figure 7, with the label above each bar
  and, on standard output, the numbers quoted in Section 6.2.

    python3 lever_gains.py
"""
import json
import sys

import geometry as geo

SRC = "results/loops_buffered.json"
BUDGETS = (0, 5, 10)
L0 = geo.loop_length(4)                      # 104 slots, the buildable minimum

rows = json.load(open(SRC))
mu = {(r["dual"], r["B"], r["stagger"], r["L"]): r["mu"] for r in rows}
loops = sorted({r["L"] for r in rows})
if L0 not in loops:
    sys.exit("%s has no cells at the base loop length %d" % (SRC, L0))


def gain(new, old):
    """Relative gain in percent."""
    return 100.0 * (new / old - 1.0)


out = {"L": L0, "budgets": list(BUDGETS), "gain_percent": {}, "staggering_range_over_L": {}}
series = {"second_drop": [], "stagger_single": [], "stagger_dual": []}
print("loop length %d slots, gains in percent of the system without the device" % L0)
print("%4s %14s %16s %14s" % ("B", "second drop", "staggering 1-drop", "staggering 2-drop"))
for B in BUDGETS:
    g2 = gain(mu[(True, B, False, L0)], mu[(False, B, False, L0)])
    gs = gain(mu[(False, B, True, L0)], mu[(False, B, False, L0)])
    gd = gain(mu[(True, B, True, L0)], mu[(True, B, False, L0)])
    series["second_drop"].append((B, g2))
    series["stagger_single"].append((B, gs))
    series["stagger_dual"].append((B, gd))
    out["gain_percent"][str(B)] = {"second_drop": g2, "stagger_single": gs, "stagger_dual": gd}
    print("%4d %14.3f %16.4f %17.4f" % (B, g2, gs, gd))

print()
for B in BUDGETS:
    rng = {}
    for dual in (False, True):
        g = [gain(mu[(dual, B, True, L)], mu[(dual, B, False, L)]) for L in loops]
        rng["dual" if dual else "single"] = [min(g), max(g)]
    out["staggering_range_over_L"][str(B)] = rng
    print("B=%-2d staggering gain over the %d loop lengths: single %.3f..%.3f, dual %.3f..%.3f"
          % (B, len(loops), rng["single"][0], rng["single"][1], rng["dual"][0], rng["dual"][1]))
d0 = mu[(True, 0, False, L0)] - mu[(False, 0, False, L0)]
d10 = mu[(True, 10, False, L0)] - mu[(False, 10, False, L0)]
print("second drop point in items per time step: %.4f without buffers, %.4f at B=10 (factor %.1f)"
      % (d0, d10, d0 / d10))

json.dump(out, open("results/lever_gains.json", "w"), indent=1)


def label(v):
    """The value as printed above its bar."""
    return "%.1f" % v if v >= 1 else "%.2g" % v


with open("results/lever_gains.tex", "w") as f:
    f.write("% pgfplots coordinates of Figure 7, with the label printed above each bar:\n")
    f.write("% gain in percent against the per-primary-belt budget B\n")
    for name, pts in series.items():
        f.write("%% %s\n" % name)
        f.write("coordinates {" + " ".join("(%d,%.4g) [%s]" % (b, v, label(v)) for b, v in pts) + "}\n")
print("wrote results/lever_gains.json and results/lever_gains.tex")
