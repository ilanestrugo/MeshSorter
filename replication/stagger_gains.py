#!/usr/bin/env python3
"""The value of staggering the loop lengths, unbuffered, over the whole grid.

Figure 5 of the manuscript (Section 4.5).  For every system with n = 3..9
primary belts and m = 3..9 feeder loops it plots the gain of the staggered
loop lengths over identical ones, in percent, against the number of feeder
loops m, one line per n:

    panel (a)   single-drop, against the exact throughput of Table 2(a)
    panel (b)   dual-drop, added slots at the turnaround, against Table 2(b)

These are the percentages printed in parentheses in Table S5.  Nothing is
simulated here: the inputs are the cells that stagger_grid.py writes to
results/table2.csv, which carries the gain of each cell, so this script reads
what Table S5 prints and cannot disagree with it.

Output
    results/stagger_gains.json   the gain of every cell, by panel, n and m
    results/stagger_gains.tex    the pgfplots coordinates of Figure 5, one line
                                 per n and panel
  and, on standard output, the range of each panel quoted in Section 4.5.

    python3 stagger_gains.py
"""
import csv
import json

SRC = "results/table2.csv"
PANELS = {"single": "table2a_single_staggered",
          "dual": "table2c_dual_staggered_turnaround"}
BELTS = range(3, 10)
FEEDERS = range(3, 10)

rows = list(csv.DictReader(open(SRC)))
gain = {name: {} for name in PANELS}
for name, panel in PANELS.items():
    for r in rows:
        if r["panel"] == panel:
            gain[name][(int(r["n"]), int(r["m"]))] = float(r["gain_percent"])
    missing = [(n, m) for n in BELTS for m in FEEDERS if (n, m) not in gain[name]]
    if missing:
        raise SystemExit("%s: %s has no cells %s" % (SRC, panel, missing[:5]))

out = {name: {str(n): [gain[name][(n, m)] for m in FEEDERS] for n in BELTS} for name in PANELS}
out["feeder_loops"] = list(FEEDERS)
json.dump(out, open("results/stagger_gains.json", "w"), indent=1)

with open("results/stagger_gains.tex", "w") as f:
    f.write("% pgfplots coordinates of Figure 5: gain of staggering in percent against the number of feeder loops m\n")
    for name in PANELS:
        f.write("%% panel: %s-drop\n" % name)
        for n in BELTS:
            pts = " ".join("(%d,%.2f)" % (m, gain[name][(n, m)]) for m in FEEDERS)
            f.write("%% n = %d\ncoordinates {%s}\n" % (n, pts))

for name in PANELS:
    g = list(gain[name].values())
    print("%-6s %d cells, gain %.1f%% to %.1f%%, positive in all: %s"
          % (name, len(g), min(g), max(g), all(v > 0 for v in g)))
print("wrote results/stagger_gains.json and results/stagger_gains.tex")
