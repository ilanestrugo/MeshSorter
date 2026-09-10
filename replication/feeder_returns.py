#!/usr/bin/env python3
"""Diminishing returns from additional feeder loops, four primary belts.

Replaces the old Figure 5.  Both mechanisms, staggered loop lengths, m = 1..10.
The single-drop identical-length values are added for reference; they are exact.
"""
import sys, json
import common
from exact_single_drop import exact

N = 4
MS = list(range(1, 11))
common.require_binary()
rows = []
for dual in (False, True):
    for m in MS:
        r = common.run(N, m, dual=dual, stagger=True, extra="turnaround")
        rows.append(dict(dual=dual, m=m, mu=r["throughput"], hw=r["halfwidth"]))
        print(f"{'dual' if dual else 'single':6s} m={m:2d}  {r['throughput']:.5f} "
              f"+-{r['halfwidth']:.5f}", file=sys.stderr, flush=True)
json.dump(rows, open("results/feeder_returns.json", "w"), indent=1)

print(f"\n{'m':>3} {'single stag':>12} {'single exact':>13} {'dual stag':>11} "
      f"{'util single':>12} {'util dual':>10} {'marginal s':>11} {'marginal d':>11}")
prev_s = prev_d = 0.0
for m in MS:
    s = next(r["mu"] for r in rows if not r["dual"] and r["m"] == m)
    d = next(r["mu"] for r in rows if r["dual"] and r["m"] == m)
    e = float(exact(N, m))
    print(f"{m:>3} {s:12.5f} {e:13.5f} {d:11.5f} {100*s/N:11.1f}% {100*d/N:9.1f}% "
          f"{s-prev_s:11.4f} {d-prev_d:11.4f}")
    prev_s, prev_d = s, d
