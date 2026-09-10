#!/usr/bin/env python3
"""Section 5.1: how feeder-loop length interacts with buffering, 4 x 4 system.

Two questions.
  1. With buffers in place, does staggering the loop lengths still help?
  2. Does making all the loops longer help, with and without buffers?

The buffered runs use the near-optimal allocation reported for n=4, m=4, B=10:
c = (0,1,2,7) along each primary belt, at the backward drop points in the
dual-drop system and at the single drop points in the single-drop system.
"""
import sys, json, math
import common

BELTS, FEEDERS = 4, 4
LOOPS = [24, 32, 40, 56, 72, 100, 140, 200]
B_ALLOC = [0, 1, 2, 7]          # per-primary-belt budget B = 10

def buffers(dual, on):
    if not on: return None
    if not dual: return B_ALLOC                       # m values, one per feeder
    v = []
    for c in B_ALLOC:                                 # 2*m*n, forward then backward
        v += [0]*BELTS + [c]*BELTS
    return v

common.require_binary()
rows = []
for dual in (False, True):
    for on in (False, True):
        for stag in (False, True):
            for L in LOOPS:
                r = common.run(BELTS, FEEDERS, dual=dual, stagger=stag,
                               extra="turnaround", buffers=buffers(dual, on),
                               loop=L)
                rows.append(dict(dual=dual, buffered=on, stagger=stag, L=L,
                                 mu=r["throughput"], hw=r["halfwidth"],
                                 sd=r["sd"], warmup=r["warmup"],
                                 loops=r["loops"]))
                print(f"{'dual' if dual else 'single':6s} "
                      f"{'B=10' if on else 'B=0 ':4s} "
                      f"{'staggered' if stag else 'common   '} "
                      f"L={L:4d}  {r['throughput']:.5f} +-{r['halfwidth']:.5f}",
                      file=sys.stderr, flush=True)
json.dump(rows, open("results/loops_buffered.json", "w"), indent=1)

print(f"\n{'mechanism':10s} {'buffers':8s} {'L':>5s} {'common':>18s} {'staggered':>18s} {'gain %':>8s}")
for dual in (False, True):
    for on in (False, True):
        for L in LOOPS:
            c = next(r for r in rows if r["dual"]==dual and r["buffered"]==on
                     and not r["stagger"] and r["L"]==L)
            s = next(r for r in rows if r["dual"]==dual and r["buffered"]==on
                     and r["stagger"] and r["L"]==L)
            se = math.sqrt((c["sd"]**2 + s["sd"]**2)/30)
            star = "*" if abs(s["mu"]-c["mu"]) > 2.045*se else " "
            print(f"{'dual' if dual else 'single':10s} {'B=10' if on else 'B=0':8s} "
                  f"{L:5d} {c['mu']:11.5f}+-{c['hw']:.5f} {s['mu']:11.5f}+-{s['hw']:.5f} "
                  f"{100*(s['mu']/c['mu']-1):7.2f}{star}")
print("\n* the difference between common and staggered exceeds a 95% interval")
