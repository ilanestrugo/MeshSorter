#!/usr/bin/env python3
"""Every number that appears in the prose of Section 7, taken from the data.

The manuscript keeps the prose with placeholders in place of the numbers, so
that a rerun of the grid regenerates the sentences rather than inviting a hand
edit.  Run this after approx_eval.py.

    python3 sec7_numbers.py            readable, with the phrase each belongs to
    python3 sec7_numbers.py --map      the substitutions, as JSON
    python3 sec7_numbers.py --tag large   read the fifteen-belt results instead
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
suffix = ""
if "--tag" in sys.argv:
    suffix = "_" + sys.argv[sys.argv.index("--tag") + 1]
path = os.path.join(HERE, "results", f"approx_eval{suffix}.json")
if not os.path.exists(path):
    sys.exit(f"{path} not found; run approx_eval.py first")
C = json.load(open(path))

pos = [c for c in C if c["B"] > 0]
if not pos:
    sys.exit("no classes with a positive budget")
tot = sum(c["full"]["count"] for c in pos)
offered = sum(c["offered"] for c in pos)
rhos = [c["full"]["rho"] for c in pos
        if c["full"]["rho"] == c["full"]["rho"] and c["full"]["count"] >= 10]
worst = max(pos, key=lambda c: c["full"]["loss"])


def vec(v):
    return ",".join(str(x) for x in v)


def mapping():
    w = worst["full"]
    bym, spread = {}, {}
    for m in sorted({c["m"] for c in pos if c["mech"] == "single"}):
        g = [c for c in pos if c["mech"] == "single" and c["m"] == m]
        n = sum(x["full"]["count"] for x in g)
        bym[m] = sum(x["full"]["gap_mean"] * x["full"]["count"] for x in g) / n
        spread[m] = (max(x["full"]["gap_mean"] for x in g)
                     - min(x["full"]["gap_mean"] for x in g))
    #  The classes the model orders worst, which are the ones the prose
    #  explains.  Restricted to classes big enough for a correlation to mean
    #  anything, and to those where that correlation is poor.
    sp = sorted((c["spread"], c["biasrange"]) for c in pos
                if c["full"]["count"] >= 10 and c["full"]["rho"] < 0.8)
    if not sp:
        sp = sorted((c["spread"], c["biasrange"]) for c in pos
                    if c["full"]["count"] >= 10)
    return {
        "@SPREADLO@": f"{min(x for x, _ in sp):.3f}" if sp else "--",
        "@SPREADHI@": f"{max(x for x, _ in sp):.3f}" if sp else "--",
        "@BIASLO@": f"{min(y for _, y in sp):.3f}" if sp else "--",
        "@BIASHI@": f"{max(y for _, y in sp):.3f}" if sp else "--",
        "@RHOMAX@": f"{max(rhos):.3f}" if rhos else "--",
        "@NWORST@": f"{len(sp)}",
        "@NPOS@": f"{tot:,}",
        "@NCLASS@": f"{len(pos)}",
        #  A number like 0.000068 is unreadable in running text.
        "@NOISE@": "$%.1f\\times10^{-5}$" % (
            1e5 * sum(c["full"]["noise"] for c in pos) / len(pos)),
        "@GAPMEAN@": f"{sum(c['full']['gap_mean'] * c['full']['count'] for c in pos) / tot:.3f}",
        "@GAPMAX@": f"{max(c['full']['gap_max'] for c in pos):.3f}",
        "@NNEG@": f"{sum(c['full']['negative'] for c in pos)}",
        "@BYM@": ", ".join(f"{g:.3f} at {m}" for m, g in sorted(bym.items())),
        "@SPREAD@": f"{max(spread.values()):.3f}",
        "@HITS@": f"{sum(1 for c in pos if c['full']['rank'] == 1)}",
        "@WORSTRANK@": f"{max(c['full']['rank'] for c in pos)}",
        "@LOSSMAX@": f"{w['loss']:.4f}",
        "@LOSSB@": f"{worst['B']}",
        "@LOSSM@": f"{worst['m']}",
        "@LOSSN@": f"{w['count']}",
        "@LOSSPICK@": vec(w["pick_cb"] if worst["mech"] == "dual" else w["pick_cf"]),
        "@LOSSBEST@": vec(w["best_cb"] if worst["mech"] == "dual" else w["best_cf"]),
        "@MODELDIFF@": f"{(w['pick_approx'] - w['best_approx']) / w['best_approx']:.4f}",
        "@LOSSSMALL@": f"{sum(1 for c in pos if c['full']['loss'] <= 0.001)}",
        "@RHOMIN@": f"{min(rhos):.3f}" if rhos else "--",
        "@NRANK@": f"{len(rhos)}",
    }


if "--map" in sys.argv:
    #  --prefix renames every key, so the two grids can be merged into one
    #  substitution table without their placeholders colliding.
    pre = ""
    if "--prefix" in sys.argv:
        pre = sys.argv[sys.argv.index("--prefix") + 1]
    m = {f"@{pre}{k[1:-1]}@": v for k, v in mapping().items()}
    print(json.dumps(m, indent=1))
    raise SystemExit


def show(label, value):
    print(f"  {label:<50s} {value}")


print(f"\nSection 7 prose, from {len(pos)} classes with a positive budget\n")
show("allocations in the structured class", f"{tot:,}")
show("allocations they were drawn from", f"{offered:,}")
show("mean |pass 1 - pass 2|, relative",
     f"{sum(c['full']['noise'] for c in pos) / len(pos):.6f}")
print()
show("mean gap, relative",
     f"{sum(c['full']['gap_mean'] * c['full']['count'] for c in pos) / tot:.4f}")
show("largest gap, relative", f"{max(c['full']['gap_max'] for c in pos):.4f}")
show("allocations underestimated", f"{sum(c['full']['negative'] for c in pos)}")
print()
show("classes where the model picks the simulated best",
     f"{sum(1 for c in pos if c['full']['rank'] == 1)} of {len(pos)}")
show("worst rank of the model's choice",
     f"{max(c['full']['rank'] for c in pos)}")
show("largest selection loss, relative", f"{worst['full']['loss']:.5f}")
show("  attained at",
     f"{worst['mech']} m={worst['m']} B={worst['B']} "
     f"({worst['full']['count']} allocations in S)")
show("  the model's choice against the simulated best",
     f"{vec(worst['full']['pick_cf'])} vs {vec(worst['full']['best_cf'])}")
show("selection loss within the indifference zone of 0.001",
     f"{sum(1 for c in pos if c['full']['loss'] <= 0.001)} of {len(pos)}")
print()
show("smallest Spearman over classes of ten or more",
     f"{min(rhos):.4f} ({len(rhos)} classes)" if rhos else "no class is large enough")
print()
