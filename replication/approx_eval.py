#!/usr/bin/env python3
"""Section 7: accuracy of the approximation model on the structured class.

The design rules of Section 5.1 define a small class S of buffer allocations,
and the certification reported there establishes that S contains an allocation
that is optimal up to the indifference zone.  A designer therefore searches S,
not the whole allocation space.  This script asks how well the approximation
model of Section 6 serves that search: whether it predicts the throughput of an
allocation in S, whether it ranks the members of S as the simulation does, and
what is lost by building the one it likes best.

Input
    <dir>/*_allocs.csv    written by sweep_structured.py, one row per
                          allocation, with the two independent simulation
                          estimates of its throughput.  certify_all.py writes
                          the same format for the whole allocation space.

Output
    results/approx_eval<tag>.csv        one row per class
    results/approx_eval<tag>.tex        LaTeX body, budgets collapsed
    results/approx_eval_full<tag>.tex   LaTeX body, class by class
    results/approx_eval<tag>.json       everything, for further analysis
    results/approx_scatter<tag>_*.dat   pgfplots data

Usage
    python3 approx_eval.py                       the four-belt grid of Table 7
    python3 approx_eval.py --dir results/large --tag large --scatter 15
    python3 approx_eval.py --all                 include allocations outside S

A class is one combination of drop mechanism, number of feeder loops and buffer
budget.  The unbuffered system is left out of the tables: it holds a single
allocation, so nothing in it can be ranked or selected, and Section 4 solves it
exactly, so the approximation is never used there.  Its numbers stay in the CSV.
"""

import csv, glob, json, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import approx_model

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")

#  Smallest class for which a reported rank correlation means anything.
RANKABLE = 10


def finite(values):
    """Drop the undefined entries, which come from classes too small to rank."""
    return [x for x in values if x == x]


# ------------------------------------------------------------------ ranking --
def ranks(v):
    """Ranks of v, smallest first, ties sharing their average rank."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    """Spearman rank correlation, as Pearson on the ranks so that ties are
    handled correctly."""
    if len(x) < 3:
        return float("nan")
    a, b = ranks(x), ranks(y)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    sa = sum((v - ma) ** 2 for v in a)
    sb = sum((v - mb) ** 2 for v in b)
    if sa <= 0 or sb <= 0:
        return float("nan")
    cov = sum((a[i] - ma) * (b[i] - mb) for i in range(len(a)))
    return cov / math.sqrt(sa * sb)


# --------------------------------------------------------------- one class ---
def read_cell(path):
    """Read one allocation dump and attach the approximation to every row."""
    rows = []
    with open(path, newline="") as f:
        for d in csv.DictReader(f):
            n, m = int(d["n"]), int(d["m"])
            dual = d["dual"] == "1"
            cf = [int(x) for x in d["c_forward"].split("|")]
            cb = [int(x) for x in d["c_backward"].split("|")] if d["c_backward"] \
                else None
            r1, r3 = int(d["r1"]), int(d["r3"])
            m1, m3 = float(d["mean1"]), float(d["mean3"])
            s1, s3 = float(d["sd1"]), float(d["sd3"])

            #  Pooled over the two independent passes.  Their difference is the
            #  only measure of simulation noise that costs nothing to obtain.
            sim = (r1 * m1 + r3 * m3) / (r1 + r3)
            var = ((r1 - 1) * s1 * s1 + (r3 - 1) * s3 * s3) / (r1 + r3 - 2)

            rows.append(dict(
                n=n, m=m, B=int(d["B"]), dual=dual, cf=cf, cb=cb,
                structured=d["structured"] == "1",
                sim=sim, se=math.sqrt(var / (r1 + r3)), half_diff=abs(m1 - m3),
                approx=approx_model.throughput(n, m, cf, cb if dual else None)))
    return rows


def summarize(rows):
    """Accuracy, ranking and selection statistics over one class."""
    if not rows:
        return None
    #  Relative throughout, on the scale of the indifference zone eps = 0.001
    #  that the rest of the paper uses, never as a percentage.
    gaps = [(x["approx"] - x["sim"]) / x["sim"] for x in rows]
    best = max(rows, key=lambda x: x["sim"])
    pick = max(rows, key=lambda x: x["approx"])
    order = sorted(rows, key=lambda x: -x["sim"])
    return dict(
        count=len(rows),
        gap_mean=sum(gaps) / len(gaps),
        gap_max=max(gaps, key=abs),
        gap_min=min(gaps),
        negative=sum(1 for g in gaps if g < 0),
        rho=spearman([x["approx"] for x in rows], [x["sim"] for x in rows]),
        best_sim=best["sim"], pick_sim=pick["sim"],
        best_approx=best["approx"], pick_approx=pick["approx"],
        loss=(best["sim"] - pick["sim"]) / best["sim"],
        rank=next(i for i, x in enumerate(order) if x is pick) + 1,
        pick_cf=pick["cf"], pick_cb=pick["cb"],
        best_cf=best["cf"], best_cb=best["cb"],
        se=sum(x["se"] for x in rows) / len(rows),
        noise=sum(x["half_diff"] for x in rows) / len(rows)
        / (sum(x["sim"] for x in rows) / len(rows)))


# ------------------------------------------------------------------- driver --
def main():
    src = os.path.join(OUT, "structured")
    if "--dir" in sys.argv:
        src = sys.argv[sys.argv.index("--dir") + 1]
        if not os.path.isabs(src):
            src = os.path.join(HERE, src)
    #  Not named "tag": that name is taken by a loop variable below, and the
    #  two silently collided once already.
    suffix = ""
    if "--tag" in sys.argv:
        suffix = "_" + sys.argv[sys.argv.index("--tag") + 1]
    scatter_at = [5, 10]
    if "--scatter" in sys.argv:
        scatter_at = [int(x) for x in
                      sys.argv[sys.argv.index("--scatter") + 1].split(",")]
    everything = "--all" in sys.argv

    files = sorted(glob.glob(os.path.join(src, "*_allocs.csv")))
    if not files:
        sys.exit(f"no allocation dumps in {src}.\n"
                 "Run sweep_structured.py first; it writes one next to the JSON\n"
                 "of every cell.  certify_all.py writes them too, for the whole\n"
                 "allocation space rather than the structured class: read those\n"
                 "with --dir results/certify.")

    classes, cells = [], {}
    for path in files:
        rows = read_cell(path)
        if not rows:
            continue
        head = rows[0]
        kept = rows if everything else [r for r in rows if r["structured"]]
        if not kept:
            continue
        cells[(head["dual"], head["m"], head["B"])] = kept
        #  How widely the class is spread, and how much the model's bias
        #  varies across it.  The first has to exceed the second for the model
        #  to order the class, which is the condition Section 7 reports.
        mean = sum(x["sim"] for x in kept) / len(kept)
        spread = (max(x["sim"] for x in kept)
                  - min(x["sim"] for x in kept)) / mean
        g = [(x["approx"] - x["sim"]) / x["sim"] for x in kept]
        classes.append(dict(mech="dual" if head["dual"] else "single",
                            n=head["n"], m=head["m"], B=head["B"],
                            full=summarize(kept), offered=len(rows),
                            spread=spread, biasrange=max(g) - min(g)))
    classes.sort(key=lambda c: (c["mech"], c["m"], c["B"]))

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"approx_eval{suffix}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mechanism", "n", "m", "B", "allocations", "gap_mean",
                    "gap_max", "spearman", "rank_of_choice",
                    "selection_loss", "sim_noise"])
        for c in classes:
            a = c["full"]
            w.writerow([c["mech"], c["n"], c["m"], c["B"], a["count"],
                        f"{a['gap_mean']:.6f}", f"{a['gap_max']:.6f}",
                        "" if a["rho"] != a["rho"] else f"{a['rho']:.4f}",
                        a["rank"], f"{a['loss']:.6f}", f"{a['noise']:.6f}"])
    with open(os.path.join(OUT, f"approx_eval{suffix}.json"), "w") as f:
        json.dump(classes, f, indent=1)

    # --------------------------------------------------------------- LaTeX --
    lines, full = [], []
    for mech in ("single", "dual"):
        for m in sorted({c["m"] for c in classes if c["mech"] == mech}):
            sel = [c for c in classes if c["mech"] == mech and c["m"] == m
                   and c["B"] > 0]
            if not sel:
                continue
            tot = sum(c["full"]["count"] for c in sel)
            gm = sum(c["full"]["gap_mean"] * c["full"]["count"] for c in sel) / tot
            #  The collapsed table carries no rank correlation.  At four
            #  primary belts a class of S holds between one and 23
            #  allocations, and a correlation over a handful of points says
            #  nothing; the class-by-class body reports it where the class is
            #  large enough, which is everywhere in the fifteen-belt grid.
            #  What is meaningful at any size is where the model's choice
            #  stands in the simulated ordering.
            hits = sum(1 for c in sel if c["full"]["rank"] == 1)
            lines.append(
                f"{mech.capitalize()} & {m} & {tot} & {gm:.4f} & "
                f"{max(c['full']['gap_max'] for c in sel):.4f} & "
                f"{hits}/{len(sel)} & "
                f"{max(c['full']['rank'] for c in sel)} & "
                f"{max(c['full']['loss'] for c in sel):.5f} " + r"\\")
            for c in sel:
                a = c["full"]
                r = a["rho"]
                rtx = "--" if r != r else f"{r:.4f}"
                full.append(
                    f"{mech.capitalize()} & {m} & {c['B']} & {a['count']} & "
                    f"{a['gap_mean']:.4f} & {a['gap_max']:.4f} & {rtx} & "
                    f"{a['rank']} & {a['loss']:.5f} " + r"\\")
        lines.append(r"\hline")
        full.append(r"\hline")
    with open(os.path.join(OUT, f"approx_eval{suffix}.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT, f"approx_eval_full{suffix}.tex"), "w") as f:
        f.write("\n".join(full) + "\n")

    # ------------------------------------------------------------- scatter --
    for (dual, m, B), rows in sorted(cells.items()):
        if B not in scatter_at:
            continue
        mech = "dual" if dual else "single"
        path = os.path.join(OUT, f"approx_scatter{suffix}_{mech}_B{B}_m{m}.dat")
        with open(path, "w") as f:
            f.write("Model_Utilization\tSim_Utilization\n")
            for r in sorted(rows, key=lambda r: r["sim"]):
                f.write(f"{r['approx'] / r['n']:.6f}\t{r['sim'] / r['n']:.6f}\n")
        print(f"  {os.path.basename(path)}: {len(rows)} allocations")

    report(classes, everything)


def report(classes, everything):
    print()
    print("over " + ("every allocation" if everything
                     else "the structured class only"))
    print(f"{'mech':7s} {'m':>2s} {'B':>3s} {'in S':>6s} {'of':>7s} {'gap':>8s} "
          f"{'max':>8s} {'rho':>8s} {'rank':>5s} {'loss':>8s}")
    for c in classes:
        a = c["full"]
        print(f"{c['mech']:7s} {c['m']:2d} {c['B']:3d} {a['count']:6d} "
              f"{c['offered']:7d} {a['gap_mean']:8.4f} {a['gap_max']:8.4f} "
              f"{a['rho']:8.4f} {a['rank']:5d} {a['loss']:8.5f}")
    pos = [c for c in classes if c["B"] > 0]
    if not pos:
        return
    tot = sum(c["full"]["count"] for c in pos)
    rhos = finite([c["full"]["rho"] for c in pos
                   if c["full"]["count"] >= RANKABLE])
    print(f"\n{len(pos)} classes with a positive budget, {tot} allocations")
    print(f"  mean gap {sum(c['full']['gap_mean'] * c['full']['count'] for c in pos) / tot:.4f}, "
          f"largest {max(c['full']['gap_max'] for c in pos):.4f}, "
          f"smallest {min(c['full']['gap_min'] for c in pos):.4f}")
    print(f"  allocations the model underestimates: "
          f"{sum(c['full']['negative'] for c in pos)}")
    if rhos:
        print(f"  Spearman between {min(rhos):.4f} and {max(rhos):.4f} "
              f"over the {len(rhos)} classes of {RANKABLE} or more")
    print(f"  the model picks the simulated best in "
          f"{sum(1 for c in pos if c['full']['rank'] == 1)} of {len(pos)} classes")
    print(f"  selection loss at most {max(c['full']['loss'] for c in pos):.5f}, "
          f"below the indifference zone of 0.001 in "
          f"{sum(1 for c in pos if c['full']['loss'] <= 0.001)}")
    print(f"  simulation noise, mean |pass 1 - pass 2| "
          f"{sum(c['full']['noise'] for c in pos) / len(pos):.6f}")


if __name__ == "__main__":
    main()
