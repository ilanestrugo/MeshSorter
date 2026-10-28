#!/usr/bin/env python3
"""Robustness of the unbuffered 4x4 dual-drop system to dependent destinations.

Section 6.4 of the main manuscript and the corresponding table of the
Supplemental Online Material.  An admitted item is bound for the primary belt
of the previous admission with probability rho, and otherwise for one of the
other three belts with equal probability (simulator option --rho).  The belt
marginals are uniform for every rho, and rho = 1/4 is the independent case, so
rho above 1/4 makes runs of items to one belt longer, with mean run length
1/(1-rho), and rho below 1/4 makes them shorter.

Two forms of dependence are run on the same grid:
  feeder   the previous admission at the same feeder; the feeders' streams
           are independent of one another;
  shared   the previous admission anywhere in the system, by time step and
           then by feeder index, so a run reaches several feeders at once.

Every cell is 30 replications of 1,330,000 measured steps after the warm-up
rule (20,000 steps here), on the base seed of the paper; the reference row is
the uniform draw of Table 2(b).

    python3 destination_dependence.py            full grid
    python3 destination_dependence.py --quick    wiring check

Output: results/destination_dependence.json and .tex.
"""
import argparse, hashlib, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
THREADS = os.environ.get("MESHSORTER_THREADS")

N, M = 4, 4
RHOS = [round(0.05 * k, 2) for k in range(20)] + [0.98, 0.99]
MODES = ("feeder", "shared")


def run(rho=None, shared=False, steps=1_330_000, reps=30, seed=20260901):
    cmd = [BIN, "-n", str(N), "-m", str(M), "--dual", "-T", str(steps),
           "-R", str(reps), "--seed", str(seed), "--json"]
    if rho is not None:
        cmd += ["--rho", f"{rho:.2f}"]
        if shared:
            cmd += ["--rho-shared"]
    key = hashlib.sha1(" ".join(cmd[1:]).encode()).hexdigest()[:16]
    cache = os.path.join(HERE, "results", "cache", key + ".json")
    if not os.environ.get("MESHSORTER_NOCACHE") and os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)
    if THREADS:
        cmd += ["-t", str(THREADS)]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"meshsorter failed: {out.stderr.strip()}")
    res = json.loads(out.stdout)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(res, f)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="short wiring check")
    a = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit(f"{BIN} not found.  Build it first:\n"
                 f"    c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp")
    kw = dict(steps=50_000, reps=4) if a.quick else {}
    rhos = [0.0, 0.25, 0.9] if a.quick else RHOS

    base = run(**kw)
    print(f"  uniform        {base['throughput']:.5f} +-{base['halfwidth']:.5f}",
          file=sys.stderr, flush=True)
    rows = []
    for rho in rhos:
        row = {"rho": rho, "mean_run": 1.0 / (1.0 - rho)}
        for mode in MODES:
            r = run(rho, mode == "shared", **kw)
            row[mode] = {"throughput": r["throughput"], "halfwidth": r["halfwidth"],
                         "change_pct": 100.0 * (r["throughput"] / base["throughput"] - 1.0),
                         "replications": r["replications"]}
        rows.append(row)
        print(f"  rho={rho:.2f}  feeder {row['feeder']['throughput']:.5f} "
              f"+-{row['feeder']['halfwidth']:.5f}   shared {row['shared']['throughput']:.5f} "
              f"+-{row['shared']['halfwidth']:.5f}", file=sys.stderr, flush=True)

    tag = "_quick" if a.quick else ""
    out = {"n": N, "m": M, "dual": True, "loops": base["loops"], "reps": base["reps"],
           "steps": base["steps"], "warmup": base["warmup"], "seed": base["seed"],
           "uniform": {"throughput": base["throughput"], "halfwidth": base["halfwidth"],
                       "replications": base["replications"]},
           "rows": rows}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", f"destination_dependence{tag}.json"), "w") as f:
        json.dump(out, f, indent=1)

    lines = [r"\begin{tabular}{rrcrcr}", r"\hline",
             r" & & \multicolumn{2}{c}{Separate runs per feeder} & \multicolumn{2}{c}{Runs shared by the feeders} \\",
             r"\cline{3-4}\cline{5-6}",
             r"$\rho$ & $1/(1-\rho)$ & Throughput & Change (\%) & Throughput & Change (\%) \\ \hline"]
    def pct(x):                     # no "-0.00"
        x = round(x, 2)
        return f"{x:+.2f}" if x else "0.00"

    for row in rows:
        f_, s_ = row["feeder"], row["shared"]
        lines.append(f"{row['rho']:.2f} & {row['mean_run']:.2f} & "
                     f"${f_['throughput']:.4f}\\pm{f_['halfwidth']:.4f}$ & ${pct(f_['change_pct'])}$ & "
                     f"${s_['throughput']:.4f}\\pm{s_['halfwidth']:.4f}$ & ${pct(s_['change_pct'])}$ \\\\")
    lines += [r"\hline", r"\end{tabular}", ""]
    with open(os.path.join(HERE, "results", f"destination_dependence{tag}.tex"), "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
