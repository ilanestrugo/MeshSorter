#!/usr/bin/env python3
"""Shuffled-sequence control for the order-derived check (Supplement S5).

The order-derived check replaces independent uniform destination draws by one
chronological sequence of 98,816 labels, read cyclically and consumed by the
feeders in index order.  That changes several things at once:

    1. the marginal share of each belt (0.2494, 0.2500, 0.2500, 0.2506 instead
       of 0.25 exactly),
    2. the fact that the destinations come from one sequence that repeats
       after 98,816 labels, about fifty times in a replication,
    3. the way the four feeders share that sequence, in index order, so that
       the labels the feeders hold at one instant are neighbours in it,
    4. the chronological dependence of the orders themselves.

To isolate (4) this script keeps (1) to (3) and removes the chronology.  Each
control sequence is a uniformly random permutation of the same labels, so it
has exactly the same counts, the same length, the same cyclic structure and is
consumed in the same way, and replications again differ in where they enter it.
Many independent permutations are used, so that what is reported is the
distribution of the throughput over shuffled sequences and not the accident of
one shuffle.

    python3 order_shuffled.py              the control, K = 30 shuffles
    python3 order_shuffled.py --shuffles K another number of shuffles
    python3 order_shuffled.py --budgets 10 only the listed budgets (comma list)

The allocations are those of results/table6.csv, the ones the order-derived
table used, so that the three columns of the comparison share them.
"""
import csv, json, os, random, statistics, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
OUT = os.path.join(HERE, "results")
SEQ = os.path.join(OUT, "order_derived_sequence.txt")
SHUF = os.path.join(OUT, "shuffled")

import geometry as _geo

BELTS, FEEDERS, DUAL = 4, 4, True
STEPS, REPS, SEED = 1_330_000, 30, 20260901
K = 30
BUDGETS = list(range(0, 11))
if "--shuffles" in sys.argv:
    K = int(sys.argv[sys.argv.index("--shuffles") + 1])
if "--budgets" in sys.argv:
    BUDGETS = [int(x) for x in sys.argv[sys.argv.index("--budgets") + 1].split(",")]
SHUFFLE_SEED = 20261004          # fixed, so the control is reproducible
TAG = "" if len(BUDGETS) == 11 else "_" + "_".join(map(str, BUDGETS))


def read_sequence():
    with open(SEQ) as f:
        return [int(l) for l in f if l.strip() and not l.startswith("#")]


def buffer_spec(alloc):
    c = [int(x) for x in alloc.split(",")]
    out = []
    for cj in c:
        out += [0] * BELTS + [cj] * BELTS
    return ",".join(map(str, out))


def simulate(alloc, path):
    cmd = [BIN, "-n", str(BELTS), "-m", str(FEEDERS), "--dual",
           "--spacing", str(_geo.SPACING), "--turn", str(_geo.TURN),
           "--df", str(_geo.DF), "--width", str(_geo.WIDTH),
           "-b", buffer_spec(alloc), "-T", str(STEPS), "-R", str(REPS),
           "--seed", str(SEED), "--sequence", path, "--json"]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("meshsorter failed: " + p.stderr.strip())
    return json.loads(p.stdout)


def write_tex(rows):
    """The body of the supplement's table: one row per budget."""
    with open(os.path.join(OUT, f"order_shuffled{TAG}.tex"), "w") as f:
        for r in rows:
            f.write(f"{r['B']} & {r['alloc']} & {float(r['th_uniform']):.4f} & "
                    f"{float(r['th_order']):.4f} & {float(r['th_shuffled_mean']):.4f} & "
                    f"{float(r['th_shuffled_sd']):.4f} & "
                    f"{float(r['th_shuffled_min']):.4f}--{float(r['th_shuffled_max']):.4f} & "
                    f"{r['shuffles_at_or_below_order']} \\\\\n")
        f.write("\\hline\n")


if "--tex-only" in sys.argv:           # rebuild the LaTeX body from the saved CSV
    write_tex(list(csv.DictReader(open(os.path.join(OUT, f"order_shuffled{TAG}.csv")))))
    sys.exit(0)

labels = read_sequence()
base = {r["B"]: r for r in csv.DictReader(open(os.path.join(OUT, "table6.csv")))}
os.makedirs(SHUF, exist_ok=True)

paths = []
for k in range(K):
    rng = random.Random(SHUFFLE_SEED + k)
    perm = labels[:]
    rng.shuffle(perm)
    assert sorted(perm) == sorted(labels)
    p = os.path.join(SHUF, f"shuffle_{k:03d}.txt")
    with open(p, "w") as f:
        f.write(f"# permutation {k} of the order-derived sequence, seed {SHUFFLE_SEED + k}\n")
        f.write("\n".join(map(str, perm)) + "\n")
    paths.append(p)

print(f"{len(labels)} labels, {K} shuffles, budgets {BUDGETS}", file=sys.stderr)
rows = []
t0 = time.time()
for B in BUDGETS:
    a = base[str(B)]["alloc"]
    th_order = float(base[str(B)]["th_order"])
    th_unif = float(base[str(B)]["th_uniform"])
    means = []
    for p in paths:
        means.append(simulate(a, p)["throughput"])
    m = statistics.mean(means)
    sd = statistics.stdev(means) if K > 1 else 0.0
    below = sum(1 for x in means if x <= th_order)     # shuffles at or below the real sequence
    rows.append(dict(B=B, alloc=a, th_uniform=th_unif, th_order=th_order,
                     th_shuffled_mean=m, th_shuffled_sd=sd, th_shuffled_min=min(means),
                     th_shuffled_max=max(means), shuffles=K,
                     shuffles_at_or_below_order=below,
                     order_minus_shuffled=th_order - m,
                     shuffled_minus_uniform=m - th_unif))
    print(f"  B={B:2d} {a:10s} uniform {th_unif:.4f}  shuffled {m:.4f} (sd {sd:.4f}, "
          f"{min(means):.4f}..{max(means):.4f})  real order {th_order:.4f}  "
          f"{below}/{K} shuffles at or below it", file=sys.stderr, flush=True)

with open(os.path.join(OUT, f"order_shuffled{TAG}.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
with open(os.path.join(OUT, f"order_shuffled{TAG}.json"), "w") as f:
    json.dump(dict(shuffles=K, shuffle_seed=SHUFFLE_SEED, steps=STEPS, reps=REPS,
                   seed=SEED, seconds=round(time.time() - t0, 1), rows=rows), f, indent=1)
write_tex([{k: str(v) for k, v in r.items()} for r in rows])
print(f"done in {time.time() - t0:.0f}s", file=sys.stderr)
