#!/usr/bin/env python3
"""Table 6: the order-derived destination-sequence robustness check.

The check asks one question.  Every other experiment in the paper draws each
item's destination uniformly and independently.  If those draws are replaced by
a single chronological sequence taken from a real order stream, do the buffering
trends survive?

For each per-primary-belt budget B the script simulates the same four-by-four
dual-drop system twice, once with uniform destinations and once with the
order-derived sequence, and reports the two throughputs and their difference.
Both columns use one protocol: the run length, warm-up rule, geometry and
replication count of Supplement S1, the same ones the certification and Section 7
use.  The earlier version of this table compared a long uniform run with a
shorter order-derived one, which is why it carried a caveat about horizons; this
one does not need it.

    python3 table6.py                 the table
    python3 table6.py --quick         a cheap pass, to check the wiring
    python3 table6.py --allocs FILE   read the allocations from FILE instead,
                                      one per line, lowest budget first

THE ALLOCATIONS
The buffer allocation at each budget is the one the certification selects, not a
free choice of this script.  Three sources are tried in turn:

  results/certify/               the certification grid itself, if it has been run
  results/table6_allocations.txt the same allocations, committed, so the short
                                 path reproduces the published table
  results/structured/            the best of the structured class, which
                                 sweep_structured.py produces in a couple of minutes

The first two agree by construction.  The third agrees at every budget except
B = 7, where it prefers 0,1,1,5 to the certification's 0,1,2,4; the two are
separated by 0.0004 items per time step, inside the half-width of either
estimate, and they agree to three decimals in both columns.  The source actually
used is printed and recorded in the output, so the table always says where its
allocations came from.

THE SEQUENCE
olist_orders_ForRun.csv, at the top of the repository, is the prepared order
stream: one row per usable order, in chronological order, with the primary belt
its destination is assigned to.  Supplement S4 describes how the raw Brazilian
e-commerce dataset was reduced to it.  This script extracts the belt column into
results/order_derived_sequence.txt, which is what the simulator reads.

The sequence is deterministic, so replications cannot differ in their draws.
They differ instead in phase: replication r begins reading the cyclic sequence at
its own position, fixed by its seed.  A replication of the reported length
consumes the 98,816 labels some fifty times over, so the spread across
replications measures how much the answer depends on where in the order stream
the day begins, and it is small.  That is a property of the check, not a
weakness of the estimate: this is one sequence, and the paper says so.
"""
import csv, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BIN  = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
OUT  = os.path.join(HERE, "results")
ORDERS = os.path.join(REPO, "olist_orders_ForRun.csv")
SEQ    = os.path.join(OUT, "order_derived_sequence.txt")

# ---------------------------------------------------------------- parameters
BELTS    = 4
FEEDERS  = 4
DUAL     = True
BUDGETS  = list(range(0, 11))
SPACING  = 4              # slots between consecutive primary belts along a loop
TURN     = 4              # slots in each end curve; loop length is 2*S*n + 2*E
STEPS    = 1_330_000
REPS     = 30
SEED     = 20260901
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, REPS, BUDGETS = 100_000, 6, [0, 1, 2]

if not os.path.exists(BIN):
    sys.exit(f"{BIN} not found.  Build it first:\n"
             f"    c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp")
os.makedirs(OUT, exist_ok=True)


def build_sequence():
    """Extract the belt column of the prepared order stream, in arrival order."""
    if not os.path.exists(ORDERS):
        sys.exit(f"{ORDERS} not found; it is the prepared order stream of Supplement S4")
    rows = []
    with open(ORDERS, newline="") as f:
        for d in csv.DictReader(f):
            rows.append((int(d["arrival_seq"]), int(d["belt"])))
    rows.sort()
    labels = [b for _, b in rows]
    bad = [b for b in labels if not 1 <= b <= BELTS]
    if bad:
        sys.exit(f"{ORDERS} carries belt labels outside 1..{BELTS}")
    share = [labels.count(i) / len(labels) for i in range(1, BELTS + 1)]
    with open(SEQ, "w") as f:
        f.write(f"# destination labels for {BELTS} primary belts, one per line,\n")
        f.write(f"# extracted from {os.path.basename(ORDERS)} in arrival order\n")
        f.write(f"# {len(labels)} orders; share per belt "
                + ", ".join(f"{s:.4f}" for s in share) + "\n")
        f.write("\n".join(str(b) for b in labels) + "\n")
    return labels, share


def allocations():
    """The certified allocation at each budget, and where it came from."""
    if "--allocs" in sys.argv:
        path = sys.argv[sys.argv.index("--allocs") + 1]
        with open(path) as f:
            got = [l.strip() for l in f if l.strip() and not l.startswith("#")]
        if len(got) != len(BUDGETS):
            sys.exit(f"{path} holds {len(got)} allocations, expected {len(BUDGETS)}")
        return dict(zip(BUDGETS, got)), os.path.basename(path)

    mech = "dual" if DUAL else "single"
    committed = os.path.join(OUT, "table6_allocations.txt")

    for folder, key in (("certify", "reference"),):
        d = os.path.join(OUT, folder)
        found, missing = {}, False
        for B in BUDGETS:
            if B == 0:
                found[B] = ",".join(["0"] * FEEDERS)
                continue
            p = os.path.join(d, f"{mech}_n{BELTS}_m{FEEDERS}_B{B}.json")
            if not os.path.exists(p):
                missing = True
                break
            v = json.load(open(p))[key]
            found[B] = v if isinstance(v, str) else ",".join(map(str, v))
        if not missing:
            return found, f"results/{folder}"

    if os.path.exists(committed):
        with open(committed) as f:
            got = [l.strip() for l in f if l.strip() and not l.startswith("#")]
        if len(got) >= max(BUDGETS) + 1:
            return {B: got[B] for B in BUDGETS}, "results/table6_allocations.txt"

    for folder, key in (("structured", "best"),):
        d = os.path.join(OUT, folder)
        found, missing = {}, False
        for B in BUDGETS:
            if B == 0:
                found[B] = ",".join(["0"] * FEEDERS)
                continue
            p = os.path.join(d, f"{mech}_n{BELTS}_m{FEEDERS}_B{B}.json")
            if not os.path.exists(p):
                missing = True
                break
            v = json.load(open(p))[key]
            found[B] = v if isinstance(v, str) else ",".join(map(str, v))
        if not missing:
            return found, f"results/{folder}"

    sys.exit("no certified allocations found.  Run either\n"
             "    python3 certify_all.py --only dual --feeders 4\n"
             "or, far more cheaply,\n"
             "    python3 sweep_structured.py")


def buffer_spec(alloc):
    """The 2*m*n form the simulator wants: nothing forward, alloc backward."""
    c = [int(x) for x in alloc.split(",")]
    if len(c) != FEEDERS:
        sys.exit(f"allocation {alloc!r} does not have {FEEDERS} entries")
    out = []
    for cj in c:
        out += [0] * BELTS + ([cj] * BELTS if DUAL else [0] * BELTS)
    if not DUAL:                       # single-drop: capacity at the forward points
        out = []
        for cj in c:
            out += [cj] * BELTS + [0] * BELTS
    return ",".join(map(str, out))


def simulate(alloc, sequence):
    cmd = [BIN, "-n", str(BELTS), "-m", str(FEEDERS),
           "--dual" if DUAL else "--single",
           "--spacing", str(SPACING), "--turn", str(TURN),
           "-b", buffer_spec(alloc),
           "-T", str(STEPS), "-R", str(REPS), "--seed", str(SEED), "--json"]
    if sequence:
        cmd += ["--sequence", SEQ]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"meshsorter failed on {alloc} "
                 f"({'sequence' if sequence else 'uniform'}): {p.stderr.strip()}")
    return json.loads(p.stdout)


labels, share = build_sequence()
allocs, source = allocations()
print(f"sequence   : {len(labels)} labels, share per belt "
      + ", ".join(f"{s:.4f}" for s in share), file=sys.stderr)
print(f"allocations: {source}", file=sys.stderr)

rows = []
t0 = time.time()
for B in BUDGETS:
    a = allocs[B]
    u = simulate(a, False)
    o = simulate(a, True)
    rows.append(dict(n=BELTS, m=FEEDERS, B=B, alloc=a,
                     th_uniform=u["throughput"], hw_uniform=u["halfwidth"],
                     th_order=o["throughput"], hw_order=o["halfwidth"],
                     gap=u["throughput"] - o["throughput"], loop=u["loops"][0]))
    print(f"  B={B:2d}  {a:12s}  uniform {u['throughput']:.4f}  "
          f"order-derived {o['throughput']:.4f}  gap {rows[-1]['gap']:+.4f}",
          file=sys.stderr, flush=True)

with open(os.path.join(OUT, "table6.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

with open(os.path.join(OUT, "table6.tex"), "w") as f:
    for r in rows:
        g = 0.0 if abs(r["gap"]) < 5e-4 else r["gap"]   # never print a signed zero
        f.write(f"{r['n']} & {r['m']} & {r['B']} & {r['alloc']} & "
                f"{r['th_uniform']:.3f} & {r['th_order']:.3f} & "
                f"{g:.3f} \\\\\n")
    f.write("\\hline\n")

with open(os.path.join(OUT, "table6.log"), "w") as f:
    json.dump(dict(belts=BELTS, feeders=FEEDERS, dual=DUAL, spacing=SPACING,
                   turn=TURN, loop=rows[0]["loop"], steps=STEPS, reps=REPS,
                   seed=SEED, allocation_source=source,
                   sequence_labels=len(labels), sequence_share=share,
                   largest_gap=max(abs(r["gap"]) for r in rows),
                   seconds=round(time.time() - t0, 1), rows=rows), f, indent=1)

print(f"\nlargest gap {max(abs(r['gap']) for r in rows):.4f} items per time step, "
      f"in {time.time() - t0:.0f}s", file=sys.stderr)
print("wrote results/table6.tex, results/table6.csv, results/table6.log",
      file=sys.stderr)
