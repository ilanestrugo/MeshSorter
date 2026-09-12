#!/usr/bin/env python3
"""Section 7: evaluate the structured class by simulation, for the approximation
model to be compared against.

The design rules of Section 5.1 define a small class of buffer allocations, and
the certification reported in that section establishes that the class contains
an allocation optimal up to the indifference zone.  A designer therefore
searches the class rather than the whole allocation space, and Section 7 asks
how well the approximation model of Section 6 carries out that search.  Only the
class has to be simulated, which is a few hundred allocations rather than the
54,120 the certification enumerates.

    python3 sweep_structured.py                 the four-belt grid of Table 7
    python3 sweep_structured.py --grid large    the fifteen-belt grid of S7
    python3 sweep_structured.py --quick         a cheap pass, to check the wiring

    main    four primary belts, 3, 4 and 5 feeder loops, budgets 1 to 10:
            388 allocations over both drop mechanisms, a couple of minutes
    large   fifteen primary belts, 10 and 15 feeder loops, budget 15:
            664 allocations, about an hour on thirty threads

The two grids answer different halves of the question.  The four-belt grid is
the one where Section 5.1 has shown, exhaustively, that the class contains a
near-optimal allocation, so an accurate ranking of the class there is an
accurate ranking of the design problem.  The fifteen-belt grid is the scale the
model exists for, where the full space holds about 1.1 billion allocations and
no exhaustive check is possible.

The run length, warm-up rule, replication count, seeds and stream families are
those of the certification, so every experiment in the paper shares one
protocol.  Results are cached per cell, so an interrupted run resumes.

Feed the output to approx_eval.py:

    python3 approx_eval.py
    python3 approx_eval.py --dir results/large --tag large --scatter 15
"""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN  = os.environ.get("SWEEP_BIN", os.path.join(HERE, "sweep_rep"))

# ---------------------------------------------------------------- parameters
GRIDS = {
    "main":  dict(belts=4,  feeders=[3, 4, 5],  budgets=list(range(1, 11)),
                  out="structured"),
    "large": dict(belts=15, feeders=[10, 15],   budgets=[15], out="large"),
}
SPACING  = 4              # slots between consecutive primary belts along a loop
TURN     = 4              # slots in each end curve; loop length is 2*S*n + 2*E
STEPS    = 1_330_000
REPS     = 15             # per pass; two disjoint stream families pool to the
                          # R = 30 replications of Section S1.4
SEED     = 20260901
# ---------------------------------------------------------------------------

name = "main"
if "--grid" in sys.argv:
    name = sys.argv[sys.argv.index("--grid") + 1]
if name not in GRIDS:
    sys.exit(f"unknown grid {name!r}; choose from {', '.join(GRIDS)}")
G = dict(GRIDS[name])
BELTS, FEEDERS, BUDGETS = G["belts"], G["feeders"], G["budgets"]
OUT = os.environ.get("SWEEP_OUT", os.path.join(HERE, "results", G["out"]))

if "--quick" in sys.argv:
    STEPS, REPS, FEEDERS, BUDGETS = 100_000, 4, FEEDERS[:1], BUDGETS[:2]
if "--feeders" in sys.argv:
    FEEDERS = [int(x) for x in sys.argv[sys.argv.index("--feeders") + 1].split(",")]
if "--budgets" in sys.argv:
    spec = sys.argv[sys.argv.index("--budgets") + 1]
    BUDGETS = (list(range(int(spec.split("-")[0]), int(spec.split("-")[1]) + 1))
               if "-" in spec else [int(x) for x in spec.split(",")])

if not os.path.exists(BIN):
    sys.exit(f"{BIN} not found.  Build it first:\n"
             f"    c++ -O2 -std=c++17 -pthread -o sweep_rep sweep_rep.cpp")
os.makedirs(OUT, exist_ok=True)


def design():
    """Everything that changes the answer.  A cached cell produced under a
    different design is ignored rather than reused."""
    return dict(belts=BELTS, spacing=SPACING, turn=TURN, steps=STEPS,
                reps=REPS, seed=SEED)


def cell(dual, m, B):
    tag = f"{'dual' if dual else 'single'}_n{BELTS}_m{m}_B{B}"
    path = os.path.join(OUT, tag + ".json")
    dump = os.path.join(OUT, tag + "_allocs.csv")
    if os.path.exists(path) and os.path.exists(dump):
        with open(path) as f:
            old = json.load(f)
        if old.get("design") == design():
            return old
        print(f"  {tag}: cached under a different design, rerunning",
              file=sys.stderr)
    cmd = [BIN, "-n", str(BELTS), "-m", str(m), "-B", str(B),
           "--dual" if dual else "--single",
           "--spacing", str(SPACING), "--turn", str(TURN),
           "-T", str(STEPS), "-R", str(REPS), "--seed", str(SEED),
           "--dump", dump, "--json"]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"sweep_rep failed on {tag}: {p.stderr.strip()}")
    r = json.loads(p.stdout)
    r["seconds"] = round(time.time() - t0, 1)
    r["design"] = design()
    with open(path, "w") as f:
        json.dump(r, f, indent=1)
    print(f"  {tag}: {r['allocations']} allocations, best {r['best']}, "
          f"TH {r['throughput']:.3f}  ({r['seconds']:.0f}s)",
          file=sys.stderr, flush=True)
    return r


res = {}
for m in FEEDERS:
    for B in BUDGETS:
        for dual in (False, True):
            res[(dual, m, B)] = cell(dual, m, B)

n = sum(r["allocations"] for r in res.values())
secs = sum(r.get("seconds", 0) for r in res.values())
print(f"\ngrid {name}: {len(res)} cells, {n} allocations, {secs / 60:.1f} minutes")
print("now run:  python3 approx_eval.py"
      + (f" --dir results/{G['out']} --tag {name} --scatter {BUDGETS[0]}"
         if name != "main" else ""))
