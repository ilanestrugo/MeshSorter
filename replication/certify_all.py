#!/usr/bin/env python3
"""Table 5: certify the structured buffer-allocation class over the reported grid.

Runs certify_rep for both drop mechanisms, m = 3, 4 and 5 feeder loops, and
per-primary-belt budgets B = 0 to 10, and writes the LaTeX body of Table 5 along
with the raw JSON of every cell and, next to it, a CSV of every allocation the
cell evaluated.  Those CSV files are the input of approx_eval.py, which is what
Section 7 rests on, so the certification and the approximation-model evaluation
come from one run and cannot drift apart.

    python3 certify_all.py                 the full grid
    python3 certify_all.py --quick         a cheap pass, to check the wiring
    python3 certify_all.py --only dual     one mechanism
    python3 certify_all.py --feeders 3,4   a subset of the feeder counts
    python3 certify_all.py --budgets 0-5   a subset of the budgets

The last two select which cells to run; they do not change how a cell is
computed, so a cell produced under them is interchangeable with the same cell
from a full grid.

The full grid is 54,120 allocations, each piloted for R0 replications and
validated again on fresh streams, so expect several hours on a machine with
thirty usable threads.  Results are cached per cell, so an interrupted run
resumes where it stopped.

Every cell also writes a CSV of every allocation it evaluated, as raw material
for anyone who wants it.  Section 7 does not read these: it needs only the
structured class, which sweep_structured.py evaluates on its own in a couple of
minutes, so there is no reason to run this grid again for it.

CERTIFY_OUT redirects the results directory, which is useful for a trial run
that should not disturb a finished grid.
"""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BIN  = os.environ.get("CERTIFY_BIN", os.path.join(HERE, "certify_rep"))
OUT  = os.environ.get("CERTIFY_OUT", os.path.join(HERE, "results", "certify"))

# ---------------------------------------------------------------- parameters
BELTS    = 4
FEEDERS  = [3, 4, 5]
BUDGETS  = list(range(0, 11))
SPACING  = 4              # slots between consecutive primary belts along a loop
TURN     = 4              # slots in each end curve; loop length is 2*S*n + 2*E
STEPS    = 1_330_000
R0       = 10
KAPPA    = 3
EPS      = 0.001
ALPHA    = 0.01
SEED     = 20260901
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, R0, FEEDERS, BUDGETS = 100_000, 4, [3], [0, 1, 2]
ONLY = None
if "--only" in sys.argv:
    ONLY = sys.argv[sys.argv.index("--only") + 1]
SUBSET = "--feeders" in sys.argv or "--budgets" in sys.argv
if "--feeders" in sys.argv:
    FEEDERS = [int(x) for x in sys.argv[sys.argv.index("--feeders") + 1].split(",")]
if "--budgets" in sys.argv:
    spec = sys.argv[sys.argv.index("--budgets") + 1]
    BUDGETS = (list(range(int(spec.split("-")[0]), int(spec.split("-")[1]) + 1))
               if "-" in spec else [int(x) for x in spec.split(",")])

if not os.path.exists(BIN):
    sys.exit(f"{BIN} not found.  Build it first:\n"
             f"    c++ -O2 -std=c++17 -pthread -o certify_rep certify_rep.cpp")
os.makedirs(OUT, exist_ok=True)

def design():
    """Everything that changes the answer.  A cached cell produced under a
    different design is ignored rather than reused, so a --quick pass cannot
    contaminate a full run."""
    return dict(belts=BELTS, spacing=SPACING, turn=TURN, steps=STEPS, R0=R0,
                kappa=KAPPA, eps=EPS, alpha=ALPHA, seed=SEED)


def cell(dual, m, B):
    tag = f"{'dual' if dual else 'single'}_n{BELTS}_m{m}_B{B}"
    path = os.path.join(OUT, tag + ".json")
    dump = os.path.join(OUT, tag + "_allocs.csv")
    if os.path.exists(path):
        with open(path) as f:
            old = json.load(f)
        if old.get("design") == design():
            return old
        print(f"  {tag}: cached under a different design, rerunning",
              file=sys.stderr)
    cmd = [BIN, "-n", str(BELTS), "-m", str(m), "-B", str(B),
           "--dual" if dual else "--single",
           "--spacing", str(SPACING), "--turn", str(TURN),
           "-T", str(STEPS), "-R", str(R0), "--kappa", str(KAPPA),
           "--eps", str(EPS), "--alpha", str(ALPHA), "--seed", str(SEED),
           "--dump", dump, "--json"]
    if os.environ.get("MESHSORTER_THREADS"):
        cmd += ["-t", os.environ["MESHSORTER_THREADS"]]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"certify_rep failed on {tag}: {p.stderr.strip()}")
    r = json.loads(p.stdout)
    r["seconds"] = round(time.time() - t0, 1)
    r["design"] = design()
    with open(path, "w") as f:
        json.dump(r, f, indent=1)
    print(f"  {tag}: ref {r['reference']}  TH {r['throughput']:.3f}  "
          f"eps_hat {r['eps_hat']:.5f}  ({r['seconds']:.0f}s)",
          file=sys.stderr, flush=True)
    return r

res = {}
for m in FEEDERS:
    for B in BUDGETS:
        for dual in (False, True):
            if ONLY == "dual" and not dual: continue
            if ONLY == "single" and dual:  continue
            res[(dual, m, B)] = cell(dual, m, B)

# ------------------------------------------------------------- LaTeX output
lines = []
for m in FEEDERS:
    for B in BUDGETS:
        s = res.get((False, m, B)); d = res.get((True, m, B))
        if not s or not d: continue
        lines.append(f"{BELTS} & {m} & {B} & {s['reference']} & {s['throughput']:.3f} & "
                     f"{s['eps_hat']:.5f} & {d['reference']} & {d['throughput']:.3f} & "
                     f"{d['eps_hat']:.5f} \\\\")
    lines.append(r"\hline")
body = "\n".join(lines) + "\n"
if ONLY is None and not SUBSET:       # a partial pass cannot fill the table
    with open(os.path.join(HERE, "results", "table5.tex"), "w") as f:
        f.write(body)
    print(body)

vals = [r for r in res.values()]
if vals:
    zero = sum(1 for r in vals if r["eps_hat"] == 0.0)
    within = sum(1 for r in vals if r["eps_hat"] <= EPS)
    print(f"{len(vals)} comparisons: eps_hat = 0 in {zero}, "
          f"eps_hat <= {EPS} in {within}, largest {max(r['eps_hat'] for r in vals):.5f}")
    bad = [k for k, r in res.items() if not r["reference_confirmed"]]
    print(f"cells where the validation did not confirm the Phase 1 reference: {bad}")
    print(f"total wall time {sum(r.get('seconds', 0) for r in vals)/3600:.1f} hours")
