#!/usr/bin/env python3
"""Loop-length variants and geometry checks, 4 x 4 system at L = 104.

Produces the numbers quoted in four places of the manuscript that no table
carries:

  1. Loop-length variants (Section 4.5 of the main manuscript and Section S4.2
     of the Supplemental Online Material).  Unbuffered, both drop mechanisms, added
     slots at the turnaround:
         (L, L, L+1, L)       one added slot, enough for gcd = 1
         (L, L+1, L+2, L+3)   every loop after the first lengthened
     Quoted: 2.7088 and 3.0903 for the first, 2.7326 and 3.1073 for the second.

  2. Spacing invariance (Section 3.4).  Single-drop, L = 104, the primary belts
     2, 4, 10 and 16 slots apart along the loop, which put the drop points at
     (1,3,5,7), (1,5,9,13), (1,11,21,31) and (1,17,33,49).  Three variants:
     unbuffered, c = (0,1,2,7), and staggered lengths (104,104,105,106).
     Quoted: bit-for-bit identical output across the four spacings.

  3. Width invariance (Section 3.4).  Dual-drop, unbuffered, L = 104, four
     geometries that share d_i - gamma = (86,66,46,26):
         gamma = 10, 8, 6, 2  with  d = (96,76,56,36), (94,..), (92,..), (88,..)
     The simulator places the backward crossings at q^b_i = L - turn + 1 -
     (i-1)*dp/2, so d_i = L - turn - (i-1)*dp, and these four are turn = 8, 10,
     12, 16 with width = 10, 8, 6, 2.  Quoted: 3.03532 in every one, to five
     decimals and with the same half-width.

  4. Offset sweep (Section 6.2).  Dual-drop, unbuffered, L = 104, turn 12, the
     width gamma swept over 1, ..., 14, so d_1 - gamma runs from 91 down to 78.
     The width cannot exceed the feeder spacing df = 14, or consecutive feeders
     would overlap on a primary belt, and width 0 would put both crossings of a
     feeder at the same point of the belt, which is not a feeder at all.
     Quoted: the throughput ranges from 3.0282 to 3.0354, and the baseline
     width 10 gives the lowest value.

Every case uses the protocol of Section S2: thirty replications of 1,330,000
measured time steps, the warm-up rule (20,000 steps for these cases), base
seed 20260901.  The script prints each quoted figure next to the value it
reproduces and exits non-zero if a check fails.

Output
    results/geometry_checks.json   one record per case

    python3 geometry_checks.py            # 34 cases, about 6 minutes on 8 cores
    python3 geometry_checks.py --quick    # low precision, to check the wiring
"""
import sys, json
import common
import geometry as _geo

# ---------------------------------------------------------------- parameters
BELTS, FEEDERS = 4, 4
L       = _geo.loop_length(BELTS)                 # 104 slots
DP      = _geo.DP
TURN    = _geo.TURN
DF      = _geo.DF
WIDTH   = _geo.WIDTH
STEPS   = 1_330_000
REPS    = 30
WARMUP  = None                                    # None = the rule of Section S2
SEED    = 20260901

VARIANTS  = [(L, L, L + 1, L), (L, L + 1, L + 2, L + 3)]
SPACINGS  = [2, 4, 10, 16]
STAGGERED = (L, L, L + 1, L + 2)
ALLOC     = [0, 1, 2, 7]                          # c at B = 10, along each belt
WIDTHS    = [(8, 10), (10, 8), (12, 6), (16, 2)]  # (turn, width), d_i - gamma = 86,...
SWEEP     = range(1, DF + 1)                      # widths 1..14 at turn 12
# ---------------------------------------------------------------------------

if "--quick" in sys.argv:
    STEPS, REPS = 100_000, 8

common.require_binary()
rows, failures = [], []


def case(label, **kw):
    kw.setdefault("dp", DP); kw.setdefault("turn", TURN)
    kw.setdefault("df", DF); kw.setdefault("width", WIDTH)
    r = common.run(BELTS, FEEDERS, steps=STEPS, reps=REPS, warmup=WARMUP,
                   seed=SEED, **kw)
    rec = dict(case=label, mu=r["throughput"], hw=r["halfwidth"], sd=r["sd"],
               warmup=r["warmup"], loops=r["loops"], reps=r["reps"],
               steps=r["steps"], replications=r["replications"],
               args={k: (list(v) if isinstance(v, tuple) else v)
                     for k, v in kw.items()})
    rows.append(rec)
    print(f"  {label:48s} {r['throughput']:.6f} +-{r['halfwidth']:.6f}",
          file=sys.stderr, flush=True)
    return rec


def check(ok, what):
    print(f"  [{'ok' if ok else 'FAIL'}] {what}")
    if not ok:
        failures.append(what)


def quoted(value, printed, places):
    """Does `value`, rounded to `places`, read as the printed figure?"""
    return f"{value:.{places}f}" == printed


# 1. loop-length variants -------------------------------------------------
var = {}
for dual in (False, True):
    for loops in VARIANTS:
        var[(dual, loops)] = case(
            f"{'dual' if dual else 'single'} loops {loops}",
            dual=dual, loop=L, loops=loops)

# 2. spacing invariance, single drop ---------------------------------------
spc = {}
for name, extra in (("unbuffered", {}), ("c=(0,1,2,7)", {"buffers": ALLOC}),
                    ("staggered", {"loops": STAGGERED})):
    for s in SPACINGS:
        spc[(name, s)] = case(f"single spacing {s:2d} {name}", dual=False,
                              loop=L, dp=2 * s, **extra)

# 3. width invariance, dual drop -------------------------------------------
wid = {}
for turn, width in WIDTHS:
    wid[(turn, width)] = case(f"dual turn {turn:2d} width {width:2d}",
                              dual=True, loop=L, turn=turn, width=width)

# 4. offset sweep, dual drop -----------------------------------------------
swp = {}
for width in SWEEP:
    swp[width] = case(f"dual turn {TURN} width {width:2d} (d1-gamma {L - TURN - width})",
                      dual=True, loop=L, turn=TURN, width=width)

with open("results/geometry_checks.json", "w", newline="\n") as f:
    json.dump(rows, f, indent=1)
print("wrote results/geometry_checks.json", file=sys.stderr)

# ---------------------------------------------------------- the quoted figures
print("\nloop-length variants (Section 4.5, Section S4.2)")
for dual, loops, printed in ((False, VARIANTS[0], "2.7088"), (True, VARIANTS[0], "3.0903"),
                             (False, VARIANTS[1], "2.7326"), (True, VARIANTS[1], "3.1073")):
    mu = var[(dual, loops)]["mu"]
    check(quoted(mu, printed, 4),
          f"{'dual' if dual else 'single'} {loops}: {mu:.6f}, quoted {printed}")

print("\nspacing invariance, single drop (Section 3.4)")
for name in ("unbuffered", "c=(0,1,2,7)", "staggered"):
    reps = [spc[(name, s)]["replications"] for s in SPACINGS]
    check(all(x == reps[0] for x in reps),
          f"{name}: spacings {SPACINGS} give identical replications "
          f"({spc[(name, SPACINGS[0])]['mu']:.6f})")

print("\nwidth invariance, dual drop (Section 3.4)")
mus = [wid[k]["mu"] for k in WIDTHS]
hws = [wid[k]["hw"] for k in WIDTHS]
check(len({f"{x:.5f}" for x in mus}) == 1 and len({f"{x:.5f}" for x in hws}) == 1,
      f"four widths agree to five decimals in value and half-width: {mus[0]:.6f}")
check(quoted(mus[0], "3.03532", 5), f"{mus[0]:.6f}, quoted 3.03532")

print("\noffset sweep, dual drop (Section 6.2)")
lo = min(swp, key=lambda w: swp[w]["mu"])
hi = max(swp, key=lambda w: swp[w]["mu"])
check(quoted(swp[lo]["mu"], "3.0282", 4) and quoted(swp[hi]["mu"], "3.0354", 4),
      f"range {swp[lo]['mu']:.6f} (width {lo}) to {swp[hi]['mu']:.6f} "
      f"(width {hi}), quoted 3.0282 to 3.0354")
check(lo == WIDTH, f"the baseline width {WIDTH} gives the lowest value")

warm = sorted({r["warmup"] for r in rows})
print(f"\n{len(rows)} cases, {REPS} replications of {STEPS} steps, warm-up {warm}")
if failures and "--quick" not in sys.argv:
    sys.exit(f"{len(failures)} check(s) failed")
