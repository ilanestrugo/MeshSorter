#!/usr/bin/env python3
"""The resource comparison of Table 1 and Table S1: three conventional sorters
against one four-by-four MeshSorter.

Every cell of those two tables, and every figure quoted in Section S1 of the
Supplemental Online Material from them, is arithmetic on a handful of stated
dimensions.  This script states the dimensions once, derives each cell from them,
and compares the result with the number printed in the paper.  It also derives
the quantities of the section's prose: the linear formulas in the number of
destinations D, the break-even points, the ratios, and the aisle-travel figures.

    python3 layout_comparison.py          print the comparison; exit status 1 if any
                                          derived value differs from the printed one
    python3 layout_comparison.py -D 640   the same layouts at another number of
                                          destinations, without the comparison with the paper

Only one input is simulated, the unbuffered throughput of the four-by-four
MeshSorter with the staggered loop lengths, which is read from
results/loops_buffered.json (the row for a base loop of 104 slots).  If that file
is missing the value 3.105747368 that it holds is used and the report says so.
Results are written to results/layout_comparison.json (not with -D).

The inputs

    pitch                    0.5 m, one slot
    belt speed               2 m/s, so a time step is 0.25 s and a belt carries
                             14,400 parcels an hour
    bins                     one every 1 m along both sides of a belt, each with
                             its funnel 1 m wide, so a belt with its two rows of
                             bins is 3 m wide
    forklift aisle           2 m, between adjacent rows of bins and along the
                             outer sides
    loading area             conventional: 6 m at the head of the belts, the belts
                             occupying the first 3 m of it.  MeshSorter: the
                             feeder loops, each 5 m across with a 2 m gap between
                             consecutive loops, and a strip 3 m across the front
                             of the loops for the stations that load them
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")

# ---------------------------------------------------------------- the inputs
PITCH = 0.5                 # m per slot
SPEED = 2.0                 # m/s
STEP = PITCH / SPEED        # seconds per time step: 0.25
REQUIRED_PER_HOUR = 43_200  # parcels per hour
DESTINATIONS = 320
BIN_PITCH = 1.0             # m of belt per bin, on each side
BIN_WIDTH = 1.0             # m, each row of bins with its funnels
BELT_WIDTH = 1.0            # m, so that belt + two rows of bins = 3 m
AISLE = 2.0                 # m
CONV_LOADING = 6.0          # m, loading area at the head of the conventional belts
CONV_BELT_IN_LOADING = 3.0  # m of it that the belts occupy
LOOP_ACROSS = 5.0           # m, one feeder loop
LOOP_GAP = 2.0              # m, between consecutive loops
LOAD_STRIP = 3.0            # m across the front of the loops
MESH_BELTS, MESH_FEEDERS = 4, 4
BIN_CAPACITY = 500          # items per bin
FORKLIFT_KMH = 12.0         # average forklift speed
CONSOLIDATION_DRAWN = 578   # m^2, the dashed box of Figure 3 for the conventional layout
FALLBACK_CEILING = 3.105747368   # items per step, the value in results/loops_buffered.json
FALLBACK_LOOPS = [104, 104, 105, 106]


def parcels_per_hour(items_per_step):
    return items_per_step * 3600.0 / STEP


def belt_per_hour():
    return parcels_per_hour(1.0)            # a belt presents one slot per step


def simulated_ceiling():
    """The unbuffered dual-drop four-by-four ceiling with the staggered loops."""
    path = os.path.join(OUT, "loops_buffered.json")
    try:
        with open(path) as f:
            rows = json.load(f)
        for r in rows:
            if r["dual"] and r["B"] == 0 and r["stagger"] and r["L"] == 104:
                return r["mu"], r["loops"], path
    except (OSError, KeyError, ValueError):
        pass
    return FALLBACK_CEILING, FALLBACK_LOOPS, None


# ---------------------------------------------------------------- the layouts
def conventional(D, sorters=None):
    """k conventional sorters side by side, each reaching all D destinations."""
    k = sorters if sorters else -(-REQUIRED_PER_HOUR // int(belt_per_hour()))   # ceil
    bin_len = D / 2 * BIN_PITCH                       # bins on both sides
    belt_len = bin_len + CONV_BELT_IN_LOADING
    width = (k + 1) * AISLE + k * (BELT_WIDTH + 2 * BIN_WIDTH)
    length = bin_len + CONV_LOADING
    return dict(
        sorters=k, bins=k * D, drops=k * D, drops_bins=k * D, drops_feeder=0,
        belt=k * belt_len, belt_primary=k * belt_len, belt_loops=0.0,
        aisle=(k + 1) * bin_len,
        length=length, width=width, footprint=length * width,
        loading=CONV_LOADING * width, sorting=bin_len * width,
        ceiling=k * belt_per_hour(), bin_len=bin_len)


def mesh(D, ceiling, loops, n=MESH_BELTS, m=MESH_FEEDERS):
    """One n-by-m MeshSorter; the destinations are divided among the n belts."""
    bin_len = D / n / 2 * BIN_PITCH                   # bins of one belt, one side
    width_bins = (n + 1) * AISLE + n * (BELT_WIDTH + 2 * BIN_WIDTH)
    loading_len = m * LOOP_ACROSS + (m - 1) * LOOP_GAP
    width_feed = width_bins + LOAD_STRIP
    length = bin_len + loading_len
    loops_m = sum(loops) * PITCH
    return dict(
        bins=D, drops=D + 2 * n * m, drops_bins=D, drops_feeder=2 * n * m,
        belt=n * length + loops_m, belt_primary=n * length, belt_loops=loops_m,
        aisle=(n + 1) * length,
        length=length, width=width_bins, width_feed=width_feed,
        footprint=bin_len * width_bins + loading_len * width_feed,
        loading=loading_len * width_feed, sorting=bin_len * width_bins,
        ceiling=parcels_per_hour(ceiling), bin_len=bin_len, loading_len=loading_len)


# ---------------------------------------------------------------- the report
failures = []


def check(label, derived, printed, digits=0):
    """Compare a derived value with the one printed, at the printed precision."""
    ok = round(derived, digits) == round(printed, digits)
    shown = f"{derived:,.{max(digits, 0)}f}"
    print(f"  {'ok  ' if ok else 'DIFF'}  {label:<50} {shown:>14}   printed {printed:,.{max(digits, 0)}f}")
    if not ok:
        failures.append(label)
    return ok


def main():
    D_arg = None
    if "-D" in sys.argv:
        D_arg = int(sys.argv[sys.argv.index("-D") + 1])
    ceiling, loops, src = simulated_ceiling()
    D = D_arg or DESTINATIONS
    conv, msh = conventional(D), mesh(D, ceiling, loops)
    paper = D_arg is None

    print(f"Sorting {REQUIRED_PER_HOUR:,} parcels an hour to {D} destinations "
          f"(a belt carries {belt_per_hour():,.0f} an hour, a time step is {STEP} s)")
    if src:
        print(f"MeshSorter ceiling {ceiling:.6f} per time step, from {os.path.relpath(src, HERE)}"
              f" (loops {loops})")
    else:
        print(f"MeshSorter ceiling {ceiling:.6f} per time step: results/loops_buffered.json "
              "not found, the value it holds is used")
    print()

    print("Table S1 (and Table 1)")
    rows = [
        ("Total conveyor belt length (m)", "belt", 489, 473.5, 1),
        ("  Primary belt length (m)", "belt_primary", 489, 264, 1),
        ("  Feeder loops length (m)", "belt_loops", 0, 209.5, 1),
        ("Dropping devices, total", "drops", 960, 352, 0),
        ("  From primary belts to bins", "drops_bins", 960, 320, 0),
        ("  From feeder to primary belt", "drops_feeder", 0, 32, 0),
        ("Total aisle length (m)", "aisle", 640, 330, 0),
        ("Total footprint (m^2)", "footprint", 2822, 1530, 0),
        ("  Loading area (m^2)", "loading", 102, 650, 0),
        ("  Sorting conveyors, bins and aisles (m^2)", "sorting", 2720, 880, 0),
        ("Facility length (m)", "length", 166, 66, 0),
        ("Facility width (m)", "width", 17, 22, 0),
        ("Throughput ceiling, no buffers (per hour)", "ceiling", 43_200, 44_700, -2),
    ]
    if paper:
        check("Destination bins, conventional", conv["bins"], 960)
        check("Destination bins, MeshSorter", msh["bins"], 320)
        for label, key, pc, pm, dg in rows:
            check(label + ", conventional", conv[key], pc, dg)
            check(label + ", MeshSorter", msh[key], pm, dg)
        check("MeshSorter facility width at the feeders (m)", msh["width_feed"], 25)
        check("Parcels an hour at 3.106 per step", parcels_per_hour(3.106), 44_700, -2)
    else:
        for label, key, pc, pm, dg in rows:
            print(f"  {label:<50} {conv[key]:>14,.1f} {msh[key]:>14,.1f}")
    print()

    cf, mf = conv["footprint"], msh["footprint"]
    cb, mb = conv["belt"], msh["belt"]
    if paper:
        print("Quoted in the text of Section S1")
        check("Footprint saved (percent)", 100 * (1 - mf / cf), 46)
        check("Conveyor saved (percent, one decimal)", 100 * (cb - mb) / cb, 3.2, 1)
        requirement = REQUIRED_PER_HOUR / belt_per_hour()
        check("Requirement in parcels per time step", requirement, 3.0, 1)
        check("Staggered loops, total slots", sum(loops), 419)
        for lab, items, printed in (("3.791 per step, single-drop at B=10", 3.791, 54_600),
                                    ("3.840 per step, dual-drop at B=10", 3.840, 55_300),
                                    ("3.106 per step, unbuffered staggered", 3.106, 44_700)):
            check("Parcels an hour at " + lab, parcels_per_hour(items), printed, -2)
        print()

    # linear formulas in D
    c0 = conventional(0)
    m0 = mesh(0, ceiling, loops)
    c1 = conventional(1)
    m1 = mesh(1, ceiling, loops)

    def line(f0, f1):
        return f1 - f0, f0           # slope per destination, constant

    cb_s, cb_c = line(c0["belt"], c1["belt"])
    mb_s, mb_c = line(m0["belt"], m1["belt"])
    cf_s, cf_c = line(c0["footprint"], c1["footprint"])
    mf_s, mf_c = line(m0["footprint"], m1["footprint"])
    print("As functions of the number of destinations D")
    print(f"  belt        conventional {cb_s:.2f} D + {cb_c:.1f}      MeshSorter {mb_s:.2f} D + {mb_c:.1f}")
    print(f"  footprint   conventional {cf_s:.2f} D + {cf_c:.0f}      MeshSorter {mf_s:.2f} D + {mf_c:.0f}")
    be_belt = (mb_c - cb_c) / (cb_s - mb_s)
    be_foot = (mf_c - cf_c) / (cf_s - mf_s)
    print(f"  belt break-even D = {be_belt:.1f}      footprint break-even D = {be_foot:.1f}")
    if paper:
        check("belt, conventional slope", cb_s, 1.5, 2)
        check("belt, conventional constant", cb_c, 9, 1)
        check("belt, MeshSorter slope", mb_s, 0.5, 2)
        check("belt, MeshSorter constant", mb_c, 313.5, 1)
        check("belt break-even", be_belt, 304.5, 1)
        for d, printed in ((640, 335.5), (1280, 975.5)):
            check(f"conveyor saved at {d} destinations (m)",
                  conventional(d)["belt"] - mesh(d, ceiling, loops)["belt"], printed, 1)
        check("footprint, conventional slope", cf_s, 8.5, 2)
        check("footprint, conventional constant", cf_c, 102, 0)
        check("footprint, MeshSorter slope", mf_s, 2.75, 2)
        check("footprint, MeshSorter constant", mf_c, 650, 0)
        check("footprint break-even", be_foot, 95.3, 1)
        for d, printed in ((160, 1.3), (320, 1.8), (640, 2.3), (1280, 2.6)):
            check(f"footprint ratio at {d} destinations",
                  conventional(d)["footprint"] / mesh(d, ceiling, loops)["footprint"], printed, 1)
        check("ratio of the slopes", cf_s / mf_s, 3.1, 1)
        check("difference of the constants (m^2)", mf_c - cf_c, 548, 0)
        print(f"  (the dashed consolidation box of the figure draws {CONSOLIDATION_DRAWN} m^2, "
              f"above the {mf_c - cf_c:.0f} m^2 that closes the gap)")
        print()

    # aisle travel
    fills = REQUIRED_PER_HOUR / BIN_CAPACITY                  # bin replacements an hour
    km_c = fills * conv["bin_len"] / 1000
    km_m = fills * msh["bin_len"] / 1000
    print("Forklift travel")
    print(f"  {fills:.1f} bin replacements an hour, a round trip averaging one aisle length")
    if paper:
        check("conventional km an hour", km_c, 13.824, 3)
        check("MeshSorter km an hour", km_m, 3.456, 3)
        check("difference km an hour", km_c - km_m, 10.368, 3)
        check("forklift-hours per hour", (km_c - km_m) / FORKLIFT_KMH, 0.86, 2)
        c640, m640 = conventional(640), mesh(640, ceiling, loops)
        d640 = fills * (c640["bin_len"] - m640["bin_len"]) / 1000
        check("difference at 640 destinations, km an hour", d640, 20.7, 1)
        check("forklift-hours per hour at 640", d640 / FORKLIFT_KMH, 1.7, 1)
    else:
        print(f"  conventional {km_c:.3f} km an hour, MeshSorter {km_m:.3f}")

    # the one asymmetry of convention that a reader may ask about
    print()
    print("Note on the aisle rows: the conventional figure counts each aisle over the bins "
          f"({conv['bin_len']:.0f} m of the {conv['length']:.0f} m facility), "
          f"the MeshSorter's over the whole facility length ({msh['length']:.0f} m). "
          f"Counted over the bins only, the MeshSorter's would be "
          f"{(MESH_BELTS + 1) * msh['bin_len']:.0f} m; counted over the whole facility, "
          f"the conventional would be {(conv['sorters'] + 1) * conv['length']:.0f} m.")

    if paper:
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, "layout_comparison.json"), "w") as f:
            json.dump(dict(destinations=D, ceiling_per_step=ceiling, loops=loops,
                           conventional=conv, meshsorter=msh,
                           differences_from_paper=failures), f, indent=1)
    print()
    if paper:
        print("all derived values agree with the paper" if not failures
              else f"{len(failures)} derived value(s) differ from the paper: " + "; ".join(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
