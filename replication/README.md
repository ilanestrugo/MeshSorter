# Replication-based MeshSorter simulator

This directory contains a standalone simulator of the MeshSorter and the scripts
that reproduce the tables of the manuscript. It shares no code with the programs
in the parent directory: `meshsorter.cpp` is a single self-contained C++17 file
that depends on nothing but the standard library.

The difference from the original `meshsim_cli` is the output analysis. Where the
older program takes one long run and forms batch means, this one runs `R`
statistically independent replications in parallel, one worker thread per
replication, discards a warm-up period from each, and forms an ordinary
Student `t` interval on the `R` replication averages. Because the replications
are independent by construction, the whole serial dependence of a run is
confined inside a single replication average, so no batch size has to be chosen
and no autocorrelation has to be tested. The protocol is described in full in
Supplement S1 of the manuscript.

## Contents

| File | Purpose |
| --- | --- |
| `meshsorter.cpp` | the simulator, one self-contained C++17 file |
| `exact_single_drop.py` | exact throughput of the unbuffered single-drop system, in rational arithmetic |
| `common.py` | helpers shared by the table scripts: run one configuration, cache the result, format LaTeX |
| `table1.py` | reproduces Table 1, panels (a) and (b) |
| `table2.py` | reproduces Table 2, panels (a), (b) and (c) |
| `run_all.sh` | builds the simulator and runs both table scripts |
| `results/` | generated output: the LaTeX bodies, the run logs, and the raw numbers |

## Build

```bash
c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp
```

Any C++17 compiler will do. The parent `CMakeLists.txt` also builds it, as the
target `meshsorter`:

```bash
cmake -S .. -B ../build -DCMAKE_BUILD_TYPE=Release
cmake --build ../build --parallel
```

Python 3.8 or newer is needed for the table scripts. They use only the standard
library.

## Quick start

```bash
./meshsorter -n 4 -m 4 --dual
```

runs the unbuffered dual-drop system with 4 primary belts and 4 feeder loops
under the manuscript's defaults, and prints

```
MeshSorter, dual-drop, n = 4 primary belts, m = 4 feeder loops
  feeder loop lengths    24 24 24 24
  crossings on a loop    forward 1 3 5 7
                         backward, feeder 1 21 19 17 15
  crossings on a belt    s^f 0 4 8 12   s^b 2 6 10 14
  buffers                none (largest capacity 0)
  design                 10 replications of 4000000 steps, warm-up 20000 each, seed 20260901

  throughput   3.02785   95% half-width 0.00021   [3.02764, 3.02806]
  replications  3.02797 3.02766 ...
```

`--json` prints the same result as one JSON object, which is what the table
scripts consume.

## The model

Time is discrete and every conveyor advances one slot per time step. There are
`n` primary belts and `m` circular feeder loops. Feeder `j` has `L_j` slots, and
slot 0 is its loading station. Each feeder crosses every primary belt once on
its forward run and, in the dual-drop design, once more on its backward run.

One time step has three phases, in this order.

1. **Loading.** For every feeder, if slot 0 is empty a new item is admitted and
   given a destination drawn uniformly from the `n` primary belts. A feeder
   whose slot 0 is occupied admits nothing; blocked items are not queued outside
   the system.
2. **Transfers.** Every crossing is examined. If the primary slot under the
   crossing is empty, a buffer at that crossing discharges into it first, and the
   feeder cell may then take the place freed in the buffer; if the buffer is
   empty the feeder cell transfers directly. If the primary slot is occupied, the
   feeder cell joins the buffer when there is room, and otherwise stays on the
   feeder for another revolution. Without buffers this is simply a transfer when
   the slot is free.
3. **Advancement.** Every conveyor moves one slot. An item that reaches the end
   of a primary belt leaves the system at once.

Throughput is reported as admissions per time step during the measurement
window. In steady state admissions and departures balance, so this is the
throughput of the system.

Crossings are examined feeder by feeder and, within a feeder, belt by belt, the
forward crossings before the backward ones. The order is immaterial under the
default geometry, where no two crossings share a primary slot, and is fixed only
so that unusual geometries remain reproducible.

## Geometry

The layout is built from two lengths measured on the feeder loop:

* `--dp D` (default 4): the loop slots consumed by each primary belt. A feeder
  crosses a belt twice, so consecutive belts sit `D/2` slots apart on the forward
  run and `D/2` apart on the backward run. `D` must be even.
* `--turn E` (default 4): the loop slots consumed by each end of the loop.

so the base loop length is

```
L = n * D + 2 * E                 (= 4n + 8 with the defaults)
```

and the crossings sit at

```
q^f_i = 1 + (i-1) * D/2                    (= 2i - 1)
q^b_i = L - E + 1 - (i-1) * D/2            (= 4n + 7 - 2i)
```

With the defaults, the last forward crossing and the first backward crossing of
a feeder are `D + E = 8` slots apart, and the last backward crossing is `E = 4`
slots from the loading station going forward around the loop.

Along a primary belt, consecutive feeders are `--df F` slots apart (default 4)
and the two crossings of one feeder are `--width W` slots apart (default 2):

```
s^f_j = (j-1) * F ,     s^b_j = s^f_j + W
```

`--width` may not exceed `--df`, or consecutive feeders would overlap.

`-L` sets the base loop length directly, overriding `n * D + 2 * E`; `--loops`
gives the `m` lengths explicitly.

### Staggered loop lengths

`--stagger` sets

```
L_1 = L_2 = L ,     L_j = L + s * (j - 2)   for j >= 3
```

with `s` from `--stagger-step` (default 1). The added slots have to go
somewhere on the loop, and where they go changes the distance a blocked item
travels between its forward and its backward crossing with the same belt.
`--extra` selects the placement:

* `turnaround` (default): the added slots extend the far end of the loop, so the
  backward crossings move away from the forward ones. This is the layout of
  Table 2 panel (c).
* `return`: the added slots extend the return run between the last backward
  crossing and the loading station, so the crossings keep the positions computed
  from the base length. This is the layout of Table 2 panel (b).

In a single-drop system the two coincide, since there are no backward crossings
to move.

### Buffers

`-b` sets the buffer capacity at the crossings, and accepts

* one value, used at every crossing;
* `m` values, one per feeder, used at all of that feeder's crossings;
* `2*m*n` values, listed feeder by feeder, the `n` forward crossings first and
  then the `n` backward ones.

For example `-b 3` puts three places at every crossing, and `-b 0,0,2,2` gives
the first two feeders of a four-feeder system no buffers and the last two three
places each at every crossing.

## Output analysis

`-R` replications differ only in their random number streams. Each discards its
first `-w` time steps and averages the throughput over the next `-T` steps. With
`Y_1..Y_R` the replication averages, the program reports

```
Ybar  +-  t_{R-1,0.975} * s / sqrt(R)
```

The default warm-up is the rule of Supplement S1,

```
W = max( 20000 , 10 * (c + 1) * L )
```

where `L` is the longest feeder loop and `c` the largest buffer capacity
anywhere in the configuration. The reasoning is that a feeder reaches each of
its crossings once per revolution, so a buffer of capacity `c` cannot fill from
empty in fewer than `c * L` steps; `(c+1) * L` is therefore the natural fill
time of the system, and the rule deletes ten of them, subject to a floor of
20,000 steps. A pilot study reported in Supplement S1 found the transient over
within 2,000 steps in every configuration in and beyond the range used in the
paper, and the estimates are insensitive to the choice: deleting nothing at all
moves them by less than 0.0006 items per step, and deleting 10,000 rather than
1,000,000 moves them by less than 0.0002.

Each `(replication, feeder)` pair draws from its own stream. The stream seed is
derived from `--seed` and the two indices, so the result does not depend on how
the replications are distributed over threads, and `--threads` changes only the
wall clock. Rerunning with the same `--seed` reproduces every number exactly.

## Reproducing the tables

```bash
./run_all.sh
```

or, one table at a time,

```bash
python3 table1.py
python3 table2.py
```

Each script carries its parameters at the top of the file, written out in full,
and prints the figures a caption needs: the largest half-width over the panel,
the warm-up actually used, and for Table 1 the pairwise-separation count. The
LaTeX bodies are written to `results/table1.tex` and `results/table2.tex`.

The raw numbers behind every published figure are kept under `results/raw`, and
are committed with the rest, so a reader can check a table without rerunning
anything:

* `<panel>_n<N>_m<M>.json`, one file per cell, holding the ten replication
  averages together with the geometry, the seed, the warm-up and the run length
  that produced them;
* `<panel>.csv`, one row per replication of that panel, with the same fields.

The panels are named `table1b_dual_identical`, `table2a_single_staggered`,
`table2b_dual_staggered_return` and `table2c_dual_staggered_turnaround`. The
mean and half-width of a cell can be recomputed from its ten values with any
tool at hand, which is the point of keeping them.

Results are also cached under `results/cache`, keyed on every argument that can
change the answer, so rerunning a table costs nothing and the two tables share
the panel they have in common. That directory is a cache and is not committed;
delete it, or set `MESHSORTER_NOCACHE=1`, to force a fresh run.

On two cores the full set takes about forty minutes. `--quick` runs the same
grids at low precision in a few minutes, which is enough to check that
everything is wired up.

`MESHSORTER_THREADS` sets the worker threads per cell, and defaults to the
number of cores. `MESHSORTER_BIN` points at the binary if it is not beside the
scripts.

## The exact single-drop values

Panel (a) of Table 1 is not simulated. The unbuffered single-drop system with a
common loop length decomposes into `L` independent copies of a chain on the
destinations the `m` feeders currently hold, and that chain is lumpable onto the
integer partitions of `m`: 30 states for `m = 9`, against `n^m = 387,420,489`
raw states. `exact_single_drop.py` builds the lumped kernel in exact rational
arithmetic and solves for the stationary distribution, so its output carries no
numerical error.

```bash
python3 exact_single_drop.py            # the 6x6 grid of Table 1(a)
python3 exact_single_drop.py 4 4        # one cell: 325/124 = 2.620967742
python3 exact_single_drop.py --verify   # check small cases against the raw chain
```

## Full option list

```
GEOMETRY
  -n, --belts N          number of primary belts                      (required)
  -m, --feeders M        number of feeder loops                       (required)
      --dual             dual-drop                                    (default)
      --single           single-drop: forward crossings only
      --dp D             loop slots per primary belt, even            (default 4)
      --turn E           loop slots per end of the loop               (default 4)
  -L, --loop LEN         base loop length                   (default N*D + 2*E)
      --loops a,b,...    explicit lengths for the M feeders
      --stagger          L, L, L+s, L+2s, ...
      --stagger-step s   the increment                                (default 1)
      --extra WHERE      turnaround or return                (default turnaround)
      --df F             feeder spacing along a primary belt          (default 4)
      --width W          forward-to-backward distance on a belt       (default 2)

BUFFERS
  -b, --buffers SPEC     1, M, or 2*M*N capacities                    (default 0)

EXPERIMENT
  -T, --steps T          measured steps per replication         (default 4000000)
  -R, --reps R           independent replications                     (default 10)
  -w, --warmup W         steps discarded per replication
                                       (default max(20000, 10*(c+1)*L))
      --seed S           base seed                             (default 20260901)
  -t, --threads K        worker threads                               (default R)

OUTPUT
      --json             machine-readable output
      --per-feeder       also report each feeder's loader utilization
  -h, --help
```

## Validation

The simulator reproduces the two independent references available for this
model:

* the exact single-drop values of `exact_single_drop.py`, for every cell of the
  6 by 6 grid, and
* the published dual-drop figures of Tables 1 and 2, which were produced by a
  separate implementation.

`validate.py` performs the first comparison:

```bash
python3 validate.py          # five representative cells
python3 validate.py --full   # the whole 6 by 6 grid
```

It prints, for each cell, the exact value, the simulated value, the half-width
and the discrepancy in standard errors.
