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
| `meshsorter_core.hpp` | the model, the geometry and one replication |
| `meshsorter.cpp` | the command-line simulator |
| `certify_rep.cpp` | replication-based certification of a structured buffer-allocation class (Supplement S3) |
| `certify_all.py` | reproduces Table 4 and Table S3, the certification grid, and writes the per-allocation dumps |
| `approx_model.py` | the analytical approximation model of Section 6, with its self-checks |
| `approx_eval.py` | reproduces Table 7 and Figure 7, and Tables S6, S7 and S8, the accuracy of the approximation against the certified grid |
| `sweep_rep.cpp` | evaluates a stated set of allocations by replications, for systems too large to enumerate |
| `sweep_structured.py` | simulates the structured class: the four-belt grid of Table 7, and the fifteen-belt grid of Section S8 |
| `sec7_numbers.py` | every number that appears in the prose of Section 7, taken from the data |
| `exact_single_drop.py` | exact throughput of the unbuffered single-drop system, in rational arithmetic |
| `common.py` | helpers shared by the table scripts: run one configuration, cache the result, format LaTeX |
| `table1.py` | reproduces Table 1, panels (a) and (b) |
| `table2.py` | reproduces Table 2, panels (a), (b) and (c) |
| `loops_buffered.py` | reproduces Table 5, how loop length and staggering act on a buffered system, at three per-primary-belt budgets |
| `tableS4S5.py` | reproduces Tables S4 and S5, throughput against load balancing over the structured class |
| `frontier.py` | reproduces Figure 6, the efficient frontier of throughput against load balancing |
| `table6.py` | reproduces Table 6, the order-derived destination-sequence robustness check |
| `feeder_returns.py` | Figure 4, the diminishing return of additional feeder loops |
| `run_all.sh` | builds the simulator and runs all four scripts |
| `results/` | generated output: the LaTeX bodies, the run logs, and the raw numbers |

## Build

```bash
c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp
c++ -O2 -std=c++17 -pthread -o certify_rep certify_rep.cpp
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
under the manuscript's defaults, and prints (the worker count shown is what a
machine with 32 hardware threads chooses)

```
MeshSorter, dual-drop, n = 4 primary belts, m = 4 feeder loops
  feeder loop lengths    24 24 24 24
  crossings on a loop    forward 1 3 5 7
                         backward, feeder 1 21 19 17 15
  crossings on a belt    s^f 0 4 8 12   s^b 2 6 10 14
  buffers                none (largest capacity 0)
  design                 30 replications of 1330000 steps, warm-up 20000 each, seed 20260901
                         30 worker threads

  throughput   3.02762   95% half-width 0.00020   [3.02742, 3.02782]
  replications  3.02699 3.02742 3.02749 3.02733 3.02809 ...
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

The defaults are 30 replications of 1,330,000 steps. Only the product of the two
fixes the width of the interval, so this is the same simulation effort as 10
replications of 4,000,000 steps, but the variance estimate carries 29 degrees of
freedom rather than 9. That narrows the interval a little, since `t` falls from
2.262 to 2.045, and, more usefully, makes the half-width itself a stable
quantity instead of one that varies by a quarter from configuration to
configuration.

The replications are handed out from a shared counter, so a worker that finishes
early takes the next one rather than waiting for the rest of its wave.
`--threads` defaults to two fewer than the hardware threads the machine reports,
capped at the number of replications: 30 workers on a machine with 32 hardware
threads, which runs the default 30 replications in a single pass and leaves two
threads for everything else. Replication `r` always draws the same stream, so
the number of workers changes the wall clock and nothing else.

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

Each `(replication, feeder)` pair draws from its own stream, seeded from
`--seed` and the two indices, so rerunning with the same `--seed` reproduces
every number exactly.

## Reproducing the tables

```bash
./run_all.sh
```

or, one table at a time,

```bash
python3 table1.py
python3 table2.py
python3 loops_buffered.py
python3 tableS4S5.py
python3 frontier.py
python3 table6.py
python3 feeder_returns.py
```

Each script carries its parameters at the top of the file, written out in full,
and prints the figures a caption needs: the largest half-width over the panel,
the warm-up actually used, and for Table 1 the pairwise-separation count.

Three kinds of output are written:

* `results/table1.tex`, `results/table2.tex`, `results/loops_buffered.tex`, the
  LaTeX bodies of the panels, ready to paste into the manuscript, and
  `results/loops_buffered_table.tex`, the whole of Table 5 with its caption, the
  allocations named and the largest half-width filled in;
* `results/table1.csv`, `results/table2.csv`, one row per cell with the
  estimate, its half-width and standard deviation, the gain over the reference
  panel where there is one, and the design and geometry that produced it;
* `results/raw/`, described below.

The CSVs are plain comma-separated text rather than a spreadsheet format, so the
scripts need nothing beyond the Python standard library; every spreadsheet
program opens them directly.

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

The full set is 240 configurations at about 4 x 10^7 simulated steps each. On
two cores that took a few hours; on a machine that can run all 30 replications at
once it is a small fraction of that. `--quick` runs the same grids at low
precision in a few minutes, which is enough to check that everything is wired up.

Every script sets the run length and the replication count explicitly at the top
of the file, at the 30 replications of 1,330,000 steps that Supplement S1
describes, so that none of them can inherit a different design from `common.py`
by accident. Every raw file records the `steps` and `replications` that produced
it, so the provenance of a number is never in doubt.

`MESHSORTER_THREADS` sets the worker threads per cell; left unset, the simulator
applies its own rule, two fewer than the machine's hardware threads, capped at
the number of replications. `MESHSORTER_BIN` points at the binary if it is not
beside the scripts.

## Throughput against load balancing

Tables S4 and S5 ask what a designer gives up by balancing the loaders instead
of maximizing throughput. At each per-primary-belt budget two allocations are
reported: `c`, the one with the largest minimum loader utilization, and `c'`,
the one with the largest throughput, together with both objectives for each and
the gap between them.

```bash
python3 tableS4S5.py             # about ten minutes
python3 tableS4S5.py --quick     # a cheap pass, to check the wiring
python3 tableS4S5.py --only dual # one drop mechanism
```

**The search is confined to the structured class.** That is a deliberate
restriction rather than an approximation of a wider search. Section 5.1 of the
manuscript certifies the class as containing a near-optimal design, so it is
what a designer following the design rules would consider; but the rules were
derived for throughput, not for balance, and the best-balanced allocation in the
whole space may well lie outside the class. What these tables report is the best
balance available to someone who follows the rules.

**Both design rules give the same set here.** With nothing on feeder 1 and
capacities nondecreasing along the feeders, the class is the partitions of the
budget into at most `m - 1` parts, whichever drop mechanism is in force. Under
the dual-drop mechanism those are the capacities at the backward drop points,
the forward ones being zero. The enumerator produces 1, 2, 3, 4, 5, 7, 8, 10,
12 and 14 allocations at budgets 1 to 10, matching the class sizes behind
Table 7.

**Loader utilization** comes from `meshsorter --per-feeder`, which reports the
fraction of rounds in which each feeder admits an item. The minimum is over all
feeders; feeder 1 is never blocked, so it is always the downstream end that
binds.

**Ties.** Where one allocation is best on both objectives, the two columns name
it once and both gaps are reported as zero rather than as a difference between
two spellings of the same design. Ties within an objective break toward the
allocation that is better on the other.

**One protocol.** Every allocation is evaluated under the design of
Supplement S1, the same geometry, run length, warm-up rule and thirty
replications the certification and Section 7 use, so these tables are directly
comparable with Table 4. The earlier version of Tables S4 and S5 was produced at
the shorter loop length and is not comparable with it: its throughputs run up to
0.034 items per time step low, and the gap widens with the budget.

**Output.** `results/tableS4.tex` and `results/tableS5.tex` are the
nine-column bodies, and `results/tableS4S5.csv` carries the per-feeder loader
utilizations behind every reported row. Results are cached per cell under
`results/loadbalance/`, so an interrupted run resumes and a rerun is free.

## The order-derived destination sequence

Every other experiment in this package draws each item's destination uniformly
and independently. Table 6 asks what happens when those draws are replaced by a
single chronological sequence taken from a real order stream, and whether the
buffering trends survive it.

```bash
python3 table6.py            # about two minutes
python3 table6.py --quick    # a cheap pass, to check the wiring
```

**The sequence.** `olist_orders_ForRun.csv`, at the top of the repository, is the
prepared order stream: one row per usable order, in chronological order, with the
primary belt its destination is assigned to. Supplement S5 describes how the raw
Brazilian e-commerce dataset was reduced to it, and the file is committed so the
reduction does not have to be repeated. `table6.py` extracts its `belt` column
into `results/order_derived_sequence.txt`, one label per line, which is what the
simulator reads. There are 98,816 labels, close to balanced across the four
belts, between 0.2494 and 0.2506 of the stream each.

**How the simulator consumes it.** `meshsorter --sequence FILE` replaces the
uniform draw at the loading step. The feeders consume the sequence in index
order, so when several admit an item in the same time step they take consecutive
labels, and the sequence repeats cyclically. That is the rule stated in
Supplement S5.

**What a replication means here.** The sequence is deterministic, so replications
cannot differ in their draws. They differ in phase: each begins reading the cycle
at its own position, fixed by its seed. A replication of the reported length
consumes the 98,816 labels some fifty times over, so the spread across
replications measures how much the answer depends on where in the order stream
the day begins. It runs about a third of the spread of the corresponding uniform
run. That is a property of the check rather than a weakness of the estimate:
this is one sequence, and the manuscript says so.

**One protocol for both columns.** The two columns of Table 6 use the same run
length, warm-up rule, geometry, seed and replication count, the ones of
Supplement S1 that the certification and Section 7 also use. An earlier version
of the table compared a long uniform run with a shorter order-derived one and
carried a caveat about horizons; this one does not need it. The uniform column
now agrees with Table 4 at `B = 10`, which is a useful check that the two
experiments really do share a protocol.

**The allocations.** The allocation at each budget is the one the certification
selects, not a free choice of the script. Three sources are tried in turn:

| Source | What it is |
| --- | --- |
| `results/certify/` | the certification grid itself, if it has been run |
| `results/table6_allocations.txt` | the same allocations, committed, so the short path reproduces the published table |
| `results/structured/` | the best of the structured class, two minutes of work |

The first two agree by construction: the committed file is the cumulative sum of
the addition order Table 4 prints. The third agrees at every budget but `B = 7`,
where it prefers `0,1,1,5` to the certification's `0,1,2,4`. Those two are
separated by 0.0004 items per time step, inside the half-width of either
estimate, and they agree to three decimals in both columns of the table, so the
choice does not move a printed figure. Whichever source was used is printed and
recorded in `results/table6.log`, so the table always says where its allocations
came from. Passing `--allocs FILE` overrides all three, one allocation per line,
lowest budget first.

**Output.** `results/table6.tex` is the LaTeX body, `results/table6.csv` has one
row per budget with both estimates and their half-widths, and
`results/table6.log` records the geometry, the protocol, the allocation source
and the largest gap.

## Certifying a buffer-allocation class

`certify_rep` answers a different question from the simulator: given a
per-primary-belt budget, does the structured class of allocations defined by the
manuscript's design rules contain the best allocation? It enumerates every
feasible allocation and splits effort three ways. A competitor is only being
screened, so it gets ten pilot replications. A member of the structured class is a
candidate for the reference, so it gets thirty. The reference itself gets at least
five hundred, because it enters every comparison in the cell and so is the one
place where precision is shared rather than spent once. Each competitor is then
sized against that reference, and the whole comparison is revalidated on a
disjoint family of random number streams. What it reports is the smallest relative
gap the data support for any allocation outside the class.

A competitor whose throughput sits at the very edge of the indifference zone needs
a reference of size proportional to the inverse square of its margin, which
diverges. No ceiling can accommodate it, so the program leaves such a competitor
at the pilot count and lets it enlarge the reported gap, which is the honest
outcome rather than a failure.

```bash
./certify_rep -n 4 -m 4 -B 10 --dual        # one cell
./certify_rep --help                        # the options
python3 certify_all.py                      # the grid of Table 4
python3 certify_all.py --verbose            # per-phase progress from each cell
python3 certify_all.py --feeders 3          # one feeder count, to split across machines
```

Three replication counts, each sized by how much work it does:

| flag | default | what it sets |
|---|---|---|
| `-R`, `--reps` | 10 | pilot replications for a competitor |
| `--reps-structured` | 30 | pilot replications for a structured allocation |
| `--reps-reference` | 500 | floor on the reference, which enters every comparison |

`--plan-only` runs the pilot and the planner, prints what validation would cost,
and exits. Use it to price a cell before committing to it. `--verbose` reports
each phase as it runs, at a granularity that adapts to the size of the cell.

The earlier program `buffer_certify` in the parent directory answers the same
question by batch means along a single long run. `certify_rep` replaces the batch
by the replication, which is what removes the need to argue that successive
observations are uncorrelated; the procedure is otherwise the same three phases,
and is documented in Supplement S1 and S2.

Work is spread over the allocations rather than over the replications of any one
of them, which is what keeps every thread busy: there are tens of thousands of
allocations and only ten replications each, so parallelizing the inner loop would
cap the program at ten active threads whatever the machine offers. The reference
is the exception. It carries hundreds of replications rather than ten, so it is
run through the pooled form instead; left in the sweep it would execute on a
single thread and become the critical path of every cell.

The grid is 54,114 allocations across both mechanisms, three feeder counts and
ten budgets, B = 1 to 10. A budget of zero admits a single allocation and leaves
no competitors, so there is nothing to certify and it is not run. each piloted and then validated, so it is an overnight run on a
machine with thirty usable threads. Cells are cached under `results/certify`,
keyed on the whole design, so an interrupted run resumes and a `--quick` pass
cannot contaminate a full one.

Every cell also writes `results/certify/<cell>_allocs.csv`, one row per
allocation, holding the capacities, whether the allocation belongs to the
structured class, and the mean and standard deviation of both passes over it.
Phases 1 and 3 draw from disjoint streams, so each row carries two independent
estimates of the same throughput and their difference measures the simulation
noise at no extra cost. These files are what the approximation-model evaluation
below consumes, which is why that evaluation and the certification cannot drift
apart: they are the same run. A cell is reused from cache only when both its
JSON and its CSV are present, so a grid certified before the dump existed is
recomputed once.

`CERTIFY_OUT` redirects the results directory, which is convenient for a trial
run that should not disturb a finished grid.

## Evaluating the approximation model

Section 6 replaces the simulation by a recursion along the feeders, in which
every buffer is a Geo/Geo/1/c queue whose birth and death probabilities are read
off the belt utilization upstream of it. `approx_model.py` is that recursion,
transcribed equation by equation, and Section 7 asks how much is lost by using
it in place of the simulation.

```bash
python3 approx_model.py        # the self-checks
python3 sweep_structured.py    # simulate the structured class, a couple of minutes
python3 approx_eval.py         # the evaluation
```

The comparison is confined to the structured class, the allocations the design
rules of Section 5.1 admit. That is the set a designer searches, since the
certification of that section establishes that it contains an allocation optimal
up to the indifference zone, and it is what makes the question tractable at
scale. Only the class has to be simulated, which is 388 allocations over both drop
mechanisms rather than the 54,120 the certification enumerates, so
`sweep_structured.py` produces it in a couple of minutes and the certification
does not have to be run again. `approx_eval.py` also reads the certification's
own dumps if you have them, with `--dir results/certify`, and `--all` then
includes the allocations outside the class.

For each class, that is each combination of drop mechanism, feeder count and
budget, it reports

* the relative gap of every allocation, by its mean and its maximum;
* where the allocation the model prefers stands in the simulated ordering, and
  the throughput given up by building it;
* the Spearman correlation of the two orderings, reported only where a class
  holds at least ten allocations, since a rank correlation over fewer says
  nothing.

At four primary belts the structured class is small, between one and 23
allocations per class, so the ranking evidence comes from the fifteen-belt grid
below.

Output lands in `results/`: `approx_eval.csv` with one row per class,
`approx_eval.tex` with the LaTeX body of the table, `approx_eval.json` with
everything, and `approx_scatter_<mechanism>_B<budget>.dat` with the data behind
the scatter panels. The scatter files are thinned by a fixed stride so that the
figure stays a reasonable size, but the extreme gaps and the two selected
allocations of every panel are always kept, so the picture cannot hide its own
worst case.

The unbuffered system is left out of the tables. It holds a single allocation,
so nothing in it can be ranked or selected, and Section 4 solves it exactly, so
the approximation is never used there. Its numbers stay in `approx_eval.csv`.

`--scatter` chooses the budgets whose allocations are written out for the
figure; the default is 5 and 10, and the fifteen-belt grid uses 15.

`sec7_numbers.py` prints every number that appears in the prose of Section 7,
and `sec7_numbers.py --map` prints them as a JSON mapping from the placeholders
the manuscript uses, so a rerun of the grid regenerates the sentences instead of
inviting a hand edit.

### A caution about the build products in a synced folder

`meshsorter`, `certify_rep` and `sweep_rep` are git-ignored, but if the clone
lives in a folder that syncs between machines, the binaries travel anyway and a
macOS build lands on the Linux box as a Mach-O file that will not execute. Build
on each machine before running, or keep the binaries outside the synced tree and
point `CERTIFY_BIN` and `SWEEP_BIN` at them.

### Systems too large to enumerate

The comparison above is exhaustive, which confines it to a four-belt system:
fifteen feeder loops with a budget of twenty admit about 1.1 billion
allocations. `sweep_rep` evaluates a stated set of allocations instead, by
default the structured class, which is what a designer following the design
rules would search. It writes the same CSV as `certify_rep --dump`, so the same
analysis reads either.

```bash
./sweep_rep -n 15 -m 15 -B 15 --dual --dump out.csv
python3 sweep_structured.py --grid large
python3 approx_eval.py --dir results/large --tag large --scatter 15
```

The reported grid is fifteen primary belts, ten and fifteen feeder loops, a
budget of fifteen, and both drop mechanisms: 664 allocations, about an hour on
thirty threads. `--wide` adds budgets of ten and twenty, about three hours. The
run length, replication count, seeds and stream families are those of the
certification, so the two experiments share one protocol and one warm-up rule.

`approx_model.py` implements all four buffer placements, including the two the
design rules reject, although nothing reported in the paper needs them: the
structured class places buffers at backward drop points only. It verifies itself
against three things: the closed form
`1 - (1 - 1/n)^m` for the unbuffered single-drop system; the identity that the
loader utilizations sum to `n * u_m`, which is flow conservation and must hold
exactly in every one of the four buffer placements; and the monotonicity the
design rules assert. One correction to the equations of the manuscript is
applied, in a case the manuscript keeps commented out and the paper never uses.
Where a feeder carries buffers at both of its drop points, the drop probability
of the backward point must carry the same factor `pi_c` of the forward buffer
that its own birth probability carries, since an item reaches the backward point
only if the forward buffer could not absorb it. Without that factor flow
conservation fails by as much as half an item per time step, and the prediction
is about three percentage points off against simulation. `MeshSim.cpp` in the
parent directory has the same error, and is likewise never exercised on that
case.

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
      --sequence FILE    destination labels from FILE, one per line,
                         instead of uniform draws

EXPERIMENT
  -T, --steps T          measured steps per replication          (default 1330000)
  -R, --reps R           independent replications                     (default 30)
  -w, --warmup W         steps discarded per replication
                                       (default max(20000, 10*(c+1)*L))
      --seed S           base seed                             (default 20260901)
  -t, --threads K        worker threads      (default: hardware threads - 2,
                                             capped at R)

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
