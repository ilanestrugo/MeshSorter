#!/bin/sh
# Build the simulator and reproduce the tables and figures of the manuscript that
# are quick to run.  Table numbers are those of the release tagged
# tre-submission; the table in README.md says which script produces which.
#
#   ./run_all.sh            the full grids
#   ./run_all.sh --quick    the same grids at low precision, a few minutes
#
# Not run here, because they take hours: the certification of Table 3 and
# Table S8 (certify_all.py enumerates 54,114 allocations), and the fifteen-belt
# systems of Section S11.  Split the certification across machines with
# --feeders and a separate CERTIFY_OUT per machine if you have more than one.
#
#   python3 certify_all.py                       Table 3 and Table S8, several hours
#   python3 sweep_structured.py --grid large      Section S11, about an hour
#   python3 approx_eval.py --dir results/large --tag large --scatter 15
#
# The geometry of the paper, loops of 20n + 24 slots with consecutive primary
# belts 10 slots apart, is set in geometry.py, which every script reads.  The
# simulator's own built-in defaults describe an earlier, shorter loop, so a
# bare "./meshsorter -n 4 -m 4 --dual" does NOT run the paper's system; add
#   --spacing 10 --turn 12 --df 14 --width 10
#
# CXX and CXXFLAGS override the compiler and its flags, for example
#   CXXFLAGS="-O3 -march=native -std=c++17 -pthread" ./run_all.sh
# MESHSORTER_THREADS overrides the worker count; left unset, the simulator uses
# two fewer than the machine's hardware threads, capped at the replication count.
# The results do not depend on the worker count or on the platform.
set -e
cd "$(dirname "$0")"
: "${CXX:=c++}"
: "${CXXFLAGS:=-O2 -std=c++17 -pthread}"
echo "building: $CXX $CXXFLAGS -o meshsorter meshsorter.cpp"
$CXX $CXXFLAGS -o meshsorter meshsorter.cpp
$CXX $CXXFLAGS -o certify_rep certify_rep.cpp
$CXX $CXXFLAGS -o sweep_rep   sweep_rep.cpp
$CXX $CXXFLAGS -o transient_pilot transient_pilot.cpp

python3 exact_single_drop.py              # Table 2(a)
python3 unbuffered_grid.py "$@"           # Table 2(b)
python3 stagger_grid.py "$@"              # Tables S5 and S6
python3 loops_buffered.py "$@"            # Table S11
python3 load_balance.py "$@"              # Tables S9 and S10
python3 frontier.py "$@"                  # Figure S4
python3 feeder_returns.py "$@"            # Figure S2
python3 three_way_comparison.py           # Table 4

python3 approx_model.py
python3 sweep_structured.py "$@"
python3 approx_eval.py                    # Table 5, Tables S13 and S14, Figure S5

# Table S12 reads the certified allocation at each budget from results/certify if
# that grid has been run, then from results/table6_allocations.txt, and from
# results/structured only as a last resort.  It comes after sweep_structured.py.
python3 order_derived.py "$@"             # Table S12
python3 sequence_dependence.py            # statistics quoted in Section S9
python3 order_shuffled.py                 # Table S17, about an hour
python3 asymmetry_check.py                # Table S16

python3 transient_pilot.py "$@"            # transient figures of Section S2.2, about half an hour
python3 buffer_relaxation.py "$@"          # buffer-content figures of Section S2.2, a few minutes
python3 warmup_sensitivity.py "$@"         # Table S3, about forty minutes

python3 gcd_verification.py               # the enumeration behind Proposition 5
python3 validate.py                       # the simulator against the exact values
echo
echo "LaTeX bodies   : results/table1.tex results/table2.tex results/loops_buffered.tex"
echo "whole tables   : results/loops_buffered_table.tex"
echo "cell-level CSV : results/table1.csv results/table2.csv"
echo "raw results    : results/raw/"
echo "approximation  : results/approx_eval.tex results/approx_eval.csv"
echo "load balancing : results/load_balance_single.tex results/load_balance_dual.tex results/load_balance.csv"
echo "frontier data  : results/m4_n4_B10_*_efficiency_front*.dat"
echo "order-derived  : results/table6.tex results/table6.csv"
echo "shuffled       : results/order_shuffled.tex results/order_shuffled.csv"
echo "asymmetry      : results/asymmetry.tex results/asymmetry.csv"
echo "sequence       : results/sequence_dependence.json"
echo "warm-up         : results/transient_pilot/pilot_summary.json results/warmup_sensitivity.tex"
echo "scatter data   : results/approx_scatter_*.dat"
