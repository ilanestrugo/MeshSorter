#!/bin/sh
# Build the simulator and reproduce every table and figure of the manuscript
# that this package is responsible for: Tables 1, 2, 5, 6 and 7, and Figures 4
# and 8.
#
#   ./run_all.sh            the full grids
#   ./run_all.sh --quick    the same grids at low precision, a few minutes
#
# Table 4 is not run here: certify_all.py enumerates 54,114 allocations and takes
# several hours.  Section 7 does not depend on it.  Split it across machines with
# --feeders and a separate CERTIFY_OUT per machine if you have more than one.
#
#   python3 certify_all.py                       Table 4, several hours
#   python3 sweep_structured.py --grid large      Section S8, about an hour
#   python3 approx_eval.py --dir results/large --tag large --scatter 15
#
# CXX and CXXFLAGS override the compiler and its flags, for example
#   CXXFLAGS="-O3 -march=native -std=c++17 -pthread" ./run_all.sh
# MESHSORTER_THREADS overrides the worker count; left unset, the simulator uses
# two fewer than the machine's hardware threads, capped at the replication count.
set -e
cd "$(dirname "$0")"
: "${CXX:=c++}"
: "${CXXFLAGS:=-O2 -std=c++17 -pthread}"
echo "building: $CXX $CXXFLAGS -o meshsorter meshsorter.cpp"
$CXX $CXXFLAGS -o meshsorter meshsorter.cpp
$CXX $CXXFLAGS -o certify_rep certify_rep.cpp
$CXX $CXXFLAGS -o sweep_rep   sweep_rep.cpp
python3 table1.py "$@"
python3 table2.py "$@"
python3 loops_buffered.py "$@"
python3 tableS4S5.py "$@"
python3 frontier.py "$@"
python3 feeder_returns.py "$@"

python3 approx_model.py
python3 sweep_structured.py "$@"
python3 approx_eval.py

# Table 6 reads the certified allocation at each budget from results/certify if
# that grid has been run and from results/structured otherwise, so it comes
# after sweep_structured.py.
python3 table6.py "$@"
echo
echo "LaTeX bodies   : results/table1.tex results/table2.tex results/loops_buffered.tex"
echo "whole tables   : results/loops_buffered_table.tex"
echo "cell-level CSV : results/table1.csv results/table2.csv"
echo "raw results    : results/raw/"
echo "approximation  : results/approx_eval.tex results/approx_eval.csv"
echo "load balancing : results/tableS4.tex results/tableS5.tex results/tableS4S5.csv"
echo "frontier data  : results/m4_n4_B10_*_efficiency_front*.dat"
echo "order-derived  : results/table6.tex results/table6.csv"
echo "scatter data   : results/approx_scatter_*.dat"
