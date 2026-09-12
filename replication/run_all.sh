#!/bin/sh
# Build the simulator and reproduce every table and figure of the manuscript
# that this package is responsible for: Tables 1, 2 and 5, and Figure 4, and,
# when the certification grid has been run, Table 7 and Figure 8.
#
#   ./run_all.sh            the full grids
#   ./run_all.sh --quick    the same grids at low precision, a few minutes
#
# Table 4 is not run here: certify_all.py enumerates 54,120 allocations and takes
# several hours.  Section 7 does not depend on it.
#
#   python3 certify_all.py                       Table 4, several hours
#   python3 sweep_structured.py --grid large      Section S7, about an hour
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
python3 feeder_returns.py "$@"

python3 approx_model.py
python3 sweep_structured.py "$@"
python3 approx_eval.py
echo
echo "LaTeX bodies   : results/table1.tex results/table2.tex results/loops_buffered.tex"
echo "whole tables   : results/loops_buffered_table.tex"
echo "cell-level CSV : results/table1.csv results/table2.csv"
echo "raw results    : results/raw/"
echo "approximation  : results/approx_eval.tex results/approx_eval.csv"
echo "scatter data   : results/approx_scatter_*.dat"
