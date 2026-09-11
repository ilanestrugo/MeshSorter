#!/bin/sh
# Build the simulator and reproduce Tables 1 and 2 of the manuscript.
#
#   ./run_all.sh            the full grids
#   ./run_all.sh --quick    the same grids at low precision, a few minutes
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
python3 table1.py "$@"
python3 table2.py "$@"
echo
echo "LaTeX bodies   : results/table1.tex results/table2.tex"
echo "cell-level CSV : results/table1.csv results/table2.csv"
echo "raw results    : results/raw/"
