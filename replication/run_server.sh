#!/bin/sh
# The expensive half of the replication, for a Linux server with ~30 cores.
#
# Everything here runs under the geometry of geometry.py: consecutive primary
# belts 5 m apart, a feeder loop of 20n + 24 slots, which is 104 slots at four
# belts, the machine drawn in Figure 3.  The earlier scripts used 4n + 8 or
# 8n + 8 slots, neither of which is buildable at the bin and aisle dimensions
# the paper costs.
#
#   ./run_server.sh              everything below, many hours
#   ./run_server.sh --quick      the same at low precision, for a smoke test
#
# Both certify_all.py and sweep_structured.py cache per cell and record the
# design that produced it, so an interrupted run resumes where it stopped and a
# cell computed under a different geometry is recomputed rather than reused.
# Killing and restarting this script is safe.
#
# CXX and CXXFLAGS override the compiler:
#   CXXFLAGS="-O3 -march=native -std=c++17 -pthread" ./run_server.sh
set -e
cd "$(dirname "$0")"

: "${CXX:=c++}"
: "${CXXFLAGS:=-O2 -std=c++17 -pthread}"
: "${MESHSORTER_THREADS:=30}"
export MESHSORTER_THREADS

echo "=============================================================="
echo " geometry in force"
echo "=============================================================="
python3 geometry.py
echo

echo "=============================================================="
echo " building"
echo "=============================================================="
echo "$CXX $CXXFLAGS"
$CXX $CXXFLAGS -o meshsorter  meshsorter.cpp
$CXX $CXXFLAGS -o certify_rep certify_rep.cpp
$CXX $CXXFLAGS -o sweep_rep   sweep_rep.cpp
echo "built meshsorter, certify_rep, sweep_rep"
echo

echo "=============================================================="
echo " 1/4  certification grid      Table 3 and Table S8"
echo "      60 cells, up to 19,448 allocations in the largest"
echo "=============================================================="
time python3 certify_all.py "$@"
echo

echo "=============================================================="
echo " 2/4  structured sweep, four belts     inputs to Table 5"
echo "      388 allocations in 60 classes"
echo "=============================================================="
time python3 sweep_structured.py "$@"
echo

echo "=============================================================="
echo " 3/4  structured sweep, fifteen belts  Section S12"
echo "      664 allocations in four designs"
echo "=============================================================="
time python3 sweep_structured.py --grid large "$@"
echo

echo "=============================================================="
echo " 4/4  approximation against both sweeps"
echo "=============================================================="
python3 approx_model.py
python3 approx_eval.py
python3 approx_eval.py --dir results/large --tag large --scatter 15
# Table S13, the order-derived check, is run by feeder_sequence.py from the
# committed allocations in results/table6_allocations.txt; it takes minutes and
# is part of run_all.sh, so it is not repeated here.
echo

echo "=============================================================="
echo " packaging"
echo "=============================================================="
STAMP=$(date +%Y%m%d-%H%M)
tar czf "results-newgeom-$STAMP.tar.gz" results geometry.py
echo "wrote results-newgeom-$STAMP.tar.gz"
echo
echo "Send back that archive.  The files that carry the numbers are:"
echo "  results/certify/*.json          the certification grid, Table 3 and S8"
echo "  results/structured/*.json       the four-belt structured sweep"
echo "  results/large/*.json            the fifteen-belt sweep, Section S12"
echo "  results/approx_eval*.csv|.tex   Table 5, Tables S14 and S15, Figure 9"
