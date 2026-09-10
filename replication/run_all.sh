#!/bin/sh
# Build the simulator and reproduce Tables 1 and 2 of the manuscript.
# On two cores this takes about forty minutes.  Add --quick to both scripts
# for a low-precision run that finishes in a few minutes.
set -e
cd "$(dirname "$0")"
c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp
python3 table1.py "$@"
python3 table2.py "$@"
echo
echo "LaTeX bodies are in results/table1.tex and results/table2.tex"
