#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_sweep.sh - buffer-allocation certification sweep
#
# Runs buffer_certify over a grid of feeder counts and per-primary-belt budgets.
# Every instance is independent, so instances are run JOBS at a time; each one is
# itself single-threaded.  Each instance writes its own summary CSV (concurrent
# appends to a shared file are not safe); the summaries are concatenated at the end.
#
# Defaults reproduce the September 2026 design:
#   pilot horizon 250,000 batches for the structured set S, 1,000 for competitors,
#   batch size 250, alpha 0.01, eps 0.001, kappa 3, minimum competitor budget 5,000.
#
# Usage:
#   ./run_sweep.sh                          # full grid, 4 jobs at a time
#   JOBS=12 ./run_sweep.sh                  # 12 instances in parallel
#   MS="5" BS="10" JOBS=1 ./run_sweep.sh    # a single instance
#   OUT=out_test BLOCKS_S=2000 BLOCKS_C=100 MIN_C=200 ./run_sweep.sh   # quick check
#
# Every parameter below can be overridden from the environment.
# ---------------------------------------------------------------------------
set -uo pipefail

BIN=${BIN:-./build/buffer_certify}
CLI=${CLI:-./build/meshsim_cli}
OUT=${OUT:-out}
JOBS=${JOBS:-4}

N=${N:-4}                       # primary belts
MS=${MS:-"3 4 5"}               # feeder belts
BS=${BS:-"10 9 8 7 6 5 4 3 2 1"}  # per-primary budgets (descending: long jobs start first)
DIRS=${DIRS:-"oneway twoway"}

BLOCKS_S=${BLOCKS_S:-250000}    # --Num_of_blocks    (pilot horizon for S)
BLOCKS_C=${BLOCKS_C:-1000}      # --Num_of_blocks_C  (pilot horizon for C)
BLOCK_SIZE=${BLOCK_SIZE:-250}
ALPHA=${ALPHA:-0.01}
EPS=${EPS:-0.001}
KAPPA=${KAPPA:-3}               # --m_mult
MIN_C=${MIN_C:-5000}            # --min_blocks_C
SEED0=${SEED0:-123}
B0=${B0:-1}                     # also produce the B=0 rows with meshsim_cli (set 0 to skip)

if [ ! -x "$BIN" ]; then
  echo "buffer_certify not found at $BIN - build first:" >&2
  echo "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel" >&2
  exit 1
fi
mkdir -p "$OUT"

run_one() {
  local m=$1 B=$2 dir=$3
  local bidir subset tag start rc
  if [ "$dir" = "oneway" ]; then bidir=false; subset=one_way_monotone
  else                           bidir=true;  subset=two_way_monotone_backward_only
  fi
  tag="m${m}_n${N}_B${B}_${dir}"
  start=$(date +%s)
  echo "[start] $tag  $(date -Is)"
  "$BIN" --m "$m" --n "$N" --B "$B" --subset "$subset" --bidirectional "$bidir" \
         --Num_of_blocks "$BLOCKS_S" --Num_of_blocks_C "$BLOCKS_C" --block_size "$BLOCK_SIZE" \
         --alpha "$ALPHA" --eps "$EPS" --seed0 "$SEED0" \
         --m_mult "$KAPPA" --min_blocks_C "$MIN_C" --fast true \
         --log_prefix "$OUT/$tag" --best_summary_csv "$OUT/summary_${tag}.csv" \
         > "$OUT/${tag}.log" 2>&1
  rc=$?
  echo "[done ] $tag  rc=$rc  $(( $(date +%s) - start ))s  $(grep -h '^epsilon_hat' "$OUT/${tag}.log" 2>/dev/null | tail -1)"
  return 0
}
export -f run_one
export BIN OUT N BLOCKS_S BLOCKS_C BLOCK_SIZE ALPHA EPS KAPPA MIN_C SEED0

echo "grid: m in [$MS] x B in [$BS] x [$DIRS], n=$N"
echo "pilot: S=$BLOCKS_S blocks, C=$BLOCKS_C blocks, block_size=$BLOCK_SIZE"
echo "plan : alpha=$ALPHA eps=$EPS kappa=$KAPPA min_blocks_C=$MIN_C seed0=$SEED0"
echo "jobs : $JOBS in parallel -> $OUT"
echo

for m in $MS; do for B in $BS; do for d in $DIRS; do printf '%s %s %s\n' "$m" "$B" "$d"; done; done; done \
  | xargs -P "$JOBS" -n 3 bash -c 'run_one "$@"' _

# ---- B = 0 rows ----
# With B=0 there is a single allocation, so the structured set has no complement and
# buffer_certify has nothing to test (epsilon_hat is 0 by construction). The throughput
# for those rows comes from a plain simulation at the same horizon, batch size and seed.
if [ "$B0" = "1" ]; then
  if [ -x "$CLI" ]; then
    echo
    echo "B=0 rows via $CLI (T = $((BLOCKS_S * BLOCK_SIZE)) time steps)"
    for m in $MS; do
      for d in $DIRS; do
        if [ "$d" = "oneway" ]; then bidir=false; else bidir=true; fi
        "$CLI" --NumberOfFeeders "$m" --NumberOfPrimary "$N" --BiDirectional "$bidir" \
               --T $((BLOCKS_S * BLOCK_SIZE)) --block_size "$BLOCK_SIZE" --seed "$SEED0" \
               --run_mode sim --fast --output_csv "$OUT/b0_rows.csv" \
               > "$OUT/m${m}_n${N}_B0_${d}.log" 2>&1
        echo "[done ] m${m}_n${N}_B0_${d}"
      done
    done
    echo "B=0 throughputs appended to $OUT/b0_rows.csv (column Throughput; epsilon_hat = 0)"
  else
    echo "meshsim_cli not found at $CLI - skipping the B=0 rows (build it or set B0=0)"
  fi
fi

# ---- collect ----
SUMMARY="$OUT/summary_all_runs.csv"
first=1
: > "$SUMMARY"
for f in "$OUT"/summary_m*_n*_B*_*.csv; do
  [ -e "$f" ] || continue
  if [ $first -eq 1 ]; then cat "$f" >> "$SUMMARY"; first=0
  else tail -n +2 "$f" >> "$SUMMARY"
  fi
done
echo
echo "sweep complete. combined summary: $SUMMARY  ($(( $(wc -l < "$SUMMARY") - 1 )) instances)"
echo "per-instance logs and phase CSVs are in $OUT/"
