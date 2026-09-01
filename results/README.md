# Reported experiments, September 2026

Output of the buffer-allocation certification sweep and of the high-precision
re-examination that follow the procedures in Sections S2 and S8 of the
supplemental material.

## Files

| File | Contents |
| --- | --- |
| `summary_all_runs.csv` | One row per certification instance (60 rows): selected structured allocation, its throughput and confidence half-widths, the worst competitor, and `epsilon_hat`. |
| `b0_rows.csv` | Throughput for the six zero-budget configurations, produced with `meshsim_cli` (with `B = 0` the structured set has no complement, so the certification procedure does not apply). |
| `m*_recheck.csv` | Re-examination of the leading competitors for the seven instances with `epsilon_hat > 0.001`: each candidate re-simulated at 250,000 batches on an independent seed range, with its gap to the reference allocation, both test statistics and `eps_c_precise`. |
| `sweep.log` | Console log of the sweep, with per-instance timings. |
| `full_run_2026-09-01.tgz` | The complete `out/` directory: per-instance Phase-1 and Phase-3 detail CSVs and logs. |

## How they were produced

Sweep (60 instances plus the six `B = 0` rows), about 14 CPU-hours:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel
JOBS=12 ./run_sweep.sh
```

which is the default grid: `n = 4`, `m` in {3,4,5}, `B` from 1 to 10, both drop
mechanisms, Phase-1 pilot of 250,000 batches for the structured set and 1,000
for the competitors, batch size 250, `alpha = 0.01`, `eps = 0.001`,
`kappa = 3` (`--m_mult`), minimum competitor budget 5,000 batches, base seed 123.

Re-examination of the seven instances with `epsilon_hat > 0.001`:

```bash
for spec in "4 9" "5 4" "5 6" "5 7" "5 8" "5 9" "5 10"; do
  set -- $spec
  ./build/recheck_candidates \
    --phase3_csv   out/m$1_n4_B$2_twoway_phase3_detail.csv \
    --phase1_S_csv out/m$1_n4_B$2_twoway_phase1_S_detail.csv \
    --m $1 --n 4 --bidirectional true --top 50 \
    --blocks 250000 --block_size 250 --alpha 0.01 --eps 0.001 --threads 0 \
    --out out/m$1_n4_B$2_twoway_recheck.csv
done
```
