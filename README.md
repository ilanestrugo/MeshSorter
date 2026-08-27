# MeshSorter

MeshSorter is a research implementation of a two-layer conveyor-sortation model. The repository contains a synchronous discrete-time simulator, an analytical throughput approximation, and a program for evaluating structured buffer-allocation rules.

The software accompanies the manuscript *MeshSorter: A Two-Layer Conveyor Architecture for High-Throughput Sortation*.

## Repository contents

| File | Purpose |
| --- | --- |
| `MeshSim.hpp`, `MeshSim.cpp` | C++ simulation library and analytical approximation |
| `MeshSim_cli.cpp` | Command-line interface for simulation and approximation runs |
| `buffer_hypothesis_certify_3phases_iut.cpp` | Three-phase buffer-allocation certification experiment |
| `CMakeLists.txt` | Portable CMake build configuration |

The build produces two executables:

- `meshsim_cli`: runs the simulator, the analytical approximation, or both.
- `buffer_certify`: enumerates buffer allocations and evaluates a specified structured subset.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler

The project has been tested on macOS with AppleClang. The same CMake commands are suitable for current Linux distributions with GCC or Clang.

On Ubuntu or Debian, install the build tools with:

```bash
sudo apt update
sudo apt install -y cmake build-essential
```

On macOS, install the Xcode command-line tools and CMake if needed:

```bash
xcode-select --install
brew install cmake
```

## Build

```bash
git clone https://github.com/ilanestrugo/MeshSorter.git
cd MeshSorter

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executables will be created as `build/meshsim_cli` and `build/buffer_certify`.

Either executable can be excluded from the build:

```bash
cmake -S . -B build \
  -DBUILD_MESHSIM_CLI=ON \
  -DBUILD_BUFFER_CERTIFY=OFF
```

## Quick start

The following command runs the simulator and analytical approximation for an unbuffered, dual-drop, 4-feeder by 4-primary system. It uses the manuscript's main simulation horizon, seed, and batch size.

```bash
mkdir -p results

./build/meshsim_cli \
  --NumberOfFeeders 4 \
  --NumberOfPrimary 4 \
  --BiDirectional true \
  --buffer_capacity '[[0,0,0,0],[0,0,0,0]]' \
  --T 250000 \
  --seed 123 \
  --block_size 250 \
  --run_mode both \
  --ci_enable \
  --ci_level 0.95 \
  --output_csv results/mesh_4x4_dual_unbuffered.csv
```

A successful run prints the detected warm-up, post-warm-up throughput and utilization, confidence-interval half-widths, and analytical approximation. It also appends one row to the requested CSV file.

Exact pseudorandom results can differ across C++ standard-library implementations. Record the compiler, platform, commit hash, command, and seed when archiving experiments.

### Analytical model only

`--T` remains a required argument in model-only mode; use `--T 0` when no simulation is requested.

```bash
./build/meshsim_cli \
  --NumberOfFeeders 4 \
  --NumberOfPrimary 4 \
  --BiDirectional false \
  --buffer_capacity '[[0,1,3,6],[0,0,0,0]]' \
  --T 0 \
  --run_mode model \
  --output_csv results/model_4x4_single_drop.csv
```

## Buffer-capacity encoding

`--buffer_capacity` describes two vectors of length `m`, where `m` is the number of feeder belts:

```text
[[forward capacities by feeder], [backward capacities by feeder]]
```

For example:

```bash
--buffer_capacity '[[0,1,2,7],[0,0,0,0]]'
```

assigns forward-drop capacities of 0, 1, 2, and 7 to feeders 1 through 4. The same feeder-level capacity is used at that feeder's intersection with every primary belt.

Accepted quoted forms include:

```text
'[[0,1,2],[0,0,0]]'
'0,1,2;0,0,0'
'0,1,2|0,0,0'
'0,1,2'
```

The last form specifies the forward vector only. With `--BiDirectional false`, only the forward vector is used. The implementation forces the first feeder's capacities to zero, fills omitted entries with zero, ignores extra entries, and clamps negative capacities to zero.

The allocation budget `B` used in the manuscript is per primary belt. Therefore, if the two vectors sum to `B` and the system has `n` primary belts, the corresponding system-wide storage capacity is `nB` units.

## `meshsim_cli` options

Flag names are case-sensitive.

| Option | Meaning | Default |
| --- | --- | --- |
| `--NumberOfFeeders <int>` | Number of feeder belts, `m` | required |
| `--NumberOfPrimary <int>` | Number of primary belts, `n` | required |
| `--T <int>` | Simulation horizon in time steps | required |
| `--BiDirectional <bool>` | `true` for dual-drop operation; `false` for single-drop | `true` |
| `--buffer_capacity <str>` | Forward/backward buffer vectors | all zero |
| `--seed <int>` | Mersenne Twister seed | `0` |
| `--block_size <int>` | Batch size in time steps | `100` |
| `--run_mode <mode>` | `sim`, `model`, or `both` | `both` |
| `--ci_enable` | Report a throughput confidence interval over post-warm-up batches | disabled |
| `--ci_level <level>` | Confidence level: `0.90`, `0.95`, or `0.99`; also enables the interval | `0.95` |
| `--no_buffer_stats`, `--fast` | Skip buffer-statistics collection; buffers remain active | disabled |
| `--Verbal` | Print detailed blocking events | disabled |
| `--output_csv <path>` | CSV output path | `MeshSim_results.csv` |
| `--help` | Print command-line help | — |

The program appends to an existing output CSV and writes a header only when the file is absent or empty. Its parent directory must already exist. Simulation statistics are computed after automatic warm-up deletion. Model-only runs leave simulation fields empty or zero; simulation-only runs leave model fields empty.

## Buffer-allocation certification

`buffer_certify` enumerates every allocation of budget `B` across feeder positions other than the first and compares a structured subset with the remaining allocations. Runtime grows combinatorially with `m` and `B`, especially for dual-drop systems.

This small smoke test verifies the executable and output pipeline. Its intentionally short horizon and loose statistical settings are **not** the manuscript's research settings.

```bash
mkdir -p results

./build/buffer_certify \
  --m 3 \
  --n 4 \
  --B 2 \
  --bidirectional false \
  --subset one_way_monotone \
  --Num_of_blocks 20 \
  --block_size 50 \
  --eps 0.10 \
  --alpha 0.10 \
  --n0 1 \
  --m_mult 0.10 \
  --min_blocks_C 10 \
  --seed0 123 \
  --fast true \
  --log_prefix results/smoke \
  --best_summary_csv results/smoke_summary.csv
```

Required arguments are `--m`, `--n`, `--B`, and `--subset`. Pair the subset and direction as follows:

| Direction | `--bidirectional` | Valid structured subset |
| --- | --- | --- |
| Single-drop | `false` | `one_way_monotone` |
| Dual-drop, backward buffers only | `true` | `two_way_monotone_backward_only` |
| Dual-drop, both directions monotone | `true` | `two_way_monotone` |

Additional options are:

| Option | Meaning | Default |
| --- | --- | --- |
| `--Num_of_blocks <int>` | Initial number of simulation batches | `150` |
| `--block_size <int>` | Time steps per batch | `200` |
| `--alpha <float>` | One-sided significance level | `0.01` |
| `--eps <float>` | Target relative indifference-zone gap | `0.001` |
| `--n0 <int>` | Number of independent pilot runs per allocation in Phase 1 | `1` |
| `--m_mult <float>` | Sampling-effort multiplier | `1.0` |
| `--min_blocks_C <int>` | Minimum retained-batch budget for outside competitors | `500` |
| `--seed0 <int>` | Base seed | `123` |
| `--fast <bool>` | Skip buffer-statistics collection | `false` |
| `--log_prefix <path>` | Prefix for phase-detail CSV files | `run` |
| `--best_summary_csv <path>` | Optional append-only run-summary CSV | none |

Each run overwrites `<log_prefix>_phase1_S_detail.csv` and `<log_prefix>_phase3_detail.csv`. If requested, `--best_summary_csv` appends a summary row and creates its parent directory.

## Model scope and current reproducibility status

The simulator represents an idealized, synchronized system:

- Destinations are sampled uniformly and independently by the public CLI.
- One belt-slot advance corresponds to one time step.
- Feeder admissions, transfers, and belt advances occur in a fixed event order.
- No external admission queues are modeled.
- Final diverts from primary belts are unconstrained.
- Transfers are instantaneous and error-free within the discrete-time model.

This repository currently contains the core simulator, analytical approximation, and buffer-allocation certification program. It is **not yet an end-to-end manuscript reproduction archive**. In particular, this snapshot does not include the processed Olist destination sequence, experiment manifests for every table and figure, stored reference outputs, plotting scripts, automated tests, or a locked software environment. The commands above demonstrate the current executables but do not reproduce every manuscript result.

## Citation

The formal article citation will be added when it becomes available. Until then, cite the repository URL and the exact commit hash used for the analysis.

## License

No software license is currently included. Until a license is added, standard copyright restrictions apply.

## Contact

Questions and reproducibility reports can be submitted through the repository's GitHub issue tracker.
