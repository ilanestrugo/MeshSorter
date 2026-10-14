#!/usr/bin/env python3
"""Finite-record destination dependence experiment, with separate feeder cursors.

Default: four belts, four feeders, dual drop, certified B=0..10 allocations,
30 paired replications, 20,000 uniform warm-up steps, then one record length
(98,816) of measurement steps. Three conditions use identical warm-up seeds:
chronological record, one new random permutation per replication, and uniform
draws. Chronological and shuffled conditions also share the four random offsets.

Record cursors advance only on admission and never revisit a position within a
feeder. Every shuffled replication uses ONE permutation shared by its feeders.
The 2,000/5,000-step exclusions re-estimate the same measurement trajectories;
they do not consume extra labels or restart the system.

All output is separate from the published table6/order_shuffled results.
Intervals are conditional on this empirical record. The shuffled condition
additionally averages over permutations; it is a control, not a new dataset.

    python3 feeder_sequence.py --build             full experiment
    python3 feeder_sequence.py --quick             wiring check
    python3 feeder_sequence.py --budgets 0,5,10     selected budgets
"""
import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shlex
import statistics
import subprocess
import sys
import tempfile
import time

import geometry

HERE = Path(__file__).resolve().parent
CONDITIONS = ("chronological", "shuffled", "uniform")


def integer_list(text):
    try:
        values = sorted(set(int(v.strip()) for v in text.split(",")))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from exc
    if not values:
        raise argparse.ArgumentTypeError("list must not be empty")
    return values


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def atomic_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(text)
    temporary.replace(path)


def write_json(path, value):
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_csv(path, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, stream.getvalue())


def read_labels(args):
    source = (args.sequence or args.orders).resolve()
    if args.sequence:
        labels = [int(line.strip()) for line in source.read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith("#")]
    else:
        with source.open(newline="", encoding="utf-8-sig") as stream:
            records = [(int(row["arrival_seq"]), int(row["belt"]))
                       for row in csv.DictReader(stream)]
        records.sort()
        if len({index for index, _ in records}) != len(records):
            raise ValueError("duplicate arrival_seq values in prepared orders")
        labels = [label for _, label in records]
    if not labels or any(label not in range(1, 5) for label in labels):
        raise ValueError("record must contain destination labels in 1..4")
    return source, labels


def read_allocations(path, budgets):
    entries = [line.strip() for line in path.read_text().splitlines()
               if line.strip() and not line.lstrip().startswith("#")]
    allocations = {}
    for budget in budgets:
        if budget < 0 or budget >= len(entries):
            raise ValueError("budget {} is missing from {}".format(budget, path))
        capacities = [int(value) for value in entries[budget].split(",")]
        if len(capacities) != 4 or min(capacities) < 0 or sum(capacities) != budget:
            raise ValueError("invalid allocation for budget {}".format(budget))
        if capacities[0] != 0:
            raise ValueError("the first feeder must have zero buffers")
        allocations[budget] = capacities
    return allocations


def buffer_spec(capacities):
    return ",".join(str(value) for capacity in capacities
                    for value in ([0] * 4 + [capacity] * 4))


def estimates(values, critical):
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    halfwidth = critical * sd / math.sqrt(len(values))
    return mean, sd, halfwidth


def observations(result, discard):
    point = result["checkpoints"].index(discard)
    return [(sum(total) - sum(prefix[point])) / (result["steps"] - discard)
            for total, prefix in zip(result["admissions_by_replication"],
                                     result["checkpoint_admissions"])]


def validate_result(result, condition, args, length):
    for field, expected in (("n", 4), ("m", 4), ("dual", True),
                            ("steps", args.steps), ("warmup", args.warmup),
                            ("reps", args.reps), ("seed", args.seed),
                            ("loops", [geometry.loop_length(4)] * 4),
                            ("checkpoints", args.discard_steps), ("warmup_input", "uniform")):
        if result.get(field) != expected:
            raise ValueError("{}: unexpected {}".format(condition, field))
    totals = result["admissions_by_replication"]
    if len(totals) != args.reps or any(len(row) != 4 for row in totals):
        raise ValueError("missing per-replication admission counts")
    if (len(result["checkpoint_admissions"]) != args.reps or
            len(result["warmup_state_hashes"]) != args.reps):
        raise ValueError("missing per-replication checkpoint or warm-up audit")
    if any(row[0] != args.steps or any(v < 0 or v > args.steps for v in row)
           for row in totals):
        raise ValueError("admission counts violate feeder limits")
    for total, prefixes in zip(totals, result["checkpoint_admissions"]):
        if len(prefixes) != len(args.discard_steps):
            raise ValueError("missing measurement checkpoints")
        previous = [0] * 4
        for discard, prefix in zip(args.discard_steps, prefixes):
            if len(prefix) != 4 or any(not (a <= b <= c and b <= discard)
                                      for a, b, c in zip(previous, prefix, total)):
                raise ValueError("invalid cumulative checkpoint counts")
            previous = prefix
    if condition == "uniform":
        if result.get("input_mode") != "uniform":
            raise ValueError("uniform control did not use uniform input")
    else:
        if (result.get("input_mode") != "separate_feeders" or
                result.get("sequence_length") != length or
                result.get("sequence_shuffled") != (condition == "shuffled") or
                result.get("sequence_consumed") != totals):
            raise ValueError("record condition used the wrong consumption protocol")
        if any(len(result[field]) != args.reps for field in
               ("feeder_start_positions", "sequence_wraps", "sequence_belt_counts")):
            raise ValueError("missing per-replication sequence audit")
        for counts, per_belt in zip(totals, result["sequence_belt_counts"]):
            if len(per_belt) != 4 or any(len(row) != 4 or min(row) < 0 or sum(row) != count
                                        for row, count in zip(per_belt, counts)):
                raise ValueError("destination counts do not match admissions")
        for starts, counts, wraps in zip(result["feeder_start_positions"], totals,
                                         result["sequence_wraps"]):
            if len(starts) != 4 or len(wraps) != 4:
                raise ValueError("missing cursor audit")
            for start, count, crossings in zip(starts, counts, wraps):
                expected = (start + count - 1) // length if count else 0
                if not (0 <= start < length and count <= length and crossings == expected):
                    raise ValueError("record position reuse or incorrect wrapping")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument("--orders", type=Path, default=HERE.parent / "olist_orders_ForRun.csv")
    sources.add_argument("--sequence", type=Path, help="use an existing destination-label file")
    parser.add_argument("--allocs", type=Path, default=HERE / "results/table6_allocations.txt")
    parser.add_argument("--budgets", type=integer_list)
    parser.add_argument("--reps", type=int)
    parser.add_argument("--steps", type=int, help="default: record length; must not exceed it")
    parser.add_argument("--warmup", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--discard-steps", type=integer_list, default=[0, 2000, 5000])
    parser.add_argument("--threads", type=int, default=int(os.environ.get("MESHSORTER_THREADS", "0")))
    parser.add_argument("--binary", type=Path, default=Path(os.environ.get("MESHSORTER_BIN", str(HERE / "meshsorter"))))
    parser.add_argument("--build", action="store_true", help="compile the simulator before running")
    parser.add_argument("--quick", action="store_true", help="4 replications, 10,000 steps, B=0,5,10")
    parser.add_argument("--out", type=Path, help="output root; each protocol has its own run directory")
    args = parser.parse_args()
    source, labels = read_labels(args)
    args.reps = args.reps if args.reps is not None else (4 if args.quick else 30)
    args.steps = args.steps if args.steps is not None else (min(10_000, len(labels)) if args.quick else len(labels))
    args.budgets = args.budgets if args.budgets is not None else ([0, 5, 10] if args.quick else list(range(11)))
    args.discard_steps = sorted(set([0] + args.discard_steps))
    if args.reps < 2 or not (1 <= args.steps <= len(labels)) or args.warmup < 0:
        parser.error("need reps >= 2, warmup >= 0, and 1 <= steps <= record length")
    if any(point < 0 or point >= args.steps for point in args.discard_steps):
        parser.error("every discard point must be between 0 and steps-1")
    if args.threads < 0 or not (0 <= args.seed <= 2**64 - 1):
        parser.error("threads must be nonnegative and seed must fit an unsigned 64-bit integer")
    allocations = read_allocations(args.allocs, args.budgets)
    args.binary = args.binary.resolve()
    if args.build:
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) +
                       ["-O2", "-std=c++17", "-pthread", "-o", str(args.binary),
                        str(HERE / "meshsorter.cpp")], check=True)
    help_text = subprocess.run([str(args.binary), "--help"], check=True,
                               capture_output=True, text=True).stdout
    if "--feeder-sequence" not in help_text or "--checkpoints" not in help_text:
        parser.error("simulator is out of date; rerun with --build")

    label_text = "".join("{}\n".format(label) for label in labels)
    protocol = dict(experiment="separate_feeder_finite_record_v1", belts=4, feeders=4,
                    dual=True, geometry=dict(spacing=geometry.SPACING, turn=geometry.TURN,
                                             df=geometry.DF, width=geometry.WIDTH),
                    warmup=args.warmup, warmup_input="uniform", steps=args.steps,
                    reps=args.reps, seed=args.seed, discard_steps=args.discard_steps,
                    allocations=allocations, sequence_length=len(labels),
                    sequence_counts=[labels.count(i) for i in range(1, 5)],
                    sequence_sha256=sha256(label_text.encode()), source=str(source),
                    source_sha256=sha256(source.read_bytes()),
                    allocation_source=str(args.allocs.resolve()),
                    allocation_sha256=sha256(args.allocs.read_bytes()),
                    code_sha256={name: sha256((HERE / name).read_bytes()) for name in
                                 ("meshsorter.cpp", "meshsorter_core.hpp", "feeder_sequence.py", "geometry.py")},
                    binary_sha256=sha256(args.binary.read_bytes()), conditions=CONDITIONS,
                    shuffled_control="one independent permutation per replication, shared by all feeders",
                    interval_scope="conditional on this record; shuffled control also varies its permutation")
    key = sha256(json.dumps(protocol, sort_keys=True).encode())
    root = (args.out or HERE / "results" / ("feeder_sequence_quick" if args.quick else "feeder_sequence")).resolve()
    out = root / key[:16]
    write_json(out / "manifest.json", dict(protocol, experiment_key=key))
    sequence = out / "inputs/chronological.txt"
    atomic_text(sequence, label_text)
    print("{} labels; {} uniform warm-up steps; {} measurement steps; {} replications".format(
          len(labels), args.warmup, args.steps, args.reps), file=sys.stderr, flush=True)
    print("Output: {}".format(out), file=sys.stderr, flush=True)
    summary, replication_rows, feeder_rows = [], [], []
    started = time.monotonic()
    for budget, capacities in allocations.items():
        results = {}
        for condition in CONDITIONS:
            command = [str(args.binary), "-n", "4", "-m", "4", "--dual",
                       "--spacing", str(geometry.SPACING), "--turn", str(geometry.TURN),
                       "--df", str(geometry.DF), "--width", str(geometry.WIDTH),
                       "-b", buffer_spec(capacities), "-T", str(args.steps), "-R", str(args.reps),
                       "-w", str(args.warmup), "--seed", str(args.seed), "--json", "--per-feeder",
                       "--checkpoints", ",".join(map(str, args.discard_steps))]
            if args.threads:
                command += ["--threads", str(args.threads)]
            if condition != "uniform":
                command += ["--feeder-sequence", str(sequence)]
            if condition == "shuffled":
                command += ["--shuffle-sequence"]
            raw_path = out / "raw/B{:02d}_{}.json".format(budget, condition)
            if raw_path.exists():
                cached = json.loads(raw_path.read_text())
                if cached["experiment_key"] != key:
                    raise ValueError("cache protocol mismatch")
                result = cached["result"]
            else:
                print("  B={:2d} {}".format(budget, condition), file=sys.stderr, flush=True)
                completed = subprocess.run(command, capture_output=True, text=True, check=True)
                result = json.loads(completed.stdout)
                write_json(raw_path, dict(experiment_key=key, command=command, result=result))
            validate_result(result, condition, args, len(labels))
            results[condition] = result
        if results["chronological"]["feeder_start_positions"] != results["shuffled"]["feeder_start_positions"]:
            raise ValueError("chronological and shuffled offsets are not paired")
        hashes = results["uniform"]["warmup_state_hashes"]
        if any(result["warmup_state_hashes"] != hashes for result in results.values()):
            raise ValueError("conditions did not start measurement in the same physical and RNG states")
        for discard in args.discard_steps:
            row = dict(B=budget, allocation=",".join(map(str, capacities)),
                       discard_steps=discard, measured_steps=args.steps-discard, replications=args.reps)
            values = {condition: observations(result, discard) for condition, result in results.items()}
            critical = results["uniform"]["critical_value"]
            for condition, data in values.items():
                mean, sd, hw = estimates(data, critical)
                row.update({condition+"_mean": mean, condition+"_sd": sd, condition+"_halfwidth": hw})
                for r, value in enumerate(data):
                    replication_rows.append(dict(B=budget, condition=condition, replication=r+1,
                                                 discard_steps=discard, throughput=value))
            for control in ("shuffled", "uniform"):
                delta = [a-b for a, b in zip(values["chronological"], values[control])]
                mean, sd, hw = estimates(delta, critical)
                row.update({"chronological_minus_"+control: mean,
                            "paired_"+control+"_sd": sd, "paired_"+control+"_halfwidth": hw})
            summary.append(row)
        for condition, result in results.items():
            for r, counts in enumerate(result["admissions_by_replication"]):
                for j, count in enumerate(counts):
                    empirical = condition != "uniform"
                    feeder_rows.append(dict(B=budget, condition=condition, replication=r+1, feeder=j+1,
                                            admissions=count, loader_utilization=count/args.steps,
                                            start_position=result["feeder_start_positions"][r][j] if empirical else "",
                                            sequence_consumed=result["sequence_consumed"][r][j] if empirical else "",
                                            end_to_start_transitions=result["sequence_wraps"][r][j] if empirical else "",
                                            belt_counts=json.dumps(result["sequence_belt_counts"][r][j]) if empirical else "",
                                            warmup_state_hash=result["warmup_state_hashes"][r]))
        print("  B={:2d}: chronological {:.6f}, shuffled {:.6f}, uniform {:.6f}".format(
              budget, *(statistics.mean(observations(results[c], 0)) for c in CONDITIONS)),
              file=sys.stderr, flush=True)
        write_csv(out / "summary.csv", summary)
        write_csv(out / "replications.csv", replication_rows)
        write_csv(out / "feeders.csv", feeder_rows)
    write_json(out / "summary.json", dict(experiment_key=key, rows=summary,
                                         elapsed_seconds=round(time.monotonic()-started, 3)))
    print("Completed. Summary: {}".format(out / "summary.csv"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print("feeder_sequence: {}".format(error), file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError) and error.stderr:
            print(error.stderr, file=sys.stderr)
        sys.exit(1)
