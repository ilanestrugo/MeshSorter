#!/usr/bin/env python3
"""Protocol checks for the separate-feeder experiment (standard library only).

Run: python3 -m unittest -v test_feeder_sequence.py
The tests compile the simulator into a temporary directory and use small
records; they never alter manuscript result files or the installed binary.
"""
from collections import Counter
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent


class FeederSequenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="meshsorter-sequence-test-")
        cls.root = Path(cls.temporary.name)
        cls.binary = cls.root / "meshsorter"
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) +
                       ["-O2", "-std=c++17", "-Wall", "-Wextra", "-Wpedantic", "-pthread",
                        "-o", str(cls.binary), str(HERE / "meshsorter.cpp")], check=True)
        cls.labels = [1 + (i // 17 + i // 113) % 4 for i in range(2048)]
        cls.record = cls.root / "record.txt"
        cls.record.write_text("".join("{}\n".format(label) for label in cls.labels))

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def run_simulator(self, extra=(), steps=2048, threads=1, checkpoints=True):
        command = [str(self.binary), "-n", "4", "-m", "4", "--dual",
                   "-b", ",".join(map(str, [0]*8 + [0]*4 + [1]*4 + [0]*4 + [2]*4 + [0]*4 + [7]*4)),
                   "-T", str(steps), "-w", "500", "-R", "4", "--seed", "20260901",
                   "--threads", str(threads), "--json"]
        if checkpoints:
            command += ["--checkpoints", "0,100,500"]
        result = subprocess.run(command + list(extra), capture_output=True, text=True, check=True)
        return json.loads(result.stdout)

    def test_each_feeder_reads_its_own_contiguous_record_once(self):
        result = self.run_simulator(["--feeder-sequence", str(self.record)])
        self.assertEqual(result["sequence_consumed"], result["admissions_by_replication"])
        self.assertNotEqual(result["feeder_start_positions"][0], result["feeder_start_positions"][1])
        for starts, consumed, belt_counts, wraps in zip(result["feeder_start_positions"],
                result["sequence_consumed"], result["sequence_belt_counts"], result["sequence_wraps"]):
            self.assertEqual(consumed[0], len(self.labels))
            self.assertTrue(any(count < len(self.labels) for count in consumed[1:]))
            for start, count, actual, crossings in zip(starts, consumed, belt_counts, wraps):
                positions = [(start + i) % len(self.labels) for i in range(count)]
                self.assertEqual(len(positions), len(set(positions)))
                expected = Counter(self.labels[position] for position in positions)
                self.assertEqual(actual, [expected[belt] for belt in range(1, 5)])
                self.assertEqual(crossings, (start + count - 1) // len(self.labels) if count else 0)

    def test_uniform_warmup_and_paired_offsets(self):
        chronological = self.run_simulator(["--feeder-sequence", str(self.record)], steps=1000)
        shuffled = self.run_simulator(["--feeder-sequence", str(self.record), "--shuffle-sequence"], steps=1000)
        uniform = self.run_simulator(steps=1000)
        self.assertEqual(chronological["warmup_state_hashes"], uniform["warmup_state_hashes"])
        self.assertEqual(shuffled["warmup_state_hashes"], uniform["warmup_state_hashes"])
        self.assertEqual(chronological["feeder_start_positions"], shuffled["feeder_start_positions"])
        # A partial pass tests that cursors did not advance during warm-up.
        for starts, counts, actual in zip(chronological["feeder_start_positions"],
                                         chronological["sequence_consumed"], chronological["sequence_belt_counts"]):
            for start, count, per_belt in zip(starts, counts, actual):
                expected = Counter(self.labels[(start+i) % len(self.labels)] for i in range(count))
                self.assertEqual(per_belt, [expected[belt] for belt in range(1, 5)])

    def test_shuffle_preserves_counts_and_is_thread_independent(self):
        options = ["--feeder-sequence", str(self.record), "--shuffle-sequence"]
        serial = self.run_simulator(options, threads=1)
        parallel = self.run_simulator(options, threads=3)
        serial.pop("threads")
        parallel.pop("threads")
        self.assertEqual(serial, parallel)
        expected = [self.labels.count(belt) for belt in range(1, 5)]
        for per_feeder in serial["sequence_belt_counts"]:
            self.assertEqual(per_feeder[0], expected)

    def test_rejects_reuse_and_defaults_to_one_record_length(self):
        base = [str(self.binary), "-n", "4", "-m", "4", "-R", "2", "-w", "0",
                "--feeder-sequence", str(self.record), "--json"]
        invalid = subprocess.run(base + ["-T", "2049"], capture_output=True, text=True)
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("may not exceed", invalid.stderr)
        valid = json.loads(subprocess.check_output(base, text=True))
        self.assertEqual(valid["steps"], len(self.labels))
        self.assertTrue(all(row[0] == len(self.labels) for row in valid["sequence_consumed"]))

    def test_checkpoints_are_cumulative_measurement_counts(self):
        result = self.run_simulator(["--feeder-sequence", str(self.record)])
        for totals, checkpoints in zip(result["admissions_by_replication"], result["checkpoint_admissions"]):
            self.assertEqual(checkpoints[0], [0]*4)
            for step, counts in zip(result["checkpoints"], checkpoints):
                self.assertEqual(counts[0], step)
                self.assertTrue(all(0 <= count <= min(step, total) for count, total in zip(counts, totals)))
                measured = (sum(totals)-sum(counts)) / (result["steps"]-step)
                self.assertGreaterEqual(measured, 1)
                self.assertLessEqual(measured, 4)

    def test_legacy_uniform_replications_are_unchanged(self):
        command = [str(self.binary), "-n", "4", "-m", "4", "--dual", "-T", "5000",
                   "-w", "500", "-R", "4", "--seed", "20260901", "--threads", "1", "--json"]
        result = json.loads(subprocess.check_output(command, text=True))
        self.assertEqual(result["replications"], [3.03, 3.0402, 3.0346, 3.0372])

    def test_runner_outputs_paired_summaries_and_resumes(self):
        output = self.root / "runner"
        command = [sys.executable, str(HERE / "feeder_sequence.py"), "--sequence", str(self.record),
                   "--binary", str(self.binary), "--out", str(output), "--budgets", "0,10",
                   "--warmup", "500", "--reps", "4", "--discard-steps", "0,100,500", "--threads", "2"]
        subprocess.run(command, capture_output=True, text=True, check=True)
        run = next(output.iterdir())
        rows = json.loads((run / "summary.json").read_text())["rows"]
        self.assertEqual(len(rows), 6)
        self.assertEqual({row["discard_steps"] for row in rows}, {0, 100, 500})
        for row in rows:
            self.assertAlmostEqual(row["chronological_minus_shuffled"],
                                   row["chronological_mean"]-row["shuffled_mean"])
            self.assertGreaterEqual(row["paired_uniform_halfwidth"], 0)
        originals = {path: path.stat().st_mtime_ns for path in (run / "raw").glob("*.json")}
        self.assertEqual(len(originals), 6)
        subprocess.run(command, capture_output=True, text=True, check=True)
        self.assertEqual(originals, {path: path.stat().st_mtime_ns for path in originals})


if __name__ == "__main__":
    unittest.main()
