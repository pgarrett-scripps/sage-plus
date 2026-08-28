import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "benchmark.py"
SPEC = importlib.util.spec_from_file_location("sage_benchmark", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
BENCHMARK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCHMARK)


class BenchmarkHarnessTests(unittest.TestCase):
    def test_parse_timing_uses_last_line(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "timing.tsv"
            path.write_text("Command exited with non-zero status 1\n1.25\t2048\t1\n")
            self.assertEqual(
                BENCHMARK.parse_timing(path),
                {"wall_seconds": 1.25, "peak_rss_kib": 2048, "exit_code": 1},
            )

    def test_synthetic_sequences_are_deterministic_and_unique(self):
        sequences = [BENCHMARK.synthetic_sequence(index) for index in range(1000)]
        self.assertEqual(len(sequences), len(set(sequences)))
        self.assertTrue(all(sequence.endswith("K") for sequence in sequences))

    def test_report_aggregates_trials_and_flags_regression(self):
        metadata = {
            "created_at_utc": "2026-08-25T00:00:00+00:00",
            "cpu": "test cpu",
            "logical_cpus": 8,
            "threads": 4,
            "baseline_commit": "baseline",
            "candidate_commit": "candidate",
            "candidate_dirty": False,
            "rustc": "rustc test",
        }
        records = [
            {
                "case": "standard",
                "version": "baseline",
                "trial": 1,
                "wall_seconds": 10.0,
                "peak_rss_kib": 1000,
                "psms_at_one_percent_fdr": 100,
                "peptides_at_one_percent_fdr": 100,
                "artifact_hashes": {},
            },
            {
                "case": "standard",
                "version": "candidate",
                "trial": 1,
                "wall_seconds": 12.0,
                "peak_rss_kib": 900,
                "psms_at_one_percent_fdr": 95,
                "peptides_at_one_percent_fdr": 95,
                "artifact_hashes": {},
            },
        ]
        report = BENCHMARK.render_report(metadata, records)
        self.assertIn("REVIEW: Baseline to candidate wall-time change", report)
        self.assertIn("REVIEW: Baseline to candidate one-percent FDR PSMs", report)

    def test_report_includes_failed_processes(self):
        metadata = {
            "created_at_utc": "2026-08-25T00:00:00+00:00",
            "cpu": "test cpu",
            "logical_cpus": 8,
            "threads": 4,
            "baseline_commit": "baseline",
            "candidate_commit": "candidate",
            "candidate_dirty": False,
            "rustc": "rustc test",
        }
        records = [
            {
                "case": "feature",
                "version": "candidate",
                "trial": "warmup-1",
                "status": "failed",
                "process_return_code": -11,
                "wall_seconds": 6.3,
                "peak_rss_kib": 2048,
            }
        ]
        report = BENCHMARK.render_report(metadata, records)
        self.assertIn("## Failures", report)
        self.assertIn("| feature | candidate | warmup-1 | -11 |", report)


if __name__ == "__main__":
    unittest.main()
