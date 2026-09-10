import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_fdrbench_validation as runner
from analyze_fdrbench_validation import completed_summaries, LIMITS
from provenance import atomic_json, file_identities


class FdrBenchTests(unittest.TestCase):
    def test_arbitrary_binary_names_and_changed_binary_invalidate_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary, spectra, fasta = [root / name for name in ("custom-engine", "input.mzML", "input.fasta")]
            for path in (binary, spectra, fasta):
                path.write_text("input")
            config = root / "config.json"
            atomic_json(config, {"database": {"fasta": str(fasta)}, "mzml_paths": [str(spectra)]})
            output, timing = root / "output", root / "timing.tsv"

            def run(*args):
                (output / "results.sage.parquet").write_text("result")
                (output / "results.json").write_text("{}")
                (output / "run-summary.json").write_text("{}")
                timing.write_text("1.0\t1024\t0\n")

            with patch.object(runner, "GNU_TIME", binary), patch.object(runner, "PRLIMIT", binary), \
                 patch.object(runner, "require_available_memory"), patch.object(runner, "run_logged", side_effect=run) as invoked:
                call = lambda: runner.run_search(binary, config, output, timing, root / "log", 1, "parquet")
                call()
                call()
                self.assertEqual(invoked.call_count, 1)
                binary.write_text("new binary")
                call()
                self.assertEqual(invoked.call_count, 2)
                binary.unlink()
                with self.assertRaises(FileNotFoundError):
                    call()

    def test_nonzero_timing_status_is_not_a_success(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "timing.tsv"
            path.write_text("1\t100\t137\n")
            with self.assertRaises(RuntimeError):
                runner.timing_values(path)
            for invalid in ("nan\t100\t0\n", "-1\t100\t0\n", "1\t0\t0\n"):
                path.write_text(invalid)
                with self.assertRaises(RuntimeError):
                    runner.timing_values(path)

    def test_complete_zero_discovery_runs_are_distinct_from_missing_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "seeds/1/summary.json"
            values = {"q_points": {str(q): None for q in LIMITS}, "matched_yield": {str(q): None for q in LIMITS}}
            summary = {"seed": 1, "status": "completed", "engines": {name: values for name in ("Sage", "Sage Plus")}}
            atomic_json(path, summary)
            manifest = {"schema_version": 2, "seeds": [1], "seed_count": 1, "summary_hashes": file_identities([path])}
            atomic_json(root / "manifest.json", manifest)
            self.assertEqual(len(completed_summaries(root)), 1)
            path.unlink()
            with self.assertRaises(FileNotFoundError):
                completed_summaries(root)
            summary["engines"]["Sage"]["q_points"].pop(str(LIMITS[0]))
            atomic_json(path, summary)
            manifest["summary_hashes"] = file_identities([path])
            atomic_json(root / "manifest.json", manifest)
            with self.assertRaisesRegex(RuntimeError, "threshold matrix"):
                completed_summaries(root)

    def test_threshold_ties_select_the_lowest_score_at_the_cutoff(self):
        base = {"q_value": "0.01", "score": "10", "n_t": "100", "n_p": "1", "paired_fdp": "0.01", "combined_fdp": "0.02", "lower_bound_fdp": "0.01"}
        later = {**base, "score": "9", "n_t": "110"}
        self.assertEqual(runner.point_at_q([base, later], 0.01)["targets"], 110)
        self.assertIsNone(runner.point_at_q([base], 0.001))


if __name__ == "__main__":
    unittest.main()
