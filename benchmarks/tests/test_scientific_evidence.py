import hashlib
import io
import json
import sys
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analyze_scientific import calibration_uncertainty, compare_public, reference_membership, quant_results, accepted_ms2_keys
from bundle_scientific import verify
from audit_scientific import audit
from scientific_diagnostics import direct_ratio_diagnostic, lfq_threshold_yields


class ScientificEvidenceTests(unittest.TestCase):
    def test_lfq_threshold_counts_distinguish_precursors_and_file_rows(self):
        row = {"peptide": "PEPTIDE", "charge": "2", "q_value": "0.043",
               "is_decoy": "false", "intensity": "100"}
        counts = lfq_threshold_yields([row, {**row, "intensity": "50"}, {**row, "is_decoy": "true"}])
        self.assertEqual(counts["0.01"]["target_precursors"], 0)
        self.assertEqual(counts["0.05"], {"target_precursors": 1, "positive_target_file_rows": 2})

    def test_ratio_diagnostic_is_explicitly_before_lfq_filter_and_requires_ms2(self):
        rows = [{"filename": name, "peptide": "PEPTIDE", "charge": "2", "proteins": "y",
                 "is_decoy": "false", "q_value": "0.9", "intensity": str(intensity)}
                for name, intensity in (("a", 100), ("b", 50))]
        design = {"a": {"condition": "A", "preparation": "Alpha"},
                  "b": {"condition": "B", "preparation": "Alpha"}}
        keys = {(name, "PEPTIDE", "2") for name in ("a", "b")}
        result = direct_ratio_diagnostic(rows, design, {"y": "yeast"}, keys, {"yeast": -1})
        self.assertEqual(result["status"], "posthoc_diagnostic")
        self.assertEqual(result["species"]["yeast"]["ratio_pairs"], 1)
        self.assertAlmostEqual(result["species"]["yeast"]["median_log2_bias"], 0)
        self.assertEqual(rows[0]["q_value"], "0.9")
        result = direct_ratio_diagnostic(rows, design, {"y": "yeast"}, set(), {"yeast": -1})
        self.assertIsNone(result["species"]["yeast"]["median_log2_bias"])

    def test_direct_ms2_evidence_requires_joint_acceptance_and_matching_charge(self):
        base = {"key": ("control.mzML", "scan=1", 2, "PEPTIDE", False),
                "is_decoy": False, "spectrum_q": 0.005, "peptide_q": 0.005}
        keys = accepted_ms2_keys([base])
        self.assertIn(("control.mzML", "PEPTIDE", "2"), keys)
        self.assertIn(("control.mzML", "PEPTIDE", ""), keys)
        self.assertNotIn(("control.mzML", "PEPTIDE", "3"), keys)
        self.assertFalse(accepted_ms2_keys([{**base, "spectrum_q": 0.5}]))
        self.assertFalse(accepted_ms2_keys([{**base, "peptide_q": 0.5}]))
        self.assertFalse(accepted_ms2_keys([{**base, "is_decoy": True}]))

    def test_quantification_ignores_partial_output_from_failed_search(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            job = root / "runs/quantification/mbr-0"
            job.mkdir(parents=True)
            (job / "lfq.parquet").write_bytes(b"partial output")
            (job / "result.json").write_text(json.dumps({"status": "failed"}))
            self.assertEqual(quant_results(root), [])

    def test_reference_variant_preserves_original_and_records_exclusions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "original.fasta", root / "defined.fasta"
            original = ">known\nACDUOK\n>ambiguous\nACDXK\n"
            source.write_text(original)
            subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "prepare_reference_variant.py"),
                            str(source), "--output", str(output)], check=True, capture_output=True)
            self.assertEqual(source.read_text(), original)
            self.assertEqual(output.read_text(), ">known\nACDUOK\n")
            manifest = json.loads(output.with_suffix(".reference.json").read_text())
            self.assertEqual(manifest["retained_proteins"], 1)
            self.assertEqual(manifest["excluded_proteins"][0]["accession"], "ambiguous")

    def test_reference_membership_excludes_isoleucine_leucine_ambiguity(self):
        references = {"human": "AAPELLDEK\x00MQRST", "yeast": "GPELLDEKW"}
        self.assertEqual(reference_membership("PEIIDEK", references), {"human", "yeast"})
        self.assertEqual(reference_membership("KMQ", references), set())

    def test_evidence_audit_detects_modified_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            job = root / "runs/suite/job"
            job.mkdir(parents=True)
            truth = root / "ptm-truth"
            truth.mkdir()
            (truth / "manifest.json").write_text("{}")
            output = job / "result.tsv"
            output.write_bytes(b"original")
            identity = {str(output): hashlib.sha256(b"original").hexdigest()}
            (job / "result.json").write_text(json.dumps({"status": "complete", "outputs": identity}))
            self.assertEqual(audit(root)["status"], "verified")
            output.write_bytes(b"changed")
            result = audit(root)
            self.assertEqual(result["status"], "failed")
            self.assertEqual(result["errors"][0]["file"], str(output))

    def test_bootstrap_preserves_engine_pairing_and_requires_every_cell(self):
        records = []
        for file_id in (0, 1):
            for seed in (20260914, 20260915, 20260916):
                base = 0.01 + file_id * 0.02 + (seed - 20260914) * 0.03
                for engine in ("upstream", "plus"):
                    records.append({"suite": f"entrapment-human-{seed}",
                                    "job": f"file-{file_id}-{engine}",
                                    "thresholds": [{"nominal_q": 0.01,
                                                    "paired_fdp_tie_max": base + (0.005 if engine == "plus" else 0)}]})
        result = calibration_uncertainty(records, replicates=200)[0]
        self.assertEqual(result["status"], "descriptive_pilot_interval")
        for endpoint in result["percentile_intervals"]["paired_difference"]:
            self.assertAlmostEqual(endpoint, 0.005)
        incomplete = calibration_uncertainty(records[:-1], replicates=200)[0]
        self.assertEqual(incomplete["status"], "insufficient_complete_evidence")
        self.assertEqual(len(incomplete["missing_cells"]), 1)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            calibration_uncertainty(records + records[:1], replicates=200)

    def test_public_missing_pair_is_not_zero_discovery(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = {"jobs": [{"id": "study-0-upstream", "engine": "upstream", "config": {}}]}
            (root / "public-comparison-plan.json").write_text(json.dumps(plan))
            result = compare_public(root)[0]
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["engine_status"]["plus"], "missing")
            self.assertNotIn("shared_target_psms", result)

    def test_archive_rejects_changed_or_missing_evidence(self):
        content = b"preserved evidence"
        digest = hashlib.sha256(content).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.tar.gz"
            with tarfile.open(path, "w:gz") as archive:
                info = tarfile.TarInfo("result.json")
                info.size = len(content)
                archive.addfile(info, io.BytesIO(content))
            verify(path, {"result.json": digest})
            with self.assertRaisesRegex(ValueError, "checksum"):
                verify(path, {"result.json": "0" * 64})
            with self.assertRaisesRegex(ValueError, "incomplete"):
                verify(path, {"result.json": digest, "missing.json": digest})


if __name__ == "__main__":
    unittest.main()
