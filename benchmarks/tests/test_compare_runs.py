import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_runs import compare


def row(peptide="PEPTIDE", q=0.005, psm_id=1, decoy=False):
    return dict(filename="file", scannr="scan=1", rank=1, peptide=peptide, charge=2,
                is_decoy=decoy, spectrum_q=q, peptide_q=q, psm_id=psm_id, hyperscore=5.0)


class ComparisonTests(unittest.TestCase):
    def test_equal_counts_do_not_hide_identity_turnover(self):
        result = compare([row("PEPTIDE")], [row("ANOTHER")], 0.01)
        self.assertEqual(result["accepted_peptidoforms"]["added"], ["ANOTHER"])
        self.assertEqual(result["accepted_peptidoforms"]["lost"], ["PEPTIDE"])
        self.assertEqual(result["accepted_target_psms"]["shared"], 0)

    def test_psm_ids_are_ignored_and_confidence_changes_are_reported(self):
        result = compare([row(psm_id=1)], [row(psm_id=400, q=0.02)], 0.01)
        self.assertEqual(result["stored_psms"]["shared"], 1)
        self.assertEqual(len(result["accepted_target_psms"]["lost"]), 1)
        self.assertAlmostEqual(result["confidence_shifts"]["spectrum_q"]["max_absolute_change"], 0.015)

    def test_invalid_or_ambiguous_rows_are_rejected(self):
        with self.assertRaises(ValueError):
            compare([row(), row()], [], 0.01)
        with self.assertRaises(ValueError):
            compare([row(q=float("nan"))], [], 0.01)
        result = compare([row(decoy=True)], [row(decoy=True)], 0.01)
        self.assertEqual(result["accepted_target_psms"]["shared"], 0)


if __name__ == "__main__":
    unittest.main()
