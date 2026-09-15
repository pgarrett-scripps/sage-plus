import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scientific_metrics import (canonical_peptide, compare_identifications,
                                identification_metrics, localization_metrics,
                                normalize_psms, peptide_representatives,
                                quantification_metrics)


def psm(**updates):
    row = dict(rank="1", peptide="PEPTIDEK", label="1", filename="one.mgf",
               scannr="scan=1", charge="2", sage_discriminant_score="8",
               spectrum_q="0.005", peptide_q="0.005", hyperscore="40", proteins="P1")
    row.update(updates)
    return row


class ScientificMetricsTests(unittest.TestCase):
    def test_format_normalization_and_decoy_exclusion(self):
        tsv = normalize_psms([psm(peptide="AC[+57.021500]K"), psm(label="-1", scannr="scan=2")])
        parquet = normalize_psms([psm(peptide="AC[57.0215]K", is_decoy="false"),
                                 psm(is_decoy="true", scannr="scan=2")])
        self.assertTrue(compare_identifications(tsv, parquet)["all_psm_identities_equal"])
        self.assertEqual(identification_metrics(tsv)["0.01"]["target_psms"], 1)
        self.assertNotEqual(canonical_peptide("AS[79.9663]TK"), canonical_peptide("AST[79.9663]K"))

    def test_representative_uses_one_real_row_and_discriminant_score(self):
        rows = normalize_psms([psm(hyperscore="100", sage_discriminant_score="1"),
                               psm(scannr="scan=2", hyperscore="10", sage_discriminant_score="4")])
        self.assertEqual(peptide_representatives(rows)[0]["hyperscore"], 10)
        rows[0]["peptide_q"] = 0.1
        with self.assertRaises(ValueError):
            peptide_representatives(rows)

    def test_missing_and_invalid_values_do_not_become_zero(self):
        with self.assertRaises(ValueError):
            normalize_psms([psm(peptide_q="nan")])
        with self.assertRaises(ValueError):
            normalize_psms([psm(), psm()])
        self.assertIsNone(compare_identifications([], [])["jaccard"])

    def test_known_ratio_quantification_and_missingness(self):
        design = {"a": {"condition": "A", "preparation": "one"},
                  "b": {"condition": "B", "preparation": "one"}}
        def quant(file, peptide, intensity, confirmed="true"):
            return dict(filename=file, peptide=peptide, intensity=str(intensity), proteins="P1",
                        is_decoy="false", q_value="0.001", ms2_confirmed=confirmed)
        rows = [quant("a", "PEPK", 100), quant("b", "PEPK", 400, "false"), quant("a", "OTHERK", 50)]
        result = quantification_metrics(rows, design, {"P1": "ecoli"}, {"ecoli": 2})
        self.assertEqual(result["species"]["ecoli"]["median_log2_bias"], 0)
        self.assertEqual(result["species"]["ecoli"]["missing_fraction_observed_union"], 0.25)
        self.assertEqual(result["quantified_rows_without_ms2_confirmation"], 1)
        self.assertIsNone(result["false_transfer_rate"])

    def test_localization_unknown_truth_is_not_correct(self):
        rows = [dict(spectrum_q="0.001", localization_q_value="0.001", modification_mass="79.966331",
                     filename="a", peptide="AS[79.966331]TK", scannr="1", position="2")]
        self.assertEqual(localization_metrics(rows, {("a", "ASTK"): {2}})["empirical_site_error_fraction"], 0)
        self.assertEqual(localization_metrics(rows, {("a", "ASTK"): {3}})["empirical_site_error_fraction"], 1)
        self.assertIsNone(localization_metrics(rows, {})["empirical_site_error_fraction"])


if __name__ == "__main__":
    unittest.main()
