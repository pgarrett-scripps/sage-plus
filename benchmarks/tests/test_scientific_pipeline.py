import base64
import csv
import io
import struct
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_scientific import mzxml_to_mgf
from scientific_entrapment import load_pairs, threshold_fdp


class ScientificPipelineTests(unittest.TestCase):
    def test_mzxml_conversion_preserves_peaks_precursor_charge_and_time(self):
        data = base64.b64encode(struct.pack(">ffff", 100.25, 42, 200.5, 84)).decode()
        xml = f'<mzXML><msRun><scan num="2" msLevel="2" peaksCount="2" retentionTime="PT2M"><precursorMz precursorCharge="3">500.2</precursorMz><peaks precision="32" byteOrder="network">{data}</peaks></scan></msRun></mzXML>'
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "input.mzXML", Path(directory) / "output.mgf"
            source.write_text(xml)
            self.assertEqual(mzxml_to_mgf(source, output), 1)
            text = output.read_text()
            for line in ("PEPMASS=500.2", "CHARGE=3+", "RTINSECONDS=120", "100.25 42", "200.5 84"):
                self.assertIn(line, text)
            source.write_text(xml.replace('peaksCount="2"', 'peaksCount="3"'))
            with self.assertRaisesRegex(ValueError, "peak count"):
                mzxml_to_mgf(source, output)

    def test_entrapment_ties_are_bounded_and_empty_sets_are_undefined(self):
        labels = {"T": "target", "E": "p_target"}
        partners = {"E": "T"}
        rows = [{"sequence": seq, "peptide_q": 0.005, "score": 10} for seq in ("T", "E")]
        result = threshold_fdp(rows, labels, partners, 0.01)
        self.assertEqual(result["paired_score_ties"], 1)
        self.assertEqual((result["paired_fdp_tie_min"], result["paired_fdp_tie_max"]), (0.5, 1.5))
        self.assertIsNone(threshold_fdp([], labels, partners, 0.01)["combined_fdp"])
        self.assertEqual(threshold_fdp(rows[1:], labels, partners, 0.01)["paired_fdp_tie_min"], 2)

    def test_pair_loader_rejects_missing_partners(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pairs.tsv"
            path.write_text("sequence\tpeptide_type\tpeptide_pair_index\nT\ttarget\t1\n")
            with self.assertRaisesRegex(ValueError, "Incomplete"):
                load_pairs(path)


if __name__ == "__main__":
    unittest.main()
