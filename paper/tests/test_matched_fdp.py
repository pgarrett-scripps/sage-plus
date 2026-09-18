"""Check incremental FDP counts against the direct paired estimator."""
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis/scripts"))
from collect_matched_fdp import curve, select, threshold_fdp


class MatchedFdpTests(unittest.TestCase):
    def test_ties_missing_partners_and_late_partner_arrivals(self):
        rng = random.Random(20260918)
        labels = {f"{kind}{i}": kind for i in range(40) for kind in ("target", "p_target")}
        partners = {f"p_target{i}": f"target{i}" for i in range(40)}
        for _ in range(30):
            rows = [dict(sequence=seq, peptide_q=rng.choice((0, .001, .01, .02, 1)),
                         score=rng.choice((0, 1, 2)))
                    for seq in labels if rng.random() > .2]
            for point in curve(rows, labels, partners):
                direct = threshold_fdp(rows, labels, partners, point["nominal_q"])
                self.assertEqual(point["targets"], direct["targets"])
                self.assertEqual(point["entrapments"], direct["entrapments"])
                self.assertEqual(point["paired_fdp"], direct["paired_fdp_tie_max"])

    def test_empty_or_unattainable_ceiling(self):
        self.assertEqual(curve([], {}, {}), [])
        self.assertIsNone(select([]))
        self.assertIsNone(select([dict(targets=5, paired_fdp=.02, nominal_q=.1)]))

    def test_missing_unobserved_partner(self):
        rows = [dict(sequence="E", peptide_q=.01, score=1)]
        self.assertEqual(curve(rows, {"E": "p_target"}, {"E": None})[0]["paired_fdp"], 2)

    def test_complete_q_step_is_never_split(self):
        rows = [dict(sequence="T", peptide_q=.01, score=1),
                dict(sequence="E", peptide_q=.01, score=1)]
        points = curve(rows, {"T": "target", "E": "p_target"}, {"E": "T"})
        self.assertEqual(len(points), 1)
        self.assertEqual(points[0]["paired_fdp"], 1.5)


if __name__ == "__main__":
    unittest.main()
