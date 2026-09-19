import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analyze_hardening import transfer_q_values


def row(score, decoy=False, file="a", eligible=True):
    return {"file_score": score, "is_decoy": decoy, "filename": file,
            "transfer_candidate": eligible, "intensity": 1}


class TransferConfidenceTests(unittest.TestCase):
    def test_ties_do_not_depend_on_decoy_order(self):
        rows = [row(2), row(1), row(1, True)]
        self.assertEqual(transfer_q_values(rows), {0: 1, 1: 1, 2: 1})
        self.assertEqual(transfer_q_values(list(reversed(rows))), {0: 1, 1: 1, 2: 1})

    def test_each_recipient_and_eligibility_is_separate(self):
        rows = [row(1) for _ in range(200)] + [row(1, file="b"), row(100, eligible=False)]
        estimates = transfer_q_values(rows)
        self.assertEqual(estimates[0], 0.005)
        self.assertEqual(estimates[200], 1)
        self.assertNotIn(201, estimates)

    def test_no_decoys_keeps_pseudocount_and_no_targets_is_conservative(self):
        self.assertEqual(transfer_q_values([row(1)]), {0: 1})
        self.assertEqual(transfer_q_values([row(1, True)]), {0: 1})
        self.assertEqual(transfer_q_values([]), {})

    def test_invalid_scores_are_rejected(self):
        with self.assertRaises(ValueError):
            transfer_q_values([row(float("nan"))])


if __name__ == "__main__":
    unittest.main()
