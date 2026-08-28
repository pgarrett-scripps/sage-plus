import importlib.util
import sys
import unittest
from pathlib import Path


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCHMARK_ROOT))
MODULE_PATH = BENCHMARK_ROOT / "charge_benchmark.py"
SPEC = importlib.util.spec_from_file_location("sage_charge_benchmark", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
CHARGE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHARGE)


class ChargeBenchmarkTests(unittest.TestCase):
    def test_parse_output_reads_typed_key_values(self):
        parsed = CHARGE.parse_output(
            "fragment_charge_array=true\npreprocess_ns_per_spectrum=7332.1\npsms=3300\n"
        )
        self.assertEqual(
            parsed,
            {
                "fragment_charge_array": True,
                "preprocess_ns_per_spectrum": 7332.1,
                "psms": 3300,
            },
        )


if __name__ == "__main__":
    unittest.main()
