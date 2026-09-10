import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from provenance import atomic_json, cached_stage


class ProvenanceTests(unittest.TestCase):
    def test_cache_requires_matching_inputs_and_intact_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output, manifest = [root / name for name in ("input", "output", "manifest.json")]
            source.write_text("first")
            for changed, expected_hit in ((False, False), (False, True), (True, False), (False, True)):
                if changed:
                    source.write_text("second")
                with cached_stage(manifest, "test", [source], [output], {"threads": 1}) as hit:
                    self.assertEqual(hit, expected_hit)
                    if not hit:
                        output.write_text(source.read_text())
                self.assertEqual(output.read_text(), source.read_text())
            output.write_text("corrupt")
            with cached_stage(manifest, "test", [source], [output], {"threads": 1}) as hit:
                self.assertFalse(hit)
                output.write_text("repaired")
            with cached_stage(manifest, "test", [source], [output], {"threads": 2}) as hit:
                self.assertFalse(hit)
                output.write_text("different settings")

    def test_missing_inputs_and_failed_stages_cannot_be_cached(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output, manifest = root / "output", root / "manifest.json"
            output.write_text("old output")
            with self.assertRaises(FileNotFoundError):
                with cached_stage(manifest, "test", [root / "missing"], [output], {}):
                    self.fail("missing input accepted")
            with self.assertRaisesRegex(RuntimeError, "failed"):
                with cached_stage(manifest, "test", [], [output], {}):
                    output.write_text("partial")
                    raise RuntimeError("failed")
            self.assertFalse(manifest.exists())
            with cached_stage(manifest, "test", [], [output], {}) as hit:
                self.assertFalse(hit)
                self.assertFalse(output.exists())
                output.write_text("complete")

    def test_atomic_json_preserves_old_manifest_on_serialization_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            atomic_json(path, {"complete": True})
            with self.assertRaises(ValueError):
                atomic_json(path, {"invalid": float("nan")})
            self.assertEqual(json.loads(path.read_text()), {"complete": True})
            self.assertEqual(list(path.parent.iterdir()), [path])

    def test_input_mutation_during_a_stage_prevents_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output, manifest = [root / name for name in ("input", "output", "manifest.json")]
            source.write_text("original")
            with self.assertRaisesRegex(RuntimeError, "inputs changed"):
                with cached_stage(manifest, "test", [source], [output], {}):
                    output.write_text("ambiguous result")
                    source.write_text("changed")
            self.assertFalse(manifest.exists())


if __name__ == "__main__":
    unittest.main()
