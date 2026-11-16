"""Guard repeated measurements against duplicate evidence and censored times."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from proteinwork.comparison import compare_runs
from proteinwork.validation import save_json


class ComparisonTests(unittest.TestCase):
    def make_run(self, root, name, *, seconds=10, passed=True, seeds=(1, 2, 3, 4), draws=512):
        directory = root / name
        directory.mkdir()
        save_json(directory / "manifest.json", {
            "fixture_id": name, "method": "direct", "seeds": list(seeds),
            "sampling_parameters": {"draws": draws, "warmup": 256},
            "target_sha256": "target", "evaluation_sha256": "evaluation", "code_sha256": "code",
            "comparison_context_sha256": "context", "sampler_metadata": {"backend": "jax"},
            "elapsed_seconds": seconds, "raw_observations": 4 * draws, "exported_candidates": 2,
            "first_time_to_accuracy_established": passed, "budget_exceeded": not passed})
        save_json(directory / "diagnostics.json", {"status": "accuracy_passed" if passed else "accuracy_not_reached"})
        return directory

    @patch("proteinwork.comparison.verify_run")
    def test_median_excludes_censored_runs_and_retains_failure(self, verify):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            runs = [self.make_run(root, "a", seconds=10),
                    self.make_run(root, "b", seconds=20, seeds=(5, 6, 7, 8)),
                    self.make_run(root, "c", seconds=900, passed=False, seeds=(9, 10, 11, 12))]
            failed = root / "failed"
            failed.mkdir()
            save_json(failed / "failure.json", {"status": "execution_failed"})
            result = compare_runs([*runs, failed], root / "comparison.json")
            summary = result["repeated_measurements"][0]
            self.assertEqual(summary["median_qualified_seconds"], 15)
            self.assertEqual(summary["unqualified_runs"], 1)
            self.assertTrue(summary["disjoint_seeds_within_method"])
            self.assertEqual(result["failed_runs"], 1)
            self.assertNotIn("observed_time_ratio_to_fastest_direct", result["runs"][2])

    @patch("proteinwork.comparison.verify_run")
    def test_duplicate_evidence_rejected_and_settings_kept_separate(self, verify):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self.make_run(root, "a")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                compare_runs([first, first], root / "duplicate.json")
            second = self.make_run(root, "b", draws=1024)
            result = compare_runs([first, second], root / "comparison.json")
            self.assertEqual(len(result["repeated_measurements"]), 2)

    @patch("proteinwork.comparison.verify_run")
    def test_shared_seeds_are_flagged(self, verify):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = compare_runs([self.make_run(root, "a"), self.make_run(root, "b")], root / "comparison.json")
            self.assertFalse(result["repeated_measurements"][0]["disjoint_seeds_within_method"])


if __name__ == "__main__":
    unittest.main()
