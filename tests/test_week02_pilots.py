"""Scientific conventions and completeness checks for retained Week 2 evidence."""

from collections import Counter
import json
import math
from pathlib import Path
import unittest

from examples.week01.model_reading import REFERENCE_ENERGIES
from examples.week02.pilot_support import exact_ising, frequencies, sha256

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "examples/week02/output"


class ScientificConventionTests(unittest.TestCase):
    def test_two_spin_energy_counts_the_edge_once_and_maps_false_to_minus_one(self):
        rows = exact_ising((0.25, -0.125), 0.5, 1.0)
        # Independent literal arithmetic from the stated b0, b1 and J.
        expected = {"00": -0.375, "01": 0.875, "10": 0.125, "11": -0.625}
        for row in rows:
            self.assertEqual(row["energy"], expected[row["state"]])
        self.assertEqual(rows[0]["spins"], [-1, -1])
        self.assertEqual(rows[-1]["spins"], [1, 1])
        self.assertAlmostEqual(math.fsum(row["probability"] for row in rows), 1.0)

    def test_energy_probability_ratio_and_inverse_temperature(self):
        cold = {row["state"]: row for row in exact_ising((0.0, 0.0), 0.5, 2.0)}
        # Aligned and opposed energies differ by 2 at beta=2; ratio is exp(2).
        self.assertAlmostEqual(cold["11"]["probability"] / cold["10"]["probability"], math.exp(2.0))
        self.assertEqual(cold["00"]["probability"], cold["11"]["probability"])
        self.assertEqual(cold["01"]["probability"], cold["10"]["probability"])

    def test_frequency_comparison_preserves_duplicates_and_zero_count_states(self):
        result = frequencies(["00", "00", "11"], {"00": 0.5, "01": 0.25, "11": 0.25})
        self.assertEqual([row["count"] for row in result["states"]], [2, 0, 1])
        self.assertEqual(result["n_returned"], 3)
        self.assertEqual(result["n_unique"], 2)
        self.assertAlmostEqual(result["total_variation_distance"], 0.25)


class RetainedArtifactTests(unittest.TestCase):
    def read_result(self, name):
        result = json.loads((OUTPUT / f"{name}.json").read_text(encoding="utf-8"))
        self.assertEqual(result["runtime"]["dependency_lock_sha256"], sha256(ROOT / result["runtime"]["dependency_lock"]))
        return result

    def check_counts(self, result, states):
        counts = Counter(states)
        comparison = result["comparison"]
        self.assertEqual(comparison["n_returned"], len(states))
        self.assertTrue(comparison["repeats_retained"])
        self.assertGreater(len(states), len(counts))
        for row in comparison["states"]:
            self.assertEqual(row["count"], counts[row["state"]])
            self.assertEqual(row["empirical_probability"], counts[row["state"]] / len(states))
        self.assertEqual(sum(row["count"] for row in comparison["states"]), len(states))
        expected_tv = 0.5 * math.fsum(abs(row["empirical_probability"] - row["exact_probability"]) for row in comparison["states"])
        self.assertAlmostEqual(comparison["total_variation_distance"], expected_tv)

    def test_adabmdca_complete_sequences_energies_and_fasta(self):
        result = self.read_result("adabmdca")
        self.assertEqual(len(result["sequences"]), result["configuration"]["n_sequences"])
        self.assertEqual(result["fixture"]["sha256"], sha256(ROOT / result["fixture"]["path"]))
        self.assertEqual(result["model_metadata"]["alphabet"], "-AC")
        for row in result["energy_checks"]:
            self.assertEqual(row["official_energy"], REFERENCE_ENERGIES[row["sequence"]])
        fasta_lines = (OUTPUT / "adabmdca.fasta").read_text(encoding="utf-8").splitlines()
        self.assertEqual(fasta_lines[1::2], result["sequences"])
        self.check_counts(result, result["sequences"])

    def test_thrml_complete_returned_boolean_states_and_configuration(self):
        result = self.read_result("thrml")
        config = result["configuration"]
        self.assertEqual(config["free_blocks"], [[0], [1]])
        self.assertEqual(result["returned_shapes"], [[config["n_samples"], 2]])
        bits = result["returned_blocks"][0]
        self.assertEqual(len(bits), config["n_samples"])
        self.assertTrue(all(len(row) == 2 and all(type(bit) is bool for bit in row) for row in bits))
        self.assertEqual(result["actual_model_arrays"]["biases"], config["biases"])
        self.assertEqual(result["actual_model_arrays"]["weights"], config["weights"])
        self.assertEqual(result["prng_key_data"]["root"], [0, config["seed"]])
        expected = exact_ising(tuple(config["biases"]), config["weights"][0], config["beta"])
        self.assertEqual(result["exact_states"], expected)
        self.check_counts(result, ["".join("1" if bit else "0" for bit in row) for row in bits])


class OfficialToolkitEnergyTests(unittest.TestCase):
    def test_installed_adabmdca_agrees_with_all_hand_reference_energies(self):
        # This test requires the pilot environment: failure is not silently skipped.
        from adabmDCA import load_model

        model = load_model(ROOT / "examples/week01/fixtures/two_site_params.dat", alphabet="-AC", device="cpu", dtype="float32")
        sequences = list(REFERENCE_ENERGIES)
        actual = model.compute_energies(sequences).tolist()
        self.assertEqual(actual, [REFERENCE_ENERGIES[sequence] for sequence in sequences])


if __name__ == "__main__":
    unittest.main()
