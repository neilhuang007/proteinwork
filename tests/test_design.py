"""Independent selection checks and held out assay leakage guards."""

import csv
import itertools
import json
import math
from pathlib import Path
import tempfile
import unittest

from proteinwork.design import design_library, validate_design_inputs, validate_design_settings


def pool(sequences, scores=None):
    if scores is None:
        scores = list(range(len(sequences), 0, -1))
    return [{"candidate_id": f"candidate_{i}", "aligned_sequence": sequence,
             "sequence": sequence.replace("-", ""), "energy": -float(score),
             "source_chain": 0, "source_draw": i}
            for i, (sequence, score) in enumerate(zip(sequences, scores))]


def run(candidates, **settings):
    return design_library(candidates, anchor="AAA", alphabet="AB-", settings={"batch_size": 2, **settings})


class DesignSelectionTests(unittest.TestCase):
    def test_milp_finds_valid_batch_where_greedy_gets_stuck(self):
        candidates = pool(["AAA", "AAB", "ABA", "BAA"], [10, 9, 8, 7])
        result = run(candidates, min_pairwise_hamming=2)
        self.assertEqual(result["status"], "completed")
        self.assertEqual([row["aligned_sequence"] for row in result["selected"]], ["AAB", "ABA"])
        self.assertTrue(result["evidence"]["solver"]["proven_optimal"])
        self.assertFalse(result["evidence"]["greedy"]["complete"])
        self.assertEqual(result["evidence"]["selected_objective"], 17.0)
        self.assertFalse(any("selected" in row for row in candidates))
        json.dumps(result, allow_nan=False)

    def test_milp_objective_matches_exhaustive_feasible_batches(self):
        candidates = pool(["AAA", "AAB", "ABA", "BAA", "ABB", "BBB"], [8, 7, 6, 4, 3, 2])
        gamma = 5.25
        result = run(candidates, diversity_weight=gamma, min_pairwise_hamming=1, batch_size=3)
        values = []
        for chosen in itertools.combinations(candidates, 3):
            utility = sum(-row["energy"] for row in chosen)
            similarity = sum(sum(a == b for a, b in zip(left["aligned_sequence"], right["aligned_sequence"])) / 3
                             for left, right in itertools.combinations(chosen, 2))
            values.append(utility - gamma * similarity)
        self.assertAlmostEqual(result["evidence"]["selected_objective"], max(values), places=10)

    def test_qubo_and_ising_equal_literal_objective_for_every_word(self):
        candidates = pool(["AAA", "AAB", "ABA", "BBB"], [4, 3, 2, 1])
        result = run(candidates, diversity_weight=1.25, min_pairwise_hamming=2)
        exported = result["evidence"]["compiled_batch_objective"]
        ordered = {row["candidate_id"]: row for row in candidates}
        ordered = [ordered[key] for key in exported["variable_candidate_ids"]]
        valid_energies = []
        all_energies = []
        for bits in itertools.product((0, 1), repeat=len(ordered)):
            literal = -sum(-row["energy"] * bit for row, bit in zip(ordered, bits))
            conflict_count = 0
            for i in range(len(ordered)):
                for j in range(i + 1, len(ordered)):
                    distance = sum(a != b for a, b in zip(ordered[i]["aligned_sequence"], ordered[j]["aligned_sequence"]))
                    literal += 1.25 * (1 - distance / 3) * bits[i] * bits[j]
                    if distance < 2:
                        conflict_count += bits[i] * bits[j]
            literal += exported["cardinality_penalty"] * (sum(bits) - 2)**2
            literal += exported["conflict_penalty"] * conflict_count
            qubo = exported["qubo_constant"] + sum(a * z for a, z in zip(exported["qubo_linear"], bits))
            qubo += sum(edge["qubo"] * bits[edge["i"]] * bits[edge["j"]] for edge in exported["pair_terms"])
            spins = [2 * z - 1 for z in bits]
            ising = exported["ising_constant"] - sum(b * s for b, s in zip(exported["thrml_biases"], spins))
            ising -= sum(edge["thrml_weight"] * spins[edge["i"]] * spins[edge["j"]] for edge in exported["pair_terms"])
            self.assertAlmostEqual(literal, qubo, places=10)
            self.assertAlmostEqual(literal, ising, places=10)
            all_energies.append(literal)
            if sum(bits) == 2 and conflict_count == 0:
                valid_energies.append(literal)
        self.assertAlmostEqual(min(all_energies), min(valid_energies))

    def test_constraints_and_duplicate_handling_are_explicit(self):
        result = run(pool(["AAA", "ABA", "ABB", "A-B", "ABA", "BBB"]),
                     batch_size=1, min_hamming=1, max_hamming=2, required_motifs=["AB"],
                     forbidden_motifs=["BB"], residue_class="B", min_residue_count=1, max_residue_count=1)
        self.assertEqual(result["status"], "completed")
        exclusions = {row["candidate_id"]: row["reasons"] for row in result["exclusions"]}
        self.assertIn("hamming_to_anchor", exclusions["candidate_0"])
        self.assertIn("forbidden_motif", exclusions["candidate_2"])
        self.assertIn("duplicate_aligned_sequence", exclusions["candidate_4"])
        self.assertEqual(len(result["selected"]), 1)
        for row in result["selected"]:
            self.assertIn("AB", row["sequence"])
            self.assertNotIn("BB", row["sequence"])
            self.assertEqual(row["sequence"].count("B"), 1)
            self.assertTrue(1 <= row["hamming_to_anchor"] <= 2)

    def test_every_selected_candidate_and_pair_satisfies_retained_constraints(self):
        candidates = pool(["AAA", "AAB", "ABA", "BAA", "ABB", "BAB", "BBA", "BBB"])
        result = run(candidates, batch_size=2, min_pairwise_hamming=2, required_motifs=["A"],
                     forbidden_motifs=["BB"], min_hamming=1, max_hamming=2, allowed_residues=["AB"] * 3)
        self.assertEqual(len(result["selected"]), 2)
        for row in result["selected"]:
            self.assertIn("A", row["sequence"])
            self.assertNotIn("BB", row["sequence"])
            self.assertEqual(row["exclusion_reasons"], [])
        for first, second in itertools.combinations(result["selected"], 2):
            self.assertGreaterEqual(sum(a != b for a, b in zip(first["aligned_sequence"], second["aligned_sequence"])), 2)
        with self.assertRaisesRegex(ValueError, "support"):
            run(candidates, allowed_residues=["A", "AB", "AB"])

    def test_insufficient_or_incompatible_pool_does_not_relax_constraints(self):
        result = run(pool(["AAA", "AAB"]), min_hamming=1)
        self.assertEqual(result["status"], "insufficient_candidates")
        self.assertEqual(result["selected"], [])
        result = run(pool(["AAA", "AAB"]), min_pairwise_hamming=2)
        self.assertEqual(result["status"], "constraints_infeasible")
        self.assertEqual(result["selected"], [])

    def test_pool_cap_and_zero_batch_are_declared(self):
        result = run(pool(["AAA", "AAB", "ABA"]), batch_size=1, max_pool=2)
        self.assertEqual(result["evidence"]["eligible_before_pool_limit"], 3)
        self.assertIn("pool_limit", result["exclusions"][-1]["reasons"])
        self.assertEqual(run([], batch_size=0)["selected"], [])

    def test_actual_thrml_optimizer_preserves_feasibility_and_truthful_incumbent(self):
        candidates = pool(["AAA", "AAB", "ABA", "BBB"], [4, 3, 2, 1])
        result = run(candidates, optimizer="thrml", optimizer_seed=729, optimizer_draws=32,
                     optimizer_warmup=8, optimizer_beta=2.0, min_pairwise_hamming=2)
        solver = result["evidence"]["solver"]
        self.assertEqual(solver["backend"], "thrml_simulator_cpu")
        self.assertEqual(solver["total_samples"], 32)
        self.assertEqual(solver["valid_samples"] + solver["invalid_samples"], 32)
        self.assertEqual(len(result["selected"]), 2)
        left, right = result["selected"]
        self.assertGreaterEqual(sum(a != b for a, b in zip(left["aligned_sequence"], right["aligned_sequence"])), 2)
        self.assertFalse(solver["physical_device_used"])
        self.assertFalse(solver["proven_optimal"])
        self.assertFalse(solver["coefficient_clipping"])
        self.assertLessEqual(solver["worst_case_absolute_gibbs_logit"], 30.0)
        self.assertGreaterEqual(result["evidence"]["selected_objective"], result["evidence"]["greedy"]["objective"])
        if solver["selection_source"] == "greedy_incumbent":
            self.assertTrue(solver["fallback_used"])
            self.assertFalse(solver["thrml_improved_greedy"])
        else:
            self.assertGreater(solver["best_feasible_sample_objective"], solver["greedy_objective"])
        encoded = result["evidence"]["compiled_batch_objective"]
        for raw, programmed in zip(encoded["thrml_biases"], solver["programmed_biases"]):
            self.assertAlmostEqual(raw / solver["coefficient_scale_divisor"], programmed)
        json.dumps(result, allow_nan=False)

    def test_thrml_all_invalid_samples_falls_back_with_explicit_source(self):
        from unittest.mock import patch
        import jax.numpy as jnp
        candidates = pool(["AAA", "AAB", "ABA", "BBB"])
        with patch("thrml.sample_states", return_value=[jnp.zeros((4, 4), dtype=jnp.bool_)]):
            result = run(candidates, optimizer="thrml", optimizer_draws=4, optimizer_warmup=0)
        solver = result["evidence"]["solver"]
        self.assertEqual(solver["invalid_fraction"], 1.0)
        self.assertEqual(solver["selection_source"], "greedy_incumbent")
        self.assertTrue(solver["fallback_used"])
        self.assertIsNone(solver["best_feasible_sample_objective"])
        self.assertEqual(len(result["selected"]), 2)
        self.assertTrue(all("greedy_incumbent" in row["selection_reason"] for row in result["selected"]))

    def test_single_candidate_does_not_claim_a_thrml_run(self):
        result = run(pool(["AAA"]), batch_size=1, optimizer="thrml")
        self.assertEqual(len(result["selected"]), 1)
        self.assertEqual(result["evidence"]["solver"]["total_samples"], 0)
        self.assertEqual(result["evidence"]["solver"]["backend"], "deterministic_no_sampling_required")

    def test_unknown_settings_and_nonfinite_candidates_fail_closed(self):
        for settings in ({"mystery": True}, {"batch_size": True}, {"max_seconds": math.inf},
                         {"forbidden_motifs": ["A.*"]}, {"min_residue_count": 1}, {"direction": "minimize"},
                         {"optimizer": "thrml", "max_pool": 65}, {"optimizer_beta": 16}, {"optimizer_draws": 0}):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                run(pool(["AAA", "BBB"]), **settings)
        for energy in (math.inf, -math.inf, math.nan):
            candidates = pool(["AAA", "BBB"])
            candidates[0]["energy"] = energy
            with self.assertRaises(ValueError):
                run(candidates)


class AssayDesignTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / "assays.csv"
        self.settings = {"batch_size": 2, "objective": "assay_ucb", "assay_id": "test_assay", "units": "arbitrary",
                         "context": "synthetic arithmetic fixture", "allowed_residues": ["AB", "AB", "AB"]}
        self.rows = [("AAA", 1, "train", 0), ("AAA", 3, "train", 0), ("BBB", 8, "train", 0),
                     ("AAB", 3, "validation", 1), ("ABB", 6, "test", 1)]
        self.candidates = pool(["AAA", "ABA", "BAA", "BAB", "BBA"])

    def tearDown(self):
        self.directory.cleanup()

    def write(self, rows=None, extra_header=None):
        fields = ["aligned_sequence", "value", "split", "assay_id", "units", "context", "round"]
        if extra_header:
            fields.append(extra_header)
        with self.path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields)
            for sequence, value, split, round_number in self.rows if rows is None else rows:
                writer.writerow([sequence, value, split, "test_assay", "arbitrary", "synthetic arithmetic fixture", round_number])

    def design(self, **overrides):
        return design_library(self.candidates, anchor="AAA", alphabet="AB", assays_path=self.path,
                              settings={**self.settings, **overrides})

    def test_test_and_validation_labels_cannot_change_acquisition(self):
        self.write()
        first = self.design()
        self.write([(sequence, value if split == "train" else -9000 * value, split, round_number)
                    for sequence, value, split, round_number in self.rows])
        second = self.design()
        for left, right in zip(first["all_candidates"], second["all_candidates"]):
            self.assertEqual(left["predicted_mean"], right["predicted_mean"])
            self.assertEqual(left["predictive_std"], right["predictive_std"])
            self.assertEqual(left["utility"], right["utility"])
        self.assertEqual([row["candidate_id"] for row in first["selected"]], [row["candidate_id"] for row in second["selected"]])
        self.assertNotEqual(first["evidence"]["predictor"]["metrics"]["test"]["rmse"],
                            second["evidence"]["predictor"]["metrics"]["test"]["rmse"])
        self.assertEqual(first["evidence"]["assays"]["raw_rows"], 5)
        self.assertEqual(first["evidence"]["predictor"]["metrics"]["train"]["sequences"], 2)
        self.assertEqual(first["all_candidates"][0]["exclusion_reasons"], ["already_measured"])
        json.dumps(first, allow_nan=False)

    def test_minimization_changes_only_mean_sign_in_ucb(self):
        self.write()
        maximum = self.design(direction="maximize", exploration_weight=0.7)
        minimum = self.design(direction="minimize", exploration_weight=0.7)
        for high, low in zip(maximum["all_candidates"], minimum["all_candidates"]):
            self.assertAlmostEqual(high["utility"] + low["utility"], 1.4 * high["predictive_std"])
            self.assertEqual(high["predicted_mean"], low["predicted_mean"])

    def test_replicates_cannot_cross_splits(self):
        self.write(self.rows + [("AAA", 10, "test", 1)])
        with self.assertRaisesRegex(ValueError, "cross"):
            self.design()

    def test_assay_identity_support_and_schema_are_strict(self):
        self.write()
        for overrides in ({"units": "other"}, {"context": "other"}, {"assay_id": "other"},
                          {"allowed_residues": None}):
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                self.design(**overrides)
        self.write(extra_header="invented")
        with self.assertRaises(ValueError):
            self.design()

    def test_assay_controls_outside_prior_support_are_retained_and_annotated(self):
        self.write()
        supported_candidates = pool(["ABA", "AAB"])
        result = design_library(supported_candidates, anchor="AAA", alphabet="AB", assays_path=self.path,
                                settings={**self.settings, "batch_size": 1, "allowed_residues": ["A", "AB", "AB"]})
        self.assertEqual(result["evidence"]["assays"]["outside_model_support_rows"], 1)
        self.assertEqual(result["evidence"]["assays"]["outside_model_support_aligned_sequences"], ["BBB"])
        self.assertIn("BBB", result["evidence"]["predictor"]["training_sequences"])
        self.assertEqual(result["selected"][0]["aligned_sequence"], "ABA")
        with self.assertRaisesRegex(ValueError, "support"):
            self.design(allowed_residues=["A", "AB", "AB"])

    def test_assay_construct_gap_pattern_must_match_anchor(self):
        self.write(self.rows + [("A-A", 4, "test", 1)])
        with self.assertRaisesRegex(ValueError, "gap pattern"):
            design_library([], anchor="AAA", alphabet="AB-", assays_path=self.path,
                           settings={**self.settings, "allowed_residues": ["AB-"] * 3})
        self.write(self.rows + [("ABC", 5, "test", 1)])
        with self.assertRaises(ValueError):
            self.design()

    def test_preflight_checks_inputs_without_model_fitting(self):
        from unittest.mock import patch
        self.write()
        with patch("proteinwork.design.BayesianRidge.fit", side_effect=AssertionError("Must not fit in preflight")):
            result = validate_design_inputs(anchor="AAA", alphabet="AB", settings=self.settings, assays_path=self.path)
        self.assertEqual(result["assays"]["unique_sequences"], 4)
        normalized = validate_design_settings({"batch_size": 1}, "AAA", "AB")
        self.assertEqual(normalized["batch_size"], 1)
        self.assertEqual(validate_design_settings(normalized, "AAA", "AB"), normalized)
        self.write(self.rows + [("AAA", 10, "test", 1)])
        with self.assertRaisesRegex(ValueError, "cross"):
            validate_design_inputs(anchor="AAA", alphabet="AB", settings=self.settings, assays_path=self.path)

    def test_batched_predictions_match_one_call_and_preserve_order(self):
        import numpy as np
        from unittest.mock import patch
        from proteinwork import design
        observations = [{"aligned_sequence": "AAA", "value": 2.0, "split": "train", "replicates": 1},
                        {"aligned_sequence": "BBB", "value": 8.0, "split": "train", "replicates": 1}]
        sequences = ["ABA", "BAA", "BAB", "BBA"] * 151
        recorded = {}
        original_fit = design.BayesianRidge.fit

        def capture(model, features, values):
            recorded["model"] = original_fit(model, features, values)
            return recorded["model"]

        original_features = design._features
        sizes = []

        def bounded_features(items, alphabet, length):
            sizes.append(len(items))
            return original_features(items, alphabet, length)

        with patch.object(design.BayesianRidge, "fit", capture), patch.object(design, "_features", bounded_features):
            means, deviations, evidence = design._fit_predictor(observations, sequences, "AB", 3)
        center = original_features(["AAA", "BBB"], "AB", 3).mean(axis=0)
        expected_mean, expected_std = recorded["model"].predict(original_features(sequences, "AB", 3) - center, return_std=True)
        np.testing.assert_allclose(means, expected_mean * 3 + 5, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(deviations, expected_std * 3, rtol=1e-14, atol=1e-14)
        self.assertLessEqual(max(sizes), 256)
        self.assertEqual(evidence["prediction_chunk_size"], 256)

    def test_training_memory_budget_rejects_before_feature_allocation(self):
        from unittest.mock import patch
        from proteinwork import design
        observations = [{"aligned_sequence": "A" * 2048, "value": 2.0, "split": "train", "replicates": 1},
                        {"aligned_sequence": "B" * 2048, "value": 8.0, "split": "train", "replicates": 1}]
        with patch.object(design, "_features", side_effect=AssertionError("No features should be allocated")):
            with self.assertRaisesRegex(ValueError, "memory budget"):
                design._fit_predictor(observations, [], "AB", 2048)

    def test_nonfinite_missing_or_constant_training_measurements_rejected(self):
        for rows in ([('AAA', math.nan, 'train', 0), ('BBB', 8, 'train', 0)],
                     [('AAA', '', 'train', 0), ('BBB', 8, 'train', 0)],
                     [('AAA', 2, 'train', 0), ('BBB', 2, 'train', 0)],
                     [('AAA', 2, 'validation', 0), ('BBB', 8, 'test', 0)]):
            with self.subTest(rows=rows):
                self.write(rows)
                with self.assertRaises(ValueError):
                    self.design()


if __name__ == "__main__":
    unittest.main()
