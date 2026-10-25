"""Falsifying tests for the bounded Week 3 arithmetic exercises."""

import itertools
import math
import unittest

from examples.week03.learning_checks import (
    ASYMMETRIC_Q, ENERGIES, decode_wall, encode_wall, exact_target, hand_score_check,
    matrix_checks, mh_acceptance, mh_matrix, scripted_wall_trace, wall_proposal_matrix,
)


class Week03LearningChecks(unittest.TestCase):
    def test_score_and_position_map_match_independent_literals(self):
        result = hand_score_check()
        self.assertEqual(result["computed_score"], 3.5)
        self.assertEqual(result["computed_energy"], -3.5)
        self.assertEqual(result["encoded_alphabet_indices"], [1, 2])
        self.assertEqual([row["model_index_zero_based"] for row in result["position_map"]], [0, 1])
        self.assertEqual([row["teaching_reference_position_one_based"] for row in result["position_map"]], [1, 2])

    def test_every_small_word_is_valid_exactly_when_prefix_ones(self):
        for domain_size in range(1, 6):
            valid = {encode_wall(category, domain_size) for category in range(domain_size)}
            for word in itertools.product((0, 1), repeat=domain_size - 1):
                if word in valid:
                    self.assertEqual(encode_wall(decode_wall(word), domain_size), word)
                else:
                    with self.assertRaises(ValueError):
                        decode_wall(word)

    def test_domain_boundary_and_malformed_words_fail(self):
        for category, size in ((-1, 3), (3, 3), (0, 0), (True, 3)):
            with self.assertRaises(ValueError):
                encode_wall(category, size)
        for word in ((0, 1), (True, 0), (1.0, 0), (2, 0)):
            with self.assertRaises(ValueError):
                decode_wall(word)

    def test_asymmetric_ratio_matches_independent_hand_arithmetic(self):
        expected = math.exp(-1.0) * 0.8 / 0.4
        self.assertAlmostEqual(mh_acceptance(-1.0, 0.0, 0.4, 0.8), expected)
        self.assertEqual(mh_acceptance(0.0, -1.0, 0.8, 0.4), 1.0)
        self.assertEqual(mh_acceptance(0.0, -1.0, 0.8, 0.0), 0.0)

    def test_bad_mh_probabilities_and_temperatures_fail(self):
        for forward, reverse in ((0.0, 0.2), (-0.1, 0.2), (1.1, 0.2), (0.5, -0.1), (0.5, math.nan)):
            with self.assertRaises(ValueError):
                mh_acceptance(0.0, 1.0, forward, reverse)
        for temperature in (0.0, -1.0, math.inf):
            with self.assertRaises(ValueError):
                mh_acceptance(0.0, 1.0, 0.5, 0.5, temperature=temperature)

    def test_generic_correction_preserves_target_and_energy_only_rule_fails(self):
        target = exact_target(ENERGIES)
        checks = matrix_checks(mh_matrix(ENERGIES, ASYMMETRIC_Q), target)
        self.assertGreaterEqual(checks["minimum_entry"], 0.0)
        for name in ("maximum_row_sum_error", "maximum_stationarity_error", "maximum_detailed_balance_error"):
            self.assertLess(checks[name], 1e-14)
        wrong = matrix_checks(mh_matrix(ENERGIES, ASYMMETRIC_Q, use_proposal_ratio=False), target)
        self.assertGreater(wrong["maximum_stationarity_error"], 0.01)

    def test_invalid_wall_moves_become_self_mass_and_kernel_is_stationary(self):
        proposal = wall_proposal_matrix(3)
        self.assertEqual(proposal, ((0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)))
        checks = matrix_checks(mh_matrix(ENERGIES, proposal), exact_target(ENERGIES))
        self.assertLess(checks["maximum_stationarity_error"], 1e-14)

    def test_trace_keeps_every_rejection_and_invalid_endpoint(self):
        result = scripted_wall_trace()
        self.assertEqual(result["states"], [[1, 0], [1, 0], [1, 1], [1, 1], [1, 0]])
        self.assertEqual([row["reason"] for row in result["events"]], ["mh_rejection", "accepted", "invalid_endpoint", "accepted"])
        self.assertEqual(len(result["states"]), len(result["events"]) + 1)
        self.assertEqual(result["events"][2]["proposed"], [0, 1])

    def test_bad_proposal_matrix_fails_before_calculation(self):
        for proposal in (((0.0, 0.8), (0.4, 0.0)), ((-0.1, 1.1), (0.0, 1.0))):
            with self.assertRaises(ValueError):
                mh_matrix((0.0, 1.0), proposal)


if __name__ == "__main__":
    unittest.main()
