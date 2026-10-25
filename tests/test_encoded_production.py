"""Independent energy and finite transition checks for encoded protein proposals."""

import importlib.util
import itertools
import unittest
from unittest.mock import patch

import numpy as np

from proteinwork.encoded import (
    EncodingBudgetError, _binarySweepMatrix, _implementedEncoding,
    compile_ising, exact_encoded_transition, sample_encoded,
)
from proteinwork.model import PottsModel, restrict_model


def makeFixture():
    fields = np.array([[0.4, -0.7, 0.9, -np.inf], [-0.2, 0.8, -0.5, 0.3], [0.1, 0.5, -0.4, 0.2]])
    couplings = np.zeros((3, 3, 4, 4))
    couplings[0, 1] = np.array([[0.1, -0.3, 0.7, 0.2], [0.4, 0.5, -0.2, -0.6], [-0.5, 0.2, 0.1, 0.8], [0.2, 0.3, 0.4, 0.5]])
    couplings[0, 2] = np.arange(16).reshape(4, 4) / 50 - 0.1
    couplings[1, 2] = np.arange(16).reshape(4, 4)[::-1] / 40 - 0.2
    for left, right in itertools.combinations(range(3), 2):
        couplings[right, left] = couplings[left, right].T
    original = PottsModel("ACDE", fields, couplings)
    return original, restrict_model(original, ((2, 0, 1), (3, 1, 0, 2), (1,)))


def literalEnergy(original, state):
    score = sum(float(original.fields[site, token]) for site, token in enumerate(state))
    score += sum(float(original.couplings[left, right, state[left], state[right]])
                 for left, right in itertools.combinations(range(len(state)), 2))
    return -score


def stationaryWeights(original, model, states, beta):
    energies = np.array([literalEnergy(original, state) for state in model.expand(states)])
    weights = np.exp(-beta * (energies - energies.min()))
    return weights / weights.sum()


class EncodedCompilerTests(unittest.TestCase):
    def test_one_hot_qubo_and_ising_match_literal_full_protein_energy(self):
        original, model = makeFixture()
        states = np.array(list(itertools.product(*model.domains)))
        expected = np.array([literalEnergy(original, row) for row in model.expand(states)])
        self.assertEqual(model.n_sites, 2)
        for penalty in (0.0, 0.3, 2.0, 7.0):
            encoded = compile_ising(model, penalty)
            codes = encoded.encode(states)
            np.testing.assert_array_equal(encoded.decode(codes), states)
            np.testing.assert_allclose(model.energy(states), expected, atol=2e-14)
            np.testing.assert_allclose(encoded.energy(codes), expected, atol=2e-14)
            pairProducts = codes[:, encoded.pairs[:, 0]] * codes[:, encoded.pairs[:, 1]]
            quboEnergy = encoded.quboConstant + codes @ encoded.quboLinear + pairProducts @ encoded.quboPairs
            np.testing.assert_allclose(quboEnergy, expected, atol=2e-14)
            words = np.array(list(itertools.product((0, 1), repeat=encoded.nSpins)))
            allQubo = encoded.quboConstant + words @ encoded.quboLinear
            allQubo += (words[:, encoded.pairs[:, 0]] * words[:, encoded.pairs[:, 1]]) @ encoded.quboPairs
            np.testing.assert_allclose(encoded.energy(words), allQubo, atol=2e-14)

    def test_forbidden_support_is_removed_and_invalid_words_fail(self):
        _, model = makeFixture()
        encoded = compile_ising(model, penalty=0.0)
        self.assertEqual(encoded.nSpins, 7)
        self.assertNotIn(3, encoded.domains[0])
        with self.assertRaisesRegex(ValueError, "forbidden"):
            encoded.encode(np.array([3, 0]))
        for word in (np.zeros(7), np.ones(7), np.array([1, 0, 0, 0, 0, 0, 0])):
            with self.assertRaisesRegex(ValueError, "invalid encoded"):
                encoded.decode(word)
        with self.assertRaises(ValueError):
            encoded.encode(np.array([0.0, 1.0]))

    def test_resource_budgets_fail_without_discarding_interactions(self):
        _, model = makeFixture()
        encoded = compile_ising(model)
        report = encoded.resource_report()
        self.assertEqual(report["logical_spins"], 7)
        self.assertEqual(report["logical_edges"], 21)
        self.assertEqual(report["logical_max_degree"], 6)
        self.assertFalse(report["physical_execution"])
        for budget in ({"max_spins": 6}, {"max_edges": 20}, {"max_bytes": 100}):
            with self.assertRaises(EncodingBudgetError):
                compile_ising(model, **budget)
        wide = restrict_model(PottsModel("ABCDEFGHIJKLMNOPQR", np.zeros((1, 18)), np.zeros((1, 1, 18, 18))), (tuple(range(18)),))
        wideReport = compile_ising(wide).resource_report()
        self.assertEqual(wideReport["logical_max_degree"], 17)
        self.assertTrue(wideReport["exceeds_published_z1_neighbor_limit"])

    def test_all_fixed_model_has_no_binary_variables(self):
        original, _ = makeFixture()
        fixed = restrict_model(original, ((0,), (2,), (1,)))
        encoded = compile_ising(fixed, max_spins=0, max_edges=0, max_bytes=0)
        self.assertEqual(encoded.nSpins, 0)
        states, transition = exact_encoded_transition(fixed, 1.0)
        self.assertEqual(states.shape, (1, 0))
        np.testing.assert_array_equal(transition, [[1.0]])
        result = sample_encoded(fixed, beta=1.0, seeds=[1, 2], draws=5, warmup=2, backend="thrml")
        self.assertEqual(result.states.shape, (2, 5, 0))
        np.testing.assert_allclose(result.energies, literalEnergy(original, [0, 2, 1]))
        self.assertIn("no sampler call", result.metadata["backend_execution"])


class EncodedTransitionTests(unittest.TestCase):
    def assertTarget(self, original, model, states, transition, beta):
        target = stationaryWeights(original, model, states, beta)
        self.assertGreaterEqual(transition.min(), -1e-15)
        np.testing.assert_allclose(transition.sum(axis=1), 1.0, atol=1e-13, rtol=0)
        np.testing.assert_allclose(target @ transition, target, atol=1e-13, rtol=0)
        flow = target[:, None] * transition
        np.testing.assert_allclose(flow, flow.T, atol=1e-13, rtol=0)
        self.assertTrue(np.all(transition > 0))

    def test_q3_q4_interactions_fixed_sites_penalties_and_scaled_quantized_correction(self):
        original, model = makeFixture()
        for penalty, sweeps, scale, quantization in ((0.0, 1, 1.0, None), (0.3, 2, 1.0, None),
                                                     (2.0, 1, 0.6, 0.125), (4.0, 2, 1.3, 0.2)):
            with self.subTest(penalty=penalty, sweeps=sweeps, scale=scale):
                states, transition = exact_encoded_transition(model, 0.9, penalty, sweeps_per_draw=sweeps,
                    proposal_beta=0.7, coefficient_scale=scale, quantization_step=quantization)
                self.assertTarget(original, model, states, transition, 0.9)

    def test_forward_only_schedule_and_uncorrected_scaling_are_falsified(self):
        original, model = makeFixture()
        beta, proposalBeta = 0.9, 0.7
        encoded = _implementedEncoding(compile_ising(model, 0.3), 0.55, 0.1)
        words = np.array(list(itertools.product((0, 1), repeat=encoded.nSpins)), dtype=bool)
        forward = _binarySweepMatrix(encoded, proposalBeta, words, range(encoded.nSpins))
        reverse = _binarySweepMatrix(encoded, proposalBeta, words, range(encoded.nSpins - 1, -1, -1))
        states, corrected = exact_encoded_transition(model, beta, 0.3, proposal_beta=proposalBeta,
                                                     coefficient_scale=0.55, quantization_step=0.1)
        codes = encoded.encode(states)
        indexes = codes @ (2 ** np.arange(encoded.nSpins - 1, -1, -1))
        target = stationaryWeights(original, model, states, beta)
        uncorrected = ((forward + reverse) / 2)[np.ix_(indexes, indexes)]
        uncorrected[np.diag_indices(len(states))] += 1 - uncorrected.sum(axis=1)
        self.assertGreater(np.max(np.abs(target @ uncorrected - target)), 0.001)
        energies, surrogate = model.energy(states), encoded.energy(codes)
        logRatio = -beta * (energies[None, :] - energies[:, None]) + proposalBeta * (surrogate[None, :] - surrogate[:, None])
        wrongSchedule = forward[np.ix_(indexes, indexes)] * np.exp(np.minimum(0.0, logRatio))
        wrongSchedule[np.diag_indices(len(states))] += 1 - wrongSchedule.sum(axis=1)
        self.assertGreater(np.max(np.abs(target @ wrongSchedule - target)), 0.0001)
        self.assertTarget(original, model, states, corrected, beta)
        self.assertGreater(np.max(np.abs(forward - forward.T)), 0.01)

    def test_full_binary_mixture_reversibility_matches_known_proposal_ratio(self):
        _, model = makeFixture()
        encoded = compile_ising(model, penalty=0.7)
        words = np.array(list(itertools.product((0, 1), repeat=encoded.nSpins)), dtype=bool)
        forward = _binarySweepMatrix(encoded, 0.8, words, range(encoded.nSpins))
        reverse = _binarySweepMatrix(encoded, 0.8, words, range(encoded.nSpins - 1, -1, -1))
        mixture = (forward @ forward + reverse @ reverse) / 2
        energy = encoded.energy(words)
        invariant = np.exp(-0.8 * (energy - energy.min()))
        invariant /= invariant.sum()
        np.testing.assert_allclose(invariant[:, None] * mixture, (invariant[:, None] * mixture).T, atol=1e-14, rtol=0)
        self.assertGreater(np.max(np.abs(mixture - mixture.T)), 0.01)
        codes = encoded.encode(np.array(list(itertools.product(*model.domains))))
        indexes = codes @ (2 ** np.arange(encoded.nSpins - 1, -1, -1))
        self.assertTrue(np.all(mixture[indexes].sum(axis=1) - mixture[np.ix_(indexes, indexes)].sum(axis=1) > 0))


class EncodedSamplingTests(unittest.TestCase):
    def test_forced_invalid_endpoint_is_single_rejection_and_retained_every_time(self):
        original, model = makeFixture()
        initial = np.array([[2, 3], [1, 0]])
        with patch("proteinwork.encoded._numpyEndpoint", return_value=np.zeros(7, dtype=bool)) as endpoint:
            result = sample_encoded(model, beta=1.0, seeds=[41, 42], draws=4, warmup=2, initial_states=initial)
        self.assertEqual(endpoint.call_count, 12)
        self.assertEqual(result.metadata["invalid_endpoint_count"], 8)
        self.assertFalse(result.accepted.any())
        np.testing.assert_array_equal(result.states, np.repeat(initial[:, None, :], 4, axis=1))
        expected = [literalEnergy(original, row) for row in model.expand(initial)]
        np.testing.assert_allclose(result.energies, np.repeat(np.array(expected)[:, None], 4, axis=1))
        self.assertEqual(np.asarray(result.metadata["warmup_states"]).shape, (2, 2, 2))

    def test_mh_rejection_keeps_previous_state_and_energy(self):
        original, model = makeFixture()
        states = np.array(list(itertools.product(*model.domains)))
        energies = model.energy(states)
        initial, endpoint = states[energies.argmin()], states[energies.argmax()]
        encoded = compile_ising(model)
        with patch("proteinwork.encoded._numpyEndpoint", return_value=encoded.encode(endpoint)):
            result = sample_encoded(model, beta=50.0, proposal_beta=0.01, coefficient_scale=0.1,
                                    seeds=[61], draws=4, warmup=0, initial_states=initial[None])
        self.assertEqual(result.metadata["mh_rejection_count"], 4)
        np.testing.assert_array_equal(result.states[0], np.repeat(initial[None], 4, axis=0))
        np.testing.assert_allclose(result.energies, literalEnergy(original, model.expand(initial)))

    def test_numpy_replay_original_energies_and_complete_attempt_counts(self):
        original, model = makeFixture()
        arguments = dict(beta=0.9, seeds=[7, 8], draws=40, warmup=4, sweeps_per_draw=2, penalty=0.5,
                         proposal_beta=0.7, coefficient_scale=0.6, quantization_step=0.125)
        first = sample_encoded(model, **arguments)
        second = sample_encoded(model, **arguments)
        np.testing.assert_array_equal(first.states, second.states)
        np.testing.assert_array_equal(first.accepted, second.accepted)
        expected = np.array([[literalEnergy(original, row) for row in chain] for chain in model.expand(first.states)])
        np.testing.assert_allclose(first.energies, expected, atol=1e-13)
        self.assertEqual(first.states.shape, (2, 40, 2))
        self.assertGreater(first.metadata["invalid_endpoint_count"], 0)
        self.assertGreater(first.metadata["mh_rejection_count"], 0)
        self.assertEqual(first.metadata["binary_updates_per_attempt"], 14)

    @unittest.skipUnless(importlib.util.find_spec("jax"), "JAX unavailable")
    def test_jax_execution_preserves_support_energy_and_replay(self):
        original, model = makeFixture()
        arguments = dict(beta=0.9, seeds=[17, 18], draws=64, warmup=4, backend="jax", penalty=0.3,
                         proposal_beta=0.7, coefficient_scale=0.6, quantization_step=0.125)
        first, second = sample_encoded(model, **arguments), sample_encoded(model, **arguments)
        np.testing.assert_array_equal(first.states, second.states)
        np.testing.assert_array_equal(first.accepted, second.accepted)
        expected = np.array([[literalEnergy(original, row) for row in chain] for chain in model.expand(first.states)])
        np.testing.assert_allclose(first.energies, expected, atol=1e-13)
        self.assertGreater(first.metadata["invalid_endpoint_count"], 0)
        self.assertTrue(first.metadata["jax_enable_x64"])

    @unittest.skipUnless(importlib.util.find_spec("thrml"), "official THRML unavailable")
    def test_real_thrml_protein_execution_uses_audited_finite_schedule(self):
        original, model = makeFixture()
        result = sample_encoded(model, beta=0.9, seeds=[27], draws=8, warmup=2, backend="thrml", penalty=0.3,
                                sweeps_per_draw=2, proposal_beta=0.7, coefficient_scale=0.6, quantization_step=0.125)
        expected = np.array([[literalEnergy(original, row) for row in chain] for chain in model.expand(result.states)])
        np.testing.assert_allclose(result.energies, expected, atol=1e-13)
        self.assertEqual(result.metadata["thrml_schedule"], {"n_warmup": 2, "n_samples": 1, "steps_per_sample": 1})
        self.assertIn("thrml.sample_states", result.metadata["official_api"])
        self.assertEqual(len(result.metadata["source_sha256"]), 3)
        self.assertTrue(result.metadata["jax_devices"])
        self.assertFalse(result.metadata["physical_execution"])
        self.assertGreater(result.timings["compilation_seconds"], 0)

    @unittest.skipUnless(importlib.util.find_spec("thrml"), "official THRML unavailable")
    def test_thrml_zero_edge_model_records_harmless_api_padding(self):
        model = restrict_model(PottsModel("ACD", np.zeros((1, 3)), np.zeros((1, 1, 3, 3))), ((0, 1, 2),))
        result = sample_encoded(model, beta=1.0, seeds=[29], draws=8, warmup=0, backend="thrml", penalty=0.0)
        self.assertTrue(result.metadata["zero_weight_api_padding_edge"])
        self.assertEqual(result.metadata["resources"]["logical_edges"], 0)
        np.testing.assert_array_equal(result.energies, 0.0)

    def test_malformed_controls_and_numerical_saturation_fail(self):
        _, model = makeFixture()
        for changes in ({"beta": 0}, {"draws": 0}, {"warmup": -1}, {"seeds": [1, 1]},
                        {"backend": "hardware"}, {"penalty": -1}, {"quantization_step": 0}, {"sweeps_per_draw": 0}):
            arguments = dict(beta=1.0, seeds=[1], draws=2, warmup=0)
            arguments.update(changes)
            with self.assertRaises(ValueError):
                sample_encoded(model, **arguments)
        with self.assertRaises(FloatingPointError):
            sample_encoded(model, beta=1.0, seeds=[1], draws=2, warmup=0, coefficient_scale=10000)
        with self.assertRaises(EncodingBudgetError):
            exact_encoded_transition(model, 1.0, max_binary_states=64)


if __name__ == "__main__":
    unittest.main()
