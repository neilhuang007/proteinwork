"""Exact kernel checks and independent sampling checks for categorical Gibbs."""

import unittest

import numpy as np

from proteinwork.direct import (
    _localFields, enumerate_target, exact_direct_transition, numpy_gibbs_sweep, sample_direct,
)
from proteinwork.model import PottsModel, restrict_model


def tinyModel():
    fields = np.array([[0.1, -0.3, 0.7], [-0.2, 0.4, 0.0], [0.0, -np.inf, 0.5]])
    couplings = np.zeros((3, 3, 3, 3))
    couplings[0, 1] = [[0.5, -0.1, 0.4], [-0.5, 0.2, -0.1], [0.4, -0.2, 0.1]]
    couplings[1, 0] = couplings[0, 1].T
    couplings[1, 2] = [[0.2, 0.1, -0.3], [0.7, -0.1, 0.4], [0.2, -0.3, -0.2]]
    couplings[2, 1] = couplings[1, 2].T
    return restrict_model(PottsModel("ABC", fields, couplings), [(0, 2), (0, 1, 2), (2,)])


class ExactDirectTests(unittest.TestCase):
    def test_stationarity_and_row_sums_at_multiple_temperatures(self):
        model = tinyModel()
        for beta in [0.2, 1.0, 2.7]:
            states, probability = enumerate_target(model, beta)
            kernelStates, transition = exact_direct_transition(model, beta)
            np.testing.assert_array_equal(states, kernelStates)
            np.testing.assert_allclose(transition.sum(axis=1), 1, atol=2e-15, rtol=0)
            np.testing.assert_allclose(probability @ transition, probability, atol=2e-15, rtol=0)
            self.assertTrue(np.all(transition >= 0))

    def test_sweep_kernel_matches_independent_full_energy_conditionals(self):
        model = tinyModel()
        states, actual = exact_direct_transition(model, 1.3)
        expected = np.eye(len(states))
        for site in range(model.n_sites):
            kernel = np.zeros_like(actual)
            for index, state in enumerate(states):
                matches = np.all(np.delete(states, site, axis=1) == np.delete(state, site), axis=1)
                energies = model.source.energy(model.expand(states[matches]))
                weights = np.exp(-1.3 * (energies - min(energies)))
                kernel[index, matches] = weights / weights.sum()
            expected = expected @ kernel
        np.testing.assert_allclose(actual, expected, atol=2e-15, rtol=0)

    def test_incremental_field_initialization_matches_scalar_terms(self):
        model = tinyModel()
        states = model.random_states(np.random.default_rng(8), 5)
        local = _localFields(model, states)
        for chain, state in enumerate(states):
            for site in range(model.n_sites):
                expected = model.fields[site].copy()
                for other in range(model.n_sites):
                    expected += model.couplings[site, other, :, state[other]]
                np.testing.assert_array_equal(local[chain, site], expected)

    def test_enumeration_limits_and_beta_are_validated(self):
        model = tinyModel()
        with self.assertRaises(ValueError):
            enumerate_target(model, 1, max_states=5)
        for beta in [0, -1, np.inf, np.nan, True]:
            with self.subTest(beta=beta), self.assertRaises(ValueError):
                exact_direct_transition(model, beta)

    def test_numpy_reference_matches_exact_endpoint_probabilities(self):
        model = tinyModel()
        states, kernel = exact_direct_transition(model, 1)
        initial = np.broadcast_to(states[0], (6000, model.n_sites))
        actual = numpy_gibbs_sweep(model, initial, beta=1, rng=np.random.default_rng(328))
        frequencies = np.array([np.mean(np.all(actual == state, axis=1)) for state in states])
        self.assertLess(float(np.max(np.abs(frequencies - kernel[0]))), 0.025)


class DirectSamplingTests(unittest.TestCase):
    def test_compiled_sampler_frequencies_and_repeated_observations(self):
        model = tinyModel()
        result = sample_direct(model, beta=1.3, seeds=[25, 26, 27, 28], draws=3500, warmup=150)
        states, probabilities = enumerate_target(model, 1.3)
        actual = result.states.reshape(-1, model.n_sites)
        frequencies = np.array([np.mean(np.all(actual == state, axis=1)) for state in states])
        self.assertLess(0.5 * float(np.abs(frequencies - probabilities).sum()), 0.035)
        np.testing.assert_array_equal(result.energies, model.energy(result.states))
        np.testing.assert_array_equal(result.accepted[:, 1:], np.any(result.states[:, 1:] != result.states[:, :-1], axis=-1))
        self.assertTrue(np.any(~result.accepted))
        self.assertTrue(result.metadata["float64"])
        self.assertEqual(result.states.shape, (4, 3500, 2))

    def test_same_seed_replays_and_chain_batching_does_not_change_stream(self):
        model = tinyModel()
        options = {"beta": 0.7, "draws": 50, "warmup": 8, "sweeps_per_draw": 2}
        first = sample_direct(model, seeds=[71, 72], **options)
        second = sample_direct(model, seeds=[71, 72], **options)
        differentBatch = sample_direct(model, seeds=[71, 73, 74], **options)
        np.testing.assert_array_equal(first.states, second.states)
        np.testing.assert_array_equal(first.states[0], differentBatch.states[0])
        np.testing.assert_array_equal(first.initial_states[0], differentBatch.initial_states[0])

    def test_zero_site_sampling_and_enumeration_are_exact(self):
        reduced = tinyModel()
        model = restrict_model(reduced.source, [(0,), (1,), (2,)])
        result = sample_direct(model, beta=1, seeds=[1, 2], draws=6, warmup=3)
        states, probability = enumerate_target(model, 1)
        _, transition = exact_direct_transition(model, 1)
        self.assertEqual(result.states.shape, (2, 6, 0))
        self.assertEqual(states.shape, (1, 0))
        np.testing.assert_array_equal(probability, [1])
        np.testing.assert_array_equal(transition, [[1]])
        np.testing.assert_array_equal(result.energies, np.full((2, 6), model.offset))
        self.assertFalse(np.any(result.accepted))

    def test_seed_controls_and_initial_support_fail_closed(self):
        model = tinyModel()
        options = {"beta": 1, "seeds": [1, 2], "draws": 2, "warmup": 0}
        cases = [{"seeds": []}, {"seeds": [1, 1]}, {"seeds": [-1]}, {"seeds": [2 ** 32]},
                 {"seeds": [True]}, {"draws": 0}, {"draws": 1.5}, {"warmup": -1},
                 {"sweeps_per_draw": 0}, {"initial_states": [[1, 0], [0, 0]]},
                 {"initial_states": [[0, 0]]}, {"initial_states": [[0.0, 0.0], [0.0, 0.0]]}]
        for override in cases:
            with self.subTest(override=override), self.assertRaises(ValueError):
                sample_direct(model, **(options | override))


if __name__ == "__main__":
    unittest.main()
