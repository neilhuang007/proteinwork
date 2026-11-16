"""Candidate review must preserve trace identity and unknown uncertainty."""

import hashlib
from types import SimpleNamespace
import unittest

import numpy as np

from proteinwork.model import PottsModel
from scripts.build_review_packet import diagnostic_audit, validate_candidate, validate_reference


class ReviewPacketTests(unittest.TestCase):
    def candidate_fixture(self):
        model = PottsModel("AC", np.array([[0.0, 0.2], [0.1, -0.1]]), np.zeros((2, 2, 2, 2)))
        inputs = SimpleNamespace(model=model, anchor="AA", mapping=[
            {"model_index_0based": site, "anchor_symbol": "A", "uniprot_position_1based": site + 7}
            for site in range(2)])
        states = np.array([[[0, 0], [0, 1], [0, 1]]])
        sequence_hash = hashlib.sha256(b"AC").hexdigest()
        row = {"candidate_id": "cm_" + sequence_hash[:16], "sequence_sha256": sequence_hash,
               "aligned_sequence": "AC", "sequence": "AC", "source_chain": "0", "source_draw": "1",
               "occurrences": "2", "energy": str(model.energy(states[0, 1])), "hamming_to_anchor": "1"}
        return inputs, states, row

    def test_candidate_links_to_repeated_trace_and_mapped_position(self):
        inputs, states, row = self.candidate_fixture()
        self.assertEqual(validate_candidate(row, inputs, states, inputs.model.energy(states), ["A", "AC"]), ["A8C"])

    def test_inconsistent_identity_count_score_and_source_are_rejected(self):
        inputs, states, row = self.candidate_fixture()
        for field, value in (("candidate_id", "wrong"), ("occurrences", "1"), ("energy", "100"),
                             ("source_draw", "0"), ("source_chain", "-1"), ("sequence", "AA")):
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_candidate({**row, field: value}, inputs, states, inputs.model.energy(states), ["A", "AC"])
        with self.assertRaises(ValueError):
            validate_candidate(row, inputs, states, inputs.model.energy(states), ["A", "A"])

    def test_zero_observations_are_unknown_uncertainty(self):
        diagnostics = {"status": "accuracy_not_reached", "checks": {"reference": False},
                       "observables": {}, "reference_qualification": {"qualified": False, "checks": {"frequency": False}},
                       "frequencies": [
                           {"key": "rare", "kind": "site", "model_indices": [0], "symbols": "C", "probability": 0.0, "mcse": None},
                           {"key": "noisy", "kind": "site", "model_indices": [1], "symbols": "A", "probability": 0.5, "mcse": 0.01},
                           {"key": "fixed", "kind": "site", "model_indices": [2], "symbols": "A", "probability": 1.0, "mcse": 0.0}]}
        audit, blockers = diagnostic_audit(diagnostics, {"broad_case_reference": {"frequency_mcse_max": 0.0025}})
        self.assertEqual(audit["unknown_frequency_mcse"], 1)
        self.assertEqual(audit["frequency_mcse_above_limit"], 1)
        self.assertEqual(len(blockers), 2)
        self.assertIsNone(blockers[0]["mcse"])

    def test_reference_requires_same_target_and_independent_seeds(self):
        source = {"target_sha256": "target", "evaluation_sha256": "evaluation", "seeds": [1, 2, 3, 4]}
        reference = {**source, "seeds": [5, 6, 7, 8], "purpose": "reference", "method": "direct"}
        validate_reference(source, reference)
        for override in ({"target_sha256": "other"}, {"seeds": [1, 6, 7, 8]}, {"purpose": "candidate_library"}):
            with self.subTest(override=override), self.assertRaises(ValueError):
                validate_reference(source, {**reference, **override})


if __name__ == "__main__":
    unittest.main()
