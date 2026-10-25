"""Production arithmetic is checked against independent frozen benchmark evidence."""

import copy
import hashlib
import io
import itertools
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from proteinwork.model import (
    PottsModel, _coordinateMap, benchmark_domains, load_benchmark, load_parameters,
    parse_parameters, restrict_model,
)
from proteinwork.direct import enumerate_target


root = Path(__file__).resolve().parents[1]


def tinyRecords(forbidden="-inf"):
    return ("J 0 1 A A 0.5\nJ 0 1 A B 0.25\nJ 0 1 B A -0.75\nJ 0 1 B B 0\n"
            f"h 0 A 1\nh 0 B {forbidden}\nh 1 A 2\nh 1 B -1\n")


def parseTiny(text=None):
    return parse_parameters(io.BytesIO((text or tinyRecords()).encode("ascii")), positions=2, alphabet="AB")


class ParameterTests(unittest.TestCase):
    def test_selected_fields_preserve_support_without_nan(self):
        model = parseTiny()
        np.testing.assert_array_equal(model.energy(np.array([[0, 0], [0, 1], [1, 0], [1, 1]])),
                                      [-3.5, -0.25, np.inf, np.inf])
        self.assertEqual(model.energy(model.encode("AA")).dtype, np.dtype("float64"))
        self.assertEqual(model.decode(model.encode("AB")), "AB")
        self.assertFalse(model.fields.flags.writeable)
        self.assertFalse(model.couplings.flags.writeable)

    def test_complete_unique_records_are_required(self):
        cases = [tinyRecords() + "h 0 A 1\n", tinyRecords() + "J 0 1 A B 0.25\n",
                 tinyRecords().replace("h 1 B -1\n", ""), tinyRecords().replace("J 0 1 B B 0\n", "")]
        for record in cases:
            with self.subTest(record=record), self.assertRaises(ValueError):
                parseTiny(record)

    def test_only_literal_negative_infinity_fields_are_admitted(self):
        for value in ["nan", "NaN", "inf", "+inf", "-1e999"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                parseTiny(tinyRecords(value))
        for value in ["-inf", "-infinity", "-Infinity"]:
            self.assertTrue(np.isneginf(parseTiny(tinyRecords(value)).fields[0, 1]))
        with self.assertRaises(ValueError):
            parseTiny(tinyRecords().replace("A A 0.5", "A A -inf"))

    def test_malformed_indices_symbols_and_empty_support_fail(self):
        badRecords = [tinyRecords().replace("J 0 1", "J 1 0", 1),
                      tinyRecords().replace("h 0 A 1", "h 2 A 1"),
                      tinyRecords().replace("h 0 A 1", "h 0 C 1"),
                      tinyRecords().replace("h 0 A 1", "h 0 A -inf"),
                      tinyRecords().replace("h 0 A 1", "h 0 A 1 extra")]
        for text in badRecords:
            with self.subTest(text=text), self.assertRaises(ValueError):
                parseTiny(text)

    def test_pinned_standalone_loader_checks_size_and_digest(self):
        data = tinyRecords().encode("ascii")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "parameters.dat"
            path.write_bytes(data)
            options = {"positions": 2, "alphabet": "AB", "expected_sha256": hashlib.sha256(data).hexdigest(), "expected_bytes": len(data)}
            self.assertEqual(load_parameters(path, **options).energy([0, 0]), -3.5)
            for changed in [{"expected_sha256": "0" * 64}, {"expected_bytes": len(data) + 1}]:
                with self.subTest(changed=changed), self.assertRaises(ValueError):
                    load_parameters(path, **(options | changed))

    def test_invalid_state_shapes_and_token_types_are_rejected(self):
        model = parseTiny()
        for states in [[0], [0, 2], [-1, 0], [0.0, 1.0], [False, True], None]:
            with self.subTest(states=states), self.assertRaises(ValueError):
                model.energy(states)
        for sequence in ["A", "AX", None]:
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                model.encode(sequence)

    def test_dense_array_symmetry_and_diagonal_are_validated(self):
        model = parseTiny()
        for location in [(0, 1, 0, 0), (0, 0, 0, 0)]:
            couplings = model.couplings.copy()
            couplings[location] += 1
            with self.assertRaises(ValueError):
                PottsModel("AB", model.fields, couplings)


class RestrictionTests(unittest.TestCase):
    def test_exact_elimination_preserves_all_original_energies(self):
        rng = np.random.default_rng(87)
        fields = rng.normal(size=(5, 3))
        couplings = np.zeros((5, 5, 3, 3))
        for first, second in itertools.combinations(range(5), 2):
            couplings[first, second] = rng.normal(size=(3, 3))
            couplings[second, first] = couplings[first, second].T
        model = PottsModel("ABC", fields, couplings)
        reduced = restrict_model(model, [(2,), (2, 0), (1,), (0, 1, 2), (0,)])
        states = np.array(list(itertools.product(*reduced.domains)))
        self.assertEqual(reduced.original_indices, (1, 3))
        np.testing.assert_allclose(reduced.energy(states), model.energy(reduced.expand(states)), atol=2e-14, rtol=0)
        self.assertTrue(np.isposinf(reduced.energy([1, 0])))

    def test_zero_active_sites_keep_fixed_energy_and_shape(self):
        model = parseTiny()
        reduced = restrict_model(model, [(0,), (1,)])
        states = reduced.random_states(np.random.default_rng(8), 3)
        self.assertEqual(states.shape, (3, 0))
        np.testing.assert_array_equal(reduced.expand(states), [[0, 1]] * 3)
        np.testing.assert_array_equal(reduced.energy(states), [-0.25] * 3)
        self.assertEqual(reduced.energy(np.empty(0, dtype=int)), -0.25)

    def test_invalid_restrictions_fail_instead_of_intersecting(self):
        model = parseTiny()
        for domains in [[(0, 1), (0,)], [(0,), ()], [(0,), (0, 0)], [(0,), (2,)], [(0,), (True,)], [(0,)]]:
            with self.subTest(domains=domains), self.assertRaises(ValueError):
                restrict_model(model, domains)

    def test_coordinate_map_uses_full_reference_positions(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "map.csv"
            path.write_text("model_index_0based,alignment_column_1based,anchor_symbol,uniprot_position_1based\n0,1,-,\n1,2,A,2\n2,3,C,3\n3,4,-,\n4,5,D,4\n")
            rows = _coordinateMap("-AC-D", "MACDE", path)
            self.assertEqual([row["uniprot_position_1based"] for row in rows], [None, 2, 3, None, 4])
            with self.assertRaises(ValueError):
                _coordinateMap("-AC-D", "ACDACD", path)
            path.write_text(path.read_text().replace("1,2,A,2", "1,2,A,1"))
            with self.assertRaises(ValueError):
                _coordinateMap("-AC-D", "MACDE", path)


class RealBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config = json.loads((root / "configs/benchmark.json").read_text())
        if not (root / config["archive"]["local_path"]).is_file():
            raise unittest.SkipTest("Pinned local archive is absent; this suite never downloads it")
        cls.inputs = load_benchmark(root / "configs/benchmark.json")
        cls.checks = json.loads((root / "data/benchmark/score_checks.json").read_text())

    def cloneInputs(self):
        result = copy.copy(self.inputs)
        result.config = copy.deepcopy(result.config)
        result.cases = copy.deepcopy(result.cases)
        result.provenance = copy.deepcopy(result.provenance)
        return result

    def test_every_saved_sequence_matches_independent_decimal_score(self):
        for name, record in self.checks["scores"].items():
            with self.subTest(name=name):
                actual = self.inputs.model.energy(self.inputs.model.encode(record["sequence"]))
                if "energy_decimal" in record:
                    self.assertAlmostEqual(float(actual), float(record["energy_decimal"]), delta=1e-10)
                else:
                    self.assertTrue(np.isposinf(actual))

    def test_all_cases_preserve_original_energy_and_support(self):
        for case in self.inputs.cases:
            with self.subTest(case=case["case_id"]):
                domains = benchmark_domains(self.inputs, case["case_id"])
                reduced = restrict_model(self.inputs.model, domains)
                states = reduced.random_states(np.random.default_rng(56), 12)
                np.testing.assert_allclose(reduced.energy(states), self.inputs.model.energy(reduced.expand(states)), atol=1e-10, rtol=0)
                self.assertEqual(reduced.n_sites, 2 if case["kind"] == "exact_conditional" else 88)

    def test_exact_case_probabilities_match_independent_saved_values(self):
        for caseId, fixture in self.checks["exact_conditional_cases"].items():
            reduced = restrict_model(self.inputs.model, benchmark_domains(self.inputs, caseId))
            states, probabilities = enumerate_target(reduced, fixture["beta"])
            actual = {self.inputs.model.decode(state): probability for state, probability in zip(reduced.expand(states), probabilities)}
            for name, expected in fixture["probabilities_binary64"].items():
                self.assertAlmostEqual(actual[self.checks["scores"][name]["sequence"]], expected, delta=5e-14)

    def test_reference_and_archive_provenance_are_preserved(self):
        self.assertEqual(self.inputs.model.n_sites, 96)
        self.assertEqual(sum(row["uniprot_position_1based"] is not None for row in self.inputs.mapping), 89)
        self.assertEqual(self.inputs.provenance["parameter_sha256"], self.inputs.config["archive"]["member_sha256"])
        self.assertEqual(int(np.isneginf(self.inputs.model.fields).sum()), 362)

    def test_forbidden_gap_duplicate_coordinate_and_semantics_requests_fail(self):
        caseId = self.inputs.cases[0]["case_id"]
        badRestrictions = [[{"model_index_0based": 38, "allowed_symbols": "Q"}],
                           [{"model_index_0based": 0, "allowed_symbols": "A"}],
                           [{"model_index_0based": 8, "allowed_symbols": "-"}],
                           [{"model_index_0based": 8, "allowed_symbols": "A"}] * 2,
                           [{"model_index_0based": True, "allowed_symbols": "A"}],
                           [{"model_index_0based": 8, "allowed_symbols": "A", "uniprot_position_1based": 8}]]
        for restrictions in badRestrictions:
            with self.subTest(restrictions=restrictions), self.assertRaises(ValueError):
                inputs = self.cloneInputs()
                inputs.cases[0]["positional_restrictions"] = restrictions
                benchmark_domains(inputs, caseId)
        inputs = self.cloneInputs()
        inputs.cases[0]["other_positions"] = "unrestricted"
        with self.assertRaises(ValueError):
            benchmark_domains(inputs, caseId)
        inputs = self.cloneInputs()
        inputs.provenance["restriction_data"]["common_constraints"]["nongap_columns_exclude_gap"] = 1
        with self.assertRaises(ValueError):
            benchmark_domains(inputs, caseId)
        with self.assertRaises(ValueError):
            benchmark_domains(self.inputs, "missing")


if __name__ == "__main__":
    unittest.main()
