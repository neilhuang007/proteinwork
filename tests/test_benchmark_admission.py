"""Admission checks use stdlib fixtures and committed evidence, with no network."""

import csv
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import unittest

from scripts.inspect_benchmark import (
    coordinate_map, exact_probabilities, inspect_parameters, parse_fasta,
    validate_case, verify_bytes,
)


ROOT = Path(__file__).resolve().parents[1]


def tiny_records(field_a="1", field_b="-inf"):
    return ("J 0 1 A A 0.5\nJ 0 1 A B 0.25\nJ 0 1 B A -0.75\nJ 0 1 B B 0\n"
            f"h 0 A {field_a}\nh 0 B {field_b}\nh 1 A 2\nh 1 B -1\n")


def inspect(text, sequences=None):
    return inspect_parameters(io.BytesIO(text.encode("ascii")), 2, "AB", sequences)


class AdmissionChecks(unittest.TestCase):
    def test_selected_field_scoring_preserves_exact_infinite_support(self):
        result = inspect(tiny_records(), {"finite": "AA", "forbidden": "BA", "other": "AB"})
        self.assertEqual(result["score_checks"]["finite"]["energy_decimal"], "-3.5")
        self.assertEqual(result["score_checks"]["other"]["energy_decimal"], "-0.25")
        self.assertEqual(result["score_checks"]["finite"]["selected_pair_count"], 1)
        self.assertEqual(result["score_checks"]["forbidden"]["energy"], "positive_infinity")
        self.assertEqual(result["support"][0]["allowed_symbols"], "A")
        self.assertFalse(math.isnan(result["score_checks"]["finite"]["energy_binary64"]))

    def test_dense_key_coverage_rejects_duplicates_and_missing_keys(self):
        for bad in [tiny_records() + "h 0 A 1\n", tiny_records().replace("J 0 1 B B 0\n", ""),
                    tiny_records() + "J 0 1 A A 0.5\n"]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                inspect(bad)

    def test_only_literal_negative_infinity_fields_are_admitted(self):
        for bad in [tiny_records(field_b="nan"), tiny_records(field_b="inf"),
                    tiny_records(field_b="-1e999"), tiny_records().replace("A A 0.5", "A A -inf")]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                inspect(bad)

    def test_position_and_symbol_errors_rejected(self):
        for bad in [tiny_records().replace("h 0 A 1", "h 2 A 1"),
                    tiny_records().replace("h 0 A 1", "h 0 C 1"),
                    tiny_records().replace("J 0 1", "J 1 0", 1)]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                inspect(bad)

    def test_checksum_size_and_altered_bytes_rejected(self):
        payload = b"immutable input"
        digest = hashlib.sha256(payload).hexdigest()
        verify_bytes(payload, digest, len(payload))
        with self.assertRaises(ValueError):
            verify_bytes(payload + b"!", digest, len(payload))
        with self.assertRaises(ValueError):
            verify_bytes(payload.replace(b"i", b"I", 1), digest, len(payload))

    def test_exact_reference_map_uses_full_accession_numbering_and_null_gaps(self):
        mapping = coordinate_map("-AC-D", "MACDE")
        self.assertEqual([row["uniprot_position_1based"] for row in mapping], [None, 2, 3, None, 4])
        self.assertEqual([row["model_index_0based"] for row in mapping], list(range(5)))
        for aligned, full in [("AAA", "AAAA"), ("AC", "ACAC"), ("AC", "DE"), ("---", "AC")]:
            with self.subTest(aligned=aligned, full=full), self.assertRaises(ValueError):
                coordinate_map(aligned, full)

    def test_fasta_has_no_sequence_before_header_or_duplicate_ids(self):
        self.assertEqual(parse_fasta(b">a\nAB\n\n>b\nBA\n"), {"a": "AB", "b": "BA"})
        for bad in [b"AB\n>a\nAB\n", b">a\nAB\n>a\nBA\n", b">a\n", b""]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                parse_fasta(bad)

    def test_requested_forbidden_or_conflicting_case_restrictions_rejected(self):
        support = [{"allowed_symbols": "AB"}, {"allowed_symbols": "A"}]
        base = {"kind": "broad_conditional", "positional_restrictions": []}
        self.assertEqual(validate_case(base, support, "AA", "AA", 1), ["AB", "A"])
        for restrictions in [[{"model_index_0based": 1, "allowed_symbols": "B"}],
                             [{"model_index_0based": 0, "allowed_symbols": ""}],
                             [{"model_index_0based": 0, "allowed_symbols": "AA"}],
                             [{"model_index_0based": 0, "allowed_symbols": "A"}] * 2,
                             [{"model_index_0based": -1, "allowed_symbols": "A"}]]:
            with self.subTest(restrictions=restrictions), self.assertRaises(ValueError):
                validate_case({**base, "positional_restrictions": restrictions}, support, "AA", "AA", 1)
        for beta in [0, -1, math.inf, math.nan]:
            with self.subTest(beta=beta), self.assertRaises(ValueError):
                validate_case(base, support, "AA", "AA", beta)

    def test_probabilities_respect_temperature(self):
        scores = {"a": {"energy_decimal": "0"}, "b": {"energy_decimal": "1"}}
        at_one = exact_probabilities(scores, 1)
        at_two = exact_probabilities(scores, 2)
        self.assertAlmostEqual(at_one["b"] / at_one["a"], math.exp(-1))
        self.assertAlmostEqual(at_two["b"] / at_two["a"], math.exp(-2))
        self.assertAlmostEqual(sum(at_two.values()), 1)
        for beta in [0, math.nan]:
            with self.assertRaises(ValueError):
                exact_probabilities(scores, beta)

    def test_frozen_case_support_mapping_and_full_model_evidence(self):
        data = ROOT / "data/benchmark"
        config = json.loads((ROOT / "configs/benchmark.json").read_text(encoding="utf-8"))
        cases = json.loads((data / "restriction_cases.json").read_text(encoding="utf-8"))
        metadata = json.loads((data / "reference_provenance.json").read_text(encoding="utf-8"))
        support = json.loads((data / "support.json").read_text(encoding="utf-8"))["positions"]
        scores = json.loads((data / "score_checks.json").read_text(encoding="utf-8"))
        anchor = metadata["aligned_sequence"]
        witness = cases["computational_witness"]["sequence"]
        mapping = coordinate_map(anchor, metadata["full_accession_sequence"])
        self.assertEqual(sum(row["uniprot_position_1based"] is not None for row in mapping), 89)
        self.assertEqual([i for i, a in enumerate(anchor) if a == "-"], config["coordinate_reference"]["gap_model_indices"])
        self.assertEqual(cases["common_constraints"]["fixed_gap_model_indices"], config["coordinate_reference"]["gap_model_indices"])
        with (data / "position_map.csv").open(encoding="utf-8", newline="") as handle:
            committed_map = list(csv.DictReader(handle))
        self.assertEqual([int(row["uniprot_position_1based"]) if row["uniprot_position_1based"] else None
                          for row in committed_map], [row["uniprot_position_1based"] for row in mapping])
        self.assertEqual(scores["scores"]["native_coordinate_anchor"]["energy"], "positive_infinity")
        self.assertEqual(scores["scores"]["native_coordinate_anchor"]["forbidden_selected_fields"],
                         [{"model_index_0based": 38, "symbol": "Q"}, {"model_index_0based": 45, "symbol": "H"}])
        for case in cases["cases"]:
            domains = validate_case(case, support, anchor, witness, config["model"]["beta"])
            self.assertTrue(all(domains))
            if case["kind"] == "exact_conditional":
                self.assertEqual(math.prod(map(len, domains)), 4)
                expected = {"".join(state) for state in itertools.product(*domains)}
                observed = {score["sequence"] for key, score in scores["scores"].items()
                            if key.startswith(case["case_id"] + "/")}
                self.assertEqual(observed, expected)
                for key, score in scores["scores"].items():
                    if key.startswith(case["case_id"] + "/"):
                        self.assertEqual(score["selected_field_count"], 96)
                        self.assertEqual(score["selected_pair_count"], 4560)
                        self.assertEqual(score["forbidden_selected_fields"], [])


if __name__ == "__main__":
    unittest.main()
