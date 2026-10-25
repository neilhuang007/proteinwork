"""Scientific-convention and teaching-input checks, with no external packages."""

import json
import math
from pathlib import Path
import unittest

from examples.week01.model_reading import (
    PROTEIN_ALPHABET,
    REFERENCE_ENERGIES,
    boltzmann_weight,
    parse_two_site_parameters,
    run_exercise,
)


FIXTURE_DIR = Path(__file__).resolve().parents[1] / "examples" / "week01" / "fixtures"
FIXTURE = (FIXTURE_DIR / "two_site_params.dat").read_text(encoding="utf-8")


class EnergyConventionTests(unittest.TestCase):
    def setUp(self):
        self.model = parse_two_site_parameters(FIXTURE, alphabet="-AC")

    def test_every_state_matches_literal_arithmetic(self):
        for sequence, expected in REFERENCE_ENERGIES.items():
            with self.subTest(sequence=sequence):
                self.assertAlmostEqual(self.model.energy(sequence), expected)
        self.assertEqual(self.model.energy("AC"), -3.5)
        self.assertEqual(self.model.energy("CA"), 2.5)

    def test_energy_score_sign_and_pair_count(self):
        self.assertEqual(self.model.terms("AC"), (1.0, 0.5, 2.0))
        self.assertEqual(self.model.score("AC"), 3.5)
        self.assertEqual(self.model.energy("AC"), -(1.0 + 0.5 + 2.0))
        self.assertNotEqual(self.model.energy("AC"), -(1.0 + 0.5 + 2 * 2.0))
        self.assertEqual(self.model.coupling_record_count, 4)

    def test_explicit_alphabet_mapping_and_reordering(self):
        self.assertEqual(PROTEIN_ALPHABET, "-ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(self.model.encode("-C"), (0, 2))
        self.assertEqual(self.model.encode("AC"), (1, 2))
        reordered = parse_two_site_parameters(FIXTURE, alphabet="CA-")
        self.assertEqual(reordered.encode("AC"), (1, 0))
        self.assertEqual(reordered.energy("AC"), -3.5)

    def test_omitted_coupling_is_zero_and_zero_field_is_explicit(self):
        self.assertEqual(self.model.terms("A-"), (1.0, 0.0, 0.0))
        self.assertEqual(self.model.energy("--"), 0.0)

    def test_weight_order_and_temperature(self):
        self.assertAlmostEqual(boltzmann_weight(-3.5, 1.0), math.exp(3.5))
        self.assertGreater(boltzmann_weight(-3.5, 1.0), boltzmann_weight(2.5, 1.0))
        self.assertAlmostEqual(boltzmann_weight(-3.5, 2.0), math.exp(1.75))

    def test_recorded_probabilities_cover_all_states(self):
        result = run_exercise()
        self.assertEqual(len(result["states"]), 9)
        self.assertEqual({row["sequence"] for row in result["states"]}, set(REFERENCE_ENERGIES))
        self.assertAlmostEqual(sum(row["probability"] for row in result["states"]), 1.0)
        best = max(result["states"], key=lambda row: row["probability"])
        self.assertEqual(best["sequence"], "AC")
        self.assertEqual(len(result["fixture"]["sha256"]), 64)
        metadata = json.loads((FIXTURE_DIR / "two_site_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["symbol_to_index"], {"-": 0, "A": 1, "C": 2})


class TeachingInputTests(unittest.TestCase):
    def test_invalid_records_are_rejected_with_line_number(self):
        cases = [
            ("J 1 0 C A 2.0", "positions"),  # Mirrored triangle would miscount a pair.
            ("J 0 0 A A 2.0", "positions"),  # No self-interaction in this exercise.
            ("J 0 2 A C 2.0", "positions"),  # Only two sites are supported.
            ("h 2 A 1.0", "zero-based"),
            ("h -1 A 1.0", "zero-based"),
            ("h first A 1.0", "zero-based"),
            ("J 0 1 1 2 2.0", "unknown symbol"),  # Numeric legacy format is unsupported.
            ("h 0 D 1.0", "unknown symbol"),
            ("h 0 AC 1.0", "unknown symbol"),
            ("J 0 1 AC A 1.0", "unknown symbol"),
            ("h 0 A", "columns"),
            ("J 0 1 A C 1.0 extra", "columns"),
            ("bias 0 A 1.0", "expected"),
            ("h 0 A nan", "finite"),
            ("h 0 A inf", "finite"),
            ("J 0 1 A C -inf", "finite"),
            ("J 0 1 A C word", "finite"),
            ("h 0 A 1.0", "duplicate h"),
            ("J 0 1 A C 2.0", "duplicate J"),
        ]
        for record, message in cases:
            with self.subTest(record=record):
                with self.assertRaisesRegex(ValueError, f"Line [0-9]+:.*{message}"):
                    parse_two_site_parameters(FIXTURE + "\n" + record, alphabet="-AC")

    def test_missing_fields_and_empty_input_are_rejected(self):
        for text in ("", FIXTURE.replace("h 0 - 0.0\n", "")):
            with self.subTest(text=text):
                with self.assertRaisesRegex(ValueError, "requires all h"):
                    parse_two_site_parameters(text, alphabet="-AC")

    def test_alphabet_must_be_explicit_unique_supported_symbols(self):
        for alphabet in ("", "A", "-AAC", "A C", "auto", "-AU", None):
            with self.subTest(alphabet=alphabet):
                with self.assertRaisesRegex(ValueError, "explicit alphabet"):
                    parse_two_site_parameters(FIXTURE, alphabet=alphabet)

    def test_sequence_length_symbols_and_case_are_strict(self):
        model = parse_two_site_parameters(FIXTURE, alphabet="-AC")
        for sequence in ("", "A", "ACA", "AD", "ac", "A."):
            with self.subTest(sequence=sequence):
                with self.assertRaises(ValueError):
                    model.energy(sequence)

    def test_nonfinite_values_and_nonpositive_temperature_are_rejected(self):
        for temperature in (0.0, -1.0, math.inf, math.nan):
            with self.subTest(temperature=temperature):
                with self.assertRaisesRegex(ValueError, "Temperature"):
                    boltzmann_weight(0.0, temperature)
        with self.assertRaisesRegex(ValueError, "Energy"):
            boltzmann_weight(math.inf, 1.0)
        for energy in (-1e308, 1e308):
            with self.subTest(energy=energy):
                with self.assertRaisesRegex(ValueError, "float range"):
                    boltzmann_weight(energy, 1.0)
        with self.assertRaisesRegex(ValueError, "float range"):
            boltzmann_weight(-1.0, 1e-320)

    def test_blank_lines_and_full_line_comments_are_accepted(self):
        model = parse_two_site_parameters("\n# fixture\n" + FIXTURE + "\n", alphabet="-AC")
        self.assertEqual(model.energy("AC"), -3.5)


if __name__ == "__main__":
    unittest.main()
