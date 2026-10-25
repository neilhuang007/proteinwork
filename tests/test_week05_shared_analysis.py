"""Check the biological coordinates and denominators used in the weekly analysis."""

import unittest

from examples.week05.shared_analysis import mapped_changes, pair_comparison, residue_profile


class SharedAnalysisTests(unittest.TestCase):
    def test_gap_denominator_and_repeated_records(self):
        profile = residue_profile(["A-", "A-", "C-", "--"], 0)
        self.assertEqual(profile["gap_records"], 1)
        self.assertEqual(profile["nongap_records"], 3)
        self.assertEqual(profile["frequencies_among_nongap"], {"A": 2 / 3, "C": 1 / 3})
        self.assertEqual(residue_profile(["A-", "C-"], 1)["frequencies_among_nongap"], {})

    def test_accession_coordinates_across_internal_gap(self):
        mapping = [{"model_index_0based": i, "anchor_symbol": residue, "uniprot_position_1based": position}
                   for i, (residue, position) in enumerate(zip("A-E", [7, None, 8]))]
        self.assertEqual(mapped_changes("A-E", "C-F", mapping), ["A7C", "E8F"])
        with self.assertRaises(ValueError):
            mapped_changes("A-E", "CAE", mapping)

    def test_unobserved_states_and_repeats_remain_in_probability_check(self):
        rows, variation = pair_comparison(["AE"] * 8 + ["CE"] * 2, [0, 1],
                                         ["AE", "AF", "CE", "CF"], [0.5, 0.1, 0.3, 0.1])
        self.assertEqual([row["count"] for row in rows], [8, 0, 2, 0])
        self.assertAlmostEqual(variation, 0.3)
        self.assertEqual(sum(row["count"] for row in rows), 10)

    def test_malformed_inputs_fail(self):
        for sequences, index in [([], 0), (["A", "AA"], 0), (["AA"], -1), (["AA"], 2)]:
            with self.subTest(sequences=sequences, index=index), self.assertRaises(ValueError):
                residue_profile(sequences, index)
        with self.assertRaises(ValueError):
            pair_comparison(["CE"], [0, 1], ["AE"], [1.0])


if __name__ == "__main__":
    unittest.main()
