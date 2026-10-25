"""Bounded Week 3 arithmetic exercises; no production sampler or model training."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from examples.week01.model_reading import parse_two_site_parameters  # noqa: E402

EXAMPLE_DIR = Path(__file__).resolve().parent
ENERGIES = (0.0, -1.0, 0.5)
ASYMMETRIC_Q = ((0.0, 0.8, 0.2), (0.4, 0.0, 0.6), (0.3, 0.7, 0.0))


def encode_wall(category: int, domain_size: int) -> tuple[int, ...]:
    if type(domain_size) is not int or domain_size < 1:
        raise ValueError("domain_size must be a positive integer")
    if type(category) is not int or not 0 <= category < domain_size:
        raise ValueError("category must index the stated domain")
    return tuple(int(category > position) for position in range(domain_size - 1))


def decode_wall(word: tuple[int, ...]) -> int:
    if any(type(bit) is not int or bit not in (0, 1) for bit in word):
        raise ValueError("wall word must contain integer bits")
    if any(left < right for left, right in zip(word, word[1:])):
        raise ValueError("reverse domain wall: invalid encoded endpoint")
    return sum(word)


def mh_acceptance(
    energy_before: float,
    energy_after: float,
    q_forward: float,
    q_reverse: float,
    *,
    temperature: float = 1.0,
) -> float:
    """General finite-energy MH ratio, including an impossible reverse move."""
    if not all(math.isfinite(value) for value in (energy_before, energy_after, temperature)):
        raise ValueError("energies and temperature must be finite")
    if temperature <= 0.0:
        raise ValueError("temperature must be positive")
    if not math.isfinite(q_forward) or not 0.0 < q_forward <= 1.0:
        raise ValueError("forward proposal probability must lie in (0, 1]")
    if not math.isfinite(q_reverse) or not 0.0 <= q_reverse <= 1.0:
        raise ValueError("reverse proposal probability must lie in [0, 1]")
    if q_reverse == 0.0:
        return 0.0
    log_ratio = -(energy_after - energy_before) / temperature + math.log(q_reverse) - math.log(q_forward)
    return math.exp(min(0.0, log_ratio))


def exact_target(energies: tuple[float, ...]) -> tuple[float, ...]:
    weights = tuple(math.exp(-(energy - min(energies))) for energy in energies)
    return tuple(weight / math.fsum(weights) for weight in weights)


def mh_matrix(energies, proposal, *, use_proposal_ratio=True):
    """Enumerate only the small, explicit matrix supplied by this exercise."""
    size = len(energies)
    if len(proposal) != size or any(len(row) != size for row in proposal):
        raise ValueError("proposal matrix shape must match energies")
    if any(
        any(not math.isfinite(p) or p < 0.0 for p in row)
        or not math.isclose(math.fsum(row), 1.0, abs_tol=1e-14, rel_tol=0.0)
        for row in proposal
    ):
        raise ValueError("proposal rows must be nonnegative and sum to one")
    matrix = []
    for source, row in enumerate(proposal):
        transition = [0.0] * size
        for destination, forward in enumerate(row):
            if source == destination or forward == 0.0:
                continue
            acceptance = mh_acceptance(
                energies[source], energies[destination],
                forward if use_proposal_ratio else 1.0,
                proposal[destination][source] if use_proposal_ratio else 1.0,
            )
            transition[destination] = forward * acceptance
        transition[source] = 1.0 - math.fsum(transition)
        matrix.append(tuple(transition))
    return tuple(matrix)


def matrix_checks(matrix, target):
    size = len(target)
    stationary = tuple(math.fsum(target[i] * matrix[i][j] for i in range(size)) for j in range(size))
    return {
        "maximum_row_sum_error": max(abs(math.fsum(row) - 1.0) for row in matrix),
        "minimum_entry": min(value for row in matrix for value in row),
        "maximum_stationarity_error": max(abs(a - b) for a, b in zip(stationary, target)),
        "maximum_detailed_balance_error": max(
            abs(target[i] * matrix[i][j] - target[j] * matrix[j][i])
            for i in range(size) for j in range(size)
        ),
    }


def wall_proposal_matrix(domain_size: int):
    """Flip one uniformly chosen bit; invalid words contribute a self proposal."""
    if domain_size < 2:
        raise ValueError("this flip exercise requires at least two categories")
    matrix = []
    for category in range(domain_size):
        word = encode_wall(category, domain_size)
        row = [0.0] * domain_size
        for position in range(len(word)):
            endpoint = list(word)
            endpoint[position] = 1 - endpoint[position]
            try:
                destination = decode_wall(tuple(endpoint))
            except ValueError:
                destination = category
            row[destination] += 1.0 / len(word)
        matrix.append(tuple(row))
    return tuple(matrix)


def scripted_wall_trace():
    """Four prescribed attempted moves, including MH and invalid-word rejection."""
    word = encode_wall(1, 3)
    raw_states = [list(word)]
    events = []
    for index, (position, uniform) in enumerate(((0, 0.9), (1, 0.01), (0, 0.5), (1, 0.5)), start=1):
        before = word
        endpoint = list(word)
        endpoint[position] = 1 - endpoint[position]
        try:
            category = decode_wall(tuple(endpoint))
            acceptance = mh_acceptance(ENERGIES[decode_wall(word)], ENERGIES[category], 0.5, 0.5)
            reason = "accepted" if uniform < acceptance else "mh_rejection"
        except ValueError:
            acceptance = 0.0
            reason = "invalid_endpoint"
        accepted = uniform < acceptance
        if accepted:
            word = tuple(endpoint)
        events.append({
            "attempt": index, "before": list(before), "proposed": endpoint,
            "uniform": uniform, "acceptance": acceptance, "accepted": accepted,
            "reason": reason, "after": list(word),
        })
        raw_states.append(list(word))
    return {"randomization": "prescribed moves and uniforms; no random seed", "states": raw_states, "events": events}


def hand_score_check():
    fixture = REPO_ROOT / "examples/week01/fixtures/two_site_params.dat"
    contents = fixture.read_bytes()
    model = parse_two_site_parameters(contents.decode("utf-8"), alphabet="-AC")
    # Literal arithmetic and mapping, written without calling the model methods.
    hand_terms = (1.0, 0.5, 2.0)
    hand_score = 1.0 + 0.5 + 2.0
    position_map = [
        {"model_index_zero_based": 0, "teaching_reference_position_one_based": 1, "residue": "A"},
        {"model_index_zero_based": 1, "teaching_reference_position_one_based": 2, "residue": "C"},
    ]
    if model.terms("AC") != hand_terms or model.score("AC") != hand_score or model.energy("AC") != -hand_score:
        raise AssertionError("independent hand arithmetic disagrees with model")
    if model.encode("AC") != (1, 2):
        raise AssertionError("independent alphabet indices disagree with model")
    for entry in position_map:
        if "AC"[entry["model_index_zero_based"]] != entry["residue"]:
            raise AssertionError("teaching coordinate map is inconsistent")
    return {
        "fixture": str(fixture.relative_to(REPO_ROOT)), "sha256": hashlib.sha256(contents).hexdigest(),
        "sequence": "AC", "hand_expression": "1.0 + 0.5 + 2.0 = 3.5; E = -3.5",
        "hand_terms": list(hand_terms), "computed_score": model.score("AC"), "computed_energy": model.energy("AC"),
        "alphabet": "-AC", "encoded_alphabet_indices": list(model.encode("AC")), "position_map": position_map,
        "status": "synthetic arithmetic and teaching numbering checked; biological map and both human lead checks pending",
    }


def run_checks():
    target = exact_target(ENERGIES)
    asymmetric = mh_matrix(ENERGIES, ASYMMETRIC_Q)
    incorrect = mh_matrix(ENERGIES, ASYMMETRIC_Q, use_proposal_ratio=False)
    wall_q = wall_proposal_matrix(3)
    wall_kernel = mh_matrix(ENERGIES, wall_q)
    return {
        "exercise": "week03_bounded_learning_checks", "production_sampler": False,
        "hand_score": hand_score_check(),
        "domain_wall_words": [
            {"word": list(word), "valid": word != (0, 1), "category": None if word == (0, 1) else decode_wall(word)}
            for word in itertools.product((0, 1), repeat=2)
        ],
        "energies": list(ENERGIES), "target": list(target),
        "asymmetric_proposal": ASYMMETRIC_Q, "corrected_transition": asymmetric,
        "corrected_checks": matrix_checks(asymmetric, target),
        "incorrect_energy_only_checks": matrix_checks(incorrect, target),
        "wall_proposal": wall_q, "wall_transition": wall_kernel,
        "wall_checks": matrix_checks(wall_kernel, target), "scripted_trace": scripted_wall_trace(),
        "limit": "three valid states, one invalid word; no finite THRML endpoint-ratio claim",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=EXAMPLE_DIR / "output/learning_checks.json")
    args = parser.parse_args()
    result = run_checks()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Hand score AC: {result['hand_score']['computed_score']}; energy: {result['hand_score']['computed_energy']}")
    print(f"Corrected stationarity error: {result['corrected_checks']['maximum_stationarity_error']:.3g}")
    print(f"Wrong energy-only stationarity error: {result['incorrect_energy_only_checks']['maximum_stationarity_error']:.3g}")
    print(f"Wall stationarity error: {result['wall_checks']['maximum_stationarity_error']:.3g}")
    print(f"Raw scripted states: {result['scripted_trace']['states']}")
    print(f"Saved {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
