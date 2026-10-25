"""Week 1: read a synthetic two-site adabmDCA text model using only stdlib.

This is a teaching exercise, not the project's future production model loader.
Run from the repository root: .venv/Scripts/python.exe examples/week01/model_reading.py
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys


PROTEIN_ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"
EXAMPLE_DIR = Path(__file__).resolve().parent

# A second calculation: literal arithmetic written from the fixture's table.
# It neither reads the file nor calls the parser, indexing or energy methods.
REFERENCE_ENERGIES = {
    "--": -(0.0 + 0.0 + 0.0),
    "-A": -(0.0 - 0.5 + 0.0),
    "-C": -(0.0 + 0.5 + 0.0),
    "A-": -(1.0 + 0.0 + 0.0),
    "AA": -(1.0 - 0.5 + 0.25),
    "AC": -(1.0 + 0.5 + 2.0),
    "C-": -(-1.0 + 0.0 + 0.0),
    "CA": -(-1.0 - 0.5 - 1.0),
    "CC": -(-1.0 + 0.5 + 0.75),
}


@dataclass(frozen=True)
class TwoSiteModel:
    """Two fields and one q-by-q coupling table; each site pair occurs once."""

    alphabet: str
    fields: tuple[tuple[float, ...], ...]
    couplings: tuple[tuple[float, ...], ...]
    coupling_record_count: int

    def encode(self, sequence: str) -> tuple[int, int]:
        if not isinstance(sequence, str) or len(sequence) != 2:
            raise ValueError("This teaching model requires a sequence of exactly two symbols.")
        if any(symbol not in self.alphabet for symbol in sequence):
            raise ValueError(f"Sequence {sequence!r} contains a symbol outside {self.alphabet!r}.")
        return self.alphabet.index(sequence[0]), self.alphabet.index(sequence[1])

    def terms(self, sequence: str) -> tuple[float, float, float]:
        a, b = self.encode(sequence)
        return self.fields[0][a], self.fields[1][b], self.couplings[a][b]

    def score(self, sequence: str) -> float:
        """Unnormalized field-plus-coupling score S = -E, not a probability."""
        return math.fsum(self.terms(sequence))

    def energy(self, sequence: str) -> float:
        return -self.score(sequence)


def parse_two_site_parameters(text: str, *, alphabet: str) -> TwoSiteModel:
    """Read a narrow teaching subset of the current symbol-based J/h format.

    Require every h value explicitly, accept only J 0 1 records, and set omitted
    J values to zero. Reject duplicate records rather than guessing semantics.
    The caller supplies the exact alphabet order; it is never inferred.
    """
    if (
        not isinstance(alphabet, str)
        or len(alphabet) < 2
        or len(set(alphabet)) != len(alphabet)
        or any(symbol not in PROTEIN_ALPHABET for symbol in alphabet)
    ):
        raise ValueError("Supply an explicit alphabet of distinct protein symbols, e.g. '-AC'.")

    fields: dict[tuple[int, str], float] = {}
    couplings: dict[tuple[str, str], float] = {}
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        parts = raw_line.split()
        if not parts or parts[0].startswith("#"):
            continue
        label = parts[0]
        if label not in {"h", "J"}:
            raise ValueError(f"Line {line_number}: expected an h or J record.")
        expected_width = 4 if label == "h" else 6
        if len(parts) != expected_width:
            raise ValueError(f"Line {line_number}: {label} record requires {expected_width} columns.")
        if label == "h":
            if parts[1] not in {"0", "1"}:
                raise ValueError(f"Line {line_number}: h position must be 0 or 1 (zero-based).")
            symbols = (parts[2],)
        else:
            if parts[1:3] != ["0", "1"]:
                raise ValueError(f"Line {line_number}: teaching J positions must be '0 1' (i < j).")
            symbols = (parts[3], parts[4])
        if any(len(symbol) != 1 or symbol not in alphabet for symbol in symbols):
            raise ValueError(f"Line {line_number}: unknown symbol; alphabet is {alphabet!r}.")
        try:
            value = float(parts[-1])
        except ValueError as error:
            raise ValueError(f"Line {line_number}: parameter value must be a finite number.") from error
        if not math.isfinite(value):
            raise ValueError(f"Line {line_number}: parameter value must be a finite number.")
        if label == "h":
            key = (int(parts[1]), parts[2])
            if key in fields:
                raise ValueError(f"Line {line_number}: duplicate h record {key!r}.")
            fields[key] = value
        else:
            pair = (parts[3], parts[4])
            if pair in couplings:
                raise ValueError(f"Line {line_number}: duplicate J record {pair!r}.")
            couplings[pair] = value

    missing = [(i, symbol) for i in (0, 1) for symbol in alphabet if (i, symbol) not in fields]
    if missing:
        raise ValueError(f"Teaching input requires all h records, including zeros; missing {missing!r}.")
    return TwoSiteModel(
        alphabet=alphabet,
        fields=tuple(tuple(fields[i, symbol] for symbol in alphabet) for i in (0, 1)),
        couplings=tuple(tuple(couplings.get((a, b), 0.0) for b in alphabet) for a in alphabet),
        coupling_record_count=len(couplings),
    )


def boltzmann_weight(energy: float, temperature: float) -> float:
    """Unnormalized w = exp(-E/T); T is dimensionless in this exercise."""
    if not math.isfinite(energy):
        raise ValueError("Energy must be finite.")
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("Temperature must be finite and strictly positive.")
    try:
        weight = math.exp(-energy / temperature)
    except OverflowError as error:
        raise ValueError("Boltzmann weight exceeds this teaching calculation's float range.") from error
    if not math.isfinite(weight):
        raise ValueError("Boltzmann weight exceeds this teaching calculation's float range.")
    if weight == 0.0:
        raise ValueError("Boltzmann weight underflows this teaching calculation's float range.")
    return weight


def run_exercise() -> dict:
    """Evaluate the fixed fixture; fail if it disagrees with the hand reference."""
    fixture_path = EXAMPLE_DIR / "fixtures" / "two_site_params.dat"
    metadata_path = EXAMPLE_DIR / "fixtures" / "two_site_metadata.json"
    fixture_bytes = fixture_path.read_bytes()
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if (
        metadata["length"] != 2
        or metadata["position_index_base"] != 0
        or metadata["coupling_records"] != "i<j, counted once"
        or metadata["alphabet"] != "-AC"
    ):
        raise ValueError("Metadata must describe the fixed two-site '-AC' teaching fixture.")
    model = parse_two_site_parameters(fixture_bytes.decode("utf-8"), alphabet=metadata["alphabet"])
    temperature = metadata["temperature"]
    rows = []
    for symbols in itertools.product(model.alphabet, repeat=2):
        sequence = "".join(symbols)
        h0, h1, coupling = model.terms(sequence)
        energy = model.energy(sequence)
        reference = REFERENCE_ENERGIES[sequence]
        if not math.isclose(energy, reference, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"Energy for {sequence!r} disagrees with the hand reference.")
        rows.append({
            "sequence": sequence,
            "encoded": list(model.encode(sequence)),
            "field_0": h0,
            "field_1": h1,
            "coupling_01": coupling,
            "score": model.score(sequence),
            "energy": energy,
            "reference_energy": reference,
            "boltzmann_weight": boltzmann_weight(energy, temperature),
        })
    partition_sum = math.fsum(row["boltzmann_weight"] for row in rows)
    for row in rows:
        row["probability"] = row["boltzmann_weight"] / partition_sum
    return {
        "exercise": "week01_synthetic_two_site_model_reading",
        "metadata": metadata,
        "fixture": {
            "path": "examples/week01/fixtures/two_site_params.dat",
            "sha256": hashlib.sha256(fixture_bytes).hexdigest(),
            "field_record_count": 2 * len(model.alphabet),
            "coupling_record_count": model.coupling_record_count,
        },
        "reference_check": {"passed": True, "states_checked": len(rows), "absolute_tolerance": 1e-12},
        "partition_sum": partition_sum,
        "probability_sum": math.fsum(row["probability"] for row in rows),
        "states": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=EXAMPLE_DIR / "output" / "results.json")
    args = parser.parse_args()
    try:
        result = run_exercise()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Exercise failed: {error}", file=sys.stderr)
        return 2
    print(f"Synthetic teaching model: L=2, alphabet=-AC, T={result['metadata']['temperature']}")
    print("sequence  score  energy  weight     probability")
    for row in result["states"]:
        print(
            f"{row['sequence']:>8}  {row['score']:5.2f}  {row['energy']:6.2f}"
            f"  {row['boltzmann_weight']:9.6f}  {row['probability']:.6f}"
        )
    print(f"Hand-reference check: {result['reference_check']['states_checked']}/9 states passed")
    print(f"Partition sum: {result['partition_sum']:.6f}; probability sum: {result['probability_sum']:.6f}")
    print(f"Saved {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
