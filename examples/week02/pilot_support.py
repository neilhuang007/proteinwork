"""Small standard-library helpers for Week 2 evidence, not a sampler."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import itertools
import json
import math
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runtime_record() -> dict:
    """Record every installed distribution as well as the interpreter/host."""
    versions = {
        distribution.metadata["Name"]: distribution.version
        for distribution in importlib.metadata.distributions()
    }
    return {
        "observed_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": dict(sorted(versions.items(), key=lambda item: item[0].lower())),
        "dependency_lock": "requirements/pilots-win-py313.txt",
        "dependency_lock_sha256": sha256(ROOT / "requirements/pilots-win-py313.txt"),
    }


def frequencies(samples: list[str], probabilities: dict[str, float]) -> dict:
    """Compare full returned states to an exact distribution; keep duplicates."""
    if not samples or set(samples) - probabilities.keys():
        raise ValueError("Samples must be nonempty and belong to the enumerated state space.")
    counts = Counter(samples)
    rows = [
        {
            "state": state,
            "count": counts[state],
            "empirical_probability": counts[state] / len(samples),
            "exact_probability": probability,
            "absolute_error": abs(counts[state] / len(samples) - probability),
        }
        for state, probability in probabilities.items()
    ]
    return {
        "n_returned": len(samples),
        "n_unique": len(counts),
        "repeats_retained": True,
        "states": rows,
        "total_variation_distance": 0.5 * math.fsum(row["absolute_error"] for row in rows),
        "max_absolute_error": max(row["absolute_error"] for row in rows),
    }


def exact_ising(biases: tuple[float, float], weight: float, beta: float) -> list[dict]:
    """Enumerate four states with False=-1, True=+1 and each edge counted once."""
    if any(not math.isfinite(value) for value in (*biases, weight, beta)) or beta <= 0:
        raise ValueError("The fixed pilot needs finite parameters and positive beta.")
    rows = []
    for bits in itertools.product((False, True), repeat=2):
        spins = tuple(1 if bit else -1 for bit in bits)
        score = biases[0] * spins[0] + biases[1] * spins[1] + weight * spins[0] * spins[1]
        rows.append({
            "state": "".join("1" if bit else "0" for bit in bits),
            "bits": list(bits),
            "spins": list(spins),
            "score": score,
            "energy": -beta * score,
        })
    shift = min(row["energy"] for row in rows)
    weights = [math.exp(-(row["energy"] - shift)) for row in rows]
    partition = math.fsum(weights)
    for row, value in zip(rows, weights, strict=True):
        row["probability"] = value / partition
    return rows


def save_json(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
