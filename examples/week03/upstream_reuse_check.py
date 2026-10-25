"""Check a pinned, licensed Gibbsiq checkout against the synthetic Week 1 model.

No upstream code is copied or installed. Supply a trusted clean source checkout
at the audited commit. This is a tiny adapter exercise, not a protein loader.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from examples.week01.model_reading import REFERENCE_ENERGIES, parse_two_site_parameters  # noqa: E402

PINNED_COMMIT = "06590b82ba509ecf721d8a7685ed27a319d6b0c1"
PINNED_REPOSITORY = "https://github.com/neilhuang007/Gibbsiq.git"
LICENSE_SHA256 = "791eecf284aeec83d27433325ac514a9ade03633e4947a69f85498bbcd709c0d"


def git_read(source: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["git", "--no-optional-locks", "-C", str(source), *arguments], text=True,
    ).strip()


def run_check(source: Path) -> dict:
    source = source.resolve(strict=True)
    if git_read(source, "rev-parse", "HEAD") != PINNED_COMMIT:
        raise ValueError("Gibbsiq checkout must match the audited commit")
    if git_read(source, "remote", "get-url", "origin") != PINNED_REPOSITORY:
        raise ValueError("Gibbsiq origin must match the verified first-party repository")
    if git_read(source, "status", "--porcelain"):
        raise ValueError("Use a clean Gibbsiq checkout; preserve existing changes")
    license_bytes = (source / "LICENSE").read_bytes()
    if hashlib.sha256(license_bytes).hexdigest() != LICENSE_SHA256:
        raise ValueError("Gibbsiq core license differs from the audited MIT license")
    # The supplied CLI path is a trusted checkout, not an imported artifact manifest.
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(source / "src"))
    import gibbsiq
    from gibbsiq.categorical import CategoricalModel
    from gibbsiq.domain_wall import compile_domain_wall

    if Path(gibbsiq.__file__).resolve() != source / "src/gibbsiq/__init__.py":
        raise ValueError("Another gibbsiq import shadows the audited checkout")
    fixture = REPO_ROOT / "examples/week01/fixtures/two_site_params.dat"
    fixture_bytes = fixture.read_bytes()
    model = parse_two_site_parameters(fixture_bytes.decode("utf-8"), alphabet="-AC")
    categorical = CategoricalModel(
        variables=(0, 1), domains={0: tuple(model.alphabet), 1: tuple(model.alphabet)},
        unary={position: {symbol: -model.fields[position][index] for index, symbol in enumerate(model.alphabet)} for position in (0, 1)},
        pairwise={(0, 1): {(left, right): -model.couplings[i][j] for i, left in enumerate(model.alphabet) for j, right in enumerate(model.alphabet)}},
        metadata={"source": "synthetic Week 1 fixture", "adapter": "E=-S; canonical i<j counted once"},
    )
    encoding = compile_domain_wall(categorical, penalty=2.0)
    valid_rows = []
    for symbols in itertools.product(model.alphabet, repeat=2):
        sequence = "".join(symbols)
        sample = dict(enumerate(symbols))
        bits = encoding.encode(sample)
        values = (model.energy(sequence), categorical.energy(sample), encoding.qubo_energy(bits), encoding.ising_model.energy(bits, vartype="BINARY"))
        if any(value != REFERENCE_ENERGIES[sequence] for value in values):
            raise AssertionError(f"Original/categorical/QUBO/Ising energy mismatch for {sequence}")
        if encoding.decode(bits) != sample:
            raise AssertionError("Gibbsiq codec failed round trip")
        valid_rows.append({"sequence": sequence, "energy": values[0], "wall_bits": [bits[label] for label in encoding.variables]})
    valid_words = invalid_words = 0
    for word in itertools.product((0, 1), repeat=len(encoding.variables)):
        bits = dict(zip(encoding.variables, word))
        if encoding.is_valid(bits):
            valid_words += 1
        else:
            invalid_words += 1
            try:
                encoding.decode(bits)
            except ValueError:
                pass
            else:
                raise AssertionError("Gibbsiq decoded an invalid endpoint")
    if (valid_words, invalid_words) != (9, 7):
        raise AssertionError("Unexpected number of tiny encoded states")
    if git_read(source, "status", "--porcelain"):
        raise AssertionError("Upstream checkout changed during read-only exercise")
    audit_files = ("LICENSE", "pyproject.toml", "src/gibbsiq/categorical.py", "src/gibbsiq/domain_wall.py", "src/gibbsiq/reference_sampler.py", "src/gibbsiq/verification.py")
    return {
        "exercise": "week03_pinned_gibbsiq_synthetic_adapter", "passed": True,
        "python_version": platform.python_version(), "upstream_version": gibbsiq.__version__,
        "source_checkout": str(source), "repository": PINNED_REPOSITORY, "commit": PINNED_COMMIT,
        "upstream_tree_unchanged": True, "vendored_code": False,
        "license_scope": "MIT core; aggregate optional package MIT AND Apache-2.0; no Z1T adapter used",
        "upstream_file_sha256": {name: hashlib.sha256((source / name).read_bytes()).hexdigest() for name in audit_files},
        "fixture_sha256": hashlib.sha256(fixture_bytes).hexdigest(),
        "states_checked": valid_rows, "valid_word_count": valid_words, "invalid_word_count": invalid_words,
        "scope": "energy adapter and codec only; no encoded sampling or biology validation",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="trusted clean checkout of the pinned Gibbsiq repository")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "output/upstream_reuse.json")
    args = parser.parse_args()
    result = run_check(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Gibbsiq {result['upstream_version']} at {PINNED_COMMIT}: 9 energy/codec checks passed")
    print("Enumerated 16 words: 9 valid, 7 rejected invalid endpoints; upstream tree unchanged")
    print(f"Saved {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
