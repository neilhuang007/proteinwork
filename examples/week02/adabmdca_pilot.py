"""Run the official documented adabmDCA model.sample API on a synthetic model."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from examples.week01.model_reading import run_exercise
from examples.week02.pilot_support import frequencies, runtime_record, save_json, sha256


def run_pilot() -> dict:
    configuration = {
        "alphabet": "-AC",
        "device": "cpu",
        "dtype": "float32",
        "n_sequences": 2048,
        "n_sweeps": 128,
        "sampler": "metropolis",
        "beta": 1.0,
        "seed": 42,
        "torch_num_threads": 1,
    }
    fixture = ROOT / "examples/week01/fixtures/two_site_params.dat"
    reference = run_exercise()
    reference_sequences = [row["sequence"] for row in reference["states"]]
    with warnings.catch_warnings(record=True) as observed_warnings:
        warnings.simplefilter("always")
        import torch
        import adabmDCA
        import adabmDCA.api.model
        import adabmDCA.api.sampling
        from adabmDCA import load_model

        torch.set_num_threads(configuration["torch_num_threads"])
        model = load_model(
            fixture,
            alphabet=configuration["alphabet"],
            device=configuration["device"],
            dtype=configuration["dtype"],
        )
        official_energies = model.compute_energies(reference_sequences).tolist()
        energy_checks = []
        for row, energy in zip(reference["states"], official_energies, strict=True):
            if not math.isclose(energy, row["reference_energy"], rel_tol=0, abs_tol=1e-6):
                raise ValueError(f"Official energy disagrees with reference for {row['sequence']}.")
            energy_checks.append({
                "sequence": row["sequence"],
                "official_energy": energy,
                "hand_reference_energy": row["reference_energy"],
                "absolute_error": abs(energy - row["reference_energy"]),
            })
        started = time.perf_counter()
        sequences = list(model.sample(
            configuration["n_sequences"],
            n_sweeps=configuration["n_sweeps"],
            sampler=configuration["sampler"],
            beta=configuration["beta"],
            seed=configuration["seed"],
        ))
        elapsed = time.perf_counter() - started
        if len(sequences) != configuration["n_sequences"]:
            raise ValueError("The official call returned an unexpected number of sequences.")
        source_modules = [adabmDCA, adabmDCA.api.model, adabmDCA.api.sampling]
        source_hashes = {module.__name__: sha256(Path(module.__file__)) for module in source_modules}
        metadata = model.metadata.to_dict()
    probabilities = {row["sequence"]: row["probability"] for row in reference["states"]}
    return {
        "pilot": "week02_official_adabmDCA_synthetic_sequence_generation",
        "scope": "Synthetic two-position -AC fixture; no training, biological benchmark or fitness claim.",
        "runtime": runtime_record(),
        "configuration": configuration,
        "fixture": {"path": str(fixture.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(fixture)},
        "official_api": "adabmDCA.load_model; DCAModel.compute_energies; DCAModel.sample",
        "official_source_sha256": source_hashes,
        "model_metadata": metadata,
        "energy_checks": energy_checks,
        "comparison": frequencies(sequences, probabilities),
        "sequences": sequences,
        "sample_call_elapsed_seconds": elapsed,
        "timing_scope": "Single CPU call after loading and energy checks; includes first-call compilation if any, not a speed benchmark.",
        "warnings": [{"category": type(item.message).__name__, "message": str(item.message)} for item in observed_warnings],
        "limitations": [
            "Finite Monte Carlo error and repeats are expected; no convergence or effective-sample-size claim.",
            "Exact reference applies only to the nine-state synthetic fixture.",
            "No reference alignment, mixing-time diagnostics, fitted protein model or biological validation.",
        ],
        "documentation": "https://spqb.github.io/adabmDCApy/high_level_api/#generate-sequences",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "examples/week02/output/adabmdca.json")
    args = parser.parse_args()
    result = run_pilot()
    save_json(args.output, result)
    fasta = args.output.with_suffix(".fasta")
    fasta.write_text("".join(f">synthetic_{index}\n{sequence}\n" for index, sequence in enumerate(result["sequences"])), encoding="utf-8")
    comparison = result["comparison"]
    print(f"Official adabmDCA synthetic pilot: {comparison['n_returned']} sequences, {comparison['n_unique']} unique")
    print(f"Energy checks: 9/9; total variation: {comparison['total_variation_distance']:.6f}")
    print(f"Sample call: {result['sample_call_elapsed_seconds']:.3f} seconds; warnings: {len(result['warnings'])}")
    print(f"Saved {args.output} and {fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
