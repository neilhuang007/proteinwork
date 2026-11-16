"""Prepare traceable candidate annotations and a sampling adequacy review."""

from __future__ import annotations

import argparse
import csv
import hashlib
from itertools import combinations
from pathlib import Path
import sys
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from examples.week05.shared_analysis import mapped_changes, residue_profile
from proteinwork.model import _parseFasta, load_benchmark
from proteinwork.validation import digest, read_json, save_json
from proteinwork.workflow import verify_run


def write_csv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def diagnostic_audit(diagnostics, evaluation):
    threshold = evaluation["broad_case_reference"]["frequency_mcse_max"]
    blockers = []
    for row in diagnostics["frequencies"]:
        if row["mcse"] is None or row["mcse"] > threshold:
            blockers.append({"key": row["key"], "kind": row["kind"],
                             "model_indices": ";".join(map(str, row["model_indices"])),
                             "symbols": row["symbols"], "probability": row["probability"],
                             "mcse": row["mcse"], "reference_mcse_limit": threshold,
                             "reason": "unknown_mcse" if row["mcse"] is None else "above_limit"})
    return {"sampling_status": diagnostics["status"],
            "failed_sampling_checks": [name for name, passed in diagnostics["checks"].items() if not passed],
            "reference_qualified": diagnostics["reference_qualification"]["qualified"],
            "failed_reference_checks": [name for name, passed in diagnostics["reference_qualification"]["checks"].items() if not passed],
            "frequency_observables": len(diagnostics["frequencies"]),
            "unknown_frequency_mcse": sum(row["reason"] == "unknown_mcse" for row in blockers),
            "frequency_mcse_above_limit": sum(row["reason"] == "above_limit" for row in blockers),
            "observables": diagnostics["observables"]}, blockers


def validate_candidate(row, inputs, full_states, energies, allowed_symbols):
    sequence = row["aligned_sequence"]
    sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
    if row["sequence_sha256"] != sequence_hash or row["candidate_id"] != "cm_" + sequence_hash[:16]:
        raise ValueError("Candidate identity does not match its sequence")
    tokens = inputs.model.encode(sequence)
    changes = mapped_changes(inputs.anchor, sequence, inputs.mapping)
    if any(symbol not in allowed for symbol, allowed in zip(sequence, allowed_symbols, strict=True)):
        raise ValueError("Candidate violates the saved residue rules")
    chain, draw = int(row["source_chain"]), int(row["source_draw"])
    if not (0 <= chain < full_states.shape[0] and 0 <= draw < full_states.shape[1]):
        raise ValueError("Candidate source observation is outside the trace")
    if not np.array_equal(full_states[chain, draw], tokens):
        raise ValueError("Candidate does not match its source observation")
    occurrences = int(np.all(full_states == tokens, axis=-1).sum())
    if occurrences != int(row["occurrences"]):
        raise ValueError("Candidate occurrence count disagrees with the complete trace")
    energy = float(inputs.model.energy(tokens))
    if not np.isclose(energy, float(row["energy"]), atol=1e-8, rtol=1e-11) or not np.isclose(
            energy, energies[chain, draw], atol=1e-8, rtol=1e-11):
        raise ValueError("Candidate energy disagrees with its original model or trace")
    if row["sequence"] != sequence.replace("-", "") or int(row["hamming_to_anchor"]) != len(changes):
        raise ValueError("Candidate sequence export or anchor distance is inconsistent")
    return changes


def validate_reference(source, reference):
    for field in ("target_sha256", "evaluation_sha256"):
        if source[field] != reference[field]:
            raise ValueError("Reference target or evaluation rules differ from the candidate run")
    if reference["purpose"] != "reference" or reference["method"] != "direct":
        raise ValueError("Supply a separate direct reference run")
    if set(source["seeds"]) & set(reference["seeds"]):
        raise ValueError("Reference and candidate run seeds must be disjoint")


def build_packet(run_directory, benchmark_path, output, reference_directory=None):
    run_directory, benchmark_path, output = map(lambda path: Path(path).resolve(),
                                               (run_directory, benchmark_path, output))
    if output.exists():
        raise ValueError("Use a new review directory to preserve previous results")
    verify_run(run_directory)
    manifest = read_json(run_directory / "manifest.json")
    evaluation = read_json(run_directory / "evaluation.json")
    inputs = load_benchmark(benchmark_path)
    if inputs.config != read_json(run_directory / "benchmark.json"):
        raise ValueError("Review benchmark differs from the source run")
    rules = read_json(run_directory / "rules.json")
    source_audit, source_blockers = diagnostic_audit(read_json(run_directory / "diagnostics.json"), evaluation)
    reference_audit = reference_manifest = None
    reference_blockers = []
    if reference_directory is not None:
        reference_directory = Path(reference_directory).resolve()
        verify_run(reference_directory)
        reference_manifest = read_json(reference_directory / "manifest.json")
        validate_reference(manifest, reference_manifest)
        reference_audit, reference_blockers = diagnostic_audit(
            read_json(reference_directory / "diagnostics.json"), evaluation)
        reference_audit.update(run_status=reference_manifest["status"],
                               raw_observations=reference_manifest["raw_observations"],
                               elapsed_seconds=reference_manifest["elapsed_seconds"],
                               budget_exceeded=reference_manifest["budget_exceeded"])

    archive_path = benchmark_path.parent.parent / inputs.config["archive"]["local_path"]
    if digest(archive_path) != inputs.config["archive"]["sha256"]:
        raise ValueError("Alignment archive changed after validation")
    with zipfile.ZipFile(archive_path) as archive:
        alignment_bytes = archive.read(inputs.config["coordinate_reference"]["archive_member"])
    alignment = list(_parseFasta(alignment_bytes).values())
    if any(len(sequence) != inputs.model.n_sites or set(sequence) - set(inputs.model.alphabet) for sequence in alignment):
        raise ValueError("Alignment differs from the modeled sequence alphabet or length")
    profiles = [residue_profile(alignment, site) for site in range(inputs.model.n_sites)]
    with (run_directory / "candidates.csv").open(encoding="utf-8", newline="") as stream:
        candidates = list(csv.DictReader(stream))
    if len(candidates) != manifest["exported_candidates"] or len({row["candidate_id"] for row in candidates}) != len(candidates):
        raise ValueError("Candidate count or identities disagree with the manifest")
    annotations, substitutions = [], []
    with np.load(run_directory / "trace.npz", allow_pickle=False) as trace:
        full_states, energies = trace["full_states"], trace["energies"]
        for candidate in candidates:
            changes = validate_candidate(candidate, inputs, full_states, energies,
                                         rules["allowed_symbols_by_model_index"])
            annotations.append({"candidate_id": candidate["candidate_id"],
                                "selection_rank": int(candidate["selection_rank"]),
                                "substitution_count": len(changes), "substitutions": "; ".join(changes),
                                "model_energy": float(candidate["energy"]),
                                "occurrences": int(candidate["occurrences"]),
                                "source_chain": int(candidate["source_chain"]),
                                "source_draw": int(candidate["source_draw"]),
                                "ungapped_length": len(candidate["sequence"]),
                                "aligned_sequence": candidate["aligned_sequence"]})
            for site, (before, after) in enumerate(zip(inputs.anchor, candidate["aligned_sequence"], strict=True)):
                if before == after:
                    continue
                profile = profiles[site]
                substitutions.append({"candidate_id": candidate["candidate_id"], "model_index_0based": site,
                                      "accession_position": inputs.mapping[site]["uniprot_position_1based"],
                                      "anchor_residue": before, "candidate_residue": after,
                                      "alignment_count": profile["counts"].get(after, 0),
                                      "alignment_nongap_count": profile["nongap_records"],
                                      "alignment_frequency": profile["frequencies_among_nongap"].get(after, 0),
                                      "allowed_symbols": rules["allowed_symbols_by_model_index"][site]})
    distances = [{"candidate_a": left["candidate_id"], "candidate_b": right["candidate_id"],
                  "aligned_hamming_distance": sum(a != b for a, b in zip(left["aligned_sequence"], right["aligned_sequence"], strict=True))}
                 for left, right in combinations(candidates, 2)]
    minimum_distance = min((row["aligned_hamming_distance"] for row in distances), default=None)
    requested_distance = rules["selection_rules"].get("min_pairwise_hamming", 0)
    if minimum_distance is not None and minimum_distance < requested_distance:
        raise ValueError("Selected candidates violate the saved pairwise distance rule")
    summary = {"schema_version": 1, "source_status": manifest["status"], "source_audit": source_audit,
               "reference_audit": reference_audit, "selected_candidates": len(candidates),
               "substitution_rows": len(substitutions), "alignment_records": len(alignment),
               "minimum_pairwise_hamming": minimum_distance,
               "requested_minimum_pairwise_hamming": requested_distance,
               "human_review_status": "pending", "new_biological_measurements": 0,
               "limits": ["Alignment frequencies are unweighted descriptions, not functional predictions.",
                          "Selection changes the distribution represented by the candidate batch.",
                          "Reviewing a new reference does not retroactively qualify the source run.",
                          "Blank review fields require a human decision before new biological restrictions are used."]}
    output.mkdir(parents=True)
    write_csv(output / "candidate_annotations.csv", annotations, ["candidate_id", "selection_rank", "substitution_count", "substitutions", "model_energy", "occurrences", "source_chain", "source_draw", "ungapped_length", "aligned_sequence"])
    write_csv(output / "substitutions.csv", substitutions, ["candidate_id", "model_index_0based", "accession_position", "anchor_residue", "candidate_residue", "alignment_count", "alignment_nongap_count", "alignment_frequency", "allowed_symbols"])
    write_csv(output / "pairwise_distances.csv", distances, ["candidate_a", "candidate_b", "aligned_hamming_distance"])
    for name, blockers in (("source", source_blockers), ("reference", reference_blockers)):
        write_csv(output / f"{name}_frequency_blockers.csv", blockers, ["key", "kind", "model_indices", "symbols", "probability", "mcse", "reference_mcse_limit", "reason"])
    review_fields = ["candidate_id", "aligned_sequence", "reviewer", "decision", "reason", "proposed_position", "proposed_allowed_symbols", "evidence_source"]
    write_csv(output / "review_worksheet.csv", [{field: row[field] if field in ("candidate_id", "aligned_sequence") else "" for field in review_fields} for row in candidates], review_fields)
    save_json(output / "summary.json", summary)
    save_json(output / "provenance.json", {
        "source_manifest_sha256": digest(run_directory / "manifest.json"),
        "reference_manifest_sha256": digest(reference_directory / "manifest.json") if reference_manifest else None,
        "script_sha256": digest(Path(__file__)),
        "analysis_helpers_sha256": digest(ROOT / "examples/week05/shared_analysis.py"),
        "alignment_sha256": hashlib.sha256(alignment_bytes).hexdigest(),
        "benchmark_sha256": digest(benchmark_path),
        "artifacts": {path.name: digest(path) for path in sorted(output.iterdir()) if path.is_file()}})
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, default=ROOT / "configs/benchmark.json")
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = build_packet(args.run, args.benchmark, args.output, args.reference)
    print(f"Prepared {summary['selected_candidates']} candidates and {summary['substitution_rows']} substitution records")


if __name__ == "__main__":
    main()
