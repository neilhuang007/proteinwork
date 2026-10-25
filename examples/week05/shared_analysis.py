"""Connect a residue counting lesson to biological interpretation of a saved run."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from proteinwork.direct import enumerate_target
from proteinwork.model import _parseFasta, benchmark_domains, load_benchmark, restrict_model
from proteinwork.validation import digest, save_json
from proteinwork.workflow import verify_run


def residue_profile(sequences, model_index):
    """Count every alignment record; residue proportions exclude gaps."""
    if not sequences or len({len(sequence) for sequence in sequences}) != 1:
        raise ValueError("Supply a nonempty alignment with equal sequence lengths")
    if not 0 <= model_index < len(sequences[0]):
        raise ValueError("Position is outside the alignment")
    counts = Counter(sequence[model_index] for sequence in sequences)
    nongap = len(sequences) - counts["-"]
    frequencies = {
        residue: count / nongap
        for residue, count in sorted(counts.items()) if residue != "-"
    }
    return {"records": len(sequences), "nongap_records": nongap,
            "gap_records": counts["-"], "counts": dict(sorted(counts.items())),
            "frequencies_among_nongap": frequencies}


def mapped_changes(anchor, sequence, mapping):
    """Name substitutions against the unchanged accession coordinate anchor."""
    if len(sequence) != len(anchor) or len(mapping) != len(anchor):
        raise ValueError("Sequence, anchor and position map must have equal lengths")
    changes = []
    for index, (before, after, row) in enumerate(zip(anchor, sequence, mapping, strict=True)):
        if row["model_index_0based"] != index or row["anchor_symbol"] != before:
            raise ValueError("Position map disagrees with the anchor")
        if before == after:
            continue
        position = row["uniprot_position_1based"]
        if position is None or before == "-" or after == "-":
            raise ValueError("This substitution exercise requires the anchor gap pattern")
        changes.append(f"{before}{position}{after}")
    return changes


def pair_comparison(sequences, positions, exact_sequences, probabilities):
    """Preserve repeat counts when comparing a two site trace with its exact target."""
    if not sequences or len(positions) != 2 or len(set(positions)) != 2:
        raise ValueError("Supply observations and two distinct positions")
    counts = Counter("".join(sequence[index] for index in positions) for sequence in sequences)
    expected = {}
    for sequence, probability in zip(exact_sequences, probabilities, strict=True):
        key = "".join(sequence[index] for index in positions)
        expected[key] = expected.get(key, 0.0) + float(probability)
    if set(counts) - set(expected) or not np.isclose(sum(expected.values()), 1):
        raise ValueError("Observed states or probabilities disagree with the exact support")
    rows = [{"residue_pair": pair, "count": counts[pair],
             "observed_probability": counts[pair] / len(sequences),
             "exact_probability": probability}
            for pair, probability in sorted(expected.items())]
    return rows, 0.5 * sum(abs(row["observed_probability"] - row["exact_probability"]) for row in rows)


def write_table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(run_directory, output):
    run_directory, output = Path(run_directory).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError("Use a new analysis directory to preserve previous results")
    receipt = verify_run(run_directory)
    request = json.loads((run_directory / "request.json").read_text(encoding="utf-8"))
    if request["case_id"] != "cm_dev_exact_8_12" or request.get("restrictions"):
        raise ValueError("This lesson uses the unchanged two position development case")
    inputs = load_benchmark(ROOT / "configs/benchmark.json")
    saved_config = json.loads((run_directory / "benchmark.json").read_text(encoding="utf-8"))
    if saved_config != inputs.config:
        raise ValueError("Saved run and analysis must use the same benchmark")
    model = inputs.model
    domains = benchmark_domains(inputs, request["case_id"])
    reduced = restrict_model(model, domains)
    positions = reduced.original_indices
    exact_states, probabilities = enumerate_target(reduced, request.get("beta", inputs.config["model"]["beta"]))
    exact_sequences = [model.decode(state) for state in reduced.expand(exact_states)]
    with np.load(run_directory / "trace.npz", allow_pickle=False) as trace:
        full_states = trace["full_states"]
        shape = full_states.shape
        full_flat = full_states.reshape(-1, model.n_sites)
        unique_states, inverse = np.unique(full_flat, axis=0, return_inverse=True)
        unique_sequences = [model.decode(state) for state in unique_states]
        if set(unique_sequences) - set(exact_sequences):
            raise ValueError("Trace contains states outside the lesson's exact support")
        energies = model.energy(unique_states)[inverse].reshape(shape[:2])
        np.testing.assert_allclose(energies, trace["energies"], atol=1e-8, rtol=1e-11)
        np.testing.assert_array_equal(reduced.expand(trace["states"]), full_states)
        sequences = [unique_sequences[index] for index in inverse]
        repeated_endpoints = int(np.all(full_states[:, 1:] == full_states[:, :-1], axis=-1).sum())

    archive_path = ROOT / inputs.config["archive"]["local_path"]
    member = inputs.config["coordinate_reference"]["archive_member"]
    if digest(archive_path) != inputs.config["archive"]["sha256"]:
        raise ValueError("Alignment archive changed after model loading")
    with zipfile.ZipFile(archive_path) as archive:
        alignment_bytes = archive.read(member)
    natural = list(_parseFasta(alignment_bytes).values())
    if any(len(sequence) != model.n_sites or set(sequence) - set(model.alphabet) for sequence in natural):
        raise ValueError("Natural sequences must use the model's aligned alphabet and length")

    profiles = []
    for row in inputs.mapping:
        index = row["model_index_0based"]
        profile = residue_profile(natural, index)
        profiles.append({**row, **profile})
    frequency_rows = []
    for index in positions:
        profile = profiles[index]
        sampled = residue_profile(sequences, index)
        for token in domains[index]:
            residue = model.alphabet[token]
            frequency_rows.append({"accession_position": profile["uniprot_position_1based"],
                "model_index": index, "residue": residue,
                "natural_count": profile["counts"].get(residue, 0),
                "natural_nongap_count": profile["nongap_records"],
                "natural_frequency": profile["frequencies_among_nongap"].get(residue, 0.0),
                "sampled_count": sampled["counts"].get(residue, 0),
                "sampled_observations": len(sequences),
                "sampled_frequency": sampled["frequencies_among_nongap"].get(residue, 0.0)})

    pair_rows, total_variation = pair_comparison(sequences, positions, exact_sequences, probabilities)
    witness_changes = set(mapped_changes(inputs.anchor, inputs.witness, inputs.mapping))
    annotation_rows = []
    counts = Counter(sequences)
    for sequence in sorted(counts):
        changes = mapped_changes(inputs.anchor, sequence, inputs.mapping)
        sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
        annotation_rows.append({"candidate_id": "cm_" + sequence_hash[:16],
            "residue_pair": "".join(sequence[index] for index in positions),
            "occurrences": counts[sequence], "accession_changes": "; ".join(changes),
            "background_changes": "; ".join(change for change in changes if change in witness_changes),
            "variable_changes": "; ".join(change for change in changes if change not in witness_changes),
            "ungapped_length": len(sequence.replace("-", "")), "aligned_sequence": sequence})

    excluded_anchor_choices = []
    for index, token in enumerate(model.encode(inputs.anchor)):
        if not np.isfinite(model.fields[index, token]):
            excluded_anchor_choices.append({"accession_position": inputs.mapping[index]["uniprot_position_1based"],
                "model_index": index, "residue": inputs.anchor[index],
                "natural_count": profiles[index]["counts"].get(inputs.anchor[index], 0),
                "witness_residue": inputs.witness[index]})

    diagnostics = json.loads((run_directory / "diagnostics.json").read_text(encoding="utf-8"))
    checks = {"source_run_verified": receipt["status"] == "verified",
              "full_model_scores_match_trace": True,
              "all_states_in_exact_support": True,
              "all_candidates_have_89_residues": all(row["ungapped_length"] == 89 for row in annotation_rows),
              "all_sampling_accuracy_checks_passed": all(diagnostics["checks"].values()),
              "pair_counts_cover_all_observations": sum(row["count"] for row in pair_rows) == len(sequences)}
    summary = {"case_id": request["case_id"], "natural_alignment_records": len(natural),
               "alignment_columns": model.n_sites, "mapped_reference_residues": sum(row["uniprot_position_1based"] is not None for row in inputs.mapping),
               "accession": inputs.config["coordinate_reference"]["uniprot_accession"],
               "focus_positions": [inputs.mapping[index] for index in positions],
               "natural_frequency_definition": "Unweighted count divided by nongap records at this aligned column",
               "sample_frequency_definition": "Every post warmup observation of the restricted model, including repeats",
               "chains": shape[0], "draws_per_chain": shape[1], "observations": len(sequences),
               "unique_sampled_sequences": len(counts), "repeated_adjacent_endpoints": repeated_endpoints,
               "sampling_status": diagnostics["status"], "total_variation": total_variation,
               "residue_frequencies": frequency_rows, "pair_frequencies": pair_rows,
               "anchor_choices_excluded_by_model": excluded_anchor_choices, "checks": checks}
    if not all(checks.values()):
        raise ValueError(f"Analysis checks failed: {checks}")
    output.mkdir(parents=True)
    save_json(output / "summary.json", summary)
    save_json(output / "alignment_profiles.json", profiles)
    write_table(output / "residue_frequencies.csv", frequency_rows)
    write_table(output / "pair_frequencies.csv", pair_rows)
    write_table(output / "candidate_annotations.csv", annotation_rows)
    save_json(output / "provenance.json", {
        "source_run": str(run_directory.relative_to(ROOT)) if run_directory.is_relative_to(ROOT) else str(run_directory),
        "run_manifest_sha256": digest(run_directory / "manifest.json"),
        "script_sha256": digest(Path(__file__)), "archive_sha256": digest(archive_path),
        "alignment_member": member, "alignment_member_sha256": hashlib.sha256(alignment_bytes).hexdigest(),
        "model_sha256": inputs.provenance["parameter_sha256"],
        "position_map_sha256": digest(ROOT / inputs.config["coordinate_reference"]["map"]),
        "artifacts": {path.name: digest(path) for path in sorted(output.iterdir()) if path.is_file()}})
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    result = analyze(options.run, options.output)
    print(json.dumps(result, indent=2, allow_nan=False))
