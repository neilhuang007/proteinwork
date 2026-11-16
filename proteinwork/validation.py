"""Strict public configuration validation and provenance utilities."""

from __future__ import annotations

import hashlib
import json
import math
from numbers import Integral
from pathlib import Path


class SourceDocument(dict):
    """Dictionary carrying the fingerprint of the bytes actually parsed."""


def _uniqueObjectPairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _integerControl(value, name, minimum):
    """Validate sampler controls, including NumPy integer scalar inputs."""
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return int(value)


def read_json(path: str | Path) -> dict:
    def constant(value):
        raise ValueError(f"Nonfinite JSON value: {value}")

    data = Path(path).read_bytes()
    value = json.loads(data.decode("utf-8"),
                       object_pairs_hook=_uniqueObjectPairs, parse_constant=constant)
    if not isinstance(value, dict):
        raise ValueError("The JSON document must be an object")
    document = SourceDocument(value)
    document.source_sha256 = hashlib.sha256(data).hexdigest()
    return document


def keys(value: dict, allowed: set[str], required: set[str], name: str) -> None:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be an object")
    unknown, missing = set(value) - allowed, required - set(value)
    if unknown or missing:
        raise ValueError(f"{name}: unknown keys {sorted(unknown)}, missing keys {sorted(missing)}")


def integer(value, name: str, minimum: int = 0, maximum: int | None = None) -> int:
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f"{name} must be an integer in the declared range")
    return value


def positive(value, name: str, allow_zero: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or (result < 0 if allow_zero else result <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")
    return result


def digest(path: str | Path) -> str:
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def object_digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def save_json(path: str | Path, value) -> None:
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_request(path: str | Path) -> tuple[dict, Path]:
    path = Path(path).resolve()
    request = read_json(path)
    keys(request, {"schema_version", "benchmark_config", "case_id", "beta", "sampling",
                   "restrictions", "selection", "evaluation_config", "reference_run",
                   "assays", "planner", "allow_unqualified_export", "execution", "purpose"},
         {"schema_version", "benchmark_config", "case_id", "sampling", "selection",
          "evaluation_config"}, "request")
    if type(request["schema_version"]) is not int or request["schema_version"] != 1:
        raise ValueError("Unsupported request schema_version")
    if request.get("execution", "simulation") != "simulation":
        raise ValueError("Physical execution requires a verified device adapter; only simulation is configured")
    if request.get("purpose", "candidate_library") not in ("candidate_library", "reference"):
        raise ValueError("purpose must be candidate_library or reference")
    for key in ("case_id", "benchmark_config", "evaluation_config"):
        if not isinstance(request[key], str) or not request[key]:
            raise ValueError(f"{key} must be a nonempty string")
    if "beta" in request:
        positive(request["beta"], "beta")
    if type(request.get("allow_unqualified_export", False)) is not bool:
        raise ValueError("allow_unqualified_export must be a boolean")
    sampling = request["sampling"]
    keys(sampling, {"method", "backend", "seeds", "draws", "warmup", "sweeps_per_draw",
                    "penalty", "proposal_beta", "coefficient_scale", "quantization_step",
                    "max_spins", "max_edges", "max_bytes", "checkpoint_policy"},
         {"method", "seeds", "draws", "warmup"}, "sampling")
    if sampling["method"] not in ("direct", "encoded", "planner"):
        raise ValueError("method must be direct, encoded or planner")
    if sampling.get("checkpoint_policy", "final") not in ("final", "first_qualified"):
        raise ValueError("checkpoint_policy must be final or first_qualified")
    if sampling.get("backend", "numpy") not in ("numpy", "thrml", "jax"):
        raise ValueError("Unknown encoded backend")
    seeds = sampling["seeds"]
    if not isinstance(seeds, list) or len(seeds) != 4 or len(set(seeds)) != 4:
        raise ValueError("Exactly four distinct chain seeds are required")
    for seed in seeds:
        integer(seed, "seed", 0, 2**32 - 1)
    integer(sampling["draws"], "draws", 4)
    integer(sampling["warmup"], "warmup", 0)
    integer(sampling.get("sweeps_per_draw", 1), "sweeps_per_draw", 1)
    for key in ("max_spins", "max_edges", "max_bytes"):
        if key in sampling:
            integer(sampling[key], key, 1)
    for key in ("penalty", "proposal_beta", "coefficient_scale", "quantization_step"):
        if key in sampling:
            positive(sampling[key], key)
    if sampling["method"] == "direct" and set(sampling) & {
            "backend", "penalty", "proposal_beta", "coefficient_scale", "quantization_step"}:
        raise ValueError("Encoded proposal settings cannot be attached to the direct sampler")
    if sampling["method"] == "planner" and not request.get("planner"):
        raise ValueError("Planner selection requires a frozen planner file")
    if request.get("purpose") == "reference" and sampling["method"] != "direct":
        raise ValueError("Independent reference runs must use direct sampling")
    if request.get("purpose") == "reference":
        selection = request["selection"]
        if (not isinstance(selection, dict) or type(selection.get("batch_size")) is not int
                or selection["batch_size"] != 0
                or selection.get("objective", "model_score") != "model_score"
                or request.get("allow_unqualified_export", False)
                or "assays" in request or "reference_run" in request):
            raise ValueError("Reference runs require batch_size zero and model_score, without assays, another reference or candidate export")
    for name in ("reference_run", "assays", "planner"):
        if name in request and (not isinstance(request[name], str) or not request[name]):
            raise ValueError(f"{name} must be a path string")
    if not isinstance(request.get("restrictions", []), list):
        raise ValueError("restrictions must be an array")
    return request, path.parent


def load_evaluation(path):
    value = read_json(path)
    required = {"schema_version", "status", "development_case_ids", "evaluation_case_ids", "methods",
                "temperature", "precision", "chains_per_condition", "development_seeds", "evaluation_seeds",
                "reference_seeds", "initialization", "selected_pair_observables", "exact_kernel_checks",
                "sampling_accuracy_targets", "broad_case_reference", "budget", "timing", "review_note"}
    keys(value, required, required, "evaluation")
    if value["schema_version"] != 1 or value["precision"] != "float64" or value["chains_per_condition"] != 4:
        raise ValueError("Only the version 1 float64 four chain protocol is supported")
    positive(value["temperature"], "evaluation temperature")
    for key in ("development_case_ids", "evaluation_case_ids", "methods"):
        values = value[key]
        if not isinstance(values, list) or not values or any(not isinstance(v, str) or not v for v in values) or len(set(values)) != len(values):
            raise ValueError(f"Invalid {key}")
    if set(value["development_case_ids"]) & set(value["evaluation_case_ids"]):
        raise ValueError("Development and evaluation case sets must be disjoint")
    seed_sets = []
    for key in ("development_seeds", "evaluation_seeds", "reference_seeds"):
        values = value[key]
        if not isinstance(values, list) or len(values) != 4:
            raise ValueError("Evaluation seed groups require four seeds")
        for v in values:
            integer(v, key, 0, 2**32 - 1)
        if len(set(values)) != 4 or any(set(values) & prior for prior in seed_sets):
            raise ValueError("Evaluation seed groups must be internally unique and disjoint")
        seed_sets.append(set(values))
    target_keys = {"exact_case_total_variation_max", "single_site_probability_absolute_error_max",
                   "selected_pair_probability_absolute_error_max", "energy_mean_error_in_reference_standard_deviations_max",
                   "hamming_mean_error_in_reference_standard_deviations_max", "rank_normalized_split_rhat_max_exclusive",
                   "bulk_ess_min_inclusive", "tail_ess_min_inclusive", "diagnostic_scope", "constant_observable_rule"}
    keys(value["sampling_accuracy_targets"], target_keys, target_keys, "accuracy targets")
    for key, val in value["sampling_accuracy_targets"].items():
        if key not in ("diagnostic_scope", "constant_observable_rule"):
            positive(val, key)
    reference_keys = {"method", "chains", "rank_normalized_split_rhat_max_exclusive", "bulk_and_tail_ess_min",
                      "frequency_mcse_max", "mean_mcse_in_reference_standard_deviations_max", "status", "uncertainty_rule"}
    keys(value["broad_case_reference"], reference_keys, reference_keys, "reference rules")
    if value["broad_case_reference"]["chains"] != 4:
        raise ValueError("Reference qualification requires four chains")
    for key in ("rank_normalized_split_rhat_max_exclusive", "bulk_and_tail_ess_min", "frequency_mcse_max",
                "mean_mcse_in_reference_standard_deviations_max"):
        positive(value["broad_case_reference"][key], key)
    budget_keys = {"tuning_seconds_per_method_per_case", "measurement_seconds_per_condition",
                   "reference_seconds_per_case_max", "failure_rule"}
    keys(value["budget"], budget_keys, budget_keys, "budget")
    for key, val in value["budget"].items():
        if key != "failure_rule":
            positive(val, key)
    checkpoints = value["timing"]["checkpoints_post_warmup_observations_per_chain"]
    if not isinstance(checkpoints, list) or not checkpoints or checkpoints != sorted(set(checkpoints)):
        raise ValueError("Accuracy checkpoints must be distinct increasing integers")
    for count in checkpoints:
        integer(count, "accuracy checkpoint", 4)
    return value


def apply_residue_requests(inputs, domains, requests: list[dict]):
    domains = list(domains)
    mapping = {int(row["uniprot_position_1based"]): int(row["model_index_0based"])
               for row in inputs.mapping if row.get("uniprot_position_1based") not in (None, "")}
    seen, evidence = set(), []
    accession = inputs.config["coordinate_reference"]["uniprot_accession"]
    for row in requests:
        keys(row, {"model_index_0based", "uniprot_position_1based", "accession", "allowed_symbols",
                   "expected_anchor_symbol", "reason"}, {"allowed_symbols", "reason"}, "restriction")
        if ("model_index_0based" in row) == ("uniprot_position_1based" in row):
            raise ValueError("Specify exactly one residue coordinate system")
        if not isinstance(row["reason"], str) or not row["reason"].strip():
            raise ValueError("Each researcher restriction requires its reason")
        if "uniprot_position_1based" in row:
            position = integer(row["uniprot_position_1based"], "UniProt position", 1)
            if row.get("accession") != accession or position not in mapping:
                raise ValueError("Unmapped biological coordinate or accession mismatch")
            site = mapping[position]
        else:
            site = integer(row["model_index_0based"], "model index", 0, len(domains) - 1)
            if "accession" in row and row["accession"] != accession:
                raise ValueError("Restriction accession mismatch")
        if site in seen:
            raise ValueError("Duplicate restriction at the same modeled residue")
        seen.add(site)
        if "expected_anchor_symbol" in row and row["expected_anchor_symbol"] != inputs.anchor[site]:
            raise ValueError("Expected reference residue does not match coordinate anchor")
        symbols = row["allowed_symbols"]
        if not isinstance(symbols, str) or not symbols or len(set(symbols)) != len(symbols):
            raise ValueError("Allowed residues must be distinct explicit symbols")
        if set(symbols) - set(inputs.model.alphabet):
            raise ValueError("Unknown residue in restriction")
        tokens = tuple(inputs.model.alphabet.index(s) for s in symbols)
        if set(tokens) - set(domains[site]):
            raise ValueError("Requested residue conflicts with original support or base case restrictions")
        domains[site] = tuple(sorted(tokens))
        evidence.append({**row, "resolved_model_index_0based": site, "status": "enforced"})
    return tuple(domains), evidence
