"""Research input, validated computation, evidence and library handoff."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import math
import os
from pathlib import Path
import platform
import sys
import time
import uuid

import numpy as np

from .diagnostics import diagnose
from .direct import SamplingResult, sample_direct
from .model import benchmark_domains, load_benchmark, restrict_model
from .planner import choose_method
from .validation import (apply_residue_requests, digest, load_evaluation, load_request, object_digest,
                         read_json, save_json)


def implementation_fingerprint():
    root = Path(__file__).parent
    return object_digest({p.name: digest(p) for p in sorted(root.glob("*.py"))})


def _csv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _runtime():
    packages = {}
    for name in ("numpy", "scipy", "jax", "jaxlib", "thrml", "adabmDCA", "arviz-stats", "scikit-learn"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    processor = os.environ.get("PROCESSOR_IDENTIFIER") or platform.processor() or platform.machine()
    return {"python": sys.version, "platform": platform.platform(), "processor": processor,
            "packages": packages, "thread_settings": {k: os.environ.get(k) for k in (
                "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "XLA_FLAGS")}}


def _target_fingerprint(inputs, domains, beta):
    return object_digest({"model": inputs.config["archive"]["member_sha256"],
                          "domains": domains, "beta": beta, "precision": "float64",
                          "anchor": inputs.anchor, "coordinate_map": inputs.mapping})


def _comparison_context(request, base):
    return object_digest({"runtime": _runtime(), "selection": request["selection"],
                          "assays_sha256": digest(base / request["assays"]) if "assays" in request else None,
                          "allow_unqualified_export": request.get("allow_unqualified_export", False)})


def _reference(path, *, target_sha256, current_seeds, evaluation_sha256):
    directory = Path(path)
    manifest = read_json(directory / "manifest.json")
    if manifest.get("status") != "reference_qualified" or manifest.get("budget_exceeded", True):
        raise ValueError("Reference did not finish as a qualified reference within its runtime budget")
    if manifest["target_sha256"] != target_sha256 or manifest["evaluation_sha256"] != evaluation_sha256:
        raise ValueError("Reference target or evaluation rules differ from the requested experiment")
    if manifest["method"] != "direct":
        raise ValueError("Broad references must use the verified direct sampler")
    if set(manifest["seeds"]) & set(current_seeds):
        raise ValueError("Reference and candidate chains require disjoint seeds")
    if digest(directory / "diagnostics.json") != manifest["artifacts"]["diagnostics.json"]:
        raise ValueError("Reference diagnostic artifact hash mismatch")
    if digest(directory / "trace.npz") != manifest["artifacts"]["trace.npz"]:
        raise ValueError("Reference trace artifact hash mismatch")
    return read_json(directory / "diagnostics.json"), {
        "path": str(directory.resolve()), "manifest_sha256": digest(directory / "manifest.json")}


def _encoding_estimate(reduced):
    sizes = [len(d) for d in reduced.domains]
    bits = sum(sizes)
    within = sum(n * (n - 1) // 2 for n in sizes)
    cross = 0
    for i, left in enumerate(reduced.domains):
        for j in range(i + 1, reduced.n_sites):
            block = reduced.couplings[i, j][np.ix_(left, reduced.domains[j])]
            cross += int(np.count_nonzero(block))
    return {"encoding": "one_hot", "mutable_sites": reduced.n_sites, "logical_spins": bits,
            "nonzero_cross_site_edges": cross, "within_site_penalty_edges": within,
            "logical_edges": cross + within, "dense_matrix_bytes_float64": bits * bits * 8,
            "physical_device_available": False, "physical_execution": False,
            "embedding_status": "Not established; a dense logical graph is not a hardware embedding"}


def validate_request(request_path):
    request, base = load_request(request_path)
    evaluation = load_evaluation(base / request["evaluation_config"])
    inputs = load_benchmark(base / request["benchmark_config"])
    domains, applied = apply_residue_requests(inputs, benchmark_domains(inputs, request["case_id"]),
                                             request.get("restrictions", []))
    beta = request.get("beta", inputs.config["model"]["beta"])
    if (evaluation.get("precision") != "float64" or evaluation.get("chains_per_condition") != 4
            or not math.isclose(beta * evaluation["temperature"], 1.0, rel_tol=1e-12)):
        raise ValueError("Request and evaluation contract differ in precision, chains or temperature")
    reduced = restrict_model(inputs.model, domains)
    settings = dict(request["selection"])
    if "allowed_residues" in settings:
        raise ValueError("Selection cannot override original model support")
    settings["allowed_residues"] = ["".join(inputs.model.alphabet[t] for t in np.flatnonzero(np.isfinite(row)))
                                   for row in inputs.model.fields]
    from .design import validate_design_inputs
    preflight = validate_design_inputs(anchor=inputs.anchor, alphabet=inputs.model.alphabet, settings=settings,
                                       assays_path=(base / request["assays"]) if "assays" in request else None)
    estimate = _encoding_estimate(reduced)
    sampling = request["sampling"]
    trace_bytes = 4 * sampling["draws"] * (inputs.model.n_sites + reduced.n_sites + 3) * 8
    if trace_bytes > sampling.get("max_bytes", 536870912):
        raise ValueError("Requested trace exceeds its memory budget")
    if sampling["method"] == "encoded":
        if estimate["logical_spins"] > sampling.get("max_spins", 2048) or estimate["logical_edges"] > sampling.get("max_edges", 2000000):
            raise ValueError("Encoded problem exceeds its declared graph budget")
        if sampling.get("backend") == "thrml" and estimate["logical_spins"] > 64:
            raise ValueError("THRML sequence adapter supports at most 64 spins; use JAX for broad encoded sampling or THRML batch selection")
    return {"status": "validated", "case_id": request["case_id"], "beta": beta,
            "mutable_sites": reduced.n_sites, "target_state_count": str(math.prod(map(len, domains))),
            "trace_memory_estimate_bytes": trace_bytes, "encoding": estimate,
            "residue_rules": applied, "selection_preflight": preflight,
            "physical_execution": False}


def _sample_with_checkpoints(reduced, *, sampling, beta, anchor, selected_pairs, evaluation, reference, started, budget_seconds):
    checkpoints = [n for n in evaluation["timing"]["checkpoints_post_warmup_observations_per_chain"]
                   if n <= sampling["draws"]]
    use_checkpoints = sampling.get("checkpoint_policy", "final") == "first_qualified"
    if not use_checkpoints:
        checkpoints = [sampling["draws"]]
    elif not checkpoints or checkpoints[-1] != sampling["draws"]:
        checkpoints.append(sampling["draws"])
    segments, records, previous, diagnosis = [], [], 0, None
    sample_seconds, diagnostic_seconds, used_seeds = 0.0, 0.0, set()
    for segment_index, endpoint in enumerate(checkpoints):
        if time.perf_counter() - started >= budget_seconds:
            if not segments:
                raise TimeoutError("Runtime budget exhausted during input preparation")
            diagnosis["status"] = "accuracy_not_reached"
            diagnosis["checks"]["elapsed_budget"] = False
            break
        if segment_index == 0:
            seeds = sampling["seeds"]
        else:
            seeds = [int(np.random.SeedSequence(seed, spawn_key=(731, segment_index)).generate_state(1)[0])
                     for seed in sampling["seeds"]]
        if len(set(seeds)) != len(seeds) or used_seeds & set(seeds):
            raise ValueError("Deterministic segment seed collision")
        used_seeds.update(seeds)
        controls = {"seeds": seeds, "draws": endpoint - previous,
                    "warmup": sampling["warmup"] if segment_index == 0 else 0,
                    "sweeps_per_draw": sampling.get("sweeps_per_draw", 1),
                    "initial_states": None if not segments else segments[-1].states[:, -1]}
        begin = time.perf_counter()
        if sampling["method"] == "direct":
            segment = sample_direct(reduced, beta=beta, **controls)
        else:
            from .encoded import sample_encoded
            extras = {k: sampling[k] for k in ("backend", "penalty", "proposal_beta", "coefficient_scale",
                       "quantization_step", "max_spins", "max_edges", "max_bytes") if k in sampling}
            segment = sample_encoded(reduced, beta=beta, **controls, **extras)
        sample_seconds += time.perf_counter() - begin
        segments.append(segment)
        result = SamplingResult(np.concatenate([s.states for s in segments], axis=1),
            np.concatenate([s.energies for s in segments], axis=1), segments[0].initial_states,
            np.concatenate([s.accepted for s in segments], axis=1), dict(segment.metadata),
            {"segments": [s.timings for s in segments]})
        begin = time.perf_counter()
        diagnosis = diagnose(reduced, result, beta=beta, anchor=anchor, selected_pairs=selected_pairs,
                             evaluation=evaluation, reference=reference)
        diagnostic_seconds += time.perf_counter() - begin
        is_declared = endpoint in evaluation["timing"]["checkpoints_post_warmup_observations_per_chain"]
        records.append({"draws_per_chain": endpoint, "segment_seeds": seeds,
                        "status": diagnosis["status"], "checks": diagnosis["checks"],
                        "elapsed_through_diagnostics_seconds": time.perf_counter() - started,
                        "predeclared_checkpoint": is_declared})
        previous = endpoint
        if time.perf_counter() - started >= budget_seconds:
            diagnosis["status"] = "accuracy_not_reached"
            diagnosis["checks"]["elapsed_budget"] = False
            break
        if use_checkpoints and diagnosis["status"] == "accuracy_passed" and is_declared:
            break
    result.metadata["checkpoint_segments"] = records
    result.metadata["segment_randomness"] = "Independent deterministic substreams; state continues across boundaries; only the first segment warms up"
    result.metadata["draws_per_chain"] = result.states.shape[1]
    return result, diagnosis, sample_seconds, diagnostic_seconds, records


def run_workflow(request_path, output_path):
    start = time.perf_counter()
    request_path, output = Path(request_path).resolve(), Path(output_path).resolve()
    if output.exists():
        raise ValueError("Output path already exists; use a new run directory to preserve evidence")
    request, base = load_request(request_path)
    request_hash = request.source_sha256
    resolve = lambda name: (base / request[name]).resolve()
    benchmark_path = resolve("benchmark_config")
    evaluation_path = resolve("evaluation_config")
    evaluation = load_evaluation(evaluation_path)
    if evaluation.get("precision") != "float64" or evaluation.get("chains_per_condition") != 4:
        raise ValueError("This release requires the recorded float64 four chain evaluation protocol")
    inputs = load_benchmark(benchmark_path)
    domains = benchmark_domains(inputs, request["case_id"])
    domains, residue_rules = apply_residue_requests(inputs, domains, request.get("restrictions", []))
    beta = request.get("beta", inputs.config["model"]["beta"])
    if not math.isclose(beta * evaluation["temperature"], 1.0, rel_tol=1e-12):
        raise ValueError("Request beta and evaluation temperature specify different targets")
    reduced = restrict_model(inputs.model, domains)
    case = next(c for c in inputs.cases if c["case_id"] == request["case_id"])
    expected_split = ("development" if case["case_id"] in evaluation["development_case_ids"] else
                      "evaluation" if case["case_id"] in evaluation["evaluation_case_ids"] else None)
    if expected_split != case["split"]:
        raise ValueError("Case split is missing or inconsistent with the evaluation contract")
    sampling = dict(request["sampling"])
    code_hash = implementation_fingerprint()
    evaluation_hash = evaluation.source_sha256
    model_hash = inputs.config["archive"]["member_sha256"]
    target_hash = _target_fingerprint(inputs, domains, beta)
    comparison_context = _comparison_context(request, base)
    plan = None
    if sampling["method"] == "planner":
        plan = choose_method(resolve("planner"), model_sha256=model_hash, evaluation_sha256=evaluation_hash,
                             code_sha256=code_hash, case_kind=case["kind"], comparison_context_sha256=comparison_context)
        sampling["method"] = plan["method"]
        for key in ("backend", "penalty", "proposal_beta", "coefficient_scale", "quantization_step",
                    "sweeps_per_draw", "warmup", "draws", "checkpoint_policy"):
            if key in plan["sampling_parameters"]:
                sampling[key] = plan["sampling_parameters"][key]
            elif key in ("backend", "penalty", "proposal_beta", "coefficient_scale", "quantization_step"):
                sampling.pop(key, None)
    # Check traces before allocating. This is independent of encoded graph storage.
    estimated_trace_bytes = 4 * sampling["draws"] * (inputs.model.n_sites + reduced.n_sites + 3) * 8
    if estimated_trace_bytes > sampling.get("max_bytes", 536870912):
        raise ValueError("Requested trace exceeds the explicit host memory budget")
    reference = reference_provenance = None
    if "reference_run" in request:
        reference, reference_provenance = _reference(resolve("reference_run"), target_sha256=target_hash,
            current_seeds=sampling["seeds"], evaluation_sha256=evaluation_hash)
    selected_pairs = evaluation["selected_pair_observables"].get(case["case_id"], [])
    settings = dict(request["selection"])
    if "allowed_residues" in settings:
        raise ValueError("Selection support is inferred from original model fields and cannot be overridden")
    original_support = ["".join(inputs.model.alphabet[t] for t in np.flatnonzero(np.isfinite(row)))
                        for row in inputs.model.fields]
    settings["allowed_residues"] = original_support
    from .design import validate_design_inputs
    preflight = validate_design_inputs(anchor=inputs.anchor, alphabet=inputs.model.alphabet, settings=settings,
                                       assays_path=resolve("assays") if "assays" in request else None)
    encoded_report = _encoding_estimate(reduced)
    prep_seconds = time.perf_counter() - start
    budget_key = "reference_seconds_per_case_max" if request.get("purpose") == "reference" else "measurement_seconds_per_condition"
    budget = evaluation["budget"][budget_key]
    result, diagnostics, sample_seconds, diagnostics_seconds, checkpoint_records = _sample_with_checkpoints(
        reduced, sampling=sampling, beta=beta, anchor=inputs.anchor, selected_pairs=selected_pairs,
        evaluation=evaluation, reference=reference, started=start, budget_seconds=budget)
    candidate_start = time.perf_counter()
    full = reduced.expand(result.states)
    unique, first, inverse, counts = np.unique(full.reshape(-1, inputs.model.n_sites), axis=0,
                                             return_index=True, return_inverse=True, return_counts=True)
    original_energies = inputs.model.energy(unique)
    if not np.allclose(original_energies[inverse], result.energies.reshape(-1), rtol=1e-11, atol=1e-8):
        raise ValueError("Original model score verification failed for sampled states")
    for site, allowed in enumerate(domains):
        if not np.all(np.isin(unique[:, site], allowed)):
            raise ValueError("A sampled sequence violates an enforced residue restriction")
    pool = []
    for state, flat_index, count, energy in zip(unique, first, counts, original_energies):
        aligned = inputs.model.decode(state)
        sequence_hash = hashlib.sha256(aligned.encode()).hexdigest()
        pool.append({"candidate_id": "cm_" + sequence_hash[:16], "sequence_sha256": sequence_hash,
                     "aligned_sequence": aligned, "sequence": aligned.replace("-", ""),
                     "energy": float(energy), "source_chain": int(flat_index // result.states.shape[1]),
                     "source_draw": int(flat_index % result.states.shape[1]), "occurrences": int(count),
                     "hamming_to_anchor": sum(a != b for a, b in zip(aligned, inputs.anchor))})
    if len({r["candidate_id"] for r in pool}) != len(pool):
        raise ValueError("Candidate identifier collision")
    from .design import design_library
    if time.perf_counter() - start >= budget:
        design = {"status": "skipped_runtime_budget", "selected": [], "all_candidates": pool, "exclusions": [], "evidence": {}}
    else:
        design = design_library(pool, anchor=inputs.anchor, alphabet=inputs.model.alphabet, settings=settings,
                                assays_path=resolve("assays") if "assays" in request else None)
    if "assays" in request and design["status"] != "skipped_runtime_budget" and preflight["assays"]["sha256"] != design["evidence"]["assays"]["sha256"]:
        raise ValueError("Assay data changed between preflight validation and candidate selection")
    selection_seconds = time.perf_counter() - candidate_start
    allow_export = diagnostics["status"] == "accuracy_passed" or request.get("allow_unqualified_export", False)
    selected = design.get("selected", []) if allow_export else []
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent / ("." + output.name + ".partial." + uuid.uuid4().hex[:10])
    work.mkdir()
    export_start = time.perf_counter()
    save_json(work / "request.json", request)
    save_json(work / "evaluation.json", evaluation)
    save_json(work / "benchmark.json", inputs.config)
    save_json(work / "diagnostics.json", diagnostics)
    save_json(work / "checkpoints.json", {"policy": sampling.get("checkpoint_policy", "final"), "records": checkpoint_records})
    save_json(work / "selection.json", design)
    save_json(work / "encoding_resources.json", encoded_report)
    np.savez_compressed(work / "trace.npz", states=result.states, full_states=full, energies=result.energies,
                        initial_states=result.initial_states, sampler_event=result.accepted)
    with gzip.open(work / "trace.csv.gz", "wt", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["chain", "draw", "aligned_sequence", "energy", "sampler_event", "state_changed_from_previous_observation"])
        for c in range(full.shape[0]):
            for d in range(full.shape[1]):
                changed = "" if d == 0 else bool(np.any(full[c, d] != full[c, d - 1]))
                writer.writerow([c, d, inputs.model.decode(full[c, d]), float(result.energies[c, d]),
                                 bool(result.accepted[c, d]), changed])
    fields = ["candidate_id", "sequence_sha256", "aligned_sequence", "sequence", "energy", "source_chain", "source_draw",
              "occurrences", "hamming_to_anchor", "utility", "predicted_mean", "predictive_std",
              "selection_rank", "selection_reason", "uncertainty_kind", "exclusion_reasons"]
    _csv(work / "candidate_pool.csv", design.get("all_candidates", pool), fields)
    _csv(work / "candidates.csv", selected, fields)
    _csv(work / "position_map.csv", inputs.mapping,
         ["model_index_0based", "alignment_column_1based", "anchor_symbol", "uniprot_position_1based"])
    for filename, field in (("candidates.fasta", "sequence"), ("candidates_aligned.fasta", "aligned_sequence")):
        with (work / filename).open("w", encoding="ascii") as stream:
            for row in selected:
                stream.write(f">{row['candidate_id']} source_chain={row['source_chain']} source_draw={row['source_draw']}\n{row[field]}\n")
    assay_fields = ["candidate_id", "aligned_sequence", "value", "split", "assay_id", "units", "context", "round"]
    _csv(work / "assay_worksheet.csv", [{"candidate_id": row["candidate_id"], "aligned_sequence": row["aligned_sequence"],
        "value": "", "split": "train", "assay_id": settings.get("assay_id", ""), "units": settings.get("units", ""),
        "context": settings.get("context", ""), "round": ""} for row in selected], assay_fields)
    rules = {"schema_version": 1, "original_parameter_hashes": "verified", "complete_parameter_coverage": "verified",
             "forbidden_original_fields": "exactly_excluded", "biological_coordinate_map": "verified",
             "base_case_constraints": "enforced", "additional_residue_rules": residue_rules,
             "allowed_symbols_by_model_index": ["".join(inputs.model.alphabet[t] for t in d) for d in domains],
             "all_sampled_original_scores": "verified", "all_sampled_residue_constraints": "verified",
             "sampling_repeats": "retained", "selection_rules": request["selection"],
             "selection_changes_sampling_distribution": True, "biological_function_verified": False,
             "independent_assay_holdout_for_potts_model": False,
             "physical_device_execution": False, "human_review_status": "Prior review records preserved; no new adviser review asserted"}
    save_json(work / "rules.json", rules)
    export_seconds = time.perf_counter() - export_start
    elapsed = time.perf_counter() - start
    budget_exceeded = elapsed > budget
    if budget_exceeded:
        diagnostics["status"] = "accuracy_not_reached"
        diagnostics["checks"]["elapsed_budget"] = False
        save_json(work / "diagnostics.json", diagnostics)
        if not request.get("allow_unqualified_export", False):
            selected = []
            _csv(work / "candidates.csv", [], fields)
            _csv(work / "assay_worksheet.csv", [], assay_fields)
            (work / "candidates.fasta").write_text("", encoding="ascii")
            (work / "candidates_aligned.fasta").write_text("", encoding="ascii")
    checkpoints = evaluation["timing"]["checkpoints_post_warmup_observations_per_chain"]
    measured_checkpoint = result.states.shape[1] in checkpoints
    eligible_time = diagnostics["status"] == "accuracy_passed" and measured_checkpoint and not budget_exceeded
    profile = {"load_validate_reduce_seconds": prep_seconds, "sampler_total_seconds": sample_seconds,
               "sampler_breakdown": result.timings, "diagnostics_seconds": diagnostics_seconds,
               "score_verify_select_seconds": selection_seconds, "export_seconds": export_seconds,
               "elapsed_before_final_manifest_seconds": elapsed,
               "optimizations": {"fixed_sites_removed": inputs.model.n_sites - reduced.n_sites,
                    "full_categorical_sites": inputs.model.n_sites, "mutable_categorical_sites": reduced.n_sites,
                    "full_pair_blocks": inputs.model.n_sites * (inputs.model.n_sites - 1) // 2,
                    "reduced_pair_blocks": reduced.n_sites * (reduced.n_sites - 1) // 2,
                    "original_scores_reused_for_repeats": full.shape[0] * full.shape[1] - len(pool)},
               "comparison_limit": "Time to accuracy is established only for first_qualified runs; neither single nor unmatched runs establish speedup."}
    save_json(work / "profile.json", profile)
    status = ("qualified_library" if diagnostics["status"] == "accuracy_passed" and selected else
              "exploratory_library" if selected else "library_not_exported")
    if request.get("purpose") == "reference":
        status = "reference_qualified" if diagnostics["reference_qualification"]["qualified"] and not budget_exceeded else "reference_inadequate"
    manifest = {"schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status, "purpose": request.get("purpose", "candidate_library"), "case": case,
        "model_sha256": model_hash, "target_sha256": target_hash,
        "evaluation_sha256": evaluation_hash, "code_sha256": code_hash, "request_sha256": request_hash,
        "comparison_context_sha256": comparison_context,
        "source_provenance": inputs.provenance, "beta": beta, "precision": "float64",
        "method": sampling["method"], "sampling_parameters": {k: v for k, v in sampling.items() if k not in ("method", "seeds")},
        "seeds": sampling["seeds"], "sampler_metadata": result.metadata, "planner_decision": plan,
        "runtime": _runtime(), "reference": reference_provenance,
        "physical_execution": False, "elapsed_seconds": elapsed, "budget_seconds": budget,
        "budget_exceeded": budget_exceeded, "predeclared_checkpoint": measured_checkpoint,
        "qualified_checkpoint_elapsed_seconds": elapsed if eligible_time else None,
        "first_time_to_accuracy_established": eligible_time and sampling.get("checkpoint_policy") == "first_qualified",
        "raw_observations": int(full.shape[0] * full.shape[1]),
        "unique_candidates": len(pool), "exported_candidates": len(selected),
        "export_permission": "explicit_exploratory" if request.get("allow_unqualified_export", False) else "qualified_only",
        "assay_source": {"path": str(resolve("assays")), "sha256": preflight["assays"]["sha256"]} if "assays" in request else None,
        "sequence_scope": {"accession": inputs.config["coordinate_reference"]["uniprot_accession"],
                           "span": inputs.config["coordinate_reference"]["reference_span_in_full_accession"],
                           "description": "Modeled domain segment, not the full accession construct"},
        "claim_limits": {"improved_function": False, "hardware_speedup": False, "independent_potts_assay_validation": False}}
    report = (f"# Protein candidate run\n\nStatus: {status}.\n\n"
              f"Case: `{case['case_id']}`. Method: `{sampling['method']}`. "
              f"Recorded {manifest['raw_observations']} observations in four chains and {len(pool)} distinct sequences.\n\n"
              f"Sampling assessment: `{diagnostics['status']}`. Exported library: {len(selected)} candidates. "
              f"Complete measured work before the final manifest took {elapsed:.3f} seconds.\n\n"
              "Every candidate is linked to its original chain and draw. Raw traces retain repeated states. "
              "The library applies separate selection rules and does not represent equilibrium frequencies.\n\n"
              "The exported sequence is the modeled CM segment corresponding to Q5SLA5 residues 2 through 90. "
              "Its model score is not measured activity or stability. No physical Extropic device was used.\n\n"
              "Inspect `rules.json`, `diagnostics.json`, `selection.json` and `profile.json` for the applied rules, "
              "failed checks, objective and timing. Fill measured values and experimental context in "
              "`assay_worksheet.csv` before ingestion. Blank measurements are never replaced by predictions.\n")
    (work / "report.md").write_text(report, encoding="utf-8")
    manifest["artifacts"] = {p.name: digest(p) for p in sorted(work.iterdir()) if p.is_file()}
    save_json(work / "manifest.json", manifest)
    work.rename(output)
    return manifest


def verify_run(directory):
    directory = Path(directory)
    manifest = read_json(directory / "manifest.json")
    failures = [name for name, expected in manifest["artifacts"].items()
                if not (directory / name).is_file() or digest(directory / name) != expected]
    if failures:
        raise ValueError(f"Run artifact verification failed: {failures}")
    with np.load(directory / "trace.npz", allow_pickle=False) as trace:
        if trace["states"].shape[:2] != trace["energies"].shape:
            raise ValueError("Trace shapes do not agree")
        if trace["energies"].size != manifest["raw_observations"] or not np.all(np.isfinite(trace["energies"])):
            raise ValueError("Trace observation count or finite energies failed validation")
    return {"status": "verified", "artifacts": len(manifest["artifacts"]),
            "raw_observations": manifest["raw_observations"]}
