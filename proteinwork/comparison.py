"""Compare matched completed runs without promoting failures into speed claims."""

from __future__ import annotations

from pathlib import Path
from collections import defaultdict
import statistics

from .validation import digest, object_digest, read_json, save_json
from .workflow import verify_run


def _summarize_repetitions(records):
    groups = defaultdict(list)
    for record in records:
        if "method" in record:
            signature = object_digest(record["sampling_parameters"])
            groups[(record["method"], signature)].append(record)
    summaries = []
    for (method, signature), runs in sorted(groups.items()):
        qualified_times = [run["elapsed_seconds"] for run in runs if run["eligible_for_timing_comparison"]]
        seed_sets = [set(run["seeds"]) for run in runs]
        independent_seeds = all(not left & right for index, left in enumerate(seed_sets)
                                for right in seed_sets[index + 1:])
        summaries.append({"method": method, "sampling_parameters_sha256": signature,
                          "sampling_parameters": runs[0]["sampling_parameters"],
                          "runs": len(runs), "qualified_runs": len(qualified_times),
                          "unqualified_runs": len(runs) - len(qualified_times),
                          "disjoint_seeds_within_method": independent_seeds,
                          "median_qualified_seconds": statistics.median(qualified_times) if qualified_times else None,
                          "minimum_qualified_seconds": min(qualified_times, default=None),
                          "maximum_qualified_seconds": max(qualified_times, default=None)})
    return summaries


def compare_runs(directories, output):
    output = Path(output)
    if output.exists():
        raise ValueError("Comparison output already exists")
    records, targets, contexts, seen_evidence = [], set(), set(), set()
    for directory in map(Path, directories):
        failure = directory / "failure.json"
        evidence_hash = digest(failure if failure.exists() else directory / "manifest.json")
        if evidence_hash in seen_evidence:
            raise ValueError("Duplicate run evidence cannot count as a repeated measurement")
        seen_evidence.add(evidence_hash)
        if failure.exists():
            records.append({"run": str(directory.resolve()), "status": "failed", "failure": read_json(failure),
                            "eligible_for_timing_comparison": False})
            continue
        verify_run(directory)
        manifest = read_json(directory / "manifest.json")
        diagnostic = read_json(directory / "diagnostics.json")
        targets.add((manifest["target_sha256"], manifest["evaluation_sha256"], manifest["code_sha256"]))
        contexts.add(manifest["comparison_context_sha256"])
        records.append({"run": str(directory.resolve()), "method": manifest["method"],
                        "manifest_sha256": evidence_hash,
                        "seeds": manifest["seeds"], "sampling_parameters": manifest["sampling_parameters"],
                        "backend": manifest["sampler_metadata"].get("backend"), "status": diagnostic["status"],
                        "elapsed_seconds": manifest["elapsed_seconds"], "raw_observations": manifest["raw_observations"],
                        "exported_candidates": manifest["exported_candidates"],
                        "eligible_for_timing_comparison": (diagnostic["status"] == "accuracy_passed" and
                            manifest["first_time_to_accuracy_established"] and not manifest["budget_exceeded"])})
    if len(records) < 2 or len(targets) > 1 or len(contexts) > 1:
        raise ValueError("Comparison requires at least two runs with matching target, evaluation, code, hardware and selection objective")
    valid = [row for row in records if row["eligible_for_timing_comparison"]]
    baseline = min((row["elapsed_seconds"] for row in valid if row["method"] == "direct"), default=None)
    for row in valid:
        row["observed_time_ratio_to_fastest_direct"] = (baseline / row["elapsed_seconds"] if baseline is not None else None)
    result = {"schema_version": 1, "runs": records,
              "timing_qualified_runs": len(valid),
              "repeated_measurements": _summarize_repetitions(records),
              "failed_runs": sum(row["status"] == "failed" for row in records),
              "claims": "Observed ratios describe these runs only. Repetition summaries group identical method settings and report qualified times only, with failures retained separately. Shared seeds are flagged. Timing order and machine contention must be controlled by the experiment. No hardware or protein function improvement is established."}
    output.parent.mkdir(parents=True, exist_ok=True)
    save_json(output, result)
    return result
