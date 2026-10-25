"""Compare matched completed runs without promoting failures into speed claims."""

from __future__ import annotations

from pathlib import Path

from .validation import object_digest, read_json, save_json
from .workflow import verify_run


def compare_runs(directories, output):
    records, targets, contexts = [], set(), set()
    for directory in map(Path, directories):
        failure = directory / "failure.json"
        if failure.exists():
            records.append({"run": str(directory.resolve()), "status": "failed", "failure": read_json(failure),
                            "eligible_for_timing_comparison": False})
            continue
        verify_run(directory)
        manifest = read_json(directory / "manifest.json")
        request = read_json(directory / "request.json")
        diagnostic = read_json(directory / "diagnostics.json")
        targets.add((manifest["target_sha256"], manifest["evaluation_sha256"], manifest["code_sha256"]))
        contexts.add(manifest["comparison_context_sha256"])
        records.append({"run": str(directory.resolve()), "method": manifest["method"],
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
              "claims": "Observed ratios describe these runs only. Repeated sequential measurements are needed for a performance conclusion; no hardware or protein function improvement is established."}
    output = Path(output)
    if output.exists():
        raise ValueError("Comparison output already exists")
    output.parent.mkdir(parents=True, exist_ok=True)
    save_json(output, result)
    return result
