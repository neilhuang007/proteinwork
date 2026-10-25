"""Freeze method choice using qualified development runs only."""

from __future__ import annotations

from collections import defaultdict
import statistics
from pathlib import Path

from .validation import digest, object_digest, read_json, save_json


def freeze_plan(run_directories, output):
    records, grouped, identities, seen_evidence = [], defaultdict(list), set(), set()
    for directory in map(Path, run_directories):
        from .workflow import verify_run
        verify_run(directory)
        evidence_hash = digest(directory / "manifest.json")
        if evidence_hash in seen_evidence:
            raise ValueError("Duplicate run evidence cannot count as a repeated measurement")
        seen_evidence.add(evidence_hash)
        manifest = read_json(directory / "manifest.json")
        diagnostics = read_json(directory / "diagnostics.json")
        if manifest["case"]["split"] != "development":
            raise ValueError("Evaluation runs cannot train the method planner")
        if digest(directory / "diagnostics.json") != manifest["artifacts"]["diagnostics.json"]:
            raise ValueError("Development diagnostic artifact hash mismatch")
        identities.add((manifest["model_sha256"], manifest["evaluation_sha256"], manifest["code_sha256"],
                        manifest["comparison_context_sha256"]))
        record = {"manifest_path": str((directory / "manifest.json").resolve()),
                  "manifest_sha256": digest(directory / "manifest.json"),
                  "case_id": manifest["case"]["case_id"], "status": diagnostics["status"],
                  "target_sha256": manifest["target_sha256"],
                  "method": manifest["method"], "sampling_parameters": manifest["sampling_parameters"],
                  "kind": manifest["case"]["kind"], "seconds": manifest["elapsed_seconds"],
                  "budget_exceeded": manifest["budget_exceeded"]}
        records.append(record)
        if diagnostics["status"] == "accuracy_passed" and not manifest["budget_exceeded"]:
            signature = object_digest(record["sampling_parameters"])
            grouped[(record["kind"], record["method"], signature)].append(record)
    if not records or len(identities) != 1:
        raise ValueError("Planner evidence requires compatible model, evaluation and implementation fingerprints")
    candidates = defaultdict(list)
    for (kind, _, _), rows in grouped.items():
        expected_targets = {r["target_sha256"] for r in records if r["kind"] == kind}
        if {r["target_sha256"] for r in rows} != expected_targets:
            continue
        per_target_seconds = [statistics.median(r["seconds"] for r in rows if r["target_sha256"] == target)
                              for target in sorted(expected_targets)]
        candidates[kind].append({"method": rows[0]["method"],
                                 "sampling_parameters": rows[0]["sampling_parameters"],
                                 "median_seconds": statistics.median(per_target_seconds),
                                 "matched_targets": len(expected_targets),
                                 "measurements": len(rows)})
    choices = {kind: min(rows, key=lambda r: (r["median_seconds"], r["method"]))
               for kind, rows in candidates.items()}
    model_hash, evaluation_hash, code_hash, context_hash = identities.pop()
    plan = {"schema_version": 1, "frozen": True, "model_sha256": model_hash,
            "evaluation_sha256": evaluation_hash, "code_sha256": code_hash,
            "comparison_context_sha256": context_hash,
            "evidence": records, "choices": choices, "fallback": "direct",
            "scope": "Selection by declared case kind from development evidence; no within run adaptation.",
            "limitations": "A frozen selection rule does not establish superiority or generalization."}
    plan["plan_sha256"] = object_digest(plan)
    output = Path(output)
    if output.exists():
        raise ValueError("A frozen plan cannot overwrite an existing artifact")
    output.parent.mkdir(parents=True, exist_ok=True)
    save_json(output, plan)
    return plan


def choose_method(path, *, model_sha256, evaluation_sha256, code_sha256, case_kind, comparison_context_sha256):
    plan = read_json(path)
    saved_hash = plan.get("plan_sha256")
    if saved_hash != object_digest({k: v for k, v in plan.items() if k != "plan_sha256"}):
        raise ValueError("Frozen planner fingerprint mismatch")
    for name, expected in (("model_sha256", model_sha256), ("evaluation_sha256", evaluation_sha256),
                           ("code_sha256", code_sha256), ("comparison_context_sha256", comparison_context_sha256)):
        if plan.get(name) != expected:
            raise ValueError(f"Planner {name} does not match the current experiment")
    if plan.get("frozen") is not True or plan.get("schema_version") != 1:
        raise ValueError("Planner is not a supported frozen plan")
    selected = plan["choices"].get(case_kind)
    if selected is None:
        return {"method": "direct", "sampling_parameters": {}, "reason": "No qualified development evidence for this case kind",
                "plan_sha256": saved_hash}
    return {**selected, "reason": "Lowest median complete runtime among qualified development measurements",
            "plan_sha256": saved_hash}
