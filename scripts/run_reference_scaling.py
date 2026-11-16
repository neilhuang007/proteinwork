"""Verify reference export equivalence, then run a fixed larger reference attempt."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from proteinwork.validation import digest, read_json, save_json
from proteinwork.workflow import _runtime, implementation_fingerprint, verify_run
from scripts.build_review_packet import diagnostic_audit, write_csv


def verify_replay(baseline, replay):
    baseline, replay = Path(baseline), Path(replay)
    baseline_manifest = read_json(baseline / "manifest.json")
    replay_manifest = read_json(replay / "manifest.json")
    for key in ("target_sha256", "evaluation_sha256", "sampling_parameters", "seeds", "raw_observations", "unique_candidates"):
        if baseline_manifest[key] != replay_manifest[key]:
            raise ValueError(f"Replay differs in {key}")
    with np.load(baseline / "trace.npz", allow_pickle=False) as before, np.load(replay / "trace.npz", allow_pickle=False) as after:
        if set(before.files) != set(after.files):
            raise ValueError("Replay trace fields differ")
        for key in before.files:
            np.testing.assert_array_equal(before[key], after[key], err_msg=f"Replay changed {key}")
    if read_json(baseline / "diagnostics.json") != read_json(replay / "diagnostics.json"):
        raise ValueError("Replay changed diagnostic evidence")
    before_bytes = sum(path.stat().st_size for path in baseline.iterdir() if path.is_file())
    after_bytes = sum(path.stat().st_size for path in replay.iterdir() if path.is_file())
    return {"all_trace_arrays_identical": True, "diagnostics_identical": True,
            "baseline_artifact_bytes": before_bytes, "replay_artifact_bytes": after_bytes,
            "artifact_bytes_saved": before_bytes - after_bytes,
            "artifact_reduction_fraction": 1 - after_bytes / before_bytes,
            "baseline_manifest_sha256": digest(baseline / "manifest.json"),
            "replay_manifest_sha256": digest(replay / "manifest.json")}


def audit_reference(directory):
    directory = Path(directory)
    verify_run(directory)
    manifest = read_json(directory / "manifest.json")
    evaluation = read_json(directory / "evaluation.json")
    summary, blockers = diagnostic_audit(read_json(directory / "diagnostics.json"), evaluation)
    with (directory / "position_map.csv").open(encoding="utf-8", newline="") as stream:
        mapping = list(csv.DictReader(stream))
    with np.load(directory / "trace.npz", allow_pickle=False) as trace:
        full_states = trace["full_states"]
        for row in blockers:
            parts = row["key"].split(":")
            if parts[0] == "site":
                sites, tokens = [int(parts[1])], [int(parts[2])]
            else:
                sites, tokens = list(map(int, parts[1:3])), list(map(int, parts[3:5]))
            indicator = np.ones(full_states.shape[:2], dtype=bool)
            for site, token in zip(sites, tokens, strict=True):
                indicator &= full_states[..., site] == token
            if not np.isclose(indicator.mean(), row["probability"], atol=1e-15, rtol=0):
                raise ValueError("Diagnostic frequency disagrees with the retained trace")
            row["accession_positions"] = ";".join(mapping[site]["uniprot_position_1based"] for site in sites)
            row.update({f"chain_{chain}_count": int(count) for chain, count in enumerate(indicator.sum(axis=1))})
    summary.update(run_status=manifest["status"], raw_observations=manifest["raw_observations"],
                   draws_per_chain=manifest["sampler_metadata"]["draws_per_chain"],
                   recorded_workflow_seconds=manifest["elapsed_seconds"],
                   budget_exceeded=manifest["budget_exceeded"],
                   artifact_bytes=sum(path.stat().st_size for path in directory.iterdir() if path.is_file()),
                   manifest_sha256=digest(directory / "manifest.json"))
    return summary, blockers


def run_study(baseline, output):
    baseline, output = Path(baseline).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError("Use a new study directory")
    verify_run(baseline)
    baseline_manifest = read_json(baseline / "manifest.json")
    if baseline_manifest["purpose"] != "reference" or baseline_manifest["method"] != "direct":
        raise ValueError("The baseline must be an independent direct reference run")
    source_request = read_json(ROOT / "configs/cm_reference.json")
    if source_request != read_json(baseline / "request.json"):
        raise ValueError("The saved baseline request differs from the declared reference configuration")
    if (baseline_manifest["sampling_parameters"]["draws"] != 16384
            or baseline_manifest["seeds"] != source_request["sampling"]["seeds"]):
        raise ValueError("This fixed study replays the declared 16,384 draw reference baseline")
    output.mkdir(parents=True)
    requests = output / "requests"
    requests.mkdir()
    runtime = _runtime()
    environment = {**os.environ, "PROCESSOR_IDENTIFIER": runtime["processor"]}
    code_hash = implementation_fingerprint()
    jobs = []
    for name, draws, seeds in (("reference_replay", 16384, baseline_manifest["seeds"]),
                              ("reference_32768", 32768, [9501, 9502, 9503, 9504])):
        request = read_json(ROOT / "configs/cm_reference.json")
        request["sampling"].update(draws=draws, seeds=seeds)
        for key in ("benchmark_config", "evaluation_config"):
            request[key] = Path(os.path.relpath(ROOT / "configs" / request[key], requests)).as_posix()
        path = requests / f"{name}.json"
        save_json(path, request)
        jobs.append({"name": name, "draws_per_chain": draws, "seeds": seeds,
                     "raw_observations": 4 * draws, "request_sha256": digest(path),
                     "trace_memory_estimate_bytes": 4 * draws * (96 + 88 + 3) * 8})
    save_json(output / "protocol.json", {
        "schema_version": 1, "declared_at_utc": datetime.now(timezone.utc).isoformat(),
        "baseline_manifest_sha256": digest(baseline / "manifest.json"),
        "code_sha256": code_hash, "runner_sha256": digest(Path(__file__)),
        "runtime_snapshot": runtime, "jobs": jobs,
        "policy": "Run sequentially. Confirm replay equality before the larger run. Preserve the model, restrictions, beta, warmup, sweep count and all accuracy thresholds. No tuning or candidate selection.",
        "resource_policy": "Keep the existing 512 MiB trace allocation guard and 3,600 second CLI cap per run. The trace estimate is not a peak process memory guarantee. Doubling observations is a bounded development attempt, not a promise of qualification.",
        "stop_rule": "Stop after the declared 32,768 draw attempt, whether it qualifies or fails. Retain failures and every raw observation. Do not invent uncertainty for unobserved indicators.",
        "timing_policy": "Record complete CLI process wall time separately from the workflow timer. Replay storage equality is the primary optimization check; single runs on different dates do not establish a runtime speedup."})
    execution = []
    for job in jobs:
        if implementation_fingerprint() != code_hash:
            raise ValueError("Implementation changed after protocol declaration")
        started = time.perf_counter()
        with (output / f"{job['name']}.log").open("w", encoding="utf-8") as log:
            result = subprocess.run([sys.executable, "-m", "proteinwork", "run", "--request",
                                     str(requests / f"{job['name']}.json"), "--output", str(output / job["name"])],
                                    cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT, check=False)
        execution.append({"name": job["name"], "exit_code": result.returncode,
                          "complete_cli_wall_seconds": time.perf_counter() - started})
        save_json(output / "execution.json", {"runs": execution})
        if result.returncode:
            raise RuntimeError(f"Run {job['name']} failed; inspect its preserved log and failure evidence")
        verify_run(output / job["name"])
        print(f"Completed {job['name']}", flush=True)
        if job["name"] == "reference_replay":
            save_json(output / "replay_equivalence.json", verify_replay(baseline, output / job["name"]))
            print("Replay arrays and diagnostics are identical", flush=True)
    records = {}
    for name, directory in (("baseline", baseline), *((job["name"], output / job["name"]) for job in jobs)):
        records[name], blockers = audit_reference(directory)
        write_csv(output / f"{name}_frequency_blockers.csv", blockers,
                  ["key", "kind", "model_indices", "symbols", "probability", "mcse", "reference_mcse_limit", "reason",
                   "accession_positions", "chain_0_count", "chain_1_count", "chain_2_count", "chain_3_count"])
    save_json(output / "summary.json", {"runs": records, "replay_equivalence": read_json(output / "replay_equivalence.json"),
              "claim_limit": "Reference qualification uses unchanged criteria. Passing general chain diagnostics alone does not qualify a reference or establish biological function."})
    print(f"Larger reference status: {records['reference_32768']['run_status']}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=ROOT / "output/week07/cm_reference")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run_study(args.baseline, args.output)
