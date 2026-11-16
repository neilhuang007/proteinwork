"""Run three fixed development repetitions per method in sequential order."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from proteinwork.comparison import compare_runs
from proteinwork.validation import digest, read_json, save_json
from proteinwork.workflow import _runtime, implementation_fingerprint


def run_repetitions(output):
    output = Path(output).resolve()
    if output.exists():
        raise ValueError("Use a new experiment directory")
    output.mkdir(parents=True)
    request_directory = output / "requests"
    request_directory.mkdir()
    code_hash = implementation_fingerprint()
    runtime = _runtime()
    # Windows processor probing can vary between fresh child processes. Capture
    # it once and pass the same observed identifier to every measured process.
    environment = {**os.environ, "PROCESSOR_IDENTIFIER": runtime["processor"]}
    schedule = []
    for repetition in range(1, 4):
        methods = ("direct", "thrml") if repetition % 2 else ("thrml", "direct")
        seeds = [3000 + 100 * repetition + chain for chain in range(1, 5)]
        for method in methods:
            name = f"repeat{repetition}_{method}"
            request = dict(read_json(ROOT / "configs" / f"cm_exact_{method}.json"))
            for key in ("benchmark_config", "evaluation_config"):
                request[key] = Path(os.path.relpath(ROOT / "configs" / request[key], request_directory)).as_posix()
            request["sampling"]["seeds"] = seeds
            request_path = request_directory / f"{name}.json"
            save_json(request_path, request)
            schedule.append({"name": name, "method": method, "repetition": repetition,
                             "seeds": seeds, "request_sha256": digest(request_path)})
    save_json(output / "protocol.json", {
        "schema_version": 1, "declared_at_utc": datetime.now(timezone.utc).isoformat(),
        "code_sha256": code_hash, "runner_sha256": digest(Path(__file__)), "schedule": schedule,
        "runtime_snapshot": runtime,
        "processor_identity_policy": "Capture the observed processor identifier once and pass it unchanged to every subprocess.",
        "case_id": "cm_dev_exact_8_12", "purpose": "development_repetition",
        "execution_policy": "Fresh sequential CLI processes. Alternate method order between repetitions. Do not run other sampling or test jobs concurrently.",
        "settings_policy": "Existing exact development requests with distinct seeds per repetition. No parameter search or threshold changes.",
        "reporting_policy": "Retain every result and failure. Report median and range of qualified complete run times, never count failures as successful capped times.",
        "scope": "Small conditional development case only. Not a final held out comparison or a general performance claim."})
    run_directories = []
    executions = []
    for job in schedule:
        if implementation_fingerprint() != code_hash:
            raise ValueError("Implementation changed after the experiment protocol was saved")
        run_directory = output / job["name"]
        with (output / f"{job['name']}.log").open("w", encoding="utf-8") as log:
            result = subprocess.run([sys.executable, "-m", "proteinwork", "run", "--request",
                                     str(request_directory / f"{job['name']}.json"), "--output", str(run_directory)],
                                    cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT, check=False)
        executions.append({"name": job["name"], "returncode": result.returncode})
        save_json(output / "execution.json", {"runs": executions})
        if not (run_directory / "manifest.json").is_file() and not (run_directory / "failure.json").is_file():
            raise RuntimeError(f"Run {job['name']} has no result artifact; inspect its log")
        run_directories.append(run_directory)
        print(f"Completed {job['name']}, exit code {result.returncode}", flush=True)
    try:
        comparison = compare_runs(run_directories, output / "comparison.json")
    except ValueError as error:
        save_json(output / "comparison_failure.json", {"status": "comparison_rejected", "reason": str(error)})
        raise
    print(f"Qualified timing measurements: {comparison['timing_qualified_runs']} of {len(schedule)}", flush=True)
    return comparison


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run_repetitions(parser.parse_args().output)
