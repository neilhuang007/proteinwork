"""Command line entry point for reproducible research runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

from .validation import digest, load_evaluation, load_request, save_json


def main():
    parser = argparse.ArgumentParser(description="Validate, sample and select protein candidate libraries")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "_run"):
        command = commands.add_parser(name, help="Execute a research request" if name == "run" else argparse.SUPPRESS)
        command.add_argument("--request", type=Path, required=True)
        command.add_argument("--output", type=Path, required=True)
    validate = commands.add_parser("validate", help="Validate a complete request before sampling")
    validate.add_argument("--request", type=Path, required=True)
    verify = commands.add_parser("verify", help="Verify saved evidence hashes and trace integrity")
    verify.add_argument("run", type=Path)
    plan = commands.add_parser("freeze-plan", help="Freeze a method choice from qualified development runs")
    plan.add_argument("--runs", type=Path, nargs="+", required=True)
    plan.add_argument("--output", type=Path, required=True)
    ingest = commands.add_parser("ingest", help="Join supplied assay values to exported candidate identities")
    ingest.add_argument("--run", type=Path, required=True)
    ingest.add_argument("--measurements", type=Path, required=True)
    ingest.add_argument("--output", type=Path, required=True)
    ingest.add_argument("--existing", type=Path)
    comparison = commands.add_parser("compare", help="Compare matched runs at their first qualified checkpoint")
    comparison.add_argument("--runs", type=Path, nargs="+", required=True)
    comparison.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "run":
            if args.output.exists():
                raise ValueError("Output path already exists; use a new run directory")
            request, base = load_request(args.request)
            evaluation = load_evaluation(base / request["evaluation_config"])
            budget_key = "reference_seconds_per_case_max" if request.get("purpose") == "reference" else "measurement_seconds_per_condition"
            cap = float(evaluation["budget"][budget_key])
            if not 0 < cap <= 86400:
                raise ValueError("Runtime budget must be positive and at most one day")
            try:
                return subprocess.run([sys.executable, "-m", "proteinwork", "_run", "--request", str(args.request.resolve()),
                                       "--output", str(args.output.resolve())], timeout=cap, check=False).returncode
            except subprocess.TimeoutExpired:
                args.output.mkdir(parents=True, exist_ok=True)
                save_json(args.output / "failure.json", {"status": "accuracy_not_reached", "reason": "wall_clock_budget_exhausted",
                    "budget_seconds": cap, "request": str(args.request.resolve()), "complete_run": False})
                print("Run reached its wall clock budget; failure evidence was saved", file=sys.stderr)
                return 3
        if args.command == "_run":
            from .workflow import run_workflow
            result = run_workflow(args.request, args.output)
            print(json.dumps({k: result[k] for k in ("status", "method", "raw_observations", "unique_candidates", "exported_candidates", "elapsed_seconds")}, indent=2))
        elif args.command == "validate":
            from .workflow import validate_request
            print(json.dumps(validate_request(args.request), indent=2))
        elif args.command == "verify":
            from .workflow import verify_run
            print(json.dumps(verify_run(args.run), indent=2))
        elif args.command == "freeze-plan":
            from .planner import freeze_plan
            result = freeze_plan(args.runs, args.output)
            print(json.dumps({"status": "frozen", "choices": result["choices"], "plan_sha256": result["plan_sha256"]}, indent=2))
        elif args.command == "ingest":
            from .ingest import ingest_measurements
            print(json.dumps(ingest_measurements(args.run, args.measurements, args.output, args.existing), indent=2))
        elif args.command == "compare":
            from .comparison import compare_runs
            print(json.dumps(compare_runs(args.runs, args.output), indent=2))
        return 0
    except (ValueError, OSError, ImportError, KeyError) as error:
        if args.command == "_run" and not args.output.exists():
            args.output.mkdir(parents=True, exist_ok=True)
            save_json(args.output / "failure.json", {"status": "execution_failed", "error_type": type(error).__name__,
                "reason": str(error), "request_sha256": digest(args.request) if args.request.is_file() else None,
                "complete_run": False})
        print(f"Input or execution error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
