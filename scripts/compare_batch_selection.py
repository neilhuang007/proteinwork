"""Compare the saved THRML library objective against MILP on the identical pool."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from proteinwork.design import design_library
from proteinwork.validation import digest, read_json, save_json
from proteinwork.workflow import verify_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Comparison output already exists")
    verify_run(args.run)
    saved = read_json(args.run / "selection.json")
    manifest = read_json(args.run / "manifest.json")
    if saved["evidence"]["objective"] != "model_score":
        raise ValueError("This replay requires model_score; assay replay also needs the original hashed assay file")
    settings = dict(saved["evidence"]["settings"])
    settings["optimizer"] = "milp"
    reference = manifest["source_provenance"]["reference"]
    benchmark = read_json(args.run / "benchmark.json")
    baseline = design_library(saved["all_candidates"], anchor=reference["aligned_sequence"],
                              alphabet=benchmark["model"]["alphabet"], settings=settings)
    if baseline["evidence"]["optimization_pool_ids"] != saved["evidence"]["optimization_pool_ids"]:
        raise ValueError("Replay changed the admitted pool")
    result = {"schema_version": 1, "source_manifest_sha256": digest(args.run / "manifest.json"),
              "same_pool": True, "pool_size": len(saved["evidence"]["optimization_pool_ids"]),
              "original_optimizer": saved["evidence"]["settings"]["optimizer"],
              "original_solver": saved["evidence"]["solver"],
              "original_selected_objective": saved["evidence"]["selected_objective"],
              "milp_status": baseline["status"], "milp_solver": baseline["evidence"].get("solver"),
              "milp_selected_objective": baseline["evidence"].get("selected_objective"),
              "milp_candidates": [row["candidate_id"] for row in baseline["selected"]],
              "claim_limit": "Model score and diversity objective only; no assay result or hardware advantage"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    save_json(args.output, result)
    print(f"Saved identical pool comparison to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
