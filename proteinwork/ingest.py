"""Join laboratory measurements to immutable exported candidate identities."""

from __future__ import annotations

import csv
import hashlib
import io
import math
from pathlib import Path

from .validation import digest, read_json, save_json
from .workflow import verify_run


ASSAY_FIELDS = ["aligned_sequence", "value", "split", "assay_id", "units", "context", "round"]


def _read_csv(path, columns):
    data = Path(path).read_bytes()
    if len(data) > 20_000_000:
        raise ValueError("Assay input exceeds the 20 MB limit")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    if reader.fieldnames is None or len(reader.fieldnames) != len(columns) or set(reader.fieldnames) != set(columns):
        raise ValueError(f"Expected exactly these CSV columns: {columns}")
    rows = list(reader)
    if len(rows) > 10000 or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("Malformed assay row or excessive row count")
    return rows, hashlib.sha256(data).hexdigest()


def ingest_measurements(run_directory, measurements, output, existing=None):
    run_directory, output = Path(run_directory).resolve(), Path(output).resolve()
    if output.exists() or output.with_suffix(output.suffix + ".provenance.json").exists():
        raise ValueError("Assay outputs cannot overwrite an existing record")
    verify_run(run_directory)
    with (run_directory / "candidates.csv").open(encoding="utf-8", newline="") as stream:
        candidates = {r["candidate_id"]: r["aligned_sequence"] for r in csv.DictReader(stream)}
    incoming, measurements_hash = _read_csv(measurements, ["candidate_id", *ASSAY_FIELDS])
    rows, existing_hash = _read_csv(existing, ASSAY_FIELDS) if existing else ([], None)
    old_count, missing = len(rows), 0
    for row in incoming:
        if candidates.get(row["candidate_id"]) != row["aligned_sequence"]:
            raise ValueError("Measurement identity does not match an exported candidate sequence")
        if not row["value"].strip():
            missing += 1
            continue
        rows.append({key: row[key] for key in ASSAY_FIELDS})
    if len(rows) == old_count:
        raise ValueError("No new measured values were supplied")
    splits, contexts, exact_records = {}, set(), set()
    for row in rows:
        if row["split"] not in ("train", "validation", "test"):
            raise ValueError("Unknown assay split")
        sequence = row["aligned_sequence"]
        if sequence in splits and splits[sequence] != row["split"]:
            raise ValueError("An exact sequence cannot cross assay evaluation splits")
        splits[sequence] = row["split"]
        for key in ("assay_id", "units", "context"):
            if not row[key].strip() or row[key] != row[key].strip():
                raise ValueError(f"Missing or ambiguous {key}")
        contexts.add((row["assay_id"], row["units"], row["context"]))
        try:
            value, round_number = float(row["value"]), int(row["round"])
        except ValueError as error:
            raise ValueError("Measurement value and experimental round must be numeric") from error
        if not math.isfinite(value) or abs(value) > 1e100 or round_number < 0:
            raise ValueError("Nonfinite assay value or negative experimental round")
        record = tuple(row[key] for key in ASSAY_FIELDS)
        if record in exact_records:
            raise ValueError("Duplicate assay record; ingestion would count an identical record twice")
        exact_records.add(record)
    if len(contexts) != 1:
        raise ValueError("Assays with different readouts, units or contexts cannot be combined")
    if len(rows) > 10000:
        raise ValueError("Combined assay records exceed the row budget")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=ASSAY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    receipt = {"schema_version": 1, "source_run_manifest_sha256": digest(run_directory / "manifest.json"),
               "measurements_sha256": measurements_hash, "existing_sha256": existing_hash,
               "output_sha256": digest(output), "new_measurements": len(rows) - old_count,
               "unmeasured_rows_skipped": missing, "total_rows": len(rows),
               "value_policy": "Only supplied finite measurements; predictions never fill missing values"}
    save_json(output.with_suffix(output.suffix + ".provenance.json"), receipt)
    return receipt
