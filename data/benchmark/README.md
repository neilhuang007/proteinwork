# Pinned CM acquisition

Run `python scripts/inspect_benchmark.py --download` from the repository root to
download the immutable upstream archive and reproduce the structural inspection.
The archive and downloaded source snapshots stay under ignored `raw/`; the
parameter member is inspected directly inside the ZIP. No model training or
parameter transformation is performed.

Small manifests, coordinate metadata, explicit support, restriction proposals,
and score evidence are retained locally here. `configs/benchmark.json` identifies the
checkpoint. The admission decision and limitations are recorded in
`docs/week04/model_admission.json` and `docs/week04/benchmark_evidence.md`.
