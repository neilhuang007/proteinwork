# Week 5 shared biology and computing exercise

The [report](../../reports/week05_report.md) follows the plan's first implementation stage. The shared exercise joins two tasks: residue counting for the biology lead's computing practice, and mapped substitution analysis for Neil's biology practice. It reuses the existing direct sampler and adds [shared_analysis.py](../../examples/week05/shared_analysis.py).

The implementation reads all 1,130 records in the pinned natural alignment. It counts residues at each column, excludes gaps from residue frequency denominators and preserves duplicate observations in the sampling analysis. Candidate annotations use the existing Q5SLA5 coordinate map and distinguish the fixed computational background from the two variable positions.

The fresh [direct run](../../output/week05/cm_direct/manifest.json) contains four chains with 512 observations each after 256 warmup sweeps. It passed all seven configured sampling checks. The [analysis summary](../../examples/week05/output/summary.json) records a total variation error of 0.0005279204513499813, three distinct observed sequences and 2,048 total observations. All three observed candidates contain 89 amino acids.

The natural frequencies are unweighted descriptions of the supplied alignment. The sampled frequencies belong to the exact conditional case with all other positions fixed. The script computes and records these denominators separately.

## Reproduce

Run these commands from the project root in the existing environment. Use fresh output directories for each new run.

```powershell
& ./.venv-pilots/Scripts/python.exe -m proteinwork run --request configs/cm_exact_direct.json --output output/week05/cm_direct_replay
& ./.venv-pilots/Scripts/python.exe examples/week05/shared_analysis.py --run output/week05/cm_direct_replay --output examples/week05/output_replay
& ./.venv-pilots/Scripts/python.exe -m unittest discover -s tests -p test_week05_shared_analysis.py -v
```

The analysis verifies the source manifest, original model scores, exact state support, expanded trace and candidate lengths before writing results. Its [provenance record](../../examples/week05/output/provenance.json) fingerprints the script, alignment, model, position map, source run and output files. The four focused tests exercise gap denominators, numbering across an internal gap, repeated and unobserved states, and malformed inputs.

## Saved results

* [Residue frequencies](../../examples/week05/output/residue_frequencies.csv)
* [Exact and observed pair frequencies](../../examples/week05/output/pair_frequencies.csv)
* [Candidate annotations](../../examples/week05/output/candidate_annotations.csv)
* [All alignment column profiles](../../examples/week05/output/alignment_profiles.json)
* [Validation record](validation.json)
