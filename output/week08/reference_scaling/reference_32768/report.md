# Protein candidate run

Status: reference_inadequate.

Case: `cm_dev_broad_8_12`. Method: `direct`. Recorded 131072 observations in four chains and 131072 distinct sequences.

Sampling assessment: `accuracy_not_reached`. Exported library: 0 candidates. Complete measured work before the final manifest took 189.901 seconds.

Raw traces retain repeated states. Reference runs skip candidate materialization and selection. Every exported candidate is linked to its original chain and draw. The library applies separate selection rules and does not represent equilibrium frequencies.

The exported sequence is the modeled CM segment corresponding to Q5SLA5 residues 2 through 90. Its model score is not measured activity or stability. No physical Extropic device was used.

Inspect `rules.json`, `diagnostics.json`, `selection.json` and `profile.json` for the applied rules, failed checks, objective and timing. Fill measured values and experimental context in `assay_worksheet.csv` before ingestion. Blank measurements are never replaced by predictions.
