# Week 8 reference sampling

This milestone improves the independent reference workflow and examines the remaining sampling limitations. It uses the existing chorismate mutase target and numerical criteria. Biological review is still pending.

## Useful software changes

Reference runs previously built and serialized a complete candidate pool despite exporting no candidates. They now skip candidate construction and optimization while retaining the full trace, original energy verification, residue support checks and diagnostic evidence. The empty candidate files retain their headers, and the manifest explicitly identifies that no candidate pool was materialized. Distinct sequence counts remain available.

The workflow also avoids copying the sole sampling segment when no checkpoint concatenation is necessary. Reference checkpoint stopping now uses reference qualification rather than the ordinary sampling assessment. Requests that combine a reference run with candidate export, assay prediction or another reference fail validation before execution.

## Controlled experiment

The [saved protocol](../../output/week08/reference_scaling/protocol.json) declares two sequential runs before execution. The first repeats the earlier 16,384 observations per chain using identical seeds and settings. The experiment compares every saved trace array and all diagnostic values against the previous result before proceeding.

The second run increases the observation count to 32,768 per chain and uses four new seeds. Warmup, sweep count, model, restrictions, temperature and all qualification thresholds remain fixed. Its 196,083,712 byte trace estimate fits the existing 512 MiB allocation guard. This estimate does not include all process memory. The existing CLI limit remains 3,600 seconds per run.

The experiment stops after the declared larger run regardless of its outcome. Unobserved indicators retain unknown uncertainty. No accuracy criterion is relaxed to obtain a passing result. Complete CLI wall time is recorded separately from the workflow timer, which excludes startup and final manifest creation.

## Measured storage improvement

The [replay result](../../output/week08/reference_scaling/replay_equivalence.json) confirms that every trace array and all diagnostic values match the original run exactly. Total saved artifacts fell from 94,465,955 bytes to 12,638,863 bytes. The change saves 81,827,092 bytes, or 86.6%, without removing observations or changing the sampling result.

The new replay completed in 155.059 seconds including CLI startup and exit. That is a separate measurement from the earlier workflow timer, so it is not presented as a runtime speedup. The primary demonstrated improvement is smaller output with identical sampling evidence.

## Larger run and next decision

The [larger reference](../../output/week08/reference_scaling/reference_32768/manifest.json) completed 131,072 observations in 202.334 seconds of complete CLI wall time. Its saved artifacts occupy 24.9 MB. No candidates were selected or exported.

Energy bulk effective sample size rose from 3,036.74 in the baseline to 5,833.08 in this independently seeded run, passing the required minimum of 4,000. Both energy and distance now pass the reference chain diagnostics. However, mean uncertainty for both observables remains above the permitted value.

Frequency observables exceeding the uncertainty limit decreased from 184 to 121. Unknown frequency uncertainties decreased from 13 to 11. The reference therefore remains `reference_inadequate`. The [summary](../../output/week08/reference_scaling/summary.json) preserves all results, including that failure.

All 11 remaining unknown values concern combinations of phenylalanine at accession position 11 with residues at position 15. None of those combinations appeared in any chain. Among observed residues, the largest uncertainty belongs to phenylalanine at position 65: 0.00683 against the 0.0025 limit. Positions 82 and 44 are the next largest contributors. The [position and chain count table](../../output/week08/reference_scaling/reference_32768_frequency_blockers.csv) identifies every remaining case.

The next sampling change should address exploration of these coupled residue choices and uncertainty in the common residues. Increasing the observation count alone improved several diagnostics but did not qualify the broad reference within this declared attempt. The existing candidate library remains exploratory. No missing frequency was converted into a known zero, no acceptance threshold was relaxed, and no biological restriction was invented.

## Reproduce

Use a fresh output directory from the project root:

```powershell
& ./.venv-pilots/Scripts/python.exe scripts/run_reference_scaling.py --output output/reference_scaling_replay
```

The runner preserves the baseline and every new run. Its output includes the protocol, raw traces, replay equivalence result, complete process times and frequency blocker tables with biological positions and counts for each chain. The baseline path can be changed with `--baseline` only when it contains the matching declared reference request.
