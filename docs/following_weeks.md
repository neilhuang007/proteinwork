# Work for the following weeks

Updated 9 November 2026. This record advances the existing project after the Week 6 work. It records completed software and experiments separately from the biological review and later evaluation that still need to happen. The original 18 week schedule remains a planning reference.

The subsequent [Week 8 milestone](week08/README.md) removes unused reference candidate records, demonstrates an 86.6% storage reduction with identical sampling evidence, and completes a larger reference attempt. Its position and chain count table now identifies the remaining broad sampling limitations.

## Completed development work

The core cleanup shares sampler integer validation, duplicate JSON key detection, target preparation, selection preflight and trace memory checks. Public interfaces and saved artifact fields remain stable. The two JSON readers retain their different contracts. Independent mathematical calculations used to verify the sampler remain separate. The changes also remove duplicate binary validation, reuse site indices during energy scoring and clarify resource counting names. No runtime improvement is claimed from this cleanup alone.

The comparison command now rejects duplicate run evidence. It groups repeated measurements by method and identical sampling settings, reports the median and range of qualified times, records unsuccessful runs separately and flags reused seeds within a method. It cannot establish that unrelated software was idle during a measurement; the experiment procedure controls the sampling and test jobs started for this project. Recorded times cover workflow preparation through export, before the final manifest. CLI startup and final manifest creation are outside this timer.

The complete Python suite passed 120 tests in 49.296 seconds. The seven added tests cover repetition accounting, duplicate evidence, source identity, mapped substitutions, residue restrictions and unknown uncertainty. See [the validation record](week07/validation.json).

## Broad sampling reference

The [independent direct reference](../output/week07/cm_reference/manifest.json) used the existing request with four distinct seeds, 2,048 warmup sweeps and 16,384 recorded observations per chain. It produced 65,536 observations with a recorded workflow time of 201.718 seconds and finished within the 3,600 second reference budget. This run used the implementation before the cleanup. Its original code fingerprint is retained in the manifest.

Its result is `reference_inadequate`. Energy bulk effective sample size was 3,036.74 against a minimum of 4,000. Energy and distance mean uncertainty also exceeded the reference thresholds. Of 1,482 frequency observables, 13 had unknown Monte Carlo standard error and 184 exceeded the permitted value of 0.0025. The unknown values remain unknown, including when a possible residue or pair was never observed.

The low energy Rhat of 1.00047 alone does not qualify the reference. No broad accuracy claim or successful time to accuracy is reported. This was one fixed endpoint attempt, not an exhaustive search or an exhausted compute budget. The [review summary](../output/week07/review_packet/summary.json) records every failed check, and the [frequency table](../output/week07/review_packet/reference_frequency_blockers.csv) identifies the observables that need attention.

The [profile](../output/week07/cm_reference/profile.json) assigns 127.13 seconds to the sampler, 17.83 seconds to diagnostics, 21.71 seconds to score verification and selection preparation, and 21.32 seconds to export. The reference also writes a 59.3 MB selection record and a 22.5 MB candidate pool even though it exports no candidates. A future reference workflow should review that reporting cost before scaling the observation budget, while retaining the raw trace, score checks and evidence needed for independent assessment.

## Candidate handoff

The new [review packet builder](../scripts/build_review_packet.py) verifies the saved source run before producing candidate annotations. It checks sequence hashes, positions, original model energies, residue support, occurrence counts and source observations. The same candidate identifier continues through every table.

The [packet](../output/week07/review_packet/summary.json) contains eight previously selected exploratory candidates, 440 substitution records and all 28 pairwise distances. The smallest pairwise distance is nine aligned positions, exceeding the requested minimum of five. Every exported sequence is the modeled 89 residue domain. Residue counts come from the pinned alignment of 1,130 records and use nongap denominators.

Daniel can use the [candidate annotations](../output/week07/review_packet/candidate_annotations.csv), [substitution table](../output/week07/review_packet/substitutions.csv) and [review worksheet](../output/week07/review_packet/review_worksheet.csv) to record decisions and their biological basis. Reviewer, decision, reason and proposed restrictions are blank. Building this packet does not establish that Daniel reviewed it. Alignment frequencies and model scores are not measurements of activity or stability.

## Repeated method measurements

The [experiment runner](../scripts/run_method_repeats.py) writes a protocol before starting six sequential CLI runs. Each method receives three repetitions of the existing exact development case. Seeds differ between repetitions and are paired across methods. The order alternates between direct first and THRML first. Each run starts in a fresh process and retains the existing thresholds, warmup, checkpoint policy and runtime cap. No parameter search is performed.

The initial six runs passed their sampling checks, but the comparison rejected inconsistent CPU metadata. One process reported the generic architecture while five reported the processor model. These [original runs and rejection](../output/week07/method_repeats/comparison_failure.json) are preserved. The runner now captures the observed processor identifier once and passes it unchanged to every child process. It does not rewrite any old manifest or relax the comparison guard.

See the replacement [declared protocol](../output/week07/method_repeats_verified/protocol.json) and [comparison](../output/week07/method_repeats_verified/comparison.json) for the execution order, settings, individual results and repetition summaries. These measurements concern the four state development target. They do not settle broad model performance or a final held out evaluation. A new code fingerprint also means earlier frozen plans cannot be reused without new compatible development evidence.

All six replacement runs passed their sampling checks at the first declared checkpoint of 512 observations per chain. Direct sampling had a median recorded time of 17.885 seconds, with a range of 17.293 to 18.251 seconds. Corrected encoded sampling through THRML had a median of 26.754 seconds, with a range of 23.233 to 27.869 seconds. These are descriptive results from three repetitions per method, without a claim of a precise confidence interval or a measured speedup from the code cleanup.

The [new frozen plan](../output/week07/frozen_plan.json) uses only these compatible development measurements and chooses direct sampling for the exact conditional case class. Its evidence covers one development target. The [evaluation request](../configs/cm_eval_planner_followup.json) supplies separate seeds and the existing evaluation target at model positions 16 and 24. This checks the refreshed planner workflow and does not replace the later final comparison.

The [evaluation run](../output/week07/planner_evaluation/manifest.json) passed all seven checks, recorded 2,048 observations and exported two candidates. Its recorded workflow time was 18.839 seconds. The method was frozen before this run; the evaluation result did not train or revise the plan.

## Next milestones and completion checks

1. Finish the broad reference investigation. Use the saved blocker table to identify rare indicators and slowly mixing observables. Declare the next sampling budget before execution and estimate peak memory and export size. Keep the existing thresholds. Completion means either a qualified independent reference or a documented limit on the broad comparison. Increasing the draw count does not guarantee qualification.
2. Complete the biological handoff. Daniel records which substitutions and positions need review, why a restriction is justified and its evidence. Neil validates each supplied rule against the position map and original support. Completion means a reviewed request and a new traceable run. No proposed rule is inferred from a blank worksheet.
3. Extend verified method comparisons. Choose additional small targets with known probabilities and explicit resource limits before considering broader encoding. Preserve the original target and every rejected observation. Completion means matched direct and corrected encoded results, including failures and logical resource costs. Physical device access remains a separate milestone.
4. Extend planner evidence beyond one small development target. Use only compatible development evidence and account for any development search cost. Freeze the choice before evaluation. Completion means matched coverage across the declared target classes and separately seeded evaluation runs, with the case size and generalization limits stated. The refreshed plan and its successful evaluation are a working starting point.

## Reproduce

Run from the repository root in the existing research environment. Each command needs a fresh output directory.

```powershell
& ./.venv-pilots/Scripts/python.exe -m unittest discover -s tests -v
& ./.venv-pilots/Scripts/python.exe -m proteinwork run --request configs/cm_reference.json --output output/reference_replay
& ./.venv-pilots/Scripts/python.exe scripts/build_review_packet.py --run output/research/cm_library_thrml --reference output/reference_replay --output output/review_replay
& ./.venv-pilots/Scripts/python.exe scripts/run_method_repeats.py --output output/method_repeats_replay
```

The reference replay uses the current implementation and records its own fingerprint. Exact timings will vary. Original saved runs are preserved. No biological assay records, adviser approval or physical hardware execution were added.
