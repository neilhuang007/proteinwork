# Week 1 scope draft

This project will develop and evaluate a tool for sampling candidate protein sequences from an existing trained model. The contribution is a measured rule for choosing a verified sampling method and update schedule. The benchmark family, model, restrictions, and numerical accuracy target remain decisions for the Week 4 benchmark contract.

This draft follows the supplied [revised project plan](C:/Users/neil_/Downloads/Protein_sequence_sampling_project_plan_revised.docx), especially its Project Scope, Timeline, Learning Resources, Implementation Plan, and Contingency Plans. The supporting Week 1 artifacts are the [candidate-family comparison](candidate_families.md) and [model-reading notes](model_reading.md). These documents record preparatory work and a synthetic teaching exercise; they do not establish an agreed benchmark or production sampling results.

## Research question

For one nonpathogenic protein family with publicly available trained model parameters and specified positional amino-acid restrictions, can a simple planner reduce the **total time to an agreed sampling accuracy**, compared with a tuned direct sampler and fixed corrected encoded alternatives, while preserving the same restricted target distribution?

Sampling means drawing sequences according to the model's probabilities. Model energy is a statistical score, not a physical energy measurement. The comparison must allow the planner to choose direct sampling when encoding and correction add more cost than they save. A reproducible result showing no speedup is a valid outcome; a speedup is a hypothesis to test.

## Work included

- Select one suitable nonpathogenic family and an existing public trained model. Use tiny artificial models for checks where all possible sequences can be enumerated.
- Define positional restrictions: fixed residues and allowed amino-acid sets at selected positions. Document the reference sequence and how model positions map to biological positions.
- Load and validate parameters, amino-acid order, indexing, score signs, pair counting, and temperature. Both leads independently check one reference score.
- Implement one direct sampler and one domain-wall encoded proposal method corrected against the original model. Reuse verified Gibbsiq code, adabmDCA's sequence-model workflow, thrml sampling tools, and jax numerical calculations where suitable.
- Build one simple planning rule from measured development costs. Select from a small documented set of verified methods and update schedules; log the decision and hold the selected plan fixed within a run.
- Retain complete raw sampling traces, including repeated states after rejected proposals, with scores, diagnostics, timings, settings, and random seeds.
- Export a separate FASTA candidate library linked to its source run and document any selection used to produce it.
- Evaluate correctness, sampling quality, runtime, diversity, and reproducibility. Compare experimental measurements only for variants whose sequences can be matched to measured records.

## Work excluded

New model training, unrestricted global constraints, a general-purpose compiler, custom accelerator code, physical-device integration, and laboratory experiments are outside the project. Hardware access must not determine whether it can finish. Candidate sequences without matching measurements must not be described as experimentally tested proteins.

## Preliminary interface

These are requirements to resolve into a benchmark contract, not an implemented API.

| Input | Required meaning or check |
| --- | --- |
| Trained model and provenance | Public parameters, version/source, layout, amino-acid order, score convention, pair counting, and temperature |
| Reference sequence and position map | Explicit mapping between model indexing and biological sequence positions |
| Positional restrictions | Fixed residues and allowed sets; reject invalid or inconsistent inputs |
| Run configuration | Method/update schedule, initialization, random seeds, software versions, numerical precision, hardware, and tuning budget |
| Benchmark cases and accuracy criteria | Separate development/evaluation cases and a documented accuracy target agreed by Week 4 |
| Experimental annotations, if available | Assay meaning and exact correspondence between measured variants and sequences |

| Output | Required record |
| --- | --- |
| Raw trace | Complete sequence states, including rejected-step repeats; linked scores and run settings |
| Candidate library | Separate FASTA export linked to its source run and documented selection procedure |
| Diagnostics | Single-site and selected pair frequencies, score distribution, distance from the reference, uncertainty, and sample dependence |
| Timing and planner record | Preparation, tuning, proposal, correction, transfers, reporting, and reason for the selected method/schedule |
| Final release | Runnable code, tests, fixed versions, model settings, example outputs, analysis, instructions, and reproducibility record |

## Acceptance principles for the eventual comparison

1. **One model interpretation.** Loading and restriction checks pass before optimization. The leads' reference-score and position interpretations agree. Disagreement triggers checks of alphabet order, indexing, signs, pair counting, and temperature.
2. **Verified target preservation.** Tiny-model results are checked against exact probabilities. Corrected-method tests check forward/reverse moves, invalid endpoints, normalized transition probabilities, and preservation of the target distribution. The direct method remains available throughout development.
3. **Fair evaluation.** Finish planner tuning before final evaluation. Use cases withheld from development, four independently initialized chains per condition, identical model/restrictions/precision/hardware, and equal tuning budgets. Compare the planner with a tuned direct sampler, a fixed encoded strategy, and correction with planning disabled.
4. **Accuracy before throughput.** The primary performance measure is total time to the agreed accuracy, including preparation and tuning. Attempted updates per second alone cannot establish improvement. Larger-model checks use separate long reference runs and agreement/uncertainty assessments; these are references, not exact probabilities.
5. **Traceable and reproducible outputs.** Preserve rejection repeats and complete run settings. Trace each exported sequence to a source run. A teammate reproduces the example in a clean environment before submission.
6. **Claims supported by evidence.** Report when the planner helps and when it does not. Biological interpretation distinguishes sampling diagnostics, sequence diversity, and measurements of actually tested variants.

Numerical thresholds, selected pair observables, chain lengths, reference-run adequacy, and library-selection rules are still open. This draft does not claim these acceptance checks have passed.

## Responsibilities

| Role | Planned responsibility |
| --- | --- |
| Systems lead | Implementation, mathematical checks, parameter interpretation, corrected sampling, planning, timing, diagnostics, and reproducibility |
| Biology lead | Family/model selection, reference sequence and position mapping, allowed residues, experimental provenance, sequence interpretation, and exported-sequence review |
| Both leads | Independent reference-score checks, benchmark agreement, development/evaluation separation, evidence review, and documented interface decisions |
| Adviser | Planned Week 4 benchmark review |
| Sampling specialist | Planned Week 8 review of encoded proposals, correction assumptions, reverse moves, invalid states, and small-model tests |
| Teammate reproducing the workflow | Planned independent clean-environment reproduction before final submission |

Names, assignments beyond these roles, feedback, and completed meetings have not been supplied. The plan budgets two 20-minute team check-ins per week; that is planned coordination, not a meeting record.

## Week 1 completion requirements

The plan's exact Week 1 work is: **review how researchers generate candidate libraries, compare two or three protein families, and begin the adabmDCA model-reading exercise.** The deliverable is **a scope draft and comparison of candidate protein families**, due at the end of Week 1. Planned effort is **8 hours per person**, including learning, documentation, and check-ins; no actual hours are claimed here.

- [x] Prepare this scope draft from the revised plan.
- [x] Document candidate-library generation practices and compare three families in [candidate_families.md](candidate_families.md), including availability and suitability of trained parameters, position mapping, and experimental evidence.
- [x] Begin the adabmDCA reading exercise and record model conventions, source references, and remaining questions in [model_reading.md](model_reading.md).
- [x] Record unresolved benchmark decisions and the next runnable examples without claiming a final family/model choice, verification of a real trained model, independent checks by both leads, or adviser feedback.

The supporting artifacts provide evidence for these Week 1 items. The synthetic exercise passed all nine arithmetic reference checks and 12 unit tests; [validation.json](validation.json) records the local run. This does not complete the later biological benchmark or production-loader milestones.

## Dependencies and decisions due Week 4

The family comparison must lead to usable public trained parameters and a reliable position map. The benchmark contract must record the reference sequence, permitted residues, parameter layout, signs, pair counting, temperature, expected outputs, accuracy criteria, and separate development/evaluation cases. Both leads must verify a reference score, and the adviser review must record decisions and changes when it occurs.

If parameters cannot be loaded or positions mapped by Week 4, switch to a smaller documented nonpathogenic example with usable trained parameters while continuing tiny-model checks. If too few measured variants match the benchmark, limit experimental analysis to matched records and retain the computational comparison and diversity analysis.

The next Week 2 work is to run a documented sequence-generation example and a tiny thrml model. The biology and systems pilot logs must preserve setup steps, software versions, outputs, failures, and runnable examples. Reading notes alone do not satisfy that milestone.

## Eighteen-week milestone alignment

Each deliverable is due at the end of its listed week. Hours are **planned hours per person**, not recorded effort. Weeks 1–3 are preparatory work; Week 4 is the plan-submission deadline. The total is 182 hours per person, or 364 team hours. Week 17 is reserved for recovery.

| Week | Deliverable and governing work | Hours/person |
| --- | --- | ---: |
| 1 | Scope draft and comparison of 2–3 candidate families; begin adabmDCA reading | 8 |
| 2 | Biology and systems pilot logs with runnable sequence-generation and tiny thrml examples | 10 |
| 3 | Draft benchmark requirements, Gibbsiq code audit, and related-work comparison; check one score by hand and numbering | 10 |
| 4 | Complete project plan and agreed benchmark requirements after adviser review | 8 |
| 5 | Passing loader, numbering, and score tests; restriction checks and jax exercises | 10 |
| 6 | Tested direct sampler and complete run records; exact tiny-model checks | 10 |
| 7 | Working reference pipeline from model input to sequence output; restrictions and ArviZ exercises | 10 |
| 8 | Reviewed encoded sampling rule, assumptions, and small-model test plan | 12 |
| 9 | Corrected sampler passes small-model checks and retains rejection repeats; direct method remains available | 12 |
| 10 | Initial prototype connecting both samplers to the protein benchmark with traceable exports | 12 |
| 11 | Full timing breakdown and results for fixed sampling strategies | 10 |
| 12 | Working planner selecting verified methods and schedules from development measurements | 10 |
| 13 | Revised prototype and fixed evaluation procedure; address peer feedback and stop adding features | 10 |
| 14 | Main runtime and quality results on held-out cases using four runs per condition and equal tuning budgets | 12 |
| 15 | Robustness checks, planning-disabled comparison, diversity analysis, and interpretation of measured variants | 12 |
| 16 | Completed testing and clean-environment reproducibility log | 10 |
| 17 | Essential recovery fixes/reruns, release candidate, final report draft, and presentation rehearsal | 8 |
| 18 | Final product, report, and demonstration with code, tests, settings, example outputs, and analysis packaged | 8 |

No calendar dates or later milestone completion are inferred from this schedule.
