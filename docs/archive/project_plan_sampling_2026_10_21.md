# Protein sequence sampling project plan

This 18-week project compares ways to sample from one existing trained protein model. The research question is whether a corrected encoded method and a simple planner can reduce total time to a specified sampling accuracy while preserving the restricted target distribution. Choosing the direct method when it is faster is a valid planner outcome.

This is the working repository version of the supplied revised plan, with technical preparation through Week 4 recorded separately from human review. The original plan is identified by [its provenance record](../week01/plan_provenance.json). The [benchmark contract](../week04/benchmark_contract.md) and [evaluation configuration](../../configs/evaluation.json) contain the concrete proposed experiment. Adviser and both-lead approval are pending, not inferred from automated tests.

## Scope and method

Use one public pretrained model, one direct categorical sampler, one domain-wall proposal method corrected against the original model, and one measured planning rule. Specify fixed residues and allowed sets at individual positions. Start correctness checks on exactly enumerable models, then connect the verified methods to the protein benchmark. Retain a direct fallback throughout.

The proposed biological context is chorismate mutase, with an explicit coordinate anchor and a separate supported computational state. Negative-infinity fields define excluded states only after their source semantics are verified. Model bytes, support, alphabet, restrictions and temperature must be identical across methods. A different encoding or proposal must not silently change that target.

Reuse verified Gibbsiq model/encoding components, adabmDCA workflow conventions, thrml sampling primitives and JAX calculations. The Week 3 audit defines the reuse boundary: encoding correctness alone does not establish the correctness of a corrected proposal or its mixing rate. Derive forward/reverse or pathwise correction for the actual update schedule before claiming target preservation.

New training, unrestricted global constraints, a general compiler, custom accelerator kernels, physical devices and laboratory experiments remain excluded. Unmeasured sequences are computational candidates, not validated functional proteins. Assays used to fit the selected model cannot serve as independent biological evaluation.

## Inputs and outputs

Inputs comprise immutable model files and hashes, alphabet order, canonical pair convention, explicit temperature, biological reference and position map, allowed residues, independent seeds and initial states, and separate development/evaluation cases. A full production loader is still a Week 5 deliverable.

Each future run must save every state after every attempted update, including rejection repeats, with original-model energy, method and schedule, random settings, acceptance events, diagnostics, software/hardware, and a timing breakdown. FASTA is a separate export linked to its source run and selection rule. Week 2 toolkit outputs preserve all states their APIs returned; they are not represented as complete per-attempt production traces.

## Evaluation and acceptance

First pass input/restriction checks, independent score calculations and exact small-model transition tests. Compare four independently initialized chains per condition with the same model, restrictions, precision, hardware and tuning budget. Reserve evaluation cases before tuning; freeze the planner before the main comparison and keep it fixed within each run.

The primary outcome is end-to-end time to the agreed accuracy, including preparation, tuning, correction and reporting. Track site and selected-pair frequencies, energy, Hamming distance, uncertainty and sample dependence. Use exact targets for small conditional cases and separately qualified long reference runs for larger cases. The concrete proposed tolerances, budgets and diagnostic rules live in [evaluation.json](../../configs/evaluation.json); they are design choices for review, not achieved results.

## Responsibilities and reviews

The systems lead owns code, mathematics, timing and reproducibility. The biology lead owns source provenance, reference mapping, restrictions and biological interpretation. Both leads independently check a real-model score and agree on the benchmark. Names, actual work hours and completed meetings are not inferred.

The plan retains two 20-minute check-ins each week, a Week 4 adviser benchmark review, and a Week 8 sampling-specialist review. The [review packet](../week04/review_packet.md) makes the Week 4 decisions reviewable. These meetings must be logged when they actually occur.

## Timeline

Hours are planned per person and include learning, documentation and check-ins. Total: 182 hours per person, 364 team hours. Calendar deadlines have not been supplied.

| Week | Work and deliverable | Hours |
| --- | --- | ---: |
| 1 | Scope draft, three-family comparison, initial model reading | 8 |
| 2 | Runnable documented sequence-generation and tiny thrml pilots; setup/output/failure logs | 10 |
| 3 | Score and numbering checks, Gibbsiq audit, related-work comparison, draft benchmark requirements | 10 |
| 4 | Updated complete plan and concrete benchmark contract; obtain both-lead and adviser review | 8 |
| 5 | Production model loader and positional restriction validation; JAX array/compilation exercises | 10 |
| 6 | Direct sampler, exact small-model tests, reproducible full run records | 10 |
| 7 | Connected reference pipeline, diagnostics and traceable FASTA; ArviZ exercises | 10 |
| 8 | Reviewed encoded proposal, correction assumptions and small-model test plan | 12 |
| 9 | Corrected sampler passes transition/target checks and retains rejections | 12 |
| 10 | Both verified methods connected to the protein benchmark | 12 |
| 11 | Timing workshop, full cost breakdown and fixed-strategy results | 10 |
| 12 | One development-driven rule selecting method and update schedule | 10 |
| 13 | Peer feedback, unusual inputs, fixed evaluation procedure and feature freeze | 10 |
| 14 | Main held-out runtime/quality comparison, four chains and equal budgets | 12 |
| 15 | Failed/uncertain reruns, planning-disabled comparison, diversity and matched assay interpretation | 12 |
| 16 | Resolve failures and document a teammate's clean-environment reproduction | 10 |
| 17 | Recovery, essential fixes/reruns, release candidate, final report and presentation rehearsal | 8 |
| 18 | Package code, data/settings, tests, raw results, instructions, final report and demonstration | 8 |

Learning stays tied to deliverables: adabmDCA and assay provenance in Weeks 1–3; thrml in Weeks 2–3; encoding/correction in Weeks 3 and 8; JAX in Weeks 5 and 11; diagnostics in Weeks 6–7. No additional training hours are added to the budget.

## Contingencies and change control

- If model support or positions cannot be verified, use another documented nonpathogenic pretrained example and assign a new benchmark identifier. Continue tiny-model checks; do not train a replacement or silently alter parameters.
- If score or numbering checks disagree, pause optimization and resolve alphabet, indexing, signs, pair counting and temperature first.
- If correction fails, isolate an enumerable case, seek the Week 8 specialist review and retain the direct baseline. Protect correctness before performance.
- If encoding does not improve total cost, cap further tuning at one week and let the planner choose direct sampling. A supported negative result is valid.
- If chains disagree, inspect initialization, autocorrelation and reference adequacy. Record capped failures; rerun every method under any revised benchmark.
- If setup blocks two sessions, restore recorded versions and reduce workload on CPU. If measured assay coverage is insufficient, limit claims to computational quality and diversity.
- If a milestone slips more than one week, remove optional sweeps/packaging and use Week 17 recovery. Protect the final verified workflow, report and demonstration.

Record each interface or benchmark change with reason, owner role, reviewer, deadline and acceptance check. Proposed Week 4 settings may change after review or development evidence; assign revisions and rerun affected checks. Do not change thresholds in response to held-out evaluation results.

