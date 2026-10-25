# Repository readiness audit — 10 October 2026

This is a code-and-artifact review of the current working directory. It does not claim that the existing experiments were rerun, that a production sampler exists, or that human approvals have occurred. No server deployment is needed for the proposed research workflow.

## What is actually available

| Area | Existing evidence | What remains |
| --- | --- | --- |
| Scope | [18-week project plan](../project_plan.md), benchmark contract, explicit role split and exclusions | Record actual human decisions; reconcile the original supplied document separately |
| Model acquisition | [Inspector](../../scripts/inspect_benchmark.py) verifies bytes, complete field/coupling keys, support, coordinates and selected scores | A reusable production loader that returns numerical arrays and model metadata |
| Biological coordinates | [Position map](../../data/benchmark/position_map.csv), provenance, separate anchor and supported witness | Biology lead review of the mapping, selected restrictions and scope |
| Direct sampling | [Week 2 pilot](../../examples/week02/adabmdca_pilot.py) uses the official toolkit on a tiny finite fixture | Support-aware, fully instrumented direct categorical sampling on the CM model |
| Encoded method | [Learning checks](../../examples/week03/learning_checks.py) demonstrate encoding, invalid words and MH correction on tiny prescribed examples | Actual proposal kernel, derivation for its schedule, full transition tests and a production implementation |
| Upstream reuse | [Gibbsiq adapter exercise](../../examples/week03/upstream_reuse_check.py) checks nine synthetic energies and all 16 encoded words | Supported finite-domain adapter, pinned install or attributable vendoring decision, actual integration |
| Evaluation | [Evaluation configuration](../../configs/evaluation.json) specifies proposed targets, caps and checkpoints | Reference generation, uncertainty calculations, timing harness, reporting and planner |
| Packaging | Python metadata and a hashed Windows/Python 3.13 pilot environment | Installable project package, production dependency lock, CLI, automated checks and clean-environment reproduction |

The saved [validation summary](../validation_summary.json) reports 37 project tests and 69 focused upstream tests on 8 October. Several project tests read retained pilot outputs; the actual source also tests input validation and small mathematical examples. These are useful preparation checks, not CM sampling, runtime or biological validation. This audit did not rerun the suite or the heavy acquisition inspector.

There is no `src/` project implementation directory. `pyproject.toml` has an empty dependency list. The pilot lock includes adabmDCA 1.0.0, thrml 0.1.4, JAX/JAXlib 0.11.2, NumPy 2.5.3 and Torch 2.14.1; these are recorded local pins, not a recommendation to update every package. ArviZ is not present in that lock. Confirm compatibility and add only the packages needed for the next deliverable.

The current `E:/projects/proteinwork` directory has no `.git` directory and `git status` reports that it is not a repository. A later implementation phase must establish the intended repository and remote before claiming versioned releases or a GitHub deployment workflow. This audit does not initialize or publish it.

## Benchmark facts that affect design

- The model has 96 aligned sites, 21 symbols, 2,016 fields and 2,010,960 canonical coupling records. The selected sequence score uses 96 fields and 4,560 pair terms, not every possible coupling value.
- The 362 negative-infinity fields define forbidden choices. Keep a separate support mask. Do not multiply a one-hot vector by negative infinity, replace infinity with a large number, or silently accept a requested forbidden residue.
- The Q5SLA5 coordinate anchor has two forbidden model states. It remains a coordinate/Hamming reference; it cannot initialize a chain. The supported witness is constructed and is not experimental evidence.
- Each exact CM case has only four states and two binary mutable sites. These are good real-model arithmetic tests, but every two-bit encoded combination is valid. Add synthetic sites with at least three categories to test invalid domain-wall words and their rejection.
- Calculated from the committed domain sizes, each broad case has 88 mutable sites. The development case requires 1,356 domain-wall bits and the evaluation case 1,366 bits after singleton sites are removed. Their categorical state spaces are about `10^105.29` and `10^105.60`. These are representation sizes, not measured costs or evidence of faster sampling.
- Both broad cases use the same model and mostly the same restrictions. One broad development case and one broad evaluation case can support a narrow case study; they cannot establish that a learned selector generalizes across protein families or constraint regimes.

Source: [benchmark configuration](../../configs/benchmark.json), [inspection](../../data/benchmark/inspection.json), [restriction cases](../../data/benchmark/restriction_cases.json) and [admission record](../week04/model_admission.json).

## Decisions to make before production benchmarking

1. **One temperature convention.** The model configuration records `beta`, while evaluation records `temperature`. Both equal one now, hiding reciprocal mistakes. Use one authoritative value and derive the other with `beta = 1 / temperature`; reject inconsistencies.
2. **Observation and budget units.** Define an attempted outer transition, internal encoded updates, one sweep, one recorded observation, warmup and checkpoint units. State whether the 900-second measurement and 3,600-second reference caps cover all four chains collectively. Vectorized execution does not make these units interchangeable.
3. **Trace size.** At 131,072 retained outer attempts per chain, four chains of 96-symbol `uint8` states require about 50 MB before events and energies. If a retained observation means an 88-site sweep and every site attempt must also be saved, states alone require about 4.4 GB. Choose chunked storage and define the units before running the broad benchmark. Keep internal proposal paths only where required for correction/audit; clearly distinguish those from post-transition chain states.
4. **Warmup and reference qualification.** Predeclare warmup, candidate checkpoints, how reference MCSE is computed and how inadequate references stop claims. Do not choose warmup by inspecting evaluation outcomes.
5. **Timing repetitions.** Four chains measure sampling disagreement; they are not four independent end-to-end timing replications. Predeclare a modest repeated-timing protocol and report failed/capped repeats as well as successful ones.
6. **Planner scope.** Start with a frozen transparent rule. Either retain the current four cases and describe a demonstration, or add a small, justified development/evaluation case grid before tuning. Avoid training a classifier from two cases.
7. **Method equivalence.** Define the difference between `fixed_corrected_encoded` and `corrected_encoded_planning_disabled`. If they execute the same schedule and settings, share one result and explicitly mark the equivalence.
8. **Human review remains human.** Pending lead/adviser review, independent score checks and the sampling-specialist review must have named reviewers and actual dated outcomes. Automated checks cannot fill those records. Code preparation can continue without falsely marking the decisions complete.

## Minimal proposed scaffold

The following paths are proposals, not implemented files. Retain the existing teaching exercises and admission inspector as independent evidence.

```text
src/proteinwork/
  model.py                 # immutable arrays, alphabet, support and canonical score
  load.py                  # hashes, parser, validation, cache provenance
  constraints.py           # biological request -> checked per-site domains
  conditional.py           # eliminate fixed sites, retain energy offset and maps
  kernels/
    direct.py              # categorical Gibbs baseline
    encoded.py             # finite-domain compiler adapter and defined proposal
    correction.py          # original-target acceptance and rejection events
  exact.py                 # enumerate bounded targets and transition matrices
  run.py                   # four chains, checkpoints, timing and resource caps
  records.py               # run manifest, chunked traces and resumability
  diagnostics.py           # observables, ESS/R-hat, MCSE and acceptance checks
  planner.py               # frozen rule and explicit direct fallback
  export.py                # trace-linked FASTA and transparent selection rules
  cli.py                   # validate, sample, reference, benchmark, summarize
tests/
  fixtures/                # q>=3, singleton, extreme and asymmetric examples
  test_model.py
  test_constraints.py
  test_kernels_exact.py
  test_records.py
  test_diagnostics.py
  test_workflow.py
docs/biology/              # short human-authored decisions and evidence
configs/                   # versioned model, cases, evaluation and frozen plan
results/<run_id>/          # ignored bulky traces; share manifests/summaries
```

Keep this as one Python package with a few clear records, not a general compiler framework. Suggested records:

| Record/interface | Required meaning |
| --- | --- |
| `Model` | Original model hash, alphabet order, `h`, canonical `J`, support, energy convention and precision |
| `RestrictedTarget` | Model identity, beta, ordered allowed domains, fixed values, original/reduced position maps, constant offset |
| `Plan` | Kernel identifier, exact proposal/schedule version, category order, relevant penalty, warmup, checkpoint policy, fixed planner version |
| `Kernel.step(state, rng)` | One documented outer transition; returns next state, next RNG state and an event, including self transitions |
| `TransitionEvent` | Attempt number, proposal validity, acceptance, reason, original energy and correction evidence when applicable |
| `RunManifest` | Target/plan/config hashes, seed and RNG, actual initial states, lock and source versions, hardware, timings, status, output hashes |

A direct Gibbs update draws a site's symbol from its exact conditional distribution over its allowed domain. Use numerically stable logits and sequential within-chain updates unless conditional independence is proved. Vectorizing separate chains is straightforward; updating interacting sites simultaneously from stale values generally changes the kernel.

When fixed sites are eliminated, absorb their interactions into each remaining site's field and preserve the fixed-fixed energy constant. Test reduced and full energies on every tiny state and selected real-model states. This both reduces work and gives Gibbsiq a finite-domain model without changing the target. Retain original arrays and hashes; the reduced representation is a recorded derived artifact.

Use a reference implementation first. Add JAX compilation and cached local fields only after agreement with that oracle. Test energy differences and cached updates against full recomputation to catch coupling orientation, update and precision errors.

## Delivery sequence and acceptance gates

| Plan stage | Deliverable | Acceptance evidence |
| --- | --- | --- |
| Before/within Week 5 | Repository/lock decision; reconciled target and run schemas | Small commands work from a clean location; pending decisions remain labelled |
| Week 5 | Production loader, restriction validator and fixed-site reduction | Reject altered hashes, incomplete keys, invalid support and numbering; reproduce independent original scores and reduced energies |
| Week 6 | Reference direct kernel plus complete run records | Tiny transition rows sum to one; nonnegative probabilities; target stationarity; connected supported states; exact rejection/retention semantics |
| Week 7 | Reference/diagnostic/export pipeline | Synthetic known-target diagnostics tested; trace length and hashes consistent; separate FASTA links to source chain/attempt; CM direct smoke run completes |
| Week 8 | Encoded proposal and correction specification | Reviewer can state the actual forward/reverse law, invalid handling, schedule, target, irreducibility and any path augmentation assumptions |
| Week 9 | Corrected encoded kernel | Exact matrices for q>=3, asymmetric proposals, fixed sites and invalid endpoints; uncorrected negative-control test fails; rejection self-mass retained |
| Week 10 | Both methods run the same CM target | Same hashes, restrictions, beta, precision and initialization policy; no claim based only on encoded energy equality |
| Week 11 | Profiled fixed-method comparison | Synchronized wall-clock costs include load, compile, tune, run, I/O and diagnostics; no speed claim from unmatched pilot timings |
| Week 12 | Simple frozen planner | Selection uses development evidence only, has a direct fallback and records its decision/cost; no within-run adaptation |
| Week 13 | Feature and evaluation freeze | Case IDs, config hashes, budgets, diagnostics and method equivalence recorded; remaining changes explicitly versioned |
| Weeks 14–15 | Main comparison and ablation | Qualified references or explicit failure; four chains per condition; repeated cost measurements; all capped/failed results retained |
| Weeks 16–18 | Reproduction, report, release and demonstration | Teammate reproduces a small workflow from instructions; raw evidence supports every claim; negative outcome remains a valid deliverable |

For a claimed reversible kernel, check detailed balance in addition to stationarity. A systematic Gibbs sweep may be stationary without being reversible. A finite encoded Gibbs trajectory cannot automatically use an energy-only acceptance rule: derive the correction for the exact kernel actually executed. The existing Gibbsiq audit explicitly does not supply a corrected protein sampler.

## Budget and likely value

The original plan allocates 182 hours per person; Weeks 1–4 account for 36 planned hours each. Weeks 5–18 therefore leave 146 planned hours for the systems lead and 146 for the biology lead, not 292 interchangeable programming hours. Actual hours are unknown. A loader, baseline, narrowly scoped corrected method, transparent selector and careful report are plausible only with substantial reuse and strict scope control. GPU-specific kernels, a web service, new model training and elaborate learned planning do not fit this budget.

The strongest near-term output is a trustworthy constrained-sampling and evaluation workflow. Whether encoding is faster is unresolved. The increase from 88 categorical sites to over 1,350 binary variables, dense couplings, correction cost and rejection behavior all need measurement. Prefer an optimized fair direct baseline; do not make encoding look useful by comparing it only with slow Python loops.

Use early stop gates: if the actual encoded correction cannot be justified, retain only the verified baseline and bounded mathematical result; if development measurements show no total-cost benefit, cap further optimization and report when the direct fallback wins. Biological impact remains limited to candidate generation and interpretation unless an independent mapped assay set is established. Existing reintegration labels were used during model fitting and cannot establish independent biological success.
