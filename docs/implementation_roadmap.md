# Implementation roadmap for protein sequence sampling

Prepared 10 October 2026. **Recommendation: proceed with a bounded computational study, with a working direct sampler as the first deliverable.** Build the encoded method only far enough to test whether it reduces total time to the same sampling accuracy. A useful negative result is possible; neither a speedup nor improved protein function has been demonstrated.

This document turns the [18-week plan](project_plan.md) and the original planning document into implementable components, acceptance checks and decisions. The accompanying [biology lead guide](biology_lead_guide.md) is the practical nonprogramming handoff. The three research notes cover [current sampling methods](research/sampling_methods_2026_10.md), [biological relevance](research/biology_and_impact_2026_10.md), and [repository readiness](research/repository_readiness_2026_10.md).

## What we are building

A researcher supplies one saved protein model and permitted amino acids at each position. The software validates these inputs, samples either directly or through a corrected binary representation, and returns complete sequence histories, accuracy diagnostics, timing and a separately selected FASTA library. A small planner chooses between verified methods using development measurements and input characteristics.

Success means obtaining trustworthy samples sooner for a defined task. The model, restrictions and target probabilities remain identical across compared methods. This is a command-line research package with reports. Model training, a general compiler, a web application, physical accelerator integration and laboratory experiments are outside the first release.

    Model files + biological map + restriction case
                         |
               Validate and condition
                         |
                Frozen method planner
                   /             \
       Direct categorical     Encoded proposal
              Gibbs           + exact correction
                   \             /
                Complete outer-chain trace
                         |
          Accuracy checks + full runtime report
                         |
           Separately selected FASTA candidates

### How the two supplied plans fit together

The original Word document proposes roughly three weeks of feasibility and proposal preparation, followed by milestones A–D: reference workflow, corrected proposal, useful execution regime/planner, and manuscript/release. The repository plan expands this into 18 weeks. Keep that 18-week schedule; do not restart the completed preparation or add an unbudgeted second project.

The original document read for this task was **planning document.docx**, at the user-supplied OneDrive path. Its SHA-256 is **036176fb764763e93fe90ea65f3fcbc3472946d6ece1cc7e706c29c51fac6dec**. It is distinct from the revised DOCX identified in [the existing provenance record](week01/plan_provenance.json). This roadmap proposes future work; it does not change the existing benchmark, settings or review status.

## Current starting point

| Available evidence | What still needs to be built or decided |
| --- | --- |
| Model acquisition inspector, immutable source hashes, full-score arithmetic and position map | Production loader and stable internal model format |
| Small adabmDCA and THRML pilots; audited Gibbsiq encoding primitives | Project direct sampler, actual corrected encoded sampler and backend integration |
| Two four-state exact CM cases and two broad restriction cases | Production transition tests, qualified broad references and runtime measurements |
| Proposed accuracy settings and review packet | Human agreement on the use case, restrictions and evaluation; specialist review of the actual proposal |
| Recorded 37 project tests and 69 focused upstream tests | Tests of the new production path; existing counts are historical preparation evidence |

The inspected workspace has no Git metadata. Establish a repository and agreed remote during implementation before claiming a commit-based reproduction. Preserve the current examples and independent inspector as reference evidence. See the [readiness audit](research/repository_readiness_2026_10.md) for the source-by-source inventory.

### Benchmark facts that implementation must preserve

- The selected CM checkpoint has **96 aligned positions**, alphabet **-ACDEFGHIKLMNPQRSTVWY**, and original negative-infinity fields defining forbidden site choices. Do not replace them with large finite penalties in the target.
- The coordinate anchor maps to **Q5SLA5 residues 2–90**, an 89-residue segment of a 355-residue protein. It is outside model support at two positions, so its energy is positive infinity. It is usable for numbering and Hamming distance, not initialization or a finite mutation-energy reference.
- The supported witness is a constructed computational state. The separate finite arithmetic fixture is also not a valid initialization for the restricted cases. Broad chains require independently randomized supported starts.
- The broad cases have **88 mutable positions**, expanding to **1,356 or 1,366 domain-wall bits**. Their dense interactions may remove much of the hoped-for parallelism.
- The existing exact cases have only two binary choices at each of two sites. Their encodings have no invalid words, so they cannot test the main encoding failure mode.
- Labels used in feedback training are ineligible as independent model evaluation. No independent mapped assay holdout currently exists. Current restriction sites are computational examples without an established functional rationale.

These facts come from the [benchmark contract](week04/benchmark_contract.md), [admission evidence](week04/benchmark_evidence.md) and [inspection record](../data/benchmark/inspection.json).

## Would this actually be impactful

**Potentially useful for sampling methodology; biological and practical impact remain unproven.** Two correct samplers targeting the same distribution have the same expected biological properties at equilibrium. A faster sampler can provide more effective exploration per hour; it cannot make the fixed model intrinsically better at predicting function. That distinction should anchor the proposal and final report.

The most credible near-term contribution is a reproducible support-aware workflow, a carefully corrected proposal, and measured evidence about when encoding helps or hurts. Encoding and Metropolis–Hastings correction are established methods. A planner is a contribution only if its choices improve outcomes relative to strong fixed strategies; choosing direct everywhere is a valid result but not evidence of useful adaptive selection.

Before investing in extensive optimization, the biology lead should document two intended-user or adviser conversations, if available: current model, requested library size, restriction patterns, job frequency, waiting time and what saved time would enable. Record unavailable evidence as unknown. The systems lead then profiles that workflow. If sampling occupies fraction f of total time and is accelerated by factor s, overall speedup is at most 1 / ((1-f) + f/s), before new overhead. For example, accelerating a 10% sampling component tenfold improves the whole workflow by only about 1.10 times. This is an illustrative calculation, not a measured project result.

Use three levels of claim:

| Claim | Evidence needed |
| --- | --- |
| Correct computational workflow | Independent energy checks, exact kernels, qualified references, complete reproducible records |
| Useful performance improvement | Same accuracy, fair optimized baseline, repeated end-to-end timings, absolute time saved on a documented task |
| Biological usefulness of candidates | Independent measured evaluation with appropriate sequence/construct coverage; outside the current release. An improvement claim additionally needs a comparator and defined outcome, such as functional yield within equal elapsed time |

For a practical go decision, propose **at least 20% less end-to-end time** on a biologically useful held-out workload, with repeated timings supporting the difference and an absolute saving the intended user considers worthwhile. This margin is a suggested project decision, not a literature standard or current configuration. Agree it before evaluation. The research notes discuss alternative ambitions; this roadmap recommends one consistent initial gate.

## Technology decisions from current research

The research was checked on 10 October 2026 and is targeted, not exhaustive.

| Technology | Role in this project |
| --- | --- |
| adabmDCA | Reuse conventions and benchmark its direct-sampling approach. Python release 1.0.0 is current in the checked source; the paper's 2.0 naming is different. Its finite-only loading path cannot directly accept this checkpoint. [Official package](https://pypi.org/project/adabmDCA/) |
| Gibbsiq domain-wall primitives | Reuse the locally audited codec/energy construction behind a narrow adapter, pinned to the audited commit. The audit does not validate its full runtime as our corrected sampler. [Local audit](week03/research_audit.md) |
| THRML and JAX | Use JAX for the direct and initial encoded kernels. Add a THRML backend only when its actual update schedule satisfies the correction contract. Start on CPU; native Windows NVIDIA GPU execution is not a supported JAX path in the checked installation matrix. [THRML](https://github.com/extropic-ai/thrml), [JAX installation](https://docs.jax.dev/en/latest/installation.html) |
| ArviZ | Compute diagnostics for declared numeric observables, including energy and Hamming distance. Categorical residue indices themselves are not meaningful continuous diagnostics. [Diagnostics guide](https://arviz-devs.github.io/EABM/Chapters/MCMC_diagnostics.html) |
| Thermalizing Stochastic Programs v2 | Relevant August 2026 prior work on approximate compilation. Reported demonstrations are software; the linked Thermalizers repository could not be accessed in this check. It is neither an exactness guarantee nor a required dependency. [Paper v2](https://arxiv.org/html/2608.01615v2) |
| New discrete samplers | The June 2026 PDHAMS paper and October 2026 categorical diffusion preprint inform related work. Their adaptation/training requirements exceed the minimal study. Add at most one stronger direct baseline if mixing evidence demands it. [Research comparison](research/sampling_methods_2026_10.md) |
| RFdiffusion3, BoltzGen and ProteinGym | Relevant design and assay context, but they solve different problems or provide different data. They are not same-target runtime baselines. Avoid a scope expansion into a second design model. [Biology research](research/biology_and_impact_2026_10.md) |

## Specify the complete mechanism before optimization

### 1 Validate and condition the original target

For sequence x and permitted domains D_i, use the existing convention:

    E(x) = -sum_i h_i(x_i) - sum_(i<j) J_ij(x_i, x_j)
    pi(x) proportional to exp(-beta * E(x)), x in product_i D_i

Validate file hashes, alphabet order, unique complete keys, pair orientation, coordinates, finite couplings and permissible fields. Reject explicit forbidden requested residues, conflicting restrictions and empty domains. Gather selected fields rather than multiplying one-hot zeros by infinity. Record beta explicitly; derive temperature as 1/beta.

Remove fixed sites by adding their interactions into mutable-site fields and keeping the fixed background energy constant. Pack the allowed categories, retaining a reversible map to original residue symbols. This is an exact conditional reduction, not new training or an approximate target. Compare every reduced score with the independent full-score oracle on tiny models and selected real states.

### 2 Implement a competitive direct sampler

For a selected mutable position, evaluate all allowed residues using its original field and incident pair terms, then sample the normalized conditional with a stable softmax. Begin with state-independent random site selection. A systematic sweep is also valid for direct Gibbs and can be benchmarked after verification.

Keep local energy differences and optionally cached conditional scores. Periodically compare against the full-score oracle. Retain repeats when Gibbs selects the existing residue. Batch independent chains, rather than independently updating coupled sites as if they were conditionally independent. A fully fixed case returns its unique supported state without attempting to select from an empty update set. For optional systematic sweeps, persist the schedule phase; test the homogeneous kernel at sweep boundaries or on the augmented sequence-and-phase state.

### 3 Add one explicit corrected encoded proposal

Use a fixed category order and a unique prefix-one domain-wall code c(x). A site with q allowed choices uses q-1 bits; fixed sites use none. Construct a finite binary surrogate energy E_s on all binary words with a fixed finite penalty. Initially use the audited exact valid-state construction.

Define R as one random-scan binary Gibbs update with fixed state-independent, strictly positive selection probability for every free bit, under exp(-beta_s E_s). Each constituent update is reversible, so R and its fixed-length power Q = R^m are reversible. For cases with free bits, m must be a positive integer. Fully fixed cases take the deterministic path above. Fix m, penalty, order and beta_s before the run.

For each outer protein-chain attempt:

1. Encode the current valid sequence x and perform exactly m inner updates.
2. Reject an invalid endpoint and retain x. Do not repair it or redraw until valid.
3. Decode a valid endpoint y and compute original and surrogate energy differences.
4. Accept with the following log probability; otherwise retain x.
5. Save the resulting outer state, original energy and event, even after rejection.

    log_alpha = min(0, -beta * (E(y) - E(x))
                      + beta_s * (E_s(c(y)) - E_s(c(x))))

This formula follows because reversibility supplies the reverse/forward proposal ratio through the surrogate energies; it does not require the inner chain to equilibrate. It applies to this specified kernel, not arbitrary THRML updates. E_s and beta_s must describe the actual effective backend coefficients and temperature, including any numerical conversion. Backend execution starts at the supplied code, with no hidden reinitialization, warmup or extra updates. Numerical conditional-probability errors still need testing; merely using float64 in the outer acceptance rule does not fix an incorrect inner kernel. The full derivation and sources are in [the sampling note](research/sampling_methods_2026_10.md).

If valid-state surrogate energy equals original energy plus a constant and beta_s equals beta, **every valid endpoint is accepted**. Rejection then only excludes invalid words. A different proposal temperature creates a nontrivial MH correction and must be tested. A high valid acceptance rate alone says nothing about decorrelation or speed.

A forward-only sweep generally cannot use this reversible-kernel shortcut. If profiling justifies an alternative, verify a fair forward/reverse sweep mixture with exact constituent block conditionals. Arbitrary adaptive or pathwise schedules require a separate derivation and are optional future work. Check connected support, positive self-transition probability and numerical underflow on small examples.

### 4 Select a method without changing the target

The first planner is a small documented rule, not a learned model. Inputs are mutable-site count, allowed-domain sizes, encoded size and development measurements. Output is a verified method and fixed schedule, with the reason and cost estimate. Select once before initialization and keep the choice fixed throughout that run. Default to direct when evidence is weak or the input is outside the calibrated regime.

The present one-broad-case-per-split design only supports a proof of concept. To study meaningful selection, propose a small revision containing two additional development and two additional evaluation restriction cases spanning few/many mutable sites and narrow/wide domains. The biology lead supplies rationale; both leads approve the split before any added timing is used. Retain one model and the time budget. If this extension is unaffordable, limit the planner claim explicitly.

## Proposed software scaffold

The following paths are **to build**, not existing commands or implemented APIs. Keep one package and one CLI; no service layer is needed.

    src/proteinwork/
      cli.py                     validate, sample, evaluate, export
      contracts.py               validated model, case, plan and run records
      model_io.py                source-specific parsing and hashes
      target.py                  support, original energy, exact conditioning
      samplers/direct.py         categorical Gibbs and local energy updates
      samplers/encoded.py        surrogate kernel and outer MH correction
      encoding.py                audited Gibbsiq codec adapter
      backends/jax_runtime.py    arrays, keys, compiled transition chunks
      backends/thrml_runtime.py  optional verified schedule adapter
      planner.py                 frozen rule and direct fallback
      traces.py                  chunked states, events and replay metadata
      diagnostics.py             observables, MCSE, ESS and R-hat
      evaluation.py              references, caps, checkpoints and comparisons
      export.py                  linked FASTA and selection ledger
    tests/
      test_target_contract.py    malformed inputs and exact reduction
      test_direct_kernel.py     independent exact transition tests
      test_encoded_kernel.py    invalid words, correction and schedule tests
      test_run_records.py       attempts, replay, caps and export provenance
    configs/                     existing benchmark plus versioned run plans
    data/benchmark/              existing frozen inputs and small fixtures
    runs/<run_id>/               immutable manifests and trace chunks
    docs/biology/                human-authored brief and review evidence

Use standard-library records/dataclasses and existing dependencies where sufficient. Pin the package, JAX, THRML, ArviZ and reused Gibbsiq version in a platform-specific lock. Do not replace the independent inspector with the same implementation it is meant to check.

| Interface | Must contain or guarantee |
| --- | --- |
| ValidatedModel | Original hashes, h/J convention, alphabet, support, beta, position map; no silent parameter changes |
| RestrictedCase | Case ID/revision, split, explicit domains, original/reduced mapping, background constant, feasible state construction |
| SamplerPlan | Method, actual schedule, encoding order, m, penalty, beta_s, effective backend coefficients, precision, backend, correctness-test version, planner version |
| AttemptEvent | Chain/attempt index, accepted or invalid/rejected event, original energy, inner-update count and replay information |
| RunManifest | Model/case/plan hashes, actual initial states, independent RNG streams, versions/hardware, warmup policy, time and memory costs, completion status |
| EvaluationResult | Observable errors/uncertainty, reference ID, first passing checkpoint or capped failure, per-condition wall time |
| ExportLedger | Source run and chain/attempt IDs, aligned sequence, exported sequence, selection rule and candidate-only label |

Store the initial state and every outer post-attempt state in bounded binary array chunks, with JSON metadata and small CSV summaries. Raw warmup remains available and is explicitly excluded from analysis. Inner surrogate states are not protein-target samples; retain inner settings/keys and auditable event data sufficient for replay. A resume checkpoint must retain current state, RNG state/counter, schedule phase if applicable, plan hash and the last committed chunk so an interrupted run cannot silently skip or duplicate attempts.

## Build sequence and acceptance gates

The existing plan budgets 182 hours per person. Weeks 1–4 account for 36 planned hours each, leaving **146 planned hours each** across Weeks 5–18. These include learning, review and writing; they are not all programming hours. Keep the existing weekly allocations unchanged: Weeks 5–18 respectively have 10, 10, 10, 12, 12, 12, 10, 10, 10, 12, 12, 10, 8 and 8 hours per person. Absorb preliminary handoff work into pending Week 4 review and Week 5 time. Actual hours and calendar dates remain unrecorded.

| When and owner | Build or supply | Done when |
| --- | --- | --- |
| Before Week 5, both leads | One-page use case, reviewed benchmark decisions, unresolved questions; agree Git/remote and artifact policy | Decisions have real names/dates, unknowns stay explicit; proposed loader work can proceed without inventing approval |
| Week 5, systems with biology review | Production parser, support, exact conditioning, validated records | Unsupported/conflicting inputs fail; reduced/full energies agree, including fixed interactions; biology verifies requested coordinates |
| Week 6, systems | Direct Gibbs, independent toy matrix builder, complete trace writer | Exact stationarity and contract tolerances pass; original-state repeats and deterministic replay are correct |
| Week 7, systems and biology | End-to-end direct run, reference runner, observables, candidate export | A run produces trace, status, timing, report and traceable FASTA; inadequate references are explicitly flagged; biology checks displayed coordinates and labels |
| Week 8, systems and sampling specialist | Written encoded kernel, correction derivation and tiny tests; biology rechecks domains | Actual schedule and reverse law are reviewed; implementation assumptions are explicit; do not infer review from test results |
| Week 9, systems | Encoded outer sampler with forced rejection cases | q=3/q=4 matrix tests pass for m=1 and m>1, equal/unequal temperatures and invalid endpoints; all rejection repeats survive |
| Week 10, both leads | Both paths on the fixed CM benchmark | Real-model score agreement, exact-case accuracy and identical support/settings; no mutation of original target |
| Week 11, systems | Fair CPU profile, incremental updates, JAX chunks and synchronized timing | Fixed-strategy cost breakdown and accuracy exist; identify whether encoding has a plausible useful regime |
| Week 12, systems with biology | Frozen simple selection rule; optional approved case extension | Aggregate development budget enforced; choices justified by development evidence; direct fallback works |
| Week 13, both leads | Resolve edge cases, review outputs, freeze protocol and package | Config revisions, method IDs, reference rules, timing repetitions and claim scope are signed off before evaluation |
| Week 14, systems | Main held-out comparison | Four independent chains per condition, equal budgets, complete logs and capped failures; report absolute seconds and ratios |
| Week 15, both leads | Failure investigation, planning-disabled comparison, diversity/assay interpretation | Equivalent methods are labelled equivalent; any assay matches have exact provenance and denominators; no held-out tuning |
| Week 16, teammate with systems support | Clean-environment reproduction | A second person reproduces a tiny complete run and one bounded benchmark from documented commands |
| Week 17, both leads | Recovery, essential fixes, report/rehearsal | Scientific failures resolved or explicitly bounded; optional work cut to protect the final result |
| Week 18, both leads | Versioned release, raw results, report and demonstration | Every claim links to data/source, a fresh user can run the reference workflow, and limitations are visible |

A release can finish with a verified direct method and a well-supported negative encoded result. It cannot claim a completed comparison if correctness or reference adequacy remains unresolved.

## Verification and experiment protocol

### Correctness before speed

Use independent literal energies and exact transition matrices, not tests that only restate the production implementation. Cover asymmetric fields/couplings, restricted and singleton sites, q=3/q=4 encodings, invalid endpoints, valid MH rejections, multiple inner lengths and finite penalties. Include a deliberately wrong forward-sweep shortcut that fails, alongside the reversible implementation that passes. Check row sums, nonnegativity, stationarity and detailed balance when claimed. The existing four-state CM references then connect this reasoning to real parameters.

For empirical comparisons, preserve the proposed [evaluation settings](../configs/evaluation.json): four independently initialized chains; exact-case total variation and site/selected-pair error at most 0.02; energy/Hamming mean error at most 0.05 reference standard deviations; R-hat below 1.01 and bulk/tail ESS at least 400 for the declared observables. These remain proposed acceptance criteria, not achieved results.

Qualify separate broad-case references with their stricter ESS/MCSE requirements before judging broad accuracy. Require absolute candidate/reference error plus 1.96 times their combined MCSE to fit the configured tolerance, including the analogous scaled rule for means. This is the contract's descriptive acceptance rule, not a simultaneous confidence guarantee. Waive variance diagnostics only for an observable provably constant under the target, requiring its exact value; a constant observed trace is insufficient. A long run is not automatically ground truth. If the reference fails its cap, report **reference_inadequate** and an inconclusive broad comparison. Diagnostics do not establish all unmeasured high-dimensional properties.

### Resolve these contract details before freezing

1. **Temperature:** make beta authoritative and validate any redundant T field with beta=1/T. Include a non-unit test; the present beta=T=1 conceals reciprocal errors.
2. **Observation unit:** propose one recorded observation per outer attempt, with the initial state separate. For direct Gibbs, one attempt is a site update; for encoding, it is m inner binary updates plus one outer decision. Record all inner work in cost counters. Compare wall time, not equal attempt counts.
3. **Budget unit:** define the 900-second measurement cap as wall time for the whole four-chain condition, not an allowance independently spent four times. Specify whether chains execute concurrently; keep that policy identical across methods. Treat the 3,600-second reference cap analogously.
4. **Warmup and checkpoints:** choose warmup policy on development data and freeze it. Evaluate configured post-warmup checkpoints without thinning the raw trace. Include checkpoint diagnostics and writing in elapsed cost.
5. **Tuning versus held-out evaluation:** use the 120-second allowance as an aggregate development search budget per method/case, including search execution and its setup. For the planner this is 120 seconds across all candidates, never 120 seconds for each candidate plus another planner allowance. Record search ID, method, development case, candidate, seconds and reused measurement IDs in a cost ledger. Shared evidence consumes the receiving method's allowance; do not give the planner uncharged access to larger baseline searches. Do not retune on held-out timing. Proposed accounting is total-to-result time = cold job time + the method's full development-search cost over its declared development cases. Also report cold job time separately. The same offline charge appears in each repeated single-use comparison; do not rerun or describe it as newly executed search. Record shared reference-generation cost separately as evaluation overhead. Agree this accounting before the freeze.
6. **Timing variability:** propose three independent repeated conditions, each with four chains and new seeds, under the same caps and machine policy. Freeze a versioned repeat-ID/seed table beyond the current single evaluation seed set. Publish all times and failures; report the median and range rather than treating three repetitions as a precise confidence interval. If compute is insufficient, reduce the case count before freezing.
7. **Failure semantics:** record accuracy_not_reached, reference_inadequate, numerical_error or invalid_input separately. A capped run has a censored time, not a successful value equal to the cap. Do not calculate a normal speedup ratio from it.

These are proposed clarifications and revisions. They have not been silently written into the existing evaluation JSON. Match the final contract, implementation and report before measurement.

### Performance priorities and stop rules

First compact supported domains and eliminate fixed positions exactly. Next use incremental energies, compile bounded JAX loops, batch independent chains, retain arrays on device and write traces in chunks. Measure cold compilation, transfer, proposal, correction, diagnostics and storage separately; synchronize asynchronous operations. JAX's official guide explains why apparent kernel speedups can disappear at application scale. [Benchmarking JAX](https://docs.jax.dev/en/latest/benchmarking.html).

At 131,072 outer observations per chain, four 96-residue uint8 state arrays occupy approximately 50 MB before energies/events. If one incorrectly interprets an observation as an 88-site sweep while retaining every attempt, states alone approach 4.4 GB. Estimate memory and disk before execution; streaming is preferable to deleting repeated states.

Do not introduce reduced precision, approximate target couplings, custom kernels or new surrogate training to rescue the first result. Proposal approximations can be future experiments only with a precisely matched correction and their own checks. A GPU is an experiment after CPU profiling, not an assumed improvement.

If encoding fails correctness, pause performance work and isolate the smallest failing case. If no credible end-to-end benefit appears, apply the existing **one-week further-tuning cap** and finish with direct sampling and the measured explanation. If the direct path is also slow to mix, consider one verified stronger mixing baseline within the remaining budget; otherwise report the unresolved regime. Protect reproduction and the report before adding optional methods.

## First concrete handoff

The biology lead uses the [full-project guide](biology_lead_guide.md) and [companion worksheet](biology/benchmark_brief_template.md) to record contributions through all 18 weeks, including benchmark preparation, output review, evaluation, reproduction and final writing. The systems lead's first implementation task is the support-aware loader plus exact conditioning and independent full-score checks. These can proceed in parallel with human review of the already documented biological questions.

At the next actual joint review, bring one supported input, one deliberately invalid input and a readable coordinate/restriction table. Completion means both leads can explain what is being sampled, why those restrictions matter and which claim the first timed experiment can resolve. Record that review when it happens.
