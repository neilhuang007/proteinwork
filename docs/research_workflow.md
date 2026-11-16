# Research workflow and engineering decisions

This release implements the complete computational path from a pinned chorismate mutase model and research constraints to a traceable candidate library and an assay data round trip. Neil confirmed chorismate mutase and THRML simulation on 23 October 2026. Physical device access and new biological measurements are not available. Existing adviser review records remain unchanged.

## Research question and outputs

A useful immediate question is: given this family model, these allowed residues, this experimental batch size and these selection criteria, which distinct CM domain sequences should a researcher review for testing? A separate methods question asks how quickly each sampler reaches the same numerical accuracy. The program keeps these questions separate so a selection heuristic cannot invalidate a sampling comparison.

The default ranking maximizes negative Potts energy. This is a model preference, not a measured activity, stability or binding prediction. An optional assay predictor learns a specified numerical readout from supplied measurements. Neither ranking guarantees a functional protein. The model covers an 89 residue domain corresponding to Q5SLA5 residues 2 through 90, not the full 355 residue accession. The exported FASTA therefore contains domain segments and does not prescribe an expression construct.

## Run the workflow

Use Python 3.13 on Windows and the complete hashed research lock:

```powershell
uv pip sync requirements/research-win-py313.txt --python .venv-pilots/Scripts/python.exe --require-hashes --link-mode=copy
uv pip install --python .venv-pilots/Scripts/python.exe -e . --no-deps --link-mode=copy
uv pip check --python .venv-pilots/Scripts/python.exe
& ./.venv-pilots/Scripts/python.exe -m proteinwork validate --request configs/cm_library_thrml.json
& ./.venv-pilots/Scripts/python.exe -m proteinwork run --request configs/cm_library_thrml.json --output output/my_cm_library
& ./.venv-pilots/Scripts/python.exe -m proteinwork verify output/my_cm_library
```

The model archive must already be cached with the pinned hash. If needed, run `scripts/inspect_benchmark.py --download` as described in the main README. Production runs never silently replace or download changed model parameters. Every run needs a new output directory.

`cm_library_thrml.json` uses the optimized direct sampler to explore all 88 mutable positions, then uses actual THRML to optimize the candidate batch. `cm_library.json` uses SciPy MILP for the batch instead. `cm_exact_direct.json` and `cm_exact_thrml.json` compare the direct and corrected encoded sequence samplers against an exactly enumerable conditional target. These small conditional cases validate the sampler; the broad cases exercise library generation.

## Input contract

The JSON request has schema version 1. All paths resolve relative to that request. Unknown top level, sampling, selection and residue rule keys fail validation. Duplicate JSON keys and nonfinite values fail rather than silently changing a rule.

| Input | Meaning and enforcement |
| --- | --- |
| `benchmark_config` | Immutable model source, archive and parameter hashes, alphabet, complete key layout, coordinate anchor and named restriction cases |
| `evaluation_config` | Four chains, float64 arithmetic, accuracy thresholds, disjoint development and evaluation cases, declared checkpoints and runtime budgets |
| `case_id` | A declared conditional target with exactly enforced original support and gap rules |
| `beta` | Inverse computational temperature, consistent with the evaluation protocol; unrelated to incubation temperature |
| `sampling` | Method, distinct seeds, warmup, draw budget, sweep count, optional encoded backend and resource limits |
| `restrictions` | Additional residue restrictions with explicit symbols, coordinate system, reason and optional expected anchor residue |
| `selection` | Batch size, objective, optimizer, diversity, motif, composition, distance and resource rules |
| `assays` | Optional strict CSV of actual supplied observations |
| `reference_run` | Optional separately qualified direct run with the identical target and evaluation rules and disjoint seeds |
| `planner` | Optional frozen selection rule learned only from eligible development evidence |
| `allow_unqualified_export` | Explicit permission to export exploratory candidates when equilibrium accuracy is unestablished |
| `purpose` | `candidate_library` by default, or an independent direct `reference` run |

For a biological residue restriction, supply `accession`, `uniprot_position_1based`, `allowed_symbols`, `expected_anchor_symbol` and `reason`. Alternatively supply `model_index_0based`. Exactly one numbering system is accepted. A request can narrow the named case but cannot widen its support. An unmapped gap, incorrect accession, duplicate site, forbidden residue or reference mismatch fails before sampling.

The sampler enforces positional support. Whole sequence motifs, mutation distances, residue composition and pairwise library distances are explicitly selection rules. They are evaluated after sampling and do not silently redefine the probability target. Motifs are literal substrings of the ungapped domain; they are not regular expressions. Hamming distances use the aligned domain.

## Every computational stage

1. **Validate source evidence.** Verify archive size and hash, parameter member size and hash, every expected field and coupling key, finite couplings and the literal meaning of forbidden fields. The unsupported coordinate anchor is never a chain initializer.
2. **Resolve biological coordinates.** Verify the exact reference substring and the supplied map. Keep gaps unmapped and preserve the distinction between the anchor and the constructed support witness.
3. **Condition the model exactly.** Remove fixed positions, add their interactions to the remaining fields and retain the fixed background energy. This changes the representation but preserves every supported full sequence energy.
4. **Initialize independent chains.** Record seeds, actual initial states, precision and backend. Randomized starts lie within the requested conditional support.
5. **Sample.** The direct method uses compiled categorical Gibbs sweeps with batched chains and cached conditional fields. The encoded method uses a one hot QUBO and an equal mixture of forward and reverse finite binary Gibbs sweeps, followed by the original target correction. Invalid endpoints are retained as rejections. See [the derivation](encoded_method.md).
6. **Measure accuracy at declared boundaries.** ArviZ computes rank normalized Rhat, bulk ESS, tail ESS and MCSE for energy and anchor distance. Site and selected pair indicators have separate uncertainty records. Small targets use exact probabilities. Broad targets require a qualified independent reference for an accuracy claim. A constant trace cannot prove a constant observable.
7. **Account for time and failure.** `checkpoint_policy=first_qualified` examines every declared checkpoint up to the draw budget, continuing each chain state with recorded independent random substreams and warming up only once. The first passing checkpoint stops the sampling experiment. `final` assesses only the requested endpoint. The CLI enforces the total wall clock cap in a child process; numerical, memory and graph failures remain failures. All available stage timings are recorded.
8. **Verify original scores.** Decode every distinct observed sequence, recompute its original full model energy and use the inverse duplicate map to check every recorded score. Check every observed residue against the final domains.
9. **Build and rank a separate pool.** Preserve source chain, draw, occurrence count and sequence hash. Deduplication occurs here, never in the trace. Exclude already measured sequences when assays are supplied.
10. **Select an experimental batch.** Enforce all individual and pairwise rules and exact batch size. MILP returns a feasible incumbent with its solver status and optimality gap. THRML samples the compiled batch objective, filters every word for exact legality and compares against a deterministic greedy incumbent. If greedy supplies the final batch, the result says so. No method relaxes constraints to fill a shortfall.
11. **Export the handoff.** Save candidate CSV, aligned and ungapped FASTA, assay worksheet, complete trace, coordinate map, applied rules, diagnostic results, timing profile and a hashed artifact manifest.
12. **Ingest measurements.** Join a completed worksheet to the immutable exported candidate identifier and aligned sequence. Preserve supplied values and experimental context. Reject split leakage and duplicate records. Reuse the resulting assay CSV in a new request for the next round.

## Assay round trip

The accepted assay CSV has exactly these columns:

```text
aligned_sequence,value,split,assay_id,units,context,round
```

`split` is `train`, `validation` or `test`. Repeated observations of the same sequence stay in one split and are averaged within that sequence; every distinct sequence has equal training weight. Identity, units and context must match the requested assay. Sequence alignment and construct gaps must match the modeled domain. A real measurement outside the Potts model's support remains measured evidence for the separate predictor and is explicitly counted; such a sequence is not admitted to the sampling target.

Only training rows fit the categorical sequence features, centering, target scale and Bayesian ridge predictor. Validation and test values are reported separately and never determine utilities. The uncertainty column is predictive standard deviation, including observation noise; it is not a claim of calibrated uncertainty. UCB uses `direction_sign * predicted_mean + exploration_weight * predictive_std`. Two distinct, nonconstant training measurements are the minimum for fitting, not evidence of a reliable predictor. More representative data and independent evaluation remain the researcher's responsibility.

The exported worksheet adds `candidate_id` for safe joining. Fill the measurement, assay identity, units, context and round, then run:

```powershell
& ./.venv-pilots/Scripts/python.exe -m proteinwork ingest --run output/my_cm_library --measurements completed_worksheet.csv --output data/assays/round1.csv
```

To append a later round, supply `--existing data/assays/round1.csv` and a new output filename. Blank measurements remain unmeasured and are omitted with a recorded count. Predictions never fill blanks. The ingestion receipt hashes the source run, supplied measurements and resulting data. No biological assay records are invented by the examples.

## Optimization audit

| Stage | Optimization and correctness boundary |
| --- | --- |
| Loading | Dense NumPy arrays and one pass over parameter records; source hashes and complete coverage remain mandatory |
| Conditioning | Exact fixed site elimination reduces positions, pair blocks and encoded spins; original energy tests guard equivalence |
| Direct sampling | JAX scan, four chain batching, cached conditional fields and periodic refresh; complete coupled sites are never updated as independent variables |
| Encoded sampling | Fixed sites and forbidden categories consume no spins; actual realized coefficients determine the correction; rounding and scale changes cannot change the target silently |
| Diagnostics | Reuse the same trace and exact target; rare unobserved indicators retain unknown uncertainty instead of false precision |
| Scoring | Score distinct full sequences once and verify repeated observations through the inverse map |
| Assay prediction | Reuse scikit learn, bounded training allocation and chunks of candidate features; no fitting to evaluation labels |
| Batch selection | Deterministic utility shortlist with an explicit pool cap, sparse MILP constraints, THRML coefficient normalization and a measured greedy baseline |
| Planning | Freeze the method and settings before evaluation; default to direct when no qualified development evidence supports a choice |
| Handoff | Compressed traces, complete provenance and separate library outputs; no implicit overwrite of prior results |

The pool cap is a declared approximation to the full selection search. MILP optimality applies only to that admitted pool. THRML does not promise an optimum. Its per call time setting is checked after the synchronized JAX call, while bounded spins, draws and warmup limit work; the overall CLI budget is hard. Current sequence THRML compilation is capped at 64 logical spins. Broad CM protein encoding is much larger, so use the JAX encoded backend for that experiment or THRML on the candidate selection problem. No interactions are silently dropped to make a graph fit.

## Artifacts and interpretation

Reference runs now skip candidate materialization and batch selection. They require `purpose=reference`, direct sampling, `batch_size=0`, the model score objective, no assays, no supplied reference and no exploratory export. All raw observations, original score verification, residue checks, diagnostics and hashes remain. Candidate CSV and worksheet files contain headers only, FASTA files are empty, and `selection.json` records `not_applicable_reference`. The manifest's `unique_candidates` still counts distinct observed sequences; `candidate_pool_materialized=false` explains why there are no candidate pool rows. Reference checkpoints stop only when the stricter reference qualification passes. Their checkpoint records identify `qualification_goal=reference`.

`trace.npz` and `trace.csv.gz` contain every retained observation, including rejections. `sampler_event` has the semantics declared by its sampler: Gibbs endpoint movement for direct sampling and MH acceptance for encoded sampling. A valid encoded self proposal can be accepted. The CSV additionally records actual movement between consecutive observations.

`candidate_pool.csv` contains candidate scores and exclusions. `candidates.csv` and both FASTA files contain only the admitted library. An insufficient pool or impossible batch produces a shortfall status and an empty library. Equilibrium qualification never implies that every possible rare sequence has appeared. `selection.json` explains the objective, finite pool, constraints, assay model, optimizer and incumbent source.

`rules.json`, `checkpoints.json`, `diagnostics.json`, `profile.json`, `encoding_resources.json` and `manifest.json` make each boundary reviewable. Output status distinguishes a qualified sampling library, an exploratory library, missing library and inadequate reference. `verify` checks saved hashes and trace counts. The manifest does not claim a Git commit when the workspace has no Git metadata.

## Comparisons, planning and hardware

Use `proteinwork compare --runs RUN1 RUN2 --output comparison.json` for matched targets, rules, code, hardware and selection settings. Runs without qualified checkpoints remain in the comparison without successful timing ratios. Use `proteinwork freeze-plan --runs DEV_RUN1 DEV_RUN2 --output plan.json` to freeze eligible development choices. Evaluation runs cannot train it. A negative result, including direct sampling winning, is retained.

Physical Extropic execution requires an identified device, supported topology, coefficient ranges, quantization, update order, initialization and authenticated runtime. The present code records logical resource needs and exact simulated behavior. It does not treat a dense Ising matrix as an embedded device program. Current source evidence and the access distinction are documented in [the Extropic engineering review](research/extropic_engineering_2026_10_23.md).
