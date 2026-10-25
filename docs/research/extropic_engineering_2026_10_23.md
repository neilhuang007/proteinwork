# Extropic engineering and a practical protein optimization workflow

Research checked 23 October 2026. This note combines primary source checks with the repository interfaces available at the start of this task. The algebra below is a project derivation. Proposed interfaces and validation gates are recommendations, not results from a device or a biological experiment.

## Recommended scope

Build a complete computational loop that accepts a protein design task, proposes valid sequences, learns from supplied assay measurements, selects a diverse experimental batch, exports a reviewable library, and ingests the resulting measurements. Extropic's useful role is a replaceable Ising sampling or optimization component. The surrounding workflow supplies the biological meaning, constraints, predictor and experimental record.

Maintain two separate result types. A sampling experiment asks whether observations follow a specified distribution. A design experiment asks whether selected sequences improve a measured objective within an assay budget. Annealing, selecting the lowest energy sample, removing duplicates, repairing a candidate, and favoring predicted fitness can be appropriate design operations, but their outputs are not a stationary trace from the original Potts target.

The present repository supplies a pinned 96 position chorismate mutase Potts model, support and coordinate evidence, a small official THRML CPU pilot, and an audited domain wall compiler reuse exercise. Its [benchmark configuration](../../configs/benchmark.json) explicitly disables experimental label evaluation because a suitable independent assay holdout has not been established. The current coordinate anchor also has infinite model energy. Preserve that evidence: use a supported computational state for sampling, and never turn an energy difference from the unsupported anchor into a fitness estimate. See the [existing readiness review](repository_readiness_2026_10.md) and [sampling derivation](sampling_methods_2026_10.md).

## What Extropic currently provides

1. **Physical prototypes exist.** Extropic reports shipping XTR-0 systems built around its X0 chip and operating an experimental remote cluster. Its August 2026 announcement describes Z1 as a sparse Ising system with 269,568 pbits, sixteen neighbors per pbit and chromatic Gibbs updates. Those are vendor specifications and plans, not measurements made by this project. [Extropic announcement](https://extropic.ai/writing/from-one-to-one-billion).
2. **The public API is a simulator.** The same announcement identifies the early access API as a GPU simulator for Torx and THRML, with later Z1 cluster access planned. The October 2026 update still calls it a simulator API. A console URL does not establish device access. [August announcement](https://extropic.ai/writing/from-one-to-one-billion), [October update](https://extropic.ai/writing/baby-thermo-rsi).
3. **Scaled products remain future access.** The hardware page lists Z1 Stick and Z1 Card early access in 2027. It describes XTR-0 as an experimental platform connecting Extropic chips with a conventional processor. [Official hardware page](https://extropic.ai/hardware).
4. **THRML is usable now.** The public library supplies JAX graphical models and block Gibbs tools. Its README positions them for software prototyping and later hardware experimentation. The project's [existing pilot](../../examples/week02/thrml_pilot.py) already uses the relevant API on CPU. Installing THRML does not establish a connection to physical Extropic silicon. [THRML repository](https://github.com/extropic-ai/thrml).
5. **Additional frameworks do not remove the integration work.** Torx represents differentiable stochastic programs. Thermalizers is described as a compilation layer using trained energy models. The published Thermalizers repository URL returned HTTP 404 during this check, so this note does not recommend depending on it. [Official announcement](https://extropic.ai/writing/from-one-to-one-billion), [checked repository URL](https://github.com/extropic-ai/thermalizers).

No verified physical device endpoint, authenticated SDK session or device execution record was supplied in the inspected project materials. The current implementation can truthfully support `thrml_simulator_cpu` and, when available and tested, `thrml_simulator_gpu`. A future `extropic_device` backend must fail with a clear unavailable status until access and capabilities are actually configured. It must not silently execute a simulator and retain the device label.

Before enabling that backend, obtain a vendor capability record covering the exact device, addressable topology, coefficient ranges, quantization, programmable inverse temperature, clamping, update order, initialization, readout, batch limits and runtime API. The sources checked here do not establish an implementable public device API or numeric coefficient limits for this workload. Record unknowns as unknowns.

## A workflow researchers can use

The following separation makes each result interpretable and lets the Extropic component change without rewriting the scientific workflow.

1. **Validate the task.** Read the reference construct, residue map, allowed residues, fixed residues, gap policy, mutation budget, objective direction, assay context and requested batch size. Reject conflicting definitions and unsupported residues before running a sampler. Distinguish the full protein sequence, aligned model sequence and actual measured construct.
2. **Build a candidate pool.** Reuse the supported Potts model as one proposal source. Also allow imported FASTA candidates and a simple legal mutation generator. Evaluate restrictions against the decoded full sequence. Keep model energy, candidate source and original sample identifiers. An existing fitted model is useful prior information; its energy is not an assay measurement.
3. **Fit an assay predictor only to available observations.** Train on accepted measured rows from earlier rounds. Store the feature transform, target scaling, split and model hash. A measured reference or control may fall outside the statistical Potts support. Accept such evidence when its aligned length, alphabet and construct gap pattern match, and annotate that relationship. The delivered predictor uses categorical sequence features, so an undefined Potts energy never enters training. With no measurements, export a prior guided batch and label it accordingly. With too little or constant data for a meaningful fitted model, expose that limitation and use an explicit fallback.
4. **Evaluate acquisition on the finite pool.** Calculate predicted mean, declared uncertainty and a reproducible utility, for example `u_i = mean_i + kappa * uncertainty_i` for a maximization task. A scalar uncertainty bonus does not automatically produce a statistically optimal batch. Store the formula, units, scaling and exclusions.
5. **Select the batch.** Exclude sequences already measured or pending unless deliberate replication is requested. Compare top utility, greedy diversity and an Ising formulation over the same pool. Keep a quota for explicit controls or replicates separate from novel design slots. The final batch must contain the requested number of valid distinct designs or return an honest shortfall.
6. **Export and ingest.** Write a candidate CSV, FASTA, assay worksheet and run manifest. Import later measurements by immutable sequence identifier, preserving raw values, missing or censored status, replicate and batch identifiers. Refit after ingestion and produce the next round. Never fill missing assay values with model predictions.

This structure follows the useful separation between models, measured sequences, exploration algorithms and assay budgets in FLEXS. Its explorer receives prior measured sequences and returns proposed sequences with model scores. [FLEXS explorer source](https://raw.githubusercontent.com/samsinai/FLEXS/master/flexs/explorer.py).

## Exact conditioning and categorical encoding

Use the repository convention, with each pair counted once:

```text
E(x) = -sum_i h_i(x_i) - sum_(i<j) J_ij(x_i, x_j)
pi(x) proportional to exp(-beta * E(x))
```

Define each allowed domain from both biological restrictions and finite original support. A field of negative infinity remains an exact excluded choice. Do not evaluate a dense one hot product containing `0 * -inf`.

For a fixed position `j` with residue `r_j`, add `J_ij(a, r_j)` to each mutable field `h_i(a)`, with the appropriate pair orientation. Fixed fields and fixed pair terms become a scalar constant. This is exact substitution. The reduced energy plus the saved constant must equal the original full energy for every test sequence. Removing fixed positions also reduces encoding size and programming cost.

For one hot variables `z_(i,a)`, the finite supported energy is quadratic:

```text
E_valid(z) = constant - sum_(i,a) h_reduced_i(a) * z_(i,a)
            - sum_(i<j,a,b) J_ij(a,b) * z_(i,a) * z_(j,b)
sum_a z_(i,a) = 1
E_penalized(z) = E_valid(z) + sum_i lambda_i * (sum_a z_(i,a) - 1)^2
```

The equality to the protein energy applies on valid words. A finite penalty still gives invalid words nonzero Boltzmann probability. Increasing the penalty can suppress invalid words while creating barriers between legal categories. If every single bit move that leaves the one hot set is rejected immediately, the chain cannot change category at all: changing one category requires two bit changes. Use a proper categorical update, a valid swap proposal, or the previously derived corrected encoded construction.

For this repository, prefer the already audited `CategoricalModel` and `compile_domain_wall` route before introducing another protein compiler. A domain with `q` choices uses `q - 1` bits. Its existing reuse test checks valid state energies and invalid word rejection, but does not establish a production sampler. Retain the [pinned reuse boundary](../week03/research_audit.md) and [reversible proposal requirements](sampling_methods_2026_10.md).

THRML also supports clamped blocks directly. `IsingSamplingProgram(model, free_blocks, clamped_blocks)` defines the blocks, and `sample_states` receives clamp state arrays separately in matching order. Use exact algebraic elimination where convenient, and test it against explicit clamping on tiny models. The free blocks used with independent spin conditionals must contain no coupled variables. A complete interaction graph generally requires singleton blocks. [THRML block API](https://docs.thrml.ai/en/latest/api-block-sampling.html), [official introductory example](https://docs.thrml.ai/en/latest/getting-started.html).

## A concrete Ising task for selecting an experimental batch

An arbitrary neural predictor, Gaussian process acquisition or uncertainty function generally cannot be transferred unchanged into a pairwise Ising model. Instead, score a finite candidate pool classically and compile the batch selection problem. This is an exact transformation of the declared selection objective, independent of how candidate utilities were obtained.

Let `z_i` indicate whether candidate `i` is selected, let `k` be the required batch size, let `u_i` be its utility, and let `r_ij` be a symmetric similarity cost. A possible similarity is sequence identity across the mutable positions, with the exact definition stored in the manifest. For nonnegative `gamma`, minimize:

```text
E_batch(z) = -sum_i u_i*z_i
             + gamma*sum_(i<j) r_ij*z_i*z_j
             + A*(sum_i z_i - k)^2
```

In the convention `E = c + sum_i a_i*z_i + sum_(i<j) b_ij*z_i*z_j`, expansion gives:

```text
c = A*k*k
a_i = -u_i + A*(1 - 2*k)
b_ij = gamma*r_ij + 2*A
```

A sufficient conservative optimization penalty is:

```text
A > max_i(abs(u_i) + gamma*sum_(j!=i) abs(r_ij))
```

This bound is derived by changing any selected bit toward cardinality `k`. The unpenalized objective can increase by at most the expression on the right, while the cardinality penalty decreases by at least `A`. Thus no invalid cardinality can be a global minimum, provided `0 <= k <= n`. This is a statement about optima, not finite temperature sample validity. Avoid applying a numerical margin smaller than the coefficient precision.

At finite temperature, inspect the cardinality of every returned word. For an optimization output, keep the best feasible result found and report valid fraction and objective. A repair or classical fallback may be useful, but label it, time it and compare it independently. It must not be presented as an unmodified device sample or as evidence of exact constrained sampling. Include the cases `k = 0`, `k = n` and an insufficient candidate pool explicitly.

This QUBO has an important hardware limitation: the squared cardinality term connects every candidate pair. Even a small candidate pool may exceed a device's native degree. Software execution is not proof of a valid Z1 placement. Reducing the pool, decomposing the problem or adding auxiliary variables are possible future strategies, each with a separate correctness and cost analysis.

## QUBO, Ising signs and temperature

With `z_i = (s_i + 1)/2`, where `s_i` is minus one or plus one:

```text
H_i = a_i/2 + sum_(j!=i) b_ij/4
K_ij = b_ij/4
c_spin = c + sum_i a_i/2 + sum_(i<j) b_ij/4
E = c_spin + sum_i H_i*s_i + sum_(i<j) K_ij*s_i*s_j
```

THRML's `IsingEBM` uses dimensionless energy `-beta * (sum bias_i*s_i + sum weight_ij*s_i*s_j)`. Therefore the adapter supplies `bias_i = -H_i` and `weight_ij = -K_ij`. The constant does not affect sampling, but retain it for arithmetic comparison. THRML stores spin state as Boolean values, with `False` corresponding to minus one and `True` to plus one. [Ising source](https://raw.githubusercontent.com/extropic-ai/thrml/main/thrml/models/ising.py), [spin factor source](https://raw.githubusercontent.com/extropic-ai/thrml/main/thrml/models/discrete_ebm.py).

Use `dimod.qubo_to_ising(Q, offset)` or `BinaryQuadraticModel.from_qubo` as an existing compiler or independent arithmetic oracle. Its Ising energy has positive coefficient sums, so negate its returned `h` and `J` when populating THRML. Store only one orientation per pair in a canonical dictionary. Passing both matrix triangles without adjusting coefficients doubles interactions. [Official conversion documentation](https://docs.dwavequantum.com/en/latest/ocean/api_ref_dimod/generated/dimod.utilities.qubo_to_ising.html), [BQM constructor](https://docs.dwavequantum.com/en/latest/ocean/api_ref_dimod/generated/dimod.binary.BinaryQuadraticModel.from_qubo.html).

Scaling all energies by `1/C` preserves the optimizer's minimizers when `C > 0`, but preserves the same sampling distribution only if inverse temperature changes from `beta` to `beta*C`. A sampler with fixed physical noise or limited programmable temperature may be unable to compensate. Clipping or quantizing individual coefficients changes the energy itself. Record the intended coefficients, realized coefficients, scaling and effective inverse temperature separately. Do not silently claim that normalization preserves a sampling target.

Large penalties can overwhelm useful energy differences. Check maximum coefficient, smallest retained nonzero coefficient, ratio of magnitudes, maximum local field, and observed numerical saturation. Exact zeros from rounded conditional probabilities can destroy transitions that exist mathematically. A prospective device adapter needs calibrated realized conditionals before applying any exact Metropolis correction argument; the intended Ising energy alone does not establish the proposal law of an imperfect physical device.

## Existing libraries to reuse

* **THRML and JAX:** Keep the installed and tested THRML API rather than inventing `chip.sample()` or assuming a remote backend. Reuse the existing pilot imports, state ordering and synchronization. Record the installed version and source hashes. [THRML repository](https://github.com/extropic-ai/thrml).
* **scikit-learn:** `BayesianRidge.fit(X, y)` and `predict(X, return_std=True)` supply a manageable first surrogate. Its returned standard deviation is predictive uncertainty, including estimated observation noise. `sigma_` stores coefficient covariance. If acquisition is intended to use uncertainty about the latent mean, calculate and label that quantity explicitly with the same feature transform. Linear feature models are a baseline, so validate them on held out sequences and report when epistasis or extrapolation defeats them. [API](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.BayesianRidge.html), [implementation](https://raw.githubusercontent.com/scikit-learn/scikit-learn/main/sklearn/linear_model/_bayes.py).
* **BoTorch:** `optimize_acqf_discrete(acq_function, q, choices, unique=True, X_avoid=...)` provides an established candidate pool route and handles sequential conditioning for batches. This is a stronger optional comparator when posterior modeling justifies the dependency. It is not an exact joint combinatorial optimizer over every possible batch. [Official optimizer API](https://botorch.readthedocs.io/en/stable/optim.html#botorch.optim.optimize.optimize_acqf_discrete).
* **FLEXS:** Reuse its model, explorer and measured round abstractions. ADALEAD, random search and other explorers can supply later comparative baselines. Direct installation needs a separate compatibility check: the current source includes `DataFrame.append`, and its setup lists TensorFlow, tf-agents, PyTorch and TAPE dependencies. Do not assume that broad historical dependency ranges make it compatible with this Python 3.13 environment. [Repository](https://github.com/samsinai/FLEXS), [explorer source](https://raw.githubusercontent.com/samsinai/FLEXS/master/flexs/explorer.py), [dependency declaration](https://raw.githubusercontent.com/samsinai/FLEXS/master/setup.py).
* **ProteinGym:** Its processed assay tables provide full mutated sequences, mutation descriptions and experimental `DMS_score`, with assay references and target constructs separately recorded. These are useful for retrospective selection experiments when the chosen assay matches the sequence task. They are not measurements of arbitrary newly generated sequences. Use only compatible assays with documented separation between training, tuning and evaluation. [Official resource schema](https://github.com/OATML-Markslab/ProteinGym).
* **Extropic codon optimization:** Its public implementation contains direct Potts and domain wall Ising paths built on THRML. It explicitly describes its current JAX execution as simulation and its energy advantage as an estimate grounded in prototype information. It optimizes synonymous codons while preserving the amino acid sequence, so reuse engineering patterns without treating its results as evidence of improved protein function or a measured speedup here. [Official companion code](https://github.com/extropic-ai/codon_opt).

## Required data and output contract

For a useful first release, the task manifest should contain:

```text
task_id, objective_name, direction, units, assay_context
reference_id, reference_sequence, measured_construct, position_map
model_source, model_hash, alphabet, support, allowed_residues, fixed_residues
mutation_budget, gap_policy, batch_size, control_and_replicate_policy
candidate_source, pool_budget, predictor_configuration, acquisition_configuration
selection_configuration, random_seeds, execution_backend, run_budget
```

The measurement table should contain:

```text
sequence_id, sequence, assay_id, value, units, value_status
replicate_id, experimental_batch, measurement_uncertainty, measured_at, round
```

Reject nonfinite measured values except as explicitly represented missing or censored records. Validate sequence identity when joining; matching an identifier alone is insufficient if the supplied sequence changed. Preserve raw records and record any aggregation or normalization. A prediction column must remain separate from `value`.

The exported candidate table should contain:

```text
sequence_id, sequence, aligned_sequence, mutations, candidate_source
model_energy, predicted_mean, uncertainty, uncertainty_kind, acquisition
selected, selection_rank, selection_reason, feasibility, source_run_id
```

Accompany it with FASTA, a blank measurement template and a run manifest containing input hashes, dependency versions, transformations, backend identity, solver settings, timings, feasibility counts, fallback use and output hashes. Keep raw sampler observations separately when making distributional claims. These proposed schemas are project interface recommendations, not a vendor SDK.

## Fair tests and useful claims

1. **Arithmetic gate:** Exhaustively compare original categorical, QUBO and Ising energies on small models, including invalid words, asymmetric pair tables, fixed sites and forbidden symbols. For batch selection, enumerate all small binary words and all feasible batches. Check both the penalty guarantee and exact conversion constants.
2. **Sampling gate:** Preserve the current fixed target contract, multiple chains, exact small distributions, rejected state repeats and uncertainty criteria. A deterministic THRML sweep preserves its declared Gibbs target when its conditionals are correct, but is generally not reversible. Apply a surrogate energy shortcut only with the proven reversible schedule described in the existing sampling note.
3. **Selection gate:** Compare the same candidate pool, utilities, similarity definition, required batch size and total budget. Include top utility, greedy diversity, uniform random feasible batches and exhaustive feasible search when tractable. Report objective gaps, valid fraction, diversity and runtime. If THRML results are repaired or improved by a classical solver, report both stages and their costs.
4. **Assay gate:** Run retrospective rounds by exposing only previously queried labels to the predictor. Match initial data and assay budgets across methods. Report best observed objective versus cumulative assays, improvement over the reference, hit rate at a declared threshold, prediction error and uncertainty calibration. A stronger selection objective without assay improvement supports a solver result only.
5. **Timing gate:** Measure input validation, fitting, candidate generation, acquisition, compilation, initialization, sampling, transfers, repair, selection and export. Report cold and cached execution separately and synchronize JAX results before stopping a timer. These requirements follow the asynchronous execution and compilation cautions in the [JAX benchmark guide](https://docs.jax.dev/en/latest/benchmarking.html).
6. **Hardware gate:** Record actual execution device and expose capability or embedding failures. Native Windows JAX supports CPU, while NVIDIA GPU execution is not supported directly and WSL2 support is listed as experimental. Use a verified environment rather than presuming GPU availability. [JAX installation support matrix](https://docs.jax.dev/en/latest/installation.html).

The immediate defensible deliverable is a reproducible sequence design loop with an actual THRML software component, explicit classical baselines and an honest device boundary. A physical chip claim requires device execution records. A speed claim requires comparable total time at comparable quality. An improved protein claim requires matching experimental evidence. None follows merely from producing a valid Ising file or a plausible amino acid sequence.

## Implementation delivered with this note

The accompanying [design module](../../proteinwork/design.py) implements strict assay ingestion, training only Bayesian ridge ranking, a deterministic greedy batch baseline, SciPy MILP selection and actual THRML CPU batch optimization. Assays outside prior support remain usable and are counted explicitly; candidates remain strictly restricted to original finite support. Predictions use chunks of 256 sequences. A conservative feature and covariance allocation estimate is checked against a 256 MiB budget before fitting. The THRML route compiles the candidate batch objective, uses singleton Gibbs blocks, and checks exact cardinality and pair conflicts on every returned word. It retains a legal greedy incumbent unless a feasible THRML sample has a strictly better objective. Returned evidence identifies the source of the selected batch, sample validity, coefficient scaling, effective inverse temperature and complete solver time. A single candidate is handled deterministically and explicitly records that no sampling was required.

Candidate optimization is capped at 128 variables generally and 64 for THRML. The THRML coefficients are divided by their largest absolute local field bound, and its configured inverse temperature is limited so the worst case absolute Gibbs logit is at most 30 in binary64. This is a declared optimization temperature, not a claim to preserve the original protein sampling target. Draw and warmup counts bound the work; the elapsed time setting is reported after a synchronous JAX call and is not an interruptible deadline. The [focused tests](../../tests/test_design.py) include exhaustive batch arithmetic, agreement with an exact feasible search, assay label leakage guards, actual THRML execution and an explicit invalid sample fallback case. These tests use synthetic arithmetic fixtures and contain no invented biological measurements.
