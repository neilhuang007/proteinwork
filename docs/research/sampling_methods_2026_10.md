# Sampling mechanism: research and implementation recommendation

Research checked on **2026-10-10** against primary papers, publisher pages, official documentation and upstream repositories. This is a proposed implementation and assessment, not evidence that a production sampler or speedup already exists. Read with the [working plan](../project_plan.md), [Week 3 audit](../week03/research_audit.md) and [benchmark contract](../week04/benchmark_contract.md). This was a targeted literature search, not an exhaustive review.

## Decision

Build the optimized direct categorical sampler first. Then implement one narrowly specified domain-wall proposal with a mathematically justified correction. Treat THRML as an optional backend behind that specification. The useful deliverable is a trustworthy comparison and a rule that chooses the faster verified method; a thermodynamic-hardware advantage is outside the evidence this project can produce.

The likely risk is that encoding increases work and slows mixing on an ordinary CPU or GPU. The opportunity is to learn precisely when a reusable encoded proposal helps after all setup, rejection and reporting costs. No retrieved source establishes that it will help this 96-position protein model.

## What the current sources actually establish

| Source and status | Verified result and project decision |
| --- | --- |
| [Chancellor, 2019; accepted journal version v4](https://arxiv.org/abs/1903.05068v4) | Domain-wall encoding represents a categorical variable using one fewer binary variable than categories and supports pair interactions. Encoding is prior work. Its optimization/embedding results do not establish better equilibrium protein sampling. Reuse the audited codec/compiler. |
| [Liu et al., cumulative self-learning Monte Carlo, 2016 preprint / 2017 paper](https://arxiv.org/pdf/1611.09364) | Many cheap reversible surrogate updates can form one proposal followed by correction against the original target. This is direct prior art for the mathematical pattern. Its benefits depend on an effective model being cheaper and useful; the original work's expensive fermion calculations differ from this small, already quadratic Potts model. Do not add surrogate training to this project. |
| [THRML upstream](https://github.com/extropic-ai/thrml) | Public JAX software provides graphical-model and block-Gibbs primitives, emphasizing sparse heterogeneous graphs. The README explicitly separates software prototyping from future hardware. It is not a ready-made corrected protein sampler. Pin a tested commit and prove the schedule used by the adapter. |
| [Energy-efficient Codon Optimization, 2026-06-15, v1](https://arxiv.org/html/2606.17327v1) | A relevant domain-wall application, but it optimizes synonymous codon choices using annealing and penalty adaptation. It does not test stationary sampling of changed amino-acid sequences under a fixed trained model. Its direct Potts implementation needs far fewer sweeps than its Ising encoding. Large chip energy advantages are estimates; the paper says categorical hardware does not yet exist. Its sparse graph and GC adaptation cannot be imported as correctness evidence for our dense fixed-target sampler. |
| [Thermalizing Stochastic Programs, v2, 2026-08-13](https://arxiv.org/html/2608.01615v2) | **The cited paper exists**, and v2 supersedes the v1 cited in Week 3. It approximately compiles stochastic kernels and studies accumulated error. All reported demonstrations execute in software. The enumerable meta-EBM example has residual stationary bias. This is related compilation work, not proof of exact target preservation or an available device speedup. |
| [Torx upstream](https://github.com/extropic-ai/torx) and [Thermalizers repository linked by the paper](https://github.com/extropic-ai/thermalizers) | Torx is publicly accessible. The Thermalizers link returned **HTTP 404 in this check**, so its public availability could not be verified; a paper citation is insufficient to make it a required dependency. This does not establish whether private or renamed code exists. Neither library is needed for the minimal sampler. |
| [Grathwohl et al., ICML 2021](https://proceedings.mlr.press/v139/grathwohl21a.html), [official code](https://github.com/wgrathwohl/GWG_release) | Gibbs-with-Gradients provides established informed discrete MH proposals, and the official repository includes Potts experiments. A sensible optional stronger baseline if direct Gibbs mixes poorly. A quadratic Potts score permits exact local changes, so blindly introducing autodiff is unnecessary. |
| [Preconditioned Discrete-HAMS, published 2026-06-10](https://link.springer.com/article/10.1007/s11222-026-10907-8) | Newer peer-reviewed work uses Gaussian auxiliary variables, second-order information and generalized detailed balance; it reports discrete benchmark improvements and a rejection-free quadratic special case. This is relevant competition, but its lattice formulation and augmented momentum rules need a separate adaptation to arbitrary categorical amino acids. Defer implementation unless the minimum project is complete. |
| [Kumar and Deb, categorical MRF discrete diffusion, 2026-10-01, v1](https://arxiv.org/abs/2610.02128v1) | A very recent preprint gives learning and sampling error bounds for a weight-sharing score learner under local dependence. It requires learned scores and discretized reverse diffusion. It is useful future-work context, not a drop-in exact replacement for one immutable pretrained Potts target. No implementation or protein advantage is asserted here. |

Newest does not mean best suited. Existing exact conditionals and established MH theory are enough to make the first experiment implementable.

## Minimal mechanism, specified precisely

The following is a project-specific construction using standard MH and cumulative reversible proposals, not a new sampling theorem. The underlying acceptance principle is in [Hastings, 1970, section 2.2](https://hedibert.org/wp-content/uploads/2013/12/1970Hastings.pdf); the surrogate-proposal pattern appears in [Liu et al., equations 4–6](https://arxiv.org/pdf/1611.09364). General auxiliary-variable formulations are discussed by [Andrieu, Lee and Livingstone](https://arxiv.org/abs/2012.14881).

### 1. One immutable target

Let `D_i` be the explicit allowed symbols with finite original fields at position `i`. Use the repository's once-counted pair convention:

```text
E(x) = -sum_i h_i(x_i) - sum_(i<j) J_ij(x_i, x_j)
pi(x) proportional to exp(-beta * E(x)), x in product_i D_i
beta = 1 for the proposed benchmark
```

Read negative-infinity fields as exact exclusions under the verified source convention. Never multiply one-hot zeros by infinity. Reject an empty domain and any explicitly requested unsupported symbol. Preserve the original parameter hashes and reference map.

Eliminate fixed positions algebraically: add their interactions to the fields of variable positions and retain their constant energy contribution. This leaves the restricted distribution unchanged. Compare the reduced score with the full 96-position score, including all fixed-variable interactions, before optimizing.

### 2. Direct categorical Gibbs baseline

Choose a mutable site with a fixed, state-independent probability. For every allowed residue `a`, calculate its exact local score and sample the normalized conditional:

```text
logit_i(a) = beta * (h_i(a) + sum_(j != i) oriented_J_ij(a, x_j))
P(x_i = a | x_-i) = softmax(logit_i)[a], a in D_i
```

Subtract the largest logit before exponentiation. The original value may be drawn again; retain that repeat. A deterministic sweep of these correct conditionals also preserves the target and may be a faster direct baseline. It need not be reversible because it is used directly, without a surrogate MH shortcut.

### 3. Encoded proposal

Give each domain a fixed recorded order. Represent a site with `q_i` choices by `q_i - 1` bits in a prefix-one word. Singletons need no bits. Let `c(x)` be this one-to-one valid code and `V` the set of valid words.

Build a finite surrogate binary energy `E_s(z)` over **all** binary words using the audited domain-wall compiler and a fixed finite penalty. Initially use its exact valid-state energy construction. Do not claim that any finite penalty makes invalid states impossible.

For the first implementation, choose a bit uniformly at random and resample it from its exact binary Gibbs conditional under `mu(z) proportional to exp(-beta_s * E_s(z))`. Call this random-scan kernel `R`. Each bit update is reversible under `mu`; a fixed state-independent mixture is reversible too. Therefore `Q = R^m` is reversible for a prechosen positive integer `m`. The chain need **not** equilibrate to `mu` before producing a proposal.

From the current supported protein sequence `x`:

1. Start the inner chain at `z = c(x)`.
2. Perform exactly `m` inner updates with fresh random keys.
3. If its endpoint `z'` is invalid, reject the entire attempt and retain `x`.
4. Otherwise decode `y` and calculate original-model and surrogate energy differences.
5. Accept `y` using the log probability below; otherwise retain `x`.

```text
log_alpha = min(0, -beta * (E(y) - E(x))
                  + beta_s * (E_s(c(y)) - E_s(c(x))))
accept iff log(U) < log_alpha, U uniform on (0, 1)
```

Why the shortcut is valid: reversibility gives `Q(c(y), c(x)) / Q(c(x), c(y)) = mu(c(x)) / mu(c(y))`. Substituting that ratio into MH gives the displayed expression. View the target on binary words as `pi(decode(z))` on `V` and zero elsewhere; invalid endpoints consequently have zero acceptance. Finite penalties affect efficiency, not this target argument.

**Crucial special case:** if `E_s(c(x)) = E(x) + constant` and `beta_s = beta`, every valid endpoint has acceptance probability one. The only rejection is for invalid endpoints. Nontrivial original-energy correction appears when the proposal uses a different temperature or a genuinely approximate surrogate. Neither acceptance one nor energy equality establishes good mixing.

Fix `m`, penalty, category ordering, bit probabilities and surrogate temperature within each run. Do not redraw until valid, repair a word, stop the inner chain when it first reaches a desirable endpoint, or retain an invalid inner state as the next protein chain state. Those changes define different proposal laws.

### 4. A faster schedule only after the simple one passes

Independent sites or independent bits in a graph color can be sampled together. Coupled variables cannot be independently resampled in parallel and called exact block Gibbs. Dense protein couplings may leave little graph-color parallelism; parallelize independent chains first.

An ordinary forward sweep `F = K_1 ... K_B` preserves `mu`, but generally is not reversible. Its adjoint is the reverse sweep `F* = K_B ... K_1`. An equal, state-independent mixture `Q = (F + F*) / 2` **is** reversible. Thus a fair coin selecting forward versus reversed block order at each proposal is a defensible THRML optimization, provided every constituent block conditional and the realized ordering are verified. A palindromic schedule is another option, at additional cost.

An arbitrary schedule requires a correct augmented/pathwise MH derivation. Its reverse must have nonzero selection probability, and the ratio must include schedule-selection probabilities and all reversed transitions. A logged forward path probability alone is not the endpoint proposal probability. Keep this out of the minimum implementation; first use the reversible construction above.

The upstream [THRML adapter audited in this repository](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/src/gibbsiq/thrml_runtime.py) is useful integration material, but its ordinary block sweeps and replica swaps do not establish this categorical correction.

## Highest-value optimizations

1. **Reduce the conditional model exactly.** Remove singleton sites, fold fixed-variable interactions into fields, and pack supported domains. Keep a scalar full-score oracle for verification. Do not threshold couplings in the target.
2. **Use incremental energies.** One changed categorical site changes only its field and incident pair terms. A local score calculation costs approximately `O(L*q)`; a full score costs `O(L^2)`. Caching all conditional scores can move repeated work into an `O(L*q)` cache update after a changed residue. Confirm periodic full-score agreement to detect drift or indexing errors.
3. **Measure encoding expansion before compiling.** The number of bits is `sum_i (q_i - 1)`. Each categorical edge may expand into up to `(q_i - 1)*(q_j - 1)` binary interactions. Record actual bit count, edge count, degree, color count and memory; fewer bits than one-hot encoding does not mean fewer variables than categorical sampling.
4. **Compile chunks of transitions.** Use an outer `jax.jit` and a bounded `jax.lax.scan` for sequential work; `scan` lowers to a loop rather than unrolling large Python loops. Batch the four independent chains with `vmap`. Avoid Python callbacks and per-update device transfers. [JAX scan documentation](https://docs.jax.dev/en/latest/_autosummary/jax.lax.scan.html).
5. **Keep all outer-chain observations.** Write chunked arrays of every post-attempt protein state, including warmup and rejections, plus energies and events. Store enough inner-event information to replay and audit the proposal. Clearly distinguish an outer proposal attempt from its `m` binary updates. Inner surrogate states are not samples from the protein target. Compression may reduce disk usage; thinning or deduplication changes the required record.
6. **Use explicit numerical and randomness settings.** Enable binary64 before array creation and record both requested and effective dtypes. Split or fold independent keys by chain and step; do not reuse a random key. Persist the RNG algorithm and package versions. [JAX randomness documentation](https://docs.jax.dev/en/latest/random-numbers.html).
7. **Time the application.** Synchronize with `block_until_ready`; report compilation, transfer, tuning and writing separately and include them in primary total time. A small model can favor CPU because GPU overhead dominates. Report cold and cached execution separately; do not silently omit compilation. [JAX benchmarking guidance](https://docs.jax.dev/en/latest/benchmarking.html).
8. **Choose a supported execution environment.** Current JAX documentation lists CPU use on Windows, no native Windows NVIDIA GPU support, and experimental GPU support through WSL2. Use a pinned CPU environment first or a verified Linux GPU environment; do not budget a native-Windows CUDA speedup. [JAX installation matrix](https://docs.jax.dev/en/latest/installation.html).

Category reordering, proposal temperature and shorter or longer inner trajectories are development-only tuning choices. A smaller penalty may increase invalid endpoints; a larger one may inhibit useful moves. There is no universal best penalty or ordering demonstrated for this model. Do not introduce approximate sparse couplings unless they are confined to a clearly identified surrogate and the correction still uses the complete original energy.

## Tests required before timing a protein result

The existing four-state CM cases are useful score fixtures, but **their two mutable binary domains have no invalid domain-wall words**. They cannot validate invalid-endpoint handling or an amino-acid ordering effect.

- Enumerate synthetic three- and four-category sites and interacting pairs with unequal fields, pair energies, restrictions and at least one singleton.
- Verify codec bijection, valid-state energy differences and fixed-site reduction against independent literal sums.
- Build the complete inner binary transition matrix and the corrected outer categorical matrix. Check nonnegativity, row sums, `pi P = pi`, and detailed balance when claimed, at the contract tolerances.
- Include `m = 1` and `m > 1`, several finite penalties, equal and unequal surrogate temperatures, and a deliberate surrogate perturbation. An always-accept test does not exercise MH correction.
- Demonstrate that a forward-only sweep plus the shortcut can fail; demonstrate that the forward/reverse mixture passes. This guards the actual schedule assumption.
- Verify reachability on each tiny model. In the ideal finite Gibbs construction, staying put has positive probability and adjacent valid category changes connect the Cartesian product. Detect numerical underflow that turns a mathematically possible transition into an implemented zero.
- Force an invalid endpoint, a valid MH rejection, and an accepted move. Verify exact repeated states, energy consistency, attempt counts and deterministic replay.
- Then run both exact CM conditional cases and independent full-score comparisons. Only after these pass should the broad-model performance comparison begin.

These are proposed tests; the current research task did not implement or run them.

## Fair baselines and planner

The essential comparison is optimized direct categorical Gibbs, a fixed corrected encoded strategy, and the frozen planner choosing between them. If “planning disabled” is identical to the fixed strategy, label it as the same configuration rather than another algorithm. Give each method the same hardware, precision, tuning allowance, initialization rules and reporting costs.

If both essential methods mix poorly, the first optional extra is an informed categorical MH baseline such as locally balanced/GWG-style moves, or small exact categorical block updates chosen on development evidence. Add only one, with its own derivation and tiny-model tests. Replica exchange is another possible reference rescue, but it introduces ladder tuning and additional chain cost; it must preserve and report the declared target-temperature marginal. Annealing and selecting the best sequence do not answer the fixed-temperature sampling question.

A small fixed lookup or threshold rule is sufficient: use the case's mutable-site count, domain sizes and actual encoded size; select only among configurations that passed correctness. Calibrate on development cases, freeze before evaluation, and default to direct when evidence is weak. Never switch methods according to the current sequence without deriving the resulting adaptive/state-dependent kernel.

The present contract has one broad development and one broad evaluation case. That supports a proof of concept, **not a generalizable learned planner**. To claim useful method selection over diverse inputs, approve a larger disjoint restriction suite before tuning, varying mutable-site count and domain size with biology-approved cases. Keep the same immutable protein model. This would be a recorded benchmark revision, not an automatic change to the existing contract.

## Stop/go and impact

Use the existing accuracy, diagnostic and runtime contract as the primary gate. A speed comparison is meaningful only when both methods meet the same accuracy requirement. A timeout or inadequate long reference is an unresolved result, not evidence that one method samples correctly.

Proceed with encoding if its exact transition checks pass and a bounded development pilot shows a plausible gain in total time to accuracy. Report invalid-endpoint rate, valid-endpoint MH acceptance, movement per accepted attempt, observable-specific ESS per second, and end-to-end time; acceptance alone is not a decision metric. Consider a predeclared practical margin, such as at least 20% lower end-to-end time with repeated timings, **only as a proposed decision threshold requiring review**, not a published standard or achieved result.

If the encoded method does not help, apply the existing one-week tuning cap and publish the negative result with the direct fallback. If the broad reference cannot be qualified, retain the exact-model correctness results and report broad performance as inconclusive. Do not relax thresholds after viewing evaluation outcomes.

Expected impact is strongest as reproducible computational methodology: a reliable target-preserving workflow, explicit failure cases for encoding/schedules, and a measured choice of execution strategy. Novelty is modest because encoding, MH correction, informed proposals and kernel compilation are established. A defensible positive research claim would need a meaningful held-out advantage beyond optimized direct sampling, preferably across a larger frozen restriction suite. Faster sampling of a trained statistical model does not itself establish better enzyme function, and this project cannot claim biological or hardware energy gains without independent evidence.
