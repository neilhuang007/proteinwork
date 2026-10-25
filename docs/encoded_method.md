# Corrected encoded protein sampling

The implementation in `proteinwork/encoded.py` compiles a supported Potts conditional distribution into an exact one hot QUBO and an Ising model. It uses a finite reversible Ising trajectory as a proposal, then applies a Metropolis Hastings correction against the original restricted Potts energy. NumPy, JAX, and the official THRML package execute the same mathematical proposal family. None of these backends uses a physical thermodynamic chip.

## Target and support

The source convention is

```text
E(x) = -sum_i h_i(x_i) - sum_(i<j) J_ij(x_i,x_j).
pi(x) is proportional to exp(-beta * E(x)) on the supported domains.
```

The model layer removes singleton sites, folds their interactions into effective fields, and retains their constant energy contribution. Every output expands back to the original protein coordinates. A symbol forbidden by an infinite source field or by the declared restriction has no binary variable. Finite encoding penalties never replace forbidden support.

Each remaining site has a recorded ordered domain `D_i`. There is one binary variable `z_ia` for each supported symbol. A word is valid exactly when each site contains one active bit. The codec is a bijection between valid words and supported protein sequences.

## Exact compiler

For a finite nonnegative penalty `lambda`, the energy on every binary word is

```text
Es(z) = offset - sum_(i,a) h_i(a) z_ia
               - sum_(i<j,a,b) J_ij(a,b) z_ia z_jb
               + lambda * sum_i (sum_a z_ia - 1)^2.
```

This defines a finite energy for invalid words as well as valid words. Every valid word has exactly the original restricted energy. The default penalty is `2.0` in source energy units; it is a fixed implementation choice, without a claim of optimal mixing.

The upper triangular QUBO uses linear coefficients `l`, pair coefficients `u`, and constant `c`. Substitution of `z = (s + 1) / 2` gives the convention used by the official THRML Ising model:

```text
Es(s) = constant - sum_k biases[k] * s[k]
                    - sum_(k<l) weights[k,l] * s[k] * s[l]
biases[k] = -l[k]/2 - sum_(edges incident on k) u[edge]/4
weights[k,l] = -u[k,l]/4
constant = c + sum(l)/2 + sum(u)/4
```

`compile_ising` retains every supported nonzero interaction. It returns the QUBO coefficients, Ising coefficients, energy constant, codec, and resource report. Exact zero couplings can be omitted because they contribute zero energy. There is no threshold based sparsification.

## Finite proposal and proof

For the implemented Ising energy `Es`, let `mu(z)` be proportional to `exp(-proposal_beta * Es(z))`. A singleton heat bath update `K_k` is reversible under `mu`. Its probability of spin `+1` is

```text
sigmoid(2 * proposal_beta * (biases[k] + sum_l weights[k,l] * s[l])).
```

Let `F = K_1 ... K_B` be one ascending singleton sweep. Its adjoint under `mu` is `F* = K_B ... K_1`, the descending sweep. With a fixed positive integer `m = sweeps_per_draw`, the proposal is

```text
Q = (F^m + (F*)^m) / 2.
```

A fair, independent coin selects one direction at the beginning of every outer attempt. That direction is repeated for exactly `m` sweeps. The adjoint of `F^m` is `(F*)^m`; therefore `Q` is reversible. No inner equilibration assumption enters this argument.

Every attempt starts at the encoding of the current supported protein state. An invalid endpoint is rejected once. The old protein state and energy are recorded again. A valid endpoint `y` is accepted with

```text
log_alpha = min(0,
    -beta * (E(y) - E(x))
    + proposal_beta * (Es(code(y)) - Es(code(x))))
```

The second term is the exact proposal ratio because reversibility gives `Q(code(y),code(x)) / Q(code(x),code(y)) = mu(code(x)) / mu(code(y))`. The target extended to invalid words has zero probability there, so their rejection is ordinary Metropolis Hastings rejection. Redrawing until valid would change this proposal and is never performed.

When the exact compiler is used without coefficient changes and both inverse temperatures agree, the two energy differences cancel for valid endpoints. All valid endpoints then have acceptance probability one, up to floating point arithmetic. This does not imply efficient movement: invalid endpoints and encoding barriers can dominate.

## Scaling, rounding, and numerical limits

The optional `coefficient_scale` multiplies Ising biases and weights. The optional `quantization_step` rounds those scaled coefficients to integer multiples of that step. These options define an explicit surrogate. The inner sampler and the correction both evaluate the resulting implemented coefficients. The original Potts target and its support remain unchanged. These options model coefficient changes in software and do not assert a calibrated physical transfer function.

All arithmetic is binary64. JAX enables and records its binary64 setting. NumPy and JAX detect any visited binary conditional that rounds to zero or one and fail explicitly. Before THRML execution, a conservative absolute field bound must keep every binary logit below `36.0` in magnitude. This can reject a feasible but numerically difficult model; it does not silently saturate it. Finite precision is still an approximation to real arithmetic, and exact transition checks use the repository's numerical tolerances.

## Official THRML integration

The adapter uses `thrml.models.IsingEBM`, `IsingSamplingProgram`, and `thrml.sample_states`. Each free block contains one spin. Two programs use opposite block orders, and the adapter checks the resulting singleton sampling order. The installed official source was inspected for the spin convention, the sigmoid conditional, the sequential block update, and the first observation semantics.

For every finite endpoint, the THRML schedule is

```python
SamplingSchedule(n_warmup=sweeps_per_draw, n_samples=1, steps_per_sample=1)
```

The official sampler records its first observation immediately after `n_warmup` sweeps. Here that parameter is the finite inner trajectory length, not an assertion of equilibrium. Outer protein warmup consists of full corrected attempts and is a separate argument. Input block states follow each program's block order; observations are returned in the common original spin order.

THRML programs are compiled before sampling. Returned arrays are synchronized before timing ends. Metadata records installed package versions, source file SHA256 hashes, actual JAX devices, effective precision, the complete finite schedule, and compilation time. The NumPy and JAX backends are named software implementations and are never presented as THRML or physical execution.

Each THRML attempt receives a separate JAX key derived from its chain seed and attempt index. A model with no logical edges receives one explicitly recorded zero weight API padding edge, because the installed library expects a nonempty edge factor. This does not change its energy or logical graph. This adapter caps official THRML compilation at 64 spins to bound its singleton program construction; larger encoded graphs require the explicitly named NumPy or JAX simulator.

## Observations and budgets

`SamplingResult.states` contains every post warmup outer attempt, including all repeated states. Warmup states, energies, acceptance flags, and invalid endpoint flags are retained in metadata. Seeds and actual initial states are returned. `accepted` means the MH decision, including accepted valid self proposals; movement must be computed by comparing successive states. Timing includes encoding, initialization, compilation, synchronized sampling, correction, transfers, and result construction within the sampler. File writing and workflow diagnostics must be added by its caller.

Compilation fails before allocating the dense interaction matrix when the requested spin count, edge count, or estimated working bytes exceed the caller's budget. Output trace size is also checked. Estimates are conservative array estimates; Python object and compiler memory can exceed them. Limits cause explicit errors rather than clipping or dropping couplings.

The resource report states logical spin count, actual logical edge count, maximum degree, mean degree, and estimated bytes. Extropic describes Z1 as a fixed sparse graph with sixteen neighbors per pbit. A dense logical graph therefore needs a separately verified physical mapping, and passing a degree check alone cannot establish compatibility. The report always states that physical access, embedding, and execution are unverified or unavailable. [Extropic's official Z1 description](https://extropic.ai/writing/from-one-to-one-billion), [official Z1T graph constraints](https://extropic.ai/writing/z1t).

## Verification

`tests/test_encoded_production.py` enumerates interacting three and four category sites with nonuniform fields, restrictions, and a fixed source residue. It compares reduced, QUBO, and Ising energies against independent literal sums over the original parameters. Complete binary transition matrices verify row sums, target stationarity, detailed balance, positive reachability, finite invalid endpoint mass, and the proposal ratio. The enumerator uses exactly the same direction mixture and repeat count as the production backends.

Deliberate incorrect kernels demonstrate failure when the correction is omitted after coefficient changes, or when a forward only sweep is substituted into the reversible shortcut. Forced invalid and valid MH rejection tests verify retained states and exact attempt counts. Repeated seeded runs check reproducibility. Actual JAX and official THRML runs check original protein energies, supported residues, metadata, and the finite schedule.

These checks establish the specified mathematical kernel and its implementation on bounded cases. Accuracy, convergence, runtime advantage, and biological utility remain separate empirical questions. A slow or poorly mixing encoded run must retain that outcome and allow the direct sampler to remain the practical fallback.
