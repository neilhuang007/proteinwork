"""Exact one hot compilation and finite, reversible encoded protein proposals."""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import importlib.metadata
import itertools
import math
from pathlib import Path
import time
from typing import TYPE_CHECKING, Sequence

import numpy as np

from .validation import _integerControl

if TYPE_CHECKING:
    from proteinwork.model import ReducedModel


class EncodingBudgetError(ValueError):
    """The requested encoding exceeds an explicit resource budget."""


def _positiveInteger(value, name, *, allowZero=False):
    return _integerControl(value, name, 0 if allowZero else 1)


def _positiveNumber(value, name, *, allowZero=False):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or result < 0 or (result == 0 and not allowZero):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allowZero else 'positive'}")
    return result


@dataclass(frozen=True)
class EncodedModel:
    """E(z) = constant - biases.s - sum(weights * s[left] * s[right])."""

    biases: np.ndarray
    pairs: np.ndarray
    weights: np.ndarray
    constant: float
    quboLinear: np.ndarray
    quboPairs: np.ndarray
    quboConstant: float
    domains: tuple[tuple[int, ...], ...]
    siteOffsets: np.ndarray
    penalty: float
    estimatedBytes: int
    coefficientScale: float = 1.0
    quantizationStep: float | None = None

    @property
    def nSpins(self):
        return len(self.biases)

    @property
    def bitSites(self):
        return np.repeat(np.arange(len(self.domains)), [len(domain) for domain in self.domains])

    @property
    def bitTokens(self):
        return np.asarray(tuple(itertools.chain.from_iterable(self.domains)), dtype=np.int64)

    def encode(self, states):
        values = np.asarray(states)
        if values.ndim < 1 or values.shape[-1] != len(self.domains) or values.dtype.kind not in "iu":
            raise ValueError("states must contain integer tokens with one column per reduced site")
        for site, domain in enumerate(self.domains):
            if not np.all(np.isin(values[..., site], domain)):
                raise ValueError(f"state contains a forbidden token at reduced site {site}")
        return values[..., self.bitSites] == self.bitTokens

    def decode(self, bits):
        states, isValid = self.decodeChecked(bits)
        if not np.all(isValid):
            raise ValueError("invalid encoded endpoint: every site must contain exactly one active bit")
        return states

    def decodeChecked(self, bits):
        values = _binaryArray(bits, self.nSpins)
        states = np.zeros((*values.shape[:-1], len(self.domains)), dtype=np.int64)
        isValid = np.ones(values.shape[:-1], dtype=bool)
        for site, domain in enumerate(self.domains):
            block = values[..., self.siteOffsets[site]:self.siteOffsets[site + 1]]
            isValid &= block.sum(axis=-1) == 1
            states[..., site] = np.asarray(domain)[block.argmax(axis=-1)]
        return states, isValid

    def energy(self, bits):
        spins = 2.0 * _binaryArray(bits, self.nSpins) - 1.0
        interaction = spins[..., self.pairs[:, 0]] * spins[..., self.pairs[:, 1]]
        return self.constant - spins @ self.biases - interaction @ self.weights

    def interactionMatrix(self):
        matrix = np.zeros((self.nSpins, self.nSpins), dtype=np.float64)
        matrix[self.pairs[:, 0], self.pairs[:, 1]] = self.weights
        matrix[self.pairs[:, 1], self.pairs[:, 0]] = self.weights
        return matrix

    def resource_report(self):
        activePairs = self.pairs[self.weights != 0]
        degrees = np.bincount(activePairs.ravel(), minlength=self.nSpins)
        return {
            "logical_spins": self.nSpins, "logical_edges": len(activePairs),
            "logical_max_degree": int(degrees.max(initial=0)),
            "logical_mean_degree": float(degrees.mean()) if len(degrees) else 0.0,
            "estimated_working_bytes": self.estimatedBytes,
            "estimate_scope": "Conservative array estimate; excludes Python objects and compiler workspace",
            "published_z1_neighbor_limit": 16,
            "exceeds_published_z1_neighbor_limit": bool(np.any(degrees > 16)),
            "physical_access": "unavailable: no physical device connection verified",
            "physical_execution": False, "physical_embedding": "not constructed or verified",
            "hardware_executability": "unverified; logical compilation does not establish physical executability",
            "hardware_source": "https://extropic.ai/writing/from-one-to-one-billion",
        }


def _binaryArray(bits, size):
    values = np.asarray(bits)
    if values.ndim < 1 or values.shape[-1] != size or not np.all((values == 0) | (values == 1)):
        raise ValueError(f"bits must be binary with final dimension {size}")
    return values.astype(bool, copy=False)


def _validateModel(model):
    domains = tuple(tuple(int(token) for token in domain) for domain in model.domains)
    if len(domains) != model.n_sites or any(len(domain) < 2 for domain in domains):
        raise ValueError("compile_ising requires a reduced model with singleton sites already removed")
    for domain in domains:
        if len(set(domain)) != len(domain) or min(domain) < 0 or max(domain) >= model.q:
            raise ValueError("domains must contain unique, supported alphabet indices")
    hasFiniteFields = all(np.isfinite(model.fields[site, domain]).all() for site, domain in enumerate(domains))
    if not hasFiniteFields or not np.isfinite(model.couplings).all() or not math.isfinite(model.offset):
        raise ValueError("all supported reduced model coefficients must be finite")
    return domains


def _edgeCount(model, domains, penalty):
    count = sum(len(domain) * (len(domain) - 1) // 2 for domain in domains) if penalty else 0
    for left, right in itertools.combinations(range(len(domains)), 2):
        count += np.count_nonzero(model.couplings[left, right][np.ix_(domains[left], domains[right])])
    return int(count)


def _checkBudget(spins, edges, maxSpins, maxEdges, maxBytes):
    limits = [_positiveInteger(value, name, allowZero=True) for value, name in (
        (maxSpins, "max_spins"), (maxEdges, "max_edges"), (maxBytes, "max_bytes"))]
    estimatedBytes = 3 * spins * spins * 8 + edges * 80 + spins * 256
    requested = (spins, edges, estimatedBytes)
    for name, count, limit in zip(("spins", "edges", "estimated bytes"), requested, limits):
        if count > limit:
            raise EncodingBudgetError(f"encoding needs {count} {name}, exceeding budget {limit}; no coefficients were discarded")
    return estimatedBytes


def _quboArrays(model, domains, offsets, penalty, edgeCount):
    linear = np.concatenate([-np.asarray(model.fields[site, domain]) - penalty
                             for site, domain in enumerate(domains)]) if domains else np.empty(0)
    pairs = np.empty((edgeCount, 2), dtype=np.int64)
    values = np.empty(edgeCount, dtype=np.float64)
    cursor = 0
    for site, domain in enumerate(domains):
        if penalty == 0:
            continue
        localPairs = np.asarray(list(itertools.combinations(range(len(domain)), 2)), dtype=np.int64)
        length = len(localPairs)
        pairs[cursor:cursor + length] = localPairs + offsets[site]
        values[cursor:cursor + length] = 2 * penalty
        cursor += length
    for left, right in itertools.combinations(range(len(domains)), 2):
        block = -model.couplings[left, right][np.ix_(domains[left], domains[right])]
        rows, columns = np.nonzero(block)
        length = len(rows)
        pairs[cursor:cursor + length] = np.column_stack((rows + offsets[left], columns + offsets[right]))
        values[cursor:cursor + length] = block[rows, columns]
        cursor += length
    return linear, pairs, values


def compile_ising(model, penalty=None, *, max_spins=2048, max_edges=2000000, max_bytes=536870912):
    """Compile every supported target coefficient without sparsification or clipping."""
    domains = _validateModel(model)
    penaltyValue = _positiveNumber(2.0 if penalty is None else penalty, "penalty", allowZero=True)
    offsets = np.concatenate(([0], np.cumsum([len(domain) for domain in domains]))).astype(np.int64)
    edgeCount = _edgeCount(model, domains, penaltyValue)
    estimatedBytes = _checkBudget(int(offsets[-1]), edgeCount, max_spins, max_edges, max_bytes)
    linear, pairs, values = _quboArrays(model, domains, offsets, penaltyValue, edgeCount)
    biases = -linear / 2
    np.add.at(biases, pairs[:, 0], -values / 4)
    np.add.at(biases, pairs[:, 1], -values / 4)
    quboConstant = float(model.offset + penaltyValue * len(domains))
    constant = float(quboConstant + linear.sum() / 2 + values.sum() / 4)
    if not all(np.isfinite(array).all() for array in (biases, values, linear)) or not math.isfinite(constant):
        raise ValueError("encoding arithmetic overflowed; rescale the model explicitly")
    return EncodedModel(biases, pairs, -values / 4, constant, linear, values,
                        quboConstant, domains, offsets, penaltyValue, estimatedBytes)


def _implementedEncoding(encoded, coefficientScale, quantizationStep):
    scale = _positiveNumber(coefficientScale, "coefficient_scale")
    step = None if quantizationStep is None else _positiveNumber(quantizationStep, "quantization_step")
    biases, weights = encoded.biases * scale, encoded.weights * scale
    if step is not None:
        biases, weights = np.rint(biases / step) * step, np.rint(weights / step) * step
    if not np.isfinite(biases).all() or not np.isfinite(weights).all() or not math.isfinite(encoded.constant * scale):
        raise ValueError("implemented coefficient arithmetic overflowed")
    return replace(encoded, biases=biases, weights=weights, constant=encoded.constant * scale,
                   coefficientScale=scale, quantizationStep=step)


def _sigmoid(logits):
    magnitude = np.exp(-np.abs(logits))
    return np.where(np.asarray(logits) >= 0, 1 / (1 + magnitude), magnitude / (1 + magnitude))


def _numpyEndpoint(bits, rng, encoded, matrix, proposalBeta, sweeps):
    spins = 2.0 * bits - 1.0
    fields = encoded.biases + matrix @ spins
    order = range(encoded.nSpins) if rng.integers(2) else range(encoded.nSpins - 1, -1, -1)
    for _ in range(sweeps):
        for bit in order:
            probability = float(_sigmoid(2 * proposalBeta * fields[bit]))
            if not 0 < probability < 1:
                raise FloatingPointError("binary conditional saturated; reduce proposal_beta or coefficient_scale")
            newSpin = 1.0 if rng.random() < probability else -1.0
            fields += matrix[:, bit] * (newSpin - spins[bit])
            spins[bit] = newSpin
    return spins > 0


def _correctEndpoint(current, currentEnergy, endpoint, rng, model, encoded, beta, proposalBeta):
    proposal, isValid = encoded.decodeChecked(endpoint)
    if not bool(isValid):
        return current, currentEnergy, False, True
    proposalEnergy = float(model.energy(proposal))
    logRatio = -beta * (proposalEnergy - currentEnergy)
    logRatio += proposalBeta * float(encoded.energy(endpoint) - encoded.energy(encoded.encode(current)))
    if not math.isfinite(proposalEnergy) or not math.isfinite(currentEnergy) or not math.isfinite(logRatio):
        raise FloatingPointError("MH energy arithmetic is nonfinite; use a numerically representable model and temperature")
    isAccepted = math.log(max(float(rng.random()), np.finfo(float).tiny)) < min(0.0, logRatio)
    if isAccepted:
        return proposal, proposalEnergy, True, False
    return current, currentEnergy, False, False


def _initialStates(model, encoded, seeds, initialStates):
    if initialStates is not None:
        states = np.asarray(initialStates)
        if states.shape != (len(seeds), model.n_sites):
            raise ValueError("initial_states shape must be (chains, reduced sites)")
        encoded.encode(states)
        return states.astype(np.int64, copy=True)
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        rows.append(np.asarray([rng.choice(domain) for domain in model.domains], dtype=np.int64))
    return np.stack(rows)


def _numpySample(model, encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, initialStates, endpointFunction=None):
    states = np.empty((len(seeds), draws + warmup, model.n_sites), dtype=np.int64)
    energies = np.empty((len(seeds), draws + warmup), dtype=np.float64)
    accepted = np.zeros((len(seeds), draws + warmup), dtype=bool)
    invalid = np.zeros_like(accepted)
    matrix = encoded.interactionMatrix()
    for chain, seed in enumerate(seeds):
        rng = np.random.default_rng(np.random.SeedSequence([int(seed), 1]))
        current = initialStates[chain].copy()
        currentEnergy = float(model.energy(current))
        for attempt in range(draws + warmup):
            bits = encoded.encode(current)
            endpoint = (_numpyEndpoint(bits, rng, encoded, matrix, proposalBeta, sweeps)
                        if endpointFunction is None else endpointFunction(bits, rng, chain, attempt))
            current, currentEnergy, accepted[chain, attempt], invalid[chain, attempt] = _correctEndpoint(
                current, currentEnergy, endpoint, rng, model, encoded, beta, proposalBeta)
            states[chain, attempt], energies[chain, attempt] = current, currentEnergy
    return states, energies, accepted, invalid


def _jaxImports():
    try:
        import jax
        jax.config.update("jax_enable_x64", True)
        import jax.numpy as jnp
    except ImportError as error:
        raise ImportError("JAX backend requires the pinned pilot environment") from error
    return jax, jnp


def _jaxEnergyFunctions(model, encoded, jnp):
    fields, couplings = jnp.asarray(model.fields), jnp.asarray(model.couplings)
    sites = jnp.arange(model.n_sites)
    left, right = np.triu_indices(model.n_sites, k=1)
    def originalEnergy(state):
        return model.offset - jnp.sum(fields[sites, state]) - jnp.sum(couplings[left, right, state[left], state[right]])
    biases, weights = jnp.asarray(encoded.biases), jnp.asarray(encoded.weights)
    pairLeft, pairRight = encoded.pairs.T
    def surrogateEnergy(bits):
        spins = bits.astype(jnp.float64) * 2 - 1
        return encoded.constant - spins @ biases - jnp.sum(weights * spins[pairLeft] * spins[pairRight])
    return originalEnergy, surrogateEnergy


def _jaxEndpointFunction(encoded, proposalBeta, sweeps, jax, jnp):
    matrix, biases = jnp.asarray(encoded.interactionMatrix()), jnp.asarray(encoded.biases)
    def endpoint(bits, key):
        directionKey, updateKey = jax.random.split(key)
        isForward = jax.random.bernoulli(directionKey)
        spins = bits.astype(jnp.float64) * 2 - 1
        fields = biases + matrix @ spins
        def update(carry, step):
            spins, fields, key, isSaturated = carry
            key, subkey = jax.random.split(key)
            index = step % encoded.nSpins
            bit = jnp.where(isForward, index, encoded.nSpins - 1 - index)
            probability = jax.nn.sigmoid(2 * proposalBeta * fields[bit])
            newSpin = jax.random.bernoulli(subkey, probability).astype(jnp.float64) * 2 - 1
            fields = fields + matrix[:, bit] * (newSpin - spins[bit])
            isSaturated |= ~jnp.isfinite(probability) | (probability <= 0) | (probability >= 1)
            return (spins.at[bit].set(newSpin), fields, key, isSaturated), None
        result, _ = jax.lax.scan(update, (spins, fields, updateKey, jnp.array(False)), jnp.arange(sweeps * encoded.nSpins))
        return result[0] > 0, result[3]
    return endpoint


def _jaxChainFunction(model, encoded, beta, proposalBeta, attempts, endpoint, jax, jnp):
    bitSites, bitTokens = jnp.asarray(encoded.bitSites), jnp.asarray(encoded.bitTokens)
    originalEnergy, surrogateEnergy = _jaxEnergyFunctions(model, encoded, jnp)
    def chain(key, initialState):
        def attempt(carry, unused):
            current, currentEnergy, key = carry
            key, endpointKey, acceptKey = jax.random.split(key, 3)
            bits = current[bitSites] == bitTokens
            proposedBits, isSaturated = endpoint(bits, endpointKey)
            counts = jax.ops.segment_sum(proposedBits.astype(jnp.int64), bitSites, model.n_sites)
            isValid = jnp.all(counts == 1)
            rawProposal = jax.ops.segment_sum(proposedBits * bitTokens, bitSites, model.n_sites)
            proposal = jnp.where(isValid, rawProposal, current)
            proposalEnergy = originalEnergy(proposal)
            logRatio = -beta * (proposalEnergy - currentEnergy) + proposalBeta * (surrogateEnergy(proposedBits) - surrogateEnergy(bits))
            isSaturated |= ~jnp.isfinite(logRatio) | ~jnp.isfinite(currentEnergy) | ~jnp.isfinite(proposalEnergy)
            logUniform = jnp.log(jax.random.uniform(acceptKey, dtype=jnp.float64, minval=jnp.finfo(jnp.float64).tiny))
            isAccepted = isValid & (logUniform < jnp.minimum(0.0, logRatio))
            current = jnp.where(isAccepted, proposal, current)
            currentEnergy = jnp.where(isAccepted, proposalEnergy, currentEnergy)
            return (current, currentEnergy, key), (current, currentEnergy, isAccepted, ~isValid, isSaturated)
        return jax.lax.scan(attempt, (initialState, originalEnergy(initialState), key), None, length=attempts)[1]
    return chain


def _jaxSample(model, encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, initialStates):
    jax, jnp = _jaxImports()
    endpoint = _jaxEndpointFunction(encoded, proposalBeta, sweeps, jax, jnp)
    chain = _jaxChainFunction(model, encoded, beta, proposalBeta, draws + warmup, endpoint, jax, jnp)
    keys = jnp.stack([jax.random.fold_in(jax.random.key(int(seed), impl="threefry2x32"), 1) for seed in seeds])
    runner = jax.jit(jax.vmap(chain))
    started = time.perf_counter()
    executable = runner.lower(keys, jnp.asarray(initialStates)).compile()
    compileSeconds = time.perf_counter() - started
    started = time.perf_counter()
    returned = executable(keys, jnp.asarray(initialStates))
    jax.block_until_ready(returned)
    sampleSeconds = time.perf_counter() - started
    arrays = tuple(np.asarray(value) for value in returned)
    if arrays[-1].any():
        raise FloatingPointError("binary conditional or MH arithmetic is nonfinite or saturated; use representable target energies and reduce proposal_beta or coefficient_scale")
    return arrays[:4], {"compilation_seconds": compileSeconds, "sampling_seconds": sampleSeconds}, {
        "jax_devices": [str(device) for device in jax.devices()], "jax_version": jax.__version__,
        "jax_enable_x64": bool(jax.config.jax_enable_x64), "rng": "JAX threefry2x32; independent chain keys"}


def _thrmlEndpointFunction(encoded, proposalBeta, sweeps, seeds):
    jax, jnp = _jaxImports()
    try:
        import thrml
        import thrml.block_sampling
        import thrml.models.ising
        from thrml import Block, SamplingSchedule, SpinNode, sample_states
        from thrml.models import IsingEBM, IsingSamplingProgram
    except ImportError as error:
        raise ImportError("THRML backend requires the official thrml package in the pinned pilot environment") from error
    maximumLogit = 2 * proposalBeta * (np.abs(encoded.biases) + np.abs(encoded.interactionMatrix()).sum(axis=1))
    if np.any(maximumLogit >= 36.0):
        raise FloatingPointError("THRML conditional bounds risk binary64 saturation; reduce proposal_beta or coefficient_scale")
    nodes = [SpinNode() for _ in range(encoded.nSpins)]
    edges = [(nodes[left], nodes[right]) for left, right in encoded.pairs]
    weights = encoded.weights
    hasPaddingEdge = not edges
    if hasPaddingEdge:
        edges, weights = [(nodes[0], nodes[1])], np.zeros(1)
    ebm = IsingEBM(nodes, edges, jnp.asarray(encoded.biases), jnp.asarray(weights), jnp.asarray(proposalBeta))
    schedule = SamplingSchedule(n_warmup=sweeps, n_samples=1, steps_per_sample=1)
    forward = list(range(encoded.nSpins))
    runners = [_compileThrmlOrder(order, nodes, ebm, schedule, jax, jnp, Block, IsingSamplingProgram, sample_states)
               for order in (forward, forward[::-1])]
    keys = [jax.random.fold_in(jax.random.key(int(seed), impl="threefry2x32"), 1) for seed in seeds]
    def endpoint(bits, rng, chain, attempt):
        direction = int(rng.integers(2))
        key = jax.random.fold_in(keys[chain], attempt)
        result = runners[direction](key, jnp.asarray(bits))
        jax.block_until_ready(result)
        return np.asarray(result, dtype=bool)
    modules = (thrml, thrml.block_sampling, thrml.models.ising)
    metadata = {"official_api": "thrml.models.IsingEBM; IsingSamplingProgram; thrml.sample_states",
                "thrml_version": importlib.metadata.version("thrml"), "jax_version": jax.__version__,
                "jax_devices": [str(device) for device in jax.devices()], "jax_enable_x64": bool(jax.config.jax_enable_x64),
                "zero_weight_api_padding_edge": hasPaddingEdge,
                "rng": "NumPy PCG64 for direction and MH; JAX threefry2x32 keys folded by chain seed and attempt",
                "thrml_schedule": {"n_warmup": sweeps, "n_samples": 1, "steps_per_sample": 1},
                "source_sha256": {module.__name__: hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() for module in modules}}
    return endpoint, metadata


def _compileThrmlOrder(order, nodes, ebm, schedule, jax, jnp, Block, IsingSamplingProgram, sampleStates):
    blocks = [Block([nodes[index]]) for index in order]
    program = IsingSamplingProgram(ebm, blocks, clamped_blocks=[])
    if program.gibbs_spec.sampling_order != [[index] for index in range(len(order))]:
        raise RuntimeError("official THRML singleton schedule differs from the audited sequential order")
    def run(key, bits):
        initial = [bits[index:index + 1] for index in order]
        return sampleStates(key, program, schedule, initial, [], [Block(nodes)])[0][0]
    return jax.jit(run).lower(jax.random.key(0), jnp.zeros(len(nodes), dtype=bool)).compile()


def _validateSampling(beta, seeds, draws, warmup, sweeps, backend):
    betaValue = _positiveNumber(beta, "beta")
    if backend not in ("numpy", "jax", "thrml"):
        raise ValueError("backend must be numpy, jax, or thrml")
    try:
        seedValues = tuple(_positiveInteger(seed, "seed", allowZero=True) for seed in seeds)
    except TypeError as error:
        raise ValueError("seeds must be a nonempty sequence of distinct integers") from error
    if not seedValues or len(set(seedValues)) != len(seedValues) or any(seed >= 2**32 for seed in seedValues):
        raise ValueError("seeds must be distinct integers in [0, 2**32)")
    return betaValue, seedValues, _positiveInteger(draws, "draws"), _positiveInteger(warmup, "warmup", allowZero=True), _positiveInteger(sweeps, "sweeps_per_draw")


def sample_encoded(model, *, beta, seeds: Sequence[int], draws, warmup, sweeps_per_draw=1,
                   initial_states=None, backend="numpy", penalty=None, max_spins=2048,
                   max_edges=2000000, proposal_beta=None, coefficient_scale=1.0,
                   quantization_step=None, max_bytes=536870912):
    """Return every outer attempt, retaining invalid endpoint and MH rejections."""
    from proteinwork.direct import SamplingResult
    started = time.perf_counter()
    beta, seeds, draws, warmup, sweeps = _validateSampling(beta, seeds, draws, warmup, sweeps_per_draw, backend)
    proposalBeta = beta if proposal_beta is None else _positiveNumber(proposal_beta, "proposal_beta")
    encoded = _implementedEncoding(compile_ising(model, penalty, max_spins=max_spins, max_edges=max_edges, max_bytes=max_bytes), coefficient_scale, quantization_step)
    if backend == "thrml" and encoded.nSpins > 64:
        raise EncodingBudgetError("official THRML singleton compilation is limited to 64 spins by this adapter; use the named JAX simulator for larger logical graphs")
    outputBytes = len(seeds) * (draws + warmup) * (8 * model.n_sites + 11)
    outputBytes += len(seeds) * warmup * (40 * model.n_sites + 128)
    if outputBytes + encoded.estimatedBytes > max_bytes:
        raise EncodingBudgetError(f"encoding and complete trace need at least {outputBytes + encoded.estimatedBytes} bytes, exceeding {max_bytes}")
    initialStates = _initialStates(model, encoded, seeds, initial_states)
    timings = {"encoding_initialization_seconds": time.perf_counter() - started, "compilation_seconds": 0.0}
    metadata = _samplingMetadata(encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, backend)
    arrays, backendTimings, backendMetadata = _runBackend(model, encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, initialStates, backend)
    timings.update(backendTimings)
    metadata.update(backendMetadata)
    states, energies, accepted, invalid = arrays
    metadata.update({"invalid_endpoint_count": int(invalid[:, warmup:].sum()), "mh_rejection_count": int((~accepted[:, warmup:] & ~invalid[:, warmup:]).sum()),
                     "warmup_states": states[:, :warmup].tolist(), "warmup_energies": energies[:, :warmup].tolist(),
                     "warmup_accepted": accepted[:, :warmup].tolist(), "invalid_endpoints": invalid[:, warmup:].tolist(),
                     "warmup_invalid_endpoints": invalid[:, :warmup].tolist()})
    timings["total_seconds"] = time.perf_counter() - started
    return SamplingResult(states=states[:, warmup:], energies=energies[:, warmup:], initial_states=initialStates,
                          accepted=accepted[:, warmup:], metadata=metadata, timings=timings)


def _samplingMetadata(encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, backend):
    return {"method": "corrected_encoded", "backend": backend, "encoding": "one_hot_potts_qubo_ising",
            "beta": beta, "proposal_beta": proposalBeta, "penalty": encoded.penalty,
            "coefficient_scale": encoded.coefficientScale, "quantization_step": encoded.quantizationStep,
            "precision": "float64", "seeds": list(seeds), "draws": draws, "warmup": warmup,
            "sweeps_per_draw": sweeps, "binary_updates_per_attempt": encoded.nSpins * sweeps,
            "proposal_schedule": "fair mixture of forward singleton sweeps and reverse singleton sweeps; one direction per attempt",
            "proposal_equilibrium_assumed": False, "invalid_endpoint_rule": "reject once and retain current state",
            "mh_log_ratio": "-beta*(E(y)-E(x)) + proposal_beta*(Es(code(y))-Es(code(x)))",
            "domains": [list(domain) for domain in encoded.domains], "resources": encoded.resource_report(),
            "rng": "NumPy PCG64; separate initialization and transition streams", "physical_execution": False,
            "accepted_semantics": "MH accepted, including valid self proposals; movement is a separate observable"}


def _runBackend(model, encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, initialStates, backend):
    if encoded.nSpins == 0:
        states = np.broadcast_to(initialStates[:, None, :], (len(seeds), draws + warmup, 0)).copy()
        energies = np.full((len(seeds), draws + warmup), float(model.offset))
        accepted = np.ones(energies.shape, dtype=bool)
        return (states, energies, accepted, ~accepted), {"sampling_seconds": 0.0}, {"backend_execution": "zero variable target; no sampler call"}
    if backend == "jax":
        return _jaxSample(model, encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, initialStates)
    endpointFunction, backendMetadata, compileSeconds = None, {}, 0.0
    if backend == "thrml":
        compileStarted = time.perf_counter()
        endpointFunction, backendMetadata = _thrmlEndpointFunction(encoded, proposalBeta, sweeps, seeds)
        compileSeconds = time.perf_counter() - compileStarted
    sampleStarted = time.perf_counter()
    arrays = _numpySample(model, encoded, beta, proposalBeta, seeds, draws, warmup, sweeps, initialStates, endpointFunction)
    return arrays, {"sampling_seconds": time.perf_counter() - sampleStarted, "compilation_seconds": compileSeconds}, backendMetadata


def _binarySweepMatrix(encoded, proposalBeta, words, order):
    count = len(words)
    transition = np.eye(count)
    matrix = encoded.interactionMatrix()
    spins = words.astype(np.float64) * 2 - 1
    indexes = np.arange(count)
    for bit in order:
        probability = _sigmoid(2 * proposalBeta * (encoded.biases[bit] + spins @ matrix[:, bit]))
        if np.any((probability <= 0) | (probability >= 1)):
            raise FloatingPointError("exact binary conditional saturated")
        flipped = indexes ^ (1 << (encoded.nSpins - bit - 1))
        stay = np.where(words[:, bit], probability, 1 - probability)
        updated = transition * stay[None, :]
        updated[:, flipped] += transition * (1 - stay)[None, :]
        transition = updated
    return transition


def exact_encoded_transition(model, beta, penalty=None, max_states=256, max_binary_states=4096, *,
                             sweeps_per_draw=1, proposal_beta=None, coefficient_scale=1.0, quantization_step=None):
    """Enumerate the same finite sweep mixture and correction used by every backend."""
    beta = _positiveNumber(beta, "beta")
    max_states = _positiveInteger(max_states, "max_states")
    max_binary_states = _positiveInteger(max_binary_states, "max_binary_states")
    proposalBeta = beta if proposal_beta is None else _positiveNumber(proposal_beta, "proposal_beta")
    sweeps = _positiveInteger(sweeps_per_draw, "sweeps_per_draw")
    encoded = _implementedEncoding(compile_ising(model, penalty), coefficient_scale, quantization_step)
    stateCount = math.prod(map(len, model.domains))
    binaryCount = 1 << encoded.nSpins
    if stateCount > max_states or binaryCount > max_binary_states:
        raise EncodingBudgetError(f"exact transition needs {stateCount} target states and {binaryCount} binary states")
    if binaryCount * binaryCount * 8 * 5 > 536870912:
        raise EncodingBudgetError("exact transition matrices exceed 512 MiB working budget")
    states = np.asarray(list(itertools.product(*model.domains)), dtype=np.int64).reshape(stateCount, model.n_sites)
    words = np.asarray(list(itertools.product((0, 1), repeat=encoded.nSpins)), dtype=bool).reshape(binaryCount, encoded.nSpins)
    forward = _binarySweepMatrix(encoded, proposalBeta, words, range(encoded.nSpins))
    reverse = _binarySweepMatrix(encoded, proposalBeta, words, range(encoded.nSpins - 1, -1, -1))
    proposal = (np.linalg.matrix_power(forward, sweeps) + np.linalg.matrix_power(reverse, sweeps)) / 2
    codes = encoded.encode(states)
    indexes = codes.astype(np.int64) @ (2 ** np.arange(encoded.nSpins - 1, -1, -1))
    validProposal = proposal[np.ix_(indexes, indexes)]
    energies, surrogate = np.asarray(model.energy(states)), encoded.energy(codes)
    logRatio = -beta * (energies[None, :] - energies[:, None]) + proposalBeta * (surrogate[None, :] - surrogate[:, None])
    transition = validProposal * np.exp(np.minimum(0.0, logRatio))
    transition[np.diag_indices(stateCount)] += 1 - transition.sum(axis=1)
    return states, transition
