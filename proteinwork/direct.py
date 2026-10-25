"""Categorical Gibbs sampling with exact support and retained observations."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import itertools
import math
from time import perf_counter

import numpy as np

from .model import ReducedModel, _validateStates


@dataclass
class SamplingResult:
    """Accepted marks a changed observation endpoint, including warmup as its predecessor."""

    states: np.ndarray
    energies: np.ndarray
    initial_states: np.ndarray
    accepted: np.ndarray
    metadata: dict
    timings: dict


def _positiveBeta(beta):
    if isinstance(beta, (bool, np.bool_)) or not isinstance(beta, (int, float, np.integer, np.floating)) or not math.isfinite(beta) or beta <= 0:
        raise ValueError("Beta must be finite and strictly positive")
    return float(beta)


def _integerControl(value, name, minimum):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return int(value)


def _samplingInputs(model, beta, seeds, draws, warmup, sweepsPerDraw, initialStates):
    beta = _positiveBeta(beta)
    draws = _integerControl(draws, "Draw count", 1)
    warmup = _integerControl(warmup, "Warmup sweep count", 0)
    sweepsPerDraw = _integerControl(sweepsPerDraw, "Sweeps per draw", 1)
    try:
        seeds = tuple(seeds)
    except TypeError as error:
        raise ValueError("Seeds must be a nonempty sequence of distinct integers") from error
    if not seeds or any(isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or not 0 <= seed < 2 ** 32 for seed in seeds):
        raise ValueError("Each chain seed must be an integer in the unsigned 32 bit range")
    if len(set(seeds)) != len(seeds):
        raise ValueError("Independent chains require distinct seeds")
    seeds = tuple(int(seed) for seed in seeds)
    if initialStates is None:
        initialStates = np.concatenate([model.random_states(np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(0,))), 1)
                                        for seed in seeds], axis=0)
    initialStates = _validateStates(initialStates, model.n_sites, model.q).copy()
    if initialStates.shape != (len(seeds), model.n_sites) or not np.all(np.isfinite(model.energy(initialStates))):
        raise ValueError("Initial states must contain one supported restricted state per chain")
    return beta, seeds, draws, warmup, sweepsPerDraw, initialStates


def _probabilities(scores, beta):
    scores = np.asarray(scores, dtype=np.float64)
    if not scores.size or not np.all(np.isfinite(scores)):
        raise ValueError("Conditional probabilities require nonempty finite scores")
    with np.errstate(over="ignore", under="ignore"):
        weights = np.exp((scores - np.max(scores, axis=-1, keepdims=True)) * beta)
    return weights / weights.sum(axis=-1, keepdims=True)


def _allStates(model, maxStates):
    maxStates = _integerControl(maxStates, "Enumeration limit", 1)
    count = math.prod(map(len, model.domains))
    if count > maxStates:
        raise ValueError(f"Enumeration requires {count} states, exceeding the limit {maxStates}")
    return np.array(list(itertools.product(*model.domains)), dtype=np.int64).reshape(count, model.n_sites)


def enumerate_target(model: ReducedModel, beta, max_states=65536):
    beta = _positiveBeta(beta)
    states = _allStates(model, max_states)
    energies = model.energy(states)
    return states, _probabilities(-energies, beta)


def _conditionalProbabilities(model, state, site, beta):
    domain = np.asarray(model.domains[site], dtype=np.int64)
    scores = model.fields[site, domain].copy()
    for other, token in enumerate(state):
        scores += model.couplings[site, other, domain, token]
    return _probabilities(scores, beta)


def exact_direct_transition(model: ReducedModel, beta, max_states=256):
    """Return one complete ascending site Gibbs sweep, matching the production schedule."""
    beta = _positiveBeta(beta)
    states = _allStates(model, max_states)
    indices = {tuple(state): index for index, state in enumerate(states)}
    transition = np.eye(len(states), dtype=np.float64)
    for site in range(model.n_sites):
        transition = transition @ _siteTransition(model, states, indices, site, beta)
    return states, transition


def _siteTransition(model, states, indices, site, beta):
    siteKernel = np.zeros((len(states), len(states)), dtype=np.float64)
    for index, state in enumerate(states):
        probabilities = _conditionalProbabilities(model, state, site, beta)
        for token, probability in zip(model.domains[site], probabilities, strict=True):
            target = state.copy()
            target[site] = token
            siteKernel[index, indices[tuple(target)]] = probability
    return siteKernel


def numpy_gibbs_sweep(model: ReducedModel, states, *, beta, rng):
    """Small transparent reference sweep using full conditional recomputation."""
    beta = _positiveBeta(beta)
    stateArray = _validateStates(states, model.n_sites, model.q).copy()
    if stateArray.ndim != 2 or not np.all(np.isfinite(model.energy(stateArray))):
        raise ValueError("Reference sweep requires a supported matrix of chain states")
    for site, domain in enumerate(model.domains):
        for chain, state in enumerate(stateArray):
            probabilities = _conditionalProbabilities(model, state, site, beta)
            stateArray[chain, site] = rng.choice(domain, p=probabilities)
    return stateArray


def _localFields(model, states):
    local = np.broadcast_to(model.fields, (len(states), *model.fields.shape)).copy()
    for site in range(model.n_sites):
        local += np.take(model.couplings[:, site], states[:, site], axis=2).transpose(2, 0, 1)
    return local


@lru_cache(maxsize=32)
def _compiledFunctions(nSites, draws, warmup, sweepsPerDraw):
    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)

    def refresh(states, fields, couplings):
        def oneState(state):
            return fields + jax.lax.fori_loop(0, nSites, lambda site, total: total + couplings[:, site, :, state[site]], jnp.zeros_like(fields))
        return jax.vmap(oneState)(states)

    def sweep(carry, parameters):
        fields, couplings, allowed, beta = parameters
        states, local, keys, count = carry
        local = jax.lax.cond(count % 128 == 0, lambda: refresh(states, fields, couplings), lambda: local)

        def update(site, current):
            currentStates, currentLocal, currentKeys = current
            splitKeys = jax.vmap(lambda key: jax.random.split(key, 2))(currentKeys)
            uniform = jax.vmap(lambda key: jax.random.uniform(key, (), dtype=jnp.float64))(splitKeys[:, 1])
            logits = jnp.where(allowed[site], currentLocal[:, site, :], -jnp.inf)
            weights = jnp.exp((logits - jnp.max(logits, axis=1, keepdims=True)) * beta)
            cumulative = jnp.cumsum(weights / jnp.sum(weights, axis=1, keepdims=True), axis=1)
            # Clamp the tail to one so accumulated rounding cannot choose an invalid last token.
            cumulative = jnp.where(jnp.arange(cumulative.shape[1]) >= jnp.max(jnp.where(allowed[site], jnp.arange(cumulative.shape[1]), 0)), 1.0, cumulative)
            chosen = jnp.sum(uniform[:, None] >= cumulative, axis=1).astype(jnp.int64)
            previous = currentStates[:, site]
            delta = jax.vmap(lambda new, old: couplings[:, site, :, new] - couplings[:, site, :, old])(chosen, previous)
            return currentStates.at[:, site].set(chosen), currentLocal + delta, splitKeys[:, 0]

        states, local, keys = jax.lax.fori_loop(0, nSites, update, (states, local, keys))
        return states, local, keys, count + 1

    def warm(carry, parameters):
        return jax.lax.fori_loop(0, warmup, lambda _, current: sweep(current, parameters), carry)

    def sample(carry, parameters):
        def observe(current, _):
            previous = current[0]
            current = jax.lax.fori_loop(0, sweepsPerDraw, lambda _, value: sweep(value, parameters), current)
            return current, (current[0], jnp.any(current[0] != previous, axis=1))
        return jax.lax.scan(observe, carry, None, length=draws)

    return jax.jit(warm), jax.jit(sample)


def _sampleJax(model, beta, seeds, draws, warmup, sweepsPerDraw, initialStates):
    import jax
    import jax.numpy as jnp

    started = perf_counter()
    warmFunction, sampleFunction = _compiledFunctions(model.n_sites, draws, warmup, sweepsPerDraw)
    allowed = np.zeros_like(model.fields, dtype=bool)
    for site, domain in enumerate(model.domains):
        allowed[site, list(domain)] = True
    keys = jax.vmap(lambda seed: jax.random.key(seed, impl="threefry2x32"))(np.asarray(seeds, dtype=np.uint32))
    carry = (jnp.asarray(initialStates), jnp.asarray(_localFields(model, initialStates)), keys, jnp.asarray(0, dtype=jnp.int64))
    parameters = (jnp.asarray(model.fields), jnp.asarray(model.couplings), jnp.asarray(allowed), jnp.asarray(beta))
    warmExecutable = warmFunction.lower(carry, parameters).compile()
    sampleExecutable = sampleFunction.lower(carry, parameters).compile()
    compileSeconds = perf_counter() - started
    started = perf_counter()
    warmed = warmExecutable(carry, parameters)
    jax.block_until_ready(warmed)
    warmupSeconds = perf_counter() - started
    started = perf_counter()
    _, observations = sampleExecutable(warmed, parameters)
    jax.block_until_ready(observations)
    states = np.asarray(observations[0]).swapaxes(0, 1).copy()
    changed = np.asarray(observations[1]).swapaxes(0, 1).copy()
    samplingSeconds = perf_counter() - started
    return states, changed, {"compile_seconds": compileSeconds, "warmup_seconds": warmupSeconds,
                             "sampling_seconds": samplingSeconds}, {"backend": "jax", "jax_version": jax.__version__,
                             "device": str(jax.devices()[0]), "prng": "threefry2x32", "float64": True,
                             "local_field_refresh_sweeps": 128}


def _sampleNumpy(model, beta, seeds, draws, warmup, sweepsPerDraw, initialStates):
    rngs = [np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(1,))) for seed in seeds]
    current = initialStates.copy()

    def sweep():
        for chain, rng in enumerate(rngs):
            current[chain] = numpy_gibbs_sweep(model, current[chain:chain + 1], beta=beta, rng=rng)[0]

    started = perf_counter()
    for _ in range(warmup):
        sweep()
    warmupSeconds = perf_counter() - started
    states = np.empty((len(seeds), draws, model.n_sites), dtype=np.int64)
    changed = np.empty((len(seeds), draws), dtype=bool)
    started = perf_counter()
    for draw in range(draws):
        previous = current.copy()
        for _ in range(sweepsPerDraw):
            sweep()
        states[:, draw] = current
        changed[:, draw] = np.any(current != previous, axis=1)
    return states, changed, {"compile_seconds": 0.0, "warmup_seconds": warmupSeconds,
                             "sampling_seconds": perf_counter() - started}, {"backend": "numpy", "float64": True,
                             "prng": "PCG64", "fallback_reason": "JAX is not installed"}


def sample_direct(model: ReducedModel, *, beta, seeds, draws, warmup, sweeps_per_draw=1,
                  initial_states=None) -> SamplingResult:
    """Record every requested endpoint, retaining repeated states and independent chain keys."""
    totalStarted = perf_counter()
    beta, seeds, draws, warmup, sweepsPerDraw, initialStates = _samplingInputs(
        model, beta, seeds, draws, warmup, sweeps_per_draw, initial_states)
    if model.n_sites == 0:
        states = np.empty((len(seeds), draws, 0), dtype=np.int64)
        changed = np.zeros((len(seeds), draws), dtype=bool)
        timings = {"compile_seconds": 0.0, "warmup_seconds": 0.0, "sampling_seconds": 0.0}
        metadata = {"backend": "deterministic", "float64": True, "prng": "none"}
    else:
        try:
            import jax
        except ModuleNotFoundError as error:
            if error.name != "jax":
                raise
            states, changed, timings, metadata = _sampleNumpy(model, beta, seeds, draws, warmup, sweepsPerDraw, initialStates)
        else:
            states, changed, timings, metadata = _sampleJax(model, beta, seeds, draws, warmup, sweepsPerDraw, initialStates)
    started = perf_counter()
    energies = model.energy(states)
    timings["energy_seconds"] = perf_counter() - started
    timings["total_seconds"] = perf_counter() - totalStarted
    timings["setup_seconds"] = max(0.0, timings["total_seconds"] - sum(value for key, value in timings.items() if key != "total_seconds"))
    metadata.update({"method": "direct_categorical_gibbs", "schedule": "systematic ascending site order",
                     "observation_unit": "complete sweeps", "sweeps_per_draw": sweepsPerDraw, "warmup_sweeps": warmup,
                     "seeds": list(seeds), "beta": beta, "repeats_retained": True,
                     "accepted_semantics": "endpoint differs from the preceding observation or final warmup state",
                     "initialization": "NumPy PCG64 SeedSequence with spawn key (0,), independent of transition keys"})
    return SamplingResult(states, energies, initialStates, changed, metadata, timings)
