"""Run official THRML block Gibbs on a hand-specified two-spin Ising model."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from examples.week02.pilot_support import exact_ising, frequencies, runtime_record, save_json, sha256


def run_pilot() -> dict:
    configuration = {
        "node_order": ["s0", "s1"],
        "edges": [[0, 1]],
        "biases": [0.25, -0.125],
        "weights": [0.5],
        "beta": 1.0,
        "dtype": "float32",
        "boolean_to_spin": {"false": -1, "true": 1},
        "free_blocks": [[0], [1]],
        "clamped_blocks": [],
        "observed_blocks": [[0, 1]],
        "batch_shape": [],
        "seed": 42,
        "prng_implementation": "threefry2x32",
        "n_warmup": 128,
        "n_samples": 4096,
        "steps_per_sample": 2,
        "requested_platform": "cpu",
    }
    with warnings.catch_warnings(record=True) as observed_warnings:
        warnings.simplefilter("always")
        import jax
        jax.config.update("jax_platforms", "cpu")
        jax.config.update("jax_enable_x64", False)
        import jax.numpy as jnp
        import thrml
        import thrml.block_sampling
        import thrml.models.ising
        from thrml import Block, SamplingSchedule, SpinNode, sample_states
        from thrml.models import IsingEBM, IsingSamplingProgram, hinton_init

        nodes = [SpinNode(), SpinNode()]
        model = IsingEBM(
            nodes,
            [(nodes[0], nodes[1])],
            jnp.asarray(configuration["biases"], dtype=jnp.float32),
            jnp.asarray(configuration["weights"], dtype=jnp.float32),
            jnp.asarray(configuration["beta"], dtype=jnp.float32),
        )
        blocks = [Block([nodes[0]]), Block([nodes[1]])]
        program = IsingSamplingProgram(model, blocks, clamped_blocks=[])
        key = jax.random.key(configuration["seed"], impl=configuration["prng_implementation"])
        init_key, sample_key = jax.random.split(key, 2)
        initial = hinton_init(init_key, model, blocks, ())
        # Synchronize before timing so initialization is excluded.
        jax.block_until_ready(initial)
        schedule = SamplingSchedule(
            n_warmup=configuration["n_warmup"],
            n_samples=configuration["n_samples"],
            steps_per_sample=configuration["steps_per_sample"],
        )
        started = time.perf_counter()
        returned = sample_states(sample_key, program, schedule, initial, [], [Block(nodes)])
        jax.block_until_ready(returned)
        elapsed = time.perf_counter() - started
        full_returned_blocks = [block.tolist() for block in returned]
        shapes = [list(block.shape) for block in returned]
        if shapes != [[configuration["n_samples"], 2]]:
            raise ValueError(f"Unexpected THRML returned shapes: {shapes}.")
        bits = full_returned_blocks[0]
        states = ["".join("1" if value else "0" for value in row) for row in bits]
        source_modules = [thrml, thrml.block_sampling, thrml.models.ising]
        source_hashes = {module.__name__: sha256(Path(module.__file__)) for module in source_modules}
        actual_model = {
            "biases": model.biases.tolist(),
            "weights": model.weights.tolist(),
            "beta": float(model.beta),
            "bias_dtype": str(model.biases.dtype),
            "weight_dtype": str(model.weights.dtype),
            "beta_dtype": str(model.beta.dtype),
        }
        devices = [str(device) for device in jax.devices()]
        keys = {"root": jax.random.key_data(key).tolist(), "initialization": jax.random.key_data(init_key).tolist(), "sampling": jax.random.key_data(sample_key).tolist()}
    exact = exact_ising(tuple(actual_model["biases"]), actual_model["weights"][0], actual_model["beta"])
    return {
        "pilot": "week02_official_thrml_two_spin_ising",
        "scope": "Hand-specified binary systems model; no training or protein encoding.",
        "runtime": runtime_record(),
        "configuration": configuration,
        "actual_model_arrays": actual_model,
        "jax_devices": devices,
        "prng_key_data": keys,
        "initial_free_block_states": [block.tolist() for block in initial],
        "official_api": "thrml.models.IsingEBM; IsingSamplingProgram; hinton_init; thrml.sample_states",
        "official_source_sha256": source_hashes,
        "returned_shapes": shapes,
        "returned_dtypes": [str(block.dtype) for block in returned],
        "returned_blocks": full_returned_blocks,
        "exact_states": exact,
        "comparison": frequencies(states, {row["state"]: row["probability"] for row in exact}),
        "sample_call_elapsed_seconds": elapsed,
        "timing_scope": "Synchronized first CPU sample call after initialization; includes JAX compilation, not a speed benchmark.",
        "warnings": [{"category": type(item.message).__name__, "message": str(item.message)} for item in observed_warnings],
        "limitations": [
            "One trajectory; observations can be correlated, even with two sweeps between samples.",
            "Finite frequency agreement is a pilot check, not proof of convergence or effective sample size.",
            "Samples include the post-warmup observation returned by the official schedule; no rows are removed or deduplicated.",
            "No hardware TSU use, categorical encoding, final sampler or performance-speedup claim.",
        ],
        "documentation": "https://docs.thrml.ai/en/latest/getting-started.html",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "examples/week02/output/thrml.json")
    args = parser.parse_args()
    result = run_pilot()
    save_json(args.output, result)
    comparison = result["comparison"]
    print(f"Official THRML CPU pilot: {comparison['n_returned']} states, repeats retained")
    print(f"Total variation: {comparison['total_variation_distance']:.6f}; maximum probability error: {comparison['max_absolute_error']:.6f}")
    print(f"Sample call: {result['sample_call_elapsed_seconds']:.3f} seconds; warnings: {len(result['warnings'])}")
    print(f"Saved {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
