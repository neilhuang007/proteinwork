"""Traceable assay ranking and constrained selection of an experimental library."""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_array
from sklearn.linear_model import BayesianRidge


MAX_POOL = 128
MAX_ASSAY_ROWS = 10_000
MAX_FEATURES = 4096
MAX_TRAINING_WORKING_BYTES = 256 * 1024 * 1024
ASSAY_COLUMNS = {"aligned_sequence", "value", "split", "assay_id", "units", "context", "round"}
SETTING_KEYS = {
    "batch_size", "min_hamming", "max_hamming", "min_pairwise_hamming",
    "forbidden_motifs", "required_motifs", "diversity_weight", "objective", "direction",
    "exploration_weight", "optimizer", "max_pool", "max_seconds", "assay_id", "units",
    "context", "allowed_residues", "residue_class", "min_residue_count", "max_residue_count",
    "optimizer_seed", "optimizer_draws", "optimizer_warmup", "optimizer_beta",
}


def _integer(value, name, lower, upper):
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError(f"{name} must be an integer from {lower} through {upper}")
    return value


def _number(value, name, lower=0.0, upper=1e9):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or not lower <= result <= upper:
        raise ValueError(f"{name} must be finite and between {lower} and {upper}")
    return result


def _text(value, name):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{name} must be nonempty text without surrounding whitespace")
    return value


def _configuration(settings, anchor, alphabet, has_assays):
    if not isinstance(settings, dict) or set(settings) - SETTING_KEYS:
        unknown = sorted(set(settings) - SETTING_KEYS) if isinstance(settings, dict) else []
        raise ValueError(f"Unknown design settings: {unknown}")
    if not isinstance(alphabet, str) or not alphabet or len(set(alphabet)) != len(alphabet):
        raise ValueError("Alphabet must contain distinct symbols")
    if any(symbol.isspace() for symbol in alphabet):
        raise ValueError("Alphabet cannot contain whitespace")
    if not isinstance(anchor, str) or not anchor or set(anchor) - set(alphabet):
        raise ValueError("Anchor must be a nonempty aligned sequence inside the alphabet")
    length = len(anchor)
    config = dict(settings)
    default_pool = 64 if config.get("optimizer") == "thrml" else MAX_POOL
    config["max_pool"] = _integer(config.get("max_pool", default_pool), "max_pool", 1, MAX_POOL)
    config["batch_size"] = _integer(config.get("batch_size", 8), "batch_size", 0, config["max_pool"])
    for name, default in (("min_hamming", 0), ("max_hamming", length), ("min_pairwise_hamming", 0)):
        config[name] = _integer(config.get(name, default), name, 0, length)
    if config["min_hamming"] > config["max_hamming"]:
        raise ValueError("Minimum Hamming distance exceeds maximum")
    for name, default in (("diversity_weight", 0.0), ("exploration_weight", 1.0)):
        config[name] = _number(config.get(name, default), name)
    config["max_seconds"] = _number(config.get("max_seconds", 30.0), "max_seconds", 0.001, 3600.0)
    for name, default, values in (
        ("objective", "model_score", {"model_score", "assay_ucb"}),
        ("direction", "maximize", {"maximize", "minimize"}),
        ("optimizer", "milp", {"milp", "greedy", "thrml"}),
    ):
        config[name] = config.get(name, default)
        if not isinstance(config[name], str) or config[name] not in values:
            raise ValueError(f"Unknown {name}: {config[name]}")
    if config["optimizer"] == "thrml" and config["max_pool"] > 64:
        raise ValueError("THRML batch optimization supports at most 64 candidate variables")
    config["optimizer_seed"] = _integer(config.get("optimizer_seed", 1729), "optimizer_seed", 0, 2**32 - 1)
    config["optimizer_draws"] = _integer(config.get("optimizer_draws", 512), "optimizer_draws", 1, 8192)
    config["optimizer_warmup"] = _integer(config.get("optimizer_warmup", 128), "optimizer_warmup", 0, 8192)
    config["optimizer_beta"] = _number(config.get("optimizer_beta", 4.0), "optimizer_beta", 0.001, 15.0)
    if config["objective"] == "assay_ucb" and not has_assays:
        raise ValueError("assay_ucb requires an assay CSV")
    if config["objective"] == "model_score" and config["direction"] != "maximize":
        raise ValueError("model_score always maximizes the negative Potts energy proxy")
    for name in ("required_motifs", "forbidden_motifs"):
        motifs = config.get(name, [])
        if not isinstance(motifs, list) or any(not isinstance(m, str) or not m or set(m) - (set(alphabet) - {"-"}) for m in motifs):
            raise ValueError(f"{name} must be literal nonempty ungapped strings inside the alphabet")
        if len(set(motifs)) != len(motifs):
            raise ValueError(f"{name} contains duplicates")
        config[name] = list(motifs)
    residue_class = config.get("residue_class", "")
    if not isinstance(residue_class, str) or len(set(residue_class)) != len(residue_class) or set(residue_class) - (set(alphabet) - {"-"}):
        raise ValueError("residue_class must contain distinct ungapped alphabet symbols")
    if not residue_class and (config.get("min_residue_count", 0) != 0 or config.get("max_residue_count", length) != length):
        raise ValueError("Residue count limits require residue_class")
    config["residue_class"] = residue_class
    config["min_residue_count"] = _integer(config.get("min_residue_count", 0), "min_residue_count", 0, length)
    config["max_residue_count"] = _integer(config.get("max_residue_count", length), "max_residue_count", 0, length)
    if config["min_residue_count"] > config["max_residue_count"]:
        raise ValueError("Minimum residue count exceeds maximum")
    domains = config.get("allowed_residues")
    if domains is None and has_assays:
        raise ValueError("Assay validation requires allowed_residues from original model support")
    if domains is not None:
        if not isinstance(domains, (list, tuple)) or len(domains) != length:
            raise ValueError("allowed_residues needs one symbol string per aligned position")
        if any(not isinstance(d, str) or not d or len(set(d)) != len(d) or set(d) - set(alphabet) for d in domains):
            raise ValueError("allowed_residues contains an empty, duplicated or unknown domain")
        config["allowed_residues"] = list(domains)
    if has_assays:
        for name in ("assay_id", "units", "context"):
            config[name] = _text(config.get(name), name)
    return config


def _validate_sequence(sequence, anchor, alphabet, domains, label):
    if not isinstance(sequence, str) or len(sequence) != len(anchor) or set(sequence) - set(alphabet):
        raise ValueError(f"{label} has incorrect aligned length or unknown residues")
    if domains is not None and any(symbol not in domain for symbol, domain in zip(sequence, domains)):
        raise ValueError(f"{label} is outside the original model support")


def _read_assays(path, anchor, alphabet, config):
    path = Path(path)
    if path.stat().st_size > 20_000_000:
        raise ValueError("Assay CSV exceeds the 20 MB input limit")
    data = path.read_bytes()
    if len(data) > 20_000_000:
        raise ValueError("Assay CSV exceeds the 20 MB input limit")
    import io
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    if reader.fieldnames is None or len(reader.fieldnames) != len(ASSAY_COLUMNS) or set(reader.fieldnames) != ASSAY_COLUMNS:
        raise ValueError(f"Assay CSV columns must be exactly {sorted(ASSAY_COLUMNS)}")
    grouped = {}
    rows, outside_support_rows = 0, 0
    for line, row in enumerate(reader, 2):
        rows += 1
        if rows > MAX_ASSAY_ROWS or None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed or oversized assay CSV at row {line}")
        sequence = row["aligned_sequence"]
        _validate_sequence(sequence, anchor, alphabet, None, f"Assay row {line}")
        if any((symbol == "-") != (anchor_symbol == "-") for symbol, anchor_symbol in zip(sequence, anchor)):
            raise ValueError(f"Assay row {line} has a different construct gap pattern from the coordinate anchor")
        outside_support = any(symbol not in domain for symbol, domain in zip(sequence, config["allowed_residues"]))
        outside_support_rows += int(outside_support)
        split = row["split"]
        if split not in {"train", "validation", "test"}:
            raise ValueError(f"Unknown assay split at row {line}")
        for field in ("assay_id", "units", "context"):
            if row[field] != config[field]:
                raise ValueError(f"Assay {field} mismatch at row {line}")
        try:
            value = float(row["value"])
            round_number = int(row["round"])
        except (ValueError, TypeError) as error:
            raise ValueError(f"Invalid assay value or round at row {line}") from error
        if not math.isfinite(value) or abs(value) > 1e100 or round_number < 0:
            raise ValueError(f"Assay value must be finite and round nonnegative at row {line}")
        if sequence in grouped and grouped[sequence]["split"] != split:
            raise ValueError("Repeated assay sequences cannot cross train, validation and test splits")
        group = grouped.setdefault(sequence, {"aligned_sequence": sequence, "split": split, "values": [], "rounds": [],
                                              "outside_model_support": outside_support})
        group["values"].append(value)
        group["rounds"].append(round_number)
    if not grouped:
        raise ValueError("Assay CSV contains no observations")
    observations = []
    for sequence in sorted(grouped):
        group = grouped[sequence]
        values = group["values"]
        observations.append({"aligned_sequence": sequence, "split": group["split"],
                             "value": math.fsum(v / len(values) for v in values),
                             "replicates": len(values), "rounds": sorted(set(group["rounds"])),
                             "outside_model_support": group["outside_model_support"]})
    return observations, {"sha256": hashlib.sha256(data).hexdigest(), "raw_rows": rows,
                          "unique_sequences": len(observations), "replicate_policy": "Arithmetic mean within exact sequence and split; equal weight per sequence",
                          "outside_model_support_rows": outside_support_rows,
                          "outside_model_support_sequences": sum(row["outside_model_support"] for row in observations),
                          "outside_model_support_aligned_sequences": [row["aligned_sequence"] for row in observations if row["outside_model_support"]],
                          "support_policy": "Measured sequences may lie outside Potts support; aligned length, alphabet and anchor gap pattern must match. Candidate support remains enforced."}


def _features(sequences, alphabet, length):
    if length * len(alphabet) > MAX_FEATURES:
        raise ValueError(f"One hot assay features exceed the {MAX_FEATURES} feature limit")
    indices = {symbol: i for i, symbol in enumerate(alphabet)}
    result = np.zeros((len(sequences), length * len(alphabet)), dtype=np.float64)
    for row, sequence in enumerate(sequences):
        result[row, [site * len(alphabet) + indices[symbol] for site, symbol in enumerate(sequence)]] = 1.0
    return result


def _training_values(observations):
    training = [row for row in observations if row["split"] == "train"]
    if len(training) < 2:
        raise ValueError("Assay prediction requires at least two distinct training sequences")
    y_train = np.asarray([row["value"] for row in training], dtype=np.float64)
    if not np.isfinite(y_train.std()) or y_train.std() <= 0:
        raise ValueError("Assay prediction requires nonconstant finite training measurements")
    return training, y_train


def _training_allocation(rows, features):
    feature_bytes = rows * features * 8
    covariance_bytes = features * features * 8
    estimate = 3 * feature_bytes + 2 * covariance_bytes
    if features > MAX_FEATURES or estimate > MAX_TRAINING_WORKING_BYTES:
        raise ValueError("Assay training exceeds the explicit feature and covariance memory budget "
                         f"({estimate} estimated bytes, limit {MAX_TRAINING_WORKING_BYTES})")
    return {"training_feature_bytes": feature_bytes, "coefficient_covariance_bytes": covariance_bytes,
            "estimated_working_bytes": estimate, "working_byte_limit": MAX_TRAINING_WORKING_BYTES,
            "estimate_scope": "Three dense feature arrays and two coefficient covariance arrays; library overhead is additional"}


def validate_design_settings(settings: dict, anchor: str, alphabet: str, has_assays: bool = False) -> dict:
    """Return normalized strict settings without fitting a model or sampling."""
    return _configuration(settings, anchor, alphabet, has_assays)


def validate_design_inputs(*, anchor: str, alphabet: str, settings: dict, assays_path: Path | None = None) -> dict:
    """Check configuration and assay evidence before an expensive sequence run."""
    config = _configuration(settings, anchor, alphabet, assays_path is not None)
    evidence = None
    if assays_path is not None:
        observations, evidence = _read_assays(assays_path, anchor, alphabet, config)
        if config["objective"] == "assay_ucb":
            training, _ = _training_values(observations)
            evidence["training_allocation"] = _training_allocation(len(training), len(anchor) * len(alphabet))
    return {"settings": config, "assays": evidence}


def _fit_predictor(observations, candidate_sequences, alphabet, length):
    training, y_train = _training_values(observations)
    allocation = _training_allocation(len(training), length * len(alphabet))
    x_train = _features([row["aligned_sequence"] for row in training], alphabet, length)
    x_center = x_train.mean(axis=0)
    y_center = float(y_train.mean())
    y_scale = float(y_train.std())
    if not math.isfinite(y_scale) or y_scale <= 0:
        raise ValueError("Assay prediction requires nonconstant finite training measurements")
    # Explicit centering keeps predictive covariance in the fitted feature coordinates.
    predictor = BayesianRidge(fit_intercept=False)
    predictor.fit(x_train - x_center, (y_train - y_center) / y_scale)

    def predict(sequences):
        if not sequences:
            return np.empty(0), np.empty(0)
        mean_parts, deviation_parts = [], []
        for start in range(0, len(sequences), 256):
            means, deviations = predictor.predict(_features(sequences[start:start + 256], alphabet, length) - x_center, return_std=True)
            mean_parts.append(means * y_scale + y_center)
            deviation_parts.append(deviations * y_scale)
        means, deviations = np.concatenate(mean_parts), np.concatenate(deviation_parts)
        if not np.all(np.isfinite(means)) or not np.all(np.isfinite(deviations)) or np.any(deviations < 0):
            raise ValueError("Assay prediction produced invalid numerical values")
        return means, deviations

    metrics = {}
    for split in ("train", "validation", "test"):
        rows = [row for row in observations if row["split"] == split]
        predicted, deviation = predict([row["aligned_sequence"] for row in rows])
        measured = np.asarray([row["value"] for row in rows])
        errors = predicted - measured
        metrics[split] = {"sequences": len(rows), "rmse": float(np.sqrt(np.mean(errors**2))) if rows else None,
                          "mae": float(np.mean(np.abs(errors))) if rows else None,
                          "rows": [{"aligned_sequence": row["aligned_sequence"], "measured_mean": row["value"],
                                    "predicted_mean": float(mean), "predictive_std": float(std), "replicates": row["replicates"],
                                    "outside_model_support": row.get("outside_model_support", False)}
                                   for row, mean, std in zip(rows, predicted, deviation)]}
    means, deviations = predict(candidate_sequences)
    fit_evidence = {"model": "sklearn.linear_model.BayesianRidge", "fit_intercept": False,
                    "feature_transform": "Full alphabet one hot; explicitly centered using training rows only",
                    "training_sequences": [row["aligned_sequence"] for row in training],
                    "training_target_mean": y_center, "training_target_scale": y_scale,
                    "allocation": allocation, "prediction_chunk_size": 256,
                    "uncertainty_kind": "Predictive standard deviation including observation noise; uncalibrated",
                    "hyperparameter_selection": "None; validation and test used only for separate reporting",
                    "metrics": metrics}
    return means, deviations, fit_evidence


def _hamming(left, right):
    return sum(a != b for a, b in zip(left, right))


def _objective(indices, utilities, similarity, weight):
    ordered = sorted(int(i) for i in indices)
    return float(math.fsum(float(utilities[i]) for i in ordered)
                 - weight * math.fsum(float(similarity[i, j]) for p, i in enumerate(ordered) for j in ordered[p + 1:]))


def _greedy(candidates, utilities, distances, similarity, config):
    chosen = []
    remaining = set(range(len(candidates)))
    while len(chosen) < config["batch_size"]:
        feasible = [i for i in remaining if all(distances[i, j] >= config["min_pairwise_hamming"] for j in chosen)]
        if not feasible:
            break
        best = min(feasible, key=lambda i: (-(utilities[i] - config["diversity_weight"] * sum(similarity[i, j] for j in chosen)), candidates[i]["candidate_id"]))
        chosen.append(best)
        remaining.remove(best)
    return chosen


def _solve_milp(utilities, distances, similarity, config):
    size = len(utilities)
    pairs = [(i, j) for i in range(size) for j in range(i + 1, size)
             if config["diversity_weight"] > 0 and similarity[i, j] > 0
             and distances[i, j] >= config["min_pairwise_hamming"]]
    conflicts = [(i, j) for i in range(size) for j in range(i + 1, size)
                 if distances[i, j] < config["min_pairwise_hamming"]]
    objective = np.asarray([-float(u) for u in utilities] + [config["diversity_weight"] * similarity[i, j] for i, j in pairs])
    row_indices, column_indices, coefficients = [], [], []
    lower, upper = [], []

    def constraint(entries, lo, hi):
        row = len(lower)
        for column, coefficient in entries:
            row_indices.append(row)
            column_indices.append(column)
            coefficients.append(coefficient)
        lower.append(lo)
        upper.append(hi)

    constraint([(i, 1.0) for i in range(size)], config["batch_size"], config["batch_size"])
    for i, j in conflicts:
        constraint([(i, 1.0), (j, 1.0)], -np.inf, 1.0)
    for p, (i, j) in enumerate(pairs, size):
        constraint([(p, 1.0), (i, -1.0), (j, -1.0)], -1.0, np.inf)
        constraint([(p, 1.0), (i, -1.0)], -np.inf, 0.0)
        constraint([(p, 1.0), (j, -1.0)], -np.inf, 0.0)
    matrix = coo_array((np.asarray(coefficients), (np.asarray(row_indices, dtype=np.int32), np.asarray(column_indices, dtype=np.int32))),
                       shape=(len(lower), len(objective))).tocsc()
    result = milp(objective, integrality=np.ones(len(objective), dtype=np.int32),
                  bounds=Bounds(np.zeros(len(objective)), np.ones(len(objective))),
                  constraints=LinearConstraint(matrix, lower, upper),
                  options={"time_limit": config["max_seconds"], "mip_rel_gap": 0.0})
    chosen = [] if result.x is None else list(np.flatnonzero(result.x[:size] > 0.5))
    feasible = (len(chosen) == config["batch_size"] and all(distances[i, j] >= config["min_pairwise_hamming"]
                for p, i in enumerate(chosen) for j in chosen[p + 1:]))
    evidence = {"solver": "scipy.optimize.milp (HiGHS)", "status_code": int(result.status), "message": str(result.message),
                "proven_optimal": bool(result.status == 0 and feasible), "feasible_incumbent": feasible,
                "binary_variables": len(objective), "pair_variables": len(pairs), "hard_conflicts": len(conflicts),
                "time_limit_seconds": config["max_seconds"],
                "mip_gap": _finite_optional(getattr(result, "mip_gap", None)),
                "mip_dual_bound": _finite_optional(getattr(result, "mip_dual_bound", None))}
    return chosen if feasible else [], evidence


def _finite_optional(value):
    return float(value) if value is not None and math.isfinite(float(value)) else None


def _solve_thrml(candidates, utilities, distances, similarity, config, compiled, greedy):
    """Run official THRML Gibbs sampling as a bounded optimization heuristic."""
    started = time.perf_counter()
    import importlib.metadata
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from thrml import Block, SamplingSchedule, SpinNode, sample_states
    from thrml.models import IsingEBM, IsingSamplingProgram

    size = len(candidates)
    biases = np.asarray(compiled["thrml_biases"], dtype=np.float64)
    pair_terms = compiled["pair_terms"]
    weights = np.asarray([edge["thrml_weight"] for edge in pair_terms], dtype=np.float64)
    local_bound = np.abs(biases).copy()
    for edge, weight in zip(pair_terms, weights):
        local_bound[edge["i"]] += abs(weight)
        local_bound[edge["j"]] += abs(weight)
    scale = max(1.0, float(local_bound.max()))
    if not math.isfinite(scale):
        raise ValueError("THRML batch objective has an invalid coefficient scale")
    logit_bound = 2 * config["optimizer_beta"] * float(local_bound.max()) / scale
    if logit_bound > 30.0 + 1e-12:
        raise ValueError("THRML batch Gibbs logits exceed the declared numerical bound")
    greedy_complete = len(greedy) == config["batch_size"]
    initial_indices = greedy if greedy_complete else list(range(config["batch_size"]))
    initial_word = np.zeros(size, dtype=bool)
    initial_word[initial_indices] = True
    cpu = jax.devices("cpu")[0]
    with jax.default_device(cpu):
        nodes = [SpinNode() for _ in range(size)]
        edges = [(nodes[edge["i"]], nodes[edge["j"]]) for edge in pair_terms]
        model = IsingEBM(nodes, edges, jnp.asarray(biases / scale, dtype=jnp.float64),
                         jnp.asarray(weights / scale, dtype=jnp.float64),
                         jnp.asarray(config["optimizer_beta"], dtype=jnp.float64))
        # Dense cardinality couplings mean that singleton blocks are required.
        blocks = [Block([node]) for node in nodes]
        program = IsingSamplingProgram(model, blocks, clamped_blocks=[])
        initial = [jnp.asarray([bit], dtype=jnp.bool_) for bit in initial_word]
        schedule = SamplingSchedule(n_warmup=config["optimizer_warmup"], n_samples=config["optimizer_draws"], steps_per_sample=1)
        key = jax.random.key(config["optimizer_seed"], impl="threefry2x32")
        returned = sample_states(key, program, schedule, initial, [], [Block(nodes)])
        jax.block_until_ready(returned)
        words = np.asarray(returned[0], dtype=bool)
    if words.shape != (config["optimizer_draws"], size):
        raise ValueError("THRML batch optimizer returned an unexpected sample shape")
    valid = words.sum(axis=1) == config["batch_size"]
    for i in range(size):
        for j in range(i + 1, size):
            if distances[i, j] < config["min_pairwise_hamming"]:
                valid &= ~(words[:, i] & words[:, j])
    best_sample, best_value = [], None
    for word in words[valid]:
        indices = list(np.flatnonzero(word))
        objective = _objective(indices, utilities, similarity, config["diversity_weight"])
        if best_value is None or objective > best_value:
            best_sample, best_value = indices, objective
    baseline_value = _objective(greedy, utilities, similarity, config["diversity_weight"]) if greedy_complete else None
    improvement = best_value is not None and (baseline_value is None or best_value > baseline_value)
    if improvement:
        chosen, source = best_sample, "thrml_best_feasible_sample"
    elif greedy_complete:
        chosen, source = greedy, "greedy_incumbent"
    else:
        chosen, source = [], "no_feasible_batch"
    elapsed = time.perf_counter() - started
    evidence = {"solver": "Official THRML IsingSamplingProgram", "backend": "thrml_simulator_cpu",
                "physical_device_used": False, "proven_optimal": False, "selection_source": source,
                "fallback_used": source == "greedy_incumbent", "thrml_improved_greedy": bool(improvement and greedy_complete),
                "best_feasible_sample_objective": best_value, "greedy_objective": baseline_value,
                "total_samples": len(words), "valid_samples": int(valid.sum()),
                "invalid_samples": int((~valid).sum()), "valid_fraction": float(valid.mean()),
                "invalid_fraction": float((~valid).mean()), "unique_valid_words": int(len(np.unique(words[valid], axis=0))),
                "seed": config["optimizer_seed"], "rng": "threefry2x32",
                "warmup_sweeps": config["optimizer_warmup"], "sweeps_per_sample": 1,
                "block_schedule": "One singleton spin block per candidate, sequential order",
                "initialization": "Legal greedy batch" if greedy_complete else "Top utility cardinality batch; pair conflicts may be present",
                "initial_word": initial_word.astype(int).tolist(),
                "coefficient_scale_divisor": scale, "programmed_beta": config["optimizer_beta"],
                "effective_beta_for_original_objective": config["optimizer_beta"] / scale,
                "worst_case_absolute_gibbs_logit": logit_bound,
                "programmed_biases": (biases / scale).tolist(), "programmed_weights": (weights / scale).tolist(),
                "coefficient_clipping": False, "dtype": "float64", "jax_device": str(cpu),
                "thrml_version": importlib.metadata.version("thrml"), "jax_version": jax.__version__,
                "elapsed_seconds": elapsed, "requested_seconds": config["max_seconds"],
                "budget_exceeded": elapsed > config["max_seconds"],
                "budget_scope": "Bounded draw and warmup counts; JAX compilation and sampling calls are synchronous and not interrupted",
                "claim_scope": "Optimization heuristic with feasible sample filtering and a separately labelled classical incumbent; no equilibrium sampling claim"}
    return chosen, evidence


def _compile_report(candidates, utilities, distances, similarity, config):
    """Export a penalized binary objective and its exactly equivalent Ising form."""
    size, k = len(candidates), config["batch_size"]
    gamma = config["diversity_weight"]
    local_bound = max((abs(float(utilities[i])) + gamma * float(similarity[i].sum() - similarity[i, i]) for i in range(size)), default=0.0)
    penalty = local_bound + max(1.0, 0.01 * local_bound)
    # More than the full unpenalized objective range excludes any conflicting optimum.
    unpenalized_range = float(np.abs(utilities).sum() + gamma * np.triu(similarity, 1).sum())
    conflict_penalty = unpenalized_range + max(1.0, 0.01 * unpenalized_range)
    linear = -utilities + penalty * (1 - 2 * k)
    pairwise = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            pairwise[i, j] = gamma * similarity[i, j] + 2 * penalty
            if distances[i, j] < config["min_pairwise_hamming"]:
                pairwise[i, j] += conflict_penalty
    # The conflict penalty changes the local repair bound. Recompute a safe size penalty.
    if config["min_pairwise_hamming"] > 0:
        conflict_counts = (distances < config["min_pairwise_hamming"]).sum(axis=1) - 1
        bound_with_conflicts = max((abs(float(utilities[i])) + gamma * float(similarity[i].sum() - similarity[i, i])
                                    + conflict_penalty * int(conflict_counts[i]) for i in range(size)), default=0.0)
        revised_penalty = bound_with_conflicts + max(1.0, 0.01 * bound_with_conflicts)
        pairwise[np.triu_indices(size, 1)] += 2 * (revised_penalty - penalty)
        penalty = revised_penalty
        linear = -utilities + penalty * (1 - 2 * k)
    ising_fields = linear / 2 + (pairwise.sum(axis=0) + pairwise.sum(axis=1)) / 4
    edges = [{"i": i, "j": j, "qubo": float(pairwise[i, j]), "thrml_weight": float(-pairwise[i, j] / 4)}
             for i in range(size) for j in range(i + 1, size)]
    constant = penalty * k * k
    return {"kind": "Finite penalty export for an optimization objective; not a physical device execution",
            "variable_candidate_ids": [row["candidate_id"] for row in candidates],
            "binary_variables": size, "edges": len(edges), "maximum_degree": max(0, size - 1),
            "quadratic_pair_convention": "Each pair i < j counted once",
            "cardinality_penalty": penalty, "conflict_penalty": conflict_penalty,
            "qubo_constant": constant, "qubo_linear": linear.tolist(), "pair_terms": edges,
            "thrml_biases": (-ising_fields).tolist(),
            "ising_constant": float(constant + linear.sum() / 2 + pairwise.sum() / 4),
            "bit_to_spin": "s = 2*z - 1; False maps to -1 and True maps to +1",
            "sampling_note": "Finite penalties permit infeasible thermal samples; decode and check every returned state",
            "hardware_status": "Simulation specification only; dense graph embedding and device coefficient limits unverified"}


def design_library(pool: list[dict], *, anchor: str, alphabet: str, settings: dict, assays_path: Path | None = None) -> dict:
    """Rank validated sequence candidates and select a constrained experimental batch.

    ``allowed_residues`` must be supplied from original model support by the caller
    for candidate validation and annotation of assay controls outside prior support.
    Assay evidence may be outside that support but must match the anchor gap pattern.
    Motifs and residue counts refer to ungapped sequences; Hamming
    distances refer to aligned sequences. No assay labels are synthesized.
    """
    started = time.perf_counter()
    config = _configuration(settings, anchor, alphabet, assays_path is not None)
    if not isinstance(pool, list):
        raise ValueError("Candidate pool must be a list")
    observations, assay_evidence = ([], None) if assays_path is None else _read_assays(assays_path, anchor, alphabet, config)
    measured = {row["aligned_sequence"] for row in observations}
    rows, seen_ids = [], set()
    for original in pool:
        if not isinstance(original, dict):
            raise ValueError("Every candidate must be an object")
        row = dict(original)
        candidate_id = _text(row.get("candidate_id"), "candidate_id")
        if candidate_id in seen_ids:
            raise ValueError("Candidate identifiers must be unique")
        seen_ids.add(candidate_id)
        _validate_sequence(row.get("aligned_sequence"), anchor, alphabet, config.get("allowed_residues"), candidate_id)
        if row.get("sequence") != row["aligned_sequence"].replace("-", ""):
            raise ValueError("Candidate ungapped sequence does not match its aligned sequence")
        row["energy"] = _number(row.get("energy"), "energy", -1e100, 1e100)
        row.update({"selected": False, "selection_rank": None, "hamming_to_anchor": _hamming(row["aligned_sequence"], anchor)})
        rows.append(row)
    evidence = {"objective": config["objective"], "direction": config["direction"], "settings": config,
                "assays": assay_evidence, "physical_device_used": False,
                "sequence_semantics": "Hamming on aligned sequences; literal motifs and composition on ungapped sequences",
                "model_support_validation": "Candidate support enforced; assays outside prior support annotated" if config.get("allowed_residues") is not None else "Alphabet and finite candidate energies only; caller validates model support"}
    if config["objective"] == "assay_ucb":
        means, deviations, fit_evidence = _fit_predictor(observations, [row["aligned_sequence"] for row in rows], alphabet, len(anchor))
        evidence["predictor"] = fit_evidence
        sign = 1.0 if config["direction"] == "maximize" else -1.0
        for row, mean, deviation in zip(rows, means, deviations):
            row.update({"predicted_mean": float(mean), "predictive_std": float(deviation),
                        "uncertainty_kind": "Predictive standard deviation including observation noise; uncalibrated",
                        "utility": float(sign * mean + config["exploration_weight"] * deviation)})
    else:
        evidence["scientific_claim"] = "Negative Potts energy is a model score proxy; no measured function claim"
        for row in rows:
            row.update({"predicted_mean": None, "predictive_std": None, "uncertainty_kind": None, "utility": -row["energy"]})
    if any(not math.isfinite(row["utility"]) for row in rows):
        raise ValueError("Acquisition utility must be finite")
    exclusions, eligible, seen_sequences = [], [], set()
    for row in rows:
        sequence, reasons = row["sequence"], []
        if row["aligned_sequence"] in seen_sequences:
            reasons.append("duplicate_aligned_sequence")
        seen_sequences.add(row["aligned_sequence"])
        if row["aligned_sequence"] in measured:
            reasons.append("already_measured")
        if not config["min_hamming"] <= row["hamming_to_anchor"] <= config["max_hamming"]:
            reasons.append("hamming_to_anchor")
        if any(motif in sequence for motif in config["forbidden_motifs"]):
            reasons.append("forbidden_motif")
        if any(motif not in sequence for motif in config["required_motifs"]):
            reasons.append("required_motif_missing")
        residue_count = sum(symbol in config["residue_class"] for symbol in sequence)
        if config["residue_class"] and not config["min_residue_count"] <= residue_count <= config["max_residue_count"]:
            reasons.append("residue_class_count")
        row["exclusion_reasons"] = reasons
        if reasons:
            exclusions.append({"candidate_id": row["candidate_id"], "reasons": reasons})
        else:
            eligible.append(row)
    eligible.sort(key=lambda row: (-row["utility"], row["candidate_id"]))
    evidence["eligible_before_pool_limit"] = len(eligible)
    for row in eligible[config["max_pool"]:]:
        row["exclusion_reasons"] = ["pool_limit"]
        exclusions.append({"candidate_id": row["candidate_id"], "reasons": ["pool_limit"]})
    eligible = eligible[:config["max_pool"]]
    evidence["optimization_scope"] = "Top utility capped candidate pool; optimality applies only to this pool"
    evidence["optimization_pool_ids"] = [row["candidate_id"] for row in eligible]
    if len(eligible) < config["batch_size"]:
        evidence["elapsed_seconds"] = time.perf_counter() - started
        return {"status": "insufficient_candidates", "selected": [], "all_candidates": rows, "exclusions": exclusions, "evidence": evidence}
    if config["batch_size"] == 0:
        evidence["elapsed_seconds"] = time.perf_counter() - started
        return {"status": "completed", "selected": [], "all_candidates": rows, "exclusions": exclusions, "evidence": evidence}
    utilities = np.asarray([row["utility"] for row in eligible], dtype=np.float64)
    distances = np.asarray([[_hamming(a["aligned_sequence"], b["aligned_sequence"]) for b in eligible] for a in eligible])
    similarity = 1.0 - distances / len(anchor)
    greedy_started = time.perf_counter()
    greedy = _greedy(eligible, utilities, distances, similarity, config)
    evidence["greedy"] = {"candidate_ids": [eligible[i]["candidate_id"] for i in greedy],
                          "complete": len(greedy) == config["batch_size"],
                          "objective": _objective(greedy, utilities, similarity, config["diversity_weight"]) if len(greedy) == config["batch_size"] else None,
                          "elapsed_seconds": time.perf_counter() - greedy_started}
    evidence["compiled_batch_objective"] = _compile_report(eligible, utilities, distances, similarity, config)
    if config["optimizer"] == "milp":
        solver_started = time.perf_counter()
        selected, evidence["solver"] = _solve_milp(utilities, distances, similarity, config)
        evidence["solver"]["elapsed_seconds"] = time.perf_counter() - solver_started
        if not selected:
            status = "constraints_infeasible" if evidence["solver"]["status_code"] == 2 else "no_feasible_solution_within_budget"
        else:
            status = "completed" if evidence["solver"]["proven_optimal"] else "feasible_budget_limited"
    elif config["optimizer"] == "thrml":
        if len(eligible) == 1:
            selected = [0]
            evidence["solver"] = {"solver": "Unique feasible candidate", "backend": "deterministic_no_sampling_required",
                                    "physical_device_used": False, "selection_source": "unique_feasible_batch",
                                    "proven_optimal": True, "total_samples": 0, "fallback_used": False,
                                    "thrml_improved_greedy": False, "budget_exceeded": False}
            status = "completed"
        else:
            selected, evidence["solver"] = _solve_thrml(eligible, utilities, distances, similarity, config,
                                                        evidence["compiled_batch_objective"], greedy)
            status = ("feasible_budget_limited" if evidence["solver"]["budget_exceeded"] else "completed") if selected else "no_feasible_solution_within_budget"
    else:
        selected = greedy if len(greedy) == config["batch_size"] else []
        evidence["solver"] = {"solver": "Deterministic greedy marginal utility", "proven_optimal": False}
        status = "completed" if selected else "greedy_incomplete"
    chosen_rows = sorted((eligible[i] for i in selected), key=lambda row: (-row["utility"], row["candidate_id"]))
    for rank, row in enumerate(chosen_rows, 1):
        row.update({"selected": True, "selection_rank": rank,
                    "selection_reason": evidence["solver"].get("selection_source", config["optimizer"]) + " constrained batch objective"})
    evidence["selected_objective"] = _objective(selected, utilities, similarity, config["diversity_weight"]) if selected else None
    evidence["objective_gain_over_greedy"] = (evidence["selected_objective"] - evidence["greedy"]["objective"]
                                               if selected and evidence["greedy"]["complete"] else None)
    evidence["elapsed_seconds"] = time.perf_counter() - started
    return {"status": status, "selected": chosen_rows, "all_candidates": rows, "exclusions": exclusions, "evidence": evidence}
