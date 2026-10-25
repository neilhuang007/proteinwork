"""Sampling diagnostics with explicit evidence and failure states."""

from __future__ import annotations

import math
import warnings

import numpy as np
from arviz_stats.base import array_stats

from .direct import enumerate_target


def _finite(value):
    value = float(value)
    return value if math.isfinite(value) else None


def summarize_observable(values, *, proven_value: float | None = None) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or not np.all(np.isfinite(values)):
        raise ValueError("Diagnostic observations must be finite arrays of chains and draws")
    base = {"mean": float(values.mean()), "sd": float(values.std(ddof=1)),
            "chains": values.shape[0], "draws_per_chain": values.shape[1]}
    if proven_value is not None:
        correct = bool(np.allclose(values, proven_value, rtol=0, atol=1e-9))
        return {**base, "constant_proven": True, "expected_value": float(proven_value),
                "constant_value_verified": correct, "rhat": None, "ess_bulk": None,
                "ess_tail": None, "mcse": 0.0 if correct else None,
                "status": "constant_verified" if correct else "constant_value_mismatch"}
    if values.shape[0] < 4 or values.shape[1] < 4 or np.ptp(values) == 0:
        return {**base, "constant_proven": False, "rhat": None, "ess_bulk": 0.0,
                "ess_tail": 0.0, "mcse": None, "status": "insufficient_or_constant_trace"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        stats = {"rhat": _finite(array_stats.rhat(values, method="rank")),
                 "ess_bulk": _finite(array_stats.ess(values, method="bulk")),
                 "ess_tail": _finite(array_stats.ess(values, method="tail", prob=(0.05, 0.95))),
                 "mcse": _finite(array_stats.mcse(values, method="mean"))}
    return {**base, **stats, "constant_proven": False,
            "status": "computed" if all(v is not None for v in stats.values()) else "undefined_diagnostic"}


def _observable_pass(summary, targets, reference=False):
    if summary["constant_proven"]:
        return summary["constant_value_verified"]
    if summary["status"] != "computed":
        return False
    ess = targets["bulk_and_tail_ess_min"] if reference else targets["bulk_ess_min_inclusive"]
    tail = ess if reference else targets["tail_ess_min_inclusive"]
    return (summary["rhat"] < targets["rank_normalized_split_rhat_max_exclusive"]
            and summary["ess_bulk"] >= ess and summary["ess_tail"] >= tail)


def _indicator_mcse(values, constant):
    if constant:
        return 0.0
    if np.all(values == values.flat[0]) or values.shape[1] < 4:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _finite(array_stats.mcse(values.astype(float), method="mean"))


def diagnose(model, result, *, beta, anchor: str, selected_pairs, evaluation: dict,
             reference: dict | None = None, exact_max_states=4096) -> dict:
    target_count = math.prod(len(d) for d in model.domains)
    full = model.expand(result.states)
    anchor_tokens = model.source.encode(anchor)
    hamming = np.count_nonzero(full != anchor_tokens, axis=-1)
    exact_states = exact_probs = exact_energies = exact_hamming = None
    if target_count <= exact_max_states:
        exact_states, exact_probs = enumerate_target(model, beta, max_states=exact_max_states)
        exact_energies = np.asarray(model.energy(exact_states))
        exact_hamming = np.count_nonzero(model.expand(exact_states) != anchor_tokens, axis=-1)
    all_domains = [(int(t),) for t in model.fixed_state]
    for i, original in enumerate(model.original_indices):
        all_domains[original] = model.domains[i]
    distance_proven = all(len(d) == 1 or int(anchor_tokens[i]) not in d
                          for i, d in enumerate(all_domains))
    energy_proven = (exact_energies is not None and np.all(exact_energies == exact_energies[0]))
    observables = {
        "energy": summarize_observable(result.energies,
                 proven_value=float(exact_energies[0]) if energy_proven else None),
        "hamming": summarize_observable(hamming,
                 proven_value=float(sum(int(anchor_tokens[i]) not in d for i, d in enumerate(all_domains)))
                 if distance_proven else None)}
    targets = evaluation["sampling_accuracy_targets"]
    checks = {f"{name}_diagnostics": _observable_pass(summary, targets)
              for name, summary in observables.items()}
    exact_full = model.expand(exact_states) if exact_states is not None else None
    frequencies = []
    for site, domain in enumerate(all_domains):
        for token in domain:
            values = full[..., site] == token
            row = {"key": f"site:{site}:{token}", "kind": "site", "model_indices": [site],
                   "symbols": model.alphabet[token], "probability": float(values.mean()),
                   "mcse": _indicator_mcse(values, len(domain) == 1),
                   "target_constant": len(domain) == 1}
            if exact_full is not None:
                row["reference_probability"] = float(exact_probs[exact_full[:, site] == token].sum())
            frequencies.append(row)
    seen = set()
    for pair in selected_pairs:
        if (not isinstance(pair, (list, tuple)) or len(pair) != 2
                or any(type(i) is not int or not 0 <= i < full.shape[-1] for i in pair)
                or pair[0] >= pair[1] or tuple(pair) in seen):
            raise ValueError("Selected pair indices must be unique ordered model index pairs")
        seen.add(tuple(pair))
        i, j = pair
        for a in all_domains[i]:
            for b in all_domains[j]:
                values = (full[..., i] == a) & (full[..., j] == b)
                constant = len(all_domains[i]) == len(all_domains[j]) == 1
                row = {"key": f"pair:{i}:{j}:{a}:{b}", "kind": "pair", "model_indices": [i, j],
                       "symbols": model.alphabet[a] + model.alphabet[b], "probability": float(values.mean()),
                       "mcse": _indicator_mcse(values, constant), "target_constant": constant}
                if exact_full is not None:
                    row["reference_probability"] = float(exact_probs[(exact_full[:, i] == a) & (exact_full[:, j] == b)].sum())
                frequencies.append(row)
    comparison = {"kind": "exact" if exact_states is not None else "unavailable"}
    if exact_states is not None:
        lookup = {tuple(state): i for i, state in enumerate(exact_states)}
        counts = np.zeros(len(exact_states))
        for row in result.states.reshape(result.states.shape[0] * result.states.shape[1], model.n_sites):
            counts[lookup[tuple(row)]] += 1
        empirical = counts / counts.sum()
        comparison["total_variation"] = float(np.abs(empirical - exact_probs).sum() / 2)
        comparison["state_count"] = len(exact_states)
        comparison["state_probabilities"] = [
            {"state": list(map(int, state)), "exact": float(p), "observed": float(e)}
            for state, p, e in zip(exact_states, exact_probs, empirical)]
        checks["exact_total_variation"] = comparison["total_variation"] <= targets["exact_case_total_variation_max"]
        for name, values in (("energy", exact_energies), ("hamming", exact_hamming)):
            mean = float(np.dot(exact_probs, values))
            sd = float(np.sqrt(np.dot(exact_probs, (values - mean) ** 2)))
            error = abs(observables[name]["mean"] - mean)
            comparison[name] = {"reference_mean": mean, "reference_sd": sd,
                                "absolute_mean_error": error, "reference_mcse": 0.0}
            bound = targets[f"{name}_mean_error_in_reference_standard_deviations_max"] * sd
            checks[f"{name}_mean_accuracy"] = error <= max(bound, 1e-9 if sd < 1e-12 else 0)
    elif reference is not None:
        if not reference.get("reference_qualification", {}).get("qualified", False):
            raise ValueError("Broad reference failed its independent adequacy checks")
        comparison["kind"] = "qualified_independent_direct"
        ref_rows = {r["key"]: r for r in reference["frequencies"]}
        for row in frequencies:
            if row["key"] not in ref_rows:
                raise ValueError("Reference lacks a requested frequency observable")
            ref = ref_rows[row["key"]]
            row["reference_probability"] = ref["probability"]
            row["reference_mcse"] = ref["mcse"]
        for name, summary in observables.items():
            ref = reference["observables"][name]
            if summary["mcse"] is None or ref["mcse"] is None:
                upper = None
            else:
                upper = abs(summary["mean"] - ref["mean"]) + 1.96 * math.hypot(summary["mcse"], ref["mcse"])
            comparison[name] = {"reference_mean": ref["mean"], "reference_sd": ref["sd"],
                                "reference_mcse": ref["mcse"], "error_upper_bound": upper}
            checks[f"{name}_mean_accuracy"] = upper is not None and upper <= max(
                targets[f"{name}_mean_error_in_reference_standard_deviations_max"] * ref["sd"],
                1e-9 if ref["constant_proven"] else 0)
    else:
        checks["independent_accuracy_reference"] = False
    for kind, limit in (("site", "single_site_probability_absolute_error_max"),
                        ("pair", "selected_pair_probability_absolute_error_max")):
        rows = [r for r in frequencies if r["kind"] == kind]
        bounds = []
        for row in rows:
            if "reference_probability" not in row:
                continue
            error = abs(row["probability"] - row["reference_probability"])
            row["absolute_error"] = error
            if exact_states is None:
                error = (error + 1.96 * math.hypot(row["mcse"], row["reference_mcse"])
                         if row["mcse"] is not None and row.get("reference_mcse") is not None else None)
            row["error_upper_bound"] = error
            bounds.append(error)
        checks[f"{kind}_frequency_accuracy"] = len(bounds) == len(rows) and all(
            b is not None and b <= targets[limit] for b in bounds)
        comparison[f"{kind}_maximum_error_bound"] = max(
            (b for b in bounds if b is not None), default=None)
    ref_targets = evaluation["broad_case_reference"]
    ref_checks = {f"{name}_diagnostics": _observable_pass(summary, ref_targets, reference=True)
                  for name, summary in observables.items()}
    ref_checks["frequency_uncertainty"] = all(r["mcse"] is not None and
        r["mcse"] <= ref_targets["frequency_mcse_max"] for r in frequencies)
    for name, summary in observables.items():
        ref_checks[f"{name}_mean_uncertainty"] = (summary["mcse"] is not None and
            (summary["constant_proven"] or summary["mcse"] <=
             ref_targets["mean_mcse_in_reference_standard_deviations_max"] * summary["sd"]))
    return {"schema_version": 1, "status": "accuracy_passed" if all(checks.values()) else "accuracy_not_reached",
            "checks": checks, "observables": observables, "frequencies": frequencies,
            "comparison": comparison, "reference_qualification": {
                "qualified": all(ref_checks.values()), "checks": ref_checks},
            "diagnostic_library": "arviz_stats.base.array_stats",
            "limitations": ["Diagnostics do not prove convergence or biological function.",
                            "An unobserved nonconstant indicator has unknown MCSE and cannot qualify a broad reference.",
                            "An observation is a complete sweep for direct sampling and a corrected attempt for encoded sampling."]}
