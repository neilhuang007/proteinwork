"""Stdlib acquisition/admission inspection, independent of the Week 5 loader.

This script never trains, changes parameters, samples chains, or asserts that a
synthetic support witness is a natural protein. The ZIP stays ignored locally.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import urllib.request
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def verify_bytes(data: bytes, expected_hash: str, expected_size: int) -> None:
    if len(data) != expected_size:
        raise ValueError(f"Byte count mismatch: expected {expected_size}, got {len(data)}")
    if hashlib.sha256(data).hexdigest() != expected_hash:
        raise ValueError("SHA-256 mismatch; do not inspect or silently accept changed input")


def parse_fasta(data: bytes) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    pieces: list[str] = []
    for line in data.decode("ascii").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(pieces)
            header, pieces = line[1:], []
            if header in records:
                raise ValueError("Duplicate FASTA header")
        elif header is None:
            raise ValueError("Sequence before FASTA header")
        else:
            pieces.append(line)
    if header is not None:
        if header in records:
            raise ValueError("Duplicate FASTA header")
        records[header] = "".join(pieces)
    if not records or any(not s for s in records.values()):
        raise ValueError("Empty FASTA input or sequence")
    return records


def coordinate_map(aligned: str, full_reference: str) -> list[dict]:
    """Map exact nongap substring to full accession positions; gaps have no residue."""
    ungapped = aligned.replace("-", "")
    first = full_reference.find(ungapped)
    if not ungapped or first < 0 or full_reference.find(ungapped, first + 1) >= 0:
        raise ValueError("Reference must have exactly one exact nongap substring match")
    residue = first + 1
    result = []
    for index, symbol in enumerate(aligned):
        result.append({"model_index_0based": index, "alignment_column_1based": index + 1,
                       "anchor_symbol": symbol,
                       "uniprot_position_1based": residue if symbol != "-" else None})
        if symbol != "-":
            residue += 1
    return result


def inspect_parameters(stream, length: int, alphabet: str,
                       sequences: dict[str, str] | None = None) -> dict:
    """Check every expected key once and sum selected terms from original records."""
    if length < 1 or not alphabet or len(set(alphabet)) != len(alphabet):
        raise ValueError("Invalid model dimensions or duplicate alphabet")
    sequences = sequences or {}
    if any(len(s) != length or set(s) - set(alphabet) for s in sequences.values()):
        raise ValueError("Invalid selected sequence length or symbol")
    q = len(alphabet)
    token_index = {a: k for k, a in enumerate(alphabet)}
    expected_j = length * (length - 1) // 2 * q * q
    seen_j = bytearray(expected_j)
    seen_h = bytearray(length * q)
    fields: dict[tuple[int, str], str] = {}
    coupling_terms = {name: [] for name in sequences}
    selected_examples = {name: [] for name in sequences}
    counts = {"J": 0, "h": 0}
    digest = hashlib.sha256()
    byte_count = 0
    for line_number, raw in enumerate(stream, 1):
        digest.update(raw)
        byte_count += len(raw)
        parts = raw.decode("ascii").split()
        if not parts or parts[0].startswith("#"):
            continue
        kind = parts[0]
        if kind not in counts or len(parts) != (6 if kind == "J" else 4):
            raise ValueError(f"Malformed parameter record at line {line_number}")
        i = int(parts[1])
        if not 0 <= i < length:
            raise ValueError("Model position outside declared dimensions")
        value = float(parts[-1])
        if kind == "J":
            j, a, b = int(parts[2]), parts[3], parts[4]
            if not i < j < length or a not in token_index or b not in token_index:
                raise ValueError("Coupling must use i<j and declared symbols")
            if not math.isfinite(value):
                raise ValueError("All couplings must be finite")
            pair = i * (2 * length - i - 1) // 2 + j - i - 1
            key = (pair * q + token_index[a]) * q + token_index[b]
            if seen_j[key]:
                raise ValueError("Duplicate coupling key")
            seen_j[key] = 1
            for name, sequence in sequences.items():
                if sequence[i] == a and sequence[j] == b:
                    coupling_terms[name].append(parts[-1])
                    if len(selected_examples[name]) < 5:
                        selected_examples[name].append({"i": i, "j": j, "a": a,
                                                        "b": b, "J": parts[-1]})
        else:
            a = parts[2]
            if a not in token_index:
                raise ValueError("Unknown field symbol")
            if not math.isfinite(value) and (value != -math.inf or parts[-1].lower() not in {"-inf", "-infinity"}):
                raise ValueError("Only negative-infinity fields have admitted semantics")
            key = i * q + token_index[a]
            if seen_h[key]:
                raise ValueError("Duplicate field key")
            seen_h[key] = 1
            fields[i, a] = parts[-1]
        counts[kind] += 1
    if counts != {"J": expected_j, "h": length * q} or not all(seen_j) or not all(seen_h):
        raise ValueError("Incomplete parameter key coverage")
    support = []
    for i in range(length):
        allowed = "".join(a for a in alphabet if math.isfinite(float(fields[i, a])))
        if not allowed:
            raise ValueError("Position has no supported state")
        support.append({"model_index_0based": i, "allowed_symbols": allowed,
                        "forbidden_symbols": "".join(a for a in alphabet if a not in allowed)})
    score_checks = {}
    for name, sequence in sequences.items():
        selected_h = [fields[i, a] for i, a in enumerate(sequence)]
        forbidden = [{"model_index_0based": i, "symbol": a}
                     for i, a in enumerate(sequence) if float(fields[i, a]) == -math.inf]
        terms = coupling_terms[name]
        if len(terms) != length * (length - 1) // 2:
            raise ValueError("Selected sequence pair-count mismatch")
        pair_sum = sum(map(Decimal, terms), Decimal(0))
        if forbidden:
            score_checks[name] = {"sequence": sequence, "energy": "positive_infinity",
                                  "forbidden_selected_fields": forbidden,
                                  "selected_field_count": length, "selected_pair_count": len(terms)}
        else:
            field_sum = sum(map(Decimal, selected_h), Decimal(0))
            energy = -field_sum - pair_sum
            float_energy = -math.fsum(map(float, selected_h)) - math.fsum(map(float, terms))
            if abs(float(energy) - float_energy) > 1e-10:
                raise ValueError("Decimal and binary64 score disagreement")
            score_checks[name] = {
                "sequence": sequence, "energy_decimal": str(energy), "energy_binary64": float_energy,
                "sum_selected_fields_decimal": str(field_sum),
                "sum_selected_upper_triangle_couplings_decimal": str(pair_sum),
                "selected_field_count": length, "selected_pair_count": len(terms),
                "first_five_pair_terms": selected_examples[name],
                "forbidden_selected_fields": []}
    return {"sha256": digest.hexdigest(), "bytes": byte_count, "record_counts": counts,
            "complete_unique_key_coverage": True, "support": support, "fields": fields,
            "score_checks": score_checks}


def validate_case(case: dict, support: list[dict], anchor: str,
                  witness: str, beta: float) -> list[str]:
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta must be finite and positive")
    if len(anchor) != len(support) or len(witness) != len(support):
        raise ValueError("Case sequences have wrong dimensions")
    domains = [s["allowed_symbols"] if anchor[i] == "-" else s["allowed_symbols"].replace("-", "")
               for i, s in enumerate(support)]
    for i, a in enumerate(anchor):
        if a == "-":
            domains[i] = "-" if "-" in domains[i] else ""
    if case["kind"] == "exact_conditional":
        domains = [witness[i] for i in range(len(witness))]
    elif case["kind"] != "broad_conditional":
        raise ValueError("Unknown restriction case kind")
    changed: set[int] = set()
    for item in case["positional_restrictions"]:
        i, allowed = item["model_index_0based"], item["allowed_symbols"]
        if type(i) is not int or not 0 <= i < len(support) or i in changed:
            raise ValueError("Invalid or duplicate restriction position")
        changed.add(i)
        original_domain = support[i]["allowed_symbols"] if anchor[i] != "-" else "-"
        if anchor[i] != "-":
            original_domain = original_domain.replace("-", "")
        if not allowed or len(set(allowed)) != len(allowed) or any(a not in original_domain for a in allowed):
            raise ValueError("Requested state is forbidden, unknown, empty, or violates gap pattern")
        domains[i] = allowed
    if any(not d or any(a not in support[i]["allowed_symbols"] for a in d) for i, d in enumerate(domains)):
        raise ValueError("Empty or unsupported case domain")
    if any((anchor[i] == "-") != (d == "-") or (anchor[i] != "-" and "-" in d)
           for i, d in enumerate(domains)):
        raise ValueError("Case violates reference gap pattern")
    return domains


def exact_probabilities(scores: dict[str, dict], beta: float = 1.0) -> dict[str, float]:
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta must be finite and positive")
    if not scores or any("energy_decimal" not in score for score in scores.values()):
        raise ValueError("Exact states must all have finite scores")
    logits = {name: -beta * float(score["energy_decimal"]) for name, score in scores.items()}
    maximum = max(logits.values())
    weights = {name: math.exp(value - maximum) for name, value in logits.items()}
    total = math.fsum(weights.values())
    return {name: value / total for name, value in weights.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="download pinned archive if missing")
    args = parser.parse_args()
    config = json.loads((ROOT / "configs/benchmark.json").read_text(encoding="utf-8"))
    archive_config = config["archive"]
    archive_path = ROOT / archive_config["local_path"]
    if not archive_path.exists():
        if not args.download:
            parser.error("Archive missing; run once with --download")
        data = urllib.request.urlopen(archive_config["url"], timeout=60).read(archive_config["bytes"] + 1)
        verify_bytes(data, archive_config["sha256"], archive_config["bytes"])
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        archive_path.write_bytes(data)
    archive_data = archive_path.read_bytes()
    verify_bytes(archive_data, archive_config["sha256"], archive_config["bytes"])
    metadata = json.loads((ROOT / config["coordinate_reference"]["metadata"]).read_text(encoding="utf-8"))
    cases = json.loads((ROOT / config["restriction_cases"]).read_text(encoding="utf-8"))
    model = config["model"]
    with zipfile.ZipFile(io.BytesIO(archive_data)) as archive:
        natural_data = archive.read(config["coordinate_reference"]["archive_member"])
        natural = parse_fasta(natural_data)
        anchor = natural[config["coordinate_reference"]["archive_header"]]
        fixture = natural[config["finite_score_fixture"]["archive_header"]]
        if anchor != metadata["aligned_sequence"]:
            raise ValueError("Anchor differs from frozen metadata")
        mapping = coordinate_map(anchor, metadata["full_accession_sequence"])
        actual_gaps = [i for i, symbol in enumerate(anchor) if symbol == "-"]
        if (actual_gaps != config["coordinate_reference"]["gap_model_indices"] or
                actual_gaps != cases["common_constraints"]["fixed_gap_model_indices"] or
                cases["common_constraints"]["fixed_gap_symbol"] != "-" or
                not cases["common_constraints"]["nongap_columns_exclude_gap"] or
                not cases["common_constraints"]["finite_original_field_support_required"]):
            raise ValueError("Declared common constraints disagree with anchor gap pattern/support")
        witness = cases["computational_witness"]["sequence"]
        sequences = {"native_coordinate_anchor": anchor, "unmodified_finite_fixture": fixture,
                     "computational_witness": witness}
        # Enumerate only the declared varying positions; retain the full 96-site background.
        for case in cases["cases"]:
            if case["kind"] == "exact_conditional":
                varying = case["positional_restrictions"]
                for n, residues in enumerate(itertools.product(*(x["allowed_symbols"] for x in varying))):
                    state = list(witness)
                    for item, symbol in zip(varying, residues):
                        state[item["model_index_0based"]] = symbol
                    sequences[f'{case["case_id"]}/state_{n}'] = "".join(state)
        with archive.open(archive_config["member"]) as stream:
            inspection = inspect_parameters(stream, model["positions"], model["alphabet"], sequences)
        if (inspection["sha256"] != archive_config["member_sha256"] or
                inspection["bytes"] != archive_config["member_bytes"]):
            raise ValueError("Parameter member differs from pinned model")
        members = [{"name": entry.filename, "bytes": entry.file_size} for entry in archive.infolist()]
    support = inspection["support"]
    domain_checks = []
    for case in cases["cases"]:
        domains = validate_case(case, support, anchor, witness, model["beta"])
        domain_checks.append({"case_id": case["case_id"], "split": case["split"],
                              "nonempty_support": True, "state_count": str(math.prod(map(len, domains))),
                              "mutable_positions": [i for i, d in enumerate(domains) if len(d) > 1],
                              "domain_sizes": list(map(len, domains))})
    development = {item["model_index_0based"] for case in cases["cases"] if case["split"] == "development"
                   for item in case["positional_restrictions"]}
    evaluation = {item["model_index_0based"] for case in cases["cases"] if case["split"] == "evaluation"
                  for item in case["positional_restrictions"]}
    if development & evaluation:
        raise ValueError("Development and evaluation selected restriction positions overlap")
    output = ROOT / "data/benchmark"
    with (output / "position_map.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping[0]))
        writer.writeheader()
        writer.writerows(mapping)
    write_json(output / "support.json", {"semantics": "finite original h is permitted; -inf original h is forbidden",
                                        "positions": support})
    scores = inspection["score_checks"]
    exact_cases = {}
    for case in cases["cases"]:
        if case["kind"] == "exact_conditional":
            selected = {name: score for name, score in scores.items() if name.startswith(case["case_id"] + "/")}
            exact_cases[case["case_id"]] = {"state_count": len(selected), "beta": model["beta"],
                                          "probabilities_binary64": exact_probabilities(selected, model["beta"])}
    write_json(output / "score_checks.json", {
        "method": "selected original h and i<j J terms; exact Decimal sums plus independent binary64 fsum",
        "pair_count_per_full_sequence": model["positions"] * (model["positions"] - 1) // 2,
        "role_limit": "one software inspection; not two independent lead attestations or production-loader verification",
        "scores": scores, "exact_conditional_cases": exact_cases})
    audit = {"inspection_date": config["inspection_date"], "archive_sha256": archive_config["sha256"],
             "archive_bytes": len(archive_data), "archive_members": members,
             "parameter_sha256": inspection["sha256"], "parameter_bytes": inspection["bytes"],
             "record_counts": inspection["record_counts"], "complete_unique_key_coverage": True,
             "negative_infinity_fields": sum(len(s["forbidden_symbols"]) for s in support),
             "positions_with_forbidden_states": sum(bool(s["forbidden_symbols"]) for s in support),
             "natural_alignment_sha256": hashlib.sha256(natural_data).hexdigest(),
             "natural_alignment_records": len(natural),
             "reference_exact_substring_match": True, "reference_mapped_nongaps": sum(x["anchor_symbol"] != "-" for x in mapping),
             "reference_residue_span_1based": [min(x["uniprot_position_1based"] for x in mapping if x["uniprot_position_1based"]),
                                              max(x["uniprot_position_1based"] for x in mapping if x["uniprot_position_1based"])],
             "restriction_checks": domain_checks, "selected_development_positions": sorted(development),
             "selected_evaluation_positions": sorted(evaluation), "selected_positions_disjoint": True,
             "production_loader_implemented": False}
    write_json(output / "inspection.json", audit)
    print(json.dumps({"model": config["benchmark_id"], "complete_unique_key_coverage": True,
                      "negative_infinity_fields": audit["negative_infinity_fields"],
                      "reference_span": audit["reference_residue_span_1based"],
                      "native_anchor_energy": scores["native_coordinate_anchor"]["energy"],
                      "finite_fixture_energy": scores["unmodified_finite_fixture"]["energy_decimal"],
                      "cases": [{"case_id": x["case_id"], "state_count": x["state_count"]} for x in domain_checks]}, indent=2))


if __name__ == "__main__":
    main()
