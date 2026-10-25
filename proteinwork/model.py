"""Validated Potts parameters and exact conditional restrictions."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import hashlib
import io
import json
import math
from pathlib import Path
from typing import BinaryIO, Sequence
import zipfile

import numpy as np


def _validateArrays(alphabet, fields, couplings):
    if not isinstance(alphabet, str) or not alphabet or len(set(alphabet)) != len(alphabet):
        raise ValueError("Alphabet must contain distinct single character symbols")
    fieldArray = np.array(fields, dtype=np.float64, copy=True)
    couplingArray = np.array(couplings, dtype=np.float64, copy=True)
    if fieldArray.ndim != 2 or fieldArray.shape[1] != len(alphabet):
        raise ValueError("Fields must have shape (positions, alphabet size)")
    nSites, q = fieldArray.shape
    if couplingArray.shape != (nSites, nSites, q, q):
        raise ValueError("Couplings must have shape (positions, positions, alphabet size, alphabet size)")
    if np.any(np.isnan(fieldArray)) or np.any(np.isposinf(fieldArray)):
        raise ValueError("Only negative infinity fields are allowed")
    if np.any(~np.any(np.isfinite(fieldArray), axis=1)) or not np.all(np.isfinite(couplingArray)):
        raise ValueError("Each site needs finite field support and every coupling must be finite")
    if not np.array_equal(couplingArray, couplingArray.transpose(1, 0, 3, 2)):
        raise ValueError("Couplings must be symmetric under simultaneous site and token exchange")
    if np.any(couplingArray[np.arange(nSites), np.arange(nSites)] != 0):
        raise ValueError("Diagonal coupling blocks must be zero")
    fieldArray.setflags(write=False)
    couplingArray.setflags(write=False)
    return fieldArray, couplingArray


def _validateStates(states, nSites, q):
    stateArray = np.asarray(states)
    if stateArray.ndim < 1 or stateArray.shape[-1] != nSites:
        raise ValueError(f"States must end in a dimension of length {nSites}")
    if stateArray.dtype.kind not in "iu" or np.any(stateArray < 0) or np.any(stateArray >= q):
        raise ValueError("States must contain integer tokens inside the alphabet")
    return stateArray.astype(np.int64, copy=False)


def _selectedEnergy(fields, couplings, states):
    """Gather selected fields so forbidden unused symbols never enter arithmetic."""
    nSites = fields.shape[0]
    outputShape = states.shape[:-1]
    if nSites == 0:
        return np.zeros(outputShape, dtype=np.float64)
    flatStates = states.reshape(-1, nSites)
    pairI, pairJ = np.triu_indices(nSites, 1)
    energies = np.empty(len(flatStates), dtype=np.float64)
    batchSize = max(1, 1_000_000 // max(1, len(pairI)))
    for start in range(0, len(flatStates), batchSize):
        selected = flatStates[start:start + batchSize]
        fieldValues = fields[np.arange(nSites), selected]
        supported = np.all(np.isfinite(fieldValues), axis=1)
        supportedStates = selected[supported]
        pairValues = couplings[pairI, pairJ, supportedStates[:, pairI], supportedStates[:, pairJ]]
        batchEnergy = np.full(len(selected), np.inf, dtype=np.float64)
        try:
            with np.errstate(over="raise", invalid="raise"):
                batchEnergy[supported] = -fieldValues[supported].sum(axis=1) - pairValues.sum(axis=1)
        except FloatingPointError as error:
            raise ValueError("Finite energy terms exceed float64 arithmetic range") from error
        energies[start:start + batchSize] = batchEnergy
    return energies.reshape(outputShape)


@dataclass(frozen=True)
class PottsModel:
    alphabet: str
    fields: np.ndarray
    couplings: np.ndarray
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        fields, couplings = _validateArrays(self.alphabet, self.fields, self.couplings)
        object.__setattr__(self, "fields", fields)
        object.__setattr__(self, "couplings", couplings)
        object.__setattr__(self, "metadata", dict(self.metadata or {}))

    @property
    def n_sites(self):
        return self.fields.shape[0]

    @property
    def q(self):
        return len(self.alphabet)

    def encode(self, sequence: str):
        if not isinstance(sequence, str) or len(sequence) != self.n_sites:
            raise ValueError("Sequence length must match the model")
        if set(sequence) - set(self.alphabet):
            raise ValueError("Sequence contains an unknown symbol")
        return np.array([self.alphabet.index(symbol) for symbol in sequence], dtype=np.int64)

    def decode(self, states):
        stateArray = _validateStates(states, self.n_sites, self.q)
        if stateArray.ndim != 1:
            raise ValueError("Decode accepts exactly one state")
        return "".join(self.alphabet[token] for token in stateArray)

    def energy(self, states):
        stateArray = _validateStates(states, self.n_sites, self.q)
        return _selectedEnergy(self.fields, self.couplings, stateArray)


def _validateDomains(model, domains):
    if not isinstance(domains, Sequence) or len(domains) != model.n_sites:
        raise ValueError("Supply one allowed token collection for every model site")
    result = []
    for site, domain in enumerate(domains):
        try:
            values = tuple(domain)
        except TypeError as error:
            raise ValueError("Each domain must be a collection of integer tokens") from error
        if not values or any(isinstance(token, (bool, np.bool_)) or not isinstance(token, (int, np.integer)) for token in values):
            raise ValueError("Domains must contain nonempty integer token collections")
        if len(set(values)) != len(values) or any(token < 0 or token >= model.q for token in values):
            raise ValueError("Domains cannot repeat tokens or contain unknown tokens")
        if not np.all(np.isfinite(model.fields[site, list(values)])):
            raise ValueError("A requested token is outside the finite original field support")
        result.append(tuple(int(token) for token in values))
    return tuple(result)


@dataclass(frozen=True)
class ReducedModel:
    alphabet: str
    fields: np.ndarray
    couplings: np.ndarray
    domains: tuple[tuple[int, ...], ...]
    original_indices: tuple[int, ...]
    fixed_state: np.ndarray
    offset: float
    source: PottsModel

    def __post_init__(self):
        fields, couplings = _validateArrays(self.alphabet, self.fields, self.couplings)
        object.__setattr__(self, "fields", fields)
        object.__setattr__(self, "couplings", couplings)
        object.__setattr__(self, "domains", _validateDomains(self, self.domains))
        rawFixed = np.asarray(self.fixed_state)
        if rawFixed.dtype.kind not in "iu":
            raise ValueError("Fixed states must contain integer tokens")
        fixedState = np.array(rawFixed, dtype=np.int64, copy=True)
        if self.alphabet != self.source.alphabet or fixedState.shape != (self.source.n_sites,):
            raise ValueError("Reduced model must retain its source alphabet and full fixed state")
        active = tuple(int(site) for site in np.flatnonzero(fixedState == -1))
        if active != self.original_indices or len(active) != self.n_sites:
            raise ValueError("Active source indices disagree with the fixed state")
        fixed = fixedState != -1
        if np.any(fixedState[fixed] < 0) or np.any(fixedState[fixed] >= self.q) or not math.isfinite(self.offset):
            raise ValueError("Fixed tokens and energy offset must be valid and finite")
        if not np.all(np.isfinite(self.source.fields[np.flatnonzero(fixed), fixedState[fixed]])):
            raise ValueError("Fixed tokens must belong to the original field support")
        fixedState.setflags(write=False)
        object.__setattr__(self, "fixed_state", fixedState)

    @property
    def n_sites(self):
        return self.fields.shape[0]

    @property
    def q(self):
        return len(self.alphabet)

    def energy(self, states):
        stateArray = _validateStates(states, self.n_sites, self.q)
        energies = _selectedEnergy(self.fields, self.couplings, stateArray) + self.offset
        allowed = np.ones(stateArray.shape[:-1], dtype=bool)
        for site, domain in enumerate(self.domains):
            allowed &= np.isin(stateArray[..., site], domain)
        return np.where(allowed, energies, np.inf)

    def expand(self, states):
        stateArray = _validateStates(states, self.n_sites, self.q)
        expanded = np.broadcast_to(self.fixed_state, (*stateArray.shape[:-1], self.source.n_sites)).copy()
        expanded[..., list(self.original_indices)] = stateArray
        return expanded

    def random_states(self, rng, chains):
        if isinstance(chains, bool) or not isinstance(chains, (int, np.integer)) or chains < 1:
            raise ValueError("Chain count must be a positive integer")
        states = np.empty((chains, self.n_sites), dtype=np.int64)
        for site, domain in enumerate(self.domains):
            states[:, site] = rng.choice(domain, size=chains)
        return states


def restrict_model(model: PottsModel, domains) -> ReducedModel:
    """Eliminate fixed sites exactly, retaining their full energy contribution."""
    validated = _validateDomains(model, domains)
    active = np.array([site for site, domain in enumerate(validated) if len(domain) > 1], dtype=np.int64)
    fixed = np.array([site for site, domain in enumerate(validated) if len(domain) == 1], dtype=np.int64)
    fixedState = np.full(model.n_sites, -1, dtype=np.int64)
    fixedState[fixed] = [validated[site][0] for site in fixed]
    effectiveFields = model.fields[active].copy()
    for site in fixed:
        effectiveFields += model.couplings[active, site, :, fixedState[site]]
    fixedTokens = fixedState[fixed]
    pairI, pairJ = np.triu_indices(len(fixed), 1)
    offset = -model.fields[fixed, fixedTokens].sum(dtype=np.float64)
    offset -= model.couplings[fixed[pairI], fixed[pairJ], fixedTokens[pairI], fixedTokens[pairJ]].sum(dtype=np.float64)
    return ReducedModel(model.alphabet, effectiveFields, model.couplings[active[:, None], active[None, :]],
                        tuple(validated[site] for site in active), tuple(map(int, active)), fixedState,
                        float(offset), model)


def _checkPin(digest, byteCount, expectedHash, expectedBytes):
    if expectedBytes is not None and (type(expectedBytes) is not int or byteCount != expectedBytes):
        raise ValueError(f"Parameter byte count mismatch: expected {expectedBytes}, received {byteCount}")
    if expectedHash is not None and digest != expectedHash:
        raise ValueError("SHA256 mismatch for pinned input")


def _validatePin(expectedHash, expectedBytes):
    if not isinstance(expectedHash, str) or len(expectedHash) != 64 or set(expectedHash) - set("0123456789abcdef"):
        raise ValueError("A lowercase hexadecimal SHA256 digest is required")
    if type(expectedBytes) is not int or expectedBytes < 1:
        raise ValueError("A positive pinned byte count is required")


def _parseRecord(parts, positions, symbolIndex, lineNumber):
    kind = parts[0]
    if kind not in (b"h", b"J") or len(parts) != (4 if kind == b"h" else 6):
        raise ValueError(f"Malformed parameter record at line {lineNumber}")
    try:
        site = int(parts[1])
        value = float(parts[-1])
        symbol = symbolIndex[parts[2] if kind == b"h" else parts[3]]
        other = None if kind == b"h" else (int(parts[2]), symbolIndex[parts[4]])
    except (ValueError, KeyError) as error:
        raise ValueError(f"Invalid parameter token at line {lineNumber}") from error
    if not 0 <= site < positions:
        raise ValueError(f"Parameter position outside model at line {lineNumber}")
    if other is not None and not site < other[0] < positions:
        raise ValueError(f"Coupling requires ordered distinct sites at line {lineNumber}")
    if not math.isfinite(value) and not (kind == b"h" and parts[-1].lower() in (b"-inf", b"-infinity")):
        raise ValueError(f"Only literal negative infinity fields are admitted at line {lineNumber}")
    return kind, site, symbol, other, value


def parse_parameters(stream: BinaryIO, *, positions: int, alphabet: str,
                     expected_sha256=None, expected_bytes=None) -> PottsModel:
    """Read complete unique h and upper triangle J records without changing values."""
    if type(positions) is not int or positions < 1 or not isinstance(alphabet, str) or not alphabet or len(set(alphabet)) != len(alphabet):
        raise ValueError("Parameter dimensions and alphabet must be explicit and valid")
    q = len(alphabet)
    fields = np.empty((positions, q), dtype=np.float64)
    couplings = np.zeros((positions, positions, q, q), dtype=np.float64)
    seenFields = np.zeros_like(fields, dtype=bool)
    seenPairs = np.zeros((positions * (positions - 1) // 2, q, q), dtype=bool)
    symbolIndex = {symbol.encode("ascii"): token for token, symbol in enumerate(alphabet)}
    digest, byteCount = hashlib.sha256(), 0
    for lineNumber, raw in enumerate(stream, 1):
        raw = raw.encode("ascii") if isinstance(raw, str) else raw
        digest.update(raw)
        byteCount += len(raw)
        parts = raw.split()
        if not parts or parts[0].startswith(b"#"):
            continue
        kind, site, symbol, other, value = _parseRecord(parts, positions, symbolIndex, lineNumber)
        if kind == b"h":
            if seenFields[site, symbol]:
                raise ValueError(f"Duplicate field key at line {lineNumber}")
            seenFields[site, symbol], fields[site, symbol] = True, value
            continue
        second, secondSymbol = other
        pair = site * (2 * positions - site - 1) // 2 + second - site - 1
        if seenPairs[pair, symbol, secondSymbol]:
            raise ValueError(f"Duplicate coupling key at line {lineNumber}")
        seenPairs[pair, symbol, secondSymbol] = True
        couplings[site, second, symbol, secondSymbol] = value
        couplings[second, site, secondSymbol, symbol] = value
    if not np.all(seenFields) or not np.all(seenPairs):
        raise ValueError("Incomplete parameter key coverage")
    _checkPin(digest.hexdigest(), byteCount, expected_sha256, expected_bytes)
    return PottsModel(alphabet, fields, couplings, {"sha256": digest.hexdigest(), "bytes": byteCount,
                     "field_records": int(seenFields.size), "coupling_records": int(seenPairs.size)})


def load_parameters(path, *, positions, alphabet, expected_sha256, expected_bytes):
    """Load a standalone parameter file only with an explicit content pin."""
    _validatePin(expected_sha256, expected_bytes)
    with Path(path).open("rb") as stream:
        return parse_parameters(stream, positions=positions, alphabet=alphabet,
                                expected_sha256=expected_sha256, expected_bytes=expected_bytes)


@dataclass
class BenchmarkInputs:
    model: PottsModel
    anchor: str
    witness: str
    mapping: list[dict]
    cases: list[dict]
    provenance: dict
    config: dict


def _readJson(path):
    def uniqueObject(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=uniqueObject)


def _parseFasta(data):
    records, current = {}, None
    for raw in data.decode("ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:]
            if not current or current in records:
                raise ValueError("FASTA headers must be nonempty and unique")
            records[current] = ""
            continue
        if current is None:
            raise ValueError("FASTA sequence appeared before a header")
        records[current] += line
    if not records or any(not sequence for sequence in records.values()):
        raise ValueError("FASTA contains no records or an empty sequence")
    return records


def _coordinateMap(anchor, reference, path):
    ungapped = anchor.replace("-", "")
    first = reference.find(ungapped)
    if not ungapped or first < 0 or reference.find(ungapped, first + 1) != -1:
        raise ValueError("Anchor must match exactly one contiguous reference segment")
    position, expected = first + 1, []
    for site, symbol in enumerate(anchor):
        expected.append({"model_index_0based": site, "alignment_column_1based": site + 1,
                         "anchor_symbol": symbol, "uniprot_position_1based": None if symbol == "-" else position})
        position += int(symbol != "-")
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    try:
        actual = [{"model_index_0based": int(row["model_index_0based"]),
                   "alignment_column_1based": int(row["alignment_column_1based"]),
                   "anchor_symbol": row["anchor_symbol"],
                   "uniprot_position_1based": int(row["uniprot_position_1based"]) if row["uniprot_position_1based"] else None}
                  for row in rows]
    except (ValueError, KeyError, TypeError) as error:
        raise ValueError("Invalid biological coordinate map") from error
    if actual != expected:
        raise ValueError("Coordinate map disagrees with the exact reference sequence mapping")
    return expected


def _validateContract(config, restrictionData):
    modelConfig = config["model"]
    required = {"model_index_base": 0, "layout": "complete_upper_triangle_J_i_lt_j_and_complete_h",
                "negative_infinity_field_semantics": "zero_probability_site_state",
                "allowed_nonfinite": "h=-inf only; reject NaN,+inf and all nonfinite J",
                "energy": "E(s) = -sum_i h_i(s_i) - sum_{i<j} J_ij(s_i,s_j)",
                "target": "p(s | case) proportional to exp(-beta*E(s)) on the case support",
                "parameter_transformations": []}
    if any(modelConfig.get(key) != value for key, value in required.items()):
        raise ValueError("Unsupported model parameter or energy semantics")
    if type(modelConfig["model_index_base"]) is not int or type(modelConfig["beta"]) not in (int, float) or not math.isfinite(modelConfig["beta"]) or modelConfig["beta"] <= 0:
        raise ValueError("Model index base and positive beta must be explicit numeric values")
    if type(config.get("schema_version")) is not int or config.get("schema_version") != 1 or type(restrictionData.get("schema_version")) is not int or restrictionData.get("schema_version") != 1:
        raise ValueError("Unsupported benchmark schema version")
    if type(restrictionData.get("model_index_base")) is not int or restrictionData.get("model_index_base") != 0 or type(restrictionData.get("reference_position_base")) is not int or restrictionData.get("reference_position_base") != 1 or restrictionData.get("model_parameter_transformations") != []:
        raise ValueError("Unsupported restriction coordinates or parameter transformations")
    if type(config["coordinate_reference"].get("reference_position_base")) is not int or config["coordinate_reference"].get("reference_position_base") != 1:
        raise ValueError("Biological coordinates must use reference position base one")
    archive = config["archive"]
    _validatePin(archive["sha256"], archive["bytes"])
    _validatePin(archive["member_sha256"], archive["member_bytes"])
    if not isinstance(restrictionData.get("cases"), list) or not restrictionData["cases"]:
        raise ValueError("Benchmark must declare a nonempty restriction case list")


def _commonDomains(inputs):
    constraints = inputs.provenance["restriction_data"]["common_constraints"]
    expectedKeys = {"fixed_gap_model_indices", "fixed_gap_symbol", "nongap_columns_exclude_gap", "finite_original_field_support_required", "role"}
    if set(constraints) - expectedKeys or constraints.get("fixed_gap_symbol") != "-" or constraints.get("nongap_columns_exclude_gap") is not True or constraints.get("finite_original_field_support_required") is not True:
        raise ValueError("Unsupported common restriction semantics")
    gaps = [site for site, symbol in enumerate(inputs.anchor) if symbol == "-"]
    declared = constraints.get("fixed_gap_model_indices")
    configured = inputs.config["coordinate_reference"]["gap_model_indices"]
    if not isinstance(declared, list) or any(type(site) is not int for site in declared) or declared != gaps or configured != gaps:
        raise ValueError("Common fixed gaps disagree with biological coordinates")
    model = inputs.model
    domains = [tuple(token for token, symbol in enumerate(model.alphabet)
                     if np.isfinite(model.fields[site, token]) and ((symbol == "-") == (site in gaps)))
               for site in range(model.n_sites)]
    if any(not domain for domain in domains):
        raise ValueError("Common restrictions leave an empty site domain")
    return domains


def benchmark_domains(inputs: BenchmarkInputs, case_id: str):
    if not isinstance(case_id, str):
        raise ValueError("Case identifier must be a string")
    if not isinstance(inputs.cases, list) or any(not isinstance(case, dict) or not isinstance(case.get("case_id"), str) for case in inputs.cases):
        raise ValueError("Cases must be objects with string identifiers")
    identifiers = [case.get("case_id") for case in inputs.cases]
    if len(set(identifiers)) != len(identifiers) or case_id not in identifiers:
        raise ValueError("Unknown or duplicate restriction case identifier")
    case = inputs.cases[identifiers.index(case_id)]
    common = _commonDomains(inputs)
    semantics = {"exact_conditional": "fixed to computational_witness",
                 "broad_conditional": "common gap pattern and finite original field support"}
    if case.get("kind") not in semantics or case.get("other_positions") != semantics[case["kind"]]:
        raise ValueError("Unsupported restriction case semantics")
    witness = inputs.model.encode(inputs.witness)
    if any(token not in common[site] for site, token in enumerate(witness)):
        raise ValueError("Computational witness violates common restrictions or original support")
    domains = [(int(token),) for token in witness] if case["kind"] == "exact_conditional" else common.copy()
    changed = set()
    restrictions = case.get("positional_restrictions")
    if not isinstance(restrictions, list):
        raise ValueError("Positional restrictions must be a list")
    for restriction in restrictions:
        if not isinstance(restriction, dict) or set(restriction) - {"model_index_0based", "allowed_symbols", "uniprot_position_1based"}:
            raise ValueError("Unsupported positional restriction semantics")
        site = restriction.get("model_index_0based")
        symbols = restriction.get("allowed_symbols")
        if type(site) is not int or not 0 <= site < inputs.model.n_sites or site in changed:
            raise ValueError("Invalid or duplicate restriction site")
        changed.add(site)
        if "uniprot_position_1based" in restriction and (type(restriction["uniprot_position_1based"]) is not int or restriction["uniprot_position_1based"] != inputs.mapping[site]["uniprot_position_1based"]):
            raise ValueError("Requested biological coordinate does not match the model site")
        if not isinstance(symbols, str) or not symbols or len(set(symbols)) != len(symbols) or set(symbols) - set(inputs.model.alphabet):
            raise ValueError("Allowed symbols must be explicit, unique and known")
        requested = tuple(inputs.model.alphabet.index(symbol) for symbol in symbols)
        if any(token not in common[site] for token in requested):
            raise ValueError("Requested state is forbidden or violates the fixed gap pattern")
        domains[site] = requested
    return tuple(domains)


def _validateSupport(model, supportData):
    positions = supportData.get("positions")
    if not isinstance(positions, list) or len(positions) != model.n_sites:
        raise ValueError("Support evidence has the wrong number of positions")
    for site, record in enumerate(positions):
        allowed = "".join(symbol for token, symbol in enumerate(model.alphabet) if np.isfinite(model.fields[site, token]))
        forbidden = "".join(symbol for symbol in model.alphabet if symbol not in allowed)
        if record != {"model_index_0based": site, "allowed_symbols": allowed, "forbidden_symbols": forbidden}:
            raise ValueError("Frozen support evidence disagrees with original parameters")


def load_benchmark(config_path) -> BenchmarkInputs:
    """Load local pinned bytes and validate coordinates, support and every case."""
    configPath = Path(config_path).resolve()
    root = configPath.parent.parent
    config = _readJson(configPath)
    restrictionData = _readJson(root / config["restriction_cases"])
    _validateContract(config, restrictionData)
    referenceConfig, archiveConfig = config["coordinate_reference"], config["archive"]
    reference = _readJson(root / referenceConfig["metadata"])
    archivePath = root / archiveConfig["local_path"]
    if not archivePath.is_file():
        raise FileNotFoundError(f"Pinned archive is missing: {archivePath}. Acquire it explicitly before loading.")
    archiveData = archivePath.read_bytes()
    _checkPin(hashlib.sha256(archiveData).hexdigest(), len(archiveData), archiveConfig["sha256"], archiveConfig["bytes"])
    with zipfile.ZipFile(io.BytesIO(archiveData)) as archive:
        names = archive.namelist()
        required = {archiveConfig["member"], referenceConfig["archive_member"], config["finite_score_fixture"]["archive_member"]}
        if any(names.count(name) != 1 for name in required):
            raise ValueError("Required archive members must be present exactly once")
        with archive.open(archiveConfig["member"]) as stream:
            model = parse_parameters(stream, positions=config["model"]["positions"], alphabet=config["model"]["alphabet"],
                                     expected_sha256=archiveConfig["member_sha256"], expected_bytes=archiveConfig["member_bytes"])
        natural = _parseFasta(archive.read(referenceConfig["archive_member"]))
        anchor = natural[referenceConfig["archive_header"]]
        finiteRecords = _parseFasta(archive.read(config["finite_score_fixture"]["archive_member"]))
        fixture = finiteRecords[config["finite_score_fixture"]["archive_header"]]
    model.encode(anchor)
    if anchor != reference["aligned_sequence"] or reference["uniprot_accession"] != referenceConfig["uniprot_accession"] or reference["organism"]["taxonId"] != referenceConfig["taxon_id"]:
        raise ValueError("Coordinate anchor or accession provenance disagrees with the frozen reference")
    mapping = _coordinateMap(anchor, reference["full_accession_sequence"], root / referenceConfig["map"])
    mapped = [row["uniprot_position_1based"] for row in mapping if row["uniprot_position_1based"] is not None]
    if [min(mapped), max(mapped)] != referenceConfig["reference_span_in_full_accession"]:
        raise ValueError("Mapped reference span disagrees with the benchmark contract")
    if not np.isfinite(model.energy(model.encode(fixture))):
        raise ValueError("Declared finite score fixture is outside original model support")
    _validateSupport(model, _readJson(root / config["support"]))
    provenance = {"archive_sha256": archiveConfig["sha256"], "parameter_sha256": model.metadata["sha256"],
                  "config_sha256": hashlib.sha256(configPath.read_bytes()).hexdigest(), "reference": reference,
                  "restriction_data": restrictionData, "finite_fixture": fixture}
    inputs = BenchmarkInputs(model, anchor, restrictionData["computational_witness"]["sequence"],
                             mapping, restrictionData["cases"], provenance, config)
    for case in inputs.cases:
        benchmark_domains(inputs, case["case_id"])
    return inputs
