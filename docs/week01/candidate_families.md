# Week 1 candidate families and library-generation research

Research note and scope draft; sources accessed **2026-10-08**. The findings below are a literature/resource review and a read-only archive audit. They are not sampling results, experimental validation, or a final family selection.

## Scope draft

Compare a direct categorical sampler with a corrected encoded method, with a planner choosing computational configurations, for **one fixed, previously trained Potts distribution**. Select one nonpathogenic reference/family context before the biological benchmark. Preserve the supplied fields, couplings, alphabet, support, and temperature across methods. Measure computational efficiency and distributional agreement; existing assays may provide limited retrospective context. No new model training or wet-lab work is included.

The provisional first choice is chorismate mutase (CM), because a public archive contains an actual trained parameter file and the family has a direct statistical-generation precedent. This is conditional: the archived model has negative-infinity fields, generic sequence identifiers, and experimental-feedback provenance. Its support semantics and reference-position mapping remain unresolved. **None of the three candidates is benchmark-ready yet.**

## How researchers construct candidate libraries

Three relevant routes have different meanings for this project:

1. **Family-model generation.** Russ et al. inferred a pairwise statistical model from an aligned CM family and generated artificial sequences by Monte Carlo sampling. Their study compared conserved-site-only and correlated models and measured activity after constructing libraries. It also used temperatures below one and a phenotype classifier to focus candidates. These choices motivate studying both correlations and assay context; changing temperature or phenotype-filtering samples changes the distribution being evaluated. Our comparison must instead keep its chosen target fixed. [Russ et al., Science (2020), author-hosted paper](https://www.phys.ens.fr/~monasson/Articles/a125.pdf).
2. **Local mutation libraries.** Sarkisyan et al. generated a local avGFP landscape through random mutagenesis and measured fluorescence. Its labels describe sampled variants near a particular reference; they do not constitute a pretrained generative model or label arbitrary distant samples. [Sarkisyan et al., Nature (2016), author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC4968632/).
3. **Systematic substitutions.** Roscoe et al. measured ubiquitin point-mutant effects on yeast growth; Mavor et al. examined environmental dependence using chemical perturbations. These are useful assay contexts, but single-mutant labels do not establish whole-family generation quality. [Roscoe et al., JMB (2013)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3615125/), [Mavor et al., eLife (2016)](https://elifesciences.org/articles/15802).

For the computational study, the candidate-library workflow is: acquire a fixed model and its coordinate metadata; generate samples under that model; retain the full sample record for distributional analysis; then optionally form a separately labelled retrospective subset that matches existing assays. Deduplication, selecting only low-energy sequences, or enforcing extra biological constraints must not silently replace the raw sampler output. A planner can select a backend, proposal, or run budget, but correctness requires the same target and appropriate correction for the proposal/encoding.

## Comparison of three candidates

ProteinGym separates processed assay tables, MSAs, and model prediction scores into different downloads. Its reference metadata gives assay target sequences and MSA spans; neither a score table nor an alignment establishes that trained Potts parameters are downloadable. Counts below are metadata observations, not newly measured quantities. [ProteinGym resource inventory](https://github.com/OATML-Markslab/ProteinGym#resources), [reference-field definitions](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/reference_files_description.md).

| Candidate | Assay and coordinate evidence | Alignment versus trained model evidence | Fit and remaining gaps |
| --- | --- | --- | --- |
| **Chorismate mutase, AroQ context** | Russ provides the relevant generation/activity precedent. The later CM archive contains natural and artificial alignments plus experimental labels. Its `Natural_CM.faa` has 1,130 records of aligned length 96; the inspected first header is a generic sequence identifier with a functional label, not an organism identifier. [Russ paper](https://www.phys.ens.fr/~monasson/Articles/a125.pdf), [pinned CM archive](https://github.com/FrancescoCalvanese/DCA_Reintegration/blob/433fda34132986d64db564bda867abb9ec9d46f2/DATA/Chorismate_Mutase.zip). | **Trained parameters verified present:** `CM_model_lambda_1.dat`, a feedback-reintegrated model. The same archive has sampled sequences for lambda zero, but no lambda-zero CM parameter member. This is a later model, not a verified copy of Russ's original checkpoint. [release record](https://zenodo.org/records/15115193), [Calvanese et al., NAR (2025)](https://doi.org/10.1093/nar/gkaf832). | Strongest acquisition lead and compact aligned representation. The presence of 362 `h=-inf` entries, missing biological position map, and unresolved nonpathogenic reference provenance prevent benchmark admission. Existing feedback labels overlap model development. |
| **Aequorea victoria GFP** | `GFP_AEQVI_Sarkisyan_2016`: target length 238; mutated region 3–237; fluorescence assay; 51,714 ProteinGym rows. [pinned assay metadata](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/DMS_substitutions.csv). | Metadata identifies `GFP_AEQVI_full_04-29-2022_b08.a2m`, spanning 1–238, with 396 sequences and reported effective depth 13.1. **No trained GFP Potts checkpoint verified** in the resources inspected. [same metadata](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/DMS_substitutions.csv). | Reporter context and useful local multi-mutant labels. Larger model than the archived CM representation. Low reported effective alignment depth is context, not a measured failure. A fixed checkpoint, its training alignment, and exact assay/model coordinate agreement are acquisition requirements. |
| **Yeast ubiquitin** | `RL40A_YEAST_Roscoe_2013`: 1,195 rows; growth assay; mutated region 2–76. ProteinGym's target is length 128, so the assay's ubiquitin region must be distinguished from the full target. [pinned metadata](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/DMS_substitutions.csv). | Metadata names `RL40A_YEAST_full_11-26-2021_b01.a2m`, span 1–128, 16,228 sequences, effective depth 3,906.9. **No trained ubiquitin Potts checkpoint verified.** [same metadata](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/DMS_substitutions.csv). | Benign housekeeping-protein alternative with a short assayed region. A supplied 128-position checkpoint cannot be cropped to ubiquitin by simply dropping interactions: that changes the target. Need a compatible existing domain model or an explicitly defined fixed-context distribution, plus a coordinate map. |

These candidates concern a metabolic enzyme, a fluorescent reporter, and a yeast housekeeping protein. This does **not** certify every homolog in a family alignment as nonpathogenic. In particular, the CM archive's generic headers do not establish taxonomic provenance. The admission decision must identify a specific reference and document the intended biological context. GB1 was considered but is not on the shortlist: its IgG-binding assay adds an unnecessary suitability question for the requested nonpathogenic scope. [Wu et al., eLife (2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4985287/).

## Verified CM acquisition lead and material obstacle

The March 2025 `v1.0.0` release is archived at [Zenodo record 15115193](https://zenodo.org/records/15115193). The inspected GitHub tag resolves to commit `433fda34132986d64db564bda867abb9ec9d46f2`. Its [notebook](https://github.com/FrancescoCalvanese/DCA_Reintegration/blob/433fda34132986d64db564bda867abb9ec9d46f2/Notebook_Implementation/Notebook_implementation.ipynb) identifies the labelled artificial CM sequences as reintegration data sourced from Russ. The publication describes updating the model using experimental feedback; its generated-CM functionality assessment uses a classifier rather than new experimental tests of those generated sequences. [Calvanese et al. (2025)](https://doi.org/10.1093/nar/gkaf832).

Read-only audit of the pinned `DATA/Chorismate_Mutase.zip` found:

| Artifact observation | Value |
| --- | --- |
| Compressed archive size | 13,145,579 bytes |
| Archive SHA-256 | `193c7dc912b394df12f0662fc6fde70d265567e7c9265b5106c3cec3d5746329` |
| Parameter member | `Chorismate_Mutase/CM_model_lambda_1.dat` |
| Uncompressed parameter size | 44,044,380 bytes |
| Parameter-member SHA-256 | `bf180a3ed43637ab082a0d398d6663ffbd2566e73b2232404519ac8e4f32a004` |
| Parameter records | 2,010,960 `J`; 2,016 `h` |
| Position/symbol observations | 96 unique indices, 0–95; 21 symbols, `-ACDEFGHIKLMNPQRSTVWY` |
| Non-finite classification | 362 negative-infinity fields (`h=-inf`), across 92 positions; no non-finite `J`, NaN, or positive infinity observed |
| Other relevant members | `Natural_CM.faa`, `P1_sample_CM.faa`, `P1_sample_exp.txt`, `CM_sequences_lambda_0.fasta`, `CM_sequences_lambda_1.fasta` |

These are file-content observations from the [pinned archive](https://github.com/FrancescoCalvanese/DCA_Reintegration/blob/433fda34132986d64db564bda867abb9ec9d46f2/DATA/Chorismate_Mutase.zip), not a claim that the parameters are validated. The counts are consistent with a dense 96-position, 21-state upper-triangle representation; duplicate-key, complete-coverage, and loader checks remain to be performed.

An independently repeated read-only audit confirmed the archive hash and counts. Its machine-readable record is [cm_parameter_audit.json](cm_parameter_audit.json), including the parameter-member hash. Both inspections used the unmodified archive in memory.

Negative-infinity fields may encode prohibited site states, giving them zero probability. That interpretation must be verified against the model's construction and intended support. **Do not replace them with finite cutoffs, delete affected symbols, or retrain:** those actions can change the target. There is also an implementation issue to resolve: the inspected adabmDCA energy routine uses a one-hot matrix product with the fields. Under ordinary floating-point arithmetic, zero times negative infinity yields NaN. It therefore needs an explicit compatibility/support check before this model is scored or sampled. This is a code-based inference, not an observed end-to-end failure. [Pinned `compute_energy` source](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/statmech.py).

The lambda-one model may still be an appropriate **fixed computational target** if its support can be represented exactly and its provenance fits the scope. It is unsuitable for claiming independent functional prediction on the labels used for reintegration. No independent held-out assay split has been established.

### Reproduce the archive audit

This Python standard-library script downloads into memory and inspects the archive without executing upstream code or modifying parameters. The archive was not retained locally during this review. It reproduces the primary counts and non-finite classification; it is not a full model validator.

```python
import collections
import hashlib
import io
import json
import math
import urllib.request
import zipfile

url = (
    "https://raw.githubusercontent.com/FrancescoCalvanese/"
    "DCA_Reintegration/433fda34132986d64db564bda867abb9ec9d46f2/"
    "DATA/Chorismate_Mutase.zip"
)
archive = urllib.request.urlopen(url, timeout=30).read()
records = collections.Counter()
nonfinite = collections.Counter()
positions, symbols, nonfinite_field_positions = set(), set(), set()
member = "Chorismate_Mutase/CM_model_lambda_1.dat"
with zipfile.ZipFile(io.BytesIO(archive)) as z:
    members = {entry.filename: entry.file_size for entry in z.infolist()}
    with z.open(member) as stream:
        for raw in stream:
            parts = raw.decode("ascii").split()
            if not parts:
                continue
            kind = parts[0]
            records[kind] += 1
            if kind == "J":
                positions.update((int(parts[1]), int(parts[2])))
                symbols.update((parts[3], parts[4]))
            elif kind == "h":
                positions.add(int(parts[1]))
                symbols.add(parts[2])
            if not math.isfinite(float(parts[-1])):
                nonfinite[f"{kind}:{parts[-1]}"] += 1
                if kind == "h":
                    nonfinite_field_positions.add(int(parts[1]))
print(json.dumps({
    "archive_bytes": len(archive),
    "archive_sha256": hashlib.sha256(archive).hexdigest(),
    "members": members,
    "records": dict(records),
    "position_count": len(positions),
    "position_range": [min(positions), max(positions)],
    "symbols": "".join(sorted(symbols)),
    "nonfinite": dict(nonfinite),
    "nonfinite_field_position_count": len(nonfinite_field_positions),
}, indent=2))
```

## Beginning adabmDCA model reading

The current Python implementation documents explicit loading of existing `params.dat` files, sequence scoring, and Metropolis/Gibbs generation. This provides a reading route without running the training workflow. Its public example directory at the inspected commit contains the RNA alignment `RF00379.fasta`; the older C++ repository's test directory contains protein alignments. These examples do not supply a verified pretrained protein benchmark. [High-level API](https://spqb.github.io/adabmDCApy/high_level_api/), [pinned Python example tree](https://github.com/spqb/adabmDCApy/tree/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/example_data), [older test tree](https://github.com/anna-pa-m/adabmDCA/tree/0a308aaae13e38a8c444bfefef3c526b8293d74a/test).

The initial reading checklist is:

- **Alphabet and positions:** the documented protein alphabet includes a gap state. Establish the original alphabet order and retained-column map; an MSA span is not sufficient when columns have been filtered. [Preprocessing documentation](https://spqb.github.io/adabmDCApy/preprocessing/), [ProteinGym mapping definitions](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/reference_files_description.md).
- **Parameter representation:** `load_params` reads `J`/`h` text into bias shape `(L, q)` and couplings shape `(L, q, L, q)` and supports one or both triangles. Check symmetrization, diagonal handling, index conventions, and forbidden states before conversion. [Import/export API](https://spqb.github.io/adabmDCApy/api/adabmDCA.io/).
- **Energy convention:** verify the negative field term and the half factor for symmetric couplings against a categorical sum over `i < j`. Fix `p(s) proportional to exp(-beta * E(s))` with the same beta across methods. The source has been read; agreement on a real model has not been tested. [Pinned energy implementation](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/statmech.py), [sampling API](https://spqb.github.io/adabmDCApy/api/adabmDCA.sampling/).
- **Generation and precision:** determine sweep semantics, proposal probabilities, seeds, initialization, burn-in, and mixing diagnostics. Default to documented floating-point precision until numerical equivalence is established; the API also offers BF16 coupling copies, which introduce a separate accuracy question. [High-level generation documentation](https://spqb.github.io/adabmDCApy/high_level_api/#generate-sequences).

## Provisional decision and follow-ups

**Prioritize the archived CM model for an acquisition/compatibility assessment.** GFP is the reporter fallback; yeast ubiquitin is the short-domain fallback. Neither fallback currently meets the fixed-pretrained-model requirement. Absence of a verified checkpoint in this review is not proof that no public checkpoint exists.

For Week 2, begin recording the computational question and admission criteria: one immutable checkpoint and checksum; defined alphabet/support; an identified nonpathogenic reference context; a model-column-to-reference-position map; and a correction method that preserves the target after encoding. Specify that the planner changes execution choices, not model parameters or the target distribution. Investigate whether the feedback-informed CM checkpoint fits the intended study; if a natural-only model is required, the lambda-one archive is insufficient. Final benchmark agreement remains scheduled for Week 4.

By the Week 4 acquisition/mapping decision, verify the CM negative-infinity semantics first, then check parameter key coverage, energy equivalence on supported states, and loader behavior without changing support. Obtain the actual reference and mapping metadata; generic headers and matching sequence lengths cannot establish it. Keep Russ/reintegration labels marked as development-overlapping unless an independent split is documented. For GFP/ubiquitin, acquire an existing compatible checkpoint before biological benchmark integration proceeds; creating one from the available MSA would be new training and outside this scope.

After admission, validate encoding/correction on small enumeratable models before making real-family efficiency claims. Agree in advance on runtime/resource accounting, acceptance or valid-state rates, autocorrelation/effective sample size, and distributional diagnostics. Existing DMS labels can only support claims for mapped assayed variants. No timing, mixing advantage, functional improvement, or planner benefit has been measured in Week 1.
