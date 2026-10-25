# Week 4 model acquisition and admission evidence

**Prepared 2026-10-08.** The pinned CM lambda-one checkpoint now has a concrete
technical benchmark proposal: immutable parameters, exact supported states,
one identified HB8 coordinate reference, and four nonempty restriction cases.
Human agreement, two independent lead checks and adviser review remain pending.
The production loader is Week 5 work. No model training or sampler run was done.
See [model_admission.json](model_admission.json) for those separate statuses and
[benchmark.json](../../configs/benchmark.json) for the input contract.

## Acquisition and model coverage

The existing checkpoint is the feedback-informed `CM_model_lambda_1.dat` in the
[CM archive at commit 433fda34132986d64db564bda867abb9ec9d46f2](https://github.com/FrancescoCalvanese/DCA_Reintegration/blob/433fda34132986d64db564bda867abb9ec9d46f2/DATA/Chorismate_Mutase.zip).
It is not a verified copy of Russ's original natural-only checkpoint. The
[upstream notebook](https://github.com/FrancescoCalvanese/DCA_Reintegration/blob/433fda34132986d64db564bda867abb9ec9d46f2/Notebook_Implementation/Notebook_implementation.ipynb)
sets reintegration strength to one and incorporates experimentally labelled P1
designs. It installs adabmDCA without pinning its version; the precise training
software environment and a training log for this parameter file are absent.

The inspection reproduced the Week 1 archive and member hashes:

| Input | Bytes | SHA-256 |
| --- | ---: | --- |
| ZIP | 13,145,579 | `193c7dc912b394df12f0662fc6fde70d265567e7c9265b5106c3cec3d5746329` |
| Parameter member | 44,044,380 | `bf180a3ed43637ab082a0d398d6663ffbd2566e73b2232404519ac8e4f32a004` |

[inspect_benchmark.py](../../scripts/inspect_benchmark.py) checked every expected
key once: 2,010,960 upper-triangle `J` records and 2,016 `h` records for 96
positions and alphabet `-ACDEFGHIKLMNPQRSTVWY`. All couplings are finite. Exactly
362 fields at 92 positions are literal negative infinity; no NaN, positive
infinity, missing keys or duplicate keys passed validation. The full results
are in [inspection.json](../../data/benchmark/inspection.json). No finite cutoff,
state merging, field shift, coupling removal or other parameter change is used.

## Negative-infinity semantics and toolkit compatibility

The `get_effective_frequencies` function inside the pinned
[Notebook_DATA.zip](https://github.com/FrancescoCalvanese/DCA_Reintegration/blob/433fda34132986d64db564bda867abb9ec9d46f2/Notebook_Implementation/Notebook_DATA.zip)
uses zero pseudocounts and clamps effective frequencies at zero. The notebook
initializes the model from those frequencies. The inspected
[adabmDCA initializer](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/utils.py#L109)
uses their logarithms as fields. This gives a source-based explanation for
negative-infinity fields, without claiming to reconstruct the training run.

With the documented energy convention, selecting such a field gives `E=+inf`
and probability zero at any finite positive beta. Leaving it unselected must
not contaminate a supported state's energy. Every site's allowed set is
recorded explicitly in [support.json](../../data/benchmark/support.json).
Finite selected fields plus finite couplings give finite energy; every admitted
case has at least one such state. A finite Cartesian support therefore has a
positive, finite normalizing constant. These are mathematical conclusions from
the archived parameters and the fixed target definition.

Current [adabmDCA loading code](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/io.py#L195)
rejects nonfinite fields. Its [one-hot energy implementation](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/statmech.py#L7)
also multiplies unselected fields by zero, which can yield NaN for zero times
negative infinity. The Week 2 installed version 1.0.0 has the same compatibility
constraint. The independent admission inspector gathers selected fields;
it does not establish a working production-loader or sampler path.

## Identified reference and numbering

The archive's `sequence_2 functional_false` is a 96-column aligned sequence.
Its 89 nongap letters match **exactly once** to residues 2–90 of the 355-residue
UniProt accession **Q5SLA5**, sequence version 1, entry version 119. The retrieved
[UniProt entry](https://rest.uniprot.org/uniprotkb/Q5SLA5.json) identifies the
organism as *Thermus thermophilus* HB8 / DSM 579 / ATCC 27634, taxon 300852.
The frozen relevant fields and complete accession sequence are retained in
[reference_provenance.json](../../data/benchmark/reference_provenance.json).
[PDB 2D8D](https://www.rcsb.org/structure/2D8D) corroborates this accession's
90-residue chorismate-mutase domain context. The model is not a 355-site full
protein model, and the missing initiator methionine is not renumbered into it.

[DSMZ's HB8 record](https://www.dsmz.de/collection/catalogue/details/culture/dsm-579)
lists DSM 579 as risk group 1 under German TRBA. This supports the selected
nonpathogenic reference context. It does not classify all homologs in the
training alignment or every generated sequence.

[position_map.csv](../../data/benchmark/position_map.csv) retains all 96 model
columns. Model indices are zero based; alignment columns and full-accession
residue positions are one based. Gap columns 0, 1, 2, 63, 64, 94, 95 have a null
biological residue position. Model index 3 maps to Q5SLA5 residue 2; index 93 maps
to residue 90. This is an explicitly recorded sequence-based coordinate
inference, rather than an upstream supplied map or a claim that generic FASTA
headers prove taxonomy. A request for a biological restriction at an unmapped
gap must be rejected; computational gap constraints name model columns.

The native coordinate anchor is outside lambda-one support: model 38 contains
Q at biological residue 37, and model 45 contains H at biological residue 44;
both fields are negative infinity. Its model energy is **positive infinity**. It can
anchor numbering and Hamming distance, but cannot initialize a supported chain.

## Finite arithmetic fixture and computational witness

The unmodified archive `Natural_CM.faa` record `sequence_1 functional_true`
serves only as a finite computational arithmetic fixture. No nonpathogenic
taxonomic assertion is made for that record. Summing the original 96 selected
fields and 4,560 selected upper-triangle couplings gives:

```text
sum h = -45.43422150
sum J(i<j) = 102.792881272664
E = -(-45.43422150) - 102.792881272664
  = -57.358659772664
```

Exact decimal addition and binary64 `fsum` agree within 1e-10. Those are two
arithmetic paths in one inspection, not two independent project-lead checks.
This fixture does not have the HB8 gap pattern and cannot initialize the
restricted cases below. [score_checks.json](../../data/benchmark/score_checks.json)
preserves the sequence, sums, selected-term counts and first five pair terms.

For nonempty-support and exact-case checks, the separately named computational
witness retains the coordinate-anchor letters except 38 Q→E and 45 H→A. Those
replacement symbols are the first permitted amino acids in the fixed alphabet.
Its finite full-model energy is 135.208916945824. It is a constructed validity
witness, not a native HB8 protein, assayed variant or biological recommendation.
No parameter is changed to create it. Future broadly supported chains require
their own documented valid initialization; this witness is not mixing evidence.

## Frozen proposed restriction inputs

All four cases fix the seven reference gap columns to `-`, prohibit gap at every
other column, and preserve the original finite-field support. The case file
[restriction_cases.json](../../data/benchmark/restriction_cases.json) records
zero-based model coordinates; the corresponding Q5SLA5 residues are 8→7,
12→11, 16→15 and 24→23. The selected development and evaluation positions are
disjoint. These are computational restriction-input splits, not assay holdouts.
Evaluation cases must not be used to fit the planner.

| Case | Split | Explicit selected domains | Remaining positions |
| --- | --- | --- | --- |
| `cm_dev_exact_8_12` | Development | 8:{A,C};12:{E,F} | Fixed computational witness |
| `cm_eval_exact_16_24` | Evaluation | 16:{V,I};24:{L,M} | Fixed computational witness |
| `cm_dev_broad_8_12` | Development | 8:{A};12:{E,F} | Full finite-field domains and common gap pattern |
| `cm_eval_broad_16_24` | Evaluation | 16:{V};24:{L,M} | Full finite-field domains and common gap pattern |

Each exact case has 4 states. Each broad case has 88 mutable columns and more
than 1e99 supported states. Every requested symbol and fixed witness symbol was
checked against original field support. The exact conditional energies retain
all 96 fields and 4,560 pairs, including all mutable-background and
background-background terms. They do not crop the model to two sites. Their
beta 1 conditional probabilities are retained in the score artifact; these are
arithmetic references, not measured sampler frequencies. Selecting these
sites defines transparent stress cases and implies no functional rationale.
Human approval of the restrictions and feedback-informed checkpoint remains
pending.

## Experimental provenance and next checks

P1 measurements informed reintegration, so they cannot serve as an independent
test of lambda-one biological prediction. The upstream notebook's generated
sequence classification is a trained classifier result. No independent mapped
assay holdout was established. The sampled lambda-one FASTA also has
`log_weight=nan` headers; those archived outputs are not accepted as validated
equilibrium chains or an accuracy baseline. Experimental evaluation is disabled
in this proposal. There is no need to substitute a new GFP/ubiquitin model after
CM's exact support and reference coordinates were established.

The remaining work before real-model sampling is: human scope/restriction
agreement, independent score/numbering checks, the adviser review, and the
Week 5 support-aware production loader with its own tests. Runtime, mixing,
functional improvement and planner benefit have not been measured.

## Reproduction

From the project root, with Python 3.13.5 and only the standard library:

```powershell
python scripts/inspect_benchmark.py --download
python -m unittest discover -s tests -p test_benchmark_admission.py -v
```

The downloader accepts only the pinned archive checksum and byte count. An
existing corrupted local archive fails rather than being silently replaced.
The parameter member is streamed inside the ignored ZIP. The inspector
reproduces support, mapping, all selected full-model scores and exact
conditional probabilities from the locally frozen reference metadata. It
does not reretrieve mutable UniProt/DSMZ pages; their snapshot hashes, versions,
source URLs and relevant metadata are retained in the provenance files.

Ten admission tests pass locally. They cover corrupted inputs, missing and
duplicate keys, invalid positions/symbols, all nonfinite values except literal
negative-infinity fields, overlapping reference matches, gap/reference mapping,
conflicting restrictions, temperature semantics, and the four frozen cases'
full-model evidence. [source_manifest.json](../../data/benchmark/source_manifest.json)
records downloaded primary-source fingerprints; this acquisition exercise does
not complete the Week 5 production-loader milestone.
