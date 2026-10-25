# Week 4 benchmark contract

**Status: technically specified proposal, pending both-lead and adviser review.** The acquisition inspector checks this model and the small conditional cases; the production loader and samplers remain future milestones. The [admission record](model_admission.json) and [source evidence](benchmark_evidence.md) distinguish verified facts from outstanding review.

## Fixed model and coordinates

The proposed benchmark is `cm_lambda1_supported_v1`. [benchmark.json](../../configs/benchmark.json) pins the CM lambda-one checkpoint from DCA_Reintegration commit `433fda34132986d64db564bda867abb9ec9d46f2`, with archive and parameter-member SHA256. The model has 96 aligned positions, alphabet `-ACDEFGHIKLMNPQRSTVWY`, complete single-site fields and canonical `i < j` couplings. No parameters are trained, clipped, replaced or otherwise changed.

Use `E(s) = -sum_i h_i(s_i) - sum_(i<j) J_ij(s_i,s_j)`, with `beta = 1/T = 1` and binary64 for the comparison. The restricted target is proportional to `exp(-E(s))` on each declared case. An original `h=-inf` excludes that site state exactly; a sequence selecting it has `E=+inf` and zero probability. All couplings are finite. The current toolkit's finite-only loader is incompatible with these fields, so Week 5 must use explicit support and selected-state scoring without `0 * infinity` arithmetic.

The biological coordinate anchor is the unmodified archive record `sequence_2 functional_false`. Its 89 nongap residues uniquely match UniProt Q5SLA5 residues 2–90, a segment of the 355-residue protein from *Thermus thermophilus* HB8. The strain context is documented in the source evidence; it does not classify every family member or generated state. The [position map](../../data/benchmark/position_map.csv) records zero-based model indices, one-based alignment columns and one-based full-accession residue positions. Gaps have no biological residue number.

This anchor lies outside model support at model positions 38 (`Q`) and 45 (`H`). Preserve its infinite energy and use it only for coordinates and Hamming distance. It is not a valid sampler initialization. The separate computational witness changes those symbols to `E` and `A`; this is an explicitly constructed supported state, not a native or experimentally tested protein. The unmodified archive `sequence_1 functional_true` is a separate finite arithmetic fixture; its taxonomy is not asserted.

## Positional cases and split

All cases fix model positions `0, 1, 2, 63, 64, 94, 95` to gap and exclude gap at every other position. Allowed choices must also have finite original fields. Reject an explicit forbidden requested residue or an empty domain; do not silently repair the restriction. Hamming distance is measured against the anchor over the 89 nongap model columns.

The authoritative domains are in [restriction_cases.json](../../data/benchmark/restriction_cases.json), together with the complete computational witness.

| Case | Selected restrictions using model indices | Other positions | Role |
| --- | --- | --- | --- |
| `cm_dev_exact_8_12` | 8 allows A/C; 12 allows E/F | Fixed to witness | Four-state development reference |
| `cm_dev_broad_8_12` | 8 fixed A; 12 allows E/F | Common gap pattern and finite field support | Development timing/quality input |
| `cm_eval_exact_16_24` | 16 allows V/I; 24 allows L/M | Fixed to witness | Four-state evaluation reference |
| `cm_eval_broad_16_24` | 16 fixed V; 24 allows L/M | Common gap pattern and finite field support | Evaluation timing/quality input |

These selected model indices map to full-accession residues 7, 11, 15 and 23, respectively. The four-state energies include **all 96 fields and all 4,560 selected pair terms**. Fixing a background does not permit dropping its interactions with variable sites.

This is a split of computational restriction inputs, not a held-out training alignment or independent assay split. Both groups share the immutable model and common constraints. Use only development cases to choose planner settings. Evaluation reference probabilities can be generated for correctness checking; no evaluation sampling/timing result may guide tuning. Existing reintegration labels are excluded from independent biological evaluation.

## Proposed accuracy and runtime targets

[evaluation.json](../../configs/evaluation.json) is the machine-readable proposal. Its concrete thresholds and budgets are for review and later implementation; the current pilots do not establish that the final methods meet them.

- Four independently initialized chains per condition; distinct development, evaluation and reference seeds; the same precision, host/device and tuning budget across methods. Record actual states, RNG algorithm, package lock and thread settings.
- Exact transition checks: maximum row-sum and stationarity error `1e-12`, with only tiny numerical roundoff tolerated in nonnegativity checks. A claim of reversibility additionally requires detailed balance. Check invalid and rejected proposals explicitly.
- Exact-case empirical total variation at most `0.02`. For every single-site and selected-pair category, absolute probability error at most `0.02`. Energy and Hamming mean errors at most `0.05` reference standard deviations.
- Across the four chains, require rank-normalized split R-hat below `1.01` and bulk/tail ESS at least `400` for energy and Hamming. The numerical diagnostic thresholds follow the general recommendations in [Vehtari et al.](https://sites.stat.columbia.edu/gelman/research/published/Vehtari_etal_2020_rhat_ess.pdf); diagnostics alone do not prove convergence. The probability/mean tolerances are project design choices.
- For broad cases, qualify separate long direct reference runs first: R-hat below `1.01`, bulk/tail ESS at least `4,000`, frequency MCSE at most `0.0025`, and mean MCSE at most `0.01` reference SD. Require candidate/reference errors plus their combined uncertainty to fit the configured tolerance. Constant observables need exact agreement; a constant observed trace is not proof that the target is constant.
- Proposed caps: 120 seconds tuning per method/case, 900 seconds measurement per condition, and 3,600 seconds reference generation per case. Record capped failures or inadequate references. Review these budgets using development evidence before the Week 13 freeze.

Measure time from input loading through preparation, encoding/compilation, tuning, initialization/warmup, proposal and correction, transfers, diagnostics and reporting. Synchronize asynchronous computations. At the predeclared observation checkpoints, record the first time all applicable targets hold across four chains; reaching a fixed number of updates is insufficient. Include preparation and tuning in the primary cost without amortizing them away. Keep the chosen plan fixed within each run.

Compare direct sampling, a fixed corrected encoded strategy, correction with planning disabled, and planner-selected execution. If two configurations become identical, record that equivalence instead of implying independent algorithms. No speedup result exists yet.

## Required run and export records

Retain the initial state and every post-attempt state, including rejection repeats, plus original energy, proposal/acceptance events, method/schedule, settings, seeds, timings and diagnostics. Preserve warmup in the raw trace and explicitly exclude it only from analysis. Do not deduplicate, redraw invalid endpoints until valid, or repair encoded words without incorporating the changed proposal law.

Export a separately named FASTA linked to its run and selection procedure. With this fixed gap pattern, an ungapped export has 89 residues representing the aligned domain segment, not a reconstructed full-length 355-residue protein. Retain the 96-symbol aligned sequence and coordinate map for traceability. No exported candidate is represented as experimentally validated.

## Admission and remaining reviews

Technical preparation comprises pinned bytes, complete unique key coverage, finite support, exact coordinate matching, nonempty restrictions, full-model score checks and reference probabilities for the two small cases. Each is reproducible from the inspector and retained evidence.

Both leads still need to record independent real-model score/position checks and approve the restrictions and accuracy proposal. The adviser review must record decisions before the benchmark is described as agreed. The actual encoded proposal and correction derivation are due for the Week 8 specialist review. The [review packet](review_packet.md) lists these decisions without fabricating approvals.
