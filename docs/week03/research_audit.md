# Week 3: source audit and mathematical preparation

This is preparatory evidence for the Week 4 benchmark contract, not a production
sampler milestone or a record of either lead's approval. Sources were checked on
2026-10-08. Local upstream code will be tested without copying it into this repo.

## Gibbsiq identity and reuse boundary

The local clean checkout `E:/projects/Gibbsiq` has remote
`https://github.com/neilhuang007/Gibbsiq.git` and commit
`06590b82ba509ecf721d8a7685ed27a319d6b0c1`. The first-party
[PyPI project](https://pypi.org/project/gibbsiq/) links that repository, and its
author matches the supplied `C:/Users/neil_/Downloads/Gibbsiq_math_fixed.md`.
That supplied paper's Stage 1 status is historical; the audited checkout is
version `0.2.0a1`. The core's [MIT license](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/LICENSE)
permits reuse with attribution and license retention. The checkout's aggregate
package expression in [pyproject.toml](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/pyproject.toml)
is `MIT AND Apache-2.0` because its optional Z1T adapter
includes Apache-licensed expressions; this audit does not reuse that adapter.

| Pinned component | Evidence and reuse decision |
| --- | --- |
| [categorical.py](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/src/gibbsiq/categorical.py) | Explicit ordered domains, finite complete energy tables, canonical pair orientation. Candidate energy representation; adabmDCA score terms must be negated to obtain this positive-sum energy convention. No protein parameter loader or biological position map. |
| [domain_wall.py](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/src/gibbsiq/domain_wall.py) | Prefix-one codec, reverse-wall rejection and exact valid-state QUBO/Ising energies. Candidate compiler; caller-supplied finite penalty is explicitly not proven sufficient. Category order changes adjacency and mixing. |
| [reference_sampler.py](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/src/gibbsiq/reference_sampler.py) | Standard-library Ising Gibbs reference with schedules, seeds and transition events. Useful oracle for small binary exercises; does not directly sample categorical amino acids. |
| [verification.py](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/src/gibbsiq/verification.py) | Bounded exact Ising transition enumeration, stationarity and ergodicity checks. Reuse as a comparison check, not an oracle for unimplemented protein correction. |
| [thrml_runtime.py](https://github.com/neilhuang007/Gibbsiq/blob/06590b82ba509ecf721d8a7685ed27a319d6b0c1/src/gibbsiq/thrml_runtime.py) | Optional Ising runtime and replica exchange. Its swap Metropolis rule is not a categorical encoded-proposal correction. Backend integration is outside this core audit. |

A search of `src/gibbsiq` for `Hastings`, `corrected`, proposal ratios and reverse
probabilities found no general corrected encoded protein sampler. Therefore
reuse needs a protein loader, restricted categorical direct baseline, defined
encoded proposal law, original-target correction, and project trace/export
adapters. No upstream code was vendored, installed or edited. The
[reuse check](../../examples/week03/upstream_reuse_check.py) verifies the remote,
commit, clean tree and MIT-license hash before importing; its
[provenance JSON](../../examples/week03/output/upstream_reuse.json) records hashes
for the audited source files. Core reuse is feasible, subject to these adapters.

## Benchmark requirements to carry forward

- Pin original trained parameters, alphabet, pair-count convention, temperature,
  reference sequence and every model-to-biological coordinate. Keep synthetic
  score arithmetic distinct from the two leads' pending real-model checks.
- Define the target on exactly the restricted sequence set:
  `pi(x) proportional to exp(-E_original(x)/T)`. Resolve gap handling, singleton
  sites and empty allowed sets before compiling; use identical restrictions in
  every method.
- Log domain/category ordering, penalty, encoded graph size and density, update
  schedule, forward/reverse proposal semantics, invalid-endpoint handling and
  numerical precision. Energy equality on valid states does not establish
  distributional correctness or useful mixing.
- Prove correction for the actual kernel. General Metropolis-Hastings uses
  `alpha(x,y)=min(1, exp(-(E(y)-E(x))/T) q(x|y)/q(y|x))`.
  This is the general proposal-ratio rule of [Hastings (1970), section 2.2](https://hedibert.org/wp-content/uploads/2013/12/1970Hastings.pdf).
  A finite encoded Gibbs run's endpoint proposal need not have a tractable
  ratio; symmetric schedules or pathwise correction need explicit derivations.
  Reject invalid endpoints and retain the current state. Do not silently repair
  or discard words, or repeatedly redraw until valid without accounting for
  the changed proposal normalization.
- Compare direct, fixed corrected encoded and planner-selected methods on the
  same target. Freeze choices after development; include tuning and correction
  in total time to the agreed accuracy. Preserve every rejected-step repeat;
  filtered FASTA libraries are separate outputs.

In particular, a systematic sweep can preserve a surrogate law without being
reversible; its reverse is the reversed sweep order. An acceptance formula that
replaces `q(x|y)/q(y|x)` with a surrogate energy difference needs a proven
reversible proposal kernel and the actual backend energy/temperature. Uniform
single-bit flips below supply a deliberately simpler, known symmetric kernel.
They do not establish a correction rule for the upstream THRML sweep.

The current Gibbsiq categorical/compiler interfaces reject nonfinite energies.
The [Week 1 CM structural audit](../week01/cm_parameter_audit.json) found
`-inf` fields. Their support semantics must be verified before any adapter
excludes them; replacing them with arbitrary finite penalties changes the model.

## Related-work comparison

| Primary source | What it establishes | Project relationship |
| --- | --- | --- |
| [adabmDCA statmech.py, pinned commit](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/statmech.py#L24-L33) | Protein-model energy is negative fields plus couplings; the symmetric tensor uses half the all-pairs sum. | Preserve the trained target and its conventions; the teaching adapter negates canonical once-counted score terms. |
| [Chancellor (2019), sections II–III, arXiv v2](https://arxiv.org/html/1903.05068v2) | A q-state variable uses q−1 spins with one domain wall between fixed boundaries; pair functions can be encoded by quadratic interactions. | Reuse the encoding idea. Its optimization/embedding results do not prove improved protein sampling accuracy or time. |
| [Jelinčič and Walker (2026), section 4.4 and Supplementary Note 1, arXiv v1](https://arxiv.org/html/2606.17327v1) | Domain-wall finite-difference construction and ordering-dependent local moves in a codon-optimization workload. | The codon choices preserve amino-acid identity, whereas this project's candidate amino-acid sequences change. Its sparse codon graph/color schedule is not a general dense protein coupling graph. |
| [Thermalizing Stochastic Programs, sections III and IV.4, arXiv:2608.01615v1](https://arxiv.org/html/2608.01615v1) | Variationally compiles DFG/PSC conditional kernels using torx, thermalizers and thrml; its meta-EBM example compiles Gibbs updates for non-native energies. Finite compilation produces residual stationary error, measured on an enumerable 12-spin example. | Closely related kernel/schedule viewpoint. Our intended method instead corrects proposals against one fixed original protein target and measures when selecting an encoded method helps. No factor training or generic compiler is proposed. |

Our planned contribution is a reproducible, modest method-selection rule for one
restricted protein benchmark, evaluated against tuned direct sampling and fixed
corrected encoding with equal development budgets and full costs. The encoding,
Metropolis-Hastings rule, diagnostics and compilation abstraction are prior
work. Any performance benefit remains a hypothesis. Thermalizing Stochastic
Programs explicitly distinguishes simulator results, hardware projections and
physical measurements; this project will report local software evidence.

## Hand score, numbering and bounded learning checks

For the **synthetic** Week 1 fixture, the literal check is
`S(AC) = h_0(A)+h_1(C)+J_01(A,C) = 1.0+0.5+2.0 = 3.5`, so `E=-3.5`.
Model positions `0,1` correspond to teaching reference positions `1,2` in `AC`.
Alphabet `-AC` encodes those residues as `1,2`; alphabet indices are not sequence
coordinates. The literal arithmetic and position table in
[learning_checks.py](../../examples/week03/learning_checks.py) are checked by
the independently implemented Week 1 parser/score methods. This is not a map to
natural protein residues, a real-model loader validation, or approval by both
leads.

For three categories, the valid two-bit words are `00`, `10`, `11`; fixed outer
bits are conceptually `1` and `0`. `01` contains a reverse wall and is rejected.
The exercise enumerates an asymmetric three-state proposal, showing why its
proposal ratio matters: corrected stationarity error is `2.78e-17`, while
omitting that ratio gives `0.107`. The symmetric bit-flip kernel redirects
invalid endpoint mass to the diagonal and has zero stationarity error in this
example. These are finite matrix checks, not sampled convergence claims.
The prescribed trace `10,10,11,11,10` retains both an MH rejection and an invalid
endpoint rejection. There is no randomness or seed because its moves and
uniforms are supplied explicitly.

The source adapter checks all nine original/categorical/QUBO/Ising energies and
round trips, then exhaustively verifies all 16 words: nine valid and seven
undecodable. It uses the original fixture's independent literal energy table;
the adabmDCA package itself is not installed by this check.

## Reproduction and results

Using `.venv/Scripts/python.exe` (Python 3.13.5), from this repository:

```powershell
& .\.venv\Scripts\python.exe .\examples\week03\learning_checks.py
& .\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_week03_*.py' -v
& .\.venv\Scripts\python.exe .\examples\week03\upstream_reuse_check.py --source E:/projects/Gibbsiq
```

The learning suite passes **9 tests**, zero skips. Results are saved in
[learning_checks.json](../../examples/week03/output/learning_checks.json), with
[test log](test_run.log), [exercise log](exercise_run.log) and
[source-adapter log](upstream_reuse_run.log).

From the verified Gibbsiq checkout, the focused upstream command was:

```powershell
$env:PYTHONDONTWRITEBYTECODE = '1'
& E:/projects/proteinwork/.venv/Scripts/python.exe -m unittest test_suite.tests.test_categorical_model test_suite.tests.test_domain_wall_encoding test_suite.tests.test_reference_sampler test_suite.tests.test_exact_distribution test_suite.tests.test_statistical_verifier -v
```

All **69 upstream tests pass**, zero skips, in
[upstream_test_run.log](upstream_test_run.log). The tree remains clean. This is a
focused core audit, not a full upstream release test, THRML adapter test,
hardware validation or evidence that a corrected protein sampler exists.
Every command exited 0; [validation.json](validation.json) records the local
checks and remaining human/model decisions.
