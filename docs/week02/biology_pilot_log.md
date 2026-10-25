# Week 2 biology pilot log

On 2026-10-08, the official released **adabmDCA 1.0.0** package loaded the Week
1 hand-authored `-AC` model and ran the documented [sequence-generation
API](https://spqb.github.io/adabmDCApy/high_level_api/#generate-sequences).
This is a synthetic API pilot. The unresolved biological benchmark choice
does not get replaced by a claim that this two-position fixture is biological
data. No alignment was trained or generated as a laboratory candidate library.

Run [adabmdca_pilot.py](../../examples/week02/adabmdca_pilot.py) using the
[setup instructions](README.md). The official calls are `load_model`,
`model.compute_energies`, and `model.sample`; the local code only supplies the
fixture, records samples, and computes an enumeration-based reference. Its
probabilities reuse Week 1's independent arithmetic reference; it does not
implement a replacement sampler.

## Actual configuration and result

- Explicit custom alphabet `-AC`; zero-based positions; length 2; `float32` CPU.
- Metropolis sampler; inverse temperature `beta=1`; fixed seed 42.
- 2,048 chains initialized by the upstream API, 128 sweeps, one terminal
  sequence per chain; Torch CPU thread count 1.
- All nine official energies exactly matched the literal hand references;
  for example `AC=-3.5` and `CA=2.5`. Each pair contributes once.
- Returned all 2,048 requested sequences, with all nine states represented;
  every repeat is preserved in both JSON and FASTA.
- State `AC`: count 1,598, empirical probability 0.780273, exact 0.771203.
- Total variation 0.015731; maximum absolute probability error 0.009070.

The [full JSON](../../examples/week02/output/adabmdca.json) records fixture and
lock SHA256, all installed versions, actual model metadata, energies, counts,
every returned sequence, installed API source hashes, and warnings. The
[terminal log](adabmdca_run_output.txt) records the successful run. The saved
single-call time includes potential compilation and excludes process/import
startup; it is not evidence of a speedup.

## Setup, warnings, and limits

Python 3.13.5 needed no fallback interpreter. `uv 0.9.29` created the isolated
environment and installed 61 packages from the hashed Windows lock. Relevant
versions: Torch 2.14.1, NumPy 2.5.3, Biopython 1.88, pandas 3.0.6. The
[installation output](install_output.txt) preserves an `uv` warning: cache and
environment reside on different drives, so hardlinks fell back to copying.
That warning concerns installation, not numerical results. The recommended
setup command uses `--link-mode=copy` explicitly.

Two `FutureWarning` messages arose from upstream use of `torch.jit.script`,
which this Torch version deprecates. Both messages are retained verbatim in
the JSON. TorchScript execution succeeded; there was no failed compilation
or unavailable compiler warning. The initial [smoke log](import_smoke_output.txt)
also retains a deprecation warning from a short Gibbs test. There were no
installation errors or failed API calls in the observed setup.

This pilot has no reference alignment, automatic mixing-time estimation,
autocorrelation analysis, ESS claim, or biological interpretation. Nine-state
frequency agreement cannot establish behavior on a large protein family. The
Week 4 family agreement and later loader/sampler milestones remain distinct.
