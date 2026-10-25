# Week 2: official-toolkit pilots

This milestone installs and runs the official adabmDCA and thrml packages on
CPU. It uses the Week 1 synthetic two-position model for sequence generation
and a hand-specified two-spin Ising model for a systems pilot. Neither pilot
is a biological benchmark, a fitted model, or the project's final sampler.

The retained runs use released adabmDCA 1.0.0 and thrml 0.1.4, with Python
3.13.5 on Windows x86-64. All 61 resolved packages and their hashes are pinned
in [the environment lock](../../requirements/pilots-win-py313.txt). The lock
targets this platform; resolve a separate lock before using another platform.
The input requirement file pins the two official toolkits. Neither notebook
extras nor a training loop is used.

## Set up and run in PowerShell

From the repository root, a clean checkout with `uv` available:

```powershell
uv venv .venv-pilots --python 3.13.5
uv pip sync requirements/pilots-win-py313.txt --python .venv-pilots/Scripts/python.exe --require-hashes --link-mode=copy
uv pip check --python .venv-pilots/Scripts/python.exe
& ./.venv-pilots/Scripts/python.exe examples/week02/adabmdca_pilot.py
& ./.venv-pilots/Scripts/python.exe examples/week02/thrml_pilot.py
& ./.venv-pilots/Scripts/python.exe -m unittest discover -s tests -p 'test_week02_*.py' -v
```

`uv` can fetch the exact interpreter if it is not installed. Here it reused an
installed Python 3.13.5. The environment directory is disposable and ignored.
No global package installs, server work, fitted models, or external submissions
were made. Package installation needs network access; both pilots run locally
without fetching model data. Each script accepts `--output <path>`; adabmDCA
also saves a FASTA beside its JSON. A rerun overwrites those output files.

To deliberately refresh dependency resolution, record and review a new lock:

```powershell
uv pip compile requirements/pilots.in --python .venv-pilots/Scripts/python.exe --python-platform windows --generate-hashes --output-file requirements/pilots-win-py313.txt
```

The completed run evidence is [adabmdca.json](../../examples/week02/output/adabmdca.json),
[adabmdca.fasta](../../examples/week02/output/adabmdca.fasta), and
[thrml.json](../../examples/week02/output/thrml.json). Full returned samples,
including repeats, appear in these artifacts, alongside exact probabilities,
seeds, actual configurations, package versions, and installed-source hashes.
The [biology log](biology_pilot_log.md) and [systems log](systems_pilot_log.md)
explain the scientific boundaries. [Validation](validation.json) records the
actual test command and results. Run logs retain terminal output.

Primary references checked on 2026-10-08:

- [adabmDCA high-level API](https://spqb.github.io/adabmDCApy/high_level_api/):
  model loading and the documented `model.sample` sequence-generation call.
- [THRML getting started](https://docs.thrml.ai/en/latest/getting-started.html):
  Ising model, graph-coloured blocks, schedule, initialization, and sampling.
- [THRML spin models](https://docs.thrml.ai/en/latest/examples/02_spin_models/):
  spin-valued Boltzmann energy and its sign convention.

## Reading the comparison

Total variation is `0.5 * sum(abs(empirical_probability - exact_probability))`
over the entire enumerated state space. It is a descriptive finite-run error.
The tests check energy sign, pair counting, spin mapping, exact normalization,
and whether saved counts/FASTA agree with every returned state. They also load
the official adabmDCA model and check its nine energies. They deliberately do
not turn a Monte Carlo error tolerance into a convergence certificate.

Observed pilot total variation: adabmDCA **0.015731** (2,048 chains after 128
sweeps), THRML **0.002280** (4,096 observations from one tiny trajectory).
These sample layouts differ, so their error values and timings do not compare
sampler performance. Weeks 5 and 6 still own the validated production loader
and project samplers.
