# Protein candidate research workflow

The project now includes a working computational workflow for chorismate mutase candidate libraries: validated model loading, exact residue restrictions, optimized direct sampling, corrected encoded sampling through THRML, diagnostics, constrained candidate selection, FASTA export and an assay data round trip. Physical Extropic device access is not configured. No new biological measurements or improved protein function are claimed.

Read the [research workflow and engineering decisions](docs/research_workflow.md) for the complete input and output contract, applied rules, optimization choices, limits and commands. The [encoded method](docs/encoded_method.md) states the sampler correction and its tests. The [Extropic review](docs/research/extropic_engineering_2026_10_23.md) explains the hardware boundary and reusable libraries.

The research environment extends the earlier pilot environment:

```powershell
uv pip sync requirements/research-win-py313.txt --python .venv-pilots/Scripts/python.exe --require-hashes --link-mode=copy
uv pip install --python .venv-pilots/Scripts/python.exe -e . --no-deps --link-mode=copy
& ./.venv-pilots/Scripts/python.exe -m proteinwork validate --request configs/cm_library_thrml.json
& ./.venv-pilots/Scripts/python.exe -m proteinwork run --request configs/cm_library_thrml.json --output output/my_cm_library
& ./.venv-pilots/Scripts/python.exe -m proteinwork verify output/my_cm_library
```

The archive must already be cached at its pinned hash. Use the existing inspector download command below if needed. The full library request explores the broad CM model, then runs actual THRML batch optimization against a greedy incumbent. `configs/cm_library.json` selects with SciPy MILP. The exact direct and THRML requests validate both sequence sampling paths against known probabilities. The library requests explicitly allow exploratory exports because the broad target has no independently qualified reference yet.

Every run produces raw traces, sampling diagnostics, applied rules, a timing profile, candidate CSV, aligned and ungapped FASTA, an assay worksheet and a hashed manifest. Existing output directories are preserved. Use `python -m proteinwork --help` for validation, execution, verification, comparison, frozen planning and measurement ingestion.

## Historical preparation through Week 4

The material below records the original preparation and milestone plan. Statements about future production implementation describe that earlier stage. The current implementation and its evidence are documented above; historical biological review and assay limitations remain in effect.

The current project uses an existing trained protein model and Extropic's THRML tools to generate candidate sequences. We compare a direct sampler with a corrected encoded method, measuring the total time needed to reach the same sampling accuracy.

Start with the [shared project plan](docs/project_plan.md), then use the [biology lead guide](docs/biology_lead_guide.md) and [computer lead guide](docs/computer_lead_guide.md) for the individual assignments. The sections below record the existing technical preparation and runnable exercises.

Reproducible technical preparation through **Week 4** of the revised 18-week plan. The project compares direct and corrected encoded protein sampling, with a small planner choosing methods by measured cost while preserving a fixed target distribution.

Weeks 1–3 technical deliverables are complete. Week 4 has a concrete model, coordinate map, restriction cases, accuracy proposal and updated plan. **Both-lead agreement and adviser review remain pending.** The production loader, project samplers and planner remain later milestones.

## Reports and deliverables

The [weekly reports](reports/README.md) are concise, meet the assignment's 300-word minimum, and are ready to paste into Canvas. They describe planned milestones without inventing calendar weeks, actual hours, feedback or approvals. No Canvas submission has been made.

| Week | Report | Technical evidence |
| --- | --- | --- |
| 1 | [Report](reports/week01_report.txt) | [Scope](docs/week01/scope_draft.md), [three-family comparison](docs/week01/candidate_families.md), [model-reading exercise](docs/week01/model_reading.md) |
| 2 | [Report](reports/week02_report.txt) | [Setup and runs](docs/week02/README.md), [biology pilot](docs/week02/biology_pilot_log.md), [systems pilot](docs/week02/systems_pilot_log.md) |
| 3 | [Report](reports/week03_report.txt) | [Gibbsiq audit and related work](docs/week03/research_audit.md), [learning results](examples/week03/output/learning_checks.json), [upstream reuse check](examples/week03/output/upstream_reuse.json) |
| 4 | [Report](reports/week04_report.txt) | [Earlier project plan](docs/archive/project_plan_sampling_2026_10_21.md), [benchmark contract](docs/week04/benchmark_contract.md), [admission evidence](docs/week04/benchmark_evidence.md), [review packet](docs/week04/review_packet.md) |
| 5 | [Report](reports/week05_report.md) | [Shared biology and computing exercise](docs/week05/README.md), [residue analysis and candidate annotations](examples/week05/output/summary.json), [fresh direct sampling run](output/week05/cm_direct/manifest.json) |

## Setup and run

The existing `.venv-pilots` environment has Python 3.13.5 and official adabmDCA 1.0.0 / thrml 0.1.4. On a fresh Windows x86-64 checkout, install the exact hashed lock with `uv`:

```powershell
uv venv .venv-pilots --python 3.13.5
uv pip sync requirements/pilots-win-py313.txt --python .venv-pilots/Scripts/python.exe --require-hashes --link-mode=copy
uv pip check --python .venv-pilots/Scripts/python.exe
```

From this directory, run the milestones:

```powershell
& ./.venv-pilots/Scripts/python.exe examples/week01/model_reading.py
& ./.venv-pilots/Scripts/python.exe examples/week02/adabmdca_pilot.py
& ./.venv-pilots/Scripts/python.exe examples/week02/thrml_pilot.py
& ./.venv-pilots/Scripts/python.exe examples/week03/learning_checks.py
& ./.venv-pilots/Scripts/python.exe examples/week03/upstream_reuse_check.py --source E:/projects/Gibbsiq
& ./.venv-pilots/Scripts/python.exe scripts/inspect_benchmark.py --download
& ./.venv-pilots/Scripts/python.exe -m unittest discover -s tests -v
```

The reuse exercise requires the trusted Gibbsiq checkout at the audited commit recorded in its script and the audit note. It verifies the remote, commit, clean tree and license before importing; it leaves that project unchanged. The main project tests do not require this external checkout.

Package installation and the first model download need network access. The pinned CM ZIP is cached under ignored `data/benchmark/raw/`; subsequent inspection is offline. Saved small artifacts allow the project tests to run without downloading the archive. Run scripts overwrite their own result files; timings and observation timestamps can differ. The setup lock targets Windows/Python 3.13; resolve and record a separate lock for another platform.

Weeks 1 and 3 learning exercises, and the Week 4 inspector, use only the standard library and can also run with the original `.venv/Scripts/python.exe`. The full test suite needs `.venv-pilots` because it checks the installed adabmDCA API.

## Verified results

- Week 1: all nine exact-state arithmetic checks passed.
- Week 2: official adabmDCA returned 2,048 sequences and matched nine reference energies; thrml returned all 4,096 requested observations. Their descriptive total variation errors are 0.015731 and 0.002280, respectively. Different sample layouts make these unsuitable as a performance comparison.
- Week 3: 69 focused upstream Gibbsiq tests passed. All nine synthetic energies matched across original, categorical, QUBO and Ising representations. Small transition examples verify the required proposal ratio and rejection retention.
- Week 4: complete model key coverage, exact forbidden-state support, biological coordinates, full-model score arithmetic and four nonempty restriction cases are verified.

[validation_summary.json](docs/validation_summary.json) records the integrated project tests, artifact checks and report word counts. Per-week logs and source fingerprints remain alongside their artifacts. Upstream TorchScript deprecation warnings are retained in the Week 2 records.

## Proposed benchmark and next work

[benchmark.json](configs/benchmark.json) pins the feedback-informed CM checkpoint. The [96-column position map](data/benchmark/position_map.csv) links an unmodified anchor to Q5SLA5 residues 2–90. That anchor is outside the model's support; a separately labelled computational witness supplies valid restricted states. The original 362 negative-infinity fields are preserved as zero-probability choices.

[restriction_cases.json](data/benchmark/restriction_cases.json) separates development and evaluation inputs. [evaluation.json](configs/evaluation.json) proposes accuracy targets, four-chain comparisons, equal budgets and end-to-end timing. These are specified future experiments, not achieved sampling or speedup results. No independent biological assay holdout is claimed.

Complete the [human review record](docs/week04/review_packet.md), then Week 5 implements the support-aware production loader, restriction checks and JAX exercises. The original plan's remaining milestones and contingency rules are retained in the updated project plan.
