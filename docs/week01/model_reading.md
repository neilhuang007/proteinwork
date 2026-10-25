# Week 1: read a tiny adabmDCA model

The goal is to read fields and couplings, encode a sequence, and explain the
sign of its energy. The checked-in fixture is **synthetic and hand-authored**.
It was never fitted to an alignment and makes no biological or fitness claim.
This exercise uses Python's standard library. It does not install, import, or
execute the full adabmDCA package.

## Run it in PowerShell

From the repository root:

```powershell
& .\.venv\Scripts\python.exe .\examples\week01\model_reading.py
& .\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_week01_model_reading.py' -v
```

There is no network access during either command. The exercise reads
[two_site_params.dat](../../examples/week01/fixtures/two_site_params.dat) and
[two_site_metadata.json](../../examples/week01/fixtures/two_site_metadata.json),
then writes [results.json](../../examples/week01/output/results.json).
Use `--output .\somewhere\results.json` to choose a different output file.

The saved JSON contains the input metadata, parameter-file SHA256, every
sequence's numerical encoding, separate terms, score, energy, hand reference,
Boltzmann weight and probability. It is local evidence of the teaching
calculation, not validation against a trained model or an independent human
review.

## Read the records

The current adabmDCA writer emits these whitespace-separated records:

```text
J i j symbol_at_i symbol_at_j value
h i symbol_at_i value
```

Positions start at **0**. Symbols are letters or `-`, rather than numerical
alphabet indices. The canonical writer stores `i < j` coupling records.
These details come from the official
[parameter reader and writer](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/io.py#L91-L219)
and [writer loop](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/io.py#L347-L364).

For example, `J 0 1 A C 2.0` connects `A` at the first position to `C` at the
second. `J 0 1 C A -1.0` is a different entry; swapping the letters while
keeping the positions changes the sequence. A mirrored tensor entry would
instead swap both positions and letters.

The official protein alphabet is `-ACDEFGHIKLMNPQRSTVWY`, with gap `-` at
index 0. RNA and DNA have different token strings. The current CLI default is
automatic alphabet detection; the model API requires an explicit resolved
alphabet. See the official [alphabet definitions](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/alphabet.py#L2-L16),
[CLI alphabet selection](https://spqb.github.io/adabmDCApy/training/#choosing-the-alphabet)
and [model-loading guidance](https://spqb.github.io/adabmDCApy/high_level_api/#load-and-reuse-a-model).

This fixture explicitly chooses the custom alphabet `-AC`, so gap maps to 0,
`A` to 1 and `C` to 2. It has two positions and nine possible sequences. The
metadata records this order; no alphabet is guessed from the file.

## Work out energy and probability

In this exercise, define the unnormalized score and energy as

```text
S(x) = h_0(x_0) + h_1(x_1) + J_01(x_0, x_1)
E(x) = -S(x)
w_T(x) = exp(-E(x)/T)
P_T(x) = w_T(x) / sum_over_all_nine_states(w_T)
```

Each position pair contributes once. The upstream implementation uses a
symmetric full coupling tensor, so its all-pairs sum needs a factor of `1/2`.
With the canonical `i < j` records used here, there is no additional `1/2`.
The sign and counting follow the official
[energy implementation](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/statmech.py#L24-L33).
The official [exact partition calculation](https://github.com/spqb/adabmDCApy/blob/c1e5e2b261444001eb48dd3a0d7c0cded9823f37/adabmDCA/statmech.py#L124-L139)
uses `exp(-E)`, corresponding to `T = 1`. The explicit positive `T` in this
exercise is a dimensionless teaching scale, not a laboratory temperature in
kelvin.

For `AC`, the independent hand calculation is
`S = 1.0 + 0.5 + 2.0 = 3.5` and `E = -3.5`.
For `CA`, it is `S = -1.0 - 0.5 - 1.0 = -2.5` and `E = 2.5`.
At the same positive temperature, lower energy gives greater Boltzmann weight.
The score and weight are not normalized probabilities. A model energy also
does not by itself establish laboratory fitness or safety.

Before running, predict the encodings and energies for `A-`, `AA`, `AC` and
`CA`. Check that omitted couplings contribute zero. Then read
`REFERENCE_ENERGIES` in [model_reading.py](../../examples/week01/model_reading.py):
its nine literal arithmetic expressions do not reuse the parser or the
energy implementation. The run fails if any computed energy disagrees.

## Expected results at T = 1

| Sequence | Encoding | Score | Energy | Probability |
| --- | --- | ---: | ---: | ---: |
| `--` | `[0, 0]` | 0.00 | 0.00 | 0.023288 |
| `-A` | `[0, 1]` | -0.50 | 0.50 | 0.014125 |
| `-C` | `[0, 2]` | 0.50 | -0.50 | 0.038396 |
| `A-` | `[1, 0]` | 1.00 | -1.00 | 0.063304 |
| `AA` | `[1, 1]` | 0.75 | -0.75 | 0.049301 |
| `AC` | `[1, 2]` | 3.50 | -3.50 | 0.771203 |
| `C-` | `[2, 0]` | -1.00 | 1.00 | 0.008567 |
| `CA` | `[2, 1]` | -2.50 | 2.50 | 0.001912 |
| `CC` | `[2, 2]` | 0.25 | -0.25 | 0.029903 |

The terminal ends with `Hand-reference check: 9/9 states passed` and
`Partition sum: 42.939976; probability sum: 1.000000`.
Display values are rounded; the JSON keeps Python float precision. There is
no random seed because this enumerates all states exactly and does no
sampling. Tests check sign, pair counting, alphabet order, reference energies,
probability normalization and teaching-input errors.

## Teaching-only limits

The parser accepts only two positions (0 and 1), an explicit alphabet of
distinct protein symbols, complete `h` records including zeros, and `J 0 1`
records. Omitted couplings are zero. Full-line `#` comments and blank lines
are accepted. It rejects duplicate records, unknown or multi-character
symbols, malformed values, nonfinite numbers, self-couplings, mirrored
couplings, legacy numerical-symbol formats and additional positions.
Sequences must have exactly two valid uppercase symbols; nothing is silently
normalized. Weight calculation rejects nonpositive or nonfinite temperature
and values outside this simple float calculation's range.

These are exercise restrictions, not the complete adabmDCA input contract.
For example, the upstream [import/export documentation](https://spqb.github.io/adabmDCApy/api/adabmDCA.io/#function-load_params)
describes support for one or both coupling triangles. A successful run here
does not validate production models, gauge transformations, alignment
coordinates, training quality or biological interpretation. The planned
validated model loader belongs to Week 5 and sampling to Week 6.

The format, alphabet and energy sources were checked on 2026-10-08 at upstream
commit `c1e5e2b261444001eb48dd3a0d7c0cded9823f37`. This pins the cited source
conventions without claiming that the upstream package was installed or run.
