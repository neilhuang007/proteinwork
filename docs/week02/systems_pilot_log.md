# Week 2 systems pilot log

On 2026-10-08, official **thrml 0.1.4** ran a CPU two-spin Ising pilot adapted
from [THRML getting started](https://docs.thrml.ai/en/latest/getting-started.html).
The package supplies `IsingEBM`, `IsingSamplingProgram`, `hinton_init`,
`SamplingSchedule`, and `sample_states`. Local code enumerates four reference
states and records the result; it does not supply the sampling kernel. No
training or TSU hardware was used.

Run [thrml_pilot.py](../../examples/week02/thrml_pilot.py) using the shared
[setup instructions](README.md). It runs locally after package installation.

## Model and conventions

The [official spin-model energy](https://docs.thrml.ai/en/latest/examples/02_spin_models/)
is `E(s) = -beta * (b0*s0 + b1*s1 + J01*s0*s1)`. Boolean `False` maps to spin
`-1`; `True` maps to `+1`. One edge `[0,1]` contributes once. The fixed arrays
are `b=[0.25,-0.125]`, `J=[0.5]`, `beta=1`, all `float32`. These values are
exactly representable in binary, so saved actual model arrays agree with the
declared values. The free blocks are `[[0],[1]]`, with no clamped nodes;
adjacent nodes occupy different blocks. Observation order is `[0,1]`.

The reference enumerates the whole state space:

| Boolean state | Spin state | Energy | Exact probability | Count / 4,096 | Empirical probability |
| --- | --- | ---: | ---: | ---: | ---: |
| `00` | `[-1,-1]` | -0.375 | 0.314756 | 1,293 | 0.315674 |
| `01` | `[-1,+1]` | 0.875 | 0.090179 | 364 | 0.088867 |
| `10` | `[+1,-1]` | 0.125 | 0.190909 | 778 | 0.189941 |
| `11` | `[+1,+1]` | -0.625 | 0.404155 | 1,661 | 0.405518 |

The sign and pair count are tested with independent literal energy values.
An additional check at zero biases confirms aligned/opposed probability
ratios and the inverse-temperature scaling.

## Sampling configuration and retained output

Fixed seed 42 uses explicit JAX `threefry2x32`, split once into initialization
and sampling keys. All three key data arrays are saved. `hinton_init` returns
initial free-block states `[[true],[true]]`; this is an initialization heuristic,
not a draw claimed to be at equilibrium. There is one unbatched trajectory.
The schedule uses 128 warmup sweeps, 4,096 observations, and 2 sweeps between
successive observations. The upstream implementation records the post-warmup
state as the first observation and then runs 4,095 intervals; all returned
states are retained. Thus there are 8,318 Gibbs sweeps in the schedule, with
both free blocks updated in each sweep.

[thrml.json](../../examples/week02/output/thrml.json) stores the complete
returned list of blocks, shape `[4096,2]`, boolean dtype, every repeated row,
actual parameters/dtypes, initial states, keys, exact probabilities, counts,
source hashes, full package versions, CPU device, and lock hash. Total
variation is **0.002280**, with maximum probability error **0.001362**.
The [terminal log](thrml_run_output.txt) records the successful run.

Python 3.13.5, JAX/jaxlib 0.11.2, equinox 0.13.8, jaxtyping 0.3.11, and NumPy
2.5.3 ran successfully on `CpuDevice(id=0)`. JAX 64-bit mode is explicitly off.
The sampling timer synchronizes before and after the call; its first call
includes compilation. No warnings or failed calls occurred in this pilot.
The shared installer hardlink/copy warning is retained in the biology log.

The observations can be correlated. Agreement with four exact probabilities
is a modest setup/scientific-convention check, not a convergence proof, ESS
estimate, or timing comparison with the adabmDCA layout. This binary model
contains no protein encoding and does not advance the Week 6 final sampler
implementation. Larger categorical models and planner cost measurements
remain future work.
