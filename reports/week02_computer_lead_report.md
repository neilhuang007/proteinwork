Protein Sequence Sampling Project  
Week 2 Report — Computer Lead

This week, I moved from studying model parameters to building small, reproducible examples with the software we plan to use. My scheduled tasks were to run a tiny THRML model, support the adabmDCA sequence-generation pilot, and document the setup and results so that both examples could be repeated. To do this, I created an isolated Python 3.13.5 environment, pinned the package versions, and wrote Python scripts around the official sampling interfaces. I completed these pilot tasks for the Week 2 milestone, giving the project a working starting point for later implementation. My contribution focused on connecting the tools, checking their outputs, and recording the evidence needed to understand each run.

For the systems pilot, I used two interacting binary variables, which gave me four possible states that I could examine individually. I explicitly mapped Boolean values to spins of minus one and plus one and placed the interacting variables in separate update blocks. With those conventions established, I configured 128 warmup sweeps, 4,096 observations, and two sweeps between observations. The following excerpt from my THRML script shows how I applied those settings and waited for the calculation to finish before recording its duration:

```python
schedule = SamplingSchedule(
    n_warmup=configuration["n_warmup"],
    n_samples=configuration["n_samples"],
    steps_per_sample=configuration["steps_per_sample"],
)
started = time.perf_counter()
returned = sample_states(sample_key, program, schedule, initial, [], [Block(nodes)])
jax.block_until_ready(returned)
elapsed = time.perf_counter() - started
```

To assess whether the output was reasonable, I also wrote a small reference calculation that enumerates all four states. After calculating each state's energy, I converted the energies into probabilities using the code below. Subtracting the smallest energy keeps the exponential calculation numerically stable while preserving the normalized probabilities. This gave me a concrete reference against which to compare the observed frequencies.

```python
shift = min(row["energy"] for row in rows)
weights = [math.exp(-(row["energy"] - shift)) for row in rows]
partition = math.fsum(weights)
for row, value in zip(rows, weights, strict=True):
    row["probability"] = value / partition
```

The THRML pilot returned all 4,096 requested observations, and its total variation distance from the exact distribution was approximately 0.00228. This measures the overall difference between the observed and expected state frequencies. Alongside that work, the adabmDCA pilot returned 2,048 sequences from the two-position synthetic model, and all nine sequence energies matched the hand-calculated references. Its total variation distance was approximately 0.0157. These results support the basic setup and numerical conventions, although they do not establish that either method will perform well on a larger protein model. In particular, the smaller THRML error does not make it the better sampler: the two pilots used different models, sample counts, and sampling layouts.

Working through that difference helped me understand why the structure of the saved output matters. adabmDCA returned one final sequence from each of many chains, whereas THRML returned observations along one trajectory, where successive observations can be related. I therefore retained every returned sample, including repeated states, and saved the parameters, random settings, package versions, and warnings alongside the results. The THRML record also includes its initial state and random keys. These records preserve what the APIs returned, while a future production sampler will need a more detailed history of every attempted update. Saving this information now gives me useful examples for designing that later output format.

To check the supporting code and saved results, I used six tests covering the energy convention, temperature scaling, duplicate retention, and consistency between the recorded samples and their summaries. All six passed, and the dependency check confirmed that the 61 installed packages were compatible. However, setup did reveal a practical issue: the package cache and project environment were on different drives, so installation fell back from hardlinks to copying files. I addressed this by making copying explicit in the reproduction instructions. I also recorded two upstream TorchScript deprecation warnings from adabmDCA; execution succeeded, so I treated these as compatibility information to preserve rather than a failed pilot.

These findings led me to keep the dependencies isolated and the examples small while retaining the overall milestone order. I also made the timing records explicit about their limits, since the first THRML call includes compilation and the pilot measurements do not establish a speed advantage. Next week, I will build on this foundation by checking score calculations and position numbering more closely, testing the Gibbsiq components we may reuse, and studying how encoding and sampling correction fit together. As I do so, I will coordinate the input requirements with the biology lead so that the eventual implementation uses the intended amino-acid choices and sequence positions consistently.

The code excerpts come from [thrml_pilot.py](../examples/week02/thrml_pilot.py) and [pilot_support.py](../examples/week02/pilot_support.py). Supporting evidence is recorded in the [THRML run log](../docs/week02/systems_pilot_log.md), [adabmDCA run log](../docs/week02/biology_pilot_log.md), and [Week 2 validation record](../docs/week02/validation.json).
