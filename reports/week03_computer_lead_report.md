Protein Sequence Sampling Project
Week 3 Report - Computer Lead

This week, I built on the software pilots from Week 2 by checking the mathematics and reusable components behind our planned sampling workflow. My scheduled tasks were to verify score calculations and position numbering, inspect Gibbsiq's implementation, test the components we may reuse, and turn the findings into draft benchmark requirements. I completed these preparatory tasks using the small synthetic model, keeping the calculations simple enough to check independently. My contribution focused on identifying what the existing software already supports and what I still need to implement before we can compare sampling methods on a protein model.

I started with the sequence AC from the two-position example. Its two field terms and one interaction term give a score of 3.5, corresponding to an energy of -3.5. I wrote the arithmetic separately from the model's scoring functions so that the comparison could expose a sign or indexing mistake. The following excerpt from my learning exercise records the reference calculation:

```python
hand_terms = (1.0, 0.5, 2.0)
hand_score = 1.0 + 0.5 + 2.0
```

I also checked that model positions 0 and 1 correspond to teaching-reference positions 1 and 2. This helped me distinguish sequence coordinates from alphabet indices: A and C have indices 1 and 2 in the alphabet "-AC", but those numbers describe residue choices, not their locations. Making that distinction explicit gives me a clearer basis for the later protein loader. The exercise verifies the synthetic example; mapping the actual protein's biological positions still requires separate work with the biology lead.

For the reuse exercise, I connected the same model to Gibbsiq's categorical representation and domain-wall encoding, which represents each three-choice position with two binary variables. All nine sequence energies matched across the original, categorical, QUBO, and Ising representations, and every valid encoded sequence decoded correctly. I also checked all 16 possible four-bit combinations: nine were valid and seven were rejected. These results support reusing the representation components. They also helped me understand why matching energies is only one requirement: a correct encoding alone does not show that a sampling procedure will visit sequences with the intended probabilities.

To investigate that second requirement, I constructed a three-state example in which a move and its reverse have different proposal probabilities. The Metropolis-Hastings acceptance calculation must account for that difference as well as the energy change. In my implementation, the final calculation uses logarithms and caps acceptance at one:

```python
log_ratio = -(energy_after - energy_before) / temperature + math.log(q_reverse) - math.log(q_forward)
return math.exp(min(0.0, log_ratio))
```

With the proposal ratio included, the maximum change in a target-state probability after one transition was about 2.78e-17, effectively zero at the calculation's numerical precision. Omitting that ratio increased the error to approximately 0.107. This was the clearest lesson of the week because it showed how a plausible energy-based rule could still produce the wrong distribution. These were exact checks of a small transition matrix, so they establish the behavior of this example rather than convergence or performance on the future protein benchmark.

The rejection checks also extended what I learned about repeated samples in Week 2. I used a prescribed sequence of moves that included both an ordinary rejection and an invalid binary endpoint. The saved trace was 10, 10, 11, 11, 10, retaining the current state after each rejection. Those repetitions record time spent in a state and must remain in the analysis. This exercise led me to specify that the future sampler should record every attempted update and preserve rejection repeats instead of silently dropping or repairing invalid outputs.

Nine learning tests and 69 focused Gibbsiq tests passed. One challenge was that the supplied Gibbsiq description reflected an earlier development stage than the implementation I inspected. I resolved this by recording the actual version, 0.2.0a1, the source commit, and hashes of the audited files. The checks left the upstream code unchanged. The audit also clarified that the reusable core does not provide the complete corrected protein sampler we need. I will still need a protein parameter loader, positional restrictions, a defined proposal procedure, and a correction justified for that procedure's actual update schedule.

The related-work review helped me describe the project's contribution more precisely. Encoding and sampling correction are established techniques; our planned contribution is a reproducible comparison and a measured rule for choosing between methods. I therefore carried forward requirements for identical targets and restrictions, separate development and evaluation cases, and timing that includes preparation and correction costs. I kept the original milestone order while making these checks explicit. Next week, I will assemble the concrete benchmark proposal and updated plan, coordinate the position map and allowed residues with the biology lead, and prepare the package for both leads and the adviser to review.
