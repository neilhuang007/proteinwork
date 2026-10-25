# Your computer lead assignment for protein sampling with Extropic

## What your part of the project is

You are building the connection between a trained protein model and Extropic's sampling tools. The biology lead explains the protein task, reference positions and allowed amino acid choices. You turn those inputs into a working sampler, verify its behavior, and investigate how to execute it efficiently.

A sampler generates candidate sequences according to probabilities. It can return the same sequence repeatedly, and some candidates should occur more often than others. Your job is to preserve those probabilities while changing the representation and execution method.

The main comparison is between a direct sampler, which works with amino acid choices themselves, and an encoded method, which represents those choices using binary variables and uses the Extropic path to propose changes.

Our question is whether the encoded method reaches the same sampling accuracy sooner. A small planner will use measured results to choose which verified method to run for a given case.

## What Extropic is doing here

THRML is Extropic's software library for constructing and sampling probabilistic models. Our repository already contains a small THRML CPU pilot. It is the starting point for learning the interface and checking the numerical behavior.

The protein model has many possible choices at each position. The planned Extropic path represents these choices in a binary model, generates proposed changes, and applies a correction against the original protein model before accepting the next sequence.

You will develop this path in software, examine which operations can run together, and record the requirements for actual device execution. Identify the Extropic hardware, runtime, model representation and access method before attaching hardware measurements to the comparison.

Extropic's [THRML repository](https://github.com/extropic-ai/thrml) describes the library and its backend interface. Use our recorded dependency versions for reproducible implementation and investigate interface changes deliberately.

## Understand the difference between possible sequences and sampled sequences

Consider four positions with choices A or V, then C, then D or G, then E. The four possible sequences are ACDE, ACGE, VCDE and VCGE.

For this teaching example, suppose the model gives them weights of 4, 3, 1 and 2. The total weight is 10. Dividing each weight by 10 gives probabilities of 40%, 30%, 10% and 20%.

The complete list describes what can occur. A sample history describes what the sampler actually generated over time. It might begin ACDE, ACGE, ACDE, VCGE, VCGE.

For observations from the target distribution, 1,000 recorded states have expected counts of 400, 300, 100 and 200. Actual counts vary. Neighboring observations can also be related because each update starts from the previous state.

Your software should reproduce the intended probabilities within the agreed statistical accuracy. Listing every possibility is useful for calculating exact answers on tiny test models. The production task is to explore a much larger sequence space by sampling.

## 1. Make the model and biology inputs precise

Your input loader turns saved model files and biological requests into a consistent representation that every method can use.

The main inputs are:

1. The trained model's numerical parameters and source version.
2. The alphabet, meaning the order of amino acid symbols used in those parameters.
3. The reference sequence and the mapping between biological residue numbers and computer positions.
4. The allowed amino acid letters at each position, including fixed residues and gaps.
5. The model temperature and sampling settings.
6. Starting states, random seeds and output requirements.

A random seed initializes the software's random number generator so a run can be reproduced with the recorded setup. Temperature is a numerical parameter controlling how strongly score differences affect sequence probabilities.

Check the indexing carefully. Biological positions normally start at 1, while computer arrays often start at 0. An alphabet index identifies a residue choice, which is a separate concept from a sequence position.

The current CM model uses 96 alignment positions to represent an 89 residue domain. Review its map with the biology lead. Preserve the model's excluded choices exactly and reject an input whose requested letters conflict with those exclusions.

Keep the biological coordinate reference distinct from the valid sequence used to initialize a run. The current CM reference has two residue choices outside the saved model's support, meaning choices that the model gives zero probability.

**Your result:** A loader that returns checked model parameters, position information and allowed choices, with readable errors for inconsistent input.

## 2. Verify how a complete sequence is scored

The model contains a contribution for each selected residue and a contribution for each selected pair of residues. The complete score combines those contributions according to the model's convention.

The current project's energy convention is recorded in the benchmark contract. Energy here is the numerical model score used to define sampling probabilities. At a fixed positive temperature, lower energy means greater statistical weight.

Before building the sampler, verify the arithmetic on a tiny example. Calculate its score by hand and compare it with your implementation. Then perform the independent score check for the saved protein model with the biology lead.

Fixed positions still contribute to the score. If a fixed residue interacts with a variable residue, changing the variable residue can change that pair contribution. Preserve those interactions when simplifying the calculation.

On tiny models, calculate the score of every possible sequence, convert scores to weights, and normalize the weights to probabilities. That gives an exact target against which you can check a sampler.

**Your result:** A trusted scoring function and small examples with independently checked probabilities.

## 3. Build the direct sampler first

A direct sampler operates on the amino acid choices themselves. It gives you a reference method for the later encoded comparison.

One basic update works like this:

1. Start from a complete allowed sequence.
2. Choose a position that can vary.
3. Hold the other positions fixed.
4. Calculate the model weight for each allowed residue at the selected position.
5. Convert those weights into probabilities and randomly choose the new residue.
6. Record the resulting complete sequence.

This type of update is called Gibbs sampling. The probability of the new residue depends on the current choices at the other positions. Repeating the updates lets the sequence explore the model.

A chain is the history of complete sequences produced by one run. Neighboring states are related because an update changes only part of the current state.

Start with the tiny examples whose target probabilities are known. Then connect the same method to a small protein case. Save the settings, initial state, complete history and model energies.

**Your result:** A working direct sampler whose representation, scores and tiny model behavior have been checked.

## 4. Encode amino acid choices for the Extropic path

Encoding means representing a choice using another set of variables while preserving its meaning.

For a position allowing A or V, one binary variable can use 0 for A and 1 for V. A binary variable has two possible values. A position with more choices needs a larger representation.

The existing project preparation investigates domain wall encoding. In that representation, ordered choices are represented by patterns of binary values. The [implementation roadmap](implementation_roadmap.md) and Week 3 exercises contain the detailed construction.

For every allowed residue, check that its encoded representation decodes back to the same residue. For a full sequence, verify the relationship between the encoded calculation and the original model score.

Also identify invalid binary patterns. A binary representation may contain patterns that do not correspond to an allowed protein choice. Define their treatment as part of the proposed update method.

Connect the encoded model to a small THRML sampling program. Begin with a problem small enough that its states and update probabilities can be inspected directly.

**Your result:** A verified encoding, decoding procedure and small executable THRML proposal program.

## 5. Correct the proposals against the original protein model

Encoding and execution details can affect how often different changes are proposed. The method must account for that behavior to preserve the intended protein distribution.

The original model defines the target distribution: the probabilities we want the final sequence history to follow. The encoded program produces a proposal, meaning a suggested next sequence.

For the planned corrected method, the host program then decides whether to accept the proposal. That decision uses the original model and the mathematical relationship between proposing a move and proposing the reverse move.

The implementation roadmap describes the proposal and correction requirements. Derive the correction for the actual update schedule you will execute and review it with the sampling specialist. An update schedule specifies which variables are updated, in what order, and how many times.

If a proposed encoded state is invalid or a proposed move is rejected, retain the current complete sequence and record it again. For example, a history might be ACDE, ACGE, ACGE. The last repeated state records that the process stayed at ACGE.

Preserve those repeats in the raw history. Changing how invalid proposals are handled or retrying until a valid proposal appears changes the proposal process and must be included in its mathematics.

**Your result:** A fully specified proposal and correction method that passes the tiny model checks and records rejected steps correctly.

## 6. Check whether the sampler works

You need several kinds of evidence.

### Check individual states and calculations

Every accepted state must decode correctly, satisfy the required restrictions and have the correct original model score. Verify these properties independently of the sampling loop.

### Check the update rule mathematically on tiny models

For a tiny model, calculate the probabilities of moving from each state to every possible next state. These probabilities include remaining at the current state.

Check that the intended target probabilities remain unchanged after applying the update. This property is called stationarity. Also check that the update can explore the intended states. The detailed tests and proposed numerical tolerances are recorded in the benchmark contract.

This checks what the algorithm is designed to do. Observed frequencies then check the implemented run.

### Compare observed frequencies with the target

For the teaching example, compare the measured frequencies with 40%, 30%, 10% and 20%. Quantify the uncertainty instead of expecting exact counts from a finite random run.

Independent runs will usually produce different sequence histories. Compare their statistical behavior. Identical sequence lists are not the acceptance criterion for two different samplers.

### Check exploration and agreement between runs

Use four independently initialized chains for each comparison condition. Examine whether they agree on energies, residue frequencies, selected pairs and reference distances.

Repeatedly staying in the same region can give you many observations with little independent information. Effective sample size estimates how much independent information the correlated observations contain. The existing evaluation configuration specifies the proposed diagnostics and thresholds.

For larger protein cases, compare against separately qualified long direct runs and include uncertainty in that reference.

**Your result:** Evidence that the implementation follows the intended target and that the measured runs explore it well enough for the agreed accuracy.

## 7. Work out the parallel execution

The simplest parallel work is running independent chains at the same time. Each chain has its own starting state, random stream and current sequence. They use the same original model and restrictions.

You can also investigate parallel numerical calculations and permitted blocks of variable updates. A block is a group of variables handled together.

The important question is whether the joint update still follows the intended probability rule. Variables connected by interactions generally influence each other's conditional probabilities. Updating them independently using old information requires its own correctness argument or a valid joint update.

Use the model's interaction structure, the actual THRML program and the correction proof to choose the schedule. Verify the resulting update on tiny cases before applying it to the protein model.

Record the cost of moving information between the host and execution device. Keep the software path reproducible and document the separate steps needed for an actual hardware backend.

**Your result:** A tested execution strategy that increases useful work performed together while preserving the sampling method.

## 8. Save the right outputs

For each production run, preserve:

* The model, reference map, restrictions and temperature.
* Method, encoding, update schedule, starting states and random settings.
* The complete sequence after each attempted outer update, including repeated states.
* Original model energies and proposal or acceptance information.
* Sampling diagnostics and the timing breakdown.
* Software versions, execution device and completion reason.

Retain the initialization period in the raw record and identify the portion used for analysis. Save large histories in manageable chunks.

Produce a separate candidate export for the biology lead. It may contain distinct sequences or a selected subset, according to the agreed selection procedure. Link each exported candidate to its source run.

The raw history answers questions about the sampler's probability distribution. The selected library answers questions about which candidates someone wants to examine.

**Your result:** A reproducible sampling record, statistical summary and traceable candidate library.

## 9. Measure whether Extropic's sampling path helps

Compare the direct and corrected encoded methods on the same model, restrictions, temperature and accuracy requirement.

Measure the time needed to reach the agreed accuracy, including preparation, compilation, initialization, proposals, correction, transfers and reporting. A higher number of attempted updates per second is useful only when those updates produce the required sampling information sooner.

Use development cases to choose settings and investigate performance. From those measurements, build a simple planner that selects a verified method and schedule before a run. Reserve separate cases for the final evaluation and keep the selection rule fixed during it.

Compare direct sampling, a fixed encoded method, and planner selected execution. Record the execution device for every result. Interpret the outcome in terms of actual time saved for the biological task.

**Your result:** A reproducible comparison showing which method reaches the agreed sampling quality sooner and under which conditions.

## What to do first

Read the existing THRML CPU pilot and run the small example in its documented environment. Inspect its model, update schedule, returned states and exact probability comparison.

Next, build the production input loader and score checks. Use a tiny protein case to implement and verify the direct sampler. In parallel, review the biological reference and allowed residues with the biology lead.

Then connect the verified encoding to THRML and derive the correction for the actual proposal program. Progress from tiny examples to the same protein case used by the direct method.

Your first complete milestone is a trustworthy direct protein sampling run. The following milestone is a trustworthy Extropic sampling path using the same inputs and target probabilities. The [project plan](project_plan.md) sets out the shared schedule and deliverables.
