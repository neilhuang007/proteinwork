# Your biology assignment for protein sampling with Extropic

## What we are trying to build

A protein is a chain of amino acids. We write that chain as a sequence of letters, with each letter representing one amino acid. The order matters: the letter at position 3 describes the third amino acid in the chain. A residue is an amino acid within a protein chain.

We want to generate candidate protein sequences using an existing trained model and Extropic's sampling tools. A trained model contains numerical scores learned from protein data. It scores choices at individual positions and combinations of choices at pairs of positions. These scores determine how likely different complete sequences are to be generated.

Sampling means drawing sequences according to those probabilities. Some candidates appear often and others appear rarely. The output is a collection of sampled sequences that we can examine and use to select candidates for further study.

Your job is to give this computation a biological purpose. You choose the protein question, explain which amino acids should be allowed, investigate the relevant biological relationships, and decide how to assess the candidates. You also explain how the computer's choices become a complete sequence.

Neil will implement the sampler, translate the model into a representation supported by Extropic's tools, check its mathematical behavior, and investigate parallel execution. Parallel execution means arranging suitable operations so that they can happen at the same time.

Work through this document in order and write your answers in the same file. Use ordinary language and include links to the papers or data that support each recommendation.

## Start with a small example

Imagine a sequence with four positions:

* Position 1 allows A or V.
* Position 2 allows C.
* Position 3 allows D or G.
* Position 4 allows E.

The letters are amino acid codes: A is alanine, V is valine, C is cysteine, D is aspartic acid, G is glycine and E is glutamic acid.

Positions 2 and 4 are fixed because each has one choice. Positions 1 and 3 can vary. Choosing A, C, G and E gives the complete sequence ACGE.

There are four possible combinations: ACDE, ACGE, VCDE and VCGE. For this teaching example, suppose the model assigns them probabilities of 40%, 30%, 10% and 20%, respectively.

A sampling run might produce ACDE, then ACGE, then ACDE again, followed by VCGE. Over a sufficiently explored run, the frequencies should approach the model's probabilities. Seeing ACDE more often than VCDE is part of the intended behavior.

We can list every combination in this tiny example to calculate the expected probabilities and check the software. For our protein, we use the model to sample from the much larger space of possible sequences.

The biology decisions determine the protein task and which choices we permit. The model supplies the probabilities among those choices. Your assignment is to explain and justify the biological inputs, then interpret the sampled candidates.

## 1. Choose the protein and explain the engineering goal

First, establish what we want the candidates to help us investigate. This gives us a reason for selecting particular positions, residue choices and summaries.

The current project materials use chorismate mutase, which we call CM. Your recent alignment work uses Kunitz, identified as PF00014. Recommend the family we should use and explain how your existing work contributes. A protein family is a group of related proteins whose sequences and biological roles help us understand our chosen example.

Describe the protein's job in words a classmate could understand. Then describe our engineering goal. Protein engineering means choosing changes to preserve or alter particular properties. We might want to explore sequence variation while preserving an activity, or investigate how different choices affect a property we care about.

Explain:

1. Which protein or family you recommend.
2. What it does.
3. Which function or property we want to preserve.
4. Which features we want to vary or investigate.
5. Who would use the candidates and what they would do with them.

Identify an existing trained model for the chosen family. Record its source and explain the data used to build it. An alignment compares related sequences; the trained model provides the numerical terms the sampler uses. Discuss the model with Neil so that your biological task and the computational input fit together.

Also describe a useful candidate collection: approximately how many sequences someone needs, how often they would generate them, and what faster sampling would let them do. Use a published workflow or a discussion with an adviser or intended user to support that description.

**What to write:** The protein choice, engineering goal, model source and practical use of the candidate collection.

## 2. Establish the reference and the allowed amino acids

The reference sequence is the sequence we use to name positions and describe changes. Record its database identifier, organism and amino acid sequence.

Decide whether our output represents a complete protein or a particular domain. A domain is a region with a recognizable structural or functional role. For a domain, identify where it begins and ends in the full protein and explain its relationship to the rest of that protein.

For the current CM proposal, the reference identifier is Q5SLA5 and the modeled region is residues 2 through 90. That region contains 89 amino acids represented across 96 alignment columns, including seven gaps.

An alignment places sequences alongside one another so corresponding positions can be compared. A gap is a placeholder in the alignment. Consequently, an alignment column number can differ from the residue number in the full protein.

Review the existing position map with Neil. It connects the computer's positions to biological residue numbers. For a position discussed in a paper, match the sequence and amino acid letter so we apply the instruction to the correct place.

Now define the allowed amino acids at each position, using the format from the teaching example. One allowed letter fixes a position. Several allowed letters let the sampler explore variation there. State what happens at positions where you request no additional restriction.

For each requested restriction, record the reference position, original residue, exact allowed letters, biological reason and supporting source. If you use a chemical category, list its letters explicitly. For example, hydrophobic describes a tendency to avoid water; an exact list of residues turns that description into a computer input.

Discuss the sequence length and gap handling with Neil. Our current CM model has a fixed alignment layout. Explain how the aligned representation becomes the domain sequence a researcher will read.

The reference serves as our numbering anchor. Neil will also prepare starting sequences that satisfy both the model's permitted choices and the requested restrictions.

**What to write:** The reference and region, the position mapping review, and the allowed residue choices with their biological reasons.

## 3. Investigate the biological engineering questions

The allowed letters are the starting point. We also need to understand which sequence features and residue combinations matter to the engineering goal.

Use the following topics to guide your reading. Connect every finding to a concrete decision or a result you want to examine.

### How the protein performs its function

Identify residues involved in the function we want to preserve. For an enzyme, this can include residues involved in catalysis, which means helping a chemical reaction occur. For a protein that binds another molecule, investigate the residues involved in that contact.

Use the evidence to recommend which positions should stay fixed and which should vary. If an experiment tested a substitution, record the original residue, replacement and measured effect.

### How its structure connects different positions

A protein chain folds into a three dimensional structure. Residues far apart in the written sequence can be close together in that structure. Their effects may therefore depend on one another.

Look for positions that contact each other, help form the protein's interior, occupy its surface or interact with another molecule. Investigate whether a change at one position is associated with a particular choice at another.

At positions you want to vary, examine relevant chemical properties such as charge, size and interaction with water. Explain how these properties inform the allowed letters or which combinations should be examined.

Determine whether the protein relies on disulfide bonds, metal or cofactor binding, or contacts between subunits or domains. A disulfide bond links two cysteine residues. A subunit is one protein chain in an assembly of chains. A cofactor is a molecule or ion that helps a protein perform its function.

Explain which features matter for our protein and how they affect the sampling task.

### What related sequences tell us

Use the alignment and conservation information you have collected. Conservation describes how consistently a position contains the same or similar amino acids across the sequences being compared.

For an important position, explain what the alignment shows and how that connects to function or structure. Record how gaps and repeated or closely related sequences are treated when calculating frequencies. Then explain what the evidence suggests about fixing or varying that position.

### How the candidates would be made and used

Investigate the properties relevant to the intended application. Expression means how much protein cells produce. Solubility describes how well it remains dispersed in solution. Stability concerns retaining the relevant structure or activity under specified conditions. Aggregation describes the tendency to form clumps.

Explain which properties matter to our task and identify the measurements, sequence features or prediction tools that could assess them.

### Turn findings into usable instructions

For each recommendation, state the positions or features involved, the condition or measurement you propose, the biological reason, the supporting source and a worked example.

Label its role clearly:

* A required restriction defines which candidates the sampler may produce.
* A model preference is represented by the trained model's numerical scores and influences how frequently candidates appear.
* An assessment or selection criterion helps us interpret or choose among sampled candidates.

The first implementation supports fixed residues and allowed letters at individual positions. For a proposed rule connecting several positions, describe the biological requirement precisely and work with Neil on how it enters the sampling method. The model already scores pairs of residues, so some relationships may be investigated through its existing interactions and the resulting pair frequencies.

Replace broad instructions such as “keep the protein stable” with the particular feature, measurement or predictor you propose and how you would interpret its result. Explain the source of any numerical cutoff.

**What to write:** The biological findings and resulting instructions, with required restrictions, model relationships and candidate assessment criteria clearly identified.

## 4. Explain how a sample becomes a protein sequence

Your procedure should let another person follow one sampled set of choices all the way to a readable candidate.

Return to the teaching example. If the sampler chooses V at position 1 and G at position 3, the fixed C and E complete the sequence VCGE. The output letters must be placed in their correct positions.

The sampler then proposes changes to the current sequence. Suppose the current candidate is ACDE and an accepted change replaces D with G. The next recorded sequence is ACGE. If a later proposal is rejected, the next recorded sequence remains ACGE. The history is ACDE, ACGE, ACGE.

Those repeated states are part of the sampling record. They help describe how long the process stays in particular states and how frequently sequences appear.

Work through this process for one real case using inputs displayed by Neil. Describe:

1. What each sampled choice represents.
2. Where its amino acid belongs in the complete sequence.
3. How fixed positions are filled.
4. How alignment gaps are handled.
5. How to check the candidate against the required restrictions.
6. What complete sequence is recorded after an accepted or rejected change.
7. Which features you would examine to interpret the candidate biologically.

For a sequence that breaks a required restriction, identify the exact position and rule. For a biologically surprising sequence that satisfies the restrictions, explain the pattern you would investigate and the evidence that would help interpret it.

Neil owns the numerical sampling and correction algorithm. Your contribution is the biological sequence construction and checking procedure, with examples precise enough for him to implement and review.

**What to write:** A numbered procedure, one annotated complete candidate, one example that fails a required restriction, and one accepted or rejected update with its recorded sequence.

## 5. Decide how to assess the sampling results

We need to assess both the sampling behavior and the usefulness of the candidates.

### Understand the sampling checks

In the small example, the expected sequence proportions are 40%, 30%, 10% and 20%. Neil will check whether the sampler approaches those probabilities and quantify the uncertainty.

You should be able to interpret a residue frequency: it is how often a particular amino acid appears at a particular position. A pair frequency describes how often two specified residue choices occur together. Pair frequencies can reveal relationships that individual position counts miss.

Choose positions and pairs whose frequencies matter to the biological task. Also identify useful summaries of differences from the reference and variation within the candidate collection. Explain the question each summary answers.

Neil will run several independently initialized sampling histories and compare their results. You will review whether the displayed sequences, numbering and summaries have the intended biological meaning.

### Assess the biological goal

Identify known variants, structures, experimental datasets and relevant prediction methods. For an experiment, record the exact sequence or construct tested, the measured property, the conditions and the source. An assay is an experimental test, such as an activity or binding measurement.

Explain how the assay relates to our goal and how its data were used previously. For the current CM proposal, some experimental labels contributed to model development. Record which evidence informed the model or our choices and which provides a separate assessment of the resulting candidates.

Describe how a selected candidate could be assessed experimentally: the property to measure, the reference or control to compare against, and the outcome that would count as success. A control provides a comparison that helps interpret the measurement.

### Specify the candidate library

A candidate library is the collection of sequences exported for further use. Specify its desired size, sequence format, treatment of duplicates and any selection or ranking procedure.

The raw sampling record and the selected library have different purposes. If a run contains 40 copies of one sequence and 10 of another, keeping one copy of each changes their relative representation. Neil will therefore save the raw history and create the requested library as a separate, traceable output.

**What to write:** The biological summaries, experimental evidence, assessment criteria and candidate library requirements.

## 6. Bring the work together for Neil

Your completed document should contain:

1. The protein choice, engineering goal and trained model source.
2. The reference sequence, modeled region and checked position mapping.
3. The allowed residues and biological reasoning behind each restriction.
4. The important structural, functional and chemical relationships to investigate.
5. The procedure for turning sampled choices into a complete sequence and checking it.
6. The relevant experiments, requested summaries and candidate selection criteria.
7. Worked examples that you and Neil can review together.

Start with the protein goal and reference. Then define the choices, investigate the biological relationships and walk through one candidate. Use those decisions to specify the summaries and evidence needed to interpret a collection of samples.

As the software becomes available, review small displayed runs with Neil. Check the sequences and restrictions, explain unexpected biological patterns, and refine the questions the final comparison should answer.

Your final contribution will connect the biological question, the sampling inputs, the candidate results and the practical usefulness of the workflow. Use those explanations for the report, figures and presentation.
