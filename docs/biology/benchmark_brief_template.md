# Biology worksheet for the full project

The [current biology guide](../biology_lead_guide.md) contains the complete assignment and can be answered directly. This optional worksheet retains the original project stage structure for recording additional detail. For each decision, explain your recommendation, the alternative you considered, the evidence, what you expect to happen and what would change your mind.

Biology lead: [name]  
Systems lead: [name]  
Last updated: [date]  
Current milestone: [week or deliverable]

## Core assignment From choices to a candidate protein sequence

Choose one CM case and work through the programmer's displayed choices. Your job is to explain the conversion and make a biological assessment of the resulting candidate.

1. What does each choice mean, and which amino acid and position does it select?
2. How do you assemble those choices into one complete aligned sequence and then the gap-free domain sequence?
3. When a proposed change is accepted or rejected, which complete sequence should be saved next?
4. Which positions should stay fixed, which should vary, and what evidence supports that choice?
5. Which residue combinations could depend on one another, and how would you investigate that relationship?
6. What features would make the sequence a plausible CM candidate for our task, and how would you assess them?
7. How would you distinguish an invalid sequence representation from a biologically surprising candidate?
8. What evidence would make you recommend this candidate for further study, change its allowed choices or investigate another candidate?

Use the four-letter ACDE example in the guide to practice, then complete these outputs for the selected case:

- Choice-to-residue table: [table or link]
- Complete aligned sequence and gap-free domain with explanatory notes: [sequence or link]
- One accepted/rejected update example and the resulting recorded states: [worked example]
- One example that fails a requested rule and an explanation: [example]
- Biological assessment, alternative interpretation and next check: [your argument]

Carry those answers into the sections below. Each section develops a decision needed to move from that one example to the complete project.

## 1 Research task and practical value

Begin with the researcher's task so that the later sequence rules and performance goals have a clear purpose.

**Questions to answer:** Which user task should we prioritize, and why is it more useful than another plausible task? Where would faster sampling change what the researcher can do? What observation would show that another part of the workflow matters more?

My recommendation, alternative and reason: [answer]

- Intended user and biological question: [answer]
- What the researcher wants to vary or preserve: [answer]
- How the candidate library fits their workflow: [answer]
- Desired library size and frequency of use: [answer and source]
- Current generation process and waiting time: [answer and source]
- What faster generation would enable: [answer]
- A useful absolute time saving: [answer to agree with the intended user]
- Paper or conversation supporting this description: [link or dated note]

A task description might begin: “A CM researcher wants several collections with different residue choices. I will establish the collection sizes and how often the choices change, then use those needs to assess the value of faster sampling.”

## 2 Model and source records

Use the task to explain why the selected family and model are appropriate. Connect each source to the claim or decision it supports.

**Questions to answer:** Does this model represent the protein region and biological question we care about? How could its training information affect the candidates we obtain? If sources disagree, which evidence should carry more weight and why? What finding would make you request a different benchmark?

My assessment and the assumption I most need to check: [answer]

- Protein family and reason for the choice: [answer]
- Model source, file name, version and citation: [answer]
- Training information, including natural sequences and experimental labels: [answer]
- Relevant earlier studies and what they contribute: [answer]
- Sequence, structure and annotation sources used: [links and versions]
- Data access or reuse terms: [source and details]
- Model-history questions and person responsible for resolving them: [answer]

For example: “This paper supports the CM design motivation, this sequence record supplies the numbering, and this experimental table explains the activity measurement.”

## 3 Reference mapping and independent score check

Confirm the sequence behind each proposed residue choice, then record the shared arithmetic review.

**Questions to answer:** How would you establish that a residue in a paper is the same position we are changing? If numbering or scores disagree, what are two possible explanations? Which specific sequence comparison or arithmetic check would distinguish them?

My explanations, proposed check and conclusion: [answer]

- Reference accession, organism, sequence version and source: [answer]
- Modeled segment and full-sequence boundaries: [answer]
- Numbering used in each biological source: [answer]
- Requested positions and their expected amino-acid letters: [table or link]
- Gap, insertion, deletion or construct differences: [answer]
- Mapping discrepancies and how they were resolved: [answer]
- Purpose of the coordinate reference and the separate valid starting states: [answer]
- Independent score fixture and source entries checked: [answer]
- My independently calculated field total, pair total and final score: [values]
- Review date, reviewer and any discrepancy: [answer]

The current reference is Q5SLA5 residues 2–90. Use the programmer's displayed map to connect source positions to this segment and record the agreed numbering convention.

## 4 Residue rules and comparison cases

Translate the mapped positions into complete sampling tasks. Copy this entry for each case, adding one row or short entry for each controlled position.

**Questions to answer:** Why should this position be fixed or allowed to vary? Which candidates would each choice remove? Could the effect of a residue change depend on surrounding residues? What alternative rule should we compare during development, and what do you predict it will reveal?

My proposed rule, alternative, prediction and supporting argument: [answer]

- Case name and biological purpose: [answer]
- Position and allowed amino-acid letters: [table or list]
- Reason for each rule and its supporting source: [answer]
- Essential requirement or optional choice to explore: [answer]
- Defaults for positions not listed: [agreed behavior]
- Numbering and model-compatibility review: [result supplied by programmer and reviewed together]
- Biological-use case or numerical-check case: [role]
- Remaining question, owner and next action: [answer]

The four-letter example would say: “Keep positions 1, 2 and 4 as A, C and E; allow D or G at position 3.” Add the biological reason and source when completing a CM case.

## 5 Experimental evidence and dataset roles

Once cases are defined, identify the measurements that help interpret their biological relevance. Copy this entry for each dataset.

**Questions to answer:** What does the assay directly measure, and how does that connect to the property we care about? What other biological explanation could produce the same observation? Which existing measurement would help distinguish the explanations? How does the dataset's previous use affect its role in our evaluation?

My interpretation, competing explanation and evidence needed: [answer]

- Paper, data link and version or access date: [answer]
- Sequence and measurement table locations: [answer]
- Exact tested sequences, construct and reference numbering: [answer]
- Measurement, units or normalization and direction of the scale: [answer]
- Relevant organism, host or assay context: [answer]
- Replicates, uncertainty, missing values and classification thresholds: [answer]
- Use in protein-model or predictor development: [answer]
- Sequence overlap and construct-matching review: [result from joint check]
- Agreed role in our report: [background, model development or independent evaluation]
- Remaining questions and next action: [answer]

Inventory outcome: [dataset chosen for evaluation, or search outcome and the evidence the report will use]

A useful note connects the experiment to its interpretation: “This growth-based test assesses CM variants in the reported host context. We will link each measurement to its tested sequence and record its role in model development or evaluation.”

## 6 Biological summaries and success criteria

Use the cases and evidence to explain which summaries answer the research question. Include the reason for each request.

**Questions to answer:** Could two collections look identical in the summaries you selected while containing biologically different combinations? Work through the paired-position example in the guide: how often do A at position 1 and D at position 3 occur together in each collection? Which CM positions deserve a similar joint check, and why?

My worked comparison and recommended measurements: [answer]

- Sites and residue frequencies to show: [positions and reasons]
- Pairs of positions to show together: [pairs and reasons]
- Reference and positions used for sequence-distance calculations: [answer]
- Diversity summaries, including distinct sequences: [answer]
- Required restriction checks: [answer]
- Assay summary, if selected in Section 5: [answer]
- Interpretation of the programmer's model-score and sampling-uncertainty summaries: [notes or questions]
- Agreed model, allowed choices and computational temperature: [joint decision]
- Biological differences that matter and agreed accuracy criteria: [joint decision]
- Practical time-saving goal, linked to Section 1: [answer]

For the four-letter example, request the frequency of D and G at position 3 and confirmation that the other letters remain fixed.

## 7 Candidate library and export requirements

Turn the summary requests into a usable output package. Agree how the exported collection will be selected from the recorded sampling runs.

**Questions to answer:** How does removing duplicates change the interpretation of a collection? What happens to the guide's 75% and 25% example after keeping one copy of each sequence? Would ranking by score help the user's task, and what information would that selection change? Which records should we preserve to understand both the sampling run and the selected library?

My recommended export policy and its consequences: [answer]

- Requested library size: [answer]
- Unique sequences or another duplicate-handling rule: [answer]
- Selection rule and which part of each run it uses: [answer]
- Additional selection or ranking criteria: [answer]
- Aligned sequences, gap-free sequences or both: [answer]
- Description of the modeled domain and reference map: [answer]
- Sequence names and accompanying information: [requested fields]
- Link from each candidate to its case and source run: [agreed format]
- Example of the desired output: [short example or annotated screenshot supplied during review]

For example: “Export the gap-free domain sequences, retain their aligned forms in the accompanying table, and attach a case name, source run and reference distance to each candidate.”

## 8 Output reviews during implementation

As the direct and corrected methods are connected, use this entry to turn observations into actionable feedback.

**Questions to answer:** What are at least two explanations for an unexpected pattern in the output? What different results would they predict? Which check would you ask the programmer to perform first, and why would it be informative?

My explanations, predicted observations and first check: [answer]

- Review date, case and run identifier: [answer]
- Sequence or report item reviewed: [identifier]
- Expected behavior: [answer]
- Observed behavior: [answer]
- Biological reason this matters: [answer]
- Investigation owner and next action: [answer]
- Resolution or follow-up question: [answer]

An example is: “ACTE contains T at position 3, where the case permits D or G. Please trace that sequence and inspect how the allowed choices were applied.”

## 9 Final experiment decisions

Before the main comparison, connect each case to its role and record the agreed procedure.

**Questions to answer:** What comparison could show that our proposed method is useful or lead us to favor direct sampling? Which factor changes between the cases, and what stays constant? Which results would support each outcome? If we can run only one additional comparison, which should we choose and why?

My comparison design, predictions and priority: [answer]

- Development cases and their purposes: [list]
- Reserved evaluation cases and their purposes: [list]
- Useful differences covered by the case set: [answer]
- Final biological summaries, accuracy criteria and impact goal: [agreed choices]
- Final library selection and export rules: [agreed choices]
- Separate development/evaluation role of any assay data: [answer]
- Case and configuration versions recorded by the programmer: [identifiers]
- Adviser or specialist feedback affecting biology: [dated decisions]
- Freeze date and people who reviewed the experiment: [answer]
- Later change, reason and affected reruns: [complete when applicable]

A decision could read: “Use the development cases to choose settings, then evaluate the reserved cases using the same agreed summaries and time measurements for every method.”

## 10 Main results and biological interpretation

Use the programmer's measurement tables to answer the original task rather than interpreting each chart in isolation.

**Questions to answer:** Which method best serves the task when accuracy, full runtime and the candidate summaries are considered together? If 18 of 20 assay-matched candidates meet the function criterion within a collection of 10,000 distinct candidates, what can you conclude about that subset and the rest? What competing explanation deserves investigation before we settle on the interpretation?

My conclusion, supporting evidence and strongest alternative: [answer]

- Methods and cases that reached the agreed accuracy: [results]
- Difficult cases and reasons investigated: [results and follow-up]
- Absolute times, differences and practical relevance: [results linked to Section 1]
- Planner choices compared with fixed methods: [interpretation]
- Residue, pair and diversity patterns: [interpretation]
- Experimental matching counts: [generated, distinct, exact construct matches, eligible measured matches and exclusions]
- Meaning of the matched experimental measurements: [interpretation]
- Further questions and agreed analyses: [answer]

Draft conclusion:

> For [task], [method] reached the agreed accuracy in [time], compared with [time] for [comparison method]. The candidate collection showed [pattern], which means [practical consequence]. The next question is [question arising from the evidence].

## 11 Reproduction and user walkthrough

Give the teammate who repeats the workflow the biological context needed to interpret its files. Record their experience so that the final instructions improve.

**Questions to answer:** Where is a reader most likely to misunderstand how choices become a sequence? If the teammate makes a mistake, what does that reveal about our explanation? How would you revise it and check whether the revision works?

My diagnosis, revision and result of the second attempt: [answer]

- Teammate and reproduction date: [answer]
- Reference, case and expected interpretation provided: [answer]
- Steps or labels that required explanation: [feedback]
- Changes made and whether they resolved the issue: [answer]
- Confirmation that a candidate can be traced to its case and run: [review result]
- Short biological walkthrough for the demonstration: [paragraph]

A walkthrough should connect the steps: “Here is the reference and question. These residue rules express the task. We use these output summaries to interpret the generated collection.”

## 12 Report, presentation and release

Use the completed sections above to assemble the biological part of the final deliverable.

**Questions to answer:** What is the central conclusion, and why does the evidence favor it over the strongest alternative? Which figure best shows that distinction? What question should we investigate next, and how would either outcome change our understanding of the project?

My final argument and recommended next investigation: [answer]

- Motivation and intended-use section: [draft or link]
- Family, model, data and related-work section: [draft or link]
- Reference mapping, restrictions and case rationale: [draft or link]
- Assay methods and evidence description: [draft or link]
- Results, practical value, limitations and next questions: [draft or link]
- Biological figures and captions with positions, units and meaning: [files or links]
- Citations, final data versions and reuse information: [list or link]
- Sequence example and biological walkthrough for the presentation: [material]
- Joint review of the connection between biology and method: [date and actions]
- Final questions and owners: [answer]

## Ongoing decision log

Use the two weekly check-ins to connect progress across these sections. Copy the entry below whenever a scientific or benchmark decision changes the work.

- Date and participants: [answer]
- Question and evidence considered: [answer]
- Alternatives considered and what each predicts: [answer]
- Decision and reason: [answer]
- Evidence that would change the decision: [answer]
- Owner and next milestone: [answer]
- Completion check and reviewer: [answer]
- Sections or experiment versions affected: [answer]

The programmer can use the [technical handoff reference](technical_reference.md) to transfer these agreed answers into structured project files.
