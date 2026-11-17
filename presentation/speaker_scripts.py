"""Read-aloud scripts for the fourteen-slide protein sampling presentation.

Slide numbers follow presentation/build.py. Paragraphs are spoken prose, in order.
"""

SCRIPTS = {
    1: {
        'speaker': 'Neil Huang',
        'paragraphs': [
            "Hello, everyone. I'm Neil Huang, and this is Daniel Sun. Today we're going to talk about our QUBO machine for protein sampling. We started with a practical question: can we explore the sequence choices in a protein model by translating them into binary choices?",
            "We'll show you what distribution we're trying to sample, how the translation and correction work, and what our software experiments actually measured. I'll take us through the model and method; then Daniel will walk through the implementation, results, and next experiments.",
        ],
    },
    2: {
        'speaker': 'Neil Huang',
        'paragraphs': [
            "Why sample protein sequences at all? Imagine a landscape containing many possible sequences. Finding its lowest point gives us one answer, but the surrounding valleys can hold many other plausible sequences. Sampling lets us explore that wider landscape: which residues appear often, which choices occur together, and which alternative candidates the model favors. Frequency estimates show us which choices recur across the model's plausible sequences. That is why we look beyond a single answer.",
            "Here the comparison begins with one fixed protein model. The direct route changes residue choices in the original sequence space. The encoded route represents those same choices as binary variables, proposes moves there, and corrects them against the original model. Both routes are meant to arrive at the same sequence probabilities, so they can be judged on their ability to explore those probabilities and on their end-to-end cost.",
            "Our question is whether the binary route can produce reliable samples sooner. That order matters: first check accuracy, then compare elapsed time. I'll start with the target, then the binary method, and Daniel will bring in the measured evidence.",
        ],
    },
    3: {
        'speaker': 'Neil Huang',
        'paragraphs': [
            "The target for our comparison is a fixed model of the chorismate mutase protein family. In this model there are ninety-six aligned positions, and each position can take one of its supported residues. Alignment means matching corresponding sites across sequences; some aligned sites can also represent gaps. Those choices define the space we're sampling.",
            "The five letters in this picture, A, C, G, T, A, are just a small schematic example. They are amino-acid single-letter codes: alanine, cysteine, glycine, threonine, and alanine. These letters simply illustrate the choices and interactions in a sequence.",
            "Why draw that connection? The curved link joins cysteine and threonine to illustrate one pairwise interaction. Protein positions do not necessarily act independently. A residue that fits well at one position may fit differently depending on the residue at another. The model therefore records preferences for individual choices as well as interactions between pairs. Our source model also specifies which choices are allowed; sampling stays within that support rather than changing the model after the fact.",
        ],
    },
    4: {
        'speaker': 'Neil Huang',
        'paragraphs': [
            "Now we need to turn those preferences into a sampling target. A Potts model gives each complete supported sequence an energy. Its individual terms favor some residues at particular positions; its pair terms capture how choices work together. We keep this map fixed for both methods.",
            "The equation on screen says that a sequence's probability is proportional to the exponential of minus beta times its energy. In ordinary language, lower-energy sequences are more likely, while beta sets how strongly the probabilities favor them. The symbol x means the entire sequence, and beta is the inverse temperature. Higher-energy sequences still appear, at model-determined frequencies.",
            "So our job is to draw from the distribution, not merely find its minimum. The direct sampler works with categorical residue choices. The alternative sampler will use binary variables, but we judge its protein sequences against exactly this original distribution. That gives us a clear standard for accuracy: the two routes can take different paths, yet they should represent the same supported sequence probabilities. Now we encode that choice.",
        ],
    },
    5: {
        'speaker': 'Neil Huang',
        'paragraphs': [
            "At a position that permits alanine, cysteine, and glycine, give each allowed residue its own on-or-off bit. Choosing cysteine becomes zero, one, zero. Exactly one bit should be on. Only supported residues get bits, keeping these choices tied to the original model. A group with no bits on, or with two on, does not describe a valid residue choice.",
            "We call the resulting binary energy a QUBO: quadratic unconstrained binary optimization. The name says that the variables are binary and that energy terms involve individual bits or pairs of bits. To encourage one choice per position, we add a finite energy penalty for the invalid patterns. For every valid encoded sequence, the binary energy equals the original protein-model energy exactly.",
            "The penalty is finite, so a short binary trajectory can visit or finish in an invalid pattern. We can also rewrite each zero-or-one bit as a minus-one-or-plus-one spin. That is the Ising form: the same encoded energy expressed with two-state spins and interactions between them. It gives our sampler a way to propose new binary states.",
        ],
    },
    6: {
        'speaker': 'Neil Huang',
        'paragraphs': [
            "Here's the sampling loop. Begin with a valid protein sequence and encode it. Run a fixed, short Ising trajectory: it updates binary spins and proposes an endpoint. For each attempt we flip a fair coin to choose either forward or reverse update order, and use that direction for the prescribed number of sweeps.",
            "Pairing those directions lets us account for forward and backward proposals. Once an endpoint represents a valid protein sequence, we use Metropolis–Hastings, or M H, to correct the move. It compares the starting and proposed states under both landscapes, then accepts the proposal or keeps the old sequence. The recorded chain follows the original protein target.",
            "If the binary endpoint has no residue at a site, or more than one, we reject it once and record the current sequence again. We do not keep drawing binary endpoints until one happens to be valid: that would quietly change the way proposals are generated. Likewise, a rejected valid endpoint still contributes a repeated current sequence to the trace. Those repeats belong in the data. An invalid endpoint is one complete attempt, even though the protein sequence stays put.",
            "Each short walk supplies a proposal for correction. Exact encoding connects the energies, the paired update directions make the proposal accountable, and the final correction keeps the sequence distribution we started with. With that method in place, the next question is how efficiently it actually moves through different sequences.",
        ],
    },
    7: {
        'speaker': 'Neil Huang',
        'paragraphs': [
            "A correct sampling rule and a useful sampling speed are two different things. Imagine walking through a building with a valid map but repeatedly arriving at the same room. Our trace sketch makes that point with a few repeated sequences.",
            "There is a particularly clear case: with the exact encoded energy and matching temperatures, the correction accepts valid binary endpoints. Even then, the endpoint may encode the same protein sequence as before, or an attempt may finish in an invalid pattern and leave us where we started. Counting acceptances alone would miss those limits. Repeated observations also make frequency estimates less precise, so uncertainty tells us where better movement is needed. We therefore look at how many genuinely different parts of the target distribution are explored and how uncertain our estimates remain.",
            "Daniel will now take over, show what we implemented in software, and explain what the measured comparisons and diagnostics tell us. Daniel, over to you.",
        ],
    },
    8: {
        'speaker': 'Daniel Sun',
        'paragraphs': [
            "Thanks, Neil. I'm Daniel Sun. Now let's look at what ran. We implemented a direct route and an encoded route in software. The encoded route can run with NumPy, JAX, or the official THRML software package. In THRML we check the order of individual spin updates, so the executed schedule follows our intended proposal.",
            "We check the method first on small systems where exact sequence probabilities can be calculated and compared. The implementation keeps the supported residue choices intact, and records the seed, starting states, attempted moves, repeated states, timings, and sampling diagnostics. Those records support checks of correctness and exploration.",
            "Our results come from software execution; physical hardware execution and interaction-graph mapping are future work. With the software routes in hand, we can ask the concrete question we promised: when both routes meet the same accuracy checkpoint, which one gets there faster?",
        ],
    },
    9: {
        'speaker': 'Daniel Sun',
        'paragraphs': [
            "We begin with a small, four-state target because we can check both routes against known probabilities. Each method was run three times; each run used four chains and reached the first declared qualification checkpoint at five hundred and twelve observations per chain. All six runs passed that checkpoint. It gives the two methods a common accuracy standard set before their times are compared.",
            "The chart shows elapsed time to that qualified result. For the direct JAX sampler, the median across three runs was seventeen point eight eight five seconds. For the encoded THRML software route, the median was twenty-six point seven five four seconds. In this small case, the encoded route took about one and a half times as long, so the direct sampler won the measured comparison.",
            "This qualified small-case result sets up our next question. The larger protein target has many more interacting residue choices, so a broader comparison needs a trustworthy reference distribution. Before asking whether any route is faster there, we have to know whether our reference estimates are precise enough.",
        ],
    },
    10: {
        'speaker': 'Daniel Sun',
        'paragraphs': [
            "We therefore ran a larger direct-sampling reference on the broad target. Four chains produced a total of one hundred thirty-one thousand and seventy-two observations, with thirty-two thousand seven hundred and sixty-eight observations in each chain. The complete command took about two hundred and two seconds. The longer software run improved several diagnostics.",
            "One diagnostic is effective sample size. Because successive samples can resemble one another, a long trace may carry less independent information than its raw length suggests. For energy, it reached five thousand eight hundred and thirty-three, above the required minimum of four thousand. The reference chain diagnostics for energy and distance also passed.",
            "Qualification, however, requires the complete set of declared checks. Mean uncertainty for both energy and distance still exceeded its limits. Across the frequency checks, one hundred and twenty-one uncertainties were over their limit and eleven remained unknown. The broader reference therefore remains inadequate under its stated criteria, even though its energy effective sample size passed. The reference exported no candidates. The next slide shows where much of the remaining uncertainty appears in residue choices.",
        ],
    },
    11: {
        'speaker': 'Daniel Sun',
        'paragraphs': [
            "A residue frequency is the fraction of sampled sequences that carry a particular amino acid at a particular position. Limited exploration leaves that estimate uncertain. Frequency uncertainty tells us how precisely our sampling estimates a residue's occurrence under the model.",
            "Look at phenylalanine, shown as F, at position sixty-five. Its estimated frequency uncertainty is zero point zero zero six eight three, while the declared limit is zero point zero zero two five. That is about two point seven times the limit and is the largest observed residue uncertainty in this report. Positions eighty-two and forty-four also contribute large values.",
            "There are also eleven pair-frequency uncertainties we cannot yet estimate. All involve phenylalanine at position eleven paired with choices at position fifteen; those combinations appeared in none of the chains. These remain unknown rather than measured zeros. These paired sites give us a concrete place to investigate coupled moves and coverage. To estimate their joint frequency, the chains need to visit those residue combinations.",
        ],
    },
    12: {
        'speaker': 'Daniel Sun',
        'paragraphs': [
            "We also made an engineering improvement while investigating the reference. The workflow used to build and save a candidate pool even for a reference run that selected no candidates. We removed that work while retaining traces and diagnostics.",
            "The chart compares saved artifacts for a replay of the earlier reference setup, with sixteen thousand three hundred and eighty-four observations per chain. Saved output went from about ninety-four point four seven megabytes to twelve point six four megabytes, a reduction of eighty-six point six percent. Every saved trace array and all diagnostic values in that controlled replay matched the earlier run exactly. The same evidence is now easier to store and inspect.",
            "This replay is separate from the larger run: it measures storage savings with matching sampling results. That leaves the scientific bottleneck in plain view: getting a qualified broad reference by exploring more of the coupled residue choices.",
        ],
    },
    13: {
        'speaker': 'Daniel Sun',
        'paragraphs': [
            "So what happens next? First, we should study the residue choices that the broad reference visits too rarely. The missing combinations at positions eleven and fifteen, together with uncertain common residues such as phenylalanine at position sixty-five, suggest tests of how our chains move between coupled choices. We can compare sampling settings and diagnostics while preserving the target, its allowed residue choices, and the declared uncertainty thresholds.",
            "Second, we need a qualified reference before making a reliable broad-target timing comparison. We already have a direct baseline and a corrected binary route. The next comparison should ask whether changes to the encoded proposal improve exploration at a competitive end-to-end cost, including the work needed to produce usable protein samples. For both routes, we should check frequency estimates and other accuracy diagnostics before comparing time. That way a faster number means we reached the same scientific goal sooner.",
            "Third, once the sampling evidence is strong enough to guide selection, candidates need biological review. A model-derived sequence is a starting point for considering sequence plausibility, practical constraints, and eventual experimental assays. Biological function awaits experimental testing. Those are successive decisions: qualify the reference, compare routes fairly, then review and test any sequences we choose to pursue.",
            "We can now locate the limiting evidence. Our small target has a qualified speed comparison; the broader target has an explicit list of uncertainty blockers; and our software produces traces we can revisit. Those results frame the next experiment: how do we explore the original distribution more efficiently?",
        ],
    },
    14: {
        'speaker': 'Daniel Sun',
        'paragraphs': [
            "Let me close by returning to the opening question. We built a software route from supported protein residue choices to binary energies, used a corrected Ising proposal, and kept the original protein distribution as the shared target. This lets us measure the value of binary sampling.",
            "The measured result today is direct sampling winning on the qualified small case. On the broader target, a larger direct run improved some diagnostics, while residue and pair-frequency uncertainty still leaves the reference inadequate. Our next milestone is better exploration and a qualified broad reference, followed by a fair end-to-end comparison. Biological review and experimental evidence come after model-based selection.",
            "Thank you for listening. Neil and I would be happy to take your questions.",
        ],
    },
}
