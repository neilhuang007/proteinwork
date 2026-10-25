# Technical handoff reference for the biology and systems leads

Use this during joint implementation and review. The [beginner guide](../biology_lead_guide.md) and [worksheet](benchmark_brief_template.md) are the starting point for the biology lead. The systems lead converts the worksheet into the structured records below and presents readable results for review.

## Current model and coordinate details

The proposed benchmark is cm_lambda1_supported_v1. The model has 96 alignment columns covering an 89-residue CM segment, Q5SLA5 residues 2–90, within a 355-residue protein. Alignment gaps have no biological residue number. Computer model indices start at zero; accession positions start at one. Use the existing [position map](../../data/benchmark/position_map.csv), not a newly guessed conversion.

The coordinate anchor is outside the model's permitted choices at two positions and has infinite model energy. Use it for numbering and distance comparisons, not sampling initialization or a finite mutation-energy baseline. The separate supported witness is a constructed computer test input, not an experimentally validated design. Current restriction examples are numerical tests without an established functional rationale. See the [benchmark contract](../week04/benchmark_contract.md).

## Turn the worksheet into structured records

- [position_review.csv](position_review.csv): the systems lead populates computer indices from the existing map. Together, verify requested residues, boundaries and gaps against the biological source. Leave accession positions blank for gaps and explain the reason. Record pending, checked or discrepancy.
- [restriction_requests.csv](restriction_requests.csv): translate each worksheet request into one row with the exact permitted letters, source and reason. Distinguish biological requests from numerical tests. A single letter fixes the site. Omitted positions follow the existing case defaults; they do not automatically become unrestricted. The systems lead checks model support and conflicts before accepting requests.
- [assay_inventory.csv](assay_inventory.csv): biology supplies sources, measured constructs, experimental readouts, units, context and known use in model development. Systems supplies technical data locations, exact sequence matches and overlap checks. Keep uncertain fields unknown and eligibility pending. A measured label belongs to its measured sequence and context; do not transfer it from a nearby sequence.
- [review_log.md](review_log.md): the systems lead may transfer dated decisions from the worksheet into this structured log. Record actual decisions and participants, preserving unresolved questions.

CSV means a plain table; the files open in a spreadsheet. Blank templates are not approved benchmark inputs. Both leads review the translated settings before they become the agreed experiment.

## Prepare the independent score check together

This is separate from reviewing a supported candidate export. Use the unmodified archive record Natural_CM.faa / sequence_1 functional_true specified in the [review packet](../week04/review_packet.md). It is a finite arithmetic fixture, not an initialization for the restricted cases. Its header alone is not an independent assay record for this project.

The systems lead prepares a spreadsheet containing its 96 selected field terms and 4,560 selected pair terms, with source locations and residue labels. Explain that a field scores a single residue and a pair term scores a combination at two positions. Show how to check selected entries against the original file and how each pair is counted once.

The biology lead independently checks selected terms and sums the two columns. The model convention is E = minus the field sum minus the pair sum. Record the result and any disagreement. This is a later guided arithmetic review, not a prerequisite for understanding the project or a request to inspect millions of unused parameters.

## Explain these terms when they arise

| Term | Meaning for this project |
| --- | --- |
| Support | Choices the model gives nonzero probability; exclusion does not prove biological impossibility. |
| Target distribution | The intended probabilities for the allowed sequences. |
| Energy | The statistical score used to set those probabilities, not measured protein stability. |
| Chain | A sequence of dependent sampling steps. Repeated states must remain in its raw record. |
| Correction | The encoded method's mathematical accept/reject step for preserving the target. |
| Effective sample size | How much independent information correlated draws provide. |
| Holdout and leakage | Evidence reserved from development; leakage reuses development information as an apparently independent test. |
| Checkpoint and provenance | A saved trained model, and the record of where files and claims came from. |
| FASTA | A text format with a sequence name followed by amino-acid letters. |

For code, mathematical assumptions, test criteria and the full schedule, use the [implementation roadmap](../implementation_roadmap.md). Human reviews remain pending until recorded; this explanatory document does not change the benchmark or its settings.
