# Week 4 adviser review packet

The technical benchmark proposal and updated plan are ready for review. No meeting, external message, signature or adviser approval has been recorded.

Read the [project plan](../project_plan.md), [benchmark contract](benchmark_contract.md), [admission evidence](benchmark_evidence.md), and [Week 3 reuse audit](../week03/research_audit.md). Runnable commands are in the [project guide](../../README.md).

## Concrete decisions

| Decision | Proposed choice | Review needed |
| --- | --- | --- |
| Biological context | CM, Q5SLA5 residues 2–90, exact coordinate inference and documented HB8 strain context | Biology lead and adviser |
| Fixed probability model | Archived feedback-informed lambda-one checkpoint, original bytes and exact forbidden-state support | Both leads and adviser |
| Coordinate versus sampled state | Preserve native anchor with infinite model energy; use a separately labelled supported computational witness | Both leads |
| Restrictions and split | Two development and two evaluation cases in the contract; no independent assay claims | Biology lead and adviser |
| Accuracy and cost | Configured tolerances, four chains, equal budgets and end-to-end timing | Both leads and adviser |
| Reuse boundary | Reuse audited encoding/model primitives; derive the actual proposal correction before Week 8 review | Systems lead |

## Independent score worksheet

Use the original parameter member and archive record `sequence_1 functional_true`. A software inspection obtained:

```text
Selected fields: 96 terms, sum = -45.43422150
Selected canonical i<j pairs: 4,560 terms, sum = 102.792881272664
E = -sum(h) - sum(J)
  = -(-45.43422150) - 102.792881272664
  = -57.358659772664
```

Decimal accumulation and binary64 selected-state accumulation agree to the recorded tolerance. This is one automated inspection, not two human attestations. Each lead should independently select the original terms, verify the sign and pair convention, and record the result. Do not substitute a score of the constructed witness or coordinate anchor.

For numbering, inspect the [96-row map](../../data/benchmark/position_map.csv): model 8 maps to accession residue 7; model 12 to 11; model 16 to 15; model 24 to 23. The seven gap columns have no accession residue. Confirm that the accession is 355 residues and the matched domain covers 2–90. The source sequence is not silently relabelled as a full-length protein.

## Review record

| Required record | Current status |
| --- | --- |
| Systems lead independent score and numbering check | Pending human record |
| Biology lead independent score and numbering check | Pending human record |
| Agreement on reference context and allowed residues | Pending |
| Adviser decision on model, targets and split | Pending |
| Requested changes, owner and due milestone | No feedback supplied |

After review, retain the decision and update the contract revision explicitly. A changed model, support, restriction case or accuracy target requires a new configuration revision and affected checks to be rerun. Week 5 can implement the proposed support-aware loader, but final agreement must not be inferred from passing automated tests.
