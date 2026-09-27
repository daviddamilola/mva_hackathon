# 09 — Clinical features, and what may be inferred from them

Status: **implemented in the report** (2026-09-27). This document records the rule that governs every
phenotype claim in this project, why it exists, and how it is enforced mechanically rather than by
good intentions. It records no genotypes and no clinical text.

## The released feature set

Eight HPO terms were released to hackathon participants as the clinical features leading to
diagnosis:

| Feature | HPO |
|---|---|
| Rhabdomyosarcoma | HP:0002859 |
| Nephrocalcinosis | HP:0000121 |
| Short stature | HP:0004322 |
| Failure to thrive | HP:0001508 |
| Skeletal muscle atrophy | HP:0003202 |
| Premature birth | HP:0001622 |
| Small for gestational age | HP:0001518 |
| Recurrent spontaneous abortion (parental) | HP:0200067 |

## The rule

> **This set is the features made available to participants. It is not the complete clinical
> picture.**
>
> - A feature's **presence** in the set is evidence that it is present in the child.
> - A feature's **absence** from the set is **not** evidence that it is absent. It may be absent, it
>   may not have been assessed, or it may simply not have been released.

So the report never writes that this child "lacks microcephaly", "has no seizures" or "has normal
development". Where a feature is not in the set, the report says it is **not among the features
released to participants**, which is a statement about the data rather than about the child.

Every inference that depends on a feature being genuinely absent is written as conditional, states
that dependency, and does not assert the antecedent.

### Why this matters more than it looks

The asymmetry is easy to lose. Two of the strongest-sounding inferences available here run through
it:

1. **Residual function.** Simmons et al. (PMID 30668728) show nearly complete loss of BubR1 is
   required for microcephaly. Microcephaly is not in the released set. It is tempting to conclude
   the missense allele retains meaningful activity — but that conclusion needs microcephaly to be
   genuinely absent, which the released set cannot establish. The report states the inference
   conditionally and flags two further cautions: microcephaly is explicitly **not required** for an
   MVA diagnosis (Callier P et al., Am J Med Genet A 2005;137:204–7, PMID 16059936), and the child
   is young, so a head-circumference trajectory may not be complete.

2. **Seizures.** They occur in 45.1% of published MVA1 cases and are not in the released set.
   Antiepileptic pharmacogenetics is therefore presented as contingency — applicable if seizures are
   present or emerge — and never as evidence that this child does not have them.

A reviewer who knows the field will check exactly these two places.

## How it is enforced

Not by remembering. Three mechanical points:

- **`assets/report/clinical_contexts.yml`** gives every clinical context a `patient_support` value of
  `released_feature`, `syndrome_only` or `mechanism_only`. The report prints the corresponding
  sentence under each heading, so a syndrome-level context cannot silently read as patient-specific.
  The summary's "drug classes this child is likely to meet" list is built from `released_feature`
  contexts only.
- **The stage 1b candidate list** carries the same field per candidate. `bin/resolve_candidates.py`
  **fails the run** on a missing or unrecognised value rather than defaulting, and the report groups
  candidates under three headings by that field, each with a sentence saying what the group's
  support amounts to. `tests/unit/test_resolve_candidates.py::test_bad_patient_support_is_rejected`
  covers the failure.
- **Methods and Limitations** both state the rule in full, so a reader who starts anywhere in the
  report meets it.

## What would change with the full record

If the complete phenotype were available, three things would resolve:

1. whether the residual-function inference above holds;
2. whether the antiepileptic axis is contingency or current care;
3. whether the nephrocalcinosis is syndromic or a consequence of the neonatal course — this child
   was born at 32 weeks at about 1 kg, and nephrocalcinosis is common in preterm, very-low-birth-
   weight infants (loop diuretics, methylxanthines, vitamin D, parenteral calcium and phosphate,
   immature tubular handling). "Present since birth" argues somewhat against a NICU-acquired
   origin, but the timing within the first days and persistence beyond infancy are what would
   distinguish the two, and neither is available here.

Each is recorded in the report as an open question rather than resolved in whichever direction
flatters the analysis.
