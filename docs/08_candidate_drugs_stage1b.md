# 08 — Stage 1b: candidate drugs, name resolution and the anti-target filter

Status: **machinery implemented and tested** (2026-09-27); the candidate list itself is not in the
repo yet. This document describes how a candidate list is consumed, how names become ChEMBL IDs, and
how the anti-target rules exclude candidates. It records no genotypes.

## The problem this stage answers

BUB1B is not a druggable target. Open Targets reports it `isEssential: true`, with no approved drug,
no advanced clinical candidate, no structure with a ligand, no high- or medium-quality pocket and no
chemical probes. In the EMBL-EBI Perturbation Catalogue it is a significant hit in 1,167 of 1,194
cell lines (~98%). Complete knockout in mouse is embryonic-lethal with full penetrance.

Worse, the drugged neighbourhood points the wrong way. Every drugged node near BubR1 in the mitotic
machinery — PLK1, AURKB, CDK1, TTK, KIF11, TOP2A, CENP-E, BIRC5 — is targeted by **inhibitors**. The
disease is loss of function; almost every available drug reduces function further. That asymmetry is
the central problem, so the pipeline encodes it as a rule rather than leaving it to prose.

## The candidate list

`assets/drugs_track2.template.yml` is the schema, with both lists empty. A real list is copied from
it and passed with `--drug_list_file`. The file is the single source for the run: the label, the
example flag, the axes, the candidates and the anti-target rules all come from it, so a report can
never show one list's verdicts under another list's name.

Candidates are named, never given as ChEMBL IDs. Names are readable and reviewable; IDs are derived
at run time against a stated Open Targets release, so the mapping is reproducible rather than pasted
in from a chat log.

Each candidate carries:

| Field | Purpose |
|---|---|
| `name` | resolved through the Open Targets search API |
| `axis` | which rationale proposes it; the report groups by this and never ranks into one list |
| `tier` | confidence, printed verbatim |
| `move_up` | what evidence would raise the tier — tiering without this is just an opinion |
| `patient_support` | how well the candidate is supported for **this child**: `released_feature`, `syndrome_only` or `mechanism_only`. Required; an unrecognised value stops the run |
| `investigational` | optional; printed next to the Open Targets maximum clinical stage |

`patient_support` is how the epistemic rule in docs/09 reaches the output. The report groups
candidates under three headings by this field, each carrying a sentence about what that group's
support amounts to, so an axis resting on the published series can never be read as one resting on
this child's own features. A candidate whose supporting feature is merely *absent from the released
HPO set* is `syndrome_only`, never `released_feature`, because absence from that set is not evidence
of absence in the child.

`pgx_annotate.nf` refuses to start when `candidates` is empty, when a candidate lacks `name`, `axis`
or `tier`, or when a candidate references an axis that is not defined; `resolve_candidates.py` fails
on a missing or misspelt `patient_support`. Nothing is assumed.

## Name resolution (`bin/resolve_candidates.py`)

| Status | Meaning |
|---|---|
| `resolved` | one drug entity matches the name case-insensitively, or one hit dominates (score ≥ 10 and ≥ 2× the runner-up) |
| `ambiguous` | several entities share the name, or no exact match and no dominant hit; the competing hits are listed |
| `not_found` | the search returned no drug entity |

An unresolved name is **reported as a row in the verdict table**, never dropped: a candidate that
silently disappears is indistinguishable from one that was never proposed.

## The anti-target filter

> **An inhibitor of a mitotic-fidelity gene is contraindicated in a mitotic-fidelity disorder.**

**Rule 1, automated.** A candidate is excluded when a mechanism-of-action row in Open Targets has an
inhibitory action type (`INHIBITOR`, `ANTAGONIST`, `BLOCKER`, `DEGRADER`, `DISRUPTING AGENT`,
`NEGATIVE ALLOSTERIC MODULATOR`, `NEGATIVE MODULATOR`) **and** one of that row's targets carries a
direct Gene Ontology biological-process annotation to one of:

| GO term | Label |
|---|---|
| GO:0007059 | chromosome segregation |
| GO:0000070 | mitotic sister chromatid segregation |
| GO:0000819 | sister chromatid segregation |
| GO:0098813 | nuclear chromosome segregation |
| GO:0007094 | mitotic spindle assembly checkpoint signaling |
| GO:0007052 | mitotic spindle organization |

The GO annotations are read from the Open Targets target record at run time, so the rule follows the
ontology rather than a hand-written gene list.

**Rule 2, curated.** Targets whose inhibition is shown in primary literature to compromise mitotic
fidelity but which carry none of those GO terms. Each entry cites a PMID that is printed in the
exclusion reason. Currently one entry: MTOR/FKBP1A, from Yamada C et al., *iScience*
2022;25:103675 (PMID 35141499) — mTORC1 inactivation induces mitotic slippage, overriding spindle
assembly checkpoint arrest, and causes chromosome instability in yeast and human cells. Verified
against the abstract before being added.

**Not assessable.** A drug with no mechanism-of-action row cannot be tested by either rule. It is
marked `not assessable`, carried forward, and the caveat is printed: absence of a mechanism
annotation is not evidence of safety.

Excluded candidates stay in every table with their exclusion reason. **The exclusions are the
evidence that the filter works**, so removing them would remove the result.

### Verified behaviour

Run against the live API (2026-09-27, Open Targets 26.09 / API 26.9.0):

| Drug | Mechanism in Open Targets | Verdict |
|---|---|---|
| volasertib | INHIBITOR of PLK1 | excluded — GO:0000070, GO:0007052, GO:0007094 |
| etoposide | INHIBITOR of TOP2A, TOP2B | excluded — GO:0000819, GO:0007059 |
| sirolimus | INHIBITOR of FKBP1A | excluded — curated rule, PMID 35141499 |
| tolvaptan | ANTAGONIST of AVPR2 | allowed |
| (a nonsense string) | — | `not_found`, reported |

`tests/unit/test_resolve_candidates.py` covers the rule logic offline from a fixture (11 tests). The
cases that would otherwise fail silently: an **activator** of the same segregation gene is allowed
(the rule is about inhibition, not about the target); a target caught only by the curated rule
carries its PMID; a drug with no mechanism is `not assessable` rather than `allowed`; ambiguous and
missing names are reported and not queried; and an empty `anti_target_rules` block excludes nothing
rather than appearing to have filtered.

## Outputs

```
results/pgx/open_targets/<sample>.track2.resolved.tsv     name, axis, tier, ChEMBL ID, matched name,
                                                          score, status, stage, mechanism targets,
                                                          anti-target verdict and reason
results/pgx/open_targets/<sample>.track2.anti_target.tsv  one row per drug x mechanism target:
                                                          action type, rule matched, verdict, reason
results/pgx/open_targets/<sample>.pgx_ot.verdicts.tsv     drug | axis | tier | patient_support |
                                                          gene | called | matched OT allele |
                                                          evidence level | flag | guideline text
results/pgx/open_targets/<sample>.pgx_ot.summary.json     axes, counts, exclusions, unresolved,
                                                          join coverage, flag counts
```

Only evidence levels 1A and 1B are asserted; level-3 rows are counted as exploratory and carried no
further. The report groups axes by `patient_support` first, then renders one subsection per axis with its
rationale, the per-candidate tier and `move_up`, and the verdict table; exclusions and unqueried
candidates follow in their own blocks.

## The anti-target rule and treatment the child is already receiving

The rule excludes topoisomerase II poisons, which are components of standard rhabdomyosarcoma
therapy. That tension is real and the report states it rather than leaving it unremarked: the rule
governs **proposing new drugs** for a mitotic-fidelity disorder. It is not applied to treatment the
child is already receiving, and it does not license withdrawing a therapy given with curative
intent. What the overlap does support is a risk consideration — alkylators and topo-II poisons are
the established causes of therapy-related myeloid neoplasia, and this is a chromosomal instability
syndrome — which the report raises for cumulative-dose decisions and long-term haematological
surveillance, not as a protocol change.

## Framing that the report states, and why

- Everything produced here is a **research hypothesis, not a treatment recommendation**.
- **Open Targets sources its pharmacogenomics from ClinPGx (formerly PharmGKB)** — the same knowledge
  base CPIC and PharmCAT draw on. Agreement between the Open Targets layer and PharmCAT is therefore
  *not* independent corroboration, and the report says so rather than letting a reader infer two
  witnesses where there is one.

## Running it

```sh
cp assets/drugs_track2.template.yml assets/drugs_track2.yml   # then fill in axes and candidates
CONDA_SUBDIR=osx-arm64 nextflow run pgx_annotate.nf -profile conda,macos_docker \
    --drug_list_file assets/drugs_track2.yml
```

Stage 1a is not re-run. The main pipeline's REPORT step then picks the stage 1b output up by file
name on its next run.
