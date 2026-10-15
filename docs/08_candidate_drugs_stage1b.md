# 08 — Stage 1b: candidate drugs, name resolution and the anti-target filter

Status: **machinery implemented and tested** (2026-09-27); **regimen list run** (2026-10-06);
**candidate list v2 run** (2026-10-10). `assets/drugs_track2.yml` is the Track 2 list the report is
built from: two mechanism-directed axes (readthrough, proteostasis), the tumour-directed axis, and
the regimen axis. `assets/drugs_regimen.yml` is the regimen axis alone and is kept
for a run without proposals. This document describes how a candidate list is consumed, how names
become ChEMBL IDs, how the anti-target rules exclude candidates, how a regimen axis differs, and
what the first candidate axis rests on. It records no genotypes.

## The problem this stage answers

Track 2 asks for mechanism and drug repurposing. For this gene the repurposing framing is forced
rather than chosen, and the report now says so with numbers fetched at run time
(`ot_query.py tractability`, rendered as report section 2) rather than asserted in prose.

Open Targets records BUB1B as `isEssential: true` with **no approved drug, no advanced clinical
candidate and no Phase 1 candidate in any modality**. The detail is more interesting than a flat
"undruggable", and the report prints the whole table:

| Modality | Evidence present | Notable absences |
|---|---|---|
| Small molecule | High-Quality Ligand, Druggable Family | Structure with Ligand, High-/Med-Quality Pocket |
| Antibody | none | all nine assessed labels negative |
| PROTAC / degrader | Literature, UniProt Ubiquitination, Database Ubiquitination | Small Molecule Binder, Half-life Data |
| Other clinical modality | none | all three assessed labels negative |

Two things follow. Small-molecule tractability is not zero — there is a high-quality ligand and it is
a druggable family — but no pocket and no structure with a ligand are recorded, and nothing has been
carried forward, so there is nothing to repurpose *against the target itself*. And the one modality
with positive evidence beyond ligandability is **degradation**, the inverse of what a loss-of-function
disease needs: a PROTAC against BubR1 would remove the protein this child is short of. That inversion
is what makes the protein-stabilisation axis worth stating, and it is why no axis proposes acting on
BUB1B directly.

Complete knockout in mouse is embryonic-lethal with full penetrance, and in the EMBL-EBI Perturbation
Catalogue BUB1B is a significant hit in 1,167 of 1,194 cell lines (~98%), so the target is not one to
inhibit further even if it could be drugged.

The drugged neighbourhood points the wrong way too. Every drugged node near BubR1 in the mitotic
machinery — PLK1, AURKB, CDK1, TTK, KIF11, TOP2A, CENP-E, BIRC5 — is targeted by **inhibitors**,
developed to reduce mitotic fidelity in dividing tumour cells. This disease is a loss of mitotic
fidelity; almost every drug near the pathway would deepen the defect. That asymmetry is the central
problem, so the pipeline encodes it as a rule rather than leaving it to prose.

**The strategy is therefore compensation, not replacement.** No molecule restores BubR1 kinase
activity. A repurposing hypothesis here acts on something other than the broken protein so the
consequences of its loss are reduced: restoring some full-length protein from the nonsense allele,
stabilising the partially functional protein the missense allele makes, reducing the cellular cost of
the resulting aneuploidy, or treating downstream consequences that are druggable where the spindle
checkpoint is not. Each axis in the candidate list is one such point of leverage.

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
or `tier` (tier is optional only on a regimen axis), or when a candidate references an axis that is
not defined; `resolve_candidates.py` fails on a missing or misspelt `patient_support`, and on an axis
`role` it does not know. Nothing is assumed.

### Axis roles: proposals versus treatment already given

An axis carries `role: repurposing` (the default) or `role: regimen`. The distinction exists because
the two kinds of drug need the join to run in opposite directions of trust:

| | Repurposing axis | Regimen axis |
|---|---|---|
| What the drugs are | new proposals for this disorder | standard therapy the child is already receiving |
| Anti-target rule | **excludes**; an excluded candidate is not queried | evaluated and reported as a **mechanistic note**; the drug is queried regardless |
| Tier, `move_up` | required | not used (they are not proposals) |
| Where the report shows them | by axis, grouped by `patient_support` | one block, "Treatment the child is already receiving" |
| What a rule match means | the filter works | the basis of the therapy-related-risk consideration |

The rule governs **proposing new drugs** for a mitotic-fidelity disorder. It is not applied to
therapy given with curative intent, and it does not license withdrawing it (see the section on
treatment the child is already receiving, below). Querying a regimen drug anyway is the whole point:
its ClinPGx annotations, including positions outside the stage 1a panel, are what the report's
"unannotated" line promises. The resolver records `role` on every row, counts `regimen` and
`regimen_anti_target_hits` separately from `excluded_by_anti_target` (which counts repurposing
candidates only), and lists regimen matches under `regimen_anti_target` rather than `excluded`. The
join adds no "excluded" placeholder row for a regimen drug, because it has real rows.

The regimen list is a stand-alone file today so it can run before the candidate set exists. When the
Track 2 candidates are written, the `regimen` axis is copied into that file, so one run produces one
report with both blocks and the report can never show one list's verdicts under another's label.

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

`tests/unit/test_resolve_candidates.py` covers the rule logic offline from a fixture (13 tests). The
cases that would otherwise fail silently: an **activator** of the same segregation gene is allowed
(the rule is about inhibition, not about the target); a target caught only by the curated rule
carries its PMID; a drug with no mechanism is `not assessable` rather than `allowed`; ambiguous and
missing names are reported and not queried; an empty `anti_target_rules` block excludes nothing
rather than appearing to have filtered; a regimen-axis drug the rule matches is still queried, with
the match recorded and counted apart from the exclusions; and an unknown axis role stops the run.
`tests/unit/test_pgx_join.py` checks that a matched regimen drug keeps its real rows and gets no
"excluded" placeholder, while a matched repurposing candidate still does.

## Outputs

```
results/pgx/open_targets/<sample>.track2.resolved.tsv     name, axis, tier, ChEMBL ID, matched name,
                                                          score, status, stage, mechanism targets,
                                                          anti-target verdict and reason
results/pgx/open_targets/<sample>.track2.anti_target.tsv  one row per drug x mechanism target:
                                                          action type, rule matched, verdict, reason
results/pgx/open_targets/<sample>.pgx_ot.verdicts.tsv     drug | axis | role | tier | patient_support |
                                                          gene | called | matched OT allele |
                                                          evidence level | flag | guideline text
results/pgx/open_targets/<sample>.pgx_ot.summary.json     axes (with role), counts, exclusions,
                                                          regimen_anti_target, unresolved,
                                                          join coverage, flag counts
```

Only evidence levels 1A and 1B are asserted; rows at any other level are counted per drug, with the
genes behind level-2 rows named in `evidence_levels.tsv` and the summary, and carried no further. The report groups axes by `patient_support` first, then renders one subsection per axis with its
rationale, the per-candidate tier and `move_up`, and the verdict table; exclusions and unqueried
candidates follow in their own blocks.

## Only approved medicines count

The brief asks for *existing, approved* medications. The candidate list therefore holds approved
medicines only. An investigational agent may be named in an axis rationale where it is the direction
the field is moving (ELX-02 on the readthrough axis, for instance), but it is not a candidate, is not
counted, and does not appear in the verdict table. Where an approval is narrower than it looks, the
candidate carries a `caution` field saying so rather than a bare "approved": ataluren's conditional
EU authorisation was not renewed in March 2025 and it has never been FDA-approved, so it is listed as
approved in Great Britain and the controversy is part of the entry.

## The two mechanism-directed axes (candidate list v2, 2026-10-10)

The brief's rigor criterion is conjunctive: the mechanism must be sound **and** the candidate must be
supported by that mechanism. A candidate supported by tumour biology meets the first half only. The
list is therefore ordered by an explicit `order` key on each axis, mechanism-directed first, and the
report follows that order rather than sorting axes alphabetically or by support level.

**Axis `readthrough` (order 1, tier 2).** One allele is a nonsense variant, and the pipeline
characterises it rather than assuming it is druggable (doc 04). Two findings point in opposite
directions and both are reported: the stop codon is UGA, the most permissive, while the +4 base is A,
the least permissive in the published hierarchy. All three near-cognate products are predicted
tolerated, and the residue sits in an inter-domain linker. The PTC is NMD-sensitive, which is the
reason to pair readthrough with NMD inhibition rather than an argument against the axis. Candidates:
ataluren, gentamicin, amikacin, amlexanox. The aminoglycosides carry two gates pointing different
ways — the MT-RNR1 pharmacogenomic gate is open because the genotype is reference, while
dose-dependent nephrotoxicity against a background of nephrocalcinosis is a patient-specific caution.

**Axis `proteostasis` (order 2, tier 3).** Acts on the consequence of BubR1 loss rather than on
BubR1: aneuploid cells carry stoichiometric imbalance that burdens chaperone and proteasome systems.
This is also why pathway enrichment over the MVA gene set finds nothing — the convergence is on the
consequence, not on an annotated pathway. Candidates: arimoclomol, glycerol phenylbutyrate. Tier 3,
because there is no evidence in MVA or any aneuploidy syndrome and no biomarker in this child.

**Axis 1, tumour-directed (order 3, all tier 3).** Demoted from tier 2. Two reasons, stated once in
the axis rationale rather than five times in the candidates: every one of them is selected on a
tumour genotype that germline sequencing cannot provide, and they propose intensifying a regimen the
report separately argues may warrant caution in this genotype. They are kept, because a demoted
candidate with its reason recorded is evidence the framework works.

## Axis 1: partners for irinotecan (candidate list v1)

The first proposed axis came from a supplied table of five drugs said to act synergistically with
irinotecan. Each claim was checked against primary literature before being written into the list,
and the tiers record what was found rather than what was claimed:

| Candidate | Claim supplied | What was found | Tier |
|---|---|---|---|
| olaparib | exploits "BUB1-induced repair defects" | No BUB1B-PARP literature. A paediatric phase I/II of olaparib + irinotecan (AcSé-ESMART arm D, 70 patients) defined a dose and found activity enriched in high-aneuploidy tumours; additive-to-synergistic with irinotecan in paediatric models including rhabdomyosarcoma | 2 |
| niraparib | as olaparib | PARP-class rationale only; no rhabdomyosarcoma data | 3 |
| abemaciclib | "clinically proven durable complete responses in RMS" with irinotecan/temozolomide | A paediatric phase 1b (JPCS Part A) with one complete response at the recommended dose, and one case report of a 22.6-month complete response in alveolar rhabdomyosarcoma. Real, but one phase 1b response and one case, not "clinically proven" | 2 |
| ceralasertib | strips "BUB1-deficient" cells of a backup checkpoint | No BUB1B-ATR literature. Adult phase 1 of an ATR inhibitor + irinotecan with activity in ATM-deficient tumours; ceralasertib's paediatric data are with olaparib, not irinotecan | 3, investigational |
| simvastatin | "distorts spindle dynamics in BUB1-deficient backgrounds" | Nothing found. Kept so the record shows the claim was checked, and because the SLCO1B1 join is informative | 3 |

Two corrections applied to the supplied rationale. The table referred to BUB1; this child's gene is
BUB1B, a different kinase, and no finding about one transfers to the other without evidence. And
every mechanism of the form "exploit the deficiency in tumour cells" meets the fact that the
variants are germline: whatever the drug does to BubR1-deficient tumour cells it does to every
other cell, so no selective window is claimed. The axis therefore stands on the published
combination evidence in paediatric tumours and on one observation that bears on this disease
directly, the aneuploidy enrichment in the olaparib trial. It treats the malignancy, not the
syndrome, and says so.

All five passed the anti-target rule (none has an inhibitory mechanism against a mitotic-fidelity
target). None has a ClinPGx annotation at an asserted level except simvastatin, whose SLCO1B1 rows
matched the PharmCAT call with no flag.

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
# the regimen list that is in the repo
CONDA_SUBDIR=osx-arm64 nextflow -log results/nextflow_stage1b.log run pgx_annotate.nf \
    -profile conda,macos_docker --drug_list_file assets/drugs_regimen.yml -ansi-log false

# the Track 2 candidate set, once written (copy the regimen axis into it)
cp assets/drugs_track2.template.yml assets/drugs_track2.yml   # then fill in axes and candidates
CONDA_SUBDIR=osx-arm64 nextflow -log results/nextflow_stage1b.log run pgx_annotate.nf \
    -profile conda,macos_docker --drug_list_file assets/drugs_track2.yml -ansi-log false
```

Use a separate `-log` file so the stage 1b session does not sit in the main pipeline's log, and
resume the main pipeline afterwards with its own session ID, never a bare `-resume`.

Stage 1a is not re-run. The main pipeline's REPORT step then picks the stage 1b output up by file
name on its next run.
