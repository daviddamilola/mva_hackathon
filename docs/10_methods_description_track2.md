# 10 — Methods description, Track 2 submission form

Status: **drafted** (2026-10-15). Paste-ready answers for the Track 2 methods-description
spreadsheet, one block per row of the template. The submitted report is
`david-oluwasusi_track2_report.pdf` (rendered by the pipeline; see docs/07). This document records
no genotypes and no clinical text beyond what docs/08 and docs/09 already record.

---

## Team name

David Oluwasusi

## Please describe your approach in detail (how you went from variant/mechanism to candidate medication(s))

The analysis is a Nextflow DSL2 pipeline (every tool version-pinned, every script digest part of its
task hash) that goes from the proband's WGS to a single report in which every number is traced to
the step that produced it.

**Variant to mechanism.** Reads were re-aligned to a rebuilt, contig-verified GRCh38. The causal
BUB1B variants were characterised against a published protein domain map; nonsense-mediated decay
was predicted from the MANE Select transcript by the 55-nt rule, with the distance reported rather
than the verdict alone. For the nonsense allele, the pipeline then characterised what a readthrough
agent would actually produce: the stop codon and its 3′ context were re-read from the genome and
required to match the transcript annotation, the near-cognate decoding products were taken from the
published set for that stop codon, and each product was scored with AlphaMissense (the
all-substitution release, streamed and filtered to the one UniProt accession), alongside every
substitution at the residue and per-domain means for context.

**Mechanism to strategy.** Open Targets records BUB1B as an essential gene with no approved drug and
no clinical candidate in any modality; the one modality with positive evidence beyond ligandability
is degradation — the inverse of what a loss-of-function disease needs. Nearly every drugged node in
the mitotic machinery is an inhibitor, developed to reduce mitotic fidelity. The repurposing framing
is therefore compensation, not replacement: no candidate acts on BUB1B directly, and each axis is a
point of leverage around the loss (readthrough of the nonsense allele, proteostasis support against
the cost of aneuploidy, tumour-directed combinations, and the regimen the child is already
receiving).

**Strategy to candidates.** Candidates live in a reviewed YAML list carrying, per candidate, the
axis, a confidence tier, what evidence would raise the tier, and an explicit `patient_support` field
(`released_feature`, `syndrome_only` or `mechanism_only`) that the report groups by, so
syndrome-level support can never read as patient-specific. Names are resolved to ChEMBL IDs at run
time against a stated Open Targets release. An anti-target rule written as code — a drug is excluded
when an Open Targets mechanism-of-action row has an inhibitory action type and a target annotated to
mitotic-fidelity GO terms, plus one literature-curated entry carrying its PMID — excludes any
proposal that would deepen the mitotic defect, and the exclusions are reported as evidence the
filter works. Surviving candidates were joined to the child's pharmacogene calls (PharmCAT with
three reconciled callers, stage 1a; ClinPGx annotations via Open Targets at evidence levels 1A/1B,
stage 1b), so each candidate carries its pharmacogenomic verdict for this child.

## (required) Generative AI

Anthropic Claude Code (Max plan), Claude Opus 5 — commercial terms, no training on customer content.

## Was candidate identification automated or primarily manual literature review and expert curation?

Both, in sequence. Automated: the target assessment (tractability, essentiality, drugged
neighbourhood), name-to-ChEMBL resolution, the GO-based anti-target filter and the pharmacogenomic
join all run as code against stated database releases, and their outputs are tables in the report.
Manual: the axes themselves and every candidate's rationale were curated, and each supplied or
AI-suggested claim was verified against primary literature before being written into the candidate
list, with tiers recording what was found rather than what was claimed.

## Please describe any manual literature review or expert curation in detail

Every mechanistic claim behind a candidate was checked against primary sources (PubMed) before
entering the list. Examples recorded in the repository: a supplied table of five irinotecan
partners was re-verified claim by claim — the olaparib claim ("BUB1-induced repair defects") had no
BUB1B literature, but a paediatric phase I/II (AcSé-ESMART arm D) with activity enriched in
high-aneuploidy tumours was found and cited instead; an abemaciclib claim of "clinically proven
durable complete responses" was downgraded to what the record shows (one phase 1b response and one
case report); two claims (ceralasertib/BUB1B, simvastatin/BUB1B) had no supporting literature and
were tiered accordingly, with the record kept to show the claim was checked. The curated anti-target
entry (MTOR/FKBP1A) was verified against the abstract of PMID 35141499 before being added.
Readthrough biology (near-cognate products per stop codon, the +4-base hierarchy) was taken from the
published literature and cited in the report. Regulatory status was checked by jurisdiction:
ataluren's EU conditional authorisation was not renewed in March 2025 and it has never been
FDA-approved, so it is listed as approved in Great Britain with the controversy stated. A supplied
reference to BUB1 was corrected to BUB1B, a different kinase.

## Publicly available databases only, or proprietary sources?

Publicly available databases and published literature only. No proprietary data source was used.
The only access-controlled input was the hackathon's own subject dataset (gated Hugging Face
download under the DUA), which is the WGS and released HPO terms, not a drug/target source.

## Public data sources used in detail

- **Open Targets Platform 26.09 (API 26.9.0)** — target tractability and essentiality, mechanisms of
  action, GO annotations for the anti-target rule, drug search for ChEMBL ID resolution, and
  pharmacogenomics (sourced from ClinPGx, formerly PharmGKB; evidence levels 1A/1B asserted, others
  counted). The report states that agreement between this layer and PharmCAT is not independent
  corroboration, since both draw on the same knowledge base.
- **PharmCAT** (CPIC, DPWG and FDA sources), **Aldy** and **PyPGx** for pharmacogene calls;
  **CPIC API** snapshot (sha256 recorded) for MT-RNR1 allele definitions; **IPD-IMGT/HLA 3.65.0**
  via T1K and OptiType for HLA.
- **Ensembl VEP** (cache 116, MANE Select), **gnomAD v4.1** joint frequencies, **ClinVar**
  (release 20260923), **AlphaMissense** (hg38 SNV file and all-substitution release), **SpliceAI**
  for the modifier panel.
- **EMBL-EBI Perturbation Catalogue** and mouse-knockout data for essentiality context; a published
  BUB1B domain map (citation in the asset header); **PubMed** primary literature for every curated
  claim.
- **ClinicalTrials.gov / trial literature** for the paediatric combination evidence cited on the
  tumour-directed axis.
- Reference and QC resources: NCBI GRCh38 analysis set, GRC exclusion regions, UCSC segmental
  duplications, ENCODE blacklist, Umap mappability, 1000 Genomes (somalier ancestry, VerifyBamID2
  contamination).

## Proprietary data sources used in detail

None.

## How did you characterize the variant's mechanism?

As loss of function, compound heterozygous, with the two alleles characterised separately rather
than assumed equivalent. The nonsense allele: premature stop placed on the MANE Select transcript
structure, predicted NMD-sensitive by the 55-nt rule (distance and exon reported); stop-codon
identity and 3′ context read from the genome and checked against the annotation; the three
near-cognate readthrough products scored with AlphaMissense at the residue, in an inter-domain
linker, with per-domain score context — so the readthrough axis rests on what the drug would
actually produce. The missense allele: placed on the published domain map, with a conditional
residual-function inference stated under the feature-absence rule (absence of a feature from the
released HPO set is not evidence of absence in the child). Downstream consequence: mosaic
aneuploidy measured genome-wide from depth and allele balance with a per-arm detection limit, and
the proteostasis axis framed on the consequence of aneuploidy rather than on an annotated pathway.
Gene-level context (essential, embryonic-lethal knockout, ~98% of cell lines dependent) is why no
axis proposes inhibiting or degrading the target.

## Estimate of time or effort

Approximately three weeks of part-time, single-person effort (25 September – 15 October 2026):
roughly one week on data verification, reference rebuild and alignment; one week on the
pharmacogenomic and mechanism branches; one week on candidate curation, literature verification and
the report. The pipeline re-runs end-to-end from cached inputs in hours on one workstation.

## Method abstract (up to 500 words)

Mosaic variegated aneuploidy from biallelic BUB1B loss has no approved therapy directed at its
cause, and its target is undruggable in the strict sense: Open Targets records an essential gene
with no clinical candidate in any modality, in a pathway whose every drugged neighbour is an
inhibitor of the process this disease already lacks. Our method takes that asymmetry as the central
design constraint. A Nextflow pipeline re-aligns the proband's WGS to a contig-verified GRCh38 and
characterises mechanism before proposing anything: domain placement of both alleles, NMD prediction
from the MANE Select transcript with distances reported, and — for the nonsense allele — an explicit
model of what readthrough would produce (stop-codon identity and 3′ context re-read from the
genome, near-cognate products from the published set, each scored with AlphaMissense against the
residue's full substitution profile). Candidate medications are then proposed only as compensation:
readthrough of the nonsense allele, proteostasis support against the cost of aneuploidy,
tumour-directed combinations with published paediatric evidence, and the regimen already being
given. Candidates are named in a reviewed list, resolved to ChEMBL at a stated Open Targets
release, and filtered by an anti-target rule written as code: any inhibitory mechanism against a
target with mitotic-fidelity GO annotations is excluded, plus one literature-curated entry with its
PMID, and the exclusions are printed as evidence the filter works. Survivors are joined to the
child's own pharmacogene calls (PharmCAT over force-called positions with depth and quality gates,
reconciled with Aldy and PyPGx; HLA from two callers passed only on agreement; MT-RNR1 called from
heteroplasmy fractions), so each proposal carries its pharmacogenomic verdict for this child at
ClinPGx evidence levels 1A/1B.

Strengths. The mechanism work is specific enough to falsify its own axis: stop codon and +4 base
point in opposite directions and both are reported. Epistemic discipline is enforced mechanically —
every candidate and clinical context carries a `patient_support` field, the report groups by it,
and a missing value fails the run, so syndrome-level evidence cannot read as patient-specific; a
stated feature-absence rule keeps the released HPO set from being mistaken for a complete
phenotype. Every claim supplied to us or suggested by AI tooling was verified against primary
literature before entering the list, and two supplied claims were corrected or downgraded in
writing. The pipeline is reproducible: versions pinned, inputs checksummed, script digests in task
hashes, and a provenance record accompanying the report.

Limitations. All candidates are research hypotheses; none has evidence in MVA itself, and the
mechanism-directed axes rest on predicted, not measured, readthrough and proteostasis effects.
The released HPO terms bound what can be claimed about this child, and two key inferences (residual
function, antiepileptic contingency) are stated conditionally for that reason. Germline sequencing
cannot supply the tumour genotypes the tumour-directed axis would need, and the pharmacogenomic
layer and PharmCAT share an upstream source, so their agreement is one witness, not two.
