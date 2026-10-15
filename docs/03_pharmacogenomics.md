# 03 — Pharmacogenomics (stage 1a calling, stage 1b Open Targets annotation)

Status: **implemented and run** (stage 1a 2026-09-26; stage 1b on the rhabdomyosarcoma regimen list
2026-10-06). Stage 1a runs inside `main.nf`. Stage 1b is a separate entry point, `pgx_annotate.nf`,
re-run whenever the drug list changes (doc 08).

This document holds methods, design decisions and **aggregate** figures only. Diplotypes,
genotypes and verdicts are subject-derived, so they stay in `results/pgx/` (gitignored).
Everything here is a research hypothesis, not a treatment recommendation.

## Methods (submission text)

**Input and reference validation.** Analyses used the merged, duplicate-marked CRAM from
re-alignment (doc 01) and the reference it was encoded against. Before any CRAM was read, a
validation step compared every `@SQ` line of the CRAM header (name, length, M5) with the
reference sequence dictionary: 2,580/2,580 contigs matched on name, length, MD5 and order. It
also decoded a 10 kb region with the explicit reference. Every CRAM-reading command passes the
reference explicitly (`--fasta-ref`, `-r`, `--reference`).

**Stage 1a — fixed-panel calling.** The 1,226 allele-definition positions shipped with PharmCAT
3.4.0 (bioconda `pharmcat3=3.4.0`) were genotyped directly from the CRAM with bcftools 1.22
(`mpileup --fasta-ref -R <±100 bp windows> -T <positions> -a FORMAT/AD,FORMAT/DP`). Calling was
done twice from the same pileup:
- (i) constrained to PharmCAT's alleles (`call -m -A -C alleles -T <CHROM POS REF,ALT>`), so every
  definition allele, including hom-ref at indel sites, is genotyped explicitly;
- (ii) unconstrained (`call -m`).

Both calls used `--ploidy GRCh38` with the sample's sex (XY; chrX haploid). At any position
where the unconstrained genotype carried a non-reference allele absent from PharmCAT's
definition, the unconstrained record replaced the constrained one. PharmCAT then reports the
undocumented variation instead of assuming reference.

**Minimum evidence, and why "no missing genotypes" proves nothing on its own.** Constrained
calling emits a genotype at every definition position *by construction*, so the absence of
missingness is a property of the method, not of the data. It cannot be used to argue that the
hom-ref calls are safe, and the clinically dangerous direction here is precisely a false reference
call producing a spurious `*1`. Every position is therefore tested against a stated threshold —
depth ≥ `pgx_min_depth` and GQ ≥ `pgx_min_gq`, GQ derived from PL because bcftools does not emit it
here — and a position failing either test is set to a no-call. PharmCAT's missing-to-reference
options are deliberately **not** used, so a no-call reaches PharmCAT as missing and is handled
explicitly rather than silently becoming reference. Per-position evidence is written to
`…pgx_forcecall.evidence.tsv` so the claim can be checked instead of taken on trust.

The VCF was renamed to `chr` contigs and preprocessed (`pharmcat_vcf_preprocessor`, with
`-refFna` set to a `chr`-named bgzipped copy of the same masked reference). PharmCAT 3.4.0
(data version 2026-07-13) then ran with CPIC, DPWG and FDA sources.

Aldy 4.8.3 (`genotype -p wgs --genome hg38 -r <reference>`; solver OR-Tools 9.15) called 19
genes directly from the CRAM: CYP2D6, CYP2C19, CYP2C9, CYP2B6, CYP3A4, CYP3A5, CYP4F2, TPMT,
NUDT15, DPYD, UGT1A1, SLCO1B1, G6PD, VKORC1, ABCG2, RYR1, CACNA1S, IFNL3 and NAT2. CYP3A4 and
NAT2 were added on 2026-10-02 (see "Gene panel" below); MT-RNR1 is excluded because Aldy
addresses the mitochondrion as `MT` while this reference names it `M`. Aldy's CYP2D6 call
(copy-number aware; tandem duplications converted to PharmCAT `xN` notation) was passed to
PharmCAT as an outside call, because PharmCAT does not call CYP2D6 from VCF. An Aldy call is
passed on only when it returns a single major solution.

**Second star-allele caller (Stargazer successor).** Stargazer itself can't be installed
non-interactively: it needs an academic registration at UW and is not on bioconda, PyPI or GitHub.
PyPGx 0.27.0 (bioconda, MIT), the same author's maintained successor that reimplements Stargazer's
genotyping and SV/CNV calling, runs in parallel with Aldy on the same gene set. The steps:
1. The PyPGx gene regions are cut from the CRAM into a BAM with the explicit reference (samtools
   1.22.1). PyPGx's depth step has no reference option, so this avoids implicit `UR:`-based CRAM
   decoding.
2. `create-input-vcf`, then `prepare-depth-of-coverage` (restricted to the panel genes, because
   PyPGx's default GSTT1 region lies on an ALT contig absent from the no-ALT reference), then
   `compute-control-statistics` (control gene VDR).
3. `run-ngs-pipeline --assembly GRCh38 --platform WGS` per gene, with pypgx-bundle 0.27.0
   (sha256-verified).

pandas is pinned to 2.3.3, because pandas 3.0 (new default string dtype) breaks pypgx 0.27.0 at
`api/genotype.py:643` and the recipe does not cap pandas. `bin/compare_star_alleles.py` reports
Aldy vs PyPGx per gene, using the same explicit nomenclature normalisation as the
reconciliation.

**HLA typing.** The reference has no HLA/ALT contigs, so the reads used were:
- the MHC region (6:28,510,120–33,480,577);
- unplaced unmapped reads;

decoded from the CRAM with the explicit reference, mates re-paired (samtools 1.22.1): 902,888
read pairs. Two independent callers typed them:
- **T1K 1.0.10:** `--preset hla-wgs`, with a reference built by `t1k-build.pl` from IPD-IMGT/HLA
  **3.65.0**. The `hla.dat.zip` was pinned by GitHub commit `5b915f2` and verified by sha256.
- **OptiType 1.5.0:** class I, DNA mode, bundled reference, RazerS3 3.5.12 mapping, CBC 2.10.13 ILP
  solver.

Calls were compared at 2-field (protein) resolution with allele order canonicalised
(`bin/hla_reconcile.py`). HLA-A and HLA-B are passed to PharmCAT as outside calls **only when both
callers agree**; a discordant gene is reported and not passed. This rule was checked with a
synthetic discordance. Class II loci are T1K-only and reported as single-caller.

**Reconciliation.** For every gene, the PharmCAT and Aldy diplotypes were compared. Neither
caller was preferred, and each comparison records its basis:
- **Star nomenclature:** Aldy suballeles (e.g. `*4.022`) were reduced to core alleles, and
  PharmCAT `xN` expanded.
- **Non-star nomenclature:**
  - reference synonyms were unified;
  - `c.`-notation was compared with its punctuation removed;
  - PharmCAT's `(heterozygous)` / `(homozygous)` labels were expanded;
  - PharmCAT's own star alias was used where it gives one (e.g. DPYD);
  - a hemizygous chrX call was compared with Aldy's diploid output (Aldy does not model
    hemizygosity).
- Where the two nomenclatures cannot be mapped without assumptions (Aldy VKORC1 H-haplotypes vs
  a PharmCAT rsID call; Aldy base-named IFNL3 alleles), the gene is reported as
  `not_comparable`, with both calls shown.

**Stage 1b — Open Targets annotation.** Open Targets Platform (API 26.9.0, data release 26.09)
was queried through GraphQL, sending only public identifiers:
- **Forward:** every gene called in stage 1a → `target(ensemblId){pharmacogenomics}`.
- **Reverse:** every candidate drug → `drug(chemblId){pharmacogenomics}`. If a ChEMBL ID
  returned no rows, its parent and then its child molecules were queried, and the resolution is
  reported per drug.
- **Panel from the drug list:** the union of positions in the returned rows (`genotypeId` /
  `variantId`, GRCh38) that were not already genotyped was force-called from the CRAM
  (unconstrained `bcftools call -m`, same ploidy). The drug list defines the panel; no gene list
  is hardcoded.

Rows were joined to the child's genotypes as follows:
- **Row type:**
  - If `haplotypeId` is populated, `genotype` is a star allele or diplotype, and the gene is
    taken from the `haplotypeId` prefix.
  - Otherwise, `genotype` is a nucleotide genotype keyed by `variantRsId`, with coordinates
    from `genotypeId`.
- **Star rows:**
  - Diplotype rows match on the canonical (naturally sorted) diplotype.
  - Single-allele rows match any allele the child carries.
  - Suballeles are never collapsed (`*3A` ≠ `*3B`).
- **Nucleotide rows:**
  - Rows match on the allele multiset.
  - A row whose REF is the reverse complement of the genome is flipped, and the flip is counted.
  - A haploid chrX call matches a homozygous row.
- **Missing gene fields:** rows with a null `target` / `targetFromSourceId` are kept, with the
  gene falling back to `UNKNOWN`.
- **Direct-target rows:** rows with `isDirectTarget = true` are excluded from verdicts and
  counted, because this is a metaboliser panel.
- **Evidence levels:** only levels 1A and 1B are asserted. Rows at any other level (2A, 2B, 3, 4)
  are counted per drug, with the genes behind level-2 rows named, and carried no further.
- **Axis roles:** a drug on a `regimen` axis (treatment already given) is queried whatever the
  anti-target rules say; the rule match is reported as a note (doc 08).

**Verdict flags** (per drug × gene) come from PharmCAT's structured guideline fields for the
child's genotype, not from matching keywords in free text:
- `contraindicated`: `alternateDrugAvailable`;
- `dose-adjust`: `dosingInformation`;
- `no flag`: a guideline exists with neither flag set.

Without a guideline:
- `review (annotation only)`: the child is non-reference at an asserted key;
- `no flag`: the child is reference at every asserted key;
- `no data`: there are no asserted rows, or nothing was genotyped.

The join rules are covered by 21 synthetic tests (`tests/unit/test_pgx_join.py`).

## Deviations from the specification, and why

| Specification | Implemented | Reason |
|---|---|---|
| `bcftools mpileup -T positions … \| bcftools call -m` | Also `-R` windows, and two calls (constrained `-C alleles -A` + unconstrained) merged | `-T` alone streams the whole 21.5 GB CRAM (killed after 10 min); `-R` index-jumps (2 s), with `-T` kept as the exact filter. Plain `call -m` writes hom-ref as `ALT=.`: PharmCAT then ignores indel sites and the preprocessor fails to merge `AD` at multiallelic sites. Constrained calling alone hid a real non-reference allele that PharmCAT does not define, hence the merge. |
| Four verdict flags | Adds `review (annotation only)` | Without a guideline, a non-reference genotype at a 1A/1B annotated position does not justify a dosing claim. Labelling it `dose-adjust` would overstate it, and `no flag` would hide it. |
| — | ChEMBL parent/child resolution | e.g. CHEMBL1425 (mercaptopurine anhydrous) has no PGx rows; its child CHEMBL1200751 has 178. The resolution is reported, never silent. |
| — | PharmCAT preprocessor given our reference | Keeps one explicit reference for every step, instead of PharmCAT downloading its own. |

## Gene panel: how it was chosen, and a correction (2026-10-02)

The original 17-gene Aldy set was adopted as a conventional core pharmacogene panel. **No rationale
was recorded for it at the time**, and it does not follow from any single rule: it is neither
everything Aldy supports (40 genes ship with 4.8.3), nor the Aldy/PharmCAT intersection (21), nor
the set with a CPIC guideline (four of the 17 — CYP4F2, IFNL3, VKORC1, and CYP3A4 once added — have
none).

Reviewing it exposed a real defect. PharmCAT was already calling three genes from the force-called
VCF that `params.pgx_genes` then discarded before the report:

| Gene | Call | Phenotype | Why it matters here |
|---|---|---|---|
| **NAT2** | `*5/*5` | **Poor Metabolizer** | acetylates the sulfamethoxazole in co-trimoxazole, standard Pneumocystis prophylaxis during myelosuppression; FDA annotation keyed to NAT2 |
| CYP3A4 | `*1/*1` | Normal Metabolizer | co-metabolises vincristine with CYP3A5, on which the neuropathy finding rests |
| CFTR | Reference/Reference | ivacaftor non-responsive | explains ivacaftor's presence among the original flags; no other relevance here |

The NAT2 omission mattered: the report showed a G6PD row for co-trimoxazole while suppressing a
Poor Metaboliser call bearing on the same drug. All three are now reported. NAT2 and CYP3A4 were
also added to Aldy and PyPGx so that a call driving a clinical statement is not single-caller;
all three callers agree (PharmCAT `*5/*5` = Aldy `*5/*5` = PyPGx `*5/*5` for NAT2) — but see
"What concordance establishes" below for what that is and is not worth.
CFTR is reported from PharmCAT alone and the table marks it `pharmcat_only`.

**Aldy cannot call MT-RNR1 on this reference.** Aldy ships an `mt-rnr1` definition, but it addresses
the mitochondrion as `MT` while this reference names it `M`, so pysam raises `invalid contig MT`.
MT-RNR1 is therefore excluded from `aldy_genes`, and PyPGx provides the cross-check instead
(`Reference/Reference`, agreeing with the dedicated caller's `Reference`). This also corrects an
earlier claim in this project that MT-RNR1 was absent from Aldy's panel: it is present, and the
reason for the bespoke caller is that neither Aldy nor PyPGx models heteroplasmy, which is the
clinically decisive quantity.

**Status-column fix.** A finding whose phenotype is abnormal but which carries no structured
CPIC/DPWG flag previously rendered as "no genotype-based change". It now renders the phenotype
itself (e.g. "Poor Metabolizer; no guideline flag"). This affected NAT2 and the CYP2C19
voriconazole row.

## Reference-allele naming between callers

Star-allele tables do not agree on what the reference allele is called, and the mismatch is
invisible until a sample carries one. PharmVar renamed **CYP2C19\*1 to \*38**: Aldy ships the new
name (its no-mutation allele is `CYP2C19*38.001`), PharmCAT reports the old one. The same genotype
would therefore be scored as a caller disagreement.

This sample escapes it — the call is `*2/*17`, so neither allele is reference and both callers
happen to print identical names. The agreement is real but accidental. Any sample carrying a
reference CYP2C19 allele, which is the commonest allele, would have produced a spurious
disagreement in the most frequently called gene, in a table the report presents as a quality
signal.

Both normalisers (`pgx_reconcile.py` and `compare_star_alleles.py`) now map the two names onto one
another for comparison only; the diplotypes shown in the report stay as each caller reported them.
The table is **gene-scoped and deliberately narrow**:

| Gene | Equated | Why |
|---|---|---|
| CYP2C19 | `*38` ≡ `*1` | PharmVar renaming of the same reference allele |
| NAT2 | **not equated** | Aldy's no-mutation allele is `*1`, PyPGx defaults to `*4`/`*12`; those are distinct haplotypes, not synonyms, and equating them would manufacture agreement |
| IFNL3, DPYD | **not equated** | wild-type naming is already handled by the existing reference and alias paths; IFNL3 remains `not_comparable`, which is the honest result |

`tests/unit/test_star_allele_normalisation.py` pins both directions — that CYP2C19 `*1`/`*38` are
equated, and that the mapping does not leak to other genes or collapse genuinely different
genotypes. Emptying the table fails 5 of the 11 tests; the gene-scoping and NAT2 guards still pass,
since those assert the mapping does *not* over-apply.

## What concordance establishes, and what it does not

Three callers agreeing is **reproducibility, not validation**, and the report should not be read as
if it were the latter.

PharmCAT, Aldy and PyPGx all read the same reads from the same alignment against the same
reference, and their allele definitions share provenance: Aldy ships PharmVar 6.2.14 tables,
PharmCAT's calls here come from ClinPGx, and PyPGx descends from the same PharmVar/Stargazer
lineage. So concordance controls for implementation bugs in any one caller. It cannot detect:

- **reference bias or mismapping**, because all three inherit the same alignment;
- **a shared definition gap**, because the definitions are not independent.

This project's own RYR1 result demonstrates the mechanism rather than refuting it. The single
discordance was definition coverage — PyPGx 0.27.0's RYR1 table does not define the variant that
PharmCAT and Aldy both call. Concordance tracks definition tables, so where the tables agree in
their blind spots, the callers will agree in being wrong.

The regions where this matters most are exactly the ones with paralogues that defeat short-read
mapping: CYP2D6/CYP2D7, NAT1/NAT2, CYP2A6/CYP2A7. All three callers can fail together there.

**What is missing is an anchor to known truth.** Nothing in this project has been run against a
sample with consensus genotypes. The decisive control is a GeT-RM/Coriell reference sample with
CDC consensus PGx genotypes put through this identical pipeline; until that exists, every accuracy
statement here rests on internal consistency. This is recorded as the largest open gap rather than
worked around.

## Ancestry and allele-definition coverage

Star-allele definition tables are markedly better populated for European haplotypes, so the
reassurance carried by a reference call is **ancestry-dependent**, and this report has no machinery
to express that. A `*1/*1` is reported identically whatever the sample's ancestry, when the prior
probability that a real decreased-function allele sits outside the definition table is not equal
across populations.

Alleles that are enriched outside European populations and are variably covered include DPYD
`c.557A>G` (rs115232898), CYP2C9 `*5`, `*6`, `*8` and `*11`, CYP2D6 `*17` and `*29`, and much of the
NAT2 haplotype distribution.

For **this** sample the genome QC branch predicts EUR at p = 0.956 (AMR 0.044), which is where
definition coverage is at its best, so the practical risk here is at the low end. That is a
statement about this sample, not about the pipeline: run on a non-European genome, the same report
would make equally confident reference calls with no qualifier attached. Treat it as a
generalisability limitation, and the ancestry estimate itself is a coarse superpopulation
assignment, not a haplotype-level statement.

## Results (aggregate; details in `results/pgx/`)

| Metric | Value |
|---|---|
| CRAM ↔ reference M5 validation | pass (2,580/2,580) |
| PharmCAT positions genotyped from CRAM | 1,226/1,226; **0 below the evidence threshold** (depth ≥ 10, GQ ≥ 20), so 0 no-calls |
| Depth at definition positions | diploid: min 26, p05 31, median 45. Haploid (chrX non-PAR, 171 positions): min 15 |
| GQ at definition positions (derived from PL) | min 50, p05 96, median 138; none below 20 |
| Heterozygote allele balance | 18 hets, range 0.26–0.76. The extremes are CYP2D6 positions, where three gene copies (`*4/*35x2`) predict 1:2 and 2:1 ratios — the skew is the duplication, not dropout |
| Hom-ref calls carrying any alt read | 42 of 1,203; maximum alt fraction 0.056, an order of magnitude below any heterozygote threshold |
| Positions with an allele outside PharmCAT's definitions | 1 (flagged by PharmCAT as undocumented variation) |
| Genes in the panel | 23 (19 callable by Aldy + HLA-A, HLA-B, MT-RNR1, CFTR) |
| PharmCAT vs Aldy | 16 agree, **0 disagree**, 2 not comparable (nomenclature), 4 outside calls (CYP2D6, HLA-A/B, MT-RNR1), 1 PharmCAT-only (CFTR) |
| Aldy ambiguous calls | 0 |
| HLA class I (T1K vs OptiType, 2-field) | 3/3 genes concordant (HLA-A, -B, -C); HLA-A and HLA-B passed to PharmCAT |
| PGx genes reported | 23/23 (16 agree, 2 not comparable, 4 outside calls: CYP2D6, HLA-A, HLA-B, MT-RNR1, 1 PharmCAT-only: CFTR) |
| Aldy vs PyPGx (20 compared) | 12 same, 4 same after normalisation, 2 not comparable (nomenclature), 1 PyPGx-only (MT-RNR1: Aldy cannot address contig `M`), **1 different** |
| CYP2D6 structure | Aldy and PyPGx independently call the same diplotype with a whole-gene duplication (PyPGx CNV: `WholeDup1`) |

The single Aldy/PyPGx difference is a **definition-coverage** difference, not a calling
difference. PyPGx's own input VCF carries the heterozygous variant (AD 33/19), but PyPGx 0.27.0's
RYR1 table (49 alleles) does not define it, so PyPGx reports Reference/Reference. PharmCAT 3.4.0
(newer CPIC RYR1 definitions) and Aldy both call it.
| Stage 1b, candidate list v1 (`assets/drugs_track2.yml`, Open Targets 26.09, 2026-10-08) | 16/16 names resolved (11 regimen, 5 proposed); 10/16 drugs had any ClinPGx row, **2/16 at an asserted level** (irinotecan, UGT1A1; simvastatin, SLCO1B1; both matched the PharmCAT call, no flag); 20/23 called genes had an asserted forward match (HLA-A, HLA-B, MT-RNR1 did not); 531 extra positions force-called; 0 strand flips; 2,792 rows below the asserted levels counted and not joined; anti-target rule matched 3 regimen agents (etoposide, doxorubicin, irinotecan), reported as notes; 0 proposed candidates excluded, 6 drugs not assessable |
| Stage 1b, regimen list alone (`assets/drugs_regimen.yml`, 2026-10-06) | 10/10 resolved; 1/10 at an asserted level (irinotecan); 503 extra positions force-called; superseded in the report by the list above |

**What the regimen run shows.** The nine agents PharmCAT does not annotate have no ClinPGx
annotation at level 1A or 1B either. So the two layers agree, and for the expected reason: both
draw on ClinPGx, and PharmCAT's drug set is the CPIC, DPWG and FDA guideline set, which these
agents are simply not in. What ClinPGx holds for them sits below the asserted levels:
methotrexate (MTHFR and SLC19A1 at 2A, ATIC at 2B), doxorubicin (SLC28A3 at 2B), and level 3 or 4
rows for the rest; dactinomycin and ifosfamide have no row at all. Those level-2 genes are two of
the report's named panel-gap genes, and their positions are now genotyped in
`results/pgx/open_targets/forcecall/`, but they are not joined, because the assertion threshold
is a design decision (`params.ot_assert_levels`, default `1A,1B`), not a data gap. Lowering it to
include 2A/2B is a one-parameter change; the report names the genes so that choice can be made
with the evidence in view.

An earlier test run on three example ChEMBL IDs (`assets/drugs_example.yml`, 2026-09-26) checked
the join logic and reproduced the evidence-level counts in the brief exactly (fluorouracil 94
level-1A / 604 level-3; mercaptopurine 53 / 119). Its outputs were removed on 2026-10-02 as stale.

## Limitations

- **HLA resolution and scope.** Typing is at 2-field resolution from short reads. Class II is
  single-caller (T1K), and so are rare or novel alleles absent from either reference.
  Clinically, HLA risk-allele testing (e.g. HLA-B\*57:01 before abacavir) is done with a validated
  assay; this research call does not replace it. OptiType's bundled reference is older than
  IPD-IMGT/HLA 3.65.0.
- **Sex for chrX ploidy is a parameter** (`sample_sex = 'M'`), inferred from genotype data in
  doc 00. The QC branch will confirm it from coverage. A wrong value would mis-genotype G6PD.
- **CYP2D6:** PharmCAT's phenotype derives from the Aldy call (outside call). PyPGx gives an
  independent SV-aware cross-check (concordant). Original Stargazer was not run (registration
  required); PyPGx is its successor, not the identical tool.
- **Allele-definition versions differ between callers** (PharmCAT 3.4.0 / Aldy 4.8.3 / PyPGx
  0.27.0). A difference can reflect definition coverage rather than genotype, as for RYR1 here.
- **Aldy licence:** free for academic/non-commercial use only.
- **Blood only.** Germline calls from blood are expected to represent the germline. Mosaic
  aneuploidy could skew allele balance in affected cells (see the aneuploidy analysis).
- **Secondary findings.** The panel includes RYR1 and CACNA1S, which are also on the ACMG
  secondary-findings list. They are reported only as pharmacogenes, using CPIC's
  susceptibility classification. No wider secondary-findings screen was done.
- **Coverage of Open Targets data.** Open Targets PGx data reflects PharmGKB curation as of
  release 26.09. A drug with no PGx rows is reported as `no data`, which is not evidence of no
  effect.
