# 03 — Pharmacogenomics (stage 1a calling, stage 1b Open Targets annotation)

Status: **implemented and run** (2026-09-26). Stage 1a runs inside `main.nf`. Stage 1b is a
separate entry point, `pgx_annotate.nf`, re-run whenever the candidate drug list changes.

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
undocumented variation instead of assuming reference. PharmCAT's missing-to-reference options
were not used, because no position was left uncalled.

The VCF was renamed to `chr` contigs and preprocessed (`pharmcat_vcf_preprocessor`, with
`-refFna` set to a `chr`-named bgzipped copy of the same masked reference). PharmCAT 3.4.0
(data version 2026-07-13) then ran with CPIC, DPWG and FDA sources.

Aldy 4.8.3 (`genotype -p wgs --genome hg38 -r <reference>`; solver OR-Tools 9.15) called 17
genes directly from the CRAM: CYP2D6, CYP2C19, CYP2C9, CYP2B6, CYP3A5, CYP4F2, TPMT, NUDT15,
DPYD, UGT1A1, SLCO1B1, G6PD, VKORC1, ABCG2, RYR1, CACNA1S and IFNL3. Aldy's CYP2D6 call
(copy-number aware; tandem duplications converted to PharmCAT `xN` notation) was passed to
PharmCAT as an outside call, because PharmCAT does not call CYP2D6 from VCF. An Aldy call is
passed on only when it returns a single major solution.

**Second star-allele caller (Stargazer successor).** Stargazer itself can't be installed
non-interactively: it needs an academic registration at UW and is not on bioconda, PyPI or GitHub.
PyPGx 0.27.0 (bioconda, MIT), the same author's maintained successor that reimplements Stargazer's
genotyping and SV/CNV calling, runs in parallel with Aldy on the same 17 genes. The steps:
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
- **Evidence levels:** only levels 1A and 1B are asserted. Level-3 rows are counted and reported
  as exploratory, and carried no further.

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

## Results (aggregate; details in `results/pgx/`)

| Metric | Value |
|---|---|
| CRAM ↔ reference M5 validation | pass (2,580/2,580) |
| PharmCAT positions genotyped from CRAM | 1,226/1,226; 0 missing genotypes; depth median 42 (range 15–90) |
| Positions with an allele outside PharmCAT's definitions | 1 (flagged by PharmCAT as undocumented variation) |
| Genes in the panel | 19 (17 callable by Aldy + HLA-A, HLA-B) |
| PharmCAT vs Aldy | 14 agree, **0 disagree**, 2 not comparable (nomenclature), 1 outside call (CYP2D6), 2 no call (HLA-A/B) |
| Aldy ambiguous calls | 0 |
| Aldy vs PyPGx (17 genes) | 10 same, 4 same after normalisation, 2 not comparable (nomenclature), **1 different** |
| CYP2D6 structure | Aldy and PyPGx independently call the same diplotype with a whole-gene duplication (PyPGx CNV: `WholeDup1`) |

The single Aldy/PyPGx difference is a **definition-coverage** difference, not a calling
difference. PyPGx's own input VCF carries the heterozygous variant (AD 33/19), but PyPGx 0.27.0's
RYR1 table (49 alleles) does not define it, so PyPGx reports Reference/Reference. PharmCAT 3.4.0
(newer CPIC RYR1 definitions) and Aldy both call it.
| Stage 1b test run (3 example ChEMBL IDs, not the candidate set) | 3/3 drugs had PGx rows (1 via child-molecule resolution); 17/17 called genes had an asserted Open Targets match; 277 extra positions force-called; 0 strand flips; 842 level-3 rows counted as exploratory |

The example drug list reproduces the evidence-level counts in the brief exactly:
fluorouracil 94 level-1A / 604 level-3, and mercaptopurine 53 / 119.

## Limitations

- **HLA-A/B not called.** PharmCAT needs HLA as an outside call, and HLA typing is not
  implemented yet. HLA-dependent recommendations (e.g. abacavir, carbamazepine, allopurinol) are
  reported as no call.
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
