# 04 — Modifier genotypes in the mechanistic panel

Status: **implemented and run** (2026-09-26). The numbers below come from
`results/modifier_panel/WGS_EX2312012.modifier_panel.summary.json`.

Subject-derived tables stay in `results/modifier_panel/` (gitignored). This document holds
methods and, once available, aggregate counts only. This is hypothesis generation, not diagnostic
classification.

## Rationale

BubR1 (BUB1B) inhibits APC/C–CDC20. With a weakened BubR1 the brake doesn't hold. Variants that
independently dampen APC/C, or stabilise the mitotic checkpoint complex, could partly compensate.
Variants in those genes therefore change both the rationale and the safe dose for any intervention
on this axis.

## Panel

32 genes, GRCh38, gene body ±10 kb (2,870,084 bp), resolved against Ensembl release 116 by
`bin/make_panel_bed.py`. The 10 Ensembl IDs given in the brief were looked up by ID and each
matched its stated symbol; the rest were looked up by symbol. Files:
`assets/mechanistic_panel_genes.tsv` (input), `assets/mechanistic_panel_grch38.bed` and
`.provenance.tsv`.

| Group | Genes |
|---|---|
| MVA disease genes | BUB1B, CEP57, TRIP13, CENATAC, MAD1L1, SLF2, SMC5 |
| Checkpoint rebalance | MAD2L1BP (p31comet), CDC20, BUB1, BUB3, MAD2L1 |
| APC/C subunits | ANAPC1, ANAPC2, ANAPC4, ANAPC5, ANAPC7, ANAPC10, ANAPC11, ANAPC13, ANAPC15, ANAPC16, CDC16, CDC23, CDC27, FZR1 |
| Error correction / cohesion | PPP2R5A, PPP2R5C, PPP2CA, SGO1, ESPL1, PTTG1 |

## Methods (submission text)

**Variants.** Subject variants within the panel were taken from the distributed single-sample
VCF (Sentieon 202308.02). It was called against the same GRCh38 reference, which was verified
contig-by-contig (doc 01). Multiallelic records were split and indels left-normalised against the
explicit reference (`bcftools norm -m -any -f`, bcftools 1.22). All FILTER values were retained
and reported.

**Annotation.** Ensembl VEP 116.2 ran fully offline (official container
`ensemblorg/ensembl-vep:release_116.2`, arm64, pinned by digest `sha256:5c57abdc…b8da`; the conda
build's perl 5.32.1 segfaults in `DB_File` on the analysis workstation's macOS) with the indexed GRCh38 cache (release 116),
the explicit reference FASTA, `--mane --mane_select --canonical --hgvs --protein --numbers
--check_existing`. External data were joined as exact-allele custom annotations:
- **Population frequency:** gnomAD v4.1 *joint* (exomes + genomes). Fields: `AF_joint`,
  `AF_grpmax_joint`/`grpmax_joint`, `nhomalt_joint`, `fafmax_faf95_max_joint`, per-group
  `AF_joint_*`, and exome/genome filters. Only the panel regions were retrieved, by tabix range
  requests to the public gnomAD bucket, so only public gene coordinates left the machine.
- **ClinVar:** GRCh38 release 2026-09-23, verified against NCBI's MD5. Fields: `CLNSIG`,
  `CLNREVSTAT`, `CLNDN`.
- **AlphaMissense:** hg38 (Cheng et al. 2023; CC BY-NC-SA 4.0). The file is published without an
  index, so it was downloaded once (resumable), verified against Google Cloud Storage's MD5, cut
  to the panel and deleted.
- **SpliceAI:** 1.3.1 run locally on the panel variants only (masked scores `-M 1`, distance
  `-D 50`, explicit reference). The gene annotation was built from the Ensembl 116 MANE Select
  transcripts of the 32 panel genes (`bin/make_spliceai_annotation.py`), the same transcripts as
  the VEP report, instead of SpliceAI's built-in GENCODE v24 table. Subject variants never went to an
  external scoring service. Software: TensorFlow 2.21.0 (conda-forge), numpy 1.26.4. The install
  was validated with a positive control: a synthetic disruption of a BUB1B exon-4 donor gave
  DS_DL = 0.99.

**Transcript and filters.** Per variant, the consequence is reported on the MANE Select
transcript of the panel gene, else the canonical transcript, else the most severe consequence; the
basis is recorded. Variants are flagged, and never removed from the full table:
- **rare:** gnomAD `AF_joint` < 0.01, or absent from gnomAD;
- **functional:** a non-synonymous or splice-region VEP consequence, or SpliceAI max Δ ≥ 0.2.

The rare-and-functional subset lists every such variant the child carries, including variants of
uncertain significance. No ACMG/AMP classification is attempted.

## Results

Aggregate figures only. Per-variant detail is in `results/modifier_panel/` (gitignored).

| Metric | Value |
|---|---|
| Panel | 32 genes, 2,870,084 bp |
| Variants in panel | 3,325 (3,253 PASS) |
| Consequence on MANE Select transcript | 2,328 |
| gnomAD v4.1 joint AF annotated | 3,324 / 3,325 |
| SpliceAI scored | 1,956 (the rest lie outside the panel genes' MANE transcript bodies, where SpliceAI does not score by design) |
| SpliceAI max Δ ≥ 0.2 | 1 (a common variant, AF_joint 0.98) |
| Missense with AlphaMissense score | 13 / 14 |
| Rare (AF_joint < 0.01 or absent) | 94 (1 absent from gnomAD) |
| Functional (non-synonymous / splice) | 34 |
| **Rare and functional** | **5**: MVA disease genes 3, error correction/cohesion 1, APC/C 1 |
| of which ClinVar-annotated | 3 |
| of which failing VCF filters | 1 (the APC/C variant: low MQ/QD, skewed allele balance; low confidence) |
| Rare and functional in checkpoint-rebalance genes (MAD2L1BP, CDC20, BUB1, BUB3, MAD2L1) | 0 |

**Positive control.** Without being targeted, the panel recovered two rare BUB1B variants matching
the publicly described Track 1 configuration: a truncating allele that ClinVar classes P/LP for
MVA1, and an ultra-rare missense predicted damaging by AlphaMissense. Short reads can't phase them
(about 11 kb apart, no parental data), so the biallelic configuration is inferred, not observed.

**Interpretation for Track 2 (hypothesis).** No confident rare, protein-altering variant was found
in the APC/C subunits or the checkpoint-rebalance genes. The only APC/C hit fails the caller's
quality filters, and the remaining non-BUB1B hits are predicted or classified benign. On this
evidence there is no germline modifier expected to dampen APC/C or stabilise the mitotic
checkpoint complex. The residual checkpoint defect should therefore be treated as set by BUB1B
alone when reasoning about interventions on this axis. This is a negative finding within the
limits below, not proof of absence.

## Limitations

- **One of the five rare, functional variants is a low-confidence call** (fails QD/MQ filters,
  3 alt reads). It should be confirmed from the CRAM before it carries any weight.
- **Frequency threshold ≠ benign.** A variant can be "rare" under the 1% cut-off and still be
  classified benign (one here: ClinVar B/LB, 23 gnomAD homozygotes).
- **The variant source is the distributed variant-only VCF.** A panel position without a record
  could be hom-ref or uncovered; coverage over the panel will be reported by the QC branch.
  Structural variants and CNVs in panel genes are not assessed here; they are covered by the
  Delly + exon-depth analysis in doc 05.
- **gnomAD comparison uses the joint AF and grpmax.** Ancestry-matched comparison needs the
  ancestry inference in the QC branch.
- **Germline, blood-derived calls.** Somatic or mosaic variants in affected tissues are not
  captured.
- **In-silico scores are not evidence of effect.** SpliceAI and AlphaMissense are predictions,
  AlphaMissense only scores missense variants, and a modifier effect on the BubR1–APC/C axis would
  need functional testing.
- **Panel scope.** Only the 32 listed genes ±10 kb were examined. No genome-wide or ACMG
  secondary-findings screen was done.
- **Licences.** AlphaMissense scores and SpliceAI model weights are licensed for non-commercial
  use.
