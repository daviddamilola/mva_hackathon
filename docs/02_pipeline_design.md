# 02 — Pipeline design

Status (updated 2026-09-27):

| Branch | State | Doc |
|---|---|---|
| Reference rebuild, M5/contig validation | implemented, run | 01 |
| Lane-wise FASTQ → CRAM | implemented, run | 01 |
| PGx stage 1a: PharmCAT + Aldy + PyPGx, reconciled | implemented, run | 03 |
| PGx stage 1b: Open Targets annotation (`pgx_annotate.nf`) | implemented; run on example drugs only | 03 |
| Modifier panel | implemented, run | 04 |
| HLA typing (T1K + OptiType) | implemented, run | 03 |
| Aneuploidy burden | implemented, run | 05 |
| Panel CNV/SV (Delly + exon depth) | implemented, run | 05 |
| Genome QC | implemented, run | 06 |
| Report + provenance | implemented, run | 07 |

The DAG below is the original proposal; changes since then are listed under
"Changes from the proposal".

## DAG

```
                       assets/lanes.csv ──► FETCH_ALIGN_LANE ×4 ──► MERGE_MARKDUP ──► CRAM
NCBI FASTA ─► PREPARE_FASTA ─► FAIDX/DICT ─► CHECK_CONTIGS ◄── VCF header          │
                    └────────► BWA_INDEX ─────────┘                                 │
VCF ──► NORMALISE (bcftools norm -m-any, PASS, reindex) ──┐                         │
                                                          │                         │
 ┌────────────────────────────────────────────────────────┼─────────────────────────┤
 │ 1. PGX                                                  │                         │
 │   PHARMCAT_PREPROCESS (VCF + gap-fill from CRAM at PharmCAT positions) ─► PHARMCAT │
 │   ALDY (CRAM, CYP2D6 + others) ──┐                                                │
 │   HLA typing (T1K/OptiType, CRAM region reads) ─┤                                 │
 │   PGX_RECONCILE ◄─ PharmCAT outside-call file ◄─┘ ─► pgx.json                     │
 │                                                                                   │
 │ 2. MODIFIER PANEL                                                                 │
 │   bcftools view -R panel±10kb ─► VEP (offline cache, --mane_select, gnomAD v4,     │
 │   ClinVar, SpliceAI, AlphaMissense) ─► FILTER (AF<0.01, non-syn/splice) ─► json    │
 │                                                                                   │
 │ 3. ANEUPLOIDY                                                                     │
 │   MOSDEPTH bins (100 kb) ─► GC + mappability (k=150) correction ─► LRR per bin    │
 │   BAF at het SNPs (VCF AD, cross-checked with CRAM pileup) ─────────┐            │
 │   ARM_SEGMENTS (cytoband arms) ◄── LRR ────────────────────────────┴─► ANEUPLOIDY_MODEL
 │   (implemented directly: cell fraction from BAF deviation + LRR, bootstrap CIs)   │
 │                                                                                   │
 │ 4. QC                                                                             │
 │   bcftools roh │ ancestry (projection onto 1KG/HGDP PCs) │ sex (X/Y depth + X het) │
 │   VerifyBamID2 │ samtools stats / mosdepth summary ─► MULTIQC                      │
 └───────────────────────────────────────────────────────────────────────────────────┘
                 all summaries (JSON) + figures ─► PROVENANCE ─► REPORT (Jinja2 → .md + .tex)
```

## Directory structure

```
main.nf  nextflow.config  nextflow_schema.json
conf/        base.config (labels) · modules.config (publishing/args) · test.config (stub)
modules/local/<tool>/main.nf + environment.yml     one process per module, versions.yml each
subworkflows/local/{prepare_reference,align_lanes,pgx,modifier_panel,aneuploidy,genome_qc,report}
workflows/mva.nf
bin/         python helpers (checksums, contig check, arm model, report render)
assets/      lanes.csv · masks · panel BED · cytoband arms · report templates (Jinja2)
tests/stub/  synthetic inputs (no subject data)
docs/        numbered analysis notes (this series)
data/ refs/ results/ work/    gitignored
```

## Changes from the proposal (2026-09-26)

- **CRAM validation is the first node** of every CRAM-consuming branch (`VALIDATE_CRAM_REF`: `@SQ`
  M5 against the reference dict, plus a decode test with the explicit reference).
- **PGx force-calling** uses both constrained and unconstrained calling, merged, so alleles outside
  PharmCAT's definitions surface (doc 03).
- **Second star-allele caller:** PyPGx 0.27.0 (Stargazer's successor; Stargazer itself needs a UW
  registration) runs in parallel with Aldy. `STAR_ALLELE_COMPARE` reports agreement.
- **Stage 1b** is a separate entry point (`pgx_annotate.nf`) driven by a ChEMBL drug list.
- **Modifier panel resources:** gnomAD v4.1 is fetched remotely for panel regions only, not
  downloaded. SpliceAI runs locally on panel variants rather than using precomputed scores or an
  external service. AlphaMissense is downloaded once, verified, cut to the panel and deleted.

## Design decisions (all four proposals accepted 2026-09-25)

1. **PharmCAT absent-site handling.** The VCF is variant-only. Proposal: genotype PharmCAT's
   allele-definition positions directly from the CRAM (`bcftools mpileup -T positions | call`)
   and merge them in, so hom-ref is observed, not assumed.
2. **CYP2D6 caller.** Aldy (reports CN and hybrids), with its call passed to PharmCAT as an
   outside call. *Implemented with PyPGx as the parallel cross-check*, because Stargazer needs
   registration.
3. **Aneuploidy method.** Proposal: implement directly (documented model), not MoChA. MoChA is
   built for array/WGS cohorts and phased BAF; a single sample with no panel-of-normals is
   better served by an explicit model whose assumptions the report can state.
4. **Ancestry reference panel.** Proposal: 1000 Genomes + HGDP (gnomAD v3.1 subset) PCA
   projection on LD-pruned common SNPs.
