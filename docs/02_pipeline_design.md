# 02 — Pipeline design (proposal for review)

Status: **approved 2026-09-25**. The reference and alignment branches are implemented. The
remaining branches are planned; build order is PGx first, then modifier panel, aneuploidy, QC and
report.

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

## Design decisions (all four proposals accepted 2026-09-25)

1. **PharmCAT absent-site handling.** The VCF is variant-only. Proposal: genotype PharmCAT's
   allele-definition positions directly from the CRAM (`bcftools mpileup -T positions | call`)
   and merge them in, so hom-ref is observed, not assumed.
2. **CYP2D6 caller.** Proposal: Aldy (actively maintained, handles WGS CRAM, reports CN and
   hybrids), with Stargazer as a secondary cross-check only if a license and GDF inputs are
   easy to get. PharmCAT then gets Aldy's call as an outside call.
3. **Aneuploidy method.** Proposal: implement directly (documented model), not MoChA. MoChA is
   built for array/WGS cohorts and phased BAF; a single sample with no panel-of-normals is
   better served by an explicit model whose assumptions the report can state.
4. **Ancestry reference panel.** Proposal: 1000 Genomes + HGDP (gnomAD v3.1 subset) PCA
   projection on LD-pruned common SNPs.
