# 00 — Data inspection (pre-pipeline)

Date: 2026-09-25
Scope: layout, integrity, reference build and file characteristics of the MVA Hackathon 2026
dataset, gathered before any pipeline code was written. This document holds only metadata and
aggregate counts. It contains **no genotypes and no clinical text**. Subject data stays in `data/`
(gitignored).

## 1. Source and access

| Item | Value |
|---|---|
| Upstream dataset | `SageBio/mva-hackathon-2026-data` (Hugging Face, gated, DUA) |
| Personal copy | `hf://buckets/Davidreal001/mva-hackathon-2026-data-bucket` |
| Local working copy | `data/` (via `hf download SageBio/mva-hackathon-2026-data <files> --repo-type dataset --local-dir data`) |
| Upstream file count / size | 11 files, ~85 GB |

### Integrity finding: do not read from the NFS mount

The bucket was first mounted at `local/` via `hf-mount-nfs`. Two problems showed up:

1. **The mount served corrupted bytes.** `WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi` read through the
   mount had sha256 `2fcbcf3ea8dd93e7346d5c5d79a894ef1478aefd41cb5a406937df37f7874632` and was
   not valid BGZF. The same file fetched with `hf buckets cp` / `hf download` matched the
   upstream LFS oid.
2. **The mount listing was incomplete.** Only lane L001 FASTQs were visible. `hf buckets ls`
   shows all four lanes.

**Decision:** all inputs are downloaded with `hf download` and checked against the upstream
LFS sha256 before use. The pipeline never reads from `local/`.

## 2. File inventory (upstream)

| File | Size | sha256 (LFS oid) | Local status |
|---|---|---|---|
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` | 315 MB | `bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175` | downloaded, verified |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi` | 2.3 MB | `6f8fed62f11c475fc63a8e2b50925ffc7be33b6930225e821cb977951761a2e0` | downloaded, verified |
| `Challenge_Clinical_Phenotype_1.docx` | 17 kB | (not LFS) | downloaded |
| `README.md` | 5 kB | (not LFS) | downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz` | 10.47 GB | `5ec244c0648552f2b23ae5b1a0b1350e8e464eb5cccd40178a109c44c174939f` | not downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L001_R2_001.fastq.gz` | 11.04 GB | `b80f2cc9f29c229a2bfe557988a4b93c74494fc0d9108b3752bc04cf57239cab` | not downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L002_R1_001.fastq.gz` | 10.38 GB | `0f60e090d9873a18a6914d374fa070a20502aa546d05b51ec520bd9e06ae7e71` | not downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L002_R2_001.fastq.gz` | 10.91 GB | `992ee50db8a89709b7adb327daabcacbf58261503f8c439bb18a3c682e43cd91` | not downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L003_R1_001.fastq.gz` | 10.20 GB | `bbeb51349d4411a661c5922e44650b357d968551ef85735ebb39e6c6071e5d46` | not downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L003_R2_001.fastq.gz` | 10.62 GB | `f3f08f04cd17b089ceb8a8598634f3eb28bc299dc1c00cbfa8e24ba2adb9017b` | not downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L004_R1_001.fastq.gz` | 10.31 GB | `ebc3a6853b9afb5b839ba1e4ff7966eb4f716fba5fee5a39c526d1468079f438` | not downloaded |
| `WGS_EX2312012_HGWCNDSX7_S16_L004_R2_001.fastq.gz` | 10.76 GB | `b16eeae1f84c57fdb889349c25d625eb96a9cb265059004081630131e1037d8b` | not downloaded |

**No BAM or CRAM is distributed.** The VCF header names the BAM it was called from
(`WGS_EX2312012_HGWCNDSX7.bam`), but that BAM is not in the dataset. Any analysis that needs
reads has to start from FASTQ.

## 3. Sequencing characteristics (from a FASTQ sample)

- Instrument `A01973` (Illumina NovaSeq 6000 series), run 164, flowcell `HGWCNDSX7`,
  4 lanes (L001–L004), sample index S16, dual 10 bp indexes.
- Paired-end, **149 bp** reads (10,000 R1 reads sampled from L001; all 149 bp).
- Read length sets the mappability track: use a 150-mer track (Umap k=150 or GenMap k=150).

## 4. VCF characteristics

| Property | Value |
|---|---|
| Format | VCFv4.2, bgzipped, tabix-indexed |
| Samples | **1** (`WGS_EX2312012`), proband only (no trio) |
| Reference | **GRCh38**: `GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta` (chr1 length 248,956,422 confirms GRCh38) |
| Contig naming | **No `chr` prefix** (`1`…`22`, `X`, `Y`, `M` [not `MT`], `EBV`, `Un_*`, `*_random`, `Un_*_decoy`); 2,580 header contigs = UCSC-named NCBI analysis set with `chr` stripped; no ALT or HLA contigs |
| Caller | Sentieon 202308.02 Haplotyper (`--emit_mode GVCF --call_conf 10`, BQSR applied, dbSNP138) → GVCFtyper (`--emit_conf 10`, `--annotation FS,QD,AC,AF,AN`) → GATK VariantFiltration (QD<2, MQ<40, FS>60, MQRankSum<−12.5, ReadPosRankSum<−8) |
| Platform | DNAnexus (header paths) |
| Type | **Variant-only VCF, not a gVCF**: no reference-confidence blocks and no `<NON_REF>` |
| INFO | AC, AF, AN, DP, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, BaseQRankSum, ExcessHet, … (no functional annotation) |
| FORMAT | GT, AD, DP, GQ, PL, PGT, PID |
| Index | Valid, but old tabix format without count metadata (`bcftools index -n/-s` fail). The pipeline should re-index. |

### Aggregate counts

| Metric | Value |
|---|---|
| Total records (all FILTER values) | 5,012,204 |
| PASS records | 4,740,790 |
| PASS SNVs / indels | 3,841,470 / 903,892 |
| PASS multiallelic sites | 66,249 |
| PASS genotype classes | 0/1: 2,847,564 · 1/1: 1,826,977 · 1/2: 66,249 |
| Most common non-PASS filter | MQ40 (177,522) |
| Records on non-primary contigs | 116,498 across 1,666 contigs |
| Ts/Tv (PASS) | 2.00 |
| FORMAT/DP at PASS SNVs, chr1–10 (1-in-50 sample) | median 44, p10 32, p90 56 |
| Autosomal het fraction (het / all calls), per chromosome | 0.59–0.72 |
| chrX non-PAR het fraction (PASS SNVs, X:2,781,480–155,701,382) | 0.012 (994 / 80,141) |
| chrY records | 13,840 |

**Inferred karyotypic sex: XY.** This rests on near-zero non-PAR X heterozygosity plus many Y
calls. The pipeline's sex-check step will confirm it from coverage. The high het fraction among Y
calls (0.78) most likely comes from X-transposed/PAR/repeat mismapping and does not mean anything
biologically.

**Depth:** the median DP at variant sites is about 44×. Genome-wide mean coverage can only be
measured from aligned reads, so treat this as provisional.

## 5. Clinical document

`Challenge_Clinical_Phenotype_1.docx` is marked confidential and "do not redistribute". Its
contents are **not** reproduced here. Two points matter for pipeline design:

- It has **no karyotype, cytogenetic or aneuploidy measurement**, so Analysis 3 has nothing
  published to cross-check against. The report must say that.
- It records an oncology history. This makes chemotherapy exposures (and the PGx genes that
  govern them) concrete, not hypothetical.

## 6. Consequences for the four analyses

| Analysis | From VCF alone? | Why |
|---|---|---|
| 1. PGx (PharmCAT) | **Partly** | This is variant-only, so a site that is absent is ambiguous: it could be hom-ref or no coverage. PharmCAT would have to assume absent = reference, and that must be flagged. A gVCF, or genotyping PharmCAT's positions from aligned reads, removes the ambiguity. |
| 1. CYP2D6 (Aldy / Stargazer) | **No** | Both need BAM/CRAM for copy-number, SV and hybrid-allele calling. |
| 1. HLA-A/B | **No** | No HLA contigs; typing needs reads (OptiType, HLA-LA, T1K). |
| 2. Modifier panel (VEP) | **Yes** | Region-restricted annotation of the VCF. |
| 3. Aneuploidy: BAF | **Yes** | FORMAT/AD at ~2.8 M het sites. |
| 3. Aneuploidy: depth/LRR | **No** | GC/mappability-corrected binned depth needs aligned reads. Variant-site DP is biased and too sparse. |
| 4. ROH | **Yes** | `bcftools roh` with population AF. |
| 4. Ancestry | **Yes** | Projection onto a reference panel. |
| 4. Sex check | **Yes** (genotype-based); coverage-based needs reads | |
| 4. Contamination | **Weak** | VerifyBamID2 needs reads. A VCF-only allele-balance estimate is less sensitive. |
| 4. Coverage summary | **No** | Needs reads. |

**Bottom line:** analysis 2, most of 4, and parts of 1 and 3 run from the VCF. CYP2D6, HLA,
depth-based aneuploidy, contamination and coverage need **re-alignment from FASTQ** to the same
reference (GRCh38 no_alt + hs38d1, no `chr`, masked GRC exclusions v2) to stay consistent with
the VCF.

## 7. Decisions (resolved 2026-09-25)

1. **Re-align from FASTQ: yes.** Process lane by lane and delete each lane's FASTQ after it
   is aligned. See [01_reference_and_alignment.md](01_reference_and_alignment.md).
2. **Reference:** rebuild the exact masked no-chr reference from the NCBI analysis set. See doc 01.
3. **Annotation resources** (VEP cache, gnomAD v4, ClinVar, SpliceAI, AlphaMissense): not yet
   downloaded.
4. **Report framework: Jinja2** templates rendered from JSON summaries. This gives Markdown and
   LaTeX from one source, and missing results are reported explicitly.
