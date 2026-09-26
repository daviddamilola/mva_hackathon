# 01 — Reference reconstruction and re-alignment

Status: **running** (relaunched 2026-09-25 with native arm64 envs). Results below are filled
in from pipeline outputs as they finish. Anything not yet measured is marked *pending*.

## Why re-align

The dataset has FASTQ and a variant-only VCF, but no BAM/CRAM (see
[00_data_inspection.md](00_data_inspection.md)). Aligned reads are needed for:

- CYP2D6 star-allele calling with copy number and hybrid detection (Aldy/Stargazer)
- HLA-A/B typing
- GC- and mappability-corrected binned depth for mosaic aneuploidy (LRR)
- Contamination (VerifyBamID2), coverage summary and coverage-based sex check
- Resolving hom-ref vs no-coverage at PharmCAT positions that are absent from the VCF

## Reference

The VCF header names
`GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta`.
That exact file isn't published, so it is rebuilt:

1. **Source:** NCBI `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz`
   (md5 `a056c57649f3c9964c68aead3849bbf8`, checked in the pipeline).
2. **Rename:** strip the `chr` prefix and FASTA description (`chrM` → `M`, `chrEBV` → `EBV`,
   `chrUn_*` → `Un_*`).
3. **Hard-mask** the GRC exclusion regions with `bedtools maskfasta`
   ([assets/grch38_grc_exclusions_v2_inferred.bed](../assets/grch38_grc_exclusions_v2_inferred.bed)).

### How the mask was inferred

NCBI's `GCA_000001405.15_GRCh38_GRC_exclusions.bed` lists 14 regions. A hard-masked region cannot
carry variant calls, so each region was checked against the subject VCF:

| Region (no-chr) | Reason | VCF calls | In mask? |
|---|---|---|---|
| `Un_KI270752v1:0-27745` | contamination | 0 | yes |
| `21:5010000-5166246` | false duplication | 0 | yes |
| `21:5966593-6161371` | false duplication | 0 | yes |
| `21:6427258-6580181` | false duplication | 0 | yes |
| `21:6789085-6934219` | false duplication | 0 | yes |
| `21:7743700-7865746` | false duplication | 0 | yes |
| `11_KI270721v1_random:0-100316` | false duplication | 0 | yes |
| `X:37085008-37099262` | false duplication | 0 | yes |
| `9:70701801-70719794` | false duplication | 0 | yes |
| `13:86252979-86269524` | false duplication | 0 | yes |
| `22_KI270734v1_random:0-165050` | false duplication | 0 | yes |
| `13:111753855-111793441` | false duplication | 4 (PASS) | **no** |
| `16:34339329-34521510` | false duplication | 11 (MQ40/LowQual) | **no** |
| `10_KI270825v1_alt:0-188315` | contamination | n/a | not in no_alt set |

**Caveat:** a region with zero calls is consistent with masking but doesn't prove it. The
inferred mask covers about 1.2 Mb. A wrong call here would only change alignment inside these
repeat regions, and none of them overlaps a PGx or modifier-panel gene.

### Verification

`CHECK_CONTIGS` compares the rebuilt `.fai` with the VCF `##contig` lines (names, lengths, order)
and fails the run on any mismatch.

**Result (2026-09-25): passed.** 2,580 / 2,580 contigs match the VCF header on name, length and
order (`results/qc/reference/contig_check.json`). Spot check with `samtools faidx`:
`21:5010001-5166246` and `X:37085009-37099262` are 100% N (masked), and the deliberately unmasked
`13:111753856-111793441` has no N (39,586 bases of sequence).

## Alignment

| Step | Tool (pinned) | Parameters |
|---|---|---|
| Fetch | huggingface_hub 0.35.3 (`hf download`) | from upstream `SageBio/mva-hackathon-2026-data`; sha256 checked against the LFS oid |
| Index | bwa 0.7.19 | `index -a bwtsw` |
| Align | bwa 0.7.19 | `mem -K 100000000 -Y -R @RG(ID=HGWCNDSX7.<lane>, SM, LB, PL=ILLUMINA, PU=HGWCNDSX7.<lane>.S16)` |
| Mate tags | samtools 1.22.1 | `fixmate -m` |
| Sort | samtools 1.22.1 | `sort` → CRAM (per lane) |
| Merge + duplicates | samtools 1.22.1 | `merge` then `markdup -d 2500` (optical distance for patterned flowcells) → CRAM + CRAI |

Why these choices:

- **bwa mem 0.7.x**, not bwa-mem2 or DRAGMAP: the workstation has 16 GB RAM, and both of those
  need 32–64 GB for a human index. bwa mem is also the algorithm Sentieon's BWA implements, so
  alignments stay comparable to the ones behind the distributed VCF.
- **`-K 100000000`** makes output independent of thread count, so runs are reproducible.
- **Lane by lane, with FASTQ deleted when each task exits:** disk is limited (~140 GB free). Each
  `FETCH_ALIGN_LANE` task downloads one lane (~21 GB), aligns it, writes a lane CRAM and deletes
  the FASTQ (`trap ... EXIT`, so failed attempts clean up too). `maxForks = 1`. `-resume` reuses
  the lane CRAMs and never needs the FASTQ again.
- **Native conda (osx-arm64)** for the local run: Docker would run amd64 biocontainers under
  emulation. The local miniconda is an **x86_64 build** (`platform: osx-64`), so envs had to be
  solved with `CONDA_SUBDIR=osx-arm64`. The first launch built Intel envs; it was stopped after
  about 5 minutes, those envs were deleted, and the run was relaunched natively. `main.nf` now
  warns if this setting is missing. The docker/singularity profiles use pinned biocontainers where one exists. Wave
  builds the fetch+align image from its pinned `environment.yml`.

### Provenance exception: lane 3 helper tools (2026-09-26)

Lane 3 (the first lane processed) ran with the **pinned bwa 0.7.19**, but `samtools fixmate`
and `samtools sort` came from **Homebrew samtools 1.21** instead of the pinned 1.22.1, and the
download used Homebrew `hf` 2.0.0 instead of huggingface_hub 0.35.3. The CRAM `@PG` header
records this.

- **Cause:** the launching shell had conda `base` active. Activating a task env on top of an
  active env replaces base's PATH entry in place. Base sits after `/opt/homebrew/bin` in PATH,
  so Homebrew binaries shadowed the env's. bwa wasn't affected because Homebrew has no bwa.
- **Impact:** alignments come from bwa, which was the pinned version. `fixmate -m` and
  coordinate `sort` behave the same in 1.21 and 1.22.1. The download was verified by sha256
  regardless of client. So the lane 3 CRAM is valid, but its helper-tool versions differ from
  the other lanes.
- **Fix:** the conda profile now sets
  `process.beforeScript = 'unset CONDA_SHLVL CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER'`,
  which was confirmed to run before activation. Lanes 1, 2 and 4 use the pinned env. Lane 3 was
  kept via `-resume`; the beforeScript change does not alter task hashes.
- **Also affected (from each task's `versions.yml`):**
  - `PREPARE_FASTA`: system curl 8.7.1 instead of 8.14.1. bedtools 2.31.1 was pinned
    correctly. The download was md5-verified.
  - `SAMTOOLS_FAIDX_DICT`: samtools 1.21 instead of 1.22.1.
  - `CHECK_CONTIGS`: Homebrew python 3.14.7 instead of 3.12.11.
  - `BWA_INDEX`: pinned bwa 0.7.19, so not affected.

  None of these change the outputs: a byte-identical download, a FASTA index, and a
  header comparison. They are recorded here for exact provenance.

## Results

| Metric | Value |
|---|---|
| Contig check | passed (2,580/2,580, identical order) |
| Lane FASTQ sha256 verification | *pending* (`results/alignment/lanes/*.fastq.sha256.tsv`) |
| Reads, mapping rate, duplicate rate | *pending* (`results/alignment/*.markdup.stats.txt`, samtools stats) |
| Mean / median coverage | *pending* (mosdepth, genome QC branch) |
| Wall time per lane | *pending* |

## Disk footprint (expected)

| Item | Size | Location |
|---|---|---|
| Masked FASTA + fai + dict | ~3.2 GB | `refs/` (hard-linked from `work/`) |
| bwa index | ~5 GB | `refs/bwa/` |
| Lane CRAMs | ~4–6 GB each | `work/` (can be removed with `nextflow clean` once the merged CRAM exists) |
| Merged CRAM | ~15–20 GB | `results/alignment/` |
| Peak transient per lane | ~21 GB FASTQ + sort temp | task dir, deleted on exit |
