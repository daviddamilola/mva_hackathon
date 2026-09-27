# 06 — Genome QC: coverage, sex, ancestry, contamination, ROH

Status: **implemented and run** (2026-09-27). The numbers below are aggregates from
`results/qc/genome_qc.json`; no genotypes are recorded here. The per-tool outputs and the MultiQC
report are under `results/qc/` and `results/multiqc/` (gitignored).

## Methods (submission text)

All tools read the CRAM or the VCF with the reference passed explicitly (`subworkflows/local/genome_qc`).

| Check | Tool | Settings |
|---|---|---|
| Mapping, duplicates | samtools 1.22.1 stats/flagstat; markdup stats from alignment | — |
| Coverage | mosdepth 0.3.14 (100 kb bins, shared with the aneuploidy branch) | MAPQ ≥ 20 |
| Sex | somalier 0.3.5 `extract` + `relate` (X/Y depth, X heterozygosity) | cross-checked with chrX non-PAR LRR from the aneuploidy model |
| Ancestry | somalier `ancestry` against 1000 Genomes labelled samples | superpopulation probabilities |
| Contamination | VerifyBamID2 2.0.3, resource `1000g.phase3.100k.b38` (bed/mu renamed to no-chr contigs in-task) | FREEMIX |
| Runs of homozygosity | bcftools 1.22 roh `-G30 --AF-dflt 0.4`, autosomes | cross-check: heterozygous-site density in 1 Mb windows |
| Variant stats | bcftools 1.22 stats on the distributed VCF | — |
| Aggregation | MultiQC 1.35 | — |

`bin/qc_summary.py` reads every metric from a tool output (a missing input is reported as
missing, never filled in) and assigns each check a status against a stated rule.

## Results

| Check | Value | Rule | Status |
|---|---|---|---|
| Sex matches declared | M (X/autosome depth 0.49, 0 X het sites, Y depth 19.3×, chrX non-PAR LRR −1.01) | X/A < 0.7 and Y > 5× ⇒ M | pass |
| Contamination (FREEMIX) | 1.3 × 10⁻⁶ | < 0.02 | pass |
| Mean autosomal depth | 42.6× (X 22.1×) | ≥ 30× | pass |
| Primary mapped | 99.61% (properly paired 98.12%) | ≥ 95% | pass |
| Duplicates | 12.91% | ≤ 20% | pass |
| Ancestry | EUR, p = 0.956 (AMR 0.044) | top probability ≥ 0.9 | pass |
| F_ROH (≥ 1.5 Mb) | 0 (one ROH ≥ 1 Mb, 1.08 Mb; 25 ≥ 0.5 Mb totalling 17 Mb) | < 0.0156 (second-cousin expectation); informational | pass |

1.076 × 10⁹ reads in total. The het-density cross-check agrees with bcftools roh: median 801 het
sites/Mb, and the longest low-het run is 1 Mb.

**Interpretation.**
- The sample is consistent with a single, uncontaminated XY individual at the expected depth.
- There is no evidence of parental relatedness. This is relevant because the recessive MVA genes
  (BUB1B, CEP57, TRIP13) are often found in consanguineous families. The modifier-panel results
  (docs/04) should be read with an outbred background in mind.

## Bug found and fixed during this step

The first QC summary reported ancestry AFR, p = 0.9999. somalier's ancestry table also lists every
labelled 1000 Genomes reference sample, and `qc_summary.py` took the first row, which was a
reference sample. The script now selects rows by sample ID in all three tables (ancestry, relate,
selfSM), and fails if there is not exactly one match.

## Other fixes

- **MultiQC 1.35** names its data directory after the report file (`multiqc_report_data`). The
  module declared `multiqc_data`, so the task failed after MultiQC had finished. `--data-dir` is a
  boolean flag in 1.35, not a path, so the module now renames the directory.
- **QC params:** the QC resource params were first missing from `nextflow.config`, and somalier
  failed with "url null". `main.nf` now checks the required resource params at startup and stops if
  any is undefined.
