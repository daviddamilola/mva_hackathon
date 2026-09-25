# mva-track2

Nextflow DSL2 pipeline for **Rare Disease, Real Kid: MVA Hackathon 2026, Track 2 (mechanism and
drug repurposing)**. It extracts four things from one proband's WGS: pharmacogenomics, modifier
genotypes in a mechanistic panel, mosaic aneuploidy burden, and genome QC. This is **not** a
variant-discovery pipeline. Every output is a research hypothesis, not a treatment recommendation.

## Data handling

The data is controlled-access under a DUA. This repo holds analysis code only.

- Subject inputs go in `data/`, outputs in `results/`, and the reference build in `refs/`. All are gitignored,
  along with `*.vcf*`, `*.bam`, `*.cram`, and `*.fastq*`.
- Download inputs with `hf download` from `SageBio/mva-hackathon-2026-data` and check the
  sha256. **Do not read from an `hf-mount-nfs` mount.** It served corrupted bytes during
  inspection (see [docs/00_data_inspection.md](docs/00_data_inspection.md)).
- The pipeline fetches FASTQ one lane at a time and deletes it after alignment.

## Status

| Branch | State |
|---|---|
| Reference rebuild + contig check | implemented |
| Lane-wise FASTQ → CRAM | implemented |
| Pharmacogenomics (PharmCAT + Aldy/Stargazer) | planned |
| Modifier panel (VEP) | planned |
| Mosaic aneuploidy (depth + BAF) | planned |
| Genome QC (ROH, ancestry, sex, contamination, MultiQC) | planned |
| Report (Jinja2 → Markdown + LaTeX) | planned |

## Running

```bash
# 1. DAG check on synthetic inputs (no data, no downloads)
nextflow run . -profile stub -stub-run

# 2. Inputs (VCF + index), verified against upstream sha256
hf download SageBio/mva-hackathon-2026-data \
  WGS_EX2312012_HGWCNDSX7.vcf.gz WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi \
  --repo-type dataset --local-dir data

# 3. Real run (Apple Silicon: force native arm64 conda envs)
CONDA_SUBDIR=osx-arm64 nextflow run . -profile conda -resume
```

Profiles: `conda`, `docker`, `singularity`, `test`/`stub`. Parameters are documented in
[nextflow_schema.json](nextflow_schema.json). You need Hugging Face access to the gated dataset,
through `hf auth login` or `HF_TOKEN`.

## Docs

- [00 — Data inspection](docs/00_data_inspection.md)
- [01 — Reference reconstruction and re-alignment](docs/01_reference_and_alignment.md)

## Acknowledgement

> This work was made possible through the Hackathon, organized by Sage Bionetworks
> in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking,
> Evaluation, and Assessment Consortium for Science), with prize sponsorship from
> AWS and Anthropic. We are deeply grateful to the child and their family who
> generously contributed their data and their story to advance research into this
> rare disease.
