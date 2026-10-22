# mva-track2

Nextflow DSL2 pipeline for **Rare Disease, Real Kid: MVA Hackathon 2026, Track 2 (mechanism and
drug repurposing)**. From one proband's germline WGS it derives the variant mechanism, then proposes
repurposing candidates that follow from it. This is **not** a variant-discovery pipeline. Every output
is a research hypothesis, not a treatment recommendation.

**What it does, in order.** Re-aligns the FASTQ against a checked reference; characterises the causal
variants against the protein's domain map; works out what a readthrough agent would actually produce
at the nonsense allele (stop codon identity, 3′ context, and the near-cognate products scored with
AlphaMissense); calls pharmacogenes three ways and reconciles them; lets a named drug list define a
pharmacogenomic panel and joins it to the child's genotypes; measures mosaic aneuploidy with a stated
detection limit; and renders one report with every figure and number traced to the step that produced
it. A rule written as code excludes any candidate that inhibits a mitotic-fidelity gene, and reports
what it rejected.

**Run it on another sample** (no data, no downloads — checks the whole DAG):

```bash
cd tests/stub && nextflow run ../.. -profile stub -stub-run && cd ../..
```

For a real sample, point `--sample_id`, `--input` and `--reference_dir` at your own data and run the
command under "Running" below. Everything is version-pinned and every bin/ script's digest is part of
its task hash, so a changed script re-runs rather than silently reusing a cached result.

| Doc | Report section it backs |
|---|---|
| [docs/00_data_inspection.md](docs/00_data_inspection.md) | data provenance, the corrupted-mount finding |
| [docs/01_reference_and_alignment.md](docs/01_reference_and_alignment.md) | methods: reference and alignment |
| [docs/02_pipeline_design.md](docs/02_pipeline_design.md) | overall structure; status of every branch |
| [docs/03_pharmacogenomics.md](docs/03_pharmacogenomics.md) | sections 3 and 4: PGx calls, stage 1b join |
| [docs/04_modifier_panel.md](docs/04_modifier_panel.md) | sections 1 and 5: variant mechanism, readthrough characterisation, modifier genotypes |
| [docs/05_aneuploidy_burden.md](docs/05_aneuploidy_burden.md) | section 6: aneuploidy burden and detection limit |
| [docs/06_genome_qc.md](docs/06_genome_qc.md) | section 7: genome QC |
| [docs/07_report_and_provenance.md](docs/07_report_and_provenance.md) | how the report, its figures and its citations are built |
| [docs/08_candidate_drugs_stage1b.md](docs/08_candidate_drugs_stage1b.md) | section 3: the candidate axes and the anti-target rule |
| [docs/09_clinical_features_and_inference.md](docs/09_clinical_features_and_inference.md) | the epistemic rule governing every phenotype claim |
| [docs/10_methods_description_track2.md](docs/10_methods_description_track2.md) | paste-ready answers for the submission's methods-description form |
| [docs/11_variant_phase_check.md](docs/11_variant_phase_check.md) | the trans-configuration check behind the "Variant phase" limitation |
| [docs/12_video_script.md](docs/12_video_script.md) | spoken script for the submission video |

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
| PGx stage 1a: PharmCAT (CRAM force-call) + Aldy, reconciled | implemented |
| PGx stage 1b: Open Targets annotation (`pgx_annotate.nf`, named drug list) | implemented; awaits Track 2 candidate list |
| PGx second star-allele caller: PyPGx (Stargazer successor) + Aldy comparison | implemented |
| HLA typing: T1K + OptiType, concordant HLA-A/B → PharmCAT | implemented, run |
| Modifier panel: VEP 116 + gnomAD v4.1 + ClinVar + AlphaMissense + SpliceAI | implemented, run |
| Mosaic aneuploidy (mosdepth LRR + BAF model, per-arm cell fraction) | implemented, run |
| Panel-gene CNV/SV (Delly + exon depth) | implemented, run |
| Genome QC (ROH, ancestry, sex, contamination, MultiQC) | implemented, run |
| Report (Jinja2 → Markdown + LaTeX + PDF) and provenance record | implemented, run |

## Running

```bash
# 1. DAG check on synthetic inputs (no data, no downloads); run from tests/stub so it has its
#    own session history and can never be picked up by a real run's -resume
cd tests/stub && nextflow run ../.. -profile stub -stub-run && cd ../..

# 2. Inputs (VCF + index), verified against upstream sha256
hf download SageBio/mva-hackathon-2026-data \
  WGS_EX2312012_HGWCNDSX7.vcf.gz WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi \
  --repo-type dataset --local-dir data

# 3. Real run (Apple Silicon: force native arm64 conda envs)
CONDA_SUBDIR=osx-arm64 nextflow run . -profile conda,macos_docker -resume <session-id>
#    (macos_docker: only VEP and Delly run in native arm64 containers, because their osx-arm64 builds crash;
#     needs Docker Desktop, and once:  docker build -t mva-track2/delly:2.6.0 containers/delly)

# 4. PGx stage 1b for a candidate drug list (re-run whenever candidates change)
CONDA_SUBDIR=osx-arm64 nextflow run pgx_annotate.nf -profile conda -params-file my_drugs.yml
#    my_drugs.yml:  drug_list_label: track2_candidates_v1
#                   drugs: [CHEMBL185, ...]   (see assets/drugs_example.yml)

# modifier panel needs the verified VEP cache once (resumable, checksum-verified, ~28 GB):
bin/fetch_vep_cache.sh homo_sapiens GRCh38 116 refs/vep

# unit tests for the Open Targets join rules (synthetic data)
python3 tests/unit/test_pgx_join.py
python3 tests/unit/test_aneuploidy_model.py   # needs numpy/scipy/matplotlib (modules/local/aneuploidy_model/environment.yml)
```

Resume by explicit session ID (`nextflow log` lists them). A bare `-resume` takes the most
recent session in the launch directory, which may not be the one holding your cached work.

Profiles: `conda`, `docker`, `singularity`, `test`/`stub`. Parameters are documented in
[nextflow_schema.json](nextflow_schema.json). You need Hugging Face access to the gated dataset,
through `hf auth login` or `HF_TOKEN`.

## Docs

Each doc records the design decisions and the aggregate findings behind one part of the report; the
table at the top of this file maps them to report sections. They hold no genotypes and no clinical
text. Numbered `docs/00` through `docs/09`.

## Acknowledgement

> This work was made possible through the Hackathon, organized by Sage Bionetworks
> in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking,
> Evaluation, and Assessment Consortium for Science), with prize sponsorship from
> AWS and Anthropic. We are deeply grateful to the child and their family who
> generously contributed their data and their story to advance research into this
> rare disease.
