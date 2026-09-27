# 07 — Report and provenance

Status: **implemented** (2026-09-27). The rendered report contains subject genotypes, so it is
written only to `results/report/` (gitignored). This document describes how it is built, not what
it says.

## Single source, two formats

`bin/render_report.py` builds one context dict from the pipeline's own outputs, then renders it
through two Jinja2 templates in `assets/report/`:

| Template | Role |
|---|---|
| `report.md.j2` | Markdown report |
| `report.tex.j2` | LaTeX report (same context and sections; the `tex` filter escapes text) |
| `methods.j2` | Methods text, shared by both formats; written to be pasted into a submission |
| `limitations.j2` | Limitations, shared by both formats |
| `acknowledgement.j2` | The data acknowledgement, verbatim |
| `references.txt` | One reference per line |

Section order: Summary → Methods → Results (PGx, modifier panel, aneuploidy, QC) → Limitations →
Provenance → References → Acknowledgement.

## Rules the renderer follows

- **Nothing is hardcoded.** Every number, count, version and threshold comes from:
  - a pipeline output (JSON/TSV);
  - the provenance record (tool versions from `versions.yml`);
  - or `params`.
  
  Two early hardcoded values (a site count and the panel size) were found and replaced.
- **Missing is explicit.** Inputs are staged into the REPORT task and resolved by file name. A
  missing input is listed at the top of the report under "Missing inputs", and each dependent
  section prints "not available". Nothing is filled in.
- **Figures are regenerated.** The genome-wide LRR/BAF plot, the per-chromosome plot and the
  per-arm burden plot are copied from this run's outputs into `figures/`. No figure is stored in
  the repo.
- **Framing.** The report opens with the research-hypothesis statement. Methods state that the
  analysis was scoped to the listed panels and that ACMG secondary findings were not screened.
- **Stage 1b (Open Targets verdicts)** is shown only when its summary records the drug list that
  produced it (`drug_list.label`, `drug_list.chembl_ids_queried`).
  - The bundled `assets/drugs_example.yml` sets `drug_list_example: true`. Its output is labelled
    in the report as an example list that checks the join logic, not the Track 2 candidate set.
  - Stage-1b output that does not record its drug list is withheld, and the report says so.
  - This fixed a real problem: the first report displayed the example-list verdicts as if they
    were the candidate drugs.
  - `pgx_annotate.nf` now stops if `drug_list_label` is unset.
  - `assets/NO_OT_SUMMARY` is an empty placeholder, so the input channel is never empty.

## PDF

The REPORT task compiles the `.tex` with `latexmk -pdf` if a LaTeX engine is on `PATH`. A failed
compile does not fail the task: the `.tex` and `.md` remain the deliverables, and `latex.log` is
kept in the work directory.

## Provenance record

`bin/provenance.py` (the PROVENANCE process) writes `results/pipeline_info/provenance.json` with
the following.

**Run metadata** (from the workflow, via `JsonOutput`):
- pipeline name and version;
- git commit and whether the tree had uncommitted changes;
- Nextflow version, session ID, profile, start time;
- all `params`.

**Reference:**
- build, name, source md5.

**Databases:**
- PharmCAT data version (read from the PharmCAT report JSON);
- VEP cache, gnomAD, ClinVar and AlphaMissense versions;
- IPD-IMGT/HLA version;
- the Open Targets release and API version with the query time, or "stage 1b not run".

**Tools:**
- every tool version collected from every task's `versions.yml`;
- a tool seen at more than one version lists all of them, e.g. samtools for the lane-3 exception in
  docs/01. Versions are never collapsed silently;
- `per_process` keeps the per-task record.

The report's Provenance section is rendered from this file, so the report and the record can't
disagree. The run that produced a submitted report should be made from a clean commit, because
the record flags uncommitted changes.

## Outputs

```
results/report/<sample>.track2_report.md
results/report/<sample>.track2_report.tex
results/report/<sample>.track2_report.pdf   (if LaTeX available)
results/report/figures/
results/pipeline_info/provenance.json
```
