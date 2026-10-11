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
- **Methods covers every results section.** `assets/report/methods.j2` is the single source for both
  templates, and it carries one paragraph per branch of the pipeline, including the variant
  mechanism: where the domain coordinates come from, how NMD sensitivity is computed, and how the
  readthrough characterisation derives the stop codon, its 3' context and the scored products. Every
  version and threshold in it is interpolated from the provenance record or the step's own output,
  so a methods claim cannot drift from what ran. When a results section is added, its methods
  paragraph is added with it. One trap: the templates render with `trim_blocks`, so a line ending
  in a block tag swallows the newline after it and merges the next paragraph into itself. Lines
  that end a paragraph with `{% endif %}` close with `{{ '' }}` to keep the break.
- **Tables carry values; quoted text is keyed.** The stage 1b verdict tables hold the genotype and
  the join result. A guideline paragraph longer than `NOTE_INLINE_MAX` moves out of the cell to a
  note quoted once below the table and keyed a, b, c; a short basis ("no PGx rows in Open Targets")
  stays in the cell. Three defects drove this: the embedded paragraphs made cells four lines deep,
  the columns overran the text width and ran together, and the same paragraph was printed once per
  allele tested.
  - `collapse_rows` merges a drug's rows that differ only in the allele tested, the evidence level
    or the annotation category, listing the alleles and levels in one cell. It is display only: the
    verdicts TSV keeps every raw row, and the test pins that no allele, level or matched allele is
    lost.
  - A note records any *other* drug that PharmCAT's own text names. PharmCAT holds DPWG's
    atorvastatin wording under simvastatin too, and printing that unremarked reads as a join error,
    so the note says the text is quoted as PharmCAT holds it.
- **A LaTeX table cannot overrun the margin.** Column widths are fractions of `\linewidth` through
  the `L{}` column type and sum to 1, instead of fixed centimetres that summed past the text width.
  `L{}` uses `ragged2e`'s `\RaggedRight`, which keeps hyphenation active, so a long single word
  (drug name, gene symbol) breaks instead of standing past the column edge. The `tex` filter inserts
  break points after `_ & / . : ;` and inside any run of 16 or more characters, because an md5, an
  accession or a URL otherwise has none. `placeins` puts a float barrier at each section, so a
  figure can no longer land inside Limitations and split a sentence across pages.
- **One section hierarchy.** The LaTeX report numbers its own sections (Results is 3, the mechanism
  subsections 3.1.1 and 3.1.2, the clinical contexts 3.4.1 onward). Cross-references go through
  `\label`/`\ref`, never hand-written numbers, which previously disagreed with the headings they
  pointed at. The Markdown report carries the same hierarchy: `report.md.j2` sets a `sec` map of
  section numbers once and both its headings and its cross-references read from it, so the two
  formats number the same section the same way and a reference cannot drift from its heading.
- **The feature-absence rule is named once.** `methods.j2` states it and names it; the results and
  limitations cite it by name. It was previously restated nearly verbatim three times, which invites
  skimming over a rule that is load-bearing.
- **Citations follow APA 7.** Literature lives in `assets/report/literature.yml` as structured
  entries (authors, year, sentence-case title, journal, volume, issue, pages, DOI or URL),
  generated from NCBI E-utilities and Crossref records rather than transcribed. Any rendered text
  cites with `[@key]` for a parenthetical citation, `[@+key]` for a narrative one and `[@-key]`
  for the bare form inside an existing parenthesis; structured literature entries (clinical
  contexts, anti-target rules) are matched by PMID. The References section lists the entries
  actually cited, alphabetically, in APA 7 form with the DOI as a link, followed by the tool and
  database list. An unknown key renders as `[?key]` and is reported on the renderer's stdout.
- **Figures are regenerated.** The genome-wide LRR/BAF plot, the per-chromosome plot, the per-arm
  burden plot and the domain map are copied from this run's outputs into `figures/`. Nine more are
  drawn by `REPORT_FIGURES` (`bin/report_figures.py`) from tables the report also prints, so a
  figure never carries a number the text does not. No figure is stored in the repo.
  - In the body: three-caller concordance and the CYP2D6 locus depth track (section 4, the latter
    from `PGX_LOCUS_DEPTH`, all reads and MAPQ >= 20 over CYP2D6/CYP2D7); the stage 1b evidence-level
    heatmap (section 3, regimen block); the ancestry PCA against 1000 Genomes and the
    runs-of-homozygosity plot (section 7).
  - In the appendix: sequencing QC (depth distribution, per-chromosome depth as the sex check,
    insert size), depth at the PharmCAT positions by genotype class, MT-RNR1 allele fractions
    against the heteroplasmy thresholds, and exon copy number across the mechanistic panel.
  - Figures are numbered in order of appearance (`fignum` in the renderer, `\ref` in LaTeX) and
    every one is referred to from the prose that it supports, so the figure is part of the
    argument rather than an illustration beside it. A figure that was not drawn takes no number
    and its reference sentence is omitted with it.
  - Every figure is optional. `<sample>.figures.json` records which were drawn and why any was
    not, and the appendix prints the gaps. Colours follow one palette (categorical slots in fixed
    order, one-hue sequential blue, reserved status colours always paired with a text label).
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
  - Axes with `role: regimen` are rendered in a separate block, "Treatment the child is already
    receiving", never as repurposing candidates, and the summary bullet counts them apart. The
    rhabdomyosarcoma context's "no PharmCAT annotation" line then points to that block.
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
