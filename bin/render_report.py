#!/usr/bin/env python3
"""Render the Track 2 report (Markdown + LaTeX) from pipeline outputs with Jinja2.

Single source: one context dict is built from the pipeline's JSON/TSV outputs and rendered through
two templates (report.md.j2, report.tex.j2). Nothing is hardcoded. A missing input becomes an
explicit 'not available' entry in the context, and the templates print that.

Usage: render_report.py --results <outdir> --sample S --templates assets/report --provenance provenance.json
                        --out-md report.md --out-tex report.tex
"""
import argparse
import csv
import glob
import gzip
import html
import json
import os
import re

import yaml

from jinja2 import Environment, FileSystemLoader, StrictUndefined

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)
ap.add_argument("--sample", required=True)
ap.add_argument("--templates", required=True)
ap.add_argument("--provenance", required=True)
ap.add_argument("--out-md", required=True)
ap.add_argument("--out-tex", required=True)
ap.add_argument("--figdir", default="figures")
a = ap.parse_args()
R, S = a.results, a.sample
missing = []


def jload(rel):
    p = resolve(rel)
    if p and os.path.getsize(p) > 0:
        return json.load(open(p))
    missing.append(rel)
    return None


def tload(rel):
    p = resolve(rel)
    if p and os.path.getsize(p) > 0:
        return list(csv.DictReader(open(p), delimiter="\t"))
    missing.append(rel)
    return None


def fig(rel, label):
    p = resolve(rel)
    if p:
        os.makedirs(a.figdir, exist_ok=True)
        dst = os.path.join(a.figdir, os.path.basename(p))
        if not os.path.exists(dst):
            import shutil
            shutil.copy(p, dst)
        return {"path": dst, "label": label}
    missing.append(rel)
    return None


# inputs are staged into the task; resolve each by file name, never by an assumed directory layout
INDEX = {}
for root, _, files in os.walk(R, followlinks=True):
    for f in files:
        INDEX.setdefault(f, os.path.join(root, f))


def resolve(rel):
    return INDEX.get(os.path.basename(rel))


prov = json.load(open(a.provenance))
unmet = jload("ot_literature_MONDO_0009759_MONDO_0000141.json")
tract = jload("ot_tractability_BUB1B.json")
# what the drugs near the mitotic machinery do: computed per gene against the Open Targets
# release, so the neighbourhood statement is evidence rather than assertion
nbhd = jload("ot_neighbourhood.json")
# small molecule first: it is the modality a repurposing candidate would most plausibly use
MODALITY_ORDER = ["SM", "AB", "PR", "OC"]
_bm = (tract or {}).get("by_modality") or {}
tract_mod = [(_bm[k]["name"], _bm[k]) for k in MODALITY_ORDER if k in _bm]
tract_mod += [(v["name"], v) for k, v in sorted(_bm.items()) if k not in MODALITY_ORDER]
unmet_rows = sorted((unmet or {}).get("diseases") or [], key=lambda d: d.get("publications") or 0)
qc = jload("qc/genome_qc.json")
pgx = jload(f"pgx/{S}.pgx_calls.json")
star = tload(f"pgx/{S}.aldy_vs_pypgx.tsv")
hla = jload(f"pgx/hla/{S}.hla_summary.json")
hla_calls = tload(f"pgx/hla/{S}.hla_calls.tsv")
pc_report = jload(f"pgx/pharmcat/{S}.report.json")
modsum = jload(f"modifier_panel/{S}.modifier_panel.summary.json")
nmd = jload(f"modifier_panel/{S}.nmd_ptc.json")
domain_map = jload(f"modifier_panel/{S}.domain_map.json")
modrf = tload(f"modifier_panel/{S}.modifier_panel.rare_functional.tsv")
aneu = jload(f"aneuploidy/{S}.aneuploidy.summary.json")
arms = jload(f"aneuploidy/{S}.aneuploidy.arms.json")
cnv = jload(f"panel_cnv_sv/{S}.panel.cnv_sv.summary.json")
genecn = tload(f"panel_cnv_sv/{S}.panel.gene_cn.tsv")
# read-level check of each exon-run CN candidate; without it every candidate stays open
cn_check = jload(f"panel_cnv_sv/{S}.panel.cn_candidates.json")
cn_verdict = {c["gene"]: c for c in (cn_check or {}).get("candidates") or []}
ot_sum = jload(f"pgx/open_targets/{S}.pgx_ot.summary.json")
forcecall = jload(f"pgx/forcecall/{S}.pgx_forcecall.summary.json")
ot_ver = tload(f"pgx/open_targets/{S}.pgx_ot.verdicts.tsv")
# stage 1b is reported only when its summary names the drug list it came from; output from before
# that was recorded (or from an unnamed list) is withheld rather than shown as the Track 2 candidates
ot_status = "not_run"
if ot_sum:
    dl = ot_sum.get("drug_list")
    if not dl or dl.get("label") in (None, "unrecorded"):
        ot_status, ot_sum, ot_ver = "unrecorded", None, None
    else:
        ot_status = "example" if dl.get("example") else "reported"
# stage 1b candidates are presented by axis (rationale), never as one ranked list
dl = (ot_sum or {}).get("drug_list") or {}
ot_axes = dl.get("axes") or {}
ot_counts = dl.get("counts") or {}
ot_cands = dl.get("candidates") or {}
SUPPORT_ORDER = ["released_feature", "syndrome_only", "mechanism_only"]
SUPPORT_HEADING = {
    "released_feature": ("Supported by this child's released clinical features",
                         "These axes are proposed because of features recorded for this child."),
    "syndrome_only": ("Supported by the MVA1 literature, not by the released features",
                      "These axes rest on the published series. The relevant features are not among those "
                      "released to participants, which is not evidence that they are absent. They are "
                      "presented as contingency, not as candidates for this child today."),
    "mechanism_only": ("Derived from this child's causal variants",
                       "These axes are proposed because of the variants found in this child's own genome: each "
                       "candidate under them acts on a consequence of those variants. No clinical feature in the "
                       "released set bears on them either way, so they stand or fall on the mechanism."),
}
# ---- stage 1b verdict rows, as printed
# The join emits one row per (allele tested, evidence level, annotation category), so a drug whose
# guideline text is the same for every allele produced three or four identical paragraphs. The rows
# are collapsed for display only: the TSV keeps every raw row, and the alleles and levels that
# differed between them are merged into one cell rather than dropped.
def collapse_rows(rows):
    out, seen = [], {}
    for r in rows:
        text = r.get("guideline_text") or r.get("flag_basis") or ""
        key = (r.get("gene"), r.get("caller"), r.get("flag"), text, r.get("tier"), r.get("patient_support"))
        if key not in seen:
            seen[key] = dict(r, calls=[], levels=[], matched=[], text=text)
            out.append(seen[key])
        g = seen[key]
        for field, bucket in (("called", "calls"), ("evidence_level", "levels"), ("matched_ot_allele", "matched")):
            v = r.get(field)
            if v and v not in g[bucket]:
                g[bucket].append(v)
    for g in out:
        g["called_merged"] = "; ".join(g["calls"]) or "-"
        g["levels_merged"] = ", ".join(sorted(g["levels"])) or "-"
        g["matched_merged"] = "; ".join(g["matched"]) or "-"
    return out


# every drug name PharmCAT holds, used to spot a guideline paragraph worded for a different drug
PC_DRUG_NAMES = sorted({name for src in ((pc_report or {}).get("drugs") or {}).values() for name in src},
                       key=len, reverse=True)


def note_key(i):
    """a, b, ... z, aa, ab, ... : a key for a guideline paragraph quoted once below its table."""
    k = ""
    i += 1
    while i:
        i, rem = divmod(i - 1, 26)
        k = chr(97 + rem) + k
    return k


# a basis this short says all it has to say in the cell; only a quoted guideline paragraph is keyed
NOTE_INLINE_MAX = 80


def keyed_notes(drugs):
    """Guideline text out of the table cells: one keyed note per unique paragraph, quoted once.

    A note also records any other drug PharmCAT's own text names, because a 'no recommendation'
    paragraph is sometimes worded for a different member of the drug class (PharmCAT holds it that
    way); printing that without comment reads as a join error.
    """
    seen, order = {}, []
    for d in drugs:
        for r in d["rows_collapsed"]:
            t = (r.get("text") or "").strip()
            r["note"], r["basis_short"] = "", ""
            if not t or t == "-":
                continue
            if len(t) <= NOTE_INLINE_MAX:     # "no PGx rows in Open Targets" and the like stay in the cell
                r["basis_short"] = t
                continue
            if t not in seen:
                seen[t] = {"key": note_key(len(order)), "text": t, "drugs": set()}
                order.append(t)
            seen[t]["drugs"].add((d["name"] or "").lower())
            r["note"] = seen[t]["key"]
    notes = []
    for t in order:
        n = seen[t]
        others = [x for x in PC_DRUG_NAMES
                  if x.lower() not in n["drugs"] and re.search(r"\b" + re.escape(x) + r"\b", t, re.I)]
        notes.append({"key": n["key"], "text": t, "for": sorted(n["drugs"]), "names_other_drug": others[:2],
                      # the join concatenates one source's statements with '||'; they are separate statements
                      "statements": [x.strip() for x in t.split("||") if x.strip()]})
    return notes


ot_by_axis = []
def axis_order(ax):
    o = (ot_axes.get(ax) or {}).get("order")
    return (0, o, str(ax)) if o is not None else (1, 0, str(ax))


for ax in sorted(ot_axes, key=axis_order):
    rows = [r for r in (ot_ver or []) if str(r.get("axis")) == str(ax)]
    drugs = []
    for name in dict.fromkeys(r.get("drug_name") for r in rows):
        drows = [r for r in rows if r.get("drug_name") == name]
        c = next((v for v in ot_cands.values() if (v.get("matched_name") or v.get("name")) == name), {})
        ba = ((ot_sum or {}).get("below_asserted_by_drug") or {}).get(drows[0].get("drug")) or {}
        ba_text = ", ".join(f"level {k}: {v}" for k, v in (ba.get("by_level") or {}).items())
        if ba.get("genes_level_2"):
            ba_text += " (level-2 genes: " + ", ".join(ba["genes_level_2"]) + ")"
        drugs.append({"name": name, "tier": drows[0].get("tier"), "rows": drows,
                      "rows_collapsed": collapse_rows(drows), "below_asserted": ba_text,
                      "support": drows[0].get("patient_support") or c.get("patient_support"),
                      "stage": c.get("max_clinical_stage") or "", "investigational": c.get("investigational"),
                      "move_up": c.get("move_up") or "", "caution": c.get("caution") or "",
                      "anti_target": c.get("anti_target"),
                      "anti_target_reason": c.get("anti_target_reason") or ""})
    support = next((d["support"] for d in drugs if d.get("support")), None)
    ot_by_axis.append({"axis": ax, "title": ot_axes[ax].get("title"), "rationale": ot_axes[ax].get("rationale"),
                       "role": ot_axes[ax].get("role") or "repurposing",
                       # an axis whose listing belongs in the appendix: it does not bear on whether a
                       # mechanism-derived candidate is valid (see the note in the drug list)
                       "supplementary": bool(ot_axes[ax].get("supplementary")),
                       "drugs": drugs, "patient_support": support, "notes": keyed_notes(drugs)})
# treatment the child is already receiving is reported in its own block, never as a repurposing axis:
# the anti-target rules do not apply to it, and its drugs were queried whatever the rules said
ot_regimen = [x for x in ot_by_axis if x["role"] == "regimen"]
ot_by_axis = [x for x in ot_by_axis if x["role"] != "regimen"]
ot_n_regimen = ot_counts.get("regimen", 0) if ot_counts else 0
ot_n_repurposing = (ot_counts.get("candidates", 0) - ot_n_regimen) if ot_counts else 0
ot_regimen_anti_target = dl.get("regimen_anti_target") or []
ot_regimen_names = {str(n).lower() for x in ot_regimen for d in x["drugs"] for n in (d["name"] or "").split("; ")}
ot_regimen_names |= {str(v.get("name") or "").lower() for v in ot_cands.values() if v.get("role") == "regimen"}
# axes grouped by how well the candidates in them are supported for THIS child
ot_support_groups = [{"key": k, "heading": SUPPORT_HEADING[k][0], "blurb": SUPPORT_HEADING[k][1],
                      "axes": [x for x in ot_by_axis if x["patient_support"] == k]}
                     for k in SUPPORT_ORDER]
ot_support_groups = [g for g in ot_support_groups if g["axes"]]
# the groups themselves follow the leading axis in each, so a mechanism-directed axis is not pushed
# below a tumour-directed one by its support level alone
ot_support_groups.sort(key=lambda g: min(axis_order(x["axis"]) for x in g["axes"]))
ot_support_ungrouped = [x for x in ot_by_axis if x["patient_support"] not in SUPPORT_ORDER]
# split the groups: the main report carries the axes the mechanism produced, the appendix the rest
def _split_groups(groups, supp):
    out = []
    for g in groups:
        axes = [x for x in g["axes"] if bool(x.get("supplementary")) == supp]
        if axes:
            out.append({**g, "axes": axes})
    return out
ot_groups_main = _split_groups(ot_support_groups, False)
ot_groups_supp = _split_groups(ot_support_groups, True)
ot_regimen_main = [x for x in ot_regimen if not x.get("supplementary")]
ot_regimen_supp = [x for x in ot_regimen if x.get("supplementary")]
# counts for the summary lines the main report keeps in place of each moved listing
ot_supp_stats = {"axes": sum(len(g["axes"]) for g in ot_groups_supp),
                 "drugs": sum(len(a["drugs"]) for g in ot_groups_supp for a in g["axes"]),
                 "regimen_axes": len(ot_regimen_supp),
                 "regimen_drugs": sum(len(a["drugs"]) for a in ot_regimen_supp)}
ot_excluded = dl.get("excluded") or []
ot_unresolved = dl.get("unresolved") or []
ot_rules = dl.get("anti_target_rules") or {}

# the recurrent-miscarriage argument needs the consanguinity check to have come back negative
roh_f = (((qc or {}).get("roh") or {}).get("F_ROH") or {}).get(">=1.5Mb", {}).get("F_ROH")
consanguinity_excluded = roh_f is not None and roh_f == 0

kv = lambda d: ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in d.items()) if d else "none"
ot_cov = kv(ot_sum["join_coverage"]) if ot_sum else ""
ot_flags = kv(ot_sum["flag_counts"]) if ot_sum else ""

# ---- PGx gene table
pgx_rows = []
if pgx:
    for r in pgx["calls"]:
        pgx_rows.append({"gene": r["gene"], "diplotype": r.get("pharmcat_diplotype") or r.get("aldy_diplotype") or "no call",
                         "phenotype": r.get("phenotype") or "n/a", "activity": r.get("activity_score") or "",
                         "source": (r.get("pharmcat_call_source") or "").replace("OUTSIDE", "outside call").replace("MATCHER", "PharmCAT matcher"),
                         "agreement": r.get("caller_agreement"), "undocumented": r.get("undocumented_variation")})

# ---- actionable CPIC/DPWG recommendations (structured PharmCAT flags, no keyword matching)
actionable, n_guidelines = [], 0
if pc_report:
    for src, drugs in (pc_report.get("drugs") or {}).items():
        if not src.startswith(("CPIC", "DPWG")):
            continue
        for name, d in drugs.items():
            for gl in d.get("guidelines") or []:
                for ann in gl.get("annotations") or []:
                    n_guidelines += 1
                    flag = "contraindicated / use alternative" if ann.get("alternateDrugAvailable") else (
                        "dose adjustment" if ann.get("dosingInformation") else None)
                    if flag:
                        genes = sorted({k for lk in (ann.get("lookupKey") or []) for k in (lk.keys() if isinstance(lk, dict) else [])})
                        actionable.append({"drug": name, "source": src.split()[0], "genes": ", ".join(genes), "flag": flag,
                                           "classification": ann.get("classification") or "",
                                           "recommendation": (ann.get("drugRecommendation") or "").strip()})
actionable.sort(key=lambda x: (x["flag"], x["drug"]))


# ---- PGx by clinical context (assets/report/clinical_contexts.yml). Genotype-dependent text is only
# ever PharmCAT's own annotation for the called genotype; the YAML supplies grouping and rationale.
def clean(t):
    t = html.unescape(t or "")
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", t).strip()


# how each context's support for THIS child is labelled; absence of a feature from the released HPO
# set is never rendered as absence of that feature in the child
# phenotype words that mean "not the reference expectation"; used so an abnormal call is never
# summarised as "no genotype-based change"
NON_NORMAL = ("poor metaboliz", "intermediate metaboliz", "rapid metaboliz", "ultrarapid",
              "decreased function", "poor function", "no function", "increased risk",
              "uncertain susceptibility", "non-expresser")

SUPPORT_LABEL = {
    "released_feature": "supported by the released clinical features",
    "syndrome_only": "supported by the MVA1 literature, not by the released features",
    "mechanism_only": "derived from this child's causal variants",
}

SOURCES = [("CPIC Guideline Annotation", "CPIC"), ("DPWG Guideline Annotation", "DPWG"),
           ("FDA Label Annotation", "FDA label"), ("FDA PGx Association", "FDA PGx association")]


def drug_annotations(drug):
    """Every matched PharmCAT annotation for a drug, all four sources; flags from structured fields."""
    out = []
    for key, short in SOURCES:
        d = ((pc_report or {}).get("drugs") or {}).get(key, {}).get(drug)
        for gl in (d or {}).get("guidelines") or []:
            for ann in gl.get("annotations") or []:
                rec = clean(ann.get("drugRecommendation"))
                if not rec:
                    continue
                genes = sorted({k for lk in (ann.get("lookupKey") or []) if isinstance(lk, dict) for k in lk})
                flag = ("alternative drug" if ann.get("alternateDrugAvailable") else
                        "dose adjustment" if ann.get("dosingInformation") else None)
                out.append({"drug": drug, "source": short, "genes": ", ".join(genes), "flag": flag,
                            "classification": ann.get("classification") or "", "recommendation": rec})
    return out


def drug_known(drug):
    return any(drug in (((pc_report or {}).get("drugs") or {}).get(k) or {}) for k, _ in SOURCES)


def gene_view(g):
    x = ((pc_report or {}).get("genes") or {}).get(g) or {}
    d = (x.get("recommendationDiplotypes") or x.get("sourceDiplotypes") or [{}])[0]
    alleles = [a for a in (d.get("allele1"), d.get("allele2")) if a]
    return {"gene": g, "diplotype": d.get("label") or "no call",
            "phenotype": "; ".join(d.get("phenotypes") or []) or "n/a",
            "activity": d.get("activityScore"),
            "alleles": "; ".join(f"{a.get('name')}: {a['function']}" for a in alleles if a.get("function")),
            "source": x.get("callSource")}


pgx_contexts = []
ctx_cfg = yaml.safe_load(open(os.path.join(a.templates, "clinical_contexts.yml")))
for c in ctx_cfg["contexts"] if pc_report else []:
    findings = []
    for f in c.get("findings") or []:
        anns = [x for d in f["drugs"] for x in drug_annotations(d)]
        guideline_flag = any(x["flag"] and x["source"] in ("CPIC", "DPWG") for x in anns)
        label_flag = any(x["flag"] and x["source"].startswith("FDA") for x in anns)
        known = [d for d in f["drugs"] if drug_known(d)]
        # an abnormal phenotype must not render as "no change" just because no structured flag is set:
        # that is the empty-flag-column failure this report exists to avoid
        abnormal = [g for g in [gene_view(x) for x in f["genes"]]
                    if any(k in (g["phenotype"] or "").lower() for k in NON_NORMAL)]
        status = ("flagged by guideline" if guideline_flag else
                  "FDA text only; no CPIC/DPWG flag" if label_flag else
                  # a cited literature risk outranks "no change": the guideline is silent, not reassuring
                  "literature risk; no guideline flag" if f.get("literature") else
                  f"{abnormal[0]['phenotype']}; no guideline flag" if anns and abnormal else
                  "no genotype-based change" if anns else
                  "no annotation matches this genotype" if known else
                  "no PharmCAT-annotated drug in this context")
        groups = {}
        for x in anns:  # identical text across drugs (e.g. the volatile anaesthetics) is shown once
            groups.setdefault((x["source"], x["classification"], x["flag"], x["recommendation"], x["genes"]), []).append(x["drug"])
        grouped = [{"drugs": ", ".join(dict.fromkeys(v)), "source": k[0], "classification": k[1], "flag": k[2],
                    "recommendation": k[3], "genes": k[4]} for k, v in groups.items()]
        findings.append({"label": f["label"], "genes": [gene_view(g) for g in f["genes"]], "status": status,
                         "annotations": grouped, "note": " ".join((f.get("note") or "").split()),
                         # curated, cited literature: shown apart so it is never read as a guideline
                         "literature": [{"statement": " ".join(x["statement"].split()),
                                         "citation": x.get("citation"), "pmid": x.get("pmid")}
                                        for x in f.get("literature") or []],
                         "missing_drugs": [d for d in f["drugs"] if not drug_known(d)]})
    gap = c.get("panel_gap") or {}
    pgx_contexts.append({"id": c["id"], "title": c["title"],
                         "patient_support": c["patient_support"],
                         # supplementary contexts are reported in full in the appendix: they bear on
                         # safety and feasibility, not on whether a repurposing candidate is valid
                         "supplementary": bool(c.get("supplementary")),
                         "support_label": SUPPORT_LABEL.get(c["patient_support"], c["patient_support"]),
                         "why": " ".join(c["why"].split()),
                         "findings": findings,
                         "narrative": [" ".join(x.split()) for x in c.get("narrative") or []],
                         "panel_gap_note": " ".join((gap.get("note") or "").split()),
                         "panel_gap_genes": gap.get("genes") or {},
                         "unannotated": [d for d in c.get("unannotated_drugs") or [] if not drug_known(d)],
                         # which of those the stage 1b regimen block has actually looked up, by name
                         "unannotated_in_stage1b": [d for d in c.get("unannotated_drugs") or []
                                                    if not drug_known(d) and d.lower() in ot_regimen_names],
                         "unannotated_annotated": [d for d in c.get("unannotated_drugs") or [] if drug_known(d)],
                         "unannotated_note": " ".join((c.get("unannotated_note") or "").split())})
pgx_contexts_main = [c for c in pgx_contexts if not c["supplementary"]]
pgx_contexts_supp = [c for c in pgx_contexts if c["supplementary"]]
# counts behind the one-line summary the main report keeps when contexts move to the appendix:
# a reported negative is a finding, so the main text has to say how many there were
pgx_ctx_stats = {"contexts": len(pgx_contexts), "supplementary": len(pgx_contexts_supp),
                 "findings": sum(len(c["findings"]) for c in pgx_contexts),
                 "no_change": sum(1 for c in pgx_contexts for f in c["findings"]
                                  if f["status"] == "no genotype-based change"),
                 "flagged": sum(1 for c in pgx_contexts for f in c["findings"]
                                if f["status"] == "flagged by guideline")}
released_ctx = [c for c in pgx_contexts if c["patient_support"] == "released_feature"]
pgx_negatives = [f"{f['label']} ({c['title'].lower()})" for c in released_ctx for f in c["findings"]
                 if f["status"] == "no genotype-based change"]
pgx_flagged = [f"{f['label']}: {f['status']}" for c in released_ctx for f in c["findings"]
               if f["status"] == "flagged by guideline"
               or f["status"].endswith("no guideline flag")
               or f["status"] == "FDA text only; no CPIC/DPWG flag"]
pgx_literature = [{"context": c["title"], "label": f["label"], "literature": f["literature"]}
                  for c in released_ctx for f in c["findings"] if f["literature"]]
context_drugs = {d for c in ctx_cfg["contexts"] for f in c.get("findings") or [] for d in f["drugs"]}
actionable_outside_contexts = [x for x in actionable if x["drug"] not in context_drugs]
counterfactual = jload(f"pgx/counterfactual/{S}.pgx_counterfactual.json")
cf_context_rows = [dict(c, observed=c["observed"] or {}, counterfactual=c["counterfactual"] or {})
                   for c in ((counterfactual or {}).get("changed") or []) if c["drug"] in context_drugs]
for r in cf_context_rows:
    for k in ("observed", "counterfactual"):
        r[k] = dict(r[k], recommendation=clean(r[k].get("recommendation")))
mt_rnr1 = jload(f"pgx/mt_rnr1/{S}.mt_rnr1.json")
# what a readthrough agent would produce at the PTC: stop codon, 3' context, and the scored products
readthrough = jload(f"modifier_panel/{S}.readthrough.json")
rt = ((readthrough or {}).get("ptcs") or [None])[0]

# ---- aneuploidy arms
arm_rows = [r for r in (arms or []) if r.get("tested")]
# a published MVA1 karyotype against this run's own detection limits, so the null result can be read
# as expected rather than as a broken pipeline (Pavone et al. 2022, PMID 35804254)
PUBLISHED_KARYOTYPE = {"citation": "[@+pavone2022]",
                       "aneuploid_metaphases_pct": 12, "distinct_chromosomes": 11,
                       "most_recurrent_event": "-18", "most_recurrent_pct": 2, "typical_event_pct": 1}
_limits = sorted(r["detection_limit_f"]["loss"] for r in arm_rows
                 if (r.get("detection_limit_f") or {}).get("loss") is not None)
karyotype_cmp = None
if _limits:
    best, median = _limits[0], _limits[len(_limits) // 2]
    karyotype_cmp = dict(PUBLISHED_KARYOTYPE, best_limit=best, median_limit=median,
                         recurrent_vs_best=PUBLISHED_KARYOTYPE["most_recurrent_pct"] / 100 / best,
                         typical_vs_median=PUBLISHED_KARYOTYPE["typical_event_pct"] / 100 / median)
events = [r for r in arm_rows if r.get("significant")]
focal = [r for r in arm_rows if (r.get("status") or "").startswith("focal")]

star_counts = {}
for r in star or []:
    star_counts[r["status"]] = star_counts.get(r["status"], 0) + 1

ctx = {
    "sample": S, "prov": prov, "unmet_rows": unmet_rows, "tract": tract, "tract_mod": tract_mod, "nbhd": nbhd, "nbhd_sum": (nbhd or {}).get("summary") or {}, "qc": qc, "qc_checks": (qc or {}).get("checks") or [],
    "pgx": pgx, "pgx_rows": pgx_rows, "pgx_summary": (pgx or {}).get("summary"),
    "star_counts": star_counts, "star_rows": [r for r in (star or []) if r["status"] not in ("same",)],
    "hla": hla, "hla_rows": [r for r in (hla_calls or []) if r["gene"] in ("HLA-A", "HLA-B", "HLA-C")],
    "actionable": actionable, "n_guidelines": n_guidelines,
    "pgx_contexts": pgx_contexts, "pgx_negatives": pgx_negatives, "pgx_flagged": pgx_flagged,
    "pgx_contexts_main": pgx_contexts_main, "pgx_contexts_supp": pgx_contexts_supp,
    "pgx_ctx_stats": pgx_ctx_stats,
    "pgx_literature": pgx_literature, "actionable_outside_contexts": actionable_outside_contexts,
    "counterfactual": counterfactual, "cf_context_rows": cf_context_rows, "mt_rnr1": mt_rnr1,
    "star_cyp2d6": next((r for r in (star or []) if r["gene"] == "CYP2D6"), None),
    "pharmcat_data_version": (pc_report or {}).get("dataVersion"),
    "modsum": modsum, "modrf": modrf or [],
    "roh_f": roh_f, "consanguinity_excluded": consanguinity_excluded,
    "nmd": nmd, "nmd_ptcs": (nmd or {}).get("ptcs") or [], "domain_map": domain_map,
    "truncating": [v for v in (domain_map or {}).get("variants") or [] if v.get("truncating")],
    "nontruncating": [v for v in (domain_map or {}).get("variants") or [] if not v.get("truncating")],
    "aneu": aneu, "arm_rows": arm_rows, "karyotype_cmp": karyotype_cmp, "events": events, "focal": focal,
    "cnv": cnv, "genecn_calls": [dict(g, check=cn_verdict.get(g["gene"])) for g in (genecn or []) if g.get("cn_call") not in (None, "", "none")],
    "cn_check": cn_check,
    "cn_open": [g["gene"] for g in (genecn or []) if g.get("cn_call") not in (None, "", "none")
                and (cn_verdict.get(g["gene"]) or {}).get("verdict") != "not confirmed"],
    "cn_resolved": [g for g, c in cn_verdict.items() if c["verdict"] == "not confirmed"],
    "ot_sum": ot_sum, "ot_status": ot_status, "ot_by_axis": ot_by_axis, "ot_counts": ot_counts,
    "ot_support_groups": ot_support_groups, "ot_support_ungrouped": ot_support_ungrouped,
    "ot_groups_main": ot_groups_main, "ot_groups_supp": ot_groups_supp,
    "ot_regimen_main": ot_regimen_main, "ot_regimen_supp": ot_regimen_supp,
    "ot_supp_stats": ot_supp_stats,
    "readthrough": readthrough, "rt": rt,
    "ot_regimen": ot_regimen, "ot_n_regimen": ot_n_regimen, "ot_n_repurposing": ot_n_repurposing,
    "ot_regimen_anti_target": ot_regimen_anti_target,
    "ot_excluded": ot_excluded, "ot_unresolved": ot_unresolved, "ot_rules": ot_rules, "ot_cov": ot_cov, "ot_flags": ot_flags, "ot_verdicts": ot_ver or [], "forcecall": forcecall,
    "figs": {
        "genome": fig(f"aneuploidy/{S}.aneuploidy.genome_lrr_baf.png", "Genome-wide depth (LRR) and B-allele frequency"),
        "perchrom": fig(f"aneuploidy/{S}.aneuploidy.per_chromosome.png", "Per-chromosome depth (LRR) and B-allele frequency"),
        "burden": fig(f"aneuploidy/{S}.aneuploidy.arm_burden.png", "Per-arm aneuploidy burden vs detection limit"),
        "domains": fig(f"modifier_panel/{S}.domain_map.png",
                       "Protein domain organisation and this child's variants in the causal gene"),
        # drawn by REPORT_FIGURES (bin/report_figures.py) from the tables this report also prints
        "pgx_concordance": fig(f"{S}.pgx_concordance.png",
                               "Pharmacogene calls by caller. The verdict column is the agreement status used in the table below; colour is never the only carrier."),
        "cyp2d6_depth": fig(f"{S}.cyp2d6_depth.png",
                            "CYP2D6 locus read depth relative to the autosomal mean, all reads and MAPQ >= 20. The elevated depth across CYP2D6 is the read-level evidence behind the duplication call."),
        "readthrough_products": fig(f"{S}.readthrough_products.png",
                                    "Predicted effect of every amino-acid substitution at the readthrough residue, with the three near-cognate products a readthrough agent would produce highlighted, and the residue set against its protein."),
        "stage1b_levels": fig(f"{S}.stage1b_levels.png",
                              "Open Targets (ClinPGx) rows per regimen drug by evidence level. Only the outlined columns are asserted; the rest are counted and carried no further."),
        "ancestry_pca": fig(f"{S}.ancestry_pca.png",
                            "Ancestry: this sample against 1000 Genomes superpopulations (somalier principal components)."),
        "roh": fig(f"{S}.roh.png", "Runs of homozygosity by chromosome against the 1.5 Mb and 5 Mb reference thresholds."),
        "qc_overview": fig(f"{S}.qc_overview.png",
                           "Sequencing QC: depth distribution, per-chromosome mean depth (sex check) and insert size."),
        "forcecall_depth": fig(f"{S}.forcecall_depth.png",
                               "Read depth at every PharmCAT allele-definition position, stacked by genotype class."),
        "mt_rnr1": fig(f"{S}.mt_rnr1.png",
                       "MT-RNR1: alternate allele fraction at each CPIC allele against the heteroplasmy thresholds, with the two homoplasmic positive controls."),
        "panel_exon_cn": fig(f"{S}.panel_exon_cn.png",
                             "Exon-level copy number across the mechanistic panel with |z| >= 3 exons flagged."),
    },
    "fig_manifest": jload(f"{S}.figures.json") or {},
    # figure numbers in order of appearance in the templates; a figure that was not drawn takes no number
    "fignum": {},
    "missing": missing,
}
FIGURE_ORDER = ["domains", "readthrough_products", "stage1b_levels", "pgx_concordance", "cyp2d6_depth", "genome", "perchrom", "burden",
                "ancestry_pca", "roh", "qc_overview", "forcecall_depth", "mt_rnr1", "panel_exon_cn"]
ctx["fignum"] = {k: n for n, k in enumerate([k for k in FIGURE_ORDER if ctx["figs"].get(k)], start=1)}


def tex(s):
    if s is None:
        return ""
    s = str(s)
    rep = {"\\": r"\textbackslash{}", "&": r"\&\allowbreak{}", "%": r"\%", "$": r"\$", "#": r"\#",
           # a break is allowed after _ and &, and below after / . : ; and inside any very long run:
           # reference names, accessions, URLs and DOIs are single words with no break opportunity
           # of their own, and without one they stand past the margin
           "_": r"\_\allowbreak{}", "{": r"\{", "}": r"\}",
           "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "<": r"\textless{}", ">": r"\textgreater{}",
           "|": r"\textbar{}"}
    # A long run with no break opportunity (md5, git commit, accession) would otherwise stand past the
    # margin. Marked before escaping, where such a run is by definition free of special characters.
    s = re.sub(r"[A-Za-z0-9]{16,}", lambda m: "\x00".join(re.findall(r".{1,9}", m.group(0))), s)
    s = re.sub(r"[\\&%$#_{}~^<>|]", lambda m: rep[m.group(0)], s)
    s = re.sub(r"(?<=[/.:;])(?=\S)", "\x00", s)
    return s.replace("\x00", r"\allowbreak{}")


def fmt(x, nd=3):
    if x is None:
        return "not available"
    if isinstance(x, float):
        return f"{x:.{nd}g}" if abs(x) < 1e-3 and x != 0 else f"{x:.{nd}f}".rstrip("0").rstrip(".")
    return str(x)


# Aldy truncates the reference-allele name to "*ference" in its own output (bin/compare_star_alleles.py
# normalises the same synonyms for comparison). The TSV keeps the caller's raw string; only the
# printed name is repaired here, and no other allele name is rewritten.
STAR_DISPLAY_SYNONYMS = {"*ference": "Reference", "*ref": "Reference", "*Ref": "Reference"}


def starname(x):
    """A caller's diplotype as printed: repair truncated reference-allele names, change nothing else."""
    if not x:
        return x
    return "/".join(STAR_DISPLAY_SYNONYMS.get(t, STAR_DISPLAY_SYNONYMS.get(t.lower(), t))
                    for t in str(x).split("/"))


def mdcell(x):
    """A Markdown table cell: pipes would end the column, newlines would end the row."""
    if x is None or x == "":
        return "-"
    return str(x).replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ").strip()


def pct(x, nd=1):
    return "not available" if x is None else f"{100 * x:.{nd}f}%"


def tool(name, default="version not recorded"):
    # default: for a tool whose process did not run in this session (e.g. bwa under
    # --skip_alignment), the caller may supply the version recorded elsewhere (the CRAM @PG line)
    return (prov.get("tools") or {}).get(name, default)


def db(name):
    # query timestamps print as dates: the time of day carries no provenance a reader acts on
    s = (prov.get("databases") or {}).get(name, "not recorded")
    return re.sub(r"(\d{4}-\d{2}-\d{2})T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z?", r"\1", str(s))


# ---- literature citations (APA 7). Any rendered text may cite [@key] (parenthetical), [@+key] (narrative)
# or [@-key] (bare, inside an existing parenthesis); structured literature cites by PMID. The References
# section (rendered last) lists the entries used, alphabetically, in APA 7 form. Rendering runs once per
# template and the used set is reset in between.
LIT = {e["key"]: e for e in yaml.safe_load(open(os.path.join(a.templates, "literature.yml"))) or []}
LIT_BY_PMID = {str(e["pmid"]): k for k, e in LIT.items() if e.get("pmid")}
cite_used, cite_unknown = {}, []


def _names(e):
    return [x["family"] for x in e.get("authors") or []]


def _who(e):
    names = _names(e)
    return names[0] if len(names) == 1 else (f"{names[0]} & {names[1]}" if len(names) == 2 else f"{names[0]} et al.")


# APA 7: entries that would cite identically (same lead author form and year) take year suffixes a, b, ...
# in title order; decided over the whole literature file so in-text citations are right before the
# reference list is rendered
_same = {}
for _k, _e in LIT.items():
    _same.setdefault((_who(_e), str(_e.get("year") or "n.d.")), []).append(_k)
YEAR_SUFFIX = {k: (chr(97 + i) if len(ks) > 1 else "")
               for ks in _same.values() for i, k in enumerate(sorted(ks, key=lambda k: LIT[k].get("title", "")))}


def apa_year(e):
    return str(e.get("year") or "n.d.") + YEAR_SUFFIX.get(e["key"], "")


def apa_intext(e, form="paren"):
    who, year = _who(e), apa_year(e)
    return {"paren": f"({who}, {year})", "narrative": f"{who} ({year})", "bare": f"{who}, {year}"}[form]


def cite_key(key, form="paren"):
    if key not in LIT:
        cite_unknown.append(key)
        return f"[?{key}]"
    cite_used[key] = True
    return apa_intext(LIT[key], form)


FORM = {"": "paren", "+": "narrative", "-": "bare"}


def cite(text):
    """Replace [@key] / [@+key] / [@-key], '(... PMID n)' and bare 'PMID n' with APA author-date citations."""
    if text is None:
        return ""
    text = str(text)
    text = re.sub(r"\[@([+-]?)([\w.-]+)\]", lambda m: cite_key(m.group(2), FORM[m.group(1)]), text)
    text = re.sub(r"\(([^()]*?),?\s*PMID[: ]*(\d{6,9})\)",
                  lambda m: cite_key(LIT_BY_PMID[m.group(2)]) if m.group(2) in LIT_BY_PMID else m.group(0), text)
    return re.sub(r"PMID[: ]*(\d{6,9})",
                  lambda m: cite_key(LIT_BY_PMID[m.group(1)]) if m.group(1) in LIT_BY_PMID else m.group(0), text)


def cite_pmid(pmid, citation=None, bare=False):
    """A structured literature entry: APA citation if its PMID is in literature.yml, else its own text."""
    pmid = str(pmid or "").strip()
    if pmid in LIT_BY_PMID:
        return cite_key(LIT_BY_PMID[pmid], "bare" if bare else "paren")
    if citation:
        return citation if bare else f"({citation})"
    return f"PMID {pmid}" if (pmid and bare) else (f"(PMID {pmid})" if pmid else "")


def apa_authors(e):
    au = e.get("authors") or []
    fmt = lambda x: (x["family"] + (", " + x["initials"] if x.get("initials") else ""))
    if not au:
        return ""
    if len(au) == 1:
        return fmt(au[0])
    if len(au) > 20:
        return ", ".join(fmt(x) for x in au[:19]) + ", ... " + fmt(au[-1])
    return ", ".join(fmt(x) for x in au[:-1]) + ", & " + fmt(au[-1])


def apa_ref(e, style):
    it = (lambda t: f"*{t}*") if style == "md" else (lambda t: "\\emph{" + tex(t) + "}")
    esc = (lambda t: t) if style == "md" else tex
    authors = apa_authors(e)
    authors = esc(authors if authors.endswith(".") else authors + ".")   # a group author ends with its own period
    year, title = esc(apa_year(e)), esc(e.get("title") or "")
    if e.get("type") == "web":
        link = e.get("url") or ""
        return f"{authors} ({year}). {it(title)}. {esc(e.get('journal') or '')}. {esc(link)}"
    src = esc(e.get("journal") or "")
    if e.get("volume"):
        src = it(f"{e['journal']}, {e['volume']}") + (f"({esc(str(e['issue']))})" if e.get("issue") else "")
    if e.get("pages"):
        src += f", {esc(str(e['pages']))}"
    link = f" https://doi.org/{e['doi']}" if e.get("doi") else (f" {e['url']}" if e.get("url") else "")
    return f"{authors} ({year}). {title}. {src}.{esc(link)}"


def refs_used():
    keys = sorted(cite_used, key=lambda k: ((_names(LIT[k]) or [""])[0].lower(), str(LIT[k].get("year")), LIT[k].get("title", "")))
    return [{"key": k, "md": apa_ref(LIT[k], "md"), "tex": apa_ref(LIT[k], "tex")} for k in keys]


ctx.update(tool=tool, db=db, params=prov.get("params", {}), cite=cite, cite_pmid=cite_pmid, refs_used=refs_used,
           references=[l.strip() for l in open(os.path.join(a.templates, "references.txt")) if l.strip()])
env = Environment(loader=FileSystemLoader(a.templates), undefined=StrictUndefined, trim_blocks=True, lstrip_blocks=True)
env.filters.update(tex=tex, fmt=fmt, pct=pct, mdcell=mdcell, cite=cite, starname=starname,
                   mdbold=lambda s: re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", str(s)))
# Table numbering. Markdown has no counter of its own and the LaTeX tables are longtables whose
# captions would repeat on every continuation page, so both formats number through one call here.
# The two templates emit the same tables in the same order, so the numbers agree between them; the
# counter is reset before each render.
class _TableCounter:
    def __init__(self):
        self.n = 0

    def __call__(self):
        self.n += 1
        return self.n

    def reset(self):
        self.n = 0


tnum = _TableCounter()
ctx["tnum"] = tnum
for tpl, out in (("report.md.j2", a.out_md), ("report.tex.j2", a.out_tex)):
    cite_used.clear()
    tnum.reset()
    open(out, "w").write(env.get_template(tpl).render(**ctx))
if cite_unknown:
    print(json.dumps({"unknown_citation_keys": sorted(set(cite_unknown))}))
print(json.dumps({"missing_inputs": missing, "actionable_recommendations": len(actionable)}))
