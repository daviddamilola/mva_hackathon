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
    "mechanism_only": ("Derived from the genotype",
                       "These axes follow from the variants themselves rather than from any observed feature."),
}
ot_by_axis = []
for ax in sorted(ot_axes, key=str):
    rows = [r for r in (ot_ver or []) if str(r.get("axis")) == str(ax)]
    drugs = []
    for name in dict.fromkeys(r.get("drug_name") for r in rows):
        drows = [r for r in rows if r.get("drug_name") == name]
        c = next((v for v in ot_cands.values() if (v.get("matched_name") or v.get("name")) == name), {})
        drugs.append({"name": name, "tier": drows[0].get("tier"), "rows": drows,
                      "support": drows[0].get("patient_support") or c.get("patient_support"),
                      "stage": c.get("max_clinical_stage") or "", "investigational": c.get("investigational"),
                      "move_up": c.get("move_up") or "", "anti_target": c.get("anti_target")})
    support = next((d["support"] for d in drugs if d.get("support")), None)
    ot_by_axis.append({"axis": ax, "title": ot_axes[ax].get("title"), "rationale": ot_axes[ax].get("rationale"),
                       "drugs": drugs, "patient_support": support})
# axes grouped by how well the candidates in them are supported for THIS child
ot_support_groups = [{"key": k, "heading": SUPPORT_HEADING[k][0], "blurb": SUPPORT_HEADING[k][1],
                      "axes": [x for x in ot_by_axis if x["patient_support"] == k]}
                     for k in SUPPORT_ORDER]
ot_support_groups = [g for g in ot_support_groups if g["axes"]]
ot_support_ungrouped = [x for x in ot_by_axis if x["patient_support"] not in SUPPORT_ORDER]
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
SUPPORT_LABEL = {
    "released_feature": "supported by the released clinical features",
    "syndrome_only": "supported by the MVA1 literature, not by the released features",
    "mechanism_only": "derived from the genotype",
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
        status = ("flagged by guideline" if guideline_flag else
                  "FDA text only; no CPIC/DPWG flag" if label_flag else
                  # a cited literature risk outranks "no change": the guideline is silent, not reassuring
                  "literature risk; no guideline flag" if f.get("literature") else
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
                         "support_label": SUPPORT_LABEL.get(c["patient_support"], c["patient_support"]),
                         "why": " ".join(c["why"].split()),
                         "findings": findings,
                         "narrative": [" ".join(x.split()) for x in c.get("narrative") or []],
                         "panel_gap_note": " ".join((gap.get("note") or "").split()),
                         "panel_gap_genes": gap.get("genes") or {},
                         "unannotated": [d for d in c.get("unannotated_drugs") or [] if not drug_known(d)],
                         "unannotated_annotated": [d for d in c.get("unannotated_drugs") or [] if drug_known(d)],
                         "unannotated_note": " ".join((c.get("unannotated_note") or "").split())})
released_ctx = [c for c in pgx_contexts if c["patient_support"] == "released_feature"]
pgx_negatives = [f"{f['label']} ({c['title'].lower()})" for c in released_ctx for f in c["findings"]
                 if f["status"] == "no genotype-based change"]
pgx_flagged = [f"{f['label']}: {f['status']}" for c in released_ctx for f in c["findings"]
               if f["status"] in ("flagged by guideline", "FDA text only; no CPIC/DPWG flag",
                                  "literature risk; no guideline flag")]
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

# ---- aneuploidy arms
arm_rows = [r for r in (arms or []) if r.get("tested")]
# a published MVA1 karyotype against this run's own detection limits, so the null result can be read
# as expected rather than as a broken pipeline (Pavone et al. 2022, PMID 35804254)
PUBLISHED_KARYOTYPE = {"citation": "Pavone P et al., Neurol Sci 2022;43:6529-38 (PMID 35804254)",
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
    "sample": S, "prov": prov, "unmet_rows": unmet_rows, "qc": qc, "qc_checks": (qc or {}).get("checks") or [],
    "pgx": pgx, "pgx_rows": pgx_rows, "pgx_summary": (pgx or {}).get("summary"),
    "star_counts": star_counts, "star_rows": [r for r in (star or []) if r["status"] not in ("same",)],
    "hla": hla, "hla_rows": [r for r in (hla_calls or []) if r["gene"] in ("HLA-A", "HLA-B", "HLA-C")],
    "actionable": actionable, "n_guidelines": n_guidelines,
    "pgx_contexts": pgx_contexts, "pgx_negatives": pgx_negatives, "pgx_flagged": pgx_flagged,
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
    "ot_excluded": ot_excluded, "ot_unresolved": ot_unresolved, "ot_rules": ot_rules, "ot_cov": ot_cov, "ot_flags": ot_flags, "ot_verdicts": ot_ver or [], "forcecall": forcecall,
    "figs": {
        "genome": fig(f"aneuploidy/{S}.aneuploidy.genome_lrr_baf.png", "Genome-wide depth (LRR) and B-allele frequency"),
        "perchrom": fig(f"aneuploidy/{S}.aneuploidy.per_chromosome.png", "Per-chromosome depth (LRR) and B-allele frequency"),
        "burden": fig(f"aneuploidy/{S}.aneuploidy.arm_burden.png", "Per-arm aneuploidy burden vs detection limit"),
        "domains": fig(f"modifier_panel/{S}.domain_map.png",
                       "Protein domain organisation and this child's variants in the causal gene"),
    },
    "missing": missing,
}


def tex(s):
    if s is None:
        return ""
    s = str(s)
    rep = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
           "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "<": r"\textless{}", ">": r"\textgreater{}", "|": r"\textbar{}"}
    return re.sub(r"[\\&%$#_{}~^<>|]", lambda m: rep[m.group(0)], s)


def fmt(x, nd=3):
    if x is None:
        return "not available"
    if isinstance(x, float):
        return f"{x:.{nd}g}" if abs(x) < 1e-3 and x != 0 else f"{x:.{nd}f}".rstrip("0").rstrip(".")
    return str(x)


def mdcell(x):
    """A Markdown table cell: pipes would end the column, newlines would end the row."""
    if x is None or x == "":
        return "-"
    return str(x).replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ").strip()


def pct(x, nd=1):
    return "not available" if x is None else f"{100 * x:.{nd}f}%"


def tool(name):
    return (prov.get("tools") or {}).get(name, "version not recorded")


def db(name):
    return (prov.get("databases") or {}).get(name, "not recorded")


ctx.update(tool=tool, db=db, params=prov.get("params", {}),
           references=[l.strip() for l in open(os.path.join(a.templates, "references.txt")) if l.strip()])
env = Environment(loader=FileSystemLoader(a.templates), undefined=StrictUndefined, trim_blocks=True, lstrip_blocks=True)
env.filters.update(tex=tex, fmt=fmt, pct=pct, mdcell=mdcell)
open(a.out_md, "w").write(env.get_template("report.md.j2").render(**ctx))
open(a.out_tex, "w").write(env.get_template("report.tex.j2").render(**ctx))
print(json.dumps({"missing_inputs": missing, "actionable_recommendations": len(actionable)}))
