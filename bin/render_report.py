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
import json
import os
import re

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
qc = jload("qc/genome_qc.json")
pgx = jload(f"pgx/{S}.pgx_calls.json")
star = tload(f"pgx/{S}.aldy_vs_pypgx.tsv")
hla = jload(f"pgx/hla/{S}.hla_summary.json")
hla_calls = tload(f"pgx/hla/{S}.hla_calls.tsv")
pc_report = jload(f"pgx/pharmcat/{S}.report.json")
modsum = jload(f"modifier_panel/{S}.modifier_panel.summary.json")
modrf = tload(f"modifier_panel/{S}.modifier_panel.rare_functional.tsv")
aneu = jload(f"aneuploidy/{S}.aneuploidy.summary.json")
arms = jload(f"aneuploidy/{S}.aneuploidy.arms.json")
cnv = jload(f"panel_cnv_sv/{S}.panel.cnv_sv.summary.json")
genecn = tload(f"panel_cnv_sv/{S}.panel.gene_cn.tsv")
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

# ---- aneuploidy arms
arm_rows = [r for r in (arms or []) if r.get("tested")]
events = [r for r in arm_rows if r.get("significant")]
focal = [r for r in arm_rows if (r.get("status") or "").startswith("focal")]

star_counts = {}
for r in star or []:
    star_counts[r["status"]] = star_counts.get(r["status"], 0) + 1

ctx = {
    "sample": S, "prov": prov, "qc": qc, "qc_checks": (qc or {}).get("checks") or [],
    "pgx": pgx, "pgx_rows": pgx_rows, "pgx_summary": (pgx or {}).get("summary"),
    "star_counts": star_counts, "star_rows": [r for r in (star or []) if r["status"] not in ("same",)],
    "hla": hla, "hla_rows": [r for r in (hla_calls or []) if r["gene"] in ("HLA-A", "HLA-B", "HLA-C")],
    "actionable": actionable, "n_guidelines": n_guidelines,
    "pharmcat_data_version": (pc_report or {}).get("dataVersion"),
    "modsum": modsum, "modrf": modrf or [],
    "aneu": aneu, "arm_rows": arm_rows, "events": events, "focal": focal,
    "cnv": cnv, "genecn_calls": [g for g in (genecn or []) if g.get("cn_call") not in (None, "", "none")],
    "ot_sum": ot_sum, "ot_status": ot_status, "ot_cov": ot_cov, "ot_flags": ot_flags, "ot_verdicts": ot_ver or [], "forcecall": forcecall,
    "figs": {
        "genome": fig(f"aneuploidy/{S}.aneuploidy.genome_lrr_baf.png", "Genome-wide depth (LRR) and B-allele frequency"),
        "perchrom": fig(f"aneuploidy/{S}.aneuploidy.per_chromosome.png", "Per-chromosome depth (LRR) and B-allele frequency"),
        "burden": fig(f"aneuploidy/{S}.aneuploidy.arm_burden.png", "Per-arm aneuploidy burden vs detection limit"),
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


def pct(x, nd=1):
    return "not available" if x is None else f"{100 * x:.{nd}f}%"


def tool(name):
    return (prov.get("tools") or {}).get(name, "version not recorded")


def db(name):
    return (prov.get("databases") or {}).get(name, "not recorded")


ctx.update(tool=tool, db=db, params=prov.get("params", {}),
           references=[l.strip() for l in open(os.path.join(a.templates, "references.txt")) if l.strip()])
env = Environment(loader=FileSystemLoader(a.templates), undefined=StrictUndefined, trim_blocks=True, lstrip_blocks=True)
env.filters.update(tex=tex, fmt=fmt, pct=pct)
open(a.out_md, "w").write(env.get_template("report.md.j2").render(**ctx))
open(a.out_tex, "w").write(env.get_template("report.tex.j2").render(**ctx))
print(json.dumps({"missing_inputs": missing, "actionable_recommendations": len(actionable)}))
