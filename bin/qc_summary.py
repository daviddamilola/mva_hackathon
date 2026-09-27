#!/usr/bin/env python3
"""Genome QC summary: coverage, mapping, duplicates, sex, ancestry, contamination, ROH.

Every metric is read from a tool output. A metric whose input is missing is reported as missing,
never filled in. Each check gets a status against a stated threshold.

Usage: qc_summary.py --flagstat f.txt --markdup markdup.stats.txt --mosdepth bins.mosdepth.summary.txt
                     --somalier-ancestry anc.tsv --somalier-samples relate.samples.tsv --selfsm s.selfSM
                     --roh roh.summary.json --aneuploidy aneuploidy.summary.json --sample S --declared-sex M --out qc.json
"""
import argparse
import csv
import json
import os
import re

ap = argparse.ArgumentParser()
for x in ("flagstat", "markdup", "mosdepth", "somalier-ancestry", "somalier-samples", "selfsm", "roh", "aneuploidy", "out"):
    ap.add_argument(f"--{x}", required=True)
ap.add_argument("--sample", required=True)
ap.add_argument("--declared-sex", required=True)
ap.add_argument("--max-freemix", type=float, default=0.02)
ap.add_argument("--min-mean-depth", type=float, default=30.0)
a = ap.parse_args()


def have(p):
    return p and os.path.exists(p) and os.path.getsize(p) > 0


def sample_row(path, id_col):
    """The row for this sample. somalier's ancestry table also lists every 1000 Genomes reference
    sample, so taking the first row would report a reference sample's ancestry."""
    rows = [x for x in csv.DictReader(open(path), delimiter="\t") if x.get(id_col) == a.sample]
    if len(rows) != 1:
        raise SystemExit(f"{path}: expected one row for {a.sample} in column {id_col}, found {len(rows)}")
    return rows[0]


qc, checks = {}, []

if have(a.flagstat):
    t = open(a.flagstat).read()
    g = lambda pat: (re.search(pat, t) or [None, None])[1]
    qc["reads_total"] = int(g(r"(\d+) \+ \d+ in total") or 0)
    qc["primary_mapped_pct"] = float(g(r"primary mapped \((\d+\.\d+)%") or "nan")
    qc["properly_paired_pct"] = float(g(r"properly paired \((\d+\.\d+)%") or "nan")
if have(a.markdup):
    t = open(a.markdup).read()
    ex = int(re.search(r"EXAMINED: (\d+)", t).group(1)); dp = int(re.search(r"DUPLICATE TOTAL: (\d+)", t).group(1))
    qc["duplicate_pct"] = round(100 * dp / ex, 2)
if have(a.mosdepth):
    rows = {r["chrom"]: r for r in csv.DictReader(open(a.mosdepth), delimiter="\t")}
    tot = rows.get("total_region") or rows.get("total")
    auto = [float(rows[c]["mean"]) for c in [str(i) for i in range(1, 23)] if c in rows]
    qc["mean_depth_bins"] = float(tot["mean"]) if tot else None
    qc["mean_depth_autosomes"] = round(sum(auto) / len(auto), 2) if auto else None
    qc["mean_depth_X"] = float(rows["X"]["mean"]) if "X" in rows else None
    qc["mean_depth_Y"] = float(rows["Y"]["mean"]) if "Y" in rows else None
if have(a.somalier_ancestry):
    r = sample_row(a.somalier_ancestry, "#sample_id")
    probs = {k.replace("_prob", ""): float(v) for k, v in r.items() if k.endswith("_prob")}
    qc["ancestry"] = {"predicted": r["predicted_ancestry"], "probabilities": probs,
                      "reference": "1000 Genomes superpopulations (somalier)"}
if have(a.somalier_samples):
    r = sample_row(a.somalier_samples, "sample_id")
    fx = lambda k: float(r[k]) if r.get(k) not in (None, "", "NA") else None
    qc["somalier"] = {k: fx(k) for k in ("depth_mean", "X_depth_mean", "X_n", "X_het", "Y_depth_mean", "p_middling_ab")}
if have(a.selfsm):
    r = sample_row(a.selfsm, "#SEQ_ID")
    qc["contamination"] = {"freemix": float(r["FREEMIX"]), "snps": int(r["#SNPS"]), "avg_depth": float(r["AVG_DP"])}
if have(a.roh):
    qc["roh"] = json.load(open(a.roh))
if have(a.aneuploidy):
    an = json.load(open(a.aneuploidy))
    qc["chrX_nonpar_lrr"] = an.get("internal_control_chrX_nonPAR_lrr")

# checks
def check(name, value, ok, rule):
    checks.append({"check": name, "value": value, "status": ("missing" if value is None else ("pass" if ok else "flag")), "rule": rule})

s = qc.get("somalier") or {}
x_ratio = (s["X_depth_mean"] / s["depth_mean"]) if s.get("X_depth_mean") and s.get("depth_mean") else None
inferred = None
if x_ratio is not None and s.get("Y_depth_mean") is not None:
    inferred = "M" if (x_ratio < 0.7 and s["Y_depth_mean"] > 5) else ("F" if (x_ratio > 0.8 and s["Y_depth_mean"] < 2) else "ambiguous")
qc["sex_inferred"] = {"sex": inferred, "x_to_autosome_depth": x_ratio, "x_het_sites": s.get("X_het"), "y_depth": s.get("Y_depth_mean"),
                      "chrX_lrr": qc.get("chrX_nonpar_lrr")}
check("sex matches declared (chrX ploidy used for calling)", inferred, inferred == a.declared_sex,
      f"somalier X/autosome depth < 0.7 and Y depth > 5 => M; declared {a.declared_sex}")
fm = (qc.get("contamination") or {}).get("freemix")
check("contamination (VerifyBamID2 FREEMIX)", fm, fm is not None and fm < a.max_freemix, f"< {a.max_freemix}")
md = qc.get("mean_depth_autosomes")
check("mean autosomal depth", md, md is not None and md >= a.min_mean_depth, f">= {a.min_mean_depth}x")
pm = qc.get("primary_mapped_pct")
check("primary mapped %", pm, pm is not None and pm >= 95, ">= 95%")
dup = qc.get("duplicate_pct")
check("duplicate %", dup, dup is not None and dup <= 20, "<= 20%")
anc = (qc.get("ancestry") or {})
top = max((anc.get("probabilities") or {}).values(), default=None)
check("ancestry assignment confidence", top, top is not None and top >= 0.9, "top superpopulation probability >= 0.9")
fr = ((qc.get("roh") or {}).get("F_ROH") or {}).get(">=1.5Mb", {}).get("F_ROH")
check("F_ROH (>= 1.5 Mb)", fr, fr is not None and fr < 0.0156, "< 0.0156 (below second-cousin expectation); informational")

qc["checks"] = checks
json.dump(qc, open(a.out, "w"), indent=2)
print(json.dumps({c["check"]: c["status"] for c in checks}, indent=1))
