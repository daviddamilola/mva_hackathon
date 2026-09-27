#!/usr/bin/env python3
"""Reconcile HLA calls from T1K and OptiType and write PharmCAT outside calls.

Class I (HLA-A, -B, -C) is typed by both callers with independent references (T1K: pinned
IPD-IMGT/HLA release; OptiType: its bundled reference). Alleles are compared at 2-field resolution
(protein level, e.g. *57:01), which is what CPIC/PharmCAT use, and allele order is canonicalised.
A gene is passed to PharmCAT as an outside call only if (1) it is in --pharmcat-genes and (2) both
callers agree at 2 fields. A discordant gene is reported and NOT passed; neither caller is
preferred. Class II and other loci typed only by T1K are reported with 'single caller'.

Usage: hla_reconcile.py --t1k WGS_genotype.tsv --optitype WGS_result.tsv --pharmcat-genes HLA-A,HLA-B
                        --min-t1k-quality 0 --prefix out
"""
import argparse
import csv
import json
import re

ap = argparse.ArgumentParser()
ap.add_argument("--t1k", required=True)
ap.add_argument("--optitype", required=True)
ap.add_argument("--pharmcat-genes", default="HLA-A,HLA-B")
ap.add_argument("--min-t1k-quality", type=float, default=0)
ap.add_argument("--prefix", required=True)
a = ap.parse_args()
pc_genes = {g.strip() for g in a.pharmcat_genes.split(",") if g.strip()}


def two_field(allele):
    """'HLA-A*24:02:01' / 'A*24:02' / '*24:02:01:02L' -> '*24:02' (None if unresolvable)."""
    if not allele or allele in (".", "NA"):
        return None
    m = re.search(r"\*(\d+:\d+)", allele)
    return f"*{m.group(1)}" if m else None


def dip(alleles):
    xs = [x for x in alleles if x]
    if not xs:
        return None
    if len(xs) == 1:           # one allele reported = homozygous at typed resolution
        xs = xs * 2
    return "/".join(sorted(xs, key=lambda x: tuple(int(v) for v in x[1:].split(":"))))


t1k = {}
for row in csv.reader(open(a.t1k), delimiter="\t"):
    if not row or not row[0].startswith("HLA-"):
        continue
    gene, n = row[0], int(row[1])
    als, quals = [], []
    for i in range(n):
        al, q = row[2 + 3 * i], float(row[4 + 3 * i])
        quals.append(q)
        if q >= a.min_t1k_quality:
            als.append(two_field(al))
    t1k[gene] = {"diplotype": dip(als), "raw": ";".join(row[2 + 3 * i] for i in range(n)), "min_quality": min(quals) if quals else None}

opti = {}
rows = list(csv.DictReader(open(a.optitype), delimiter="\t"))
if rows:
    r = rows[0]
    for g in ("A", "B", "C"):
        opti[f"HLA-{g}"] = {"diplotype": dip([two_field(r.get(f"{g}1")), two_field(r.get(f"{g}2"))]),
                            "raw": f"{r.get(g + '1')};{r.get(g + '2')}"}
    opti_meta = {"reads": r.get("Reads"), "objective": r.get("Objective")}
else:
    opti_meta = {}

out, outside = [], []
for gene in sorted(set(t1k) | set(opti), key=lambda g: (g not in ("HLA-A", "HLA-B", "HLA-C"), g)):
    d1 = (t1k.get(gene) or {}).get("diplotype")
    d2 = (opti.get(gene) or {}).get("diplotype")
    if d1 and d2:
        status = "agree" if d1 == d2 else "disagree"
    elif d1:
        status = "single caller (T1K)"
    elif d2:
        status = "single caller (OptiType)"
    else:
        status = "no call"
    passed = gene in pc_genes and status == "agree"
    if passed:
        outside.append((gene, d1))
    out.append({"gene": gene, "t1k_2field": d1 or "", "optitype_2field": d2 or "", "agreement": status,
                "passed_to_pharmcat": passed, "t1k_raw": (t1k.get(gene) or {}).get("raw", ""),
                "t1k_min_quality": (t1k.get(gene) or {}).get("min_quality", ""), "optitype_raw": (opti.get(gene) or {}).get("raw", "")})

with open(a.prefix + ".hla_calls.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(out)
with open(a.prefix + ".hla_outside_calls.tsv", "w") as fh:
    for g, d in outside:
        fh.write(f"{g}\t{d}\n")
summary = {"class_I_agreement": {r["gene"]: r["agreement"] for r in out if r["gene"] in ("HLA-A", "HLA-B", "HLA-C")},
           "passed_to_pharmcat": [g for g, _ in outside],
           "not_passed_discordant": [r["gene"] for r in out if r["gene"] in pc_genes and r["agreement"] == "disagree"],
           "optitype": opti_meta, "resolution": "2-field"}
json.dump(summary, open(a.prefix + ".hla_summary.json", "w"), indent=2)
print(json.dumps(summary))
