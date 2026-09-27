#!/usr/bin/env python3
"""Panel-gene copy number (exon read depth) and structural variants (Delly) overlapping the panel.

EXON COPY NUMBER
  Expected depth for an exon = median depth of autosomal genome bins (kept: mappable >= min_map,
  no blacklist) in the same GC stratum (0.5% strata; neighbours pooled until >= 30 bins), so exon
  CN = 2 * observed / expected. Per gene: median exon CN; runs of >= min_run consecutive exons below
  del_cn or above dup_cn are flagged. Exons with mappable fraction < min_map are reported but not
  used for gene calls. With mosaic aneuploidy, an exon CN shift matching the arm-level LRR is
  expected and is labelled as such, not as a focal event.

STRUCTURAL VARIANTS
  Delly PASS calls (all types) overlapping a panel gene +/- flank (from the panel BED) or a MANE exon,
  including breakends (BND) whose either end falls in the panel. Reported with type, size,
  genotype, paired-end/split-read support and the exons hit.

Usage: gene_cn.py --exon-depth exons.regions.bed.gz --exons exons.annot.tsv --bins bins_with_depth.tsv
                  --delly sv.vcf.gz --panel panel.bed --arms arms.json --prefix out
"""
import argparse
import csv
import gzip
import json
import re
from collections import defaultdict

import numpy as np

ap = argparse.ArgumentParser()
for x in ("exon-depth", "exons", "bins", "delly", "panel", "arms", "prefix"):
    ap.add_argument(f"--{x}", required=True)
ap.add_argument("--min-map", type=float, default=0.9)
ap.add_argument("--del-cn", type=float, default=1.5)
ap.add_argument("--dup-cn", type=float, default=2.5)
ap.add_argument("--min-run", type=int, default=2)
a = ap.parse_args()
AUTOS = [str(i) for i in range(1, 23)]

# expected depth by GC from genome bins
bins = list(csv.DictReader(open(a.bins), delimiter="\t"))
bgc = np.array([float(b["gc"]) for b in bins]); bdp = np.array([float(b["depth"]) for b in bins])
bok = np.array([b["chrom"] in AUTOS and float(b["map_frac"]) >= a.min_map and int(float(b["blacklist_bp"])) == 0
                and float(b["depth"]) > 0 for b in bins])
bstr = np.round(bgc / 0.005).astype(int)


def expected(gc):
    s, w = round(gc / 0.005), 0
    while True:
        sel = bok & (np.abs(bstr - s) <= w)
        if sel.sum() >= 30 or w > 40:
            return float(np.median(bdp[sel])) if sel.any() else float(np.median(bdp[bok]))
        w += 1


ex_annot = {(r["chrom"], r["start"], r["end"]): r for r in csv.DictReader(open(a.exons), delimiter="\t")}
exons = []
for l in gzip.open(a.exon_depth, "rt"):
    c, s, e, name, m = l.rstrip("\n").split("\t")[:5]
    an = ex_annot.get((c, s, e), {})
    gene, ex = name.split("|")
    gcv = float(an.get("gc", "nan")); mapv = float(an.get("map_frac", "nan"))
    exp = expected(gcv) if gcv == gcv else float("nan")
    exons.append({"gene": gene, "exon": ex, "chrom": c, "start": int(s), "end": int(e), "gc": gcv, "map_frac": mapv,
                  "depth": float(m), "expected": exp, "cn": 2 * float(m) / exp if exp else float("nan"),
                  "usable": mapv >= a.min_map})

arms = {r["arm"]: r for r in json.load(open(a.arms))}
panel = []
for l in open(a.panel):
    if l.startswith("#") or not l.strip():
        continue
    f = l.rstrip("\n").split("\t")
    panel.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]), "gene": f[3], "group": f[5]})

genes = []
by_gene = defaultdict(list)
for x in exons:
    by_gene[x["gene"]].append(x)
for g, xs in sorted(by_gene.items()):
    xs.sort(key=lambda x: int(re.sub(r"\D", "", x["exon"])))
    use = [x for x in xs if x["usable"]]
    cns = np.array([x["cn"] for x in use])
    runs = []
    for label, test in (("deletion-like", lambda v: v < a.del_cn), ("duplication-like", lambda v: v > a.dup_cn)):
        cur = []
        for x in use:
            if test(x["cn"]):
                cur.append(x["exon"])
            else:
                if len(cur) >= a.min_run:
                    runs.append(f"{label}: {cur[0]}-{cur[-1]} ({len(cur)} exons)")
                cur = []
        if len(cur) >= a.min_run:
            runs.append(f"{label}: {cur[0]}-{cur[-1]} ({len(cur)} exons)")
    arm_match = None
    if xs:
        ar = next((A for A, r in arms.items() if A[:-1] == xs[0]["chrom"]), None)
        pa = [p for p in panel if p["gene"] == g]
        arm_lrrs = [r.get("lrr_mean") for A, r in arms.items() if A[:-1] == xs[0]["chrom"] and r.get("significant")]
        arm_match = "gene lies on an arm with significant mosaic imbalance" if arm_lrrs else None
    genes.append({"gene": g, "group": next((p["group"] for p in panel if p["gene"] == g), ""),
                  "exons": len(xs), "exons_usable": len(use),
                  "median_cn": float(np.median(cns)) if len(cns) else None,
                  "min_exon_cn": float(cns.min()) if len(cns) else None, "max_exon_cn": float(cns.max()) if len(cns) else None,
                  "flag": "; ".join(runs) or "none", "arm_context": arm_match or ""})

# Delly SVs overlapping the panel
svs, sv_counts = [], defaultdict(int)
with gzip.open(a.delly, "rt") as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if f[6] != "PASS":
            continue
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        t = info.get("SVTYPE", "?")
        gt = fmt.get("GT", "./.")
        if gt in ("0/0", "./."):
            continue
        sv_counts[t] += 1
        c1, p1 = f[0], int(f[1])
        c2, p2 = info.get("CHR2", c1), int(info.get("END", p1))
        spans = [(c1, p1, p2)] if t in ("DEL", "DUP", "INV", "INS") and c1 == c2 else [(c1, p1, p1), (c2, p2, p2)]
        hits = {p["gene"] for p in panel for (c, s, e) in spans if p["chrom"] == c and s <= p["end"] and e >= p["start"]}
        if not hits:
            continue
        ex_hit = [f"{x['gene']}:{x['exon']}" for x in exons for (c, s, e) in spans if x["chrom"] == c and s <= x["end"] and e >= x["start"]]
        svs.append({"id": f[2], "type": t, "chrom": c1, "pos": p1, "chr2": c2, "end": p2,
                    "size": (p2 - p1) if c1 == c2 and t != "BND" else "", "genes": ";".join(sorted(hits)),
                    "exons_hit": ";".join(ex_hit[:20]) + (" ..." if len(ex_hit) > 20 else ""),
                    "genotype": gt, "PE": info.get("PE", ""), "SR": info.get("SR", ""),
                    "precise": "PRECISE" in info, "qual": f[5]})


def write(path, rows):
    if not rows:
        open(path, "w").close(); return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


write(a.prefix + ".gene_cn.tsv", genes)
write(a.prefix + ".exon_cn.tsv", [{k: v for k, v in x.items()} for x in exons])
write(a.prefix + ".sv_panel.tsv", svs)
summary = {
    "genes": len(genes), "exons": len(exons), "exons_usable": sum(1 for x in exons if x["usable"]),
    "exon_cn_median_all": float(np.median([x["cn"] for x in exons if x["usable"]])) if exons else None,
    "genes_flagged_cn": [g["gene"] for g in genes if g["flag"] != "none"],
    "delly_pass_nonref_by_type": dict(sv_counts),
    "sv_overlapping_panel": len(svs),
    "sv_hitting_panel_exons": sum(1 for s in svs if s["exons_hit"]),
    "thresholds": {"del_cn": a.del_cn, "dup_cn": a.dup_cn, "min_run_exons": a.min_run, "min_map": a.min_map},
}
json.dump(summary, open(a.prefix + ".cnv_sv.summary.json", "w"), indent=2)
print(json.dumps(summary))
