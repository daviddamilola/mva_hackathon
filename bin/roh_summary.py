#!/usr/bin/env python3
"""Summarise bcftools roh output and cross-check it with heterozygous-SNV density.

F_ROH at >= 0.5, 1, 1.5 and 5 Mb = total ROH length / autosomal length (from the .fai). The
cross-check counts PASS biallelic SNVs per 1 Mb window (windows with >= 300 SNVs) and reports windows
with < 10% of the median het count and the longest consecutive run of such windows.
Reference points: offspring of first cousins F ~ 0.0625; second cousins ~ 0.016.

Usage: roh_summary.py --roh roh.txt --vcf sample.vcf.gz --fai ref.fa.fai --out summary.json
"""
import argparse
import gzip
import json
import statistics
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--roh", required=True)
ap.add_argument("--vcf", required=True)
ap.add_argument("--fai", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
AUTOS = [str(i) for i in range(1, 23)]
auto_len = sum(int(l.split("\t")[1]) for l in open(a.fai) if l.split("\t")[0] in AUTOS)

segs = [l.rstrip("\n").split("\t") for l in open(a.roh) if l.startswith("RG")]
lens = [int(s[5]) for s in segs]
froh = {}
for thr in (500_000, 1_000_000, 1_500_000, 5_000_000):
    sel = [x for x in lens if x >= thr]
    froh[f">={thr / 1e6:g}Mb"] = {"n": len(sel), "total_mb": round(sum(sel) / 1e6, 3), "F_ROH": round(sum(sel) / auto_len, 5)}
longest = max(segs, key=lambda s: int(s[5])) if segs else None

het, tot = Counter(), Counter()
with gzip.open(a.vcf, "rt") as fh:
    for l in fh:
        if l[0] == "#":
            continue
        f = l.split("\t", 10)
        if f[0] not in AUTOS or f[6] != "PASS" or len(f[3]) != 1 or len(f[4]) != 1:
            continue
        w = (f[0], int(f[1]) // 1_000_000)
        tot[w] += 1
        if f[9].split(":", 1)[0] in ("0/1", "0|1", "1|0"):
            het[w] += 1
ws = sorted((w for w in tot if tot[w] >= 300), key=lambda w: (int(w[0]), w[1]))
med = statistics.median([het[w] for w in ws])
low = [w for w in ws if het[w] < 0.1 * med]
best = cur = 0
prev = None
for w in low:
    cur = cur + 1 if prev and w[0] == prev[0] and w[1] == prev[1] + 1 else 1
    best, prev = max(best, cur), w

summary = {
    "method": "bcftools roh -G30 --AF-dflt 0.4, autosomes; cross-check: het density in 1 Mb windows",
    "autosomal_length": auto_len, "F_ROH": froh,
    "longest_roh": {"region": f"{longest[2]}:{longest[3]}-{longest[4]}", "mb": round(int(longest[5]) / 1e6, 3)} if longest else None,
    "het_density_windows": len(ws), "het_density_median_per_mb": med,
    "low_het_windows": len(low), "longest_low_het_run_mb": best,
    "interpretation_reference": "offspring of first cousins F~0.0625; second cousins F~0.016",
}
json.dump(summary, open(a.out, "w"), indent=2)
print(json.dumps(summary["F_ROH"]))
