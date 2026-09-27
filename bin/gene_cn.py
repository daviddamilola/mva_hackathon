#!/usr/bin/env python3
"""Panel-gene copy number (exon read depth) and structural variants (Delly) overlapping the panel.

EXON COPY NUMBER
  r = log2(observed / expected), where expected = median depth of autosomal genome bins (mappable
  >= min_map, no blacklist) in the same GC stratum. GC bias is stronger at exon scale than at 100 kb,
  so r is re-centred against exon GC across all usable panel exons (median per 5% GC stratum, most
  exons being CN 2). Noise depends on exon length: z = r / (sigma0 * sqrt(L0 / L)), where sigma0 is
  estimated robustly (MAD) at reference length L0 = 150 bp. A CN call needs a run of >= min_run
  consecutive usable exons all with z <= -z_min (deletion-like) or >= z_min (duplication-like);
  CN = 2 * 2^mean(r) over the run. A single-exon outlier (|z| >= z_min) is reported as low confidence
  only. (First version used fixed CN thresholds 1.5/2.5: with per-exon SD ~0.34 that flagged 7/32
  genes by chance.)

STRUCTURAL VARIANTS
  Delly PASS, non-reference calls. A gene counts as hit when:
    DEL / DUP  the event overlaps the gene +/- flank
    INV / BND  a BREAKPOINT falls in the gene +/- flank (an inversion containing a whole gene does
               not disrupt it)
    INS        the insertion point falls in the gene +/- flank
  DEL/DUP are checked against depth: >= 100 kb by mean bin LRR over the span (het DEL expected -1,
  DUP +0.58; supported if the LRR is at least half-way there); < 100 kb by the z-scores of exons inside
  the span (supported if mean z <= -z_min / >= z_min); otherwise 'depth not assessable'.
  Evidence tier: 'high' = PRECISE with PE >= 3 and SR >= 3; else 'low'. Only high-tier calls not
  contradicted by depth count as supported. (First version counted any overlap: 105 'hits' were
  imprecise 100-223 Mb calls between repeat-rich chromosome ends, invisible in depth.)

Usage: gene_cn.py --exon-depth exons.regions.bed.gz --exons exons.annot.tsv --bins bins_with_depth.tsv
                  --delly sv.vcf.gz --panel panel.bed --arms arms.json --prefix out
"""
import argparse
import csv
import gzip
import json
import math
import re
from collections import defaultdict

import numpy as np

ap = argparse.ArgumentParser()
for x in ("exon-depth", "exons", "bins", "lrr-bins", "delly", "panel", "arms", "prefix"):
    ap.add_argument(f"--{x}", required=True)
ap.add_argument("--min-map", type=float, default=0.9)
ap.add_argument("--z-min", type=float, default=3.0)
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
    c, s_, e_, name, m = l.rstrip("\n").split("\t")[:5]
    an = ex_annot.get((c, s_, e_), {})
    gene, ex = name.split("|")
    gcv = float(an.get("gc", "nan")); mapv = float(an.get("map_frac", "nan"))
    exp = expected(gcv) if gcv == gcv else float("nan")
    dep = float(m)
    r = math.log2(dep / exp) if (exp and dep > 0) else float("nan")
    exons.append({"gene": gene, "exon": ex, "chrom": c, "start": int(s_), "end": int(e_), "length": int(e_) - int(s_),
                  "gc": gcv, "map_frac": mapv, "depth": dep, "expected": exp, "log2r_raw": r,
                  "usable": bool(mapv >= a.min_map and c in AUTOS and r == r)})

# exon-scale GC re-centring across usable panel exons
use_all = [x for x in exons if x["usable"]]
strata = defaultdict(list)
for x in use_all:
    strata[round(x["gc"] / 0.05)].append(x["log2r_raw"])
def gc_median(gc):
    k, w = round(gc / 0.05), 0
    while True:
        vals = [v for kk in range(k - w, k + w + 1) for v in strata.get(kk, [])]
        if len(vals) >= 15 or w > 6:
            return float(np.median(vals)) if vals else 0.0
        w += 1
L0 = 150.0
for x in exons:
    x["log2r"] = x["log2r_raw"] - gc_median(x["gc"]) if x["usable"] else float("nan")
scaled = np.array([x["log2r"] * math.sqrt(x["length"] / L0) for x in use_all])
sigma0 = float(1.4826 * np.median(np.abs(scaled - np.median(scaled))))
for x in exons:
    x["z"] = x["log2r"] / (sigma0 * math.sqrt(L0 / x["length"])) if x["usable"] else float("nan")
    x["cn"] = 2 * 2 ** x["log2r"] if x["usable"] else float("nan")

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
    calls, singles = [], []
    for label, test in (("deletion-like", lambda z: z <= -a.z_min), ("duplication-like", lambda z: z >= a.z_min)):
        cur = []
        for x in use + [None]:
            if x is not None and test(x["z"]):
                cur.append(x); continue
            if len(cur) >= a.min_run:
                cn_run = 2 * 2 ** float(np.mean([y["log2r"] for y in cur]))
                calls.append(f"{label}: {cur[0]['exon']}-{cur[-1]['exon']} ({len(cur)} exons, CN~{cn_run:.2f}, "
                             f"z {min(y['z'] for y in cur):.1f}..{max(y['z'] for y in cur):.1f})")
            elif len(cur) == 1:
                singles.append(f"{cur[0]['exon']} z={cur[0]['z']:.1f}")
            cur = []
    arm_sig = [A for A, r in arms.items() if xs and A[:-1] == xs[0]["chrom"] and r.get("significant")]
    lr = np.array([x["log2r"] for x in use])
    genes.append({"gene": g, "group": next((p["group"] for p in panel if p["gene"] == g), ""),
                  "exons": len(xs), "exons_usable": len(use),
                  "gene_cn": float(2 * 2 ** np.median(lr)) if len(lr) else None,
                  "cn_call": "; ".join(calls) or "none",
                  "single_exon_outliers_low_confidence": "; ".join(singles),
                  "arm_context": "on an arm with an arm-level mosaic event" if arm_sig else ""})

# depth lookup for SV confirmation
lrr_bins = defaultdict(list)
for l in gzip.open(a.lrr_bins, "rt"):
    if l.startswith("chrom"):
        continue
    f = l.rstrip("\n").split("\t")
    if f[5] == "1" and f[7] != "":
        lrr_bins[f[0]].append((int(f[1]), int(f[2]), float(f[7])))

def span_lrr(c, s_, e_):
    v = [lr for (bs, be, lr) in lrr_bins.get(c, []) if bs >= s_ and be <= e_]
    return (float(np.mean(v)), len(v)) if v else (None, 0)

# Delly SVs and the panel
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
        t, gt = info.get("SVTYPE", "?"), fmt.get("GT", "./.")
        if gt in ("0/0", "./."):
            continue
        sv_counts[t] += 1
        c1, p1 = f[0], int(f[1])
        c2, p2 = info.get("CHR2", c1), int(info.get("END", p1))
        if t in ("DEL", "DUP") and c1 == c2:
            test = [(c1, p1, p2)]
        else:  # INV, BND, INS: breakpoints only
            test = [(c1, p1, p1), (c2, p2, p2)]
        hits = sorted({p["gene"] for p in panel for (c, s_, e_) in test if p["chrom"] == c and s_ <= p["end"] and e_ >= p["start"]})
        if not hits:
            continue
        pe, sr = int(info.get("PE", 0) or 0), int(info.get("SR", 0) or 0)
        tier = "high" if ("PRECISE" in info and pe >= 3 and sr >= 3) else "low"
        size = (p2 - p1) if (c1 == c2 and t != "BND") else None
        depth = "not applicable"
        if t in ("DEL", "DUP") and c1 == c2:
            if size >= 100_000:
                m, nb = span_lrr(c1, p1, p2)
                if m is None:
                    depth = "not assessable"
                else:
                    target = -1.0 if t == "DEL" else 0.585
                    ok = (m <= target / 2) if t == "DEL" else (m >= target / 2)
                    depth = f"{'supported' if ok else 'contradicted'} (mean LRR {m:+.3f} over {nb} bins)"
            else:
                zs = [x["z"] for x in exons if x["usable"] and x["chrom"] == c1 and x["start"] < p2 and x["end"] > p1]
                if zs:
                    mz = float(np.mean(zs))
                    ok = (mz <= -a.z_min) if t == "DEL" else (mz >= a.z_min)
                    depth = f"{'supported' if ok else 'not supported'} (mean exon z {mz:+.1f}, {len(zs)} exons)"
                else:
                    depth = "not assessable (no exon inside)"
        ex_hit = [f"{x['gene']}:{x['exon']}" for x in exons for (c, s_, e_) in test
                  if x["chrom"] == c and s_ <= x["end"] and e_ >= x["start"]]
        supported = tier == "high" and not depth.startswith(("contradicted", "not supported"))
        svs.append({"id": f[2], "type": t, "chrom": c1, "pos": p1, "chr2": c2, "end": p2, "size": size or "",
                    "genes": ";".join(hits), "exons_hit": ";".join(ex_hit[:20]) + (" ..." if len(ex_hit) > 20 else ""),
                    "genotype": gt, "PE": pe, "SR": sr, "precise": "PRECISE" in info, "tier": tier,
                    "depth_check": depth, "supported": supported})

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
    "exon_noise_sigma_log2_at_150bp": sigma0,
    "gene_cn_median": float(np.median([g["gene_cn"] for g in genes if g["gene_cn"] is not None])),
    "genes_with_cn_call": [g["gene"] for g in genes if g["cn_call"] != "none"],
    "genes_with_single_exon_outliers": [g["gene"] for g in genes if g["single_exon_outliers_low_confidence"]],
    "delly_pass_nonref_by_type": dict(sv_counts),
    "sv_touching_panel": len(svs),
    "sv_touching_panel_high_tier": sum(1 for s_ in svs if s_["tier"] == "high"),
    "sv_supported": [f"{s_['type']} {s_['chrom']}:{s_['pos']}-{s_['end']} ({s_['genes']})" for s_ in svs if s_["supported"]],
    "sv_contradicted_by_depth": sum(1 for s_ in svs if s_["depth_check"].startswith(("contradicted", "not supported"))),
    "thresholds": {"z_min": a.z_min, "min_run_exons": a.min_run, "min_map": a.min_map, "ref_exon_length": L0},
}
json.dump(summary, open(a.prefix + ".cnv_sv.summary.json", "w"), indent=2)
print(json.dumps(summary))
