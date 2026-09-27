#!/usr/bin/env python3
"""Read-level check of each exon-run copy-number candidate from gene_cn.py.

An exon-run call rests on mean depth over short exons, which is sensitive to exon-scale GC and
capture of AT-rich sequence. A real heterozygous deletion spanning those exons also removes the
intron between them and produces a depth step with sharp edges. So for each candidate run:

  span   = first exon start .. last exon end of the run
  flank  = `flank` bp either side of the span (excluding the span)
  depth  = samtools depth -a, MAPQ >= min_mapq, CRAM decoded with the explicit reference,
           averaged in `window` bp windows
  ratio  = mean depth over the span / median of the flank window means

  consistent with deletion   ratio <= del_max  (het deletion expected ~0.5)
  consistent with duplication ratio >= dup_min (one extra copy expected ~1.5)
  not confirmed              ratio within [no_step_lo, no_step_hi] and no Delly call of the matching
                             type overlapping the span; classified as noise
  indeterminate              anything else

Usage: cn_candidate_check.py --gene-cn g.tsv --exon-cn e.tsv --sv-panel sv.tsv --cram x.cram
                             --fasta ref.fa --out-tsv out.tsv --out-json out.json
"""
import argparse
import csv
import json
import re
import statistics
import subprocess

ap = argparse.ArgumentParser()
for x in ("gene-cn", "exon-cn", "sv-panel", "cram", "fasta", "out-tsv", "out-json"):
    ap.add_argument(f"--{x}", required=True)
ap.add_argument("--flank", type=int, default=2000)
ap.add_argument("--window", type=int, default=250)
ap.add_argument("--min-mapq", type=int, default=20)
ap.add_argument("--del-max", type=float, default=0.65)
ap.add_argument("--dup-min", type=float, default=1.35)
ap.add_argument("--no-step-lo", type=float, default=0.8)
ap.add_argument("--no-step-hi", type=float, default=1.2)
a = ap.parse_args()


def rows(p):
    try:
        return list(csv.DictReader(open(p), delimiter="\t"))
    except (FileNotFoundError, StopIteration):
        return []


RUN = re.compile(r"(deletion-like|duplication-like): (exon\d+)-(exon\d+)")


def classify(ratio, kind, delly_hit):
    if kind == "deletion-like" and ratio <= a.del_max:
        return "consistent with deletion"
    if kind == "duplication-like" and ratio >= a.dup_min:
        return "consistent with duplication"
    if a.no_step_lo <= ratio <= a.no_step_hi and not delly_hit:
        return "not confirmed"
    return "indeterminate"


def depth(chrom, start, end):
    """Per-base depth over [start, end) (0-based half-open), zero-filled."""
    out = subprocess.run(["samtools", "depth", "-a", "-Q", str(a.min_mapq), "--reference", a.fasta,
                          "-r", f"{chrom}:{start + 1}-{end}", a.cram], check=True, capture_output=True, text=True).stdout
    d = {int(f[1]) - 1: int(f[2]) for f in (l.split("\t") for l in out.splitlines())}
    return [d.get(p, 0) for p in range(start, end)]


def windows(vals, start, lo, hi):
    """Window means over [lo, hi) of a depth vector that begins at `start`."""
    out = []
    for w in range(lo, hi, a.window):
        seg = vals[w - start:min(w + a.window, hi) - start]
        if seg:
            out.append(sum(seg) / len(seg))
    return out


exons = {(r["gene"], r["exon"]): r for r in rows(a.exon_cn)}
svs = rows(a.sv_panel)
results = []
for g in rows(a.gene_cn):
    for m in RUN.finditer(g.get("cn_call") or ""):
        kind, e1, e2 = m.groups()
        x1, x2 = exons[(g["gene"], e1)], exons[(g["gene"], e2)]
        chrom, s, e = x1["chrom"], int(x1["start"]), int(x2["end"])
        lo, hi = max(0, s - a.flank), e + a.flank
        d = depth(chrom, lo, hi)
        span = d[s - lo:e - lo]
        span_win = windows(d, lo, s, e)
        flank_win = windows(d, lo, lo, s) + windows(d, lo, e, hi)
        flank_med = statistics.median(flank_win)
        ratio = (sum(span) / len(span)) / flank_med if flank_med else float("nan")
        want = "DEL" if kind == "deletion-like" else "DUP"
        delly = [v["id"] for v in svs if v["type"] == want and v["chrom"] == chrom
                 and int(v["pos"]) <= e and int(v["end"]) >= s]
        results.append({
            "gene": g["gene"], "call": m.group(0), "region": f"{chrom}:{s + 1}-{e}",
            "span_mean_depth": round(sum(span) / len(span), 2),
            "span_window_range": f"{min(span_win):.0f}-{max(span_win):.0f}",
            "flank_window_median": round(flank_med, 2),
            "flank_window_range": f"{min(flank_win):.0f}-{max(flank_win):.0f}",
            "span_to_flank_ratio": round(ratio, 3),
            "delly_overlapping": ";".join(delly),
            "verdict": classify(ratio, kind, bool(delly)),
        })

fields = ["gene", "call", "region", "span_mean_depth", "span_window_range", "flank_window_median",
          "flank_window_range", "span_to_flank_ratio", "delly_overlapping", "verdict"]
with open(a.out_tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(results)
summary = {
    "method": f"samtools depth -a -Q {a.min_mapq}; {a.window} bp windows; +/-{a.flank} bp flanks",
    "thresholds": {"deletion_ratio_max": a.del_max, "duplication_ratio_min": a.dup_min,
                   "no_step_ratio": [a.no_step_lo, a.no_step_hi]},
    "candidates": results,
    "open_candidates": [r["gene"] for r in results if r["verdict"] != "not confirmed"],
    "resolved_as_noise": [r["gene"] for r in results if r["verdict"] == "not confirmed"],
}
json.dump(summary, open(a.out_json, "w"), indent=2)
print(json.dumps(summary, indent=1))
