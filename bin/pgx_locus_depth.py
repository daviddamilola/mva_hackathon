#!/usr/bin/env python3
"""Windowed read depth over a pharmacogene locus, for the report's copy-number figure.

Two tracks are written for every window: all reads, and reads with MAPQ >= 20. The locus (CYP2D6 with
its CYP2D7 paralog) is highly homologous, so the MAPQ-filtered track dips wherever reads cannot be placed
uniquely; showing both is the honest view, and the difference between them is itself informative.
Depth is reported raw; the figure normalises to the sample's autosomal mean from genome QC.

Usage: pgx_locus_depth.py --cram S.cram --fasta ref.fa --region 22:42116498-42155810 --window 250 --out S.tsv
"""
import argparse
import subprocess

ap = argparse.ArgumentParser()
ap.add_argument("--cram", required=True)
ap.add_argument("--fasta", required=True)
ap.add_argument("--region", required=True)
ap.add_argument("--window", type=int, default=250)
ap.add_argument("--out", required=True)
a = ap.parse_args()

chrom, span = a.region.split(":")
start, end = (int(x) for x in span.split("-"))


def depth(min_mq):
    cmd = ["samtools", "depth", "-a", "-Q", str(min_mq), "-r", a.region, "--reference", a.fasta, a.cram]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    d = {}
    for line in out.splitlines():
        c, p, v = line.split("\t")
        d[int(p)] = int(v)
    return d


all_reads, mq20 = depth(0), depth(20)
with open(a.out, "w") as fh:
    fh.write("chrom\tstart\tend\tdepth_all\tdepth_mq20\n")
    for ws in range(start, end + 1, a.window):
        we = min(ws + a.window - 1, end)
        n = we - ws + 1
        fh.write(f"{chrom}\t{ws}\t{we}\t{sum(all_reads.get(p, 0) for p in range(ws, we + 1)) / n:.2f}\t"
                 f"{sum(mq20.get(p, 0) for p in range(ws, we + 1)) / n:.2f}\n")
