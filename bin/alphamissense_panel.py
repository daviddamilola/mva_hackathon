#!/usr/bin/env python3
"""Filter the AlphaMissense hg38 TSV (BGZF, no index published) to rows inside the panel BED and
write a minimal VCF (no-chr contigs) for VEP --custom. The MD5 of the compressed input is checked
against the value Google Cloud Storage publishes (x-goog-hash, base64); a mismatch fails, so a
truncated or corrupted download cannot yield a silently incomplete subset.

Usage: alphamissense_panel.py --input AlphaMissense_hg38.tsv.gz --bed panel.bed --md5-b64 <b64> --out am_panel.vcf
       (--input - reads stdin)
"""
import argparse
import base64
import bisect
import hashlib
import sys
import zlib

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="-")
ap.add_argument("--bed", required=True)
ap.add_argument("--md5-b64", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

regions = {}
for line in open(a.bed):
    if line.startswith("#") or not line.strip():
        continue
    c, s, e = line.split("\t")[:3]
    regions.setdefault(c, []).append((int(s), int(e)))
starts = {}
for c, v in regions.items():   # merge overlapping windows so the bisect lookup is exact
    merged = []
    for s, e in sorted(v):
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    starts[c] = merged


def in_panel(chrom, pos):
    iv = starts.get(chrom)
    if not iv:
        return False
    i = bisect.bisect_right(iv, (pos, float("inf"))) - 1
    return i >= 0 and iv[i][0] < pos <= iv[i][1]


src = sys.stdin.buffer if a.input == "-" else open(a.input, "rb")
md5 = hashlib.md5()
decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
buf, kept, total = b"", [], 0
while True:
    chunk = src.read(1 << 22)
    if not chunk:
        break
    md5.update(chunk)
    data = chunk
    while data:
        out = decomp.decompress(data)
        buf += out
        data = decomp.unused_data
        if decomp.eof:  # BGZF = concatenated gzip members
            decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
        else:
            break
    *lines, buf = buf.split(b"\n")
    for ln in lines:
        if not ln or ln.startswith(b"#"):
            continue
        f = ln.decode().split("\t")
        # CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class
        chrom = f[0].replace("chr", "")
        total += 1
        if in_panel(chrom, int(f[1])):
            kept.append((chrom, int(f[1]), f[2], f[3], f[6], f[7], f[8], f[9]))

observed = base64.b64encode(md5.digest()).decode()
if observed != a.md5_b64:
    sys.exit(f"AlphaMissense stream MD5 mismatch: expected {a.md5_b64}, observed {observed}")

order = {c: i for i, c in enumerate([str(i) for i in range(1, 23)] + ["X", "Y"])}
kept.sort(key=lambda r: (order.get(r[0], 99), r[1], r[3]))
with open(a.out, "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    out.write("##source=AlphaMissense_hg38.tsv.gz (Cheng et al. 2023), panel subset\n")
    for c in sorted({r[0] for r in kept}, key=lambda c: order.get(c, 99)):
        out.write(f"##contig=<ID={c}>\n")
    for k, t, d in [("AM_PATHOGENICITY", "Float", "AlphaMissense pathogenicity"), ("AM_CLASS", "String", "AlphaMissense class"),
                    ("AM_TRANSCRIPT", "String", "AlphaMissense transcript"), ("AM_PROTEIN_VARIANT", "String", "AlphaMissense protein variant")]:
        out.write(f'##INFO=<ID={k},Number=1,Type={t},Description="{d}">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for c, p, r, al, tx, pv, score, cls in kept:
        out.write(f"{c}\t{p}\t.\t{r}\t{al}\t.\t.\tAM_PATHOGENICITY={score};AM_CLASS={cls};AM_TRANSCRIPT={tx};AM_PROTEIN_VARIANT={pv}\n")
print(f"AlphaMissense rows scanned: {total}; kept in panel: {len(kept)}; md5 ok ({observed})")
