#!/usr/bin/env python3
"""Collect every position named by Open Targets PGx rows (genotypeId / variantId, GRCh38 no-chr)
that is not already genotyped in the stage-1a force-call VCF, for force-calling from the CRAM.

The drug list defines the panel: positions come from all returned rows (all evidence levels), so
exploratory level-3 rows are genotyped too, even though they are not asserted.

Usage: ot_sites.py --genotyped stage1a.vcf.gz --out-tsv sites.tsv --out-bed regions.bed <ot json>...
"""
import argparse
import gzip
import json
import re

ap = argparse.ArgumentParser()
ap.add_argument("--genotyped", required=True)
ap.add_argument("--out-tsv", required=True)
ap.add_argument("--out-bed", required=True)
ap.add_argument("--primary-contigs", default=",".join([str(i) for i in range(1, 23)] + ["X", "Y"]))
ap.add_argument("json", nargs="+")
a = ap.parse_args()
primary = set(a.primary_contigs.split(","))

have = set()
with gzip.open(a.genotyped, "rt") as fh:
    for line in fh:
        if not line.startswith("#"):
            f = line.split("\t", 2)
            have.add((f[0], int(f[1])))


def rows_of(doc):
    if "rows" in doc:
        yield from doc["rows"]
    for g in (doc.get("genes") or {}).values():
        yield from g.get("rows") or []


want, skipped = set(), 0
for path in a.json:
    for r in rows_of(json.load(open(path))):
        for key in (r.get("genotypeId"), r.get("variantId")):
            if not key:
                continue
            m = re.match(r"^([0-9XYM]+)_(\d+)_", key)
            if not m or m.group(1) not in primary:
                skipped += 1
                continue
            pos = (m.group(1), int(m.group(2)))
            if pos not in have:
                want.add(pos)

order = {c: i for i, c in enumerate([str(i) for i in range(1, 23)] + ["X", "Y"])}
sites = sorted(want, key=lambda p: (order.get(p[0], 99), p[1]))
with open(a.out_tsv, "w") as out:
    for c, p in sites:
        out.write(f"{c}\t{p}\n")
with open(a.out_bed, "w") as out:
    for c, p in sites:
        out.write(f"{c}\t{max(0, p - 101)}\t{p + 100}\n")
print(json.dumps({"new_positions": len(sites), "already_genotyped": len(have), "unparsed_or_nonprimary": skipped}))
