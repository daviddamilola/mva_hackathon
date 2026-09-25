#!/usr/bin/env python3
"""Compare VCF header contigs with a FASTA index (names, lengths, order).
Reads only '##contig' header lines; never reads genotype records.
Usage: check_contigs.py <vcf|vcf.gz|header> <fasta.fai> <reference_id> <out.json> [--no-strict]"""
import gzip
import json
import re
import sys

vcf_path, fai_path, ref_id, out_path = sys.argv[1:5]
strict = "--no-strict" not in sys.argv

opener = gzip.open if vcf_path.endswith(".gz") else open
vcf = []
with opener(vcf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM") or not line.startswith("##"):
            break
        m = re.match(r"##contig=<ID=([^,>]+),length=(\d+)", line)
        if m:
            vcf.append((m.group(1), int(m.group(2))))

with open(fai_path) as fh:
    fai = [(f[0], int(f[1])) for f in (l.rstrip("\n").split("\t") for l in fh) if f[0]]

only_vcf, only_fai = sorted(set(vcf) - set(fai)), sorted(set(fai) - set(vcf))
res = {
    "reference": ref_id,
    "n_contigs_vcf": len(vcf),
    "n_contigs_fasta": len(fai),
    "identical_names_lengths": not only_vcf and not only_fai,
    "identical_order": vcf == fai,
    "only_in_vcf": [list(x) for x in only_vcf[:50]],
    "only_in_fasta": [list(x) for x in only_fai[:50]],
}
with open(out_path, "w") as out:
    json.dump(res, out, indent=2)
if strict and not res["identical_names_lengths"]:
    sys.exit(f"Reference contigs do not match the VCF header; see {out_path}")
print(json.dumps(res))
