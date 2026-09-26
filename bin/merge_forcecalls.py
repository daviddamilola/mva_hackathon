#!/usr/bin/env python3
"""Merge constrained and unconstrained force-calls at PharmCAT positions.

Constrained calls (bcftools call -C alleles -A) genotype every PharmCAT definition allele
explicitly, including hom-ref at indel sites. However, they can only test the alleles PharmCAT
defines: a real non-reference allele that PharmCAT does not define is written as 0/0.

The unconstrained call (bcftools call -m) sees the observed allele. At any position where the
unconstrained genotype carries a non-reference allele absent from the constrained ALT list, the
unconstrained record replaces the constrained one. PharmCAT then sees the real allele and raises
its own undocumented-variation warning instead of reporting reference. Every such site is written
to off_definition.tsv.

Usage: merge_forcecalls.py <constrained.vcf.gz> <unconstrained.vcf.gz> <merged.vcf> <off_definition.tsv>
"""
import gzip
import sys


def read_vcf(path):
    header, records = [], {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            records.setdefault((f[0], int(f[1])), []).append(f)
    return header, records


def gt_alleles(f):
    fmt = f[8].split(":")
    gt = dict(zip(fmt, f[9].split(":"))).get("GT", ".")
    return [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]


def field(f, key):
    return dict(zip(f[8].split(":"), f[9].split(":"))).get(key, ".")


con_path, unc_path, out_vcf, out_tsv = sys.argv[1:5]
con_header, con = read_vcf(con_path)
_, unc = read_vcf(unc_path)

off_rows = []
with open(out_vcf, "w") as out:
    out.writelines(con_header)
    for key in con:  # insertion order = input order, already sorted by header contig order
        for cf in con[key]:
            replacement = None
            con_alts = set(cf[4].split(",")) - {"."}
            for uf in unc.get(key, []):
                alleles = gt_alleles(uf)
                if not alleles or all(a == "0" for a in alleles):
                    continue
                unc_alts = uf[4].split(",")
                observed = {unc_alts[int(a) - 1] for a in alleles if a != "0"}
                # compare on the allele difference relative to REF, so representation differences
                # (same event, different REF padding) are not flagged
                def event(ref, alt):
                    return (ref, alt) if len(ref) == len(alt) == 1 else (ref[1:], alt[1:], len(ref) - len(alt))
                con_events = {event(cf[3], a) for a in con_alts}
                novel = [a for a in observed if event(uf[3], a) not in con_events]
                if novel:
                    replacement = uf
                    off_rows.append([uf[0], uf[1], cf[2], cf[3], cf[4], uf[3], ",".join(novel),
                                     field(uf, "GT"), field(uf, "AD"), field(uf, "DP"), field(cf, "GT")])
                    break
            out.write("\t".join(replacement or cf) + "\n")

with open(out_tsv, "w") as out:
    out.write("chrom\tpos\tid\tdefinition_ref\tdefinition_alts\tobserved_ref\tobserved_undefined_alt\t"
              "observed_gt\tobserved_ad\tobserved_dp\tconstrained_gt\n")
    for r in off_rows:
        out.write("\t".join(r) + "\n")
print(f"{len(off_rows)} PharmCAT position(s) carry an allele outside the PharmCAT definition")
