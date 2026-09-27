#!/usr/bin/env python3
"""MT-RNR1 (12S rRNA) genotype for the CPIC aminoglycoside guideline, with heteroplasmy.

Input is a bcftools mpileup over the whole MT-RNR1 gene (chrM, rCRS numbering) with FORMAT/AD, so
every position has observed allele depths. Calls are made from allele fractions, not from a diploid
genotype model: mtDNA is haploid per molecule, and heteroplasmy is a fraction, not a genotype.

  AF >= homoplasmic_min                 homoplasmic
  heteroplasmic_min <= AF < homo_min    heteroplasmic (reported with the fraction)
  low_level_min <= AF < het_min         low-level; reported, not called (NUMT reads and error)
  otherwise                             not detected

Each CPIC allele is matched by position and allele; the CPIC REF must equal the FASTA base (checked).
Deletions in the m.956-965 C-tract are matched by position window, since their representation in a
homopolymer is not unique. All non-reference sites in the gene are listed, CPIC or not: the
near-universal rCRS differences m.750A>G and m.1438A>G act as a positive control that the region was
genotyped.

The PharmCAT outside call is the highest-risk CPIC allele called (homoplasmic or heteroplasmic), else
'Reference'. A heteroplasmic call is passed on and flagged, because the guideline's risk statements are
made for carriers and give no safe heteroplasmy threshold.

Usage: mt_rnr1_call.py --pileup mt.vcf.gz --alleles cpic.tsv --fasta ref.fa --contig M --sample S
                       --out-json out.json --out-outside outside.tsv
"""
import argparse
import csv
import gzip
import json
import subprocess

ap = argparse.ArgumentParser()
for x in ("pileup", "alleles", "fasta", "contig", "sample", "out-json", "out-outside"):
    ap.add_argument(f"--{x}", required=True)
ap.add_argument("--homoplasmic-min", type=float, default=0.95)
ap.add_argument("--heteroplasmic-min", type=float, default=0.05)
ap.add_argument("--low-level-min", type=float, default=0.01)
ap.add_argument("--min-depth", type=int, default=100)
ap.add_argument("--gene-start", type=int, default=648)   # MT-RNR1, rCRS
ap.add_argument("--gene-end", type=int, default=1601)
a = ap.parse_args()

RISK_ORDER = {"increased": 0, "uncertain": 1, "normal": 2}
alleles = [r for r in csv.DictReader((l for l in open(a.alleles) if not l.startswith("#")), delimiter="\t")
           if r["name"] != "Reference"]


def fasta_base(pos):
    out = subprocess.run(["samtools", "faidx", a.fasta, f"{a.contig}:{pos}-{pos}"], check=True,
                         capture_output=True, text=True).stdout.splitlines()
    return out[1].upper()


for r in alleles:
    if r["ref"] and fasta_base(int(r["pos"])) != r["ref"]:
        raise SystemExit(f"CPIC allele {r['name']}: REF {r['ref']} != FASTA base at {a.contig}:{r['pos']}")


def state(af):
    if af >= a.homoplasmic_min:
        return "homoplasmic"
    if af >= a.heteroplasmic_min:
        return "heteroplasmic"
    if af >= a.low_level_min:
        return "low-level (not called)"
    return "not detected"


# per-position allele depths; indel records kept separately
snv, indel, depth = {}, [], {}
with gzip.open(a.pileup, "rt") as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        pos, ref, alts = int(f[1]), f[3], [x for x in f[4].split(",") if x != "<*>"]
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        ad = [int(x) for x in fmt.get("AD", "0").split(",") if x != "."]
        tot = sum(ad)
        is_indel = "INDEL" in f[7].split(";")
        rec = {"pos": pos, "ref": ref, "depth": tot,
               "alts": {alt: (ad[i + 1] if i + 1 < len(ad) else 0) for i, alt in enumerate(alts)}}
        if is_indel:
            indel.append(rec)
        else:
            snv[pos] = rec
            depth[pos] = tot

gene_depths = sorted(depth.get(p, 0) for p in range(a.gene_start, a.gene_end + 1))
low_cov = [p for p in range(a.gene_start, a.gene_end + 1) if depth.get(p, 0) < a.min_depth]

# every non-reference observation in the gene
observed = []
for rec in list(snv.values()) + indel:
    for alt, n in rec["alts"].items():
        af = n / rec["depth"] if rec["depth"] else 0.0
        if af >= a.low_level_min:
            observed.append({"pos": rec["pos"], "ref": rec["ref"], "alt": alt, "alt_depth": n, "depth": rec["depth"],
                             "af": round(af, 4), "state": state(af), "type": "indel" if rec in indel else "snv"})
observed.sort(key=lambda x: (x["pos"], x["alt"]))

# CPIC alleles
cpic = []
for r in alleles:
    pos = int(r["pos"])
    risk = r["clinical_function"].split()[0]
    if r["alt"] in ("A", "C", "G", "T"):
        rec = snv.get(pos)
        n = (rec or {}).get("alts", {}).get(r["alt"], 0)
        dp = (rec or {}).get("depth", 0)
    else:  # C-tract deletions: any deletion record starting in 955-965
        hits = [(x["alt_depth"], x["depth"]) for x in observed if x["type"] == "indel" and 955 <= x["pos"] <= 965
                and len(x["alt"]) < len(x["ref"])]
        n, dp = (max(hits) if hits else (0, depth.get(pos, 0)))
    af = n / dp if dp else None
    cpic.append({"allele": r["name"], "clinical_function": r["clinical_function"], "risk": risk, "depth": dp,
                 "alt_depth": n, "af": round(af, 4) if af is not None else None,
                 "state": "no coverage" if not dp else state(af)})

called = [c for c in cpic if c["state"] in ("homoplasmic", "heteroplasmic")]
called.sort(key=lambda c: (RISK_ORDER.get(c["risk"], 9), -(c["af"] or 0)))
outside = called[0]["allele"] if called else "Reference"
flags = []
if any(c["state"] == "heteroplasmic" for c in called):
    flags.append("heteroplasmic CPIC allele: risk statements are made for carriers; no safe heteroplasmy threshold is defined")
if low_cov:
    flags.append(f"{len(low_cov)} positions in MT-RNR1 below {a.min_depth}x")
controls = [x for x in observed if (x["pos"], x["alt"]) in ((750, "G"), (1438, "G"))]

summary = {
    "sample": a.sample, "contig": a.contig, "gene_region": f"{a.contig}:{a.gene_start}-{a.gene_end}",
    "depth_median": gene_depths[len(gene_depths) // 2] if gene_depths else 0,
    "depth_min": gene_depths[0] if gene_depths else 0,
    "positions_below_min_depth": len(low_cov),
    "thresholds": {"homoplasmic_min": a.homoplasmic_min, "heteroplasmic_min": a.heteroplasmic_min,
                   "low_level_min": a.low_level_min, "min_depth": a.min_depth},
    "cpic_alleles_tested": len(cpic),
    "cpic_alleles_called": called,
    "cpic_low_level": [c for c in cpic if c["state"].startswith("low-level")],
    "pharmcat_outside_call": outside,
    "phenotype": (called[0]["clinical_function"] if called else "normal risk of aminoglycoside-induced hearing loss"),
    "flags": flags,
    "positive_controls_observed": controls,
    "non_reference_observations": observed,
    "cpic_alleles": cpic,
}
json.dump(summary, open(a.out_json, "w"), indent=2)
open(a.out_outside, "w").write(f"MT-RNR1\t{outside}\n")
print(json.dumps({k: summary[k] for k in ("depth_median", "cpic_alleles_called", "pharmcat_outside_call", "flags",
                                          "positive_controls_observed")}, indent=1))
