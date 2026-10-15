#!/usr/bin/env python3
"""Parse Aldy .aldy outputs into (a) a JSON summary of calls and (b) a PharmCAT outside-call file.

Aldy reports one or more solutions per gene; the 'Major' column holds the star-allele diplotype,
with '+' joining tandem copies on one haplotype (e.g. '*4/*35+*35'). Copy-number duplications of
the same allele are converted to PharmCAT notation ('*35+*35' -> '*35x2'); hybrid tandems
(e.g. '*36+*10') are kept as-is, which PharmCAT accepts for CYP2D6.

A gene is passed to PharmCAT as an outside call only if (1) it is listed in --outside-genes and
(2) Aldy returned exactly one distinct major solution. Ambiguous calls are reported, not passed.

Usage: aldy_parse.py --outside-genes CYP2D6 --json aldy_calls.json --outside outside_calls.tsv <gene.aldy>...
"""
import argparse
import json
import re
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--outside-genes", default="CYP2D6")
ap.add_argument("--json", required=True)
ap.add_argument("--outside", required=True)
ap.add_argument("aldy", nargs="+")
args = ap.parse_args()
outside_genes = {g.strip().upper() for g in args.outside_genes.split(",") if g.strip()}


def to_pharmcat_haplotype(hap):
    copies = [c for c in hap.split("+") if c]
    counts = Counter(copies)
    if len(counts) == 1 and len(copies) > 1:
        return f"{copies[0]}x{len(copies)}"
    return "+".join(copies)


def to_pharmcat_diplotype(major):
    return "/".join(to_pharmcat_haplotype(h) for h in major.split("/"))


calls = {}
for path in args.aldy:
    gene = None
    solutions = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            gene = f[1].upper()
            sol_id, major, minor = f[2], f[3], f[4]
            solutions.setdefault(sol_id, {"major": major, "minor": minor})
    if gene is None:
        # Aldy writes no records when it cannot call the gene
        gene = re.sub(r"\.aldy$", "", path.split("/")[-1]).upper()
    majors = sorted({s["major"] for s in solutions.values()})
    copy_number = None
    if len(majors) == 1:
        copy_number = sum(len([c for c in h.split("+") if c]) for h in majors[0].split("/"))
    calls[gene] = {
        "gene": gene,
        "n_solutions": len(solutions),
        "major_solutions": majors,
        "minor_solutions": sorted({s["minor"] for s in solutions.values()}),
        "status": "no_call" if not majors else ("ambiguous" if len(majors) > 1 else "called"),
        "diplotype": majors[0] if len(majors) == 1 else None,
        "pharmcat_diplotype": to_pharmcat_diplotype(majors[0]) if len(majors) == 1 else None,
        "gene_copies_total": copy_number,
        "passed_to_pharmcat": gene in outside_genes and len(majors) == 1,
    }

with open(args.json, "w") as out:
    json.dump(calls, out, indent=2, sort_keys=True)
with open(args.outside, "w") as out:
    for gene in sorted(calls):
        c = calls[gene]
        if c["passed_to_pharmcat"]:
            out.write(f"{gene}\t{c['pharmcat_diplotype']}\n")
print(json.dumps({g: c["status"] for g, c in calls.items()}))
