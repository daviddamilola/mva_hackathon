#!/usr/bin/env python3
"""NMD sensitivity of premature termination codons (PTCs) on MANE Select transcripts.

For every stop_gained variant in the VEP-annotated VCF, on the MANE Select transcript:

  exon structure   from the SpliceAI-format MANE annotation (EXON_START 0-based, EXON_END)
  consistency      the cDNA position recomputed from the exon structure must equal VEP's
                   cDNA_position, else the variant is reported as 'structure mismatch' and not called
  PTC              the premature stop codon: CDS codon (Protein_position) -> cDNA positions
  last junction    the last exon-exon junction = cDNA position of the last base of the penultimate exon
  distance         last junction - last nucleotide of the PTC (nt, spliced mRNA)

Rules (canonical EJC model; Lindeboom et al. 2016 Nat Genet; Nagy & Maquat 1998):
  PTC in the last exon                              -> escapes NMD
  PTC <= `rule_nt` nt upstream of the last junction -> escapes NMD (50-55 nt rule)
  PTC start-proximal (CDS position < `start_nt`)    -> may escape (translation re-initiation)
  PTC in an exon longer than `long_exon_nt`         -> may escape (long-exon rule)
  otherwise                                         -> predicted NMD-sensitive

Usage: nmd_ptc.py --vcf panel.vep.vcf.gz --annotation mane.txt --out-tsv nmd.tsv --out-json nmd.json
"""
import argparse
import csv
import gzip
import json

ap = argparse.ArgumentParser()
for x in ("vcf", "annotation", "out-tsv", "out-json"):
    ap.add_argument(f"--{x}", required=True)
ap.add_argument("--rule-nt", type=int, default=55)
ap.add_argument("--start-nt", type=int, default=150)
ap.add_argument("--long-exon-nt", type=int, default=407)
a = ap.parse_args()

tx = {}
for l in open(a.annotation):
    if l.startswith("#") or not l.strip():
        continue
    name, chrom, strand, _, _, es, ee = l.rstrip("\n").split("\t")[:7]
    starts = [int(x) for x in es.strip(",").split(",")]
    ends = [int(x) for x in ee.strip(",").split(",")]
    ex = sorted(zip(starts, ends), reverse=(strand == "-"))  # transcript order
    tx[name] = {"chrom": chrom, "strand": strand, "exons": ex}


def cdna_of(t, pos):
    """1-based cDNA position of 1-based genomic `pos`, and the exon index (0-based), or (None, None)."""
    acc = 0
    for i, (s, e) in enumerate(t["exons"]):
        if s < pos <= e:
            return acc + ((pos - s) if t["strand"] == "+" else (e - pos + 1)), i
        acc += e - s
    return None, None


fields = None
rows = []
with gzip.open(a.vcf, "rt") as fh:
    for l in fh:
        if l.startswith("##INFO=<ID=CSQ"):
            fields = l.split("Format: ")[1].rstrip('">\n').split("|")
            continue
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        for c in (info.get("CSQ") or "").split(","):
            d = dict(zip(fields, c.split("|")))
            if not d.get("MANE_SELECT") or "stop_gained" not in d["Consequence"].split("&"):
                continue
            t = tx.get(d["SYMBOL"])
            if not t:
                continue
            pos = int(f[1])
            cdna, ei = cdna_of(t, pos)
            vep_cdna = int(d["cDNA_position"].split("-")[0])
            lens = [e - s for s, e in t["exons"]]
            n = len(lens)
            last_junction = sum(lens[:-1])
            aa = int(d["Protein_position"].split("-")[0])
            cds = int(d["CDS_position"].split("-")[0])
            utr5 = vep_cdna - cds
            ptc_first, ptc_last = utr5 + 3 * aa - 2, utr5 + 3 * aa
            dist = last_junction - ptc_last
            if cdna != vep_cdna:
                verdict = "structure mismatch (recomputed cDNA != VEP); not called"
            elif ei == n - 1:
                verdict = "escapes NMD (PTC in last exon)"
            elif dist <= a.rule_nt:
                verdict = f"escapes NMD (PTC within {a.rule_nt} nt of the last exon-exon junction)"
            elif 3 * aa < a.start_nt:
                verdict = "may escape NMD (start-proximal; re-initiation possible)"
            elif lens[ei] > a.long_exon_nt:
                verdict = f"may escape NMD (PTC in an exon > {a.long_exon_nt} nt)"
            else:
                verdict = "predicted NMD-sensitive"
            rows.append({
                "gene": d["SYMBOL"], "transcript": d["Feature"], "mane_select": d["MANE_SELECT"],
                "hgvsc": d["HGVSc"].split(":")[-1], "hgvsp": d["HGVSp"].split(":")[-1],
                "ptc_exon": f"{ei + 1}/{n}" if ei is not None else "", "vep_exon": d["EXON"],
                "ptc_exon_length": lens[ei] if ei is not None else "",
                "cdna_variant": vep_cdna, "cdna_recomputed": cdna,
                "ptc_cdna": f"{ptc_first}-{ptc_last}", "last_junction_cdna": last_junction,
                "ptc_to_last_junction_nt": dist,
                "exons_downstream_of_ptc": (n - 1 - ei) if ei is not None else "",
                "verdict": verdict,
            })

cols = list(rows[0]) if rows else ["gene"]
with open(a.out_tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
json.dump({"rules": {"last_junction_nt": a.rule_nt, "start_proximal_nt": a.start_nt, "long_exon_nt": a.long_exon_nt},
           "ptcs": rows}, open(a.out_json, "w"), indent=2)
print(json.dumps(rows, indent=1))
