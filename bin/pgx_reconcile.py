#!/usr/bin/env python3
"""Reconcile PharmCAT and Aldy diplotypes per gene. Neither caller is silently preferred.

Inputs : PharmCAT <sample>.phenotype.json, aldy_calls.json (aldy_parse.py), off_definition.tsv
Outputs: pgx_calls.tsv / pgx_calls.json with, per gene:
  diplotype(s) from each caller, PharmCAT phenotype and activity score, call source
  (PharmCAT matcher vs outside call), agreement flag, undocumented-variation flag, notes.

Agreement flags
  agree            same allele multiset after normalisation (Aldy suballeles *4.022 -> *4;
                   PharmCAT xN duplications expanded)
  disagree         both called with star nomenclature and the allele multisets differ
  not_comparable   the nomenclatures cannot be mapped without assumptions (e.g. Aldy VKORC1
                   H-haplotypes vs a PharmCAT rsID call; Aldy base-named IFNL3 alleles); both calls
                   are reported side by side
Non-star names are compared after explicit normalisation, recorded in 'comparison_basis':
reference synonyms (Reference/*ref/*Ref/Aldy's truncated '*ference'/'rsN reference (X)'),
c.-notation with punctuation removed, PharmCAT '(heterozygous)/(homozygous)' expanded,
PharmCAT's own star alias used where given (DPYD 'c.2194G>A (*6)'), and a hemizygous chrX call
compared against Aldy's diploid output.
  outside_call     PharmCAT's call for this gene *is* the Aldy call (e.g. CYP2D6); agreement
                   cannot be assessed independently
  pharmcat_only / aldy_only / no_call
PharmCAT ambiguity (several source diplotypes) is reported, and 'agree' requires Aldy to match
one of them.

Usage: pgx_reconcile.py <phenotype.json> <aldy_calls.json> <off_definition.tsv> <genes> <out_prefix>
"""
import csv
import json
import re
import sys

pheno_path, aldy_path, offdef_path, genes_arg, prefix = sys.argv[1:6]
genes = [g.strip() for g in genes_arg.split(",") if g.strip()]
reports = json.load(open(pheno_path))["geneReports"]
aldy = json.load(open(aldy_path))
offdef = {}
with open(offdef_path) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        offdef.setdefault(row["chrom"], []).append(row)

STAR = re.compile(r"^\*\d+[A-Z]?(x\d+)?$")


def expand(allele):
    m = re.match(r"^(\*[^x]+)x(\d+)$", allele)
    return [m.group(1)] * int(m.group(2)) if m else [allele]


def norm_aldy(dip):
    alleles = []
    for hap in dip.split("/"):
        for a in hap.split("+"):
            a = re.sub(r"\.\d+$", "", a.strip())      # suballele -> core allele
            alleles.append(a)
    return sorted(alleles)


def norm_pharmcat(label):
    # PharmCAT labels: '*1/*17', '*4/*35x2', 'c.6178G>T (heterozygous)', 'Reference/Reference'
    alleles = []
    for a in label.split("/"):
        alleles.extend(expand(a.strip()))
    return sorted(alleles)


REF_SYNONYMS = {"reference", "ref", "ference"}          # 'ference' = Aldy's truncation of 'Reference'
UNMAPPABLE = [re.compile(r"^[acgt]+$"),                  # Aldy base-named alleles (e.g. IFNL3 *C)
              re.compile(r"^h\d+$")]                      # VKORC1 H-haplotypes vs PharmCAT rsID calls
ZYGOSITY = re.compile(r"^(.*?)\s*\((heterozygous|homozygous)\)$", re.I)
ALIAS = re.compile(r"^(.*?)\s*\((\*\d+[A-Z]?)\)$")


def sem_token(x):
    low = re.sub(r"^\*", "", x.strip()).lower()
    if low in REF_SYNONYMS or re.match(r"^rs\d+ reference \([acgt]+\)$", low):
        return "reference"
    low = re.sub(r"\s*\(reference\)$", "", low)          # G6PD 'B (reference)' -> 'b'
    tok = re.sub(r"[^a-z0-9]", "", re.sub(r"^c\.", "", low))
    return None if any(p.match(tok) for p in UNMAPPABLE) else tok


def pharmcat_alias_star(label):
    """'c.2194G>A (*6) (heterozygous)' -> ['*1', '*6'] using PharmCAT's own star alias."""
    m = ZYGOSITY.match(label.strip())
    if not m:
        return None
    a = ALIAS.match(m.group(1))
    if not a:
        return None
    return sorted(["*1", a.group(2)] if m.group(2).lower() == "heterozygous" else [a.group(2)] * 2)


def semantic_pharmcat(label):
    m = ZYGOSITY.match(label.strip())
    if m:
        t = sem_token(m.group(1))
        toks = ["reference", t] if m.group(2).lower() == "heterozygous" else [t, t]
    else:
        toks = [sem_token(p) for p in label.split("/")]
    return None if any(t is None for t in toks) else sorted(toks)


def semantic_aldy(dip):
    toks = [sem_token(p) for h in dip.split("/") for p in h.split("+")]
    return None if any(t is None for t in toks) else sorted(toks)


def is_star(alleles):
    return all(STAR.match(a) or re.match(r"^\*\d+[A-Z]?$", a) for a in alleles)


rows = []
for gene in genes:
    r = reports.get(gene)
    a = aldy.get(gene)
    pc_labels, pc_pheno, pc_as, source, undoc, msgs = [], [], [], None, None, []
    if r:
        source = r.get("callSource")
        undoc = r.get("hasUndocumentedVariations")
        for d in r.get("sourceDiplotypes") or []:
            if d.get("label") and d["label"] != "Unknown/Unknown":
                pc_labels.append(d["label"])
                pc_pheno.extend(d.get("phenotypes") or [])
                if d.get("activityScore") is not None:
                    pc_as.append(str(d["activityScore"]))
        msgs = [m.get("message", "") for m in (r.get("messages") or []) if m.get("exception_type") != "note"]
    aldy_dip = a["diplotype"] if a and a["status"] == "called" else None

    basis = None
    if not pc_labels and not aldy_dip:
        flag = "no_call"
    elif source == "OUTSIDE":
        flag = "outside_call"
    elif pc_labels and not aldy_dip:
        flag = "pharmcat_only"
    elif aldy_dip and not pc_labels:
        flag = "aldy_only"
    else:
        pc_sets = [pharmcat_alias_star(l) or norm_pharmcat(l) for l in pc_labels]
        al_set = norm_aldy(aldy_dip)
        if all(is_star(s) for s in pc_sets) and is_star(al_set):
            flag = "agree" if al_set in pc_sets else "disagree"
            basis = "star alleles" + (" (PharmCAT star alias)" if any(pharmcat_alias_star(l) for l in pc_labels) else "")
        else:
            # explicit semantic normalisation for non-star nomenclature (see module docstring)
            pc_sem = [semantic_pharmcat(l) for l in pc_labels]
            al_sem = semantic_aldy(aldy_dip)
            if any(s is None for s in pc_sem) or al_sem is None:
                flag, basis = "not_comparable", "nomenclatures cannot be mapped without assumptions"
            elif al_sem in pc_sem:
                flag, basis = "agree", "normalised names"
            elif any(len(s) == 1 and len(set(al_sem)) == 1 and s[0] == al_sem[0] for s in pc_sem):
                flag, basis = "agree", "normalised names; hemizygous call vs Aldy diploid output (Aldy does not model chrX hemizygosity)"
            else:
                flag, basis = "disagree", "normalised names"

    rows.append({
        "gene": gene,
        "pharmcat_diplotype": " | ".join(pc_labels) or None,
        "pharmcat_ambiguous": len(pc_labels) > 1,
        "pharmcat_call_source": source,
        "aldy_diplotype": aldy_dip,
        "aldy_status": a["status"] if a else "not_run",
        "aldy_gene_copies": a.get("gene_copies_total") if a else None,
        "aldy_all_solutions": " | ".join(a["major_solutions"]) if a else None,
        "phenotype": " | ".join(dict.fromkeys(pc_pheno)) or None,
        "activity_score": " | ".join(dict.fromkeys(pc_as)) or None,
        "caller_agreement": flag,
        "comparison_basis": basis,
        "undocumented_variation": bool(undoc),
        "pharmcat_messages": " || ".join(msgs) or None,
    })

fields = list(rows[0].keys())
with open(prefix + ".tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for row in rows:
        w.writerow({k: ("" if v is None else v) for k, v in row.items()})
summary = {
    "genes": len(rows),
    "agreement_counts": {f: sum(1 for r in rows if r["caller_agreement"] == f)
                         for f in sorted({r["caller_agreement"] for r in rows})},
    "disagreements": [r["gene"] for r in rows if r["caller_agreement"] == "disagree"],
    "undocumented_variation_genes": [r["gene"] for r in rows if r["undocumented_variation"]],
    "off_definition_sites": sum(len(v) for v in offdef.values()),
}
json.dump({"summary": summary, "calls": rows}, open(prefix + ".json", "w"), indent=2)
print(json.dumps(summary))
