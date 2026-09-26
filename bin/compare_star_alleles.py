#!/usr/bin/env python3
"""Quick check: do Aldy and a second star-allele caller (Stargazer or PyPGx) report the same
diplotype per gene?

Inputs
  --aldy     aldy_calls.json from aldy_parse.py (or one or more raw .aldy files)
  --other    Stargazer genotype-calls.tsv (columns: name/gene, hap1_main, hap2_main, ...)
             or PyPGx data.tsv / results TSV (column 'Genotype', one file per gene, or
             a 'Gene' column); format auto-detected
  --gene     gene name when the other-caller file has no gene column (PyPGx per-gene output)

Comparison: alleles are normalised to core star alleles (Aldy '*4.022' -> '*4'; Stargazer
'*4.001' -> '*4'; PyPGx '*4' unchanged), copy-number notation is expanded ('*2x2' -> '*2','*2';
Aldy '*35+*35' -> '*35','*35'), and the allele multiset is compared (order-independent).
Suballele letters are kept (*3A != *3B).

Non-star nomenclature is normalised explicitly before comparing: reference synonyms
(Reference, *ref, *Ref, Aldy's truncated *ference) -> REF; a variant name carrying a star alias
('c.2194G>A (*6)') -> the alias, with Reference -> *1 in that case; c.-notation punctuation removed
('c.6178G>T' == Aldy '*6178GT'); PyPGx 'MALE' hemizygous marker dropped and compared on the single
allele. Alleles that cannot be mapped without assumptions (Aldy base-named '*C', VKORC1 '*H6')
give not_comparable. The result per gene is one of:
  same | same_after_normalisation | different | not_comparable | only_aldy | only_other | not_called

Usage: compare_star_alleles.py --aldy aldy_calls.json --other genotype-calls.tsv [--gene CYP2D6] [--out cmp.tsv]
Exit status is 0 even when calls differ; this is a report, not a gate.
"""
import argparse
import csv
import json
import re
import sys


def core(allele):
    a = allele.strip()
    a = re.sub(r"\.\d+$", "", a)                     # suballele numbering, not letter suffixes
    return a


def expand(token):
    token = core(token)
    m = re.match(r"^(\*[^x]+)x(\d+|N)$", token, re.I)
    if m:
        n = 2 if m.group(2).upper() == "N" else int(m.group(2))
        return [m.group(1)] * n
    return [token]


def alleles_of(diplotype):
    if not diplotype or diplotype.strip() in ("", ".", "Indeterminate", "NA", "None"):
        return None
    out = []
    for hap in diplotype.split("/"):
        for part in hap.split("+"):
            if part.strip():
                out.extend(expand(part))
    return sorted(out, key=lambda x: (re.sub(r"\D", "", x).zfill(4), x))


def read_aldy(paths):
    calls = {}
    if len(paths) == 1 and paths[0].endswith(".json"):
        for gene, c in json.load(open(paths[0])).items():
            calls[gene.upper()] = c.get("diplotype") if c.get("status") == "called" else None
        return calls
    for path in paths:
        majors, gene = set(), None
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            gene = f[1].upper()
            majors.add(f[3])
        if gene:
            calls[gene] = majors.pop() if len(majors) == 1 else None
    return calls


def read_other(path, gene_arg):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    if not rows:
        return {}, "empty"
    cols = {c.lower(): c for c in rows[0]}
    calls = {}
    if "hap1_main" in cols:                               # Stargazer genotype-calls.tsv
        gcol = cols.get("gene") or cols.get("name")
        for r in rows:
            g = (gene_arg or r.get(gcol, "")).upper()
            h1, h2 = r[cols["hap1_main"]], r[cols["hap2_main"]]
            calls[g] = f"{h1}/{h2}" if h1 and h2 and h1 != "." else None
        return calls, "stargazer"
    if "genotype" in cols:                                # PyPGx data.tsv / results
        gcol = cols.get("gene")
        for r in rows:
            g = (r[gcol] if gcol else gene_arg or "").upper()
            if not g:
                sys.exit("PyPGx output has no gene column; pass --gene")
            calls[g] = r[cols["genotype"]] or None
        return calls, "pypgx"
    sys.exit(f"Unrecognised format for {path}: columns {list(rows[0])}")


REF_SYN = {"reference", "ref", "ference"}
UNMAPPABLE = re.compile(r"^(?:[acgt]+|h\d+)$")


def normalise(diplotype):
    """Semantic allele multiset, or None if any allele cannot be mapped without assumptions."""
    if not diplotype:
        return None
    toks = [t for t in re.split(r"[/+]", diplotype) if t.strip()]
    toks = [t for t in toks if t.strip().upper() != "MALE"]           # PyPGx hemizygous marker
    aliases = [re.search(r"\((\*\d+[A-Z]?)\)", t) for t in toks]
    use_star = any(aliases)
    out = []
    for t, al in zip(toks, aliases):
        low = re.sub(r"^\*", "", core(t.strip())).lower()
        low = re.sub(r"\s*\(reference\)$", "", low)
        if al:
            out.append(al.group(1))
        elif low in REF_SYN:
            out.append("*1" if use_star else "REF")
        elif re.match(r"^\d+[a-z]?$", low):
            out.append("*" + low.upper())
        else:
            tok = re.sub(r"[^a-z0-9]", "", re.sub(r"^c\.", "", low))
            if UNMAPPABLE.match(tok):
                return None
            out.append(tok)
    if len(set(out)) == 1 and len(out) == 2 and any(t.strip().upper() == "MALE" for t in re.split(r"[/+]", diplotype)):
        out = out[:1]
    return sorted(out)


ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--aldy", nargs="+", required=True)
ap.add_argument("--other", nargs="+", required=True)
ap.add_argument("--gene")
ap.add_argument("--out")
a = ap.parse_args()

aldy = read_aldy(a.aldy)
other, fmt = {}, None
for p in a.other:
    c, fmt = read_other(p, a.gene)
    other.update(c)

results = []
for gene in sorted(set(aldy) | set(other)):
    x, y = alleles_of(aldy.get(gene)), alleles_of(other.get(gene))
    if x is None and y is None:
        status = "not_called"
    elif y is None:
        status = "only_aldy"
    elif x is None:
        status = "only_other"
    elif x == y:
        status = "same"
    else:
        nx, ny = normalise(aldy.get(gene)), normalise(other.get(gene))
        if nx is None or ny is None:
            status = "not_comparable"
        elif nx == ny or (len(ny) == 1 and len(set(nx)) == 1 and nx[0] == ny[0]) or (len(nx) == 1 and len(set(ny)) == 1 and ny[0] == nx[0]):
            status = "same_after_normalisation"
        else:
            status = "different"
    results.append({"gene": gene, "aldy": aldy.get(gene) or "", f"{fmt}": other.get(gene) or "",
                    "aldy_core": "/".join(x or []), "other_core": "/".join(y or []), "status": status})

fields = list(results[0].keys()) if results else ["gene", "status"]
out = open(a.out, "w") if a.out else sys.stdout
w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
w.writeheader()
w.writerows(results)
counts = {s: sum(1 for r in results if r["status"] == s)
          for s in ("same", "same_after_normalisation", "different", "not_comparable", "only_aldy", "only_other", "not_called")}
counts["different_genes"] = [r["gene"] for r in results if r["status"] == "different"]
print(json.dumps({"other_caller": fmt, **counts}), file=sys.stderr)
