#!/usr/bin/env python3
"""Characterise a premature termination codon and the protein it would yield if read through.

A nonsense allele is a readthrough target only if three things hold, and this derives all three from
the pipeline's own outputs rather than asserting them:

  1. WHICH STOP CODON. Taken from VEP's codon field, not inferred from the amino-acid change. UGA is
     the most readthrough-permissive of the three stop codons and the one ataluren was characterised
     against; UAA is the least.
  2. WHAT CONTEXT. The base immediately 3' of the stop (+4) exerts the strongest single-nucleotide
     influence on readthrough efficiency, in the order C > U > G > A. It is read from the reference
     at the derived codon coordinates, so a disagreement with the transcript annotation would show up
     as a codon mismatch rather than pass silently.
  3. WHAT PROTEIN. Readthrough does not restore the reference residue. Near-cognate decoding inserts
     a small set of amino acids determined by the stop codon (at UGA: Trp, Arg, Cys). Each resulting
     full-length protein is scored with AlphaMissense, so "readthrough might work" becomes a specific
     prediction about specific proteins. All 19 substitutions at the residue are reported as well,
     because a residue that tolerates everything is a different claim from one that tolerates three
     things, and the regional means say whether the residue is unusual for its part of the protein.

AlphaMissense's published SNV file cannot answer point 3: the near-cognate products are reachable by
one nucleotide change from the mutant stop codon but not from the reference codon, so they are absent
from a reference-anchored SNV list. This reads the all-amino-acid-substitution file instead, streamed
and filtered to one accession so nothing large is written.

Nothing here is specific to this gene or this disease; --gene, --uniprot and the stop codon table are
the only inputs that change.

Usage: readthrough_analysis.py --vcf panel.vep.vcf.gz --fasta ref.fa --nmd nmd_ptc.json
                               --domains domain_map.json --gene BUB1B --uniprot O60566
                               --am-url <AlphaMissense_aa_substitutions.tsv.gz> --prefix out
"""
import argparse
import gzip
import json
import re
import subprocess
import sys
import urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("--vcf", required=True)
ap.add_argument("--fasta", required=True)
ap.add_argument("--nmd", required=True)
ap.add_argument("--domains", required=True)
ap.add_argument("--gene", required=True)
ap.add_argument("--uniprot", required=True)
ap.add_argument("--am-url", required=True)
ap.add_argument("--am-file", help="local AlphaMissense substitutions file instead of --am-url (tests)")
ap.add_argument("--clinvar", help="ClinVar VCF: how many pathogenic alleles in this gene are truncating")
ap.add_argument("--prefix", required=True)
a = ap.parse_args()

# Near-cognate decoding products per stop codon, and the relative permissiveness of each stop.
# Sources are carried through to the report so the claim can be checked, not taken on trust.
STOP_TABLE = {
    "TGA": {"rna": "UGA", "products": ["W", "R", "C"], "permissiveness": "most permissive"},
    "TAG": {"rna": "UAG", "products": ["Q", "Y", "K"], "permissiveness": "intermediate"},
    "TAA": {"rna": "UAA", "products": ["Q", "Y", "K"], "permissiveness": "least permissive"},
}
PLUS4_ORDER = ["C", "T", "G", "A"]   # C > U > G > A, strongest first
AA3 = {"W": "Trp", "R": "Arg", "C": "Cys", "Q": "Gln", "Y": "Tyr", "K": "Lys"}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    return s.translate(COMP)[::-1]


def faidx(region):
    out = subprocess.run(["samtools", "faidx", a.fasta, region], check=True, capture_output=True, text=True).stdout
    return "".join(l.strip() for l in out.splitlines() if not l.startswith(">")).upper()


# ---------------------------------------------------------------- the PTC, from VEP's own fields
def csq_format(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                return re.search(r"Format: ([^\"]+)", line).group(1).split("|")
            if not line.startswith("#"):
                break
    raise SystemExit(f"{path}: no CSQ header")


FIELDS = csq_format(a.vcf)
nmd = json.load(open(a.nmd))
ptcs = [p for p in nmd.get("ptcs") or [] if p["gene"] == a.gene]
if not ptcs:
    json.dump({"gene": a.gene, "status": "no PTC in this gene", "ptcs": []}, open(a.prefix + ".json", "w"), indent=2)
    print(json.dumps({"status": "no PTC"}))
    sys.exit(0)

op = gzip.open if a.vcf.endswith(".gz") else open
records = []
with op(a.vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        m = re.search(r"(?:^|;)CSQ=([^;\t]+)", f[7])
        if not m:
            continue
        for block in m.group(1).split(","):
            c = dict(zip(FIELDS, block.split("|")))
            if c.get("SYMBOL") != a.gene or "stop_gained" not in (c.get("Consequence") or ""):
                continue
            if "NMD_transcript_variant" in (c.get("Consequence") or ""):
                continue   # the NMD-biotype copies of the same call
            records.append({"chrom": f[0], "pos": int(f[1]), "ref": f[3], "alt": f[4], "csq": c})

out_ptcs = []
for p in ptcs:
    rec = next((r for r in records if r["csq"].get("HGVSc", "").endswith(p["hgvsc"])), records[0] if records else None)
    if rec is None:
        continue
    c = rec["csq"]
    codons = c.get("Codons") or ""                      # e.g. tTa/tGa: the changed base is upper case
    ref_codon, _, alt_codon = codons.partition("/")
    strand = int(c.get("STRAND") or 1)
    idx = next((i for i, ch in enumerate(ref_codon) if ch.isupper()), None)
    stop = alt_codon.upper()
    entry = {"gene": a.gene, "transcript": c.get("Feature"), "hgvsc": p["hgvsc"], "hgvsp": p["hgvsp"],
             "residue": int(re.sub(r"\D", "", (c.get("Protein_position") or "0").split("-")[0]) or 0),
             "strand": strand, "codon_change": codons, "stop_codon_dna": stop,
             "stop_codon": STOP_TABLE.get(stop, {}).get("rna", stop),
             "stop_permissiveness": STOP_TABLE.get(stop, {}).get("permissiveness", "unknown"),
             "nmd_verdict": p["verdict"], "nmd_distance_to_last_junction_nt": p["ptc_to_last_junction_nt"],
             "ptc_exon": p["ptc_exon"]}

    # genomic span of the codon, derived from which base VEP says changed; then the +4 base
    if idx is not None and strand in (1, -1):
        start = rec["pos"] - idx if strand == 1 else rec["pos"] - (2 - idx)
        obs = faidx(f"{rec['chrom']}:{start}-{start + 2}")
        obs_codon = obs if strand == 1 else rc(obs)
        entry["codon_region"] = f"{rec['chrom']}:{start}-{start + 2}"
        entry["reference_codon_confirmed"] = (obs_codon == ref_codon.upper())
        plus4_pos = start + 3 if strand == 1 else start - 1
        base = faidx(f"{rec['chrom']}:{plus4_pos}-{plus4_pos}")
        entry["plus4_base"] = base if strand == 1 else rc(base)
        entry["plus4_base_rna"] = entry["plus4_base"].replace("T", "U")
        entry["plus4_rank"] = (PLUS4_ORDER.index(entry["plus4_base"]) + 1) if entry["plus4_base"] in PLUS4_ORDER else None
        entry["plus4_of"] = len(PLUS4_ORDER)
        entry["plus4_order"] = "C > U > G > A"
        entry["stop_context"] = f"{entry['stop_codon']}-{entry['plus4_base_rna']}"
    out_ptcs.append(entry)

# ---------------------------------------------------------------- the readthrough products
dom = json.load(open(a.domains))
domains = dom.get("domains") or []


def region_of(pos):
    for d in domains:
        if d["start"] <= pos <= d["end"]:
            return d["domain"]
    prev = max((d for d in domains if d["end"] < pos), key=lambda d: d["end"], default=None)
    nxt = min((d for d in domains if d["start"] > pos), key=lambda d: d["start"], default=None)
    if prev and nxt:
        return f"inter-domain linker between {prev['domain']} (ends {prev['end']}) and {nxt['domain']} (begins {nxt['start']})"
    return "outside the annotated domains"


def am_rows():
    """All AlphaMissense substitutions for one accession, streamed so the 1.2 GB file is never stored."""
    if a.am_file:
        op2 = gzip.open if a.am_file.endswith(".gz") else open
        fh = op2(a.am_file, "rt")
    else:
        fh = gzip.open(urllib.request.urlopen(a.am_url), "rt")
    with fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] == a.uniprot:
                yield f[1], float(f[2]), f[3]


by_res = {}
for sub, score, cls in am_rows():
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", sub)
    if m:
        by_res.setdefault(int(m.group(2)), {})[m.group(3)] = {"score": round(score, 4), "class": cls}

if not by_res:
    raise SystemExit(f"AlphaMissense: no rows for {a.uniprot}")


def region_mean(lo, hi):
    vals = [v["score"] for r, subs in by_res.items() if lo <= r <= hi for v in subs.values()]
    return round(sum(vals) / len(vals), 4) if vals else None


for e in out_ptcs:
    res = e["residue"]
    subs = by_res.get(res) or {}
    e["residue_region"] = region_of(res)
    prods = STOP_TABLE.get(e["stop_codon_dna"], {}).get("products", [])
    ref3 = (re.match(r"p\.([A-Z][a-z]{2})\d+", e["hgvsp"]) or [None, "Xaa"])[1]
    e["readthrough_products"] = [
        {"aa1": p, "aa3": AA3.get(p, p), "label": f"p.{ref3}{res}{AA3.get(p, p)}",
         "alphamissense": subs.get(p, {}).get("score"), "class": subs.get(p, {}).get("class")}
        for p in prods]
    e["residue_tolerance"] = {
        "n_substitutions_scored": len(subs),
        "n_benign": sum(1 for v in subs.values() if v["class"] == "benign"),
        "n_ambiguous": sum(1 for v in subs.values() if v["class"] == "ambiguous"),
        "n_pathogenic": sum(1 for v in subs.values() if v["class"] == "pathogenic"),
        "mean_score": round(sum(v["score"] for v in subs.values()) / len(subs), 4) if subs else None,
        "max_score": max((v["score"] for v in subs.values()), default=None),
        "all": {k: v for k, v in sorted(subs.items(), key=lambda kv: kv[1]["score"])},
    }
    e["regional_means"] = {d["domain"]: region_mean(d["start"], d["end"]) for d in domains}
    e["regional_means"]["whole protein"] = region_mean(1, dom.get("length") or max(by_res))
    prod_scores = [p["alphamissense"] for p in e["readthrough_products"] if p["alphamissense"] is not None]
    prod_classes = {p["class"] for p in e["readthrough_products"] if p["class"]}
    e["verdict"] = ("all predicted tolerated" if prod_classes and prod_classes <= {"benign", "likely_benign"}
                    else ("mixed" if prod_scores else "not scored"))

# ---------------------------------------------------------------- how far the axis reaches
# A readthrough strategy applies to whoever carries a nonsense allele, so the share of this gene's
# pathogenic alleles that are truncating says how much of the disease the axis could address.
clinvar = None
if a.clinvar:
    counts = {"pathogenic_total": 0, "nonsense": 0, "frameshift": 0, "splice": 0, "other": 0}
    op3 = gzip.open if a.clinvar.endswith(".gz") else open
    with op3(a.clinvar, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            info = line.split("\t")[7] if line.count("\t") >= 7 else ""
            gi = re.search(r"(?:^|;)GENEINFO=([^;]+)", info)
            if not gi or a.gene not in re.split(r"[|:]", gi.group(1)):
                continue
            sig = (re.search(r"(?:^|;)CLNSIG=([^;]+)", info) or [None, ""])[1]
            if sig.split("/")[0] not in ("Pathogenic", "Likely_pathogenic") and sig != "Pathogenic/Likely_pathogenic":
                continue
            counts["pathogenic_total"] += 1
            mc = (re.search(r"(?:^|;)MC=([^;]+)", info) or [None, ""])[1]
            terms = {t.split("|")[-1] for t in mc.split(",") if t}
            if "nonsense" in terms:
                counts["nonsense"] += 1
            elif "frameshift_variant" in terms:
                counts["frameshift"] += 1
            elif any("splice" in t for t in terms):
                counts["splice"] += 1
            else:
                counts["other"] += 1
    counts["truncating"] = counts["nonsense"] + counts["frameshift"]
    counts["truncating_pct"] = (round(100 * counts["truncating"] / counts["pathogenic_total"], 1)
                                if counts["pathogenic_total"] else None)
    clinvar = counts

summary = {"gene": a.gene, "uniprot": a.uniprot, "alphamissense_source": a.am_file or a.am_url,
           "clinvar_pathogenic_alleles": clinvar,
           "plus4_order_note": "the +4 base influences readthrough efficiency in the order C > U > G > A",
           "ptcs": out_ptcs}
json.dump(summary, open(a.prefix + ".json", "w"), indent=2)

with open(a.prefix + ".tsv", "w") as fh:
    fh.write("gene\thgvsp\tstop_codon\tstop_permissiveness\tplus4_base\tplus4_rank\tstop_context\t"
             "nmd_verdict\tresidue\tresidue_region\tproduct\talphamissense\tclass\n")
    for e in out_ptcs:
        for p in e["readthrough_products"]:
            fh.write(f"{e['gene']}\t{e['hgvsp']}\t{e['stop_codon']}\t{e['stop_permissiveness']}\t"
                     f"{e.get('plus4_base','')}\t{e.get('plus4_rank','')}\t{e.get('stop_context','')}\t"
                     f"{e['nmd_verdict']}\t{e['residue']}\t{e['residue_region']}\t"
                     f"{p['aa3']}\t{p['alphamissense']}\t{p['class']}\n")

print(json.dumps({"ptcs": len(out_ptcs),
                  "stop_codons": [e["stop_codon"] for e in out_ptcs],
                  "plus4": [e.get("plus4_base") for e in out_ptcs],
                  "verdicts": [e["verdict"] for e in out_ptcs]}))
