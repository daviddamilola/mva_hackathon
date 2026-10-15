#!/usr/bin/env python3
"""Stage 1b join: child genotypes (stage 1a + extra force-calls) x Open Targets PGx annotations.

JOIN RULES (each fails silently if done naively, so each is explicit here)
 1. Key type: if haplotypeId is populated, 'genotype' is a star allele or diplotype ('*1/*3A',
    '*12') and the gene comes from the haplotypeId prefix ('TPMT*1/*3A' -> TPMT; 'genotype'
    carries no prefix). Otherwise 'genotype' is a nucleotide genotype ('TT') keyed by
    variantRsId, with coordinates from genotypeId ('16_69711242_G_G,G').
 2. Star rows: a diplotype row ('*1/*3A') matches if its canonical form equals the called
    diplotype. A single-allele row ('*12') matches if the child carries that allele. Allele order
    is canonicalised (natural sort) before comparing. Suballeles are NOT collapsed (*3A != *3B != *3).
    CYP2D6 copy-number notation is expanded (*35x2 -> *35,*35).
 3. Nucleotide rows: the child's allele multiset at the position must equal the row's. The row's
    REF must match the genome; a reverse-strand row (complement matches) is flipped and counted.
    For a haploid call (chrX in XY), a homozygous row matches the single allele.
 4. Null target / targetFromSourceId: gene = haplotypeId prefix, then approvedSymbol, then
    targetFromSourceId, else 'UNKNOWN'. Such rows are kept and counted.
 5. isDirectTarget: rows where the gene is the drug's target are excluded from verdicts by default
    (this is a metaboliser panel), counted, and listed separately.
 6. Evidence: only levels in --assert-levels (1A,1B) contribute to verdicts. Level-3 rows are counted
    and reported as exploratory and carried no further. Other levels are counted.

VERDICT FLAG (per drug x gene)
  contraindicated         PharmCAT guideline annotation for the child's genotype with
                          alternateDrugAvailable = true
  dose-adjust             guideline annotation with dosingInformation = true
  no flag                 guideline annotation exists with neither flag, OR (no guideline) the
                          child is reference at every asserted key
  review (annotation only) no guideline, but the child carries a non-reference allele at an
                          asserted (1A/1B) OT key. The OT annotation text is reported and no
                          dosing claim is made.
  no data                 no asserted OT rows, or gene/positions not genotyped
Guideline flags come from PharmCAT's structured fields, not free-text keyword matching.
"""
import argparse
import csv
import gzip
import json
import os
import re
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--calls", required=True)
ap.add_argument("--pharmcat-report", required=True)
ap.add_argument("--vcf", nargs="+", required=True)
ap.add_argument("--drugs", nargs="+", required=True)
ap.add_argument("--targets", required=True)
ap.add_argument("--meta", required=True)
ap.add_argument("--assert-levels", default="1A,1B")
ap.add_argument("--include-direct-targets", action="store_true")
ap.add_argument("--prefix", required=True)
ap.add_argument("--candidates", help="candidates.json from resolve_candidates.py (axis, tier, anti-target verdict)")
ap.add_argument("--mt-rnr1", help="MT_RNR1_CALL json: per-allele observations on the mitochondrial contig")
ap.add_argument("--mt-alleles", help="CPIC MT-RNR1 allele table with the rsid column")
ap.add_argument("--drug-list-label", default="unrecorded")
ap.add_argument("--drug-list-example", action="store_true")
a = ap.parse_args()
ASSERT = set(a.assert_levels.split(","))
COMP = str.maketrans("ACGT", "TGCA")

# ---------------------------------------------------------------- child genotypes
pos_gt = {}
for path in a.vcf:
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]))
            if key in pos_gt:
                continue  # first VCF (stage 1a merged) wins
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fmt.get("GT", ".").replace("|", "/")
            alleles_all = [f[3]] + [x for x in f[4].split(",") if x != "."]
            idx = [x for x in gt.split("/") if x not in (".", "")]
            pos_gt[key] = {"ref": f[3], "gt": gt, "dp": fmt.get("DP"),
                           "alleles": [alleles_all[int(i)] for i in idx] if idx else None}

calls = {r["gene"]: r for r in json.load(open(a.calls))["calls"]}


def natkey(allele):
    m = re.match(r"^\*(\d+)(.*)$", allele)
    return (0, int(m.group(1)), m.group(2)) if m else (1, 0, allele)


def expand(allele):
    m = re.match(r"^(\*[^x]+)x(\d+)$", allele)
    return [m.group(1)] * int(m.group(2)) if m else [allele]


def canon_dip(label):
    # haplotype level: '*4/*35x2' -> ('*35x2','*4') sorted; also allele multiset
    haps = [h.strip() for h in label.split("/")]
    alleles = []
    for h in haps:
        for part in h.split("+"):
            alleles.extend(expand(re.sub(r"\.\d+$", "", part.strip())))
    return tuple(sorted(haps, key=natkey)), sorted(alleles, key=natkey)


def called_star(gene):
    """Return list of (label, haps_canon, allele_list, caller) for the gene."""
    r = calls.get(gene)
    if not r:
        return []
    out = []
    if r.get("pharmcat_diplotype"):
        for lab in r["pharmcat_diplotype"].split(" | "):
            h, al = canon_dip(lab)
            out.append((lab, h, al, "PharmCAT" + (" (outside call: Aldy)" if r.get("pharmcat_call_source") == "OUTSIDE" else "")))
    elif r.get("aldy_diplotype"):
        h, al = canon_dip(r["aldy_diplotype"])
        out.append((r["aldy_diplotype"], h, al, "Aldy"))
    return out


# ---------------------------------------------------------------- OT rows
# Open Targets' aminoglycoside rows carry an rsID but no target and no coordinates, so a
# coordinate-based join cannot place them. The mitochondrial caller already genotyped every CPIC
# allele directly from the CRAM, so those rows are matched by rsID against its own result instead.
MT_BY_RS, MT_OBS = {}, {}
if a.mt_alleles:
    with open(a.mt_alleles) as fh:
        rd = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in rd:
            if (r.get("rsid") or "").strip():
                MT_BY_RS[r["rsid"].strip()] = {"allele": r["name"], "pos": int(r["pos"]),
                                               "ref": r["ref"], "alt": r["alt"]}
if a.mt_rnr1:
    _mt = json.load(open(a.mt_rnr1))
    _thr = (_mt.get("thresholds") or {}).get("heteroplasmic_min", 0.05)
    for x in _mt.get("cpic_alleles") or []:
        MT_OBS[x["allele"]] = {"af": float(x.get("af") or 0), "depth": x.get("depth"),
                               "state": x.get("state"), "detected": float(x.get("af") or 0) >= _thr}


def mt_row(row):
    """The MT-RNR1 allele an rsID-only row refers to, if the mitochondrial caller tested it."""
    rs = (row.get("variantRsId") or "").strip()
    if not rs or row.get("genotypeId") or row.get("haplotypeId"):
        return None
    hit = MT_BY_RS.get(rs)
    return dict(hit, rsid=rs, observed=MT_OBS.get(hit["allele"])) if hit else None


def gene_of(row):
    if mt_row(row):
        return "MT-RNR1", "rsid"
    if row.get("haplotypeId") and "*" in row["haplotypeId"]:
        return row["haplotypeId"].split("*")[0], "haplotypeId"
    if row.get("target") and row["target"].get("approvedSymbol"):
        return row["target"]["approvedSymbol"], "approvedSymbol"
    if row.get("targetFromSourceId"):
        return row["targetFromSourceId"], "targetFromSourceId"
    return "UNKNOWN", "none"


def match_row(row):
    """-> dict(status in matched/not_matched/not_genotyped/not_called/unparsed, child, key, nonref)"""
    gene, _ = gene_of(row)
    if row.get("haplotypeId"):
        star = called_star(gene)
        g = (row.get("genotype") or "").strip()
        if not star:
            return {"status": "not_called", "key": row["haplotypeId"], "child": None, "nonref": None}
        child_label = " | ".join(s[0] for s in star)
        nonref = any(any(x not in ("*1", "Reference") for x in s[2]) for s in star)
        if "/" in g:
            rh, _ = canon_dip(g)
            ok = any(rh == s[1] for s in star)
            kind = "diplotype"
        else:
            ok = any(g in s[2] or g in s[1] for s in star)
            kind = "allele"
        return {"status": "matched" if ok else "not_matched", "key": f"{row['haplotypeId']} ({kind})",
                "child": child_label, "nonref": nonref, "caller": star[0][3]}
    mt = mt_row(row)
    if mt:
        obs = mt.get("observed")
        if obs is None:
            return {"status": "not_genotyped", "key": f"{mt['rsid']} ({mt['allele']})", "child": None, "nonref": None}
        # haploid mitochondrial call: the child carries the alternate only above the heteroplasmy threshold
        child_base = mt["alt"] if obs["detected"] else mt["ref"]
        g = (row.get("genotype") or "").strip()
        return {"status": "matched" if g == child_base else "not_matched",
                "key": f"{mt['rsid']} ({mt['allele']})",
                "child": f"{child_base} ({mt['allele']} {obs['state']}, AF {obs['af']:.4f}, depth {obs['depth']}x)",
                "nonref": obs["detected"], "caller": "MT_RNR1_CALL (CRAM)"}
    gid = row.get("genotypeId") or ""
    m = re.match(r"^([0-9XYM]+)_(\d+)_([^_]+)_(.+)$", gid)
    if not m:
        return {"status": "unparsed", "key": row.get("variantRsId") or gid, "child": None, "nonref": None}
    chrom, pos, ref, alls = m.group(1), int(m.group(2)), m.group(3), m.group(4).split(",")
    rs = row.get("variantRsId") or f"{chrom}:{pos}"
    c = pos_gt.get((chrom, pos))
    if not c or not c["alleles"]:
        return {"status": "not_genotyped", "key": rs, "child": None, "nonref": None}
    flipped = False
    if c["ref"][: len(ref)] != ref and c["ref"][: len(ref)] == ref.translate(COMP):
        alls, ref, flipped = [x.translate(COMP) for x in alls], ref.translate(COMP), True
    child = sorted(c["alleles"])
    row_alleles = sorted(alls)
    if len(child) == 1 and len(set(row_alleles)) == 1:
        ok = child[0] == row_alleles[0]            # haploid call vs homozygous row
    else:
        ok = child == row_alleles
    return {"status": "matched" if ok else "not_matched", "key": rs, "child": "/".join(c["alleles"]) + f" (GT {c['gt']}, DP {c['dp']})",
            "nonref": any(x != c["ref"] for x in c["alleles"]), "flipped": flipped}


# ---------------------------------------------------------------- PharmCAT guideline lookup
report = json.load(open(a.pharmcat_report))


def guideline(drug_names, gene):
    hits = []
    for source, drugs in (report.get("drugs") or {}).items():
        for name in drug_names:
            d = drugs.get(name.lower())
            if not d:
                continue
            for gl in d.get("guidelines") or []:
                for ann in gl.get("annotations") or []:
                    genes = {k for lk in (ann.get("lookupKey") or []) for k in (lk.keys() if isinstance(lk, dict) else [])}
                    if gene not in genes:
                        continue
                    hits.append({"source": source, "classification": ann.get("classification"),
                                 "recommendation": ann.get("drugRecommendation"),
                                 "implications": " ".join(ann.get("implications") or []),
                                 "alternate": bool(ann.get("alternateDrugAvailable")),
                                 "dosing": bool(ann.get("dosingInformation")),
                                 "other": bool(ann.get("otherPrescribingGuidance"))})
    return hits


# ---------------------------------------------------------------- reverse: drug -> child
meta = json.load(open(a.meta))
verdict_rows, level_counts, cov = [], [], Counter()
below_asserted = {}
flips = Counter()

# axis, tier and anti-target verdict per candidate, keyed by ChEMBL ID (empty for a plain ID list)
cand_spec = json.load(open(a.candidates)) if a.candidates else {}
CAND = cand_spec.get("candidates") or {}


def axis_tier(chembl):
    c = CAND.get(chembl) or {}
    return {"axis": c.get("axis"), "role": c.get("role"), "tier": c.get("tier"),
            "patient_support": c.get("patient_support")}
for path in a.drugs:
    d = json.load(open(path))
    qid = d["query_id"]
    rows = d.get("rows") or []
    names = sorted({n for n in [d.get("query_name")] + [r["_source_name"] for r in rows] if n})
    resolved = "; ".join(f"{x['id']} ({x['relation']})" for x in d.get("resolution") or []) or "none"
    lv = Counter(r.get("evidenceLevel") for r in rows)
    # rows below the asserted levels are counted, never joined; the genes behind level-2 rows are named
    # so a reader can see what a lower assertion threshold (params.ot_assert_levels) would bring in
    below = {str(k): v for k, v in lv.items() if k not in ASSERT}
    genes_2 = sorted({gene_of(r)[0] for r in rows if str(r.get("evidenceLevel")) in ("2A", "2B")})
    level_counts.append({"drug": qid, "drug_name": d.get("query_name"), "resolved_to": resolved, "rows_total": len(rows),
                         **{f"level_{k}": v for k, v in sorted(lv.items(), key=lambda x: str(x[0]))},
                         "direct_target_rows": sum(1 for r in rows if r.get("isDirectTarget")),
                         "rows_below_asserted": sum(below.values()),
                         "genes_level_2": "; ".join(genes_2)})
    below_asserted[qid] = {"by_level": dict(sorted(below.items())), "genes_level_2": genes_2}
    cov["drugs_total"] += 1
    if rows:
        cov["drugs_with_any_pgx"] += 1
    if any(r.get("evidenceLevel") in ASSERT for r in rows):
        cov["drugs_with_asserted_pgx"] += 1

    by_gene = defaultdict(list)
    for r in rows:
        if r.get("evidenceLevel") not in ASSERT:
            continue
        if r.get("isDirectTarget") and not a.include_direct_targets:
            cov["asserted_rows_excluded_direct_target"] += 1
            continue
        by_gene[gene_of(r)[0]].append(r)

    if not by_gene:
        verdict_rows.append({"drug": qid, "drug_name": "; ".join(names) or None, **axis_tier(qid),
                             "resolved_to": resolved, "gene": None, "called": None, "caller": None, "matched_ot_allele": None,
                             "evidence_level": None, "pgx_category": None, "flag": "no data",
                             "flag_basis": "no PGx rows at asserted levels" if rows else "no PGx rows in Open Targets",
                             "guideline_text": None, "ot_phenotype_text": None})
        continue

    for gene, grows in sorted(by_gene.items()):
        results = [(r, match_row(r)) for r in grows]
        for _, m in results:
            if m.get("flipped"):
                flips[qid] += 1
        gl = guideline(names, gene)
        matched = [(r, m) for r, m in results if m["status"] == "matched"]
        genotyped = [m for _, m in results if m["status"] in ("matched", "not_matched")]
        if gl:
            flag = "contraindicated" if any(g["alternate"] for g in gl) else ("dose-adjust" if any(g["dosing"] for g in gl) else "no flag")
            basis = "PharmCAT guideline (" + ", ".join(sorted({g["source"].split()[0] for g in gl})) + ")"
            gtext = " || ".join(f"[{g['source'].split()[0]}; {g['classification']}] {g['recommendation']}" for g in gl)
        elif not genotyped:
            flag, basis, gtext = "no data", "asserted rows present but gene/positions not called", None
        elif any(m["nonref"] for m in genotyped):
            flag, basis, gtext = "review (annotation only)", "non-reference at an asserted OT key; no guideline annotation", None
        else:
            flag, basis, gtext = "no flag", "reference at all asserted OT keys; no guideline annotation", None
        star_call = " | ".join(s[0] for s in called_star(gene)) or None
        out_rows = matched or [(None, {"key": None})]
        for r, m in out_rows:
            # the genotype reported on each row is the child's genotype AT THAT ROW'S KEY
            # (star diplotype for haplotype rows, the position's genotype for rsID rows)
            if r is None:
                called_here = star_call or "; ".join(sorted({f"{mm['key']}: {mm['child']}" for _, mm in results if mm.get("child")})) or None
                caller_here = called_star(gene)[0][3] if called_star(gene) else ("force-call" if genotyped else None)
            else:
                called_here = m.get("child")
                caller_here = m.get("caller") or "force-call (CRAM)"
            verdict_rows.append({
                "drug": qid, "drug_name": "; ".join(names), **axis_tier(qid),
                "resolved_to": resolved, "gene": gene,
                "called": called_here, "caller": caller_here,
                "matched_ot_allele": (f"{m['key']} = {r.get('genotype')}" if r else None),
                "evidence_level": r.get("evidenceLevel") if r else None,
                "pgx_category": r.get("pgxCategory") if r else None,
                "flag": flag, "flag_basis": basis, "guideline_text": gtext,
                "ot_phenotype_text": (r.get("phenotypeText") if r else None)})

# candidates the pipeline deliberately did not query: the anti-target exclusions are the evidence
# that the filter works, and an unresolved name is a gap in the list, so both are reported as rows.
# A regimen candidate the rules matched was queried anyway (the rules do not apply to treatment
# already given), so it has real rows above and no placeholder here
for chembl, c in sorted(CAND.items(), key=lambda kv: (str(kv[1].get("axis")), kv[1].get("name") or "")):
    if c.get("anti_target") != "excluded" or c.get("role") == "regimen":
        continue
    verdict_rows.append({"drug": chembl, "drug_name": c.get("matched_name") or c.get("name"),
                         "axis": c.get("axis"), "role": c.get("role"), "tier": c.get("tier"),
                         "patient_support": c.get("patient_support"), "resolved_to": f"{chembl} (as_given)",
                         "gene": None, "called": None, "caller": None, "matched_ot_allele": None,
                         "evidence_level": None, "pgx_category": None,
                         "flag": "excluded (anti-target)", "flag_basis": c.get("anti_target_reason"),
                         "guideline_text": None, "ot_phenotype_text": None})
for u in cand_spec.get("unresolved") or []:
    verdict_rows.append({"drug": None, "drug_name": u.get("name"), "axis": None, "role": None, "tier": None,
                         "patient_support": None,
                         "resolved_to": "none", "gene": None, "called": None, "caller": None,
                         "matched_ot_allele": None, "evidence_level": None, "pgx_category": None,
                         "flag": f"not queried ({u.get('status')})", "flag_basis": u.get("note"),
                         "guideline_text": None, "ot_phenotype_text": None})

# ---------------------------------------------------------------- forward: child -> annotations
targets = json.load(open(a.targets))["genes"]
fwd_rows = []
for gene, t in sorted(targets.items()):
    call = calls.get(gene, {})
    has_call = bool(call.get("pharmcat_diplotype") or call.get("aldy_diplotype"))
    cov["called_genes"] += int(has_call)
    matched_any = False
    for r in t.get("rows") or []:
        if r.get("evidenceLevel") not in ASSERT:
            continue
        m = match_row(r)
        if m["status"] != "matched":
            continue
        matched_any = True
        drugs = "; ".join(sorted({(x.get("drug") or {}).get("name") or x.get("drugFromSource") or "?" for x in r.get("drugs") or []}))
        fwd_rows.append({"gene": gene, "ensembl_id": t["ensembl_id"], "called": m["child"], "matched_ot_allele": f"{m['key']} = {r.get('genotype')}",
                         "evidence_level": r.get("evidenceLevel"), "pgx_category": r.get("pgxCategory"), "drugs": drugs,
                         "ot_phenotype_text": r.get("phenotypeText")})
    if has_call and matched_any:
        cov["called_genes_with_asserted_ot_match"] += 1

# ---------------------------------------------------------------- write
def write_tsv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in fields})


write_tsv(a.prefix + ".verdicts.tsv", verdict_rows)
write_tsv(a.prefix + ".forward_annotations.tsv", fwd_rows)
write_tsv(a.prefix + ".evidence_levels.tsv", level_counts)
summary = {
    # which candidate list produced these verdicts, so a report never shows an example list as the Track 2 set
    "drug_list": {"label": a.drug_list_label, "example": a.drug_list_example,
                  "chembl_ids_queried": sorted({json.load(open(f)).get("query_id") for f in a.drugs}),
                  **{k: v for k, v in cand_spec.items() if k in ("axes", "counts", "unresolved", "excluded",
                                                                 "regimen_anti_target", "anti_target_rules",
                                                                 "candidates")}},
    "open_targets": meta,
    "asserted_levels": sorted(ASSERT),
    "include_direct_targets": a.include_direct_targets,
    "join_coverage": dict(cov),
    "flag_counts": dict(Counter(r["flag"] for r in verdict_rows)),
    "candidates_by_axis": {str(ax): sorted({(r.get("drug_name") or "?") for r in verdict_rows if r.get("axis") == ax})
                           for ax in sorted({r.get("axis") for r in verdict_rows if r.get("axis") is not None},
                                            key=str)},
    "strand_flipped_rows": dict(flips),
    "level3_rows_exploratory_total": sum(r.get("level_3", 0) for r in level_counts),
    # every row at a level outside ASSERT (2A, 2B, 3, 4, ...): counted, carried no further
    "rows_below_asserted_total": sum(r.get("rows_below_asserted", 0) for r in level_counts),
    "below_asserted_by_drug": below_asserted,
}
json.dump(summary, open(a.prefix + ".summary.json", "w"), indent=2)
print(json.dumps(summary["join_coverage"]))
