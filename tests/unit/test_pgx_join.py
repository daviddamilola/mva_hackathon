#!/usr/bin/env python3
"""Synthetic tests for bin/pgx_join.py join rules (no subject data).

Each case targets a rule that would otherwise fail silently. Run:  python3 tests/unit/test_pgx_join.py
"""
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
JOIN = os.path.join(ROOT, "bin", "pgx_join.py")


def row(**kw):
    base = {"haplotypeId": None, "genotype": None, "genotypeId": None, "variantRsId": None, "variantId": None,
            "evidenceLevel": "1A", "pgxCategory": "toxicity", "phenotypeText": "synthetic", "isDirectTarget": False,
            "targetFromSourceId": None, "target": None, "drugs": [], "_source_chembl": "CHEMBL1", "_source_name": "DRUGA"}
    base.update(kw)
    return base


def vcf(path, records):
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for c, p, ref, alt, gt in records:
            fh.write(f"{c}\t{p}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:DP\t{gt}:30\n")


def run(tmp, calls, report, drugs, records, extra_args=()):
    json.dump({"calls": calls}, open(f"{tmp}/calls.json", "w"))
    json.dump(report, open(f"{tmp}/report.json", "w"))
    json.dump({"api_version": "t", "data_release": "t"}, open(f"{tmp}/meta.json", "w"))
    json.dump({"genes": {}}, open(f"{tmp}/targets.json", "w"))
    paths = []
    for i, d in enumerate(drugs):
        p = f"{tmp}/drug{i}.json"
        json.dump(d, open(p, "w"))
        paths.append(p)
    vcf(f"{tmp}/a.vcf.gz", records)
    vcf(f"{tmp}/b.vcf.gz", [])
    subprocess.run([sys.executable, JOIN, "--calls", f"{tmp}/calls.json", "--pharmcat-report", f"{tmp}/report.json",
                    "--vcf", f"{tmp}/a.vcf.gz", f"{tmp}/b.vcf.gz", "--drugs", *paths, "--targets", f"{tmp}/targets.json",
                    "--meta", f"{tmp}/meta.json", "--prefix", f"{tmp}/out", *extra_args], check=True, capture_output=True)
    rows = list(csv.DictReader(open(f"{tmp}/out.verdicts.tsv"), delimiter="\t"))
    return rows, json.load(open(f"{tmp}/out.summary.json"))


def ann(gene, alternate=False, dosing=False):
    return {"classification": "Strong", "drugRecommendation": "synthetic rec", "implications": [],
            "alternateDrugAvailable": alternate, "dosingInformation": dosing, "otherPrescribingGuidance": False,
            "lookupKey": [{gene: "x"}]}


def report_for(drug_name, gene, **kw):
    return {"drugs": {"CPIC Guideline Annotation": {drug_name: {"guidelines": [{"annotations": [ann(gene, **kw)]}]}}}}


failures = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        failures.append(name)


star_calls = [{"gene": "TPMT", "pharmcat_diplotype": "*3A/*1", "pharmcat_call_source": "MATCHER"},
              {"gene": "CYP2D6", "pharmcat_diplotype": "*4/*35x2", "pharmcat_call_source": "OUTSIDE"}]

with tempfile.TemporaryDirectory() as tmp:
    drug = {"query_id": "CHEMBL1", "query_name": "DRUGA", "resolution": [{"id": "CHEMBL1", "relation": "as_given"}], "rows": [
        row(haplotypeId="TPMT*1/*3A", genotype="*1/*3A"),                                  # diplotype, order differs from call
        row(haplotypeId="TPMT*3A", genotype="*3A"),                                        # single allele
        row(haplotypeId="TPMT*3B", genotype="*3B"),                                        # suballele: must NOT match *3A
        row(haplotypeId="CYP2D6*35", genotype="*35"),                                      # allele inside xN duplication
        row(haplotypeId="CYP2D6*4/*35x2", genotype="*4/*35x2"),                            # duplication diplotype
        row(variantRsId="rs1", genotype="AG", genotypeId="1_100_G_A,G", target={"approvedSymbol": "GENEA"}),
        row(variantRsId="rs2", genotype="GG", genotypeId="1_200_C_C,C", target={"approvedSymbol": "GENEB"}),  # reverse strand
        row(variantRsId="rs3", genotype="GG", genotypeId="X_300_G_G,G", target={"approvedSymbol": "GENEX"}),  # haploid chrX
        row(variantRsId="rs4", genotype="TT", genotypeId="1_400_T_T,T", target={"approvedSymbol": "GENEC"}, isDirectTarget=True),
        row(variantRsId="rs5", genotype="AA", genotypeId="1_500_A_A,A", target={"approvedSymbol": "GENED"}, evidenceLevel="3"),
    ]}
    records = [("1", 100, "G", "A", "0/1"), ("1", 200, "G", ".", "0/0"), ("X", 300, "G", ".", "0"),
               ("1", 400, "T", ".", "0/0"), ("1", 500, "A", ".", "0/0")]
    rows, summ = run(tmp, star_calls, report_for("druga", "TPMT", dosing=True), [drug], records)
    m = {(r["gene"], r["matched_ot_allele"]) for r in rows}
    check("diplotype matched despite allele order", ("TPMT", "TPMT*1/*3A (diplotype) = *1/*3A") in m)
    check("single allele *3A matched", ("TPMT", "TPMT*3A (allele) = *3A") in m)
    check("suballele *3B not collapsed into *3A", not any("*3B" in k[1] for k in m))
    check("allele inside CYP2D6 xN duplication matched", ("CYP2D6", "CYP2D6*35 (allele) = *35") in m)
    check("CYP2D6 duplication diplotype matched", ("CYP2D6", "CYP2D6*4/*35x2 (diplotype) = *4/*35x2") in m)
    check("rsID het matched regardless of allele order", ("GENEA", "rs1 = AG") in m)
    check("reverse-strand row flipped and matched", ("GENEB", "rs2 = GG") in m)
    check("strand flip counted", summ["strand_flipped_rows"].get("CHEMBL1") == 1)
    check("haploid chrX call matched homozygous row", ("GENEX", "rs3 = GG") in m)
    check("isDirectTarget row excluded", not any(r["gene"] == "GENEC" for r in rows))
    check("isDirectTarget exclusion counted", summ["join_coverage"].get("asserted_rows_excluded_direct_target") == 1)
    check("level-3 row not asserted", not any(r["gene"] == "GENED" for r in rows))
    check("level-3 row counted as exploratory", summ["level3_rows_exploratory_total"] == 1)
    tp = [r for r in rows if r["gene"] == "TPMT"]
    check("guideline dosingInformation -> dose-adjust", all(r["flag"] == "dose-adjust" for r in tp))
    ga = [r for r in rows if r["gene"] == "GENEA"]
    check("no guideline + non-ref at asserted key -> review (annotation only)", ga and ga[0]["flag"] == "review (annotation only)")
    gb = [r for r in rows if r["gene"] == "GENEB"]
    check("no guideline + reference -> no flag", gb and gb[0]["flag"] == "no flag")
    check("per-row 'called' is the genotype at that row's key",
          ga and ga[0]["called"].startswith("G/A") and gb and gb[0]["called"].startswith("G/G"))

with tempfile.TemporaryDirectory() as tmp:
    drug = {"query_id": "CHEMBL1", "query_name": "DRUGA", "resolution": [], "rows": [row(haplotypeId="TPMT*3A", genotype="*3A")]}
    rows, _ = run(tmp, star_calls, report_for("druga", "TPMT", alternate=True), [drug], [])
    check("guideline alternateDrugAvailable -> contraindicated", rows[0]["flag"] == "contraindicated")
    empty = {"query_id": "CHEMBL2", "query_name": "DRUGB", "resolution": [], "rows": []}
    rows, summ = run(tmp, star_calls, {"drugs": {}}, [empty], [])
    check("drug with no PGx rows -> no data", rows[0]["flag"] == "no data")
    check("coverage counts drugs without PGx", summ["join_coverage"]["drugs_total"] == 1 and "drugs_with_any_pgx" not in summ["join_coverage"])
    nullt = {"query_id": "CHEMBL3", "query_name": "DRUGC", "resolution": [], "rows": [row(variantRsId="rs9", genotype="AA", genotypeId="1_900_A_A,A")]}
    rows, _ = run(tmp, star_calls, {"drugs": {}}, [nullt], [("1", 900, "A", ".", "0/0")])
    check("null target and targetFromSourceId handled (gene UNKNOWN)", rows[0]["gene"] == "UNKNOWN" and rows[0]["flag"] == "no flag")

print(f"\n{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
