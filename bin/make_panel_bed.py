#!/usr/bin/env python3
"""Build the mechanistic modifier-panel BED (GRCh38, gene +/- flank) from Ensembl REST.

Genes given with an Ensembl ID are looked up by ID, and the returned symbol is checked against the
expected one (a mismatch fails the run). Genes given by symbol only are looked up by symbol. Only
public gene annotation is requested. Output contigs use the no-chr naming of the pipeline
reference. The Ensembl release is recorded in the BED header and in the provenance TSV.

Usage: make_panel_bed.py assets/mechanistic_panel_genes.tsv 10000 assets/mechanistic_panel_grch38.bed
  genes TSV columns: group, symbol, ensembl_id (may be empty)
"""
import csv
import json
import sys
import urllib.request

SERVER = "https://rest.ensembl.org"


def get(path):
    req = urllib.request.Request(SERVER + path, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def post(path, payload):
    req = urllib.request.Request(SERVER + path, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


genes_tsv, flank, out_bed = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rows = list(csv.DictReader(open(genes_tsv), delimiter="\t"))
release = get("/info/data")["releases"][0]
assembly = get("/info/assembly/homo_sapiens")["assembly_name"]

by_id = post("/lookup/id", {"ids": [r["ensembl_id"] for r in rows if r["ensembl_id"]]})
by_sym = post("/lookup/symbol/homo_sapiens", {"symbols": [r["symbol"] for r in rows if not r["ensembl_id"]]})

records, problems = [], []
for r in rows:
    hit = by_id.get(r["ensembl_id"]) if r["ensembl_id"] else by_sym.get(r["symbol"])
    if not hit:
        problems.append(f"{r['symbol']} {r['ensembl_id']}: not found in Ensembl {release}")
        continue
    if r["ensembl_id"] and hit.get("display_name") != r["symbol"]:
        problems.append(f"{r['ensembl_id']}: Ensembl symbol is {hit.get('display_name')}, panel says {r['symbol']}")
    if hit.get("assembly_name") != "GRCh38" or hit["seq_region_name"] not in [str(i) for i in range(1, 23)] + ["X", "Y"]:
        problems.append(f"{r['symbol']}: not on a GRCh38 primary chromosome ({hit.get('seq_region_name')})")
        continue
    records.append({"chrom": hit["seq_region_name"], "start0": max(0, hit["start"] - 1 - flank), "end": hit["end"] + flank,
                    "symbol": r["symbol"], "ensembl_id": hit["id"], "group": r["group"], "strand": "+" if hit["strand"] == 1 else "-",
                    "gene_start": hit["start"], "gene_end": hit["end"], "biotype": hit.get("biotype")})

if problems:
    sys.exit("Panel lookup problems:\n  " + "\n  ".join(problems))

order = {c: i for i, c in enumerate([str(i) for i in range(1, 23)] + ["X", "Y"])}
records.sort(key=lambda x: (order[x["chrom"]], x["start0"]))
with open(out_bed, "w") as out:
    out.write(f"#ensembl_release={release} assembly={assembly} flank_bp={flank} source={SERVER}\n")
    for x in records:
        out.write(f"{x['chrom']}\t{x['start0']}\t{x['end']}\t{x['symbol']}\t{x['ensembl_id']}\t{x['group']}\t{x['strand']}\n")
with open(out_bed.replace(".bed", ".provenance.tsv"), "w") as out:
    w = csv.DictWriter(out, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(records)
print(json.dumps({"ensembl_release": release, "assembly": assembly, "genes": len(records),
                  "total_bp": sum(x["end"] - x["start0"] for x in records)}))
