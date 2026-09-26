#!/usr/bin/env python3
"""Build a SpliceAI gene-annotation file for the panel genes from Ensembl MANE Select transcripts.

Why: SpliceAI 1.3.1's built-in 'grch38' annotation is GENCODE v24 (2016). For this panel it uses
outdated symbols (SGO1 appears as SGOL1) and shorter gene spans than current Ensembl (e.g. ANAPC1,
BUB3, ANAPC5, ANAPC10, PPP2R5C), so 641 intronic/splice variants in the panel went unscored. The
MANE Select transcript per gene matches the transcript the VEP report uses.

Format (SpliceAI utils.Annotator): tab-separated with header
  #NAME CHROM STRAND TX_START TX_END EXON_START EXON_END
TX_START and EXON_START are 0-based (SpliceAI adds 1), ends are 1-based, and exon lists end with ','.
CHROM uses no-chr names, matching the reference.

Usage: make_spliceai_annotation.py assets/mechanistic_panel_grch38.provenance.tsv out.txt
"""
import csv
import json
import sys
import urllib.request

SERVER = "https://rest.ensembl.org"


def post(path, payload):
    req = urllib.request.Request(SERVER + path, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.load(r)


prov, out_path = sys.argv[1:3]
genes = list(csv.DictReader(open(prov), delimiter="\t"))
req = urllib.request.Request(SERVER + "/info/data", headers={"Content-Type": "application/json"})
release = json.load(urllib.request.urlopen(req, timeout=60))["releases"][0]
looked = post("/lookup/id?expand=1;mane=1", {"ids": [g["ensembl_id"] for g in genes]})

rows, notes = [], []
for g in genes:
    hit = looked.get(g["ensembl_id"])
    txs = (hit or {}).get("Transcript") or []
    mane = [t for t in txs if any(m.get("type") == "MANE_Select" for m in (t.get("MANE") or []))]
    chosen = mane[0] if mane else next((t for t in txs if t.get("is_canonical")), None)
    if not chosen:
        sys.exit(f"{g['symbol']}: no MANE Select or canonical transcript in Ensembl {release}")
    basis = "MANE_Select" if mane else "Ensembl canonical"
    exons = sorted(chosen["Exon"], key=lambda e: e["start"])
    rows.append([g["symbol"], hit["seq_region_name"], "+" if hit["strand"] == 1 else "-",
                 chosen["start"] - 1, chosen["end"],
                 "".join(f"{e['start'] - 1}," for e in exons), "".join(f"{e['end']}," for e in exons)])
    notes.append(f"{g['symbol']}\t{chosen['id']}\t{basis}\t{len(exons)} exons")

with open(out_path, "w") as out:
    out.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
with open(out_path.replace(".txt", ".provenance.tsv"), "w") as out:
    out.write(f"# Ensembl release {release}; transcripts used for the SpliceAI annotation\n")
    out.write("symbol\ttranscript\tbasis\texons\n" + "\n".join(notes) + "\n")
print(json.dumps({"ensembl_release": release, "genes": len(rows),
                  "mane_select": sum(1 for n in notes if "MANE_Select" in n)}))
