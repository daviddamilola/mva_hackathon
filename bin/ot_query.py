#!/usr/bin/env python3
"""Open Targets Platform GraphQL queries for pharmacogenomics (stage 1b).

Subcommands
  meta                        release/API version -> JSON
  literature --diseases IDs   publication count per disease (an unmet-need measure) -> JSON
  tractability --gene SYMBOL  can the gene itself be drugged? essentiality + tractability -> JSON
  drug   --chembl ID          drug(chemblId){pharmacogenomics}. If the ID has no PGx rows, its parent
                              molecule and then child molecules are queried; the resolution is recorded,
                              never applied silently.
  target --genes A,B,...      symbol -> Ensembl via mapIds, then target(ensemblId){pharmacogenomics}

Only public reference data is requested. No subject data is sent to the API.
"""
import argparse
import json
import sys
import time
import urllib.request

API = "https://api.platform.opentargets.org/api/v4/graphql"

PGX_FIELDS = """
  haplotypeId haplotypeFromSourceId genotype genotypeId genotypeAnnotationText
  variantRsId variantId evidenceLevel pgxCategory phenotypeText phenotypeFromSourceId
  isDirectTarget targetFromSourceId datasourceId studyId
  target { id approvedSymbol }
  drugs { drugId drugFromSource drug { id name } }
"""


def gql(query, variables=None, retries=4):
    body = json.dumps({"query": query, "variables": variables or {}}).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(API, data=body, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                out = json.load(r)
            if out.get("errors"):
                raise RuntimeError(json.dumps(out["errors"])[:500])
            return out["data"]
        except Exception as e:  # noqa: BLE001 - network retries
            if attempt == retries - 1:
                raise
            time.sleep(5 * (attempt + 1))
            print(f"retry {attempt + 1}: {e}", file=sys.stderr)


def cmd_meta(_args):
    d = gql("{ meta { name apiVersion { x y z } dataVersion { year month iteration } } }")["meta"]
    av, dv = d["apiVersion"], d["dataVersion"]
    return {"api": d["name"], "api_version": f"{av['x']}.{av['y']}.{av['z']}",
            "data_release": f"{dv['year']}.{dv['month']}" + (f".{dv['iteration']}" if dv.get("iteration") else ""),
            "endpoint": API, "queried_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


DRUG_Q = "query($id: String!) { drug(chemblId: $id) { id name parentMolecule { id name } childMolecules { id name } pharmacogenomics { %s } } }" % PGX_FIELDS


def cmd_drug(args):
    q = gql(DRUG_Q, {"id": args.chembl})["drug"]
    result = {"query_id": args.chembl, "found": q is not None, "resolution": [], "rows": []}
    if q is None:
        return result
    result["query_name"] = q["name"]
    tried = [(q, "as_given")]
    if not q["pharmacogenomics"]:
        if q.get("parentMolecule"):
            tried.append((gql(DRUG_Q, {"id": q["parentMolecule"]["id"]})["drug"], "parent"))
        if not any(t[0] and t[0]["pharmacogenomics"] for t in tried):
            for c in q.get("childMolecules") or []:
                tried.append((gql(DRUG_Q, {"id": c["id"]})["drug"], "child"))
    for d, relation in tried:
        if d and d["pharmacogenomics"]:
            result["resolution"].append({"id": d["id"], "name": d["name"], "relation": relation,
                                         "n_rows": len(d["pharmacogenomics"])})
            for row in d["pharmacogenomics"]:
                row["_source_chembl"] = d["id"]
                row["_source_name"] = d["name"]
                result["rows"].append(row)
    return result


LIT_Q = """query($ids: [String!]!) {
  diseases(efoIds: $ids) { id name literatureOcurrences { count } }
}"""


def cmd_literature(args):
    ids = [x.strip() for x in args.diseases.split(",") if x.strip()]
    rows = gql(LIT_Q, {"ids": ids})["diseases"] or []
    return {"queried_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "diseases": [{"id": d["id"], "name": d["name"],
                          "publications": (d.get("literatureOcurrences") or {}).get("count")}
                         for d in rows]}


TRACT_Q = """query($e: String!) {
  target(ensemblId: $e) {
    id approvedSymbol biotype
    isEssential
    tractability { label modality value }
  }
}"""

MODALITY = {"SM": "small molecule", "AB": "antibody", "PR": "PROTAC / degrader",
            "OC": "other clinical modality"}


def cmd_tractability(args):
    """Can the causal gene itself be drugged? The answer decides whether a repurposing strategy has
    to act on something other than the target, which is the premise of the whole candidate list."""
    ens = symbol_to_ensembl([args.gene])[args.gene]
    if not ens:
        return {"gene": args.gene, "ensembl_id": None, "found": False}
    t = gql(TRACT_Q, {"e": ens})["target"] or {}
    tr = t.get("tractability") or []
    by_modality = {}
    for m in sorted({x["modality"] for x in tr}):
        rows = [x for x in tr if x["modality"] == m]
        by_modality[m] = {"name": MODALITY.get(m, m),
                          "assessed": len(rows),
                          "positive": [x["label"] for x in rows if x["value"]],
                          "negative": [x["label"] for x in rows if not x["value"]]}
    clinical = [lbl for lbl in ("Approved Drug", "Advanced Clinical", "Phase 1 Clinical")
                if any(x["label"] == lbl and x["value"] for x in tr)]
    return {"gene": t.get("approvedSymbol") or args.gene, "ensembl_id": ens, "found": True,
            "biotype": t.get("biotype"), "is_essential": t.get("isEssential"),
            "clinically_tractable_labels": clinical,
            "has_clinical_stage_drug": bool(clinical),
            "by_modality": by_modality,
            "queried_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def symbol_to_ensembl(genes):
    m = gql("query($t: [String!]!) { mapIds(queryTerms: $t, entityNames: [\"target\"]) { mappings { term hits { id name entity } } } }",
            {"t": genes})["mapIds"]["mappings"]
    out = {}
    for mp in m:
        hits = [h for h in mp["hits"] if h["entity"] == "target" and h["name"] == mp["term"]]
        out[mp["term"]] = hits[0]["id"] if hits else None
    return out


def cmd_target(args):
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    m = gql("query($t: [String!]!) { mapIds(queryTerms: $t, entityNames: [\"target\"]) { mappings { term hits { id name entity } } } }",
            {"t": genes})["mapIds"]["mappings"]
    out = {"genes": {}}
    for mp in m:
        hits = [h for h in mp["hits"] if h["entity"] == "target" and h["name"] == mp["term"]]
        ens = hits[0]["id"] if hits else None
        entry = {"symbol": mp["term"], "ensembl_id": ens, "rows": []}
        if ens:
            t = gql("query($e: String!) { target(ensemblId: $e) { id approvedSymbol pharmacogenomics { %s } } }" % PGX_FIELDS,
                    {"e": ens})["target"]
            entry["rows"] = (t or {}).get("pharmacogenomics") or []
        out["genes"][mp["term"]] = entry
    return out


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
sub.add_parser("meta")
p = sub.add_parser("drug"); p.add_argument("--chembl", required=True)
p = sub.add_parser("target"); p.add_argument("--genes", required=True)
p = sub.add_parser("literature"); p.add_argument("--diseases", required=True)
p = sub.add_parser("tractability"); p.add_argument("--gene", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
res = {"meta": cmd_meta, "drug": cmd_drug, "target": cmd_target, "literature": cmd_literature,
       "tractability": cmd_tractability}[args.cmd](args)
with open(args.out, "w") as fh:
    json.dump(res, fh, indent=1)
print(f"{args.cmd}: written {args.out}")
