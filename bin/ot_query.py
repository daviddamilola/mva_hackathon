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
  neighbourhood --genes-tsv F for each gene in a curated TSV: the mitotic-fidelity GO terms it carries
                              and the action type of every clinical candidate against it. The gene set is
                              curated; the drugs and their directions are fetched, never asserted.

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


# the six GO biological-process terms the anti-target rule treats as mitotic fidelity; kept here so the
# neighbourhood statement and the exclusion rule answer to the same definition
MITOTIC_GO = {
    "GO:0007059": "chromosome segregation",
    "GO:0000070": "mitotic sister chromatid segregation",
    "GO:0000819": "sister chromatid segregation",
    "GO:0098813": "nuclear chromosome segregation",
    "GO:0007094": "mitotic spindle assembly checkpoint signaling",
    "GO:0007052": "mitotic spindle organization",
}
# action types the anti-target rule counts as inhibitory, plus the antisense form, which reduces the
# gene product just as directly
INHIBITORY = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "DEGRADER", "DISRUPTING AGENT",
              "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR", "ANTISENSE INHIBITOR"}

NBHD_Q = """query($e: String!) { target(ensemblId: $e) { approvedSymbol
  geneOntology { aspect term { id } }
  drugAndClinicalCandidates { count rows { maxClinicalStage
    drug { id name mechanismsOfAction { rows { actionType targets { approvedSymbol } } } } } } } }"""


def cmd_neighbourhood(args):
    genes = []
    for line in open(args.genes_tsv):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "symbol":
            continue
        genes.append({"symbol": f[0], "ensembl_id": f[1], "role": f[2] if len(f) > 2 else ""})
    out = {"source_file": args.genes_tsv.split("/")[-1], "mitotic_go_terms": MITOTIC_GO,
           "inhibitory_action_types": sorted(INHIBITORY), "genes": []}
    for g in genes:
        t = gql(NBHD_Q, {"e": g["ensembl_id"]})["target"]
        if not t:
            out["genes"].append({**g, "found": False})
            continue
        # a symbol that no longer matches the recorded Ensembl id is reported, never absorbed
        sym = t.get("approvedSymbol")
        go = sorted({x["term"]["id"] for x in (t.get("geneOntology") or []) if x["aspect"] == "P"} & set(MITOTIC_GO))
        dac = t.get("drugAndClinicalCandidates") or {}
        actions, drugs = {}, []
        for r in dac.get("rows") or []:
            d = r.get("drug") or {}
            for m in ((d.get("mechanismsOfAction") or {}).get("rows") or []):
                # only mechanism rows that act on THIS gene; a multi-target drug contributes once per
                # row that names it, and its unrelated rows are not counted here
                if any((x or {}).get("approvedSymbol") == sym for x in m.get("targets") or []):
                    at = (m.get("actionType") or "UNSPECIFIED").upper()
                    actions[at] = actions.get(at, 0) + 1
                    drugs.append({"name": d.get("name"), "chembl_id": d.get("id"),
                                  "action_type": at, "max_clinical_stage": r.get("maxClinicalStage")})
        n = sum(actions.values())
        out["genes"].append({**g, "found": True, "approved_symbol": sym,
                             "symbol_matches": sym == g["symbol"],
                             "mitotic_go": go, "mitotic_go_labels": [MITOTIC_GO[x] for x in go],
                             "n_clinical_candidates": dac.get("count", 0),
                             "action_types": actions, "n_mechanism_rows": n,
                             "n_inhibitory": sum(v for k, v in actions.items() if k in INHIBITORY),
                             "drugs": drugs})
    g = [x for x in out["genes"] if x.get("found")]
    out["summary"] = {
        "genes": len(g),
        "genes_with_candidates": sum(1 for x in g if x["n_clinical_candidates"]),
        "genes_without_candidates": sorted(x["symbol"] for x in g if not x["n_clinical_candidates"]),
        "genes_without_mitotic_go": sorted(x["symbol"] for x in g if not x["mitotic_go"]),
        "mechanism_rows": sum(x["n_mechanism_rows"] for x in g),
        "inhibitory_rows": sum(x["n_inhibitory"] for x in g),
        "non_inhibitory": [{"gene": x["symbol"], "drug": d["name"], "action_type": d["action_type"]}
                           for x in g for d in x["drugs"] if d["action_type"] not in INHIBITORY],
    }
    return out


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
sub.add_parser("meta")
p = sub.add_parser("drug"); p.add_argument("--chembl", required=True)
p = sub.add_parser("target"); p.add_argument("--genes", required=True)
p = sub.add_parser("literature"); p.add_argument("--diseases", required=True)
p = sub.add_parser("tractability"); p.add_argument("--gene", required=True)
p = sub.add_parser("neighbourhood"); p.add_argument("--genes-tsv", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
res = {"meta": cmd_meta, "drug": cmd_drug, "target": cmd_target, "literature": cmd_literature,
       "tractability": cmd_tractability, "neighbourhood": cmd_neighbourhood}[args.cmd](args)
with open(args.out, "w") as fh:
    json.dump(res, fh, indent=1)
print(f"{args.cmd}: written {args.out}")
