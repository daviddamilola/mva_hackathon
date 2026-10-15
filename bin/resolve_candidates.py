#!/usr/bin/env python3
"""Resolve a candidate drug list to ChEMBL IDs and apply the anti-target rules (stage 1b input).

The candidate list names drugs; it never carries ChEMBL IDs. Every name is resolved here through the
Open Targets search API, so the list stays readable and the IDs stay reproducible. A name that does
not resolve is REPORTED with a status, never silently dropped:

  resolved    exactly one drug entity matched the name (case-insensitive) or one hit dominates
  ambiguous   several hits match the name equally well; the candidate is not carried forward
  not_found   the search returned no drug entity

ANTI-TARGET RULES (from the candidate file's `anti_target_rules`)
  Rule 1, automated. A drug is excluded when a mechanism-of-action row has an inhibitory action type
    AND one of that row's targets carries a direct Gene Ontology biological-process annotation to one
    of the listed mitotic-fidelity terms. Encodes: an inhibitor of a mitotic-fidelity gene is
    contraindicated in a mitotic-fidelity disorder.
  Rule 2, curated. Targets whose inhibition is shown in primary literature to compromise mitotic
    fidelity but which carry none of those GO terms. Each entry cites a PMID, which is printed in the
    exclusion reason so a reader can check it.

A drug with no mechanism-of-action row in Open Targets cannot be tested by either rule; it is marked
'not assessable' and carried forward with that stated, because absence of a mechanism annotation is
not evidence of safety.

AXIS ROLES. An axis is `repurposing` (the default) or `regimen`. A regimen axis holds treatment the
child is already receiving. The anti-target rules govern proposing NEW drugs for a mitotic-fidelity
disorder; they are not applied to therapy given with curative intent and do not license withdrawing
it (docs/08). So a regimen candidate is queried for its pharmacogenomic annotation whatever the rules
say. The verdict is still computed and recorded for it, as a mechanistic note (it is the basis of the
therapy-related-risk consideration), never as an exclusion.

Usage: resolve_candidates.py --candidates drugs.yml --prefix out
"""
import argparse
import csv
import json
import sys
import time
import urllib.request

try:
    import yaml
except ModuleNotFoundError:  # the task environment pins PyYAML; JSON input still works without it
    yaml = None

API = "https://api.platform.opentargets.org/api/v4/graphql"

ap = argparse.ArgumentParser()
ap.add_argument("--candidates", required=True)
ap.add_argument("--prefix", required=True)
ap.add_argument("--min-score", type=float, default=10.0)
ap.add_argument("--fixture", help="canned API responses {search, drug, go} instead of calling the API (tests)")
a = ap.parse_args()


FIX = json.load(open(a.fixture)) if a.fixture else None


def gql(query, variables=None, retries=4):
    if FIX is not None:  # offline: serve the canned response for this query shape
        if "search(" in query:
            return {"search": {"hits": FIX["search"].get(variables["q"], [])}}
        if "drug(chemblId" in query:
            return {"drug": FIX["drug"].get(variables["id"])}
        return {"target": {"approvedSymbol": variables["e"],
                           "geneOntology": [{"aspect": "P", "term": {"id": t, "label": t}}
                                            for t in FIX["go"].get(variables["e"], [])]}}
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


SEARCH_Q = """query($q: String!) {
  search(queryString: $q, entityNames: ["drug"], page: {index: 0, size: 10}) {
    hits { id name entity score }
  }
}"""

MOA_Q = """query($id: String!) {
  drug(chemblId: $id) {
    id name maximumClinicalStage drugType
    parentMolecule { id name }
    mechanismsOfAction { rows { mechanismOfAction actionType targets { id approvedSymbol } } }
  }
}"""

GO_Q = """query($e: String!) {
  target(ensemblId: $e) { approvedSymbol geneOntology { aspect term { id label } } }
}"""

text = open(a.candidates).read()
if yaml is not None:
    spec = yaml.safe_load(text)
else:
    try:
        spec = json.loads(text)
    except json.JSONDecodeError:
        raise SystemExit(f"{a.candidates}: PyYAML is not installed, so only JSON input can be parsed. "
                         "Install pyyaml, or run this through the pipeline, whose environment pins it.")
rules = spec.get("anti_target_rules") or {}
GO_TERMS = {str(k): v for k, v in (rules.get("go_terms") or {}).items()}
INHIBITORY = {str(x).upper() for x in (rules.get("inhibitory_action_types") or [])}
LIT = rules.get("literature") or []
LIT_TARGETS = {t.upper(): r for r in LIT for t in r.get("targets") or []}

# how well each candidate is supported by evidence about THIS child; the report groups by this so a
# syndrome-level or purely mechanistic proposal is never presented as a patient-specific one
SUPPORT = {"released_feature", "syndrome_only", "mechanism_only"}

# axis roles: a `regimen` axis is treatment already given, queried regardless of the anti-target verdict
ROLES = {"repurposing", "regimen"}
AXES_SPEC = {str(k): (v or {}) for k, v in (spec.get("axes") or {}).items()}


def axis_role(axis):
    role = str(AXES_SPEC.get(str(axis), {}).get("role") or "repurposing")
    if role not in ROLES:
        raise SystemExit(f"axis {axis!r}: role must be one of {sorted(ROLES)}, got {role!r}")
    return role

go_cache = {}


def go_bp(ensembl_id):
    """Direct GO biological-process term ids for a target."""
    if ensembl_id not in go_cache:
        t = gql(GO_Q, {"e": ensembl_id})["target"] or {}
        go_cache[ensembl_id] = {g["term"]["id"] for g in (t.get("geneOntology") or []) if g["aspect"] == "P"}
    return go_cache[ensembl_id]


def resolve(name):
    hits = gql(SEARCH_Q, {"q": name})["search"]["hits"]
    hits = [h for h in hits if h["entity"] == "drug"]
    if not hits:
        return None, "not_found", "Open Targets search returned no drug entity", []
    exact = [h for h in hits if h["name"].strip().lower() == name.strip().lower()]
    if len(exact) == 1:
        return exact[0], "resolved", "exact name match", hits[:3]
    if len(exact) > 1:
        return None, "ambiguous", f"{len(exact)} drug entities share this name: " + \
            ", ".join(f"{h['id']} ({h['name']})" for h in exact[:5]), hits[:5]
    top = hits[0]
    second = hits[1]["score"] if len(hits) > 1 else 0.0
    if top["score"] >= a.min_score and top["score"] >= 2 * second:
        return top, "resolved", f"no exact match; top hit {top['name']} dominates (score {top['score']:.0f} vs {second:.0f})", hits[:3]
    return None, "ambiguous", "no exact name match and no dominant hit: " + \
        ", ".join(f"{h['id']} ({h['name']}, {h['score']:.0f})" for h in hits[:5]), hits[:5]


resolved_rows, anti_rows, candidates = [], [], {}
for c in spec.get("candidates") or []:
    name = c["name"]
    hit, status, note, alts = resolve(name)
    support = c.get("patient_support")
    if support not in SUPPORT:
        raise SystemExit(f"candidate {name!r}: patient_support must be one of {sorted(SUPPORT)}, got {support!r}")
    role = axis_role(c.get("axis"))
    row = {"name": name, "axis": c.get("axis"), "role": role, "tier": c.get("tier"), "patient_support": support,
           "caution": " ".join((c.get("caution") or "").split()),
           "chembl_id": (hit or {}).get("id"), "matched_name": (hit or {}).get("name"),
           "match_score": round((hit or {}).get("score", 0.0), 1) if hit else "",
           "status": status, "note": note,
           "alternatives": "; ".join(f"{h['id']} ({h['name']})" for h in alts) if status != "resolved" else "",
           "investigational": bool(c.get("investigational")),
           "max_clinical_stage": "", "action_types": "", "mechanism_targets": "",
           "anti_target": "", "anti_target_reason": ""}

    if hit:
        d = gql(MOA_Q, {"id": hit["id"]})["drug"] or {}
        row["max_clinical_stage"] = d.get("maximumClinicalStage") or ""
        moa_rows = ((d.get("mechanismsOfAction") or {}).get("rows")) or []
        row["action_types"] = "; ".join(sorted({r.get("actionType") or "?" for r in moa_rows}))
        targets = sorted({t["approvedSymbol"] for r in moa_rows for t in r.get("targets") or []})
        row["mechanism_targets"] = "; ".join(targets[:12]) + (" ..." if len(targets) > 12 else "")
        excluded = []
        for r in moa_rows:
            action = (r.get("actionType") or "").upper()
            inhibitory = action in INHIBITORY
            for t in r.get("targets") or []:
                sym, ens = t["approvedSymbol"], t["id"]
                hit_go = sorted(go_bp(ens) & set(GO_TERMS)) if inhibitory else []
                lit = LIT_TARGETS.get(sym.upper()) if inhibitory else None
                rule = ("GO mitotic-fidelity" if hit_go else ("curated literature" if lit else ""))
                reason = ""
                if hit_go:
                    reason = (f"{action} of {sym}, annotated to " +
                              "; ".join(f"{g} {GO_TERMS[g]}" for g in hit_go))
                elif lit:
                    reason = f"{action} of {sym}: {lit['statement'].strip()} ({lit['citation']}, PMID {lit['pmid']})"
                anti_rows.append({"name": name, "chembl_id": hit["id"], "target": sym, "ensembl_id": ens,
                                  "action_type": r.get("actionType"), "mechanism": r.get("mechanismOfAction"),
                                  "rule": rule or "none matched", "verdict": "excluded" if rule else "allowed",
                                  "reason": reason})
                if rule:
                    excluded.append(reason)
        if not moa_rows:
            row["anti_target"] = "not assessable"
            row["anti_target_reason"] = "no mechanism-of-action row in Open Targets; neither rule can be applied"
        elif excluded:
            row["anti_target"] = "excluded"
            row["anti_target_reason"] = " || ".join(dict.fromkeys(excluded))
        else:
            row["anti_target"] = "allowed"
            row["anti_target_reason"] = "no inhibitory mechanism against a mitotic-fidelity target"
        candidates[hit["id"]] = {"name": name, "matched_name": hit["name"], "axis": c.get("axis"), "role": role,
                                 "tier": c.get("tier"), "patient_support": support, "anti_target": row["anti_target"],
                                 "anti_target_reason": row["anti_target_reason"],
                                 "max_clinical_stage": row["max_clinical_stage"],
                                 "investigational": bool(c.get("investigational")),
                                 "caution": " ".join((c.get("caution") or "").split()),
                                 "move_up": " ".join((c.get("move_up") or "").split())}
    resolved_rows.append(row)


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


# only resolved candidates that the anti-target rules did not exclude are queried in stage 1b; every
# other candidate stays in the tables with its status. Regimen candidates are the exception: the rules
# do not apply to treatment already given, so they are queried whatever the verdict
queried = [r["chembl_id"] for r in resolved_rows
           if r["status"] == "resolved" and (r["anti_target"] != "excluded" or r["role"] == "regimen")]
open(a.prefix + ".query_ids.txt", "w").write("\n".join(dict.fromkeys(queried)) + ("\n" if queried else ""))
write_tsv(a.prefix + ".resolved.tsv", resolved_rows)
write_tsv(a.prefix + ".anti_target.tsv", anti_rows)

# `order` sets the sequence the report presents axes in: the mechanism-directed axes lead, because the
# brief scores whether a candidate is supported by the variant mechanism, not merely by plausible biology
axes = {str(k): {"title": v.get("title"), "role": axis_role(k), "order": v.get("order"),
                 "rationale": " ".join((v.get("rationale") or "").split())}
        for k, v in AXES_SPEC.items()}
repurposing = [r for r in resolved_rows if r["role"] == "repurposing"]
regimen = [r for r in resolved_rows if r["role"] == "regimen"]
summary = {
    "drug_list": {"label": spec.get("drug_list_label"), "example": bool(spec.get("drug_list_example"))},
    "axes": axes,
    "counts": {
        "candidates": len(resolved_rows),
        "resolved": sum(1 for r in resolved_rows if r["status"] == "resolved"),
        "ambiguous": sum(1 for r in resolved_rows if r["status"] == "ambiguous"),
        "not_found": sum(1 for r in resolved_rows if r["status"] == "not_found"),
        "excluded_by_anti_target": sum(1 for r in repurposing if r["anti_target"] == "excluded"),
        "anti_target_not_assessable": sum(1 for r in resolved_rows if r["anti_target"] == "not assessable"),
        "queried_in_stage_1b": len(queried),
        "regimen": len(regimen),
        "regimen_anti_target_hits": sum(1 for r in regimen if r["anti_target"] == "excluded"),
    },
    "by_patient_support": {k: sorted(r["name"] for r in resolved_rows if r["patient_support"] == k)
                           for k in sorted(SUPPORT)},
    "unresolved": [{"name": r["name"], "status": r["status"], "note": r["note"]}
                   for r in resolved_rows if r["status"] != "resolved"],
    "excluded": [{"name": r["name"], "axis": r["axis"], "reason": r["anti_target_reason"]}
                 for r in repurposing if r["anti_target"] == "excluded"],
    # regimen agents the rules would have matched: reported as a mechanistic note, and still queried
    "regimen_anti_target": [{"name": r["name"], "axis": r["axis"], "reason": r["anti_target_reason"]}
                            for r in regimen if r["anti_target"] == "excluded"],
    "anti_target_rules": {"go_terms": GO_TERMS, "inhibitory_action_types": sorted(INHIBITORY),
                          "literature": [{"targets": r["targets"], "pmid": r["pmid"], "citation": r["citation"]}
                                         for r in LIT]},
    "candidates": candidates,
}
json.dump(summary, open(a.prefix + ".candidates.json", "w"), indent=2)
print(json.dumps({k: summary[k] for k in ("counts", "unresolved", "excluded")}, indent=1))
