#!/usr/bin/env python3
"""Copy-number counterfactual for outside star-allele calls.

What would the recommendations have been if copy number had not been detected? For every outside
call whose diplotype carries a gene duplication (e.g. *35x2, *2xN, *36+*10), the copies are collapsed
to one ('make'); PharmCAT is re-run on the same preprocessed VCF with that call, and the two reports
are compared ('diff'). Only the collapsed genes differ between the two runs.

  make:  pgx_counterfactual.py make <outside_calls.tsv> <out_counterfactual.tsv> <out_changes.json>
  diff:  pgx_counterfactual.py diff <observed.report.json> <counterfactual.report.json> <changes.json> <out.json>

A recommendation 'changes' when the CPIC or DPWG annotation text, classification, or the structured
flags (alternateDrugAvailable, dosingInformation) differ. Flag counts use the structured fields only.
"""
import json
import re
import sys

DUP = re.compile(r"x(\d+|N)$")


def collapse(allele):
    """*35x2 -> *35; *36+*10 -> *10 (tandem: the functional copy PharmCAT would see alone)."""
    a = DUP.sub("", allele)
    return a.split("+")[-1] if "+" in a else a


def make(inp, out, changes):
    rows, changed = [], {}
    for l in open(inp):
        if not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        gene, dip = f[0], f[1]
        if "/" in dip and not gene.startswith("HLA"):
            alleles = dip.split("/")
            cf = "/".join(collapse(x) for x in alleles)
            if cf != dip:
                changed[gene] = {"observed": dip, "counterfactual": cf}
                f[1] = cf
        rows.append("\t".join(f))
    open(out, "w").write("\n".join(rows) + "\n")
    json.dump(changed, open(changes, "w"), indent=2)
    print(json.dumps(changed))


def annotations(report, genes):
    out = {}
    for src, drugs in report["drugs"].items():
        if not src.startswith(("CPIC", "DPWG")):
            continue
        for name, d in drugs.items():
            for gl in d.get("guidelines") or []:
                for a in gl.get("annotations") or []:
                    keys = {k for lk in (a.get("lookupKey") or []) if isinstance(lk, dict) for k in lk}
                    if not keys & genes:
                        continue
                    flag = ("alternative drug" if a.get("alternateDrugAvailable") else
                            "dose adjustment" if a.get("dosingInformation") else None)
                    out[(src.split()[0], name, json.dumps(sorted(keys)))] = {
                        "classification": a.get("classification"), "flag": flag,
                        "recommendation": (a.get("drugRecommendation") or "").strip()}
    return out


def gene_call(report, gene):
    g = report["genes"].get(gene) or {}
    d = (g.get("recommendationDiplotypes") or g.get("sourceDiplotypes") or [{}])[0]
    return {"diplotype": d.get("label"), "phenotype": " | ".join(d.get("phenotypes") or []),
            "activity_score": d.get("activityScore")}


def diff(obs_p, cf_p, changes_p, out):
    changes = json.load(open(changes_p))
    if not changes:
        json.dump({"applicable": False, "reason": "no outside call carries a gene duplication"}, open(out, "w"), indent=2)
        return
    obs, cf = json.load(open(obs_p)), json.load(open(cf_p))
    genes = set(changes)
    ao, ac = annotations(obs, genes), annotations(cf, genes)
    changed = []
    for k in sorted(set(ao) | set(ac)):
        o, c = ao.get(k), ac.get(k)
        if o != c:
            changed.append({"source": k[0], "drug": k[1], "genes": json.loads(k[2]),
                            "observed": o, "counterfactual": c})
    res = {
        "applicable": True,
        "genes": {g: {"observed": gene_call(obs, g), "counterfactual": gene_call(cf, g)} for g in sorted(genes)},
        "annotations_compared": len(set(ao) | set(ac)),
        "annotations_changed": len(changed),
        "flagged_observed": sum(1 for v in ao.values() if v["flag"]),
        "flagged_counterfactual": sum(1 for v in ac.values() if v["flag"]),
        "drugs_flagged_only_in_counterfactual": sorted({c["drug"] for c in changed
                                                        if (c["counterfactual"] or {}).get("flag") and not (c["observed"] or {}).get("flag")}),
        "changed": changed,
    }
    json.dump(res, open(out, "w"), indent=2)
    print(json.dumps({k: v for k, v in res.items() if k != "changed"}, indent=1))


if __name__ == "__main__":
    {"make": make, "diff": diff}[sys.argv[1]](*sys.argv[2:])
