#!/usr/bin/env python3
"""Synthetic tests for bin/resolve_candidates.py: name resolution and the anti-target rules.

No network and no subject data: API responses come from a fixture, so the rule logic is what is
under test. Each case targets a rule whose failure would be silent.
Run:  python3 tests/unit/test_resolve_candidates.py
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESOLVE = os.path.join(ROOT, "bin", "resolve_candidates.py")

# GO ids the rules treat as mitotic-fidelity, and one that is not
SEG = "GO:0007059"
SAC = "GO:0007094"
OTHER = "GO:0006915"

RULES = {
    "go_terms": {SEG: "chromosome segregation", SAC: "mitotic spindle assembly checkpoint signaling"},
    "inhibitory_action_types": ["INHIBITOR", "ANTAGONIST", "DEGRADER"],
    "literature": [{"targets": ["MTOR", "FKBP1A"], "pmid": "35141499",
                    "citation": "Yamada C et al., iScience 2022;25:103675",
                    "statement": "mTORC1 inactivation induces mitotic slippage."}],
}


def drug(cid, name, action, targets, stage="APPROVAL"):
    return {"id": cid, "name": name, "maximumClinicalStage": stage, "drugType": "Small molecule",
            "mechanismsOfAction": {"rows": [{"mechanismOfAction": f"{action} of {targets[0]}",
                                             "actionType": action,
                                             "targets": [{"id": f"ENS_{t}", "approvedSymbol": t} for t in targets]}]}}


FIXTURE = {
    "search": {
        "inhibits_a_segregation_gene": [{"id": "CHEMBL_A", "name": "inhibits_a_segregation_gene", "entity": "drug", "score": 900.0}],
        "activates_a_segregation_gene": [{"id": "CHEMBL_B", "name": "activates_a_segregation_gene", "entity": "drug", "score": 900.0}],
        "inhibits_an_unrelated_gene": [{"id": "CHEMBL_C", "name": "inhibits_an_unrelated_gene", "entity": "drug", "score": 900.0}],
        "hits_the_literature_rule": [{"id": "CHEMBL_D", "name": "hits_the_literature_rule", "entity": "drug", "score": 900.0}],
        "no_mechanism_recorded": [{"id": "CHEMBL_E", "name": "no_mechanism_recorded", "entity": "drug", "score": 900.0}],
        "case insensitive NAME": [{"id": "CHEMBL_F", "name": "Case Insensitive Name", "entity": "drug", "score": 900.0}],
        "two_drugs_share_this_name": [{"id": "CHEMBL_G1", "name": "two_drugs_share_this_name", "entity": "drug", "score": 900.0},
                                      {"id": "CHEMBL_G2", "name": "two_drugs_share_this_name", "entity": "drug", "score": 880.0}],
        "weak_and_contested": [{"id": "CHEMBL_H1", "name": "something else", "entity": "drug", "score": 30.0},
                               {"id": "CHEMBL_H2", "name": "another thing", "entity": "drug", "score": 28.0}],
        "dominant_near_match": [{"id": "CHEMBL_I", "name": "dominant near match", "entity": "drug", "score": 800.0},
                                {"id": "CHEMBL_I2", "name": "distant", "entity": "drug", "score": 10.0}],
        "absent_everywhere": [],
    },
    "drug": {
        "CHEMBL_A": drug("CHEMBL_A", "inhibits_a_segregation_gene", "INHIBITOR", ["SEGGENE"]),
        "CHEMBL_B": drug("CHEMBL_B", "activates_a_segregation_gene", "ACTIVATOR", ["SEGGENE"]),
        "CHEMBL_C": drug("CHEMBL_C", "inhibits_an_unrelated_gene", "INHIBITOR", ["PLAINGENE"]),
        "CHEMBL_D": drug("CHEMBL_D", "hits_the_literature_rule", "INHIBITOR", ["FKBP1A"]),
        "CHEMBL_E": {"id": "CHEMBL_E", "name": "no_mechanism_recorded", "maximumClinicalStage": "PHASE_2",
                     "drugType": "Small molecule", "mechanismsOfAction": {"rows": []}},
        "CHEMBL_F": drug("CHEMBL_F", "Case Insensitive Name", "INHIBITOR", ["PLAINGENE"]),
        "CHEMBL_I": drug("CHEMBL_I", "dominant near match", "INHIBITOR", ["PLAINGENE"]),
    },
    "go": {"ENS_SEGGENE": [SEG, OTHER], "ENS_PLAINGENE": [OTHER], "ENS_FKBP1A": [OTHER]},
}

def cand(name, axis, tier, support="mechanism_only"):
    return {"name": name, "axis": axis, "tier": tier, "move_up": "n/a", "patient_support": support}


CANDIDATES = [
    cand("inhibits_a_segregation_gene", 7, "anti"),
    cand("activates_a_segregation_gene", 1, 2, "released_feature"),
    cand("inhibits_an_unrelated_gene", 1, 2, "syndrome_only"),
    cand("hits_the_literature_rule", 7, "anti"),
    cand("no_mechanism_recorded", 1, 3),
    cand("case insensitive NAME", 1, 3),
    cand("two_drugs_share_this_name", 1, 3),
    cand("weak_and_contested", 1, 3),
    cand("dominant_near_match", 1, 3),
    cand("absent_everywhere", 1, 3),
]


def run(tmp, candidates=None, rules=None):
    spec = {"drug_list_label": "unit_test", "drug_list_example": True,
            "axes": {1: {"title": "t", "rationale": "r"}, 7: {"title": "anti", "rationale": "r"}},
            "candidates": candidates if candidates is not None else CANDIDATES,
            "anti_target_rules": rules if rules is not None else RULES}
    json.dump(spec, open(f"{tmp}/spec.yml", "w"))          # JSON is valid YAML
    json.dump(FIXTURE, open(f"{tmp}/fix.json", "w"))
    p = subprocess.run(["python3", RESOLVE, "--candidates", f"{tmp}/spec.yml", "--fixture", f"{tmp}/fix.json",
                        "--prefix", f"{tmp}/out"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    summary = json.load(open(f"{tmp}/out.candidates.json"))
    rows = {r.split("\t")[0]: r.split("\t") for r in open(f"{tmp}/out.resolved.tsv").read().splitlines()[1:]}
    header = open(f"{tmp}/out.resolved.tsv").read().splitlines()[0].split("\t")
    queried = [x for x in open(f"{tmp}/out.query_ids.txt").read().splitlines() if x]
    return summary, {k: dict(zip(header, v)) for k, v in rows.items()}, queried


def test_inhibitor_of_mitotic_fidelity_gene_excluded():
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, queried = run(tmp)
        r = rows["inhibits_a_segregation_gene"]
        assert r["anti_target"] == "excluded", r
        assert "chromosome segregation" in r["anti_target_reason"], r["anti_target_reason"]
        assert "CHEMBL_A" not in queried, "an excluded candidate must not be queried"


def test_non_inhibitory_action_on_same_gene_allowed():
    """The rule is about inhibition, not about the target: an activator must survive."""
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, _ = run(tmp)
        assert rows["activates_a_segregation_gene"]["anti_target"] == "allowed", rows["activates_a_segregation_gene"]


def test_inhibitor_of_unrelated_gene_allowed():
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, queried = run(tmp)
        assert rows["inhibits_an_unrelated_gene"]["anti_target"] == "allowed"
        assert "CHEMBL_C" in queried


def test_literature_rule_excludes_and_cites():
    """A target with no mitotic-fidelity GO term is caught only by the curated rule, with its PMID."""
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, _ = run(tmp)
        r = rows["hits_the_literature_rule"]
        assert r["anti_target"] == "excluded", r
        assert "35141499" in r["anti_target_reason"], r["anti_target_reason"]


def test_no_mechanism_is_not_assessable_not_allowed():
    """Absence of a mechanism annotation is not evidence of safety."""
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, queried = run(tmp)
        assert rows["no_mechanism_recorded"]["anti_target"] == "not assessable", rows["no_mechanism_recorded"]
        assert "CHEMBL_E" in queried, "not-assessable candidates are still queried, with the caveat recorded"


def test_name_match_is_case_insensitive():
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, _ = run(tmp)
        assert rows["case insensitive NAME"]["status"] == "resolved"
        assert rows["case insensitive NAME"]["chembl_id"] == "CHEMBL_F"


def test_ambiguous_and_missing_names_are_reported_not_dropped():
    with tempfile.TemporaryDirectory() as tmp:
        summary, rows, queried = run(tmp)
        assert rows["two_drugs_share_this_name"]["status"] == "ambiguous"
        assert rows["weak_and_contested"]["status"] == "ambiguous"
        assert rows["absent_everywhere"]["status"] == "not_found"
        reported = {u["name"] for u in summary["unresolved"]}
        assert {"two_drugs_share_this_name", "weak_and_contested", "absent_everywhere"} <= reported, reported
        assert summary["counts"]["candidates"] == len(CANDIDATES)
        for cid in ("CHEMBL_G1", "CHEMBL_G2", "CHEMBL_H1"):
            assert cid not in queried, f"{cid} was not resolved and must not be queried"


def test_dominant_hit_resolves_without_exact_name():
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, _ = run(tmp)
        r = rows["dominant_near_match"]
        assert r["status"] == "resolved" and r["chembl_id"] == "CHEMBL_I", r


def test_axis_tier_and_support_are_carried_through():
    with tempfile.TemporaryDirectory() as tmp:
        summary, _, _ = run(tmp)
        c = summary["candidates"]["CHEMBL_A"]
        assert (c["axis"], c["tier"], c["patient_support"]) == (7, "anti", "mechanism_only"), c
        assert summary["counts"]["excluded_by_anti_target"] == 2, summary["counts"]
        by = summary["by_patient_support"]
        assert by["released_feature"] == ["activates_a_segregation_gene"], by
        assert by["syndrome_only"] == ["inhibits_an_unrelated_gene"], by


def test_bad_patient_support_is_rejected():
    """A candidate whose support level is missing or misspelt must stop the run, not default."""
    import subprocess as sp
    with tempfile.TemporaryDirectory() as tmp:
        spec = {"drug_list_label": "unit_test", "drug_list_example": True,
                "axes": {1: {"title": "t", "rationale": "r"}},
                "candidates": [{"name": "inhibits_an_unrelated_gene", "axis": 1, "tier": 3,
                                "patient_support": "patient"}],
                "anti_target_rules": RULES}
        json.dump(spec, open(f"{tmp}/spec.yml", "w"))
        json.dump(FIXTURE, open(f"{tmp}/fix.json", "w"))
        p = sp.run(["python3", RESOLVE, "--candidates", f"{tmp}/spec.yml", "--fixture", f"{tmp}/fix.json",
                    "--prefix", f"{tmp}/out"], capture_output=True, text=True)
        assert p.returncode != 0, "an invalid patient_support must fail the run"
        assert "patient_support" in (p.stderr + p.stdout)


def test_regimen_axis_is_queried_despite_anti_target_match():
    """Treatment already given is annotated whatever the rules say; the match is a note, not an exclusion."""
    with tempfile.TemporaryDirectory() as tmp:
        spec = {"drug_list_label": "unit_test", "drug_list_example": True,
                "axes": {1: {"title": "t", "rationale": "r"},
                         "regimen": {"title": "given", "role": "regimen", "rationale": "r"}},
                "candidates": [{"name": "inhibits_a_segregation_gene", "axis": "regimen",
                                "patient_support": "released_feature"},          # no tier: allowed on a regimen axis
                               {"name": "inhibits_a_segregation_gene", "axis": 1, "tier": "anti",
                                "patient_support": "mechanism_only"}],
                "anti_target_rules": RULES}
        json.dump(spec, open(f"{tmp}/spec.yml", "w"))
        json.dump(FIXTURE, open(f"{tmp}/fix.json", "w"))
        p = subprocess.run(["python3", RESOLVE, "--candidates", f"{tmp}/spec.yml", "--fixture", f"{tmp}/fix.json",
                            "--prefix", f"{tmp}/out"], capture_output=True, text=True)
        assert p.returncode == 0, p.stderr
        summary = json.load(open(f"{tmp}/out.candidates.json"))
        queried = [x for x in open(f"{tmp}/out.query_ids.txt").read().splitlines() if x]
        rows = list(csv.DictReader(open(f"{tmp}/out.resolved.tsv"), delimiter="\t"))
        reg = next(r for r in rows if r["role"] == "regimen")
        assert reg["anti_target"] == "excluded", "the verdict is still computed and recorded"
        assert "CHEMBL_A" in queried, "a regimen drug is queried even when the rules match"
        c = summary["counts"]
        assert c["regimen"] == 1 and c["regimen_anti_target_hits"] == 1, c
        assert c["excluded_by_anti_target"] == 1, "only the repurposing entry counts as excluded"
        assert [x["name"] for x in summary["regimen_anti_target"]] == ["inhibits_a_segregation_gene"]
        assert summary["axes"]["regimen"]["role"] == "regimen" and summary["axes"]["1"]["role"] == "repurposing"


def test_unknown_axis_role_is_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        spec = {"drug_list_label": "unit_test", "drug_list_example": True,
                "axes": {1: {"title": "t", "rationale": "r", "role": "current"}},
                "candidates": [cand("inhibits_an_unrelated_gene", 1, 3)], "anti_target_rules": RULES}
        json.dump(spec, open(f"{tmp}/spec.yml", "w"))
        json.dump(FIXTURE, open(f"{tmp}/fix.json", "w"))
        p = subprocess.run(["python3", RESOLVE, "--candidates", f"{tmp}/spec.yml", "--fixture", f"{tmp}/fix.json",
                            "--prefix", f"{tmp}/out"], capture_output=True, text=True)
        assert p.returncode != 0 and "role" in (p.stderr + p.stdout)


def test_empty_rules_exclude_nothing():
    """A list with no anti_target_rules must not silently behave as if the filter ran."""
    with tempfile.TemporaryDirectory() as tmp:
        _, rows, _ = run(tmp, rules={})
        assert rows["inhibits_a_segregation_gene"]["anti_target"] == "allowed"


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print(f"PASS {name}")
            except AssertionError as e:
                fails += 1; print(f"FAIL {name}: {e}")
    print(f"\n{fails} failure(s)")
    sys.exit(1 if fails else 0)
