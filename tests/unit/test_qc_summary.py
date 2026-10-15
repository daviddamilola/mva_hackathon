#!/usr/bin/env python3
"""Synthetic tests for bin/qc_summary.py row selection (no subject data).

somalier's ancestry table lists every labelled 1000 Genomes reference sample alongside the query
sample. The first QC summary took the first row and reported a reference sample's ancestry
(AFR, p = 0.9999). These tests fail if the selected row is anything other than the sample's own.
Run:  python3 tests/unit/test_qc_summary.py
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
QC = os.path.join(ROOT, "bin", "qc_summary.py")
SAMPLE = "QUERY1"
ANC_COLS = ["#sample_id", "predicted_ancestry", "given_ancestry", "AFR_prob", "EUR_prob", "EAS_prob", "AMR_prob", "SAS_prob"]


def tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols); w.writerows(rows)


def run(tmp, anc_rows, relate_ids=(SAMPLE,), selfsm_ids=(SAMPLE,)):
    tsv(f"{tmp}/anc.tsv", ANC_COLS, anc_rows)
    tsv(f"{tmp}/relate.tsv", ["sample_id", "depth_mean", "X_depth_mean", "X_n", "X_het", "Y_depth_mean", "p_middling_ab"],
        [[s, 40, 20, 300, 0, 18, 0.01] for s in relate_ids])
    tsv(f"{tmp}/selfsm.tsv", ["#SEQ_ID", "#SNPS", "FREEMIX", "AVG_DP"], [[s, 100000, 0.001, 40] for s in selfsm_ids])
    args = ["python3", QC, "--sample", SAMPLE, "--declared-sex", "M", "--out", f"{tmp}/qc.json",
            "--somalier-ancestry", f"{tmp}/anc.tsv", "--somalier-samples", f"{tmp}/relate.tsv", "--selfsm", f"{tmp}/selfsm.tsv"]
    for x in ("flagstat", "markdup", "mosdepth", "roh", "aneuploidy"):
        args += [f"--{x}", f"{tmp}/absent.{x}"]
    p = subprocess.run(args, capture_output=True, text=True)
    return p, (json.load(open(f"{tmp}/qc.json")) if p.returncode == 0 else None)


# reference rows first and last, each confidently assigned to a different superpopulation; the
# query row sits between them with a distinctive, non-saturated probability
REFS_BEFORE = [["HG03082", "AFR", "AFR", 0.9999, 0, 0, 0.0001, 0], ["HG03914", "SAS", "SAS", 0, 0, 0, 0, 1.0]]
REFS_AFTER = [["HG01843", "EAS", "EAS", 0, 0, 1.0, 0, 0]]
QUERY = [SAMPLE, "EUR", "", 0, 0.956, 0, 0.044, 0]


def test_selects_own_row():
    with tempfile.TemporaryDirectory() as tmp:
        p, qc = run(tmp, REFS_BEFORE + [QUERY] + REFS_AFTER)
        assert p.returncode == 0, p.stderr
        anc = qc["ancestry"]
        assert anc["predicted"] == "EUR", anc
        assert anc["probabilities"]["EUR"] == 0.956, anc
        conf = next(c for c in qc["checks"] if c["check"] == "ancestry assignment confidence")
        assert conf["value"] == 0.956, conf


def test_row_order_does_not_matter():
    with tempfile.TemporaryDirectory() as tmp:
        for rows in ([QUERY] + REFS_BEFORE + REFS_AFTER, REFS_BEFORE + REFS_AFTER + [QUERY]):
            p, qc = run(tmp, rows)
            assert p.returncode == 0, p.stderr
            assert qc["ancestry"]["probabilities"]["EUR"] == 0.956, qc["ancestry"]


def test_missing_sample_fails():
    with tempfile.TemporaryDirectory() as tmp:
        p, _ = run(tmp, REFS_BEFORE + REFS_AFTER)
        assert p.returncode != 0 and "expected one row" in (p.stderr + p.stdout), p.stderr


def test_duplicate_sample_fails():
    with tempfile.TemporaryDirectory() as tmp:
        p, _ = run(tmp, REFS_BEFORE + [QUERY, QUERY])
        assert p.returncode != 0, "two rows for the sample must not be resolved silently"


def test_relate_and_selfsm_select_own_row():
    with tempfile.TemporaryDirectory() as tmp:
        p, _ = run(tmp, [QUERY], relate_ids=("HG03082",), selfsm_ids=("HG03082",))
        assert p.returncode != 0, "relate/selfSM without the sample's row must fail"


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print(f"PASS {name}")
            except AssertionError as e:
                fails += 1; print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
