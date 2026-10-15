#!/usr/bin/env python3
"""Simulation test for bin/aneuploidy_model.py (synthetic data only).

A synthetic XY genome (7 chromosomes, 100 kb bins, ~45x, overdispersed allele counts, reference bias)
carries planted mosaic events. The model must:
  - call the planted events significant with the correct model, and a cell-fraction CI containing the truth
  - not call the unaltered arms
  - report chrX non-PAR LRR near -1 (internal control)
  - NOT call an arm-level event on 2p, which carries a 3 Mb hotspot of paralogue-like skewed 'het'
    sites (40% of sites at 0.25/0.75) outside any segdup annotation: the artefact that produced 13
    false arm calls in the first real run
Run: python3 tests/unit/test_aneuploidy_model.py
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODEL = os.path.join(ROOT, "bin", "aneuploidy_model.py")
rng = np.random.default_rng(42)

CHR = {"1": 60_000_000, "2": 60_000_000, "3": 60_000_000, "5": 60_000_000, "9": 60_000_000, "18": 60_000_000, "X": 60_000_000}
CEN = {c: (28_000_000, 30_000_000) for c in CHR}
# planted events: arm -> (model, cell fraction)
TRUTH = {"5q": ("gain", 0.20), "18p": ("loss", 0.10), "9p": ("cnloh", 0.30), "3q": ("loss", 0.08)}
HOTSPOT = ("2", 10_000_000, 13_000_000)   # PSV-like artefact on an otherwise normal arm
BIN, DEPTH, RHO, MBIAS = 100_000, 45.0, 0.004, 0.485


def state(c, pos):
    arm = f"{c}{'p' if pos < CEN[c][0] else 'q'}"
    return TRUTH.get(arm, (None, 0.0))


def main():
    with tempfile.TemporaryDirectory() as t:
        with open(f"{t}/cyto.txt", "w") as fh:
            for c, L in CHR.items():
                fh.write(f"{c}\t0\t{CEN[c][0]}\tp11\tgneg\n{c}\t{CEN[c][0]}\t{CEN[c][1]}\tp11.1\tacen\n{c}\t{CEN[c][1]}\t{L}\tq11\tgneg\n")
        with open(f"{t}/bins.tsv", "w") as fh:
            fh.write("chrom\tstart\tend\tgc\tmap_frac\tblacklist_bp\tdepth\n")
            for c, L in CHR.items():
                for s in range(0, L, BIN):
                    if CEN[c][0] <= s < CEN[c][1]:
                        continue
                    gc = float(np.clip(rng.normal(0.41, 0.05), 0.3, 0.6))
                    gc_effect = 1 + 0.8 * (gc - 0.41)             # GC bias to be corrected
                    mdl, f = state(c, s)
                    cn = 2 + (f if mdl == "gain" else -f if mdl == "loss" else 0.0)
                    if c == "X" and s > 2_800_000:
                        cn = 1.0                                   # XY sample
                    dp = rng.gamma(400, (DEPTH * cn / 2 * gc_effect) / 400)
                    fh.write(f"{c}\t{s}\t{s + BIN}\t{gc:.4f}\t{rng.uniform(0.92, 1):.3f}\t0\t{dp:.2f}\n")
        with gzip.open(f"{t}/s.vcf.gz", "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
            for c, L in CHR.items():
                if c == "X":
                    continue
                for p in np.sort(rng.choice(np.arange(1, L), size=L // 1500, replace=False)):
                    if CEN[c][0] <= p < CEN[c][1]:
                        continue
                    mdl, f = state(c, p)
                    if mdl == "gain":
                        mu = (1 + f) / (2 + f)
                    elif mdl == "loss":
                        mu = 1 / (2 - f)
                    elif mdl == "cnloh":
                        mu = (1 + f) / 2
                    else:
                        mu = 0.5
                    if c == HOTSPOT[0] and HOTSPOT[1] <= p < HOTSPOT[2] and rng.random() < 0.4:
                        mu = 0.25                                   # paralogous sequence variant
                    mu = mu if rng.random() < 0.5 else 1 - mu      # unknown phase
                    mu = mu + (MBIAS - 0.5)                         # reference bias
                    cn = 2 + (f if mdl == "gain" else -f if mdl == "loss" else 0.0)
                    n = rng.poisson(DEPTH * cn / 2)
                    ab = (1 - RHO) / RHO
                    k = rng.binomial(n, rng.beta(mu * ab, (1 - mu) * ab))
                    fh.write(f"{c}\t{p}\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP\t0/1:{n - k},{k}:{n}\n")
        open(f"{t}/segdup.bed", "w").close()     # hotspot deliberately NOT annotated as segdup
        res = subprocess.run([sys.executable, MODEL, "--bins", f"{t}/bins.tsv", "--vcf", f"{t}/s.vcf.gz", "--cytoband", f"{t}/cyto.txt",
                              "--segdup", f"{t}/segdup.bed", "--sex", "M", "--prefix", f"{t}/out"], capture_output=True, text=True)
        if res.returncode:
            print(res.stderr[-2000:]); sys.exit(1)
        arms = {r["arm"]: r for r in json.load(open(f"{t}/out.arms.json"))}
        summ = json.load(open(f"{t}/out.summary.json"))

    fails = []

    def check(name, ok):
        print(("PASS " if ok else "FAIL ") + name)
        if not ok:
            fails.append(name)

    for A, (mdl, f) in TRUTH.items():
        r = arms[A]
        check(f"{A}: planted {mdl} {f:.0%} called significant", bool(r.get("significant")))
        check(f"{A}: model = {mdl} (got {r.get('best_model')})", r.get("best_model") == mdl)
        ci = r.get("cell_fraction_ci95") or [None, None]
        check(f"{A}: CI {ci} contains {f}", ci[0] is not None and ci[0] <= f <= ci[1])
    normal = [A for A, r in arms.items() if r.get("tested") and A not in TRUTH]
    fp = [A for A in normal if arms[A].get("significant")]
    check(f"no false positives among {len(normal)} unaltered arms (got {fp})", not fp)
    check(f"2p hotspot artefact not called an arm-level event (status: {arms['2p'].get('status')})",
          arms["2p"].get("status") != "arm-level event")
    xc = summ["internal_control_chrX_nonPAR_lrr"]
    check(f"chrX non-PAR LRR near -1 (got {xc:.3f})", abs(xc + 1) < 0.1)
    print(f"\n{len(fails)} failure(s)")
    sys.exit(1 if fails else 0)


main()
