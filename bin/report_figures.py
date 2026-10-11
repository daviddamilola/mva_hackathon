#!/usr/bin/env python3
"""Figures for the Track 2 report, drawn from the pipeline's own outputs.

Every figure is a view of a table that already exists in results/; nothing is computed here that the
report does not also print as numbers. Each figure is independent: a missing or malformed input marks
that figure 'not available' in the manifest and the rest are still drawn. Inputs are resolved by file
name, never by directory layout, like bin/render_report.py.

Palette: the reference data-viz palette (categorical in fixed slot order, one-hue sequential blue,
reserved status colours with a text label beside every status-coloured mark).

Usage: report_figures.py --results res --sample S --outdir figures --locus-genes assets/pgx_locus_genes.tsv
"""
import argparse
import csv
import gzip
import json
import os
import re
import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)
ap.add_argument("--sample", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--locus-genes", required=True)
a = ap.parse_args()
R, S = a.results, a.sample
os.makedirs(a.outdir, exist_ok=True)

# ---------------------------------------------------------------- palette (reference instance)
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
       "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
INK, INK2, MUTED, GRID, AXIS, SURFACE, NEUTRAL = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#f0efec"
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", SEQ)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8.5, "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "axes.labelsize": 8.5, "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": SURFACE, "figure.facecolor": SURFACE,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "grid.color": GRID, "grid.linewidth": 0.6, "legend.fontsize": 8, "legend.frameon": False,
    "savefig.dpi": 200, "savefig.facecolor": SURFACE,
})

# ---------------------------------------------------------------- input resolution (by file name)
INDEX = {}
for root, _, files in os.walk(R, followlinks=True):
    for f in files:
        INDEX.setdefault(f, os.path.join(root, f))


def path(name):
    p = INDEX.get(name)
    if not p or os.path.getsize(p) == 0:
        raise FileNotFoundError(name)
    return p


def jload(name):
    return json.load(open(path(name)))


def tload(name):
    return list(csv.DictReader(open(path(name)), delimiter="\t"))


def save(fig, name):
    out = os.path.join(a.outdir, f"{S}.{name}.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def hgrid(ax):
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)


manifest = {}


def figure(name):
    def wrap(fn):
        try:
            out = fn()
            manifest[name] = {"file": os.path.basename(out), "status": "drawn"}
        except FileNotFoundError as e:
            manifest[name] = {"status": "not available", "reason": f"missing input {e}"}
        except Exception as e:  # noqa: BLE001 - one figure must not take the others down
            manifest[name] = {"status": "failed", "reason": f"{type(e).__name__}: {e}"}
            traceback.print_exc()
        return fn
    return wrap


# ================================================================ 1. three-caller concordance
def norm_call(s):
    """Allele set of a call, so 'Reference' and 'Reference/Reference' compare equal."""
    s = (s or "").strip().lower()
    s = s.replace("reference", "ref").replace("*ref", "ref").replace("(ref)", "ref")
    s = re.sub(r"\s*\((heterozygous|homozygous)\)", "", s)
    s = re.sub(r"[^a-z0-9/;:+x>.-]", "", s)
    return frozenset(x for x in re.split(r"[/;]", s) if x)


@figure("pgx_concordance")
def pgx_concordance():
    calls = tload(f"{S}.pgx_calls.tsv")
    star = {r["gene"]: r for r in tload(f"{S}.aldy_vs_pypgx.tsv")}
    try:
        hla = {r["gene"]: r for r in tload(f"{S}.hla_calls.tsv")}
    except FileNotFoundError:
        hla = {}
    rows = []
    for c in calls:
        g = c["gene"]
        pc, al = c.get("pharmcat_diplotype") or "", c.get("aldy_diplotype") or ""
        py = (star.get(g) or {}).get("pypgx") or ""
        pa, ap_ = c.get("caller_agreement") or "", (star.get(g) or {}).get("status") or ""
        cols = ["PharmCAT", "Aldy", "PyPGx"]
        vals = [pc, al, py]
        if g in hla:
            h = hla[g]
            cols, vals = ["PharmCAT (outside)", "T1K", "OptiType"], [pc, h.get("t1k_2field") or "", h.get("optitype_2field") or ""]
            verdict, status = ("2 callers agree", "good") if h.get("agreement") == "agree" else ("callers differ", "serious")
        elif ap_ == "different":
            verdict, status = "Aldy vs PyPGx differ", "serious"
        elif ap_ in ("same", "same_after_normalisation") and pa == "agree":
            verdict, status = "3 callers agree", "good"
        elif ap_ in ("same", "same_after_normalisation") and pa == "outside_call":
            verdict, status = "Aldy = PyPGx (PharmCAT outside call)", "good"
        elif pa == "not_comparable" or ap_ == "not_comparable":
            verdict, status = "not comparable (nomenclature)", "warning"
        elif pa == "pharmcat_only":
            verdict, status = "PharmCAT only", "neutral"
        elif ap_ == "only_other" and norm_call(pc) == norm_call(py):
            verdict, status = "dedicated caller = PyPGx", "good"
        elif ap_ == "only_other":
            verdict, status = "PyPGx only", "neutral"
        else:
            verdict, status = f"{pa or ap_ or 'n/a'}", "neutral"
        rows.append((g, cols, vals, verdict, status))

    n = len(rows)
    fig, ax = plt.subplots(figsize=(11.5, 0.34 * n + 1.1))
    ax.set_xlim(0, 11.5); ax.set_ylim(0, n + 1); ax.axis("off")
    xcols = [0.0, 1.5, 4.2, 6.9, 9.2]
    for x, h in zip(xcols, ["Gene", "Caller 1", "Caller 2", "Caller 3", "Verdict"]):
        ax.text(x + 0.05, n + 0.45, h, fontsize=8.5, fontweight="bold", color=INK, va="center")
    fill = {"good": STATUS["good"], "warning": STATUS["warning"], "serious": STATUS["serious"], "neutral": NEUTRAL}
    for i, (g, cols, vals, verdict, status) in enumerate(rows):
        y = n - 1 - i
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y), 11.5, 1, facecolor="#f5f5f2", edgecolor="none", zorder=0))
        ax.text(xcols[0] + 0.05, y + 0.5, g, fontsize=8.5, color=INK, va="center", fontweight="bold")
        for x, col, val in zip(xcols[1:4], cols, vals):
            txt = (val if len(val) <= 30 else val[:28] + "…") or "—"
            ax.text(x + 0.05, y + 0.5, txt, fontsize=7.5, color=INK if val else MUTED, va="center", family="monospace")
            ax.text(x + 0.05, y + 0.88, col, fontsize=5.5, color=MUTED, va="center")
        ax.add_patch(Rectangle((xcols[4], y + 0.18), 0.22, 0.64, facecolor=fill[status], edgecolor="none"))
        ax.text(xcols[4] + 0.32, y + 0.5, verdict, fontsize=7.5, color=INK, va="center")
    ax.set_title(f"Pharmacogene calls by caller ({n} genes): the verdict column is the report's agreement status",
                 loc="left")
    handles = [Patch(facecolor=fill[k], label=l) for k, l in
               [("good", "agree"), ("warning", "not comparable"), ("serious", "differ"), ("neutral", "single caller")]]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, -0.06 * 23 / max(n, 1)), ncol=4)
    return save(fig, "pgx_concordance")


# ================================================================ 2. CYP2D6 locus depth
@figure("cyp2d6_depth")
def cyp2d6_depth():
    d = tload(f"{S}.pgx_locus_depth.tsv")
    if not d:
        raise FileNotFoundError("pgx_locus_depth.tsv (empty)")
    qc = jload("genome_qc.json")
    mean = float(qc["mean_depth_autosomes"])
    genes = list(csv.DictReader(open(a.locus_genes), delimiter="\t"))
    calls = {r["gene"]: r for r in tload(f"{S}.pgx_calls.tsv")}
    mid = np.array([(int(r["start"]) + int(r["end"])) / 2 for r in d]) / 1e6
    all_r = np.array([float(r["depth_all"]) for r in d]) / mean
    mq_r = np.array([float(r["depth_mq20"]) for r in d]) / mean
    fig, ax = plt.subplots(figsize=(11, 3.4))
    hgrid(ax)
    for g in genes:
        s, e = int(g["start"]) / 1e6, int(g["end"]) / 1e6
        ax.axvspan(s, e, color=NEUTRAL, zorder=0)
        ax.text((s + e) / 2, 2.35, g["gene"], ha="center", va="top", fontsize=8, color=INK2, style="italic")
    for ratio, lbl in [(1.0, "2 copies"), (1.5, "3 copies"), (0.5, "1 copy")]:
        ax.axhline(ratio, color=AXIS, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.text(mid.min(), ratio + 0.03, lbl, fontsize=7, color=MUTED, va="bottom")
    ax.plot(mid, all_r, color=CAT[0], lw=1.6, label="all reads", zorder=3)
    ax.plot(mid, mq_r, color=CAT[1], lw=1.6, label="MAPQ ≥ 20 only", zorder=3)
    ax.set_ylim(0, 2.4)
    ax.set_xlabel(f"chr{d[0]['chrom']} position (Mb, GRCh38)")
    ax.set_ylabel("depth / autosomal mean")
    call = calls.get("CYP2D6", {})
    ax.set_title(f"CYP2D6 locus read depth in {int(d[0]['end']) - int(d[0]['start']) + 1} bp windows; "
                 f"called {call.get('pharmcat_diplotype') or call.get('aldy_diplotype') or 'n/a'} "
                 f"(autosomal mean {mean:.1f}x)", loc="left")
    ax.legend(loc="upper right", ncol=2)
    ax.text(0.0, -0.3, "CYP2D6 and CYP2D7 are near-identical paralogs: the MAPQ-filtered track drops wherever reads "
            "cannot be placed uniquely, so the two tracks bracket the true copy number rather than fixing it.",
            transform=ax.transAxes, fontsize=7, color=INK2, va="top")
    return save(fig, "cyp2d6_depth")


# ================================================================ 2b. readthrough products
@figure("readthrough_products")
def readthrough_products():
    d = jload(f"{S}.readthrough.json")
    e = (d.get("ptcs") or [None])[0]
    if not e or not e.get("residue_tolerance", {}).get("all"):
        raise FileNotFoundError("readthrough.json (no scored residue)")
    tol = e["residue_tolerance"]["all"]
    prods = {p["aa1"]: p for p in e["readthrough_products"]}
    order = sorted(tol, key=lambda k: tol[k]["score"])
    vals = [tol[k]["score"] for k in order]
    # AlphaMissense class boundaries, so a reader can place every bar without consulting the table
    BENIGN_MAX, PATH_MIN = 0.34, 0.564
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3), gridspec_kw={"width_ratios": [2.5, 1]})

    ax = axes[0]
    ax.grid(axis="x"); ax.set_axisbelow(True)
    cols = [STATUS["serious"] if k in prods else CAT[0] for k in order]
    ax.barh(range(len(order)), vals, color=cols, height=0.68, zorder=3)
    ax.axvline(BENIGN_MAX, color=AXIS, lw=0.9, ls=(0, (4, 3)), zorder=2)
    ax.axvline(PATH_MIN, color=AXIS, lw=0.9, ls=(0, (4, 3)), zorder=2)
    ax.text(BENIGN_MAX, len(order) - 0.2, " ambiguous ≥ 0.34", fontsize=7, color=MUTED, va="top")
    ax.text(PATH_MIN, len(order) - 0.2, " pathogenic ≥ 0.564", fontsize=7, color=MUTED, va="top")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{k}  (readthrough)" if k in prods else k for k in order], fontsize=8)
    for i, k in enumerate(order):
        ax.text(vals[i] + 0.008, i, f"{vals[i]:.3f}", va="center", fontsize=7,
                color=INK if k in prods else INK2, fontweight="bold" if k in prods else "normal")
    ax.set_xlim(0, max(max(vals) * 1.35, PATH_MIN * 1.1))
    ax.set_xlabel("AlphaMissense pathogenicity")
    ax.set_title(f"Every substitution AlphaMissense scores at residue {e['residue']},\n"
                 f"with the three near-cognate {e['stop_codon']} readthrough products highlighted",
                 loc="left", fontsize=9)
    ax.legend(handles=[Patch(color=STATUS["serious"], label=f"inserted at {e['stop_codon']} (Trp, Arg, Cys)"),
                       Patch(color=CAT[0], label="other substitutions at this residue")],
              loc="lower right")

    ax = axes[1]
    hgrid(ax)
    rm = {k: v for k, v in (e.get("regional_means") or {}).items() if v is not None}
    labels = list(rm) + [f"residue {e['residue']}"]
    heights = list(rm.values()) + [e["residue_tolerance"]["mean_score"]]
    colsr = [NEUTRAL] * len(rm) + [CAT[0]]
    bars = ax.bar(range(len(labels)), heights, color=colsr, edgecolor=AXIS, lw=0.6, width=0.7, zorder=3)
    bars[-1].set_edgecolor(INK)
    for i, h in enumerate(heights):
        ax.text(i, h + 0.012, f"{h:.2f}", ha="center", fontsize=7, color=INK2)
    SHORT = {"Serine/threonine kinase": "kinase", "CDC20 binding": "CDC20", "BUB3 binding": "BUB3",
             "BUB1 N-terminal": "BUB1 N-term", "whole protein": "whole\nprotein"}
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([SHORT.get(l, l) for l in labels], fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("mean AlphaMissense")
    ax.set_ylim(0, max(heights) * 1.3)
    ax.set_title("This residue against its protein", loc="left", fontsize=9)
    fig.tight_layout()
    return save(fig, "readthrough_products")


# ================================================================ 3. ancestry PCA
def somalier_js(html, name):
    i = html.find(f"var {name} = ")
    if i < 0:
        raise ValueError(f"{name} not in somalier html")
    i += len(f"var {name} = ")
    j = html.find("\nvar ", i)
    return json.loads(html[i:j].strip().rstrip(";"))


@figure("ancestry_pca")
def ancestry_pca():
    html = open(path(f"{S}.somalier-ancestry.html")).read()
    bg, q = somalier_js(html, "background_data"), somalier_js(html, "query_data")
    qc = jload("genome_qc.json")
    probs = qc["ancestry"]["probabilities"]
    pops = sorted(bg)
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    qs = list(q.values())[0]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, (px, py) in zip(axes, [(0, 1), (2, 3)]):
        hgrid(ax); ax.grid(axis="x")
        for k, pop in enumerate(pops):
            pcs = bg[pop]["pcs"]
            ax.scatter(pcs[px], pcs[py], s=9, marker=markers[k], color=CAT[k], alpha=0.45, lw=0, zorder=2)
            ax.text(np.median(pcs[px]), np.median(pcs[py]), pop, fontsize=8, color=INK, ha="center", va="center",
                    fontweight="bold", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor="none", alpha=0.75))
        ax.scatter([qs["pcs"][px][0]], [qs["pcs"][py][0]], s=120, marker="*", color=INK, edgecolor=SURFACE,
                   lw=1.2, zorder=5)
        ax.set_xlabel(f"PC{px + 1}"); ax.set_ylabel(f"PC{py + 1}")
    handles = [Line2D([], [], marker=markers[k], color=CAT[k], ls="", ms=6, label=pop) for k, pop in enumerate(pops)]
    handles.append(Line2D([], [], marker="*", color=INK, ls="", ms=11, label="this sample"))
    axes[0].legend(handles=handles, loc="best")
    top = max(probs, key=probs.get)
    fig.suptitle(f"Ancestry: this sample against 1000 Genomes superpopulations (somalier); "
                 f"predicted {top} (p = {probs[top]:.3f}"
                 + "".join(f", {k} {v:.3f}" for k, v in sorted(probs.items(), key=lambda x: -x[1]) if k != top and v > 0)
                 + ")", fontsize=9.5, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout()
    return save(fig, "ancestry_pca")


# ================================================================ 4. stage 1b evidence levels
@figure("stage1b_levels")
def stage1b_levels():
    lv = tload(f"{S}.pgx_ot.evidence_levels.tsv")
    summ = jload(f"{S}.pgx_ot.summary.json")
    asserted = set(summ.get("asserted_levels") or ["1A", "1B"])
    levels = ["1A", "1B", "2A", "2B", "3", "4"]
    rows = sorted(lv, key=lambda r: (r.get("drug_name") or r["drug"]).lower())
    M = np.array([[int(r.get(f"level_{l}") or 0) for l in levels] for r in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 0.42 * len(rows) + 1.4))
    shown = np.log1p(M)
    ax.imshow(shown, cmap=SEQ_CMAP, aspect="auto", vmin=0, vmax=max(shown.max(), 1))
    for i in range(len(rows)):
        for j in range(len(levels)):
            v = int(M[i, j])
            ax.text(j, i, str(v) if v else "·", ha="center", va="center", fontsize=8,
                    color=SURFACE if shown[i, j] > 0.6 * shown.max() else INK)
    ax.set_xticks(range(len(levels))); ax.set_xticklabels([f"level {l}" for l in levels])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([(r.get("drug_name") or r["drug"]).title() for r in rows])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for j, l in enumerate(levels):
        if l in asserted:
            ax.add_patch(Rectangle((j - 0.5, -0.5), 1, len(rows), fill=False, edgecolor=INK, lw=1.4, zorder=5))
    ax.text(len([l for l in levels if l in asserted]) / 2 - 0.5, -0.9, "asserted", ha="center", fontsize=7.5, color=INK)
    ax.text(len(levels) - 1, -0.9, "counted, not joined", ha="right", fontsize=7.5, color=INK2)
    ax.set_title(f"Open Targets (ClinPGx) rows per drug by evidence level, list '{summ['drug_list']['label']}'",
                 loc="left", pad=22)
    fig.tight_layout()
    return save(fig, "stage1b_levels")


# ================================================================ 5. runs of homozygosity
@figure("roh")
def roh():
    segs = []
    for line in open(path(f"{S}.roh.txt")):
        if line.startswith("RG\t"):
            f = line.rstrip("\n").split("\t")
            segs.append((f[2], int(f[3]), int(f[4]), int(f[5]), int(f[6])))
    qc = jload("genome_qc.json")
    r = qc["roh"]
    chroms = [str(i) for i in range(1, 23)]
    fig, ax = plt.subplots(figsize=(11, 3.4))
    hgrid(ax)
    for thr, lbl in [(1.5, "1.5 Mb (F_ROH threshold)"), (5.0, "5 Mb (consanguinity-scale)")]:
        ax.axhline(thr, color=AXIS, lw=0.9, ls=(0, (4, 3)), zorder=1)
        ax.text(21.6, thr * 1.08, lbl, fontsize=7, color=MUTED, ha="right", va="bottom")
    x = [chroms.index(c) for c, *_ in segs if c in chroms]
    y = [l / 1e6 for c, s, e, l, n in segs if c in chroms]
    ax.scatter(x, y, s=22, color=CAT[0], edgecolor=SURFACE, lw=0.6, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(0.03, 20)
    ax.set_xticks(range(22)); ax.set_xticklabels(chroms)
    ax.set_xlabel("chromosome"); ax.set_ylabel("run length (Mb, log)")
    f15 = r["F_ROH"].get(">=1.5Mb", {})
    longest = r.get("longest_roh") or {}
    n_half = sum(1 for *_, l, n in segs if l >= 5e5)
    n_small = sum(1 for *_, l, n in segs if l < 3e4)
    ax.set_title(f"Runs of homozygosity (bcftools roh, autosomes): {len(segs)} runs, {n_half} of them ≥ 0.5 Mb, "
                 f"longest {longest.get('mb', 0):.2f} Mb ({longest.get('region', '')}); "
                 f"F_ROH at ≥ 1.5 Mb = {f15.get('F_ROH', 0):.4f} (n = {f15.get('n', 0)})", loc="left")
    ax.text(0.0, -0.25, f"Reference points: {r.get('interpretation_reference', '')}. Every run is below the 1.5 Mb "
            f"threshold, which is why the report treats consanguinity as excluded. {n_small} runs shorter than "
            f"30 kb fall below the axis.", transform=ax.transAxes, fontsize=7, color=INK2, va="top")
    return save(fig, "roh")


# ================================================================ 6. QC overview (appendix)
@figure("qc_overview")
def qc_overview():
    stats = [l.rstrip("\n").split("\t") for l in open(path(f"{S}.samtools.stats")) if not l.startswith("#")]
    cov = [(int(f[2]), int(f[3])) for f in stats if f[0] == "COV"]
    ins = [(int(f[1]), int(f[2])) for f in stats if f[0] == "IS"]
    qc = jload("genome_qc.json")
    mean = float(qc["mean_depth_autosomes"])
    # per-chromosome mean depth from the aneuploidy bins (kept bins only)
    per = {}
    with gzip.open(path(f"{S}.aneuploidy.bins.tsv.gz"), "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r.get("kept") in ("1", "True", "true") and r.get("depth"):
                per.setdefault(r["chrom"], []).append(float(r["depth"]))
    order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    order = [c for c in order if c in per]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.4), gridspec_kw={"width_ratios": [1.1, 1.6, 1.1]})
    ax = axes[0]; hgrid(ax)
    xs = [c for c, _ in cov if c <= 100]; ys = [b / 1e6 for c, b in cov if c <= 100]
    ax.fill_between(xs, ys, color=CAT[0], alpha=0.25, lw=0); ax.plot(xs, ys, color=CAT[0], lw=1.6)
    ax.axvline(mean, color=INK, lw=0.9, ls=(0, (4, 3))); ax.text(mean + 1.5, max(ys) * 0.95, f"mean {mean:.1f}x", fontsize=7.5, color=INK)
    ax.set_xlabel("depth (x)"); ax.set_ylabel("bases (Mb)"); ax.set_title("Depth distribution", loc="left")
    ax = axes[1]; hgrid(ax)
    vals = [np.mean(per[c]) for c in order]
    cols = [CAT[0] if c not in ("X", "Y") else CAT[1] for c in order]
    ax.bar(range(len(order)), vals, color=cols, width=0.72, zorder=3)
    ax.axhline(mean, color=INK, lw=0.9, ls=(0, (4, 3)))
    ax.axhline(mean / 2, color=AXIS, lw=0.9, ls=(0, (4, 3))); ax.text(-0.6, mean / 2 + 0.8, "half", fontsize=7, color=MUTED, ha="left")
    ax.text(-0.6, mean + 0.8, "autosomal mean", fontsize=7, color=INK, ha="left")
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=7)
    ax.set_ylim(0, max(vals) * 1.22)
    ax.set_ylabel("mean depth (x)")
    sex = (qc.get("sex_inferred") or {}).get("sex", "?")
    ax.set_title(f"Per-chromosome depth (sex inferred {sex}; X {qc.get('mean_depth_X', 0):.1f}x, Y {qc.get('mean_depth_Y', 0):.1f}x)", loc="left")
    ax.legend(handles=[Patch(color=CAT[0], label="autosome"), Patch(color=CAT[1], label="sex chromosome")],
              loc="upper left", ncol=2)
    ax = axes[2]; hgrid(ax)
    xs = [i for i, _ in ins if i <= 1000]; ys = [n / 1e6 for i, n in ins if i <= 1000]
    ax.fill_between(xs, ys, color=CAT[0], alpha=0.25, lw=0); ax.plot(xs, ys, color=CAT[0], lw=1.6)
    mode = xs[int(np.argmax(ys))] if ys else 0
    ax.axvline(mode, color=INK, lw=0.9, ls=(0, (4, 3))); ax.text(mode + 15, max(ys) * 0.95, f"mode {mode} bp", fontsize=7.5, color=INK)
    ax.set_xlabel("insert size (bp)"); ax.set_ylabel("read pairs (M)"); ax.set_title("Insert size", loc="left")
    fig.tight_layout()
    return save(fig, "qc_overview")


# ================================================================ 7. force-call depth (appendix)
@figure("forcecall_depth")
def forcecall_depth():
    dp = {"0/0": [], "0/1": [], "1/1": [], "haploid": []}
    with gzip.open(path(f"{S}.pgx_forcecall.vcf.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fmt, smp = f[8].split(":"), f[9].split(":")
            rec = dict(zip(fmt, smp))
            gt = rec.get("GT", "./.")
            key = gt if gt in dp else ("haploid" if "/" not in gt and "|" not in gt else "0/0")
            if gt in ("./.", "."):
                continue
            dp[key].append(int(rec.get("DP") or 0))
    summ = jload(f"{S}.pgx_forcecall.summary.json")
    keys = [k for k in ["0/0", "0/1", "1/1", "haploid"] if dp[k]]
    hi = max(max(v) for v in dp.values() if v)
    bins = np.arange(0, hi + 5, 5)
    fig, ax = plt.subplots(figsize=(8.5, 3.4))
    hgrid(ax)
    ax.hist([dp[k] for k in keys], bins=bins, stacked=True, color=CAT[:len(keys)],
            label=[f"{k} (n = {len(dp[k])})" for k in keys], edgecolor=SURFACE, lw=0.6, zorder=3)
    ax.axvline(summ["dp_min"], color=INK, lw=0.9, ls=(0, (4, 3)))
    ax.text(summ["dp_min"] + 1, ax.get_ylim()[1] * 0.92, f"min {summ['dp_min']}x", fontsize=7.5, color=INK)
    ax.set_xlabel("read depth at position (x)"); ax.set_ylabel("PharmCAT positions")
    ax.set_title(f"Depth at the {summ['pharmcat_positions']} PharmCAT allele-definition positions "
                 f"(median {summ['dp_median']}x, {summ['missing_gt']} missing genotypes)", loc="left")
    ax.legend(loc="upper right", title="genotype")
    return save(fig, "forcecall_depth")


# ================================================================ 8. MT-RNR1 allele fractions (appendix)
@figure("mt_rnr1")
def mt_rnr1():
    d = jload(f"{S}.mt_rnr1.json")
    alleles = d.get("cpic_alleles") or []
    if not alleles:
        raise FileNotFoundError("mt_rnr1.json (no cpic_alleles)")
    thr = d["thresholds"]
    floor = 1e-4
    fig, ax = plt.subplots(figsize=(9, 0.26 * len(alleles) + 1.8))
    ax.grid(axis="x"); ax.set_axisbelow(True)
    names = [x["allele"] for x in alleles]
    afs = [max(float(x.get("af") or 0), floor) for x in alleles]
    ax.scatter(afs, range(len(alleles)), s=28, color=CAT[0], edgecolor=SURFACE, lw=0.6, zorder=3, label="CPIC allele (tested)")
    ctrl = d.get("positive_controls_observed") or []
    if ctrl:
        ax.scatter([max(float(c["af"]), floor) for c in ctrl], [len(alleles) + 0.6 + i * 0.6 for i in range(len(ctrl))],
                   s=28, marker="D", color=CAT[1], edgecolor=SURFACE, lw=0.6, zorder=3, label="positive control (expected homoplasmic)")
        for i, c in enumerate(ctrl):
            ax.text(floor * 1.15, len(alleles) + 0.6 + i * 0.6, f"m.{c['pos']}{c['ref']}>{c['alt']}", fontsize=7, va="center", color=INK2)
    for key, lbl in [("low_level_min", "low-level"), ("heteroplasmic_min", "heteroplasmic"), ("homoplasmic_min", "homoplasmic")]:
        ax.axvline(thr[key], color=AXIS, lw=0.9, ls=(0, (4, 3)), zorder=1)
        ax.text(thr[key], 1.0, f"{lbl} ≥ {thr[key]}", fontsize=7, color=MUTED, ha="center", va="bottom",
                transform=ax.get_xaxis_transform())
    ax.set_xscale("log"); ax.set_xlim(floor * 0.8, 1.5)
    ax.set_yticks(range(len(alleles))); ax.set_yticklabels(names, fontsize=7)
    ax.set_ylim(-0.8, len(alleles) + 0.6 * max(len(ctrl), 1) + 0.6)
    ax.set_xlabel("alternate allele fraction (log; values at 0 drawn at the floor)")
    ax.set_title(f"MT-RNR1: {d['cpic_alleles_tested']} CPIC alleles at median depth {d['depth_median']}x; "
                 f"called {d.get('pharmcat_outside_call', 'n/a')}", loc="left", pad=16)
    ax.legend(loc="lower right")
    return save(fig, "mt_rnr1")


# ================================================================ 9. panel exon copy number (appendix)
@figure("panel_exon_cn")
def panel_exon_cn():
    ex = tload(f"{S}.panel.exon_cn.tsv")
    genes_tbl = tload(f"{S}.panel.gene_cn.tsv")
    cn_check = jload(f"{S}.panel.cn_candidates.json")
    noise = set(cn_check.get("resolved_as_noise") or [])
    open_c = {c["gene"] for c in cn_check.get("open_candidates") or []}
    order = [g["gene"] for g in sorted(genes_tbl, key=lambda g: (g.get("group") or "", g["gene"]))]
    gi = {g: i for i, g in enumerate(order)}
    gene_cn = {g["gene"]: float(g["gene_cn"]) for g in genes_tbl if g.get("gene_cn")}
    rng = np.random.default_rng(0)
    xs, ys, flagged = [], [], []
    for r in ex:
        if r.get("usable") not in ("True", "true", "1") or r["gene"] not in gi:
            continue
        x = gi[r["gene"]] + rng.uniform(-0.22, 0.22)
        xs.append(x); ys.append(float(r["cn"]))
        if abs(float(r["z"])) >= 3:
            flagged.append((x, float(r["cn"]), r["gene"], r["exon"]))
    fig, ax = plt.subplots(figsize=(13, 3.8))
    hgrid(ax)
    for cn, lbl in [(1, "1 copy"), (2, "2 copies"), (3, "3 copies")]:
        ax.axhline(cn, color=AXIS, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.text(-0.6, cn + 0.04, lbl, fontsize=7, color=MUTED, va="bottom")
    ax.scatter(xs, ys, s=8, color=CAT[0], alpha=0.5, lw=0, zorder=3, label="exon (usable)")
    ax.scatter([gi[g] for g in order if g in gene_cn], [gene_cn[g] for g in order if g in gene_cn], s=34, marker="_",
               color=INK, lw=1.6, zorder=4, label="gene-level CN")
    if flagged:
        ax.scatter([f[0] for f in flagged], [f[1] for f in flagged], s=30, color=STATUS["serious"], edgecolor=SURFACE,
                   lw=0.6, zorder=5, label="exon |z| ≥ 3")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([g + (" †" if g in noise else "") + (" *" if g in open_c else "") for g in order],
                       rotation=90, fontsize=7)
    ax.set_xlim(-0.8, len(order) - 0.2); ax.set_ylim(0, 4)
    ax.set_ylabel("copy number estimate")
    ax.set_title("Exon copy number across the mechanistic panel (mosdepth, GC-corrected; |z| ≥ 3 flagged; "
                 "† run resolved as noise at read level; * open candidate)", loc="left")
    ax.legend(loc="upper right", ncol=3)
    return save(fig, "panel_exon_cn")


json.dump(manifest, open(os.path.join(a.outdir, f"{S}.figures.json"), "w"), indent=2)
print(json.dumps(manifest, indent=1))
