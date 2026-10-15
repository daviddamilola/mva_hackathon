#!/usr/bin/env python3
"""Map this child's coding variants onto the protein's domain organisation.

Domain coordinates come from a literature asset with its citation in the header; variant positions
come from the pipeline's own rare-and-functional panel output, parsed from HGVSp. Nothing about the
variants is hardcoded: if the panel output changes, the figure changes.

For a stop-gained variant the figure also shows what is left of the protein (residues 1..pos-1) and
which domains that truncation keeps or removes, which is the point the report makes in words.

Usage: domain_figure.py --domains bubr1_domains.tsv --rare-functional rf.tsv --gene BUB1B
                        --out-png fig.png --out-json fig.json
"""
import argparse
import csv
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

ap = argparse.ArgumentParser()
for x in ("domains", "rare-functional", "gene", "out-png", "out-json"):
    ap.add_argument(f"--{x}", required=True)
a = ap.parse_args()

AA3 = ("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter Sec "
       "Xaa").split()
HGVSP = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*|=|fs.*)?")

domain_source = next((l.split("Source:", 1)[1].strip() for l in open(a.domains)
                      if l.startswith("# Source")), "source not recorded")

rows = [r for r in csv.DictReader((l for l in open(a.domains) if not l.startswith("#")), delimiter="\t")
        if r["protein"]]
if not rows:
    raise SystemExit(f"{a.domains}: no domain rows")
length = int(rows[0]["length"])
protein = rows[0]["protein"]
domains = [{"domain": r["domain"], "start": int(r["start"]), "end": int(r["end"]), "role": r["role"]}
           for r in rows]

variants = []
for r in csv.DictReader(open(a.rare_functional), delimiter="\t"):
    if r["gene"] != a.gene or not r.get("hgvsp"):
        continue
    m = HGVSP.search(r["hgvsp"])
    if not m:
        continue
    ref_aa, pos, alt_aa = m.group(1), int(m.group(2)), m.group(3) or ""
    if ref_aa not in AA3:
        continue
    truncating = "stop_gained" in r["consequence"] or alt_aa in ("Ter", "*") or alt_aa.startswith("fs")
    variants.append({"hgvsp": r["hgvsp"].split(":")[-1], "position": pos, "consequence": r["consequence"],
                     "zygosity": r.get("zygosity"), "truncating": truncating,
                     "in_domains": [d["domain"] for d in domains if d["start"] <= pos <= d["end"]]})
variants.sort(key=lambda v: v["position"])


def domains_kept(pos):
    """Domains fully present in a protein truncated immediately before `pos`."""
    return [d["domain"] for d in domains if d["end"] < pos]


for v in variants:
    if v["truncating"]:
        v["retained_domains"] = domains_kept(v["position"])
        v["lost_domains"] = [d["domain"] for d in domains if d["domain"] not in v["retained_domains"]]

PALETTE = ["#4a6fa5", "#2e8b57", "#b07d2b", "#7d5ba6"]
fig, ax = plt.subplots(figsize=(11, 2.9))
ax.add_patch(Rectangle((1, 0.44), length - 1, 0.12, facecolor="#e6e6e6", edgecolor="#999", lw=0.6))
for i, d in enumerate(domains):
    ax.add_patch(Rectangle((d["start"], 0.36), d["end"] - d["start"], 0.28,
                           facecolor=PALETTE[i % len(PALETTE)], edgecolor="none", alpha=0.85))
    ax.text((d["start"] + d["end"]) / 2, 0.50, d["domain"], ha="center", va="center", fontsize=7,
            color="white", fontweight="bold")
    ax.text((d["start"] + d["end"]) / 2, 0.28, f"{d['start']}-{d['end']}", ha="center", va="center", fontsize=6,
            color="#444")

for i, v in enumerate(variants):
    y = 0.78 + 0.13 * (i % 2)
    colour = "#c0392b" if v["truncating"] else "#1f6f8b"
    ax.plot([v["position"], v["position"]], [0.62, y - 0.04], color=colour, lw=1.1)
    ax.plot([v["position"]], [0.62], marker="v", color=colour, markersize=6)
    ax.text(v["position"], y, v["hgvsp"], ha="center", va="bottom", fontsize=7.5, color=colour)
    if v["truncating"]:
        ax.add_patch(Rectangle((v["position"], 0.44), length - v["position"], 0.12,
                               facecolor="#c0392b", edgecolor="none", alpha=0.18))
        ax.text((v["position"] + length) / 2, 0.15, "not translated on this allele", ha="center",
                va="center", fontsize=6.5, color="#c0392b", style="italic")

ax.set_xlim(-20, length + 20)
ax.set_ylim(0.05, 1.02)
ax.set_yticks([])
ax.set_xlabel(f"{protein} residue")
for sp in ("top", "left", "right"):
    ax.spines[sp].set_visible(False)
ax.set_title(f"{protein} domain organisation and this child's {a.gene} variants", fontsize=9)
fig.tight_layout()
fig.savefig(a.out_png, dpi=200)
plt.close(fig)

json.dump({"protein": protein, "length": length, "gene": a.gene, "domains": domains,
           "variants": variants, "domain_source": domain_source},
          open(a.out_json, "w"), indent=2)
print(json.dumps({"variants": variants}, indent=1))
