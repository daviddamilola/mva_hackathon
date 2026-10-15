#!/usr/bin/env python3
"""Modifier-panel variant tables from a VEP-annotated VCF (CSQ in INFO).

For each variant: the child's genotype (GT, AD, DP, FILTER), the consequence on the MANE Select
transcript (else the canonical, else the most severe), HGVS, gnomAD v4.1 joint AF / grpmax AF /
homozygote count / per-group AFs, ClinVar significance and review status, SpliceAI delta scores
(max), and AlphaMissense score and class. Panel group and gene come from the panel BED.

Filters, reported and never hidden:
  rare        gnomAD AF_joint < --max-af, or absent from gnomAD
  functional  non-synonymous or splice-affecting: a VEP consequence in FUNCTIONAL (MANE/selected
              transcript), OR SpliceAI max delta >= --spliceai-min
Everything found is written to *.all.tsv; rare & functional to *.rare_functional.tsv.
Classification is not attempted: this is hypothesis generation, not diagnostic interpretation.

Usage: modifier_panel_report.py --vcf vep.vcf.gz --bed panel.bed --sample S --prefix out [--max-af 0.01] [--spliceai-min 0.2]
"""
import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict

FUNCTIONAL = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
    "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", "feature_truncation",
    "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
    "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant",
}
SEVERITY = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
            "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", "feature_truncation",
            "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
            "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
            "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant",
            "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
            "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant",
            "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
            "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
            "regulatory_region_variant", "intergenic_variant", "sequence_variant"]
RANK = {c: i for i, c in enumerate(SEVERITY)}

ap = argparse.ArgumentParser()
ap.add_argument("--vcf", required=True)
ap.add_argument("--bed", required=True)
ap.add_argument("--sample", required=True)
ap.add_argument("--prefix", required=True)
ap.add_argument("--max-af", type=float, default=0.01)
ap.add_argument("--spliceai-min", type=float, default=0.2)
a = ap.parse_args()

panel = []
for line in open(a.bed):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    panel.append((f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]))


def panel_gene(chrom, pos):
    hits = [p for p in panel if p[0] == chrom and p[1] < pos <= p[2]]
    return hits


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


csq_fields, rows = None, []
with gzip.open(a.vcf, "rt") as fh:
    for line in fh:
        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, filt = f[0], int(f[1]), f[3], f[4], f[6]
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
        gt = fmt.get("GT", "./.")
        idx = [x for x in re.split(r"[/|]", gt) if x not in (".", "")]
        zyg = ("hom_alt" if len(set(idx)) == 1 and idx[0] != "0" else "het" if len(set(idx)) > 1 and "0" in idx
               else "hom_ref" if idx and set(idx) == {"0"} else "multi_alt_het" if idx else "no_call")
        csqs = [dict(zip(csq_fields, c.split("|"))) for c in info.get("CSQ", "").split(",") if c] if csq_fields else []
        in_panel = panel_gene(chrom, pos)
        panel_symbols = {p[3] for p in in_panel}

        def pick(cs):
            for key in ("MANE_SELECT", "CANONICAL"):
                sel = [c for c in cs if c.get(key) not in (None, "", "-") and (not panel_symbols or c.get("SYMBOL") in panel_symbols)]
                if sel:
                    return sel[0], key.lower()
            if cs:
                return min(cs, key=lambda c: min(RANK.get(x, 99) for x in c.get("Consequence", "").split("&"))), "most_severe"
            return {}, "none"

        chosen, basis = pick(csqs)
        cons = chosen.get("Consequence", "")
        # SpliceAI custom annotation lives in CSQ key 'SpliceAI_SpliceAI' (key 'SpliceAI' alone is the
        # custom file's ID column). VEP replaces the '|' inside the value with '&', so each entry is
        # 10 '&'-separated fields: ALLELE, SYMBOL, DS_AG, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL;
        # several entries (one per overlapping gene) are concatenated.
        sai_raw = chosen.get("SpliceAI_SpliceAI", "") or ""
        toks = [t for t in re.split(r"[&|,]", sai_raw) if t != ""]
        entries = [toks[i:i + 10] for i in range(0, len(toks) - len(toks) % 10, 10)]
        sai_scores = [fnum(x) for e in entries for x in e[2:6]]
        sai_max = max([x for x in sai_scores if x is not None], default=None)
        sai_raw = " ; ".join("|".join(e) for e in entries)
        af = fnum(chosen.get("gnomAD_joint_AF_joint"))
        functional_csq = any(c in FUNCTIONAL for c in cons.split("&"))
        functional_splice = sai_max is not None and sai_max >= a.spliceai_min
        rows.append({
            "gene": ";".join(sorted(panel_symbols)) or chosen.get("SYMBOL", ""),
            "panel_group": ";".join(sorted({p[5] for p in in_panel})),
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "filter": filt,
            "genotype": gt, "zygosity": zyg, "AD": fmt.get("AD", ""), "DP": fmt.get("DP", ""), "GQ": fmt.get("GQ", ""),
            "consequence": cons, "impact": chosen.get("IMPACT", ""), "transcript": chosen.get("Feature", ""),
            "transcript_basis": basis, "mane_select": chosen.get("MANE_SELECT", ""),
            "hgvsc": chosen.get("HGVSc", ""), "hgvsp": chosen.get("HGVSp", ""), "existing_variation": chosen.get("Existing_variation", ""),
            "gnomad_af_joint": "" if af is None else af,
            "gnomad_af_grpmax_joint": chosen.get("gnomAD_joint_AF_grpmax_joint", ""),
            "gnomad_grpmax_joint": chosen.get("gnomAD_joint_grpmax_joint", ""),
            "gnomad_nhomalt_joint": chosen.get("gnomAD_joint_nhomalt_joint", ""),
            "gnomad_faf95_max_joint": chosen.get("gnomAD_joint_fafmax_faf95_max_joint", ""),
            "gnomad_filters": ";".join(x for x in (chosen.get("gnomAD_joint_exomes_filters", ""), chosen.get("gnomAD_joint_genomes_filters", "")) if x),
            "clinvar_significance": chosen.get("ClinVar_CLNSIG", ""), "clinvar_review": chosen.get("ClinVar_CLNREVSTAT", ""),
            "clinvar_condition": chosen.get("ClinVar_CLNDN", ""), "clinvar_id": chosen.get("ClinVar", ""),
            "spliceai_max_delta": "" if sai_max is None else sai_max, "spliceai_raw": sai_raw,
            "alphamissense_score": chosen.get("AlphaMissense_AM_PATHOGENICITY", ""), "alphamissense_class": chosen.get("AlphaMissense_AM_CLASS", ""),
            "rare": af is None or af < a.max_af, "in_gnomad": af is not None,
            "functional": functional_csq or functional_splice,
            "functional_basis": ";".join(x for x, ok in (("consequence", functional_csq), (f"spliceai>={a.spliceai_min}", functional_splice)) if ok),
        })

fields = list(rows[0].keys()) if rows else ["gene"]
for suffix, subset in (("all", rows), ("rare_functional", [r for r in rows if r["rare"] and r["functional"] and r["zygosity"] not in ("hom_ref", "no_call")])):
    with open(f"{a.prefix}.{suffix}.tsv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(subset)

rf = [r for r in rows if r["rare"] and r["functional"] and r["zygosity"] not in ("hom_ref", "no_call")]
summary = {
    "sample": a.sample, "max_af": a.max_af, "spliceai_min": a.spliceai_min,
    "panel_genes": len({p[3] for p in panel}), "panel_bp": sum(p[2] - p[1] for p in panel),
    "variants_total": len(rows), "variants_pass": sum(1 for r in rows if r["filter"] == "PASS"),
    "rare": sum(1 for r in rows if r["rare"]), "not_in_gnomad": sum(1 for r in rows if not r["in_gnomad"]),
    "functional": sum(1 for r in rows if r["functional"]),
    "rare_functional": len(rf),
    "rare_functional_by_group": dict(Counter(r["panel_group"] for r in rf)),
    "rare_functional_genes": sorted({r["gene"] for r in rf}),
    "clinvar_annotated_rare_functional": sum(1 for r in rf if r["clinvar_significance"]),
    "variants_with_mane_consequence": sum(1 for r in rows if r["transcript_basis"] == "mane_select"),
}
json.dump(summary, open(f"{a.prefix}.summary.json", "w"), indent=2)
print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, (list, dict))}))
