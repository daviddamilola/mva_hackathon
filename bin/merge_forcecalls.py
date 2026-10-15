#!/usr/bin/env python3
"""Merge constrained and unconstrained force-calls at PharmCAT positions, and gate on evidence.

MINIMUM EVIDENCE (why this exists)
  Constrained calling (-C alleles -A) emits a genotype at every definition position by
  construction, so "no position was left uncalled" is a property of the method, not of the data,
  and cannot be used to argue that hom-ref calls are safe. The clinically dangerous direction is a
  false reference call producing a spurious *1.
  So every position is tested against a stated threshold: depth >= min_dp and GQ >= min_gq, where
  GQ is derived from PL (second-smallest minus smallest) because bcftools call does not emit GQ
  here. A position failing either test is set to ./. and left for PharmCAT to handle explicitly;
  the missing-to-reference options are deliberately NOT used, so a no-call stays a no-call.
  Per-position evidence is written out so the claim can be checked rather than taken on trust.
  Thresholds are applied to the merged record, after the unconstrained replacement below.

Constrained calls (bcftools call -C alleles -A) genotype every PharmCAT definition allele
explicitly, including hom-ref at indel sites. However, they can only test the alleles PharmCAT
defines: a real non-reference allele that PharmCAT does not define is written as 0/0.

The unconstrained call (bcftools call -m) sees the observed allele. At any position where the
unconstrained genotype carries a non-reference allele absent from the constrained ALT list, the
unconstrained record replaces the constrained one. PharmCAT then sees the real allele and raises
its own undocumented-variation warning instead of reporting reference. Every such site is written
to off_definition.tsv.

Usage: merge_forcecalls.py <constrained.vcf.gz> <unconstrained.vcf.gz> <merged.vcf>
                           <off_definition.tsv> <evidence.tsv> [min_dp] [min_gq]
"""
import gzip
import json
import sys


def read_vcf(path):
    header, records = [], {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            records.setdefault((f[0], int(f[1])), []).append(f)
    return header, records


def gt_alleles(f):
    fmt = f[8].split(":")
    gt = dict(zip(fmt, f[9].split(":"))).get("GT", ".")
    return [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]


def field(f, key):
    return dict(zip(f[8].split(":"), f[9].split(":"))).get(key, ".")


con_path, unc_path, out_vcf, out_tsv, out_evidence = sys.argv[1:6]
MIN_DP = int(sys.argv[6]) if len(sys.argv) > 6 else 10
MIN_GQ = int(sys.argv[7]) if len(sys.argv) > 7 else 20


def evidence(f):
    """Depth, alt fraction and GQ (derived from PL) for the sample column."""
    d = dict(zip(f[8].split(":"), f[9].split(":")))
    ad = [int(x) for x in d.get("AD", "").split(",") if x not in (".", "")]
    pl = sorted(int(x) for x in d.get("PL", "").split(",") if x not in (".", ""))
    dp_f = d.get("DP", ".")
    tot = sum(ad) if ad else (int(dp_f) if dp_f.isdigit() else 0)
    alt = sum(ad[1:]) if len(ad) > 1 else 0
    return {"gt": d.get("GT", "."), "depth": tot,
            "alt_reads": alt, "af": round(alt / tot, 4) if tot else None,
            "gq": (pl[1] - pl[0]) if len(pl) > 1 else None,
            "haploid": "/" not in d.get("GT", "") and "|" not in d.get("GT", "")}


def no_call(f):
    """Blank the genotype, keeping the evidence fields so the reason stays inspectable."""
    keys = f[8].split(":")
    vals = f[9].split(":")
    d = dict(zip(keys, vals))
    d["GT"] = "./." if "/" in d.get("GT", "") else "."
    f = list(f)
    f[9] = ":".join(d[k] for k in keys)
    return f


con_header, con = read_vcf(con_path)
_, unc = read_vcf(unc_path)

off_rows, ev_rows = [], []
with open(out_vcf, "w") as out:
    out.writelines(con_header)
    for key in con:  # insertion order = input order, already sorted by header contig order
        for cf in con[key]:
            replacement = None
            con_alts = set(cf[4].split(",")) - {"."}
            for uf in unc.get(key, []):
                alleles = gt_alleles(uf)
                if not alleles or all(a == "0" for a in alleles):
                    continue
                unc_alts = uf[4].split(",")
                observed = {unc_alts[int(a) - 1] for a in alleles if a != "0"}
                # compare on the allele difference relative to REF, so representation differences
                # (same event, different REF padding) are not flagged
                def event(ref, alt):
                    return (ref, alt) if len(ref) == len(alt) == 1 else (ref[1:], alt[1:], len(ref) - len(alt))
                con_events = {event(cf[3], a) for a in con_alts}
                novel = [a for a in observed if event(uf[3], a) not in con_events]
                if novel:
                    replacement = uf
                    off_rows.append([uf[0], uf[1], cf[2], cf[3], cf[4], uf[3], ",".join(novel),
                                     field(uf, "GT"), field(uf, "AD"), field(uf, "DP"), field(cf, "GT")])
                    break
            rec = list(replacement or cf)
            ev = evidence(rec)
            fails = ([] if ev["depth"] >= MIN_DP else [f"depth {ev['depth']} < {MIN_DP}"]) + \
                    ([] if (ev["gq"] is None or ev["gq"] >= MIN_GQ) else [f"GQ {ev['gq']} < {MIN_GQ}"])
            if fails:
                rec = no_call(rec)
            ev_rows.append({"chrom": rec[0], "pos": rec[1], "id": rec[2], "ref": rec[3],
                            "alt_allele": rec[4], **ev,
                            "no_called": bool(fails), "reason": "; ".join(fails)})
            out.write("\t".join(rec) + "\n")

with open(out_tsv, "w") as out:
    out.write("chrom\tpos\tid\tdefinition_ref\tdefinition_alts\tobserved_ref\tobserved_undefined_alt\t"
              "observed_gt\tobserved_ad\tobserved_dp\tconstrained_gt\n")
    for r in off_rows:
        out.write("\t".join(r) + "\n")
cols = ["chrom", "pos", "id", "ref", "alt_allele", "gt", "depth", "alt_reads", "af", "gq",
        "haploid", "no_called", "reason"]
with open(out_evidence, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in ev_rows:
        out.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

dip = [r for r in ev_rows if not r["haploid"]]
hap = [r for r in ev_rows if r["haploid"]]
gqs = sorted(r["gq"] for r in ev_rows if r["gq"] is not None)
hets = [r for r in ev_rows if r["gt"] in ("0/1", "1/0")]
homref_alt = [r for r in ev_rows if r["gt"] in ("0/0", "0") and r["alt_reads"] > 0]
pct = lambda xs, q: xs[min(len(xs) - 1, max(0, int(len(xs) * q)))] if xs else None
summary = {
    "thresholds": {"min_depth": MIN_DP, "min_gq": MIN_GQ,
                   "note": "GQ derived from PL; a failing position is set to no-call and "
                           "PharmCAT's missing-to-reference options are not used"},
    "positions": len(ev_rows),
    "no_called": sum(1 for r in ev_rows if r["no_called"]),
    "no_call_reasons": [f"{r['chrom']}:{r['pos']} {r['reason']}" for r in ev_rows if r["no_called"]],
    "depth_diploid": {"min": min((r["depth"] for r in dip), default=None),
                      "p05": pct([r["depth"] for r in sorted(dip, key=lambda x: x["depth"])], 0.05),
                      "median": pct([r["depth"] for r in sorted(dip, key=lambda x: x["depth"])], 0.5)},
    "depth_haploid": {"n": len(hap), "min": min((r["depth"] for r in hap), default=None)},
    "gq": {"min": gqs[0] if gqs else None, "p05": pct(gqs, 0.05),
           "median": pct(gqs, 0.5), "below_min_gq": sum(1 for g in gqs if g < MIN_GQ)},
    "het_allele_balance": {"n": len(hets),
                           "min": min((r["af"] for r in hets), default=None),
                           "max": max((r["af"] for r in hets), default=None)},
    "hom_ref_with_alt_reads": {"n": len(homref_alt),
                               "max_alt_fraction": max((r["af"] for r in homref_alt), default=0.0)},
}
json.dump(summary, open(out_evidence.replace(".tsv", ".summary.json"), "w"), indent=2)
print(f"{len(off_rows)} PharmCAT position(s) carry an allele outside the PharmCAT definition; "
      f"{summary['no_called']} position(s) below the evidence threshold")
