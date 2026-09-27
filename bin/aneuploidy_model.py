#!/usr/bin/env python3
"""Mosaic aneuploidy burden per chromosome arm from WGS depth (LRR) and allelic imbalance (BAF).

Implemented directly (MoChA-style logic, single sample, no panel of normals), with every assumption
stated so the report can repeat it.

DEPTH (LRR)
  mosdepth mean depth in fixed bins (MAPQ >= 20). Bins kept: mappable fraction (Umap k100) >= min_map,
  no blacklist overlap, GC within [0.25, 0.75], depth > 0. Normalised to the autosomal median, then GC-corrected by
  dividing by the median ratio of bins in the same GC stratum (0.5% strata, >= 30 bins, neighbouring
  strata pooled when sparse). LRR = log2(corrected ratio / autosomal median). chrX non-PAR in an XY
  sample is expected at LRR -1 and is reported as an internal control.

ALLELIC IMBALANCE (BAF)
  Heterozygous PASS biallelic SNVs from the VCF (GT 0/1, allele depths), in kept bins, outside
  segmental duplications (UCSC genomicSuperDups: paralogous sequence variants there masquerade as
  skewed hets), with depth within [0.5, 2] x the median het depth; chrX non-PAR excluded in XY
  samples. Alt-allele counts are modelled as (1 - eps) x [50:50 mixture of beta-binomials centred at
  m +/- d] + eps x uniform, where m is the genome-wide mean alt fraction (reference bias), d the
  per-arm imbalance and eps an outlier fraction profiled out per arm over {0, .005, .01, .02, .05}
  (so clean arms are unbiased and a few extreme sites cannot drive d). Overdispersion rho is estimated per arm
  under d = 0 and summarised by the autosomal median (robust to a minority of aneuploid arms). d is
  estimated by maximum likelihood (grid), tested against d = 0 by LRT (boundary: p = 0.5 * P(chi2_1 > T)),
  with a 95% profile-likelihood CI. Arms are tested with Bonferroni correction.
  SPATIAL CONSISTENCY (window-deletion test): an arm-level mosaic event must be present along the
  arm. The arm is split into 5 Mb windows (>= 200 sites) and each window gets its own LRT. The
  ceil(20%) most imbalanced windows (at least one) are removed and the rest of the arm is re-tested.
  The arm is called only if it stays significant at the Bonferroni level. A genuine arm-wide event
  survives losing a fifth of its windows; a signal driven by one hotspot (paralogous sequence
  variants, a residual artefact or a focal event) does not. Such arms are reported as 'focal
  imbalance' with the responsible windows, and get no cell-fraction estimate.

CELL FRACTION
  For a clone fraction f:  gain (3 copies)  d = f / (2 (2 + f)),  LRR = log2(1 + f/2)
                           loss (1 copy)    d = f / (2 (2 - f)),  LRR = log2(1 - f/2)
                           CN-LOH           d = f / 2,            LRR = 0
  Each model converts d (and its CI) to f. The reported model is the one whose predicted LRR is closest
  to the observed arm LRR in SE units; the alternatives are reported too. The LRR-only estimate of f is
  also given. Detection limit: smallest f detectable with 80% power at the Bonferroni alpha, from an
  analytic variance-test approximation using the arm's site count, depth distribution, m and rho.

Usage: aneuploidy_model.py --bins bins.tsv --vcf sample.vcf.gz --cytoband cytoBand.nochr.txt --sex M --prefix out
  bins.tsv columns: chrom start end gc map_frac blacklist_bp depth
  --segdup: UCSC genomicSuperDups (bin chrom start end ...) or BED
"""
import argparse
import gzip
import json
import math
import sys

import numpy as np
from scipy import stats
from scipy.special import betaln, gammaln

ap = argparse.ArgumentParser()
ap.add_argument("--bins", required=True)
ap.add_argument("--vcf", required=True)
ap.add_argument("--cytoband", required=True)
ap.add_argument("--segdup", required=True)
ap.add_argument("--outlier-eps", default="0,0.005,0.01,0.02,0.05", help="outlier fractions profiled over")
ap.add_argument("--window", type=int, default=5_000_000)
ap.add_argument("--sex", choices=["M", "F"], required=True)
ap.add_argument("--prefix", required=True)
ap.add_argument("--min-map", type=float, default=0.9)
ap.add_argument("--min-bins-arm", type=int, default=20)
ap.add_argument("--min-sites-arm", type=int, default=200)
ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args()

AUTOS = [str(i) for i in range(1, 23)]
CHROMS = AUTOS + ["X", "Y"]
PAR = {"X": [(10001, 2781479), (155701383, 156030895)]}
ACRO_P = {"13", "14", "15", "21", "22"}


def in_par(c, p):
    return any(s <= p <= e for s, e in PAR.get(c, []))


# ------------------------------------------------------------------ arms from cytoband
cen, chrom_len = {}, {}
for line in open(a.cytoband):
    c, s, e, band, stain = line.rstrip("\n").split("\t")[:5]
    s, e = int(s), int(e)
    chrom_len[c] = max(chrom_len.get(c, 0), e)
    if stain == "acen":
        lo, hi = cen.get(c, (s, e))
        cen[c] = (min(lo, s), max(hi, e))


def arm_of(c, pos):
    if c not in cen:
        return None
    lo, hi = cen[c]
    return f"{c}p" if pos < lo else (f"{c}q" if pos >= hi else None)


ARMS = [f"{c}{x}" for c in CHROMS for x in ("p", "q") if c in cen]

# ------------------------------------------------------------------ segmental duplications (merged)
sd_raw = {}
for line in open(a.segdup):
    f = line.rstrip("\n").split("\t")
    if len(f) < 3 or line.startswith("#"):
        continue
    c, s_, e_ = (f[1], f[2], f[3]) if (f[0].isdigit() and len(f) > 4 and not f[1].isdigit()) else (f[0], f[1], f[2])
    c = c.replace("chr", "")
    sd_raw.setdefault(c, []).append((int(s_), int(e_)))
segdup = {}
for c, iv in sd_raw.items():
    merged = []
    for s_, e_ in sorted(iv):
        if merged and s_ <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e_))
        else:
            merged.append((s_, e_))
    segdup[c] = (np.array([m[0] for m in merged]), np.array([m[1] for m in merged]))


def in_segdup(c, pos):
    if c not in segdup:
        return False
    st, en = segdup[c]
    i = np.searchsorted(st, pos, side="right") - 1
    return i >= 0 and pos < en[i]

# ------------------------------------------------------------------ bins -> LRR
B = []
for line in open(a.bins):
    if line.startswith("#") or line.startswith("chrom"):
        continue
    c, s, e, gc, mp, bl, dp = line.rstrip("\n").split("\t")[:7]
    if c not in CHROMS:
        continue
    B.append((c, int(s), int(e), float(gc), float(mp), int(float(bl)), float(dp)))
chrom = np.array([b[0] for b in B]); start = np.array([b[1] for b in B]); end = np.array([b[2] for b in B])
gc = np.array([b[3] for b in B]); mp = np.array([b[4] for b in B]); bl = np.array([b[5] for b in B]); dp = np.array([b[6] for b in B])
mid = (start + end) // 2
arm = np.array([arm_of(c, m) or "" for c, m in zip(chrom, mid)])
keep = (mp >= a.min_map) & (bl == 0) & (dp > 0) & (gc >= 0.25) & (gc <= 0.75) & (arm != "")
auto = np.isin(chrom, AUTOS)
ratio = dp / np.median(dp[keep & auto])

strata = np.round(gc / 0.005).astype(int)
corr = np.full(len(B), np.nan)
ref_mask = keep & auto
for s in np.unique(strata[keep]):
    w, sel = 0, None
    while True:  # pool neighbouring strata until >= 30 autosomal reference bins
        sel = ref_mask & (np.abs(strata - s) <= w)
        if sel.sum() >= 30 or w > 20:
            break
        w += 1
    corr[keep & (strata == s)] = ratio[keep & (strata == s)] / np.median(ratio[sel])
lrr = np.full(len(B), np.nan)
lrr[keep] = np.log2(corr[keep] / np.median(corr[keep & auto]))
lrr_sd_auto = float(np.nanstd(lrr[keep & auto]))

# ------------------------------------------------------------------ het sites -> BAF
bin_index = {}
for i, (c, s) in enumerate(zip(chrom, start)):
    bin_index[(c, s)] = i
binsize = int(np.median(end - start))
sites_c, sites_p, sites_k, sites_n = [], [], [], []
n_segdup_excluded = [0]
with gzip.open(a.vcf, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.split("\t", 10)
        c, p, ref, alt, filt = f[0], int(f[1]), f[3], f[4], f[6]
        if filt != "PASS" or len(ref) != 1 or len(alt) != 1 or c not in CHROMS:
            continue
        if a.sex == "M" and c in ("X", "Y") and not in_par(c, p):
            continue
        fmt = f[8].split(":"); val = f[9].rstrip("\n").split(":")
        d = dict(zip(fmt, val))
        if d.get("GT", "").replace("|", "/") != "0/1" or "AD" not in d:
            continue
        ad = d["AD"].split(",")
        if len(ad) != 2 or "." in ad:
            continue
        r_, k_ = int(ad[0]), int(ad[1])
        i = bin_index.get((c, (p - 1) // binsize * binsize))
        if i is None or not keep[i] or in_segdup(c, p - 1):
            n_segdup_excluded[0] += in_segdup(c, p - 1) and i is not None and bool(keep[i])
            continue
        sites_c.append(c); sites_p.append(p); sites_k.append(k_); sites_n.append(r_ + k_)
sc = np.array(sites_c); sp = np.array(sites_p); sk = np.array(sites_k); sn = np.array(sites_n)
med_het_dp = float(np.median(sn[np.isin(sc, AUTOS)]))
okd = (sn >= 0.5 * med_het_dp) & (sn <= 2 * med_het_dp)
sc, sp, sk, sn = sc[okd], sp[okd], sk[okd], sn[okd]
sarm = np.array([arm_of(c, p) or "" for c, p in zip(sc, sp)])
baf = sk / sn
m_bias = float(np.mean(baf[np.isin(sc, AUTOS)]))


def bb_logpmf(k, n, mu, rho):
    # beta-binomial with mean mu and intra-class correlation rho
    ab = (1 - rho) / rho
    al, be = mu * ab, (1 - mu) * ab
    # the binomial coefficient is constant in d and omitted (it cancels in every likelihood ratio)
    return betaln(k + al, n - k + be) - betaln(al, be)


def rho_mom(k, n, mu):
    # method of moments for overdispersion under d = 0
    x = k / n
    s2 = np.mean((x - mu) ** 2)
    base = np.mean(mu * (1 - mu) / n)
    nbar = np.mean(n)
    return float(np.clip((s2 / base - 1) / (nbar - 1), 1e-5, 0.2))


rho_by_arm = {A: rho_mom(sk[sarm == A], sn[sarm == A], m_bias)
              for A in ARMS if (sarm == A).sum() >= a.min_sites_arm and A[:-1] in AUTOS}
rho = float(np.median(list(rho_by_arm.values())))
GRID = np.arange(0.0, 0.30001, 0.001)
EPS = [float(x) for x in a.outlier_eps.split(",")]


def arm_likelihood(k, n, grid=None):
    grid = GRID if grid is None else grid
    ll = np.full(len(grid), -np.inf)
    # outlier component: uniform over k in 0..n; the omitted binomial coefficient is added back to
    # the beta-binomial terms so both components are on the same (full pmf) scale. The outlier
    # fraction eps is profiled out: ll(d) = max over eps.
    lcoef = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    lunif = -np.log(n + 1)
    for j, d in enumerate(grid):
        lo, hi = np.clip(m_bias - d, 1e-4, 1 - 1e-4), np.clip(m_bias + d, 1e-4, 1 - 1e-4)
        lmix = np.logaddexp(bb_logpmf(k, n, lo, rho), bb_logpmf(k, n, hi, rho)) + math.log(0.5) + lcoef
        for eps in EPS:
            v = np.sum(lmix) if eps == 0 else np.sum(np.logaddexp(math.log(1 - eps) + lmix, math.log(eps) + lunif))
            ll[j] = max(ll[j], v)
    return ll


WGRID = np.arange(0.0, 0.15001, 0.0025)


def window_profile(A):
    """d estimate per window along the arm (>= 200 sites)."""
    idx = np.where(sarm == A)[0]
    out = []
    for w in np.unique(sp[idx] // a.window):
        sel = idx[(sp[idx] // a.window) == w]
        if len(sel) < 200:
            continue
        ll = arm_likelihood(sk[sel], sn[sel], WGRID)
        j = int(np.argmax(ll))
        T = max(0.0, 2 * (ll[j] - ll[0]))
        out.append({"start": int(w * a.window), "end": int((w + 1) * a.window), "n_sites": int(len(sel)),
                    "d": float(WGRID[j]), "p": float(0.5 * stats.chi2.sf(T, 1)) if T > 0 else 1.0})
    return out


def d_to_f(d, model):
    if model == "gain":
        return 4 * d / (1 - 2 * d) if d < 0.5 else float("nan")
    if model == "loss":
        return 4 * d / (1 + 2 * d)
    return 2 * d


def f_to_lrr(f, model):
    if model == "gain":
        return math.log2(1 + f / 2)
    if model == "loss":
        return math.log2(max(1e-6, 1 - f / 2))
    return 0.0


tested = [A for A in ARMS if (A[:-1] in AUTOS or a.sex == "F") and (sarm == A).sum() >= a.min_sites_arm
          and ((arm == A) & keep).sum() >= a.min_bins_arm]
alpha_b = a.alpha / max(1, len(tested))
z_a, z_b = stats.norm.isf(alpha_b), stats.norm.isf(0.2)

rows = []
for A in ARMS:
    bm = (arm == A) & keep
    k, n = sk[sarm == A], sn[sarm == A]
    r = {"arm": A, "chrom": A[:-1], "n_bins": int(bm.sum()), "n_het_sites": int(len(k)),
         "tested": A in tested}
    if bm.sum() > 0:
        v = lrr[bm]
        r.update(lrr_mean=float(np.mean(v)), lrr_se=float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else None)
    if A in tested:
        ll = arm_likelihood(k, n)
        j = int(np.argmax(ll)); d_hat = float(GRID[j])
        T = max(0.0, 2 * (ll[j] - ll[0]))
        p = 0.5 * stats.chi2.sf(T, 1) if T > 0 else 1.0
        inside = GRID[ll >= ll[j] - 1.92]
        d_lo, d_hi = float(inside.min()), float(inside.max())
        sig2 = m_bias * (1 - m_bias) / n * (1 + (n - 1) * rho)
        d_min = math.sqrt((z_a + z_b) * math.sqrt(2 * np.mean(sig2 ** 2) / len(n)))
        models = {}
        for mdl in ("gain", "loss", "cnloh"):
            f, flo, fhi = d_to_f(d_hat, mdl), d_to_f(d_lo, mdl), d_to_f(d_hi, mdl)
            pred = f_to_lrr(f, mdl)
            dev = abs(r["lrr_mean"] - pred) / (r["lrr_se"] or 1e9)
            models[mdl] = {"f": f, "f_ci95": [flo, fhi], "pred_lrr": pred, "lrr_dev_se": dev}
        best = min(models, key=lambda x: models[x]["lrr_dev_se"])
        wins = window_profile(A)
        sig_raw = p < alpha_b
        consistent_share = (sum(1 for w in wins if w["p"] < 0.05) / len(wins)) if wins else None
        p_drop = None
        if sig_raw and wins:
            n_drop = max(1, math.ceil(0.2 * len(wins)))
            dropped = sorted(wins, key=lambda w: w["p"])[:n_drop]
            idx = np.where(sarm == A)[0]
            keep_s = np.ones(len(idx), bool)
            for w in dropped:
                keep_s &= ~((sp[idx] >= w["start"]) & (sp[idx] < w["end"]))
            ll2 = arm_likelihood(sk[idx][keep_s], sn[idx][keep_s])
            j2 = int(np.argmax(ll2)); T2 = max(0.0, 2 * (ll2[j2] - ll2[0]))
            p_drop = 0.5 * stats.chi2.sf(T2, 1) if T2 > 0 else 1.0
        spatially_consistent = bool(p_drop is not None and p_drop < alpha_b)
        focal = [w for w in wins if w["p"] < 0.05 / max(1, len(wins)) and w["d"] > 0]
        lrr_gain_f = 2 * (2 ** r["lrr_mean"] - 1)
        lrr_loss_f = 2 * (1 - 2 ** r["lrr_mean"])
        called = sig_raw and spatially_consistent
        r.update(baf_d=d_hat, baf_d_ci95=[d_lo, d_hi], baf_lrt=T, baf_p=p, significant_raw=bool(sig_raw),
                 window_consistency=consistent_share, n_windows=len(wins), p_after_window_deletion=p_drop,
                 significant=bool(called),
                 status=("arm-level event" if called else "focal imbalance (not arm-wide)" if sig_raw else "no evidence"),
                 focal_windows=[f"{A[:-1]}:{w['start']}-{w['end']} d={w['d']:.4f}" for w in focal] if sig_raw and not called else [],
                 best_model=best if called else None, models=models,
                 cell_fraction=models[best]["f"] if called else None,
                 cell_fraction_ci95=models[best]["f_ci95"] if called else None,
                 cell_fraction_upper_bound_if_null={m: models[m]["f_ci95"][1] for m in models} if not called else None,
                 lrr_only_f={"gain": lrr_gain_f, "loss": lrr_loss_f},
                 detection_limit_d=d_min,
                 detection_limit_f={"gain": d_to_f(d_min, "gain"), "loss": d_to_f(d_min, "loss"), "cnloh": d_to_f(d_min, "cnloh")},
                 het_depth_median=float(np.median(n)))
    rows.append(r)

x_nonpar = (chrom == "X") & keep & np.array([not in_par("X", m) for m in mid])
summary = {
    "sex": a.sex, "bin_size": binsize, "bins_total": int(len(B)), "bins_kept": int(keep.sum()),
    "min_mappability": a.min_map, "lrr_sd_autosomal_bins": lrr_sd_auto,
    "het_sites_used": int(len(sk)), "het_depth_median": med_het_dp,
    "reference_bias_mean_alt_fraction": m_bias, "overdispersion_rho": rho,
    "arms_tested": len(tested), "bonferroni_alpha": alpha_b,
    "arms_significant": [r["arm"] for r in rows if r.get("significant")],
    "arms_focal_imbalance": [r["arm"] for r in rows if r.get("status", "").startswith("focal")],
    "het_sites_excluded_segdup": int(n_segdup_excluded[0]),
    "outlier_eps_profiled": EPS, "window_bp": a.window,
    "internal_control_chrX_nonPAR_lrr": float(np.median(lrr[x_nonpar])) if x_nonpar.any() else None,
    "internal_control_expected": -1.0 if a.sex == "M" else 0.0,
    "detection_limit_f_median_across_arms": {
        m: float(np.median([r["detection_limit_f"][m] for r in rows if r.get("tested")])) for m in ("gain", "loss", "cnloh")},
    "caveats": [
        "Blood-derived DNA: the aneuploidy burden in blood may not reflect the tissues driving the phenotype.",
        "In MVA, aneuploidy levels vary between tissues and over time; a single blood sample is one time point in one tissue.",
        f"Sensitivity depends on depth (median het-site depth {med_het_dp:.0f}x) and the number of informative sites per arm; "
        "the per-arm detection limit is reported, and a non-significant arm means below that limit, not absent.",
        "Random mosaic aneuploidy spread across many chromosomes at low level produces little imbalance on any single arm; "
        "per-arm tests are insensitive to it, so the genome-wide LRR/BAF dispersion is reported as well.",
        "Single sample, no panel of normals: GC/mappability correction is internal to this sample.",
    ],
}
# genome-wide dispersion (diffuse random aneuploidy leaves a variance signature, not an arm shift)
summary["genomewide_baf_excess_variance"] = float(np.mean((baf[np.isin(sc, AUTOS)] - m_bias) ** 2)
                                                  - np.mean(m_bias * (1 - m_bias) / sn[np.isin(sc, AUTOS)] * (1 + (sn[np.isin(sc, AUTOS)] - 1) * rho)))

with open(a.prefix + ".arms.json", "w") as fh:
    json.dump(rows, fh, indent=1)
with open(a.prefix + ".summary.json", "w") as fh:
    json.dump(summary, fh, indent=1)
with open(a.prefix + ".arms.tsv", "w") as fh:
    cols = ["arm", "n_bins", "n_het_sites", "lrr_mean", "lrr_se", "baf_d", "baf_p", "window_consistency", "status",
            "best_model", "cell_fraction", "cell_fraction_ci95", "detection_limit_f", "focal_windows"]
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join("" if r.get(c) is None else (json.dumps(r[c]) if isinstance(r.get(c), (list, dict)) else str(r[c])) for c in cols) + "\n")
with gzip.open(a.prefix + ".bins.tsv.gz", "wt") as fh:
    fh.write("chrom\tstart\tend\tgc\tmap_frac\tkept\tdepth\tlrr\n")
    for i in range(len(B)):
        fh.write(f"{chrom[i]}\t{start[i]}\t{end[i]}\t{gc[i]:.4f}\t{mp[i]:.3f}\t{int(keep[i])}\t{dp[i]:.2f}\t{'' if np.isnan(lrr[i]) else f'{lrr[i]:.4f}'}\n")

# ------------------------------------------------------------------ figures
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

offs, x0 = {}, 0
for c in CHROMS:
    offs[c] = x0; x0 += chrom_len.get(c, 0)
rng = np.random.default_rng(1)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 7), sharex=True)
kk = keep & ~np.isnan(lrr)
ax1.scatter([offs[c] + m for c, m in zip(chrom[kk], mid[kk])], lrr[kk], s=1, c="#4a6fa5", alpha=0.4, rasterized=True)
for r in rows:
    if r.get("lrr_mean") is None:
        continue
    c = r["chrom"]; lo, hi = cen.get(c, (0, 0))
    xs = (offs[c], offs[c] + lo) if r["arm"].endswith("p") else (offs[c] + hi, offs[c] + chrom_len[c])
    ax1.plot(xs, [r["lrr_mean"]] * 2, color="#c0392b" if r.get("significant") else "#222", lw=2)
ax1.set_ylim(-1.5, 1.0); ax1.axhline(0, color="#999", lw=0.5); ax1.set_ylabel("LRR (log2 depth ratio)")
sub = rng.choice(len(baf), size=min(60000, len(baf)), replace=False)
ax2.scatter([offs[c] + p for c, p in zip(sc[sub], sp[sub])], baf[sub], s=1, c="#2e8b57", alpha=0.3, rasterized=True)
ax2.set_ylim(0, 1); ax2.axhline(0.5, color="#999", lw=0.5); ax2.set_ylabel("BAF at het SNPs")
for c in CHROMS:
    for ax in (ax1, ax2):
        ax.axvline(offs[c], color="#ddd", lw=0.5)
ax2.set_xticks([offs[c] + chrom_len.get(c, 0) / 2 for c in CHROMS]); ax2.set_xticklabels(CHROMS, fontsize=8)
ax1.set_title("Genome-wide depth (GC/mappability-corrected LRR) and B-allele frequency; red = arm with significant imbalance")
fig.tight_layout(); fig.savefig(a.prefix + ".genome_lrr_baf.png", dpi=150); plt.close(fig)

tr = [r for r in rows if r.get("tested")]
fig, ax = plt.subplots(figsize=(16, 4))
xs = np.arange(len(tr))
lim = [r["detection_limit_f"]["loss"] for r in tr]
ax.bar(xs, lim, color="#e8e8e8", label="detection limit (loss model, 80% power)")
for i, r in enumerate(tr):
    if r.get("significant"):
        lo, hi = r["cell_fraction_ci95"]
        ax.errorbar(i, r["cell_fraction"], yerr=[[r["cell_fraction"] - lo], [hi - r["cell_fraction"]]], fmt="o", color="#c0392b")
ax.set_xticks(xs); ax.set_xticklabels([r["arm"] for r in tr], rotation=90, fontsize=7)
ax.set_ylabel("fraction of cells"); ax.set_title("Mosaic aneuploidy burden per arm: estimated cell fraction (red, 95% CI) vs detection limit (grey)")
ax.legend(loc="upper right"); fig.tight_layout(); fig.savefig(a.prefix + ".arm_burden.png", dpi=150); plt.close(fig)

# per-chromosome depth + BAF panels (one row per chromosome: LRR on top of BAF)
fig, axes = plt.subplots(12, 4, figsize=(18, 30))
for idx, c in enumerate(CHROMS):
    r_, col = divmod(idx, 2)
    axl, axb = axes[r_][col * 2], axes[r_][col * 2 + 1]
    m = kk & (chrom == c)
    axl.scatter(mid[m] / 1e6, lrr[m], s=1, c="#4a6fa5", alpha=0.5, rasterized=True)
    axl.set_ylim(-1.5, 1.0); axl.axhline(0, color="#999", lw=0.5); axl.set_title(f"chr{c} LRR", fontsize=9)
    sm = sc == c
    if sm.any():
        pick = np.where(sm)[0]
        pick = pick if len(pick) <= 8000 else rng.choice(pick, 8000, replace=False)
        axb.scatter(sp[pick] / 1e6, baf[pick], s=1, c="#2e8b57", alpha=0.3, rasterized=True)
    axb.set_ylim(0, 1); axb.axhline(0.5, color="#999", lw=0.5); axb.set_title(f"chr{c} BAF", fontsize=9)
    if c in cen:
        for ax in (axl, axb):
            ax.axvspan(cen[c][0] / 1e6, cen[c][1] / 1e6, color="#eee")
    for ax in (axl, axb):
        ax.tick_params(labelsize=7)
fig.suptitle("Per-chromosome GC/mappability-corrected depth (LRR) and B-allele frequency (x axis: Mb; grey = centromere)")
fig.tight_layout(rect=(0, 0, 1, 0.99)); fig.savefig(a.prefix + ".per_chromosome.png", dpi=110); plt.close(fig)

print(json.dumps({k: v for k, v in summary.items() if k != "caveats"}, indent=1))
