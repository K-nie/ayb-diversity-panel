#!/usr/bin/env python3
"""
Weir & Cockerham (1984) per-locus and global F_ST between PCA clusters.

The k-means partition from script 01 splits the panel into Cluster 1 (n = 94,
PC1 low) and Cluster 2 (n = 11, PC1 high). This script:

1. Computes per-locus W&C 1984 estimators of a (between-pop), b (between-ind
   within-pop) and c (within-ind) variance components.
2. Global F_ST = sum_l(a_l) / sum_l(a_l + b_l + c_l) — the ratio-of-sums form
   that W&C recommend (not the mean of per-locus ratios).
3. Bootstrap 95 % CI on the global F_ST by resampling loci with replacement.
4. Per-marker F_ST Manhattan with the empirical 99th percentile flagged as the
   "outlier" tail (puts a defensible label on differentiation outliers without
   over-claiming significance at n = 105).

Note: the cluster assignment was learned from PCA on the same SNPs, so the
F_ST values here are descriptive, not inferential. This caveat is reported in
the summary text.

Outputs (results/10_fst/)
-------------------------
tables/
    fst_per_locus.csv                rs, chrom, pos, a, b, c, fst
    fst_global_summary.csv           global F_ST + bootstrap CI
figures/
    fig25_fst_distribution.png/.pdf  per-locus F_ST histogram
    fig26_fst_manhattan.png/.pdf     per-locus F_ST Manhattan (anchored SNPs)

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "10_fst"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_BOOT = 1000
RNG = np.random.default_rng(20250516)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage (same cascade as scripts 01-03)
# ---------------------------------------------------------------------------
hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META]
calls = hm[sample_cols].astype(str)

ref_alt = hm["alleles"].str.split("/", expand=True)
ref_alt.columns = ["ref", "alt"]
dosage = np.full((len(hm), len(sample_cols)), np.nan)
for i in range(len(hm)):
    ref, alt = ref_alt.iloc[i]
    arr = calls.iloc[i].values
    dosage[i, arr == ref + ref] = 0.0
    dosage[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
    dosage[i, arr == alt + alt] = 2.0

dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values

dos_f = dosage.loc[keep_m, keep_s]
samples = dos_f.columns.tolist()
markers = dos_f.index.tolist()
chrom_raw = hm.loc[keep_m, "chrom"].values
pos_raw = pd.to_numeric(hm.loc[keep_m, "pos"], errors="coerce").values
M = dos_f.values   # markers x samples, may have NaN
print(f"[load] {len(markers)} markers x {len(samples)} samples")


# Cluster labels
pca = pd.read_csv(PCA_TAB).set_index("sample").loc[samples]
grp = pca["cluster"].astype(int).values
pop_labels = sorted(set(grp))
print(f"[grp] populations = {pop_labels}, sizes = "
      f"{[int((grp == k).sum()) for k in pop_labels]}")


# ---------------------------------------------------------------------------
# Weir & Cockerham 1984 per-locus, vectorised over loci
# ---------------------------------------------------------------------------
# Reference: equations 5-7 of Weir & Cockerham (1984), Evolution 38:1358-1370.
# For diploid biallelic data, with r populations of sample size n_i:
#   p_i  = allele frequency in pop i
#   h_i  = observed heterozygote frequency in pop i
#   n_bar = mean sample size
#   p_bar = sum(n_i p_i) / sum(n_i)  (weighted mean freq)
#   S^2 = (1 / ((r-1) n_bar)) * sum_i n_i (p_i - p_bar)^2
#   h_bar = sum(n_i h_i) / sum(n_i)
#   n_c = ((sum n_i) - sum(n_i^2)/sum(n_i)) / (r - 1)
#   a = (n_bar / n_c) * (S^2 - (1/(n_bar-1)) * (p_bar(1-p_bar) - ((r-1)/r) S^2 - h_bar/4))
#   b = (n_bar/(n_bar-1)) * (p_bar(1-p_bar) - ((r-1)/r) S^2 - ((2*n_bar - 1)/(4*n_bar)) h_bar)
#   c = h_bar / 2
def wc_fst_per_locus(M_, grp_):
    L, N = M_.shape
    pops = sorted(set(grp_))
    r = len(pops)
    a = np.zeros(L); b = np.zeros(L); c = np.zeros(L)
    # accumulate per-pop n, p, h
    n_i = np.zeros((r, L))     # non-missing per pop per locus
    p_i = np.zeros((r, L))     # alt allele freq
    h_i = np.zeros((r, L))     # heterozygote freq
    for ki, k in enumerate(pops):
        cols = np.where(grp_ == k)[0]
        sub = M_[:, cols]
        nm = ~np.isnan(sub)
        n_i[ki] = nm.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            p_i[ki] = np.nansum(sub, axis=1) / (2.0 * n_i[ki])
            h_i[ki] = (sub == 1).sum(axis=1) / n_i[ki]
    n_tot = n_i.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        n_bar = n_tot / r
        p_bar = (n_i * p_i).sum(axis=0) / n_tot
        S2 = (1.0 / ((r - 1) * n_bar)) * (n_i * (p_i - p_bar) ** 2).sum(axis=0)
        h_bar = (n_i * h_i).sum(axis=0) / n_tot
        n_c = (n_tot - (n_i ** 2).sum(axis=0) / n_tot) / (r - 1)
        common = p_bar * (1.0 - p_bar) - ((r - 1) / r) * S2 - h_bar / 4.0
        a = (n_bar / n_c) * (S2 - (1.0 / (n_bar - 1)) * common)
        b = (n_bar / (n_bar - 1)) * (
            p_bar * (1.0 - p_bar) - ((r - 1) / r) * S2
            - ((2 * n_bar - 1) / (4 * n_bar)) * h_bar
        )
        c = h_bar / 2.0
    return a, b, c, p_bar, n_tot


print("[wc] computing per-locus W&C variance components...")
a, b, c, p_bar, n_tot = wc_fst_per_locus(M, grp)
denom = a + b + c
with np.errstate(invalid="ignore", divide="ignore"):
    fst_per = np.where(denom > 0, a / denom, np.nan)
valid = np.isfinite(fst_per)
print(f"[wc] {valid.sum():,} of {len(markers):,} loci with finite F_ST")

global_fst = np.nansum(a[valid]) / np.nansum(denom[valid])
print(f"[wc] global F_ST (ratio of sums) = {global_fst:.4f}")

# Bootstrap 95 % CI by resampling loci
print(f"[wc] bootstrapping global F_ST ({N_BOOT} reps)...")
idx = np.where(valid)[0]
boot = np.empty(N_BOOT)
for i in range(N_BOOT):
    s = RNG.choice(idx, size=len(idx), replace=True)
    boot[i] = a[s].sum() / denom[s].sum()
lo, hi = np.percentile(boot, [2.5, 97.5])
print(f"[wc] bootstrap 95 % CI = [{lo:.4f}, {hi:.4f}]")


# Save tables
per_locus = pd.DataFrame({
    "rs": markers,
    "chrom": chrom_raw,
    "pos": pos_raw,
    "n_tot": n_tot,
    "p_bar": p_bar,
    "a": a, "b": b, "c": c,
    "fst": fst_per,
})
per_locus.to_csv(TAB / "fst_per_locus.csv", index=False)

global_tab = pd.DataFrame([{
    "estimator": "Weir & Cockerham 1984 (ratio of sums)",
    "global_fst": global_fst,
    "boot_ci95_low": lo,
    "boot_ci95_high": hi,
    "n_loci_used": int(valid.sum()),
    "n_pops": len(pop_labels),
    "n_total": int(len(samples)),
    "n_per_pop": ";".join(str(int((grp == k).sum())) for k in pop_labels),
    "caveat": ("Cluster labels were learned by k-means on PCA of the same SNPs; "
               "F_ST is therefore descriptive of the partition, not an "
               "independent test of differentiation."),
}])
global_tab.to_csv(TAB / "fst_global_summary.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 25. F_ST distribution
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
ax.hist(fst_per[valid], bins=50, color=WONG["blue"],
        edgecolor="white", linewidth=0.4)
ax.axvline(global_fst, color=WONG["vermillion"], linewidth=1.5,
           label=f"global F$_{{ST}}$ = {global_fst:.3f}\n"
                 f"95 % CI [{lo:.3f}, {hi:.3f}]")
q99 = np.nanpercentile(fst_per, 99)
ax.axvline(q99, color="black", linestyle="--", linewidth=1.0,
           label=f"99th pct = {q99:.3f}")
ax.set_xlabel("per-locus F$_{ST}$ (W&C 1984)")
ax.set_ylabel("count")
ax.set_title(f"F$_{{ST}}$ distribution — Cluster 1 vs Cluster 2 "
             f"(n = {(grp == 0).sum()} vs {(grp == 1).sum()})")
ax.legend()
ax.grid(True, axis="y")
save(fig, "fig25_fst_distribution")
print("[fig] fig25_fst_distribution")


# ---------------------------------------------------------------------------
# Fig 26. Per-locus F_ST Manhattan (anchored SNPs only)
# ---------------------------------------------------------------------------
per_locus["chr"] = per_locus["chrom"].astype(str).str.extract(
    r"(Vu\d+)", expand=False)
anc = per_locus.dropna(subset=["chr", "pos", "fst"]).copy()
anc["pos"] = anc["pos"].astype(int)
anc = anc.sort_values(["chr", "pos"]).reset_index(drop=True)
chrs = sorted(anc["chr"].unique(),
              key=lambda s: int(re.search(r"\d+", s).group()))

offsets, mids, xs = {}, [], []
x_cursor = 0
for c in chrs:
    sub = anc[anc["chr"] == c]
    offsets[c] = x_cursor
    xs.extend((sub["pos"].values + x_cursor).tolist())
    mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
    x_cursor += sub["pos"].max() + 5_000_000
anc["x"] = xs

fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
palette = [WONG["blue"], "#7c7c7c"]
for i, c in enumerate(chrs):
    sub = anc[anc["chr"] == c]
    ax.scatter(sub["x"], sub["fst"], s=14, color=palette[i % 2],
               alpha=0.85, edgecolor="none")
ax.axhline(global_fst, color=WONG["vermillion"], linewidth=1.0,
           label=f"global F$_{{ST}}$ = {global_fst:.3f}")
ax.axhline(q99, color="black", linestyle="--", linewidth=0.8,
           label=f"99th pct = {q99:.3f}")
ax.set_xticks(mids)
ax.set_xticklabels(chrs, rotation=0, fontsize=8)
ax.set_ylabel(r"F$_{ST}$")
ax.set_title(f"Per-locus F$_{{ST}}$ across the cowpea-anchored "
             f"{len(anc):,} SNPs")
ax.set_ylim(bottom=min(0, np.nanmin(anc["fst"]) - 0.02))
ax.grid(True, axis="y")
ax.legend(loc="upper right")
save(fig, "fig26_fst_manhattan")
print("[fig] fig26_fst_manhattan")


print(f"\nOutputs in: {OUT}")
print(f"  global F_ST = {global_fst:.4f}  [{lo:.4f}, {hi:.4f}]")
