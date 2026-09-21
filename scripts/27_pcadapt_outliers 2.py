#!/usr/bin/env python3
"""
PCAdapt-style F_ST outlier scan for candidate-selection loci.

Reference: Luu, Bazin & Blum (2017) Molecular Ecology Resources 17:67-77;
also implemented in pcadapt R package. Method:

  1. Run PCA on the standardised dosage matrix (n samples x m markers).
  2. For each marker j, regress its standardised genotype on the K leading
     PCs (here K = 2, matching the marker-PCA cluster solution).
  3. The vector of regression slopes z_j across PCs forms a per-marker
     "loading vector".
  4. Compute the Mahalanobis distance of each z_j from the multivariate
     mean using the inverse covariance of the loadings.
  5. Mahalanobis^2 follows chi^2_K under the null of neutral drift.
  6. Convert to p-values + BH-FDR; flag selection candidates.

This is the principled within-panel selection-genomics analogue of BayeScan
or OutFLANK that's been picked up in plant-breeding Q1 papers (e.g. Pavan
et al. 2020 common bean, Cantelmo et al. 2019 rice).

Outputs (results/27_pcadapt_outliers/)
--------------------------------------
tables/
    pcadapt_per_locus.csv       per-marker Mahalanobis stat, p, q
    pcadapt_outlier_top50.csv   top-50 by q value (AYB-anchored only)
figures/
    fig64_pcadapt_manhattan.png/.pdf
    fig65_pcadapt_qq.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

from _plotstyle import apply, WONG
from _figstyle import adjust_labels, nearest_named_gene, publishable_axes
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
OUT = ROOT / "results" / "27_pcadapt_outliers"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_PCs = 2


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage
# ---------------------------------------------------------------------------
print("[load] HapMap...")
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
M = dosage.loc[keep_m, keep_s].T.values
M = pd.DataFrame(M).fillna(pd.DataFrame(M).mean()).values  # samples x markers
samples = dosage.columns[keep_s].tolist()
markers = dosage.index[keep_m].tolist()
print(f"[load] M shape = {M.shape}  (samples x markers)")


# ---------------------------------------------------------------------------
# PCAdapt
# ---------------------------------------------------------------------------
print(f"[pcadapt] running PCA + Mahalanobis with K = {N_PCs} PCs...")
# standardise across samples (per-marker centering, scale by sqrt(2 p (1-p)))
p = M.mean(axis=0) / 2.0
p = np.clip(p, 1e-6, 1 - 1e-6)
denom = np.sqrt(2.0 * p * (1.0 - p))
Mz = (M - 2.0 * p[None, :]) / denom[None, :]

pca = PCA(n_components=N_PCs)
scores = pca.fit_transform(Mz)
# z_j = regression of column j of Mz on the K scores -> shape (n_markers, K)
# == Mz.T @ scores / scores' scores (one regression per marker; since scores
# are orthonormal up to scaling we just project Mz onto each PC).
scores_norm = scores / scores.std(axis=0, ddof=1)
z_load = (Mz.T @ scores_norm) / (len(samples) - 1)
print(f"[pcadapt] loadings shape: {z_load.shape}")
print(f"[pcadapt] variance explained by {N_PCs} PCs: "
      f"{pca.explained_variance_ratio_}")

# Mahalanobis distance per marker
cov_z = np.cov(z_load.T)
inv_cov = np.linalg.inv(cov_z)
mu = z_load.mean(axis=0)
delta = z_load - mu
mahal = np.einsum("ij,jk,ik->i", delta, inv_cov, delta)

# chi^2_K p-values
pvals = 1.0 - stats.chi2.cdf(mahal, df=N_PCs)
q_bh = stats.false_discovery_control(pvals)
print(f"[pcadapt] BH q < 0.05 outliers: {(q_bh < 0.05).sum()}")
print(f"[pcadapt] BH q < 0.10 outliers: {(q_bh < 0.10).sum()}")

# Join to AYB anchoring + cowpea fallback
anc = pd.read_csv(ANCHOR)
res = pd.DataFrame({
    "rs": markers,
    "z_PC1": z_load[:, 0],
    "z_PC2": z_load[:, 1],
    "mahal": mahal,
    "p": pvals,
    "q_bh": q_bh,
}).merge(anc[["rs", "chr_cowpea", "pos_cowpea", "chr_ayb", "snp_pos_ayb"]],
         on="rs", how="left")
res.to_csv(TAB / "pcadapt_per_locus.csv", index=False)

top = res.sort_values("p").head(50)
top.to_csv(TAB / "pcadapt_outlier_top50.csv", index=False)
print("[pcadapt] top-10 by p:")
print(top.head(10)[["rs", "chr_ayb", "snp_pos_ayb",
                     "z_PC1", "z_PC2", "mahal", "p", "q_bh"]].to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 64. PCAdapt Manhattan (AYB-anchored)
# ---------------------------------------------------------------------------
anc_ok = res.dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc_ok["pos"] = anc_ok["snp_pos_ayb"].astype(int)
chrs = sorted(anc_ok["chr_ayb"].unique(),
              key=lambda s: int(re.search(r"\d+", s).group()))
fig, ax = plt.subplots(figsize=(11, 4.0), constrained_layout=True)
offsets, mids, x_cursor = {}, [], 0
xs = []
for c in chrs:
    sub = anc_ok[anc_ok["chr_ayb"] == c]
    offsets[c] = x_cursor
    xs.extend((sub["pos"].values + x_cursor).tolist())
    mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
    x_cursor += sub["pos"].max() + 5_000_000
anc_ok["x"] = xs
palette = [WONG["blue"], "#7c7c7c"]
for i, c in enumerate(chrs):
    sub = anc_ok[anc_ok["chr_ayb"] == c]
    ax.scatter(sub["x"], -np.log10(sub["p"]), s=14,
               color=palette[i % 2], alpha=0.85, edgecolor="none")
# highlight FDR < 0.10
sig = anc_ok[anc_ok["q_bh"] < 0.10]
ax.scatter(sig["x"], -np.log10(sig["p"]), s=42,
           facecolor="none", edgecolor=WONG["vermillion"], linewidth=1.2,
           label=f"BH q < 0.10 (n = {len(sig)})")
alpha_bonf = 0.05 / len(anc_ok)
ax.axhline(-np.log10(alpha_bonf), ls="--", color="black", linewidth=0.7,
           label=f"Bonferroni $\\alpha$ = {alpha_bonf:.2e}")

# Funannotate gene callouts on the top-15 PCAdapt outliers. Use BH-q
# when there is any FDR-significant hit; otherwise fall back to the
# smallest raw p-values so the panel still communicates which markers the
# scan flagged as most differentiated.
gff_path = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
if len(sig) >= 5:
    sig_top = sig.nsmallest(15, "q_bh")
else:
    sig_top = anc_ok.nsmallest(15, "p")
callout_texts = []
seen: set[str] = set()
for _, r in sig_top.iterrows():
    gene = nearest_named_gene(r["chr_ayb"], r["pos"], gff_path,
                                window_bp=200_000)
    if not gene or gene in seen:
        continue
    seen.add(gene)
    callout_texts.append(ax.text(r["x"], -np.log10(r["p"]), gene,
                                  fontsize=7.5, fontstyle="italic",
                                  color=WONG["vermillion"], zorder=10))
adjust_labels(callout_texts, ax=ax,
               expand=(1.2, 1.5), force_text=(0.8, 1.2))

ax.set_xticks(mids); ax.set_xticklabels(chrs, fontsize=8)
ax.set_ylabel(r"$-\log_{10}\,p$" "\n(PCAdapt Mahalanobis$^2$)")
ax.set_title(f"PCAdapt-style F$_{{ST}}$ outlier scan, K = {N_PCs} PCs -- "
              f"top-15 outliers labelled with nearest named gene")
ax.legend(loc="upper right", fontsize=8, framealpha=0.9, edgecolor="none")
publishable_axes(ax, grid="y")
save(fig, "fig64_pcadapt_manhattan")
print("[fig] fig64_pcadapt_manhattan")


# ---------------------------------------------------------------------------
# Fig 65. PCAdapt QQ
# ---------------------------------------------------------------------------
p_sorted = np.sort(pvals)
n = len(p_sorted)
expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
observed = -np.log10(p_sorted)
fig, ax = plt.subplots(figsize=(5.5, 5), constrained_layout=True)
ax.scatter(expected, observed, s=8, color=WONG["blue"], edgecolor="none")
lim = max(expected.max(), observed.max()) + 0.3
ax.plot([0, lim], [0, lim], "k--", lw=1)
# Beta CI band
k = np.arange(1, n + 1)
lo = -np.log10(stats.beta.ppf(0.975, k, n - k + 1))
hi = -np.log10(stats.beta.ppf(0.025, k, n - k + 1))
ax.fill_between(expected, lo, hi, color="grey", alpha=0.18, linewidth=0)
ax.set_xlabel(r"Expected $-\log_{10}\,p$")
ax.set_ylabel(r"Observed $-\log_{10}\,p$")
lambda_gc = float(np.median(mahal) / stats.chi2.ppf(0.5, df=N_PCs))
ax.set_title(f"PCAdapt QQ — n = {n} markers; "
             f"$\\lambda$$_{{GC}}$ = {lambda_gc:.2f}", fontsize=11)
ax.grid(True, alpha=0.4)
save(fig, "fig65_pcadapt_qq")
print("[fig] fig65_pcadapt_qq")

print(f"\nOutputs in: {OUT}")
