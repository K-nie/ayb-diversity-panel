#!/usr/bin/env python3
"""
Discriminant Analysis of Principal Components (DAPC; Jombart, Devillard &
Balloux 2010, BMC Genetics 11:94).

DAPC complements ADMIXTURE / PCA by finding the linear combinations of
markers that maximise BETWEEN-group separation given the assigned cluster
labels. The two-step structure is:
  1. PCA on the standardised dosage to compress dimensionality and avoid
     small-n LDA failure (n samples << m markers).
  2. Linear Discriminant Analysis on the retained PC scores using cluster
     labels from script 01 (k = 2 k-means partition).
The resulting LD1 axis projects samples onto the maximum-discrimination
direction; each marker has a DAPC loading (the contribution of that marker
to LD1, computed as the dot product of its centred-scaled genotype vector
with the LD1 score divided by sample count).

Outputs (results/29_dapc/)
--------------------------
tables/
    dapc_sample_scores.csv     sample LD1 score + cluster
    dapc_marker_loadings.csv   per-marker contribution to LD1, sorted
figures/
    fig67_dapc_density.png/.pdf       LD1 density by cluster
    fig68_dapc_marker_manhattan.png/.pdf  Manhattan of |loading| on LD1

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from _plotstyle import apply, WONG, CLUSTER_PAL
from _figstyle import adjust_labels, nearest_named_gene, publishable_axes
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
OUT = ROOT / "results" / "29_dapc"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_RETAINED_PCS = 30   # standard "retain enough PCs to capture ~80% variance"


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage + cluster labels
# ---------------------------------------------------------------------------
print("[load] HapMap + cluster labels...")
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
M = dosage.loc[keep_m, keep_s].T
M = M.fillna(M.mean(axis=0)).values
samples = list(dosage.columns[keep_s])
markers = list(dosage.index[keep_m])
print(f"[load] M shape = {M.shape}")

pca_tab = pd.read_csv(PCA_CSV).set_index("sample").reindex(samples)
clu = pca_tab["cluster"].astype(int).values
print(f"[load] clusters: {np.unique(clu, return_counts=True)}")


# ---------------------------------------------------------------------------
# Step 1: PCA on standardised dosage
# ---------------------------------------------------------------------------
print(f"[dapc] PCA on dosage; retaining {N_RETAINED_PCS} PCs...")
p = M.mean(axis=0) / 2.0
p = np.clip(p, 1e-6, 1 - 1e-6)
denom = np.sqrt(2.0 * p * (1.0 - p))
Mz = (M - 2.0 * p[None, :]) / denom[None, :]

pca = PCA(n_components=min(N_RETAINED_PCS, len(samples) - 1, M.shape[1]))
PCs = pca.fit_transform(Mz)
cum_var = np.cumsum(pca.explained_variance_ratio_)
print(f"[dapc] {PCs.shape[1]} PCs retained, cumulative variance = "
      f"{cum_var[-1]:.3f} ({cum_var[len(cum_var)//2]:.3f} at midpoint)")


# ---------------------------------------------------------------------------
# Step 2: LDA on PCs using cluster labels
# ---------------------------------------------------------------------------
lda = LinearDiscriminantAnalysis(n_components=1)
LD1 = lda.fit_transform(PCs, clu).ravel()
# project LDA back to original marker space:
#   loading_marker = Mz_z (centred-scaled) @ scaling vector across PCs
# scaling_marker = pca.components_.T @ lda.scalings_
# Each marker's LD1 contribution = how much its standardised value moves
# samples along LD1.
scal_PC = lda.scalings_[:, 0]                  # (n_PCs,)
loadings_marker = pca.components_.T @ scal_PC  # (n_markers,)

# Save sample-level LD1 scores
sdf = pd.DataFrame({
    "sample": samples,
    "cluster": clu + 1,
    "LD1": LD1,
})
sdf.to_csv(TAB / "dapc_sample_scores.csv", index=False)

# Save marker loadings + AYB anchoring
anc = pd.read_csv(ANCHOR)
mdf = pd.DataFrame({"rs": markers,
                    "loading_LD1": loadings_marker,
                    "abs_loading": np.abs(loadings_marker)})
mdf = mdf.merge(anc[["rs", "chr_cowpea", "pos_cowpea",
                     "chr_ayb", "snp_pos_ayb"]], on="rs", how="left")
mdf = mdf.sort_values("abs_loading", ascending=False)
mdf.to_csv(TAB / "dapc_marker_loadings.csv", index=False)
print(f"[dapc] top-5 LD1 markers (by |loading|):")
print(mdf.head(5)[["rs", "chr_ayb", "snp_pos_ayb",
                   "loading_LD1", "abs_loading"]].to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 67. LD1 density by cluster
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
for k in sorted(set(clu)):
    pts = LD1[clu == k]
    ax.hist(pts, bins=20, alpha=0.7,
            color=CLUSTER_PAL[k], edgecolor="white", linewidth=0.4,
            label=f"Cluster {k + 1} (n = {len(pts)})")
ax.set_xlabel("DAPC LD1 score")
ax.set_ylabel("count of samples")
ax.set_title("DAPC LD1 separating marker-PCA k = 2 clusters")
ax.legend()
ax.grid(True, axis="y", alpha=0.4)
save(fig, "fig67_dapc_density")
print("[fig] fig67_dapc_density")


# ---------------------------------------------------------------------------
# Fig 68. Manhattan of |LD1 loading| (AYB-anchored)
# ---------------------------------------------------------------------------
anc_ok = mdf.dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc_ok["pos"] = anc_ok["snp_pos_ayb"].astype(int)
chrs = sorted(anc_ok["chr_ayb"].unique(),
              key=lambda s: int(re.search(r"\d+", s).group()))
offsets, mids, x_cursor = {}, [], 0
xs = []
anc_ok = anc_ok.sort_values(["chr_ayb", "pos"]).reset_index(drop=True)
for c in chrs:
    sub = anc_ok[anc_ok["chr_ayb"] == c]
    offsets[c] = x_cursor
    xs.extend((sub["pos"].values + x_cursor).tolist())
    mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
    x_cursor += sub["pos"].max() + 5_000_000
anc_ok["x"] = xs
fig, ax = plt.subplots(figsize=(11, 4.0), constrained_layout=True)
palette = [WONG["blue"], "#7c7c7c"]
for i, c in enumerate(chrs):
    sub = anc_ok[anc_ok["chr_ayb"] == c]
    ax.scatter(sub["x"], sub["abs_loading"], s=14,
               color=palette[i % 2], alpha=0.85, edgecolor="none")
# 99th-pct threshold for visual flag
top_thresh = anc_ok["abs_loading"].quantile(0.99)
ax.axhline(top_thresh, color=WONG["vermillion"], ls="--", lw=0.8,
           label=f"99th pct = {top_thresh:.3f}")

# Funannotate nearest-named-gene callouts on the top-15 DAPC LD1 outliers.
gff_path = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
top15 = anc_ok.nlargest(15, "abs_loading")
callout_texts = []
seen: set[str] = set()
for _, r in top15.iterrows():
    gene = nearest_named_gene(r["chr_ayb"], r["pos"], gff_path,
                                window_bp=200_000)
    if not gene or gene in seen:
        continue
    seen.add(gene)
    callout_texts.append(ax.text(r["x"], r["abs_loading"], gene,
                                  fontsize=7.5, fontstyle="italic",
                                  color=WONG["vermillion"], zorder=10))
adjust_labels(callout_texts, ax=ax,
               expand=(1.2, 1.5), force_text=(0.8, 1.2))

ax.set_xticks(mids); ax.set_xticklabels(chrs, fontsize=8)
ax.set_ylabel("|DAPC LD1 marker loading|")
ax.set_title(f"DAPC marker contributions to LD1 "
             f"({len(anc_ok)} AYB-anchored SNPs) -- "
             f"top-15 markers labelled with nearest named gene")
ax.legend(framealpha=0.9, edgecolor="none")
publishable_axes(ax, grid="y")
save(fig, "fig68_dapc_marker_manhattan")
print("[fig] fig68_dapc_marker_manhattan")

print(f"\nOutputs in: {OUT}")
