#!/usr/bin/env python3
"""97 - Is PC1 driven by heterozygosity itself? (Reviewer Item 5, part B)

The panel splits on PC1 into 11 near-homozygous minor-cluster lines and 84 more
heterozygous main-cluster lines. A reviewer asked us to rule out the trivial
explanation that PC1 is a heterozygosity axis rather than a structure axis, by
re-deriving the split from data in which heterozygosity cannot contribute:
  (i)  PCA on haploidised calls (every heterozygous site resolved to its
       major-allele homozygote -> a strictly homozygous 0/2 matrix), and
  (ii) PCoA on an identity-by-state (IBS) distance matrix.
If the 11/84 split and PC1 ordering survive both, PC1 indexes structure, not het.

Inputs
------
- AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv
  (same file, same QC cascade as script 01: sample call-rate >= 0.90,
   marker call-rate >= 0.90 & MAF >= 0.05 -> 1,625 markers x 95 samples)

Outputs (results/97_pc1_heterozygosity_robustness/)
--------------------------------------------------
tables/pc1_het_summary.csv         per-sample PC1 (baseline/haploid/IBS), obs het, cluster
tables/concordance_summary.csv     PC1-het correlation + cross-method ARI + split recovery
figures/fig_pc1_heterozygosity.{png,pdf}   4-panel supplementary figure

Author: Benjamin Narh-Madey
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT = ROOT / "results" / "97_pc1_heterozygosity_robustness"
FIG = OUT / "figures"; TAB = OUT / "tables"
for d in (FIG, TAB): d.mkdir(parents=True, exist_ok=True)

SAMP_CR = MARK_CR = 0.90; MIN_MAF = 0.05
RNG = 0

# ---- load + dosage (identical logic to script 01) ----
hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#","alleles","chrom","pos","strand","assembly#","center",
        "protLSID","assayLSID","panelLSID","QCcode"]
scols = [c for c in hm.columns if c not in META]
calls = hm[scols].astype(str)

def dosage_row(row, alleles):
    ref, alt = alleles.split("/")
    out = np.full(len(row), np.nan); arr = row.values
    out[arr == ref+ref] = 0.0
    out[(arr == ref+alt) | (arr == alt+ref)] = 1.0
    out[arr == alt+alt] = 2.0
    return out

dos = np.vstack([dosage_row(calls.iloc[i], hm["alleles"].iloc[i]) for i in range(len(hm))])
dos = pd.DataFrame(dos, index=hm["rs#"].values, columns=scols)

# ---- QC ----
ns = dos.shape[1]
called_m = dos.notna().sum(axis=1); cr_m = called_m / ns
mean_dose = dos.mean(axis=1, skipna=True); maf = np.minimum(mean_dose/2, 1-mean_dose/2)
cr_s = dos.notna().sum(axis=0) / dos.shape[0]
het_s = (dos == 1).sum(axis=0) / dos.notna().sum(axis=0)  # per-sample observed het (all markers)

keep_s = cr_s[cr_s >= SAMP_CR].index.tolist()
keep_m = maf[(cr_m >= MARK_CR) & (maf >= MIN_MAF)].index.tolist()
X = dos.loc[keep_m, keep_s].T.copy()      # samples x markers
print(f"[qc] samples {len(keep_s)}/{ns}  markers {len(keep_m)}/{dos.shape[0]}")

# per-sample observed het on the post-QC marker set (the axis the reviewer worries about)
obs_het = (X == 1).sum(axis=1) / X.notna().sum(axis=1)

# ---- baseline: mean-impute, standardise, PCA (reproduces manuscript Fig 2A) ----
Xi = X.fillna(X.mean(axis=0))
Xs = StandardScaler().fit_transform(Xi.values)
pca0 = PCA(n_components=5, random_state=RNG).fit(Xs)
pc0 = pca0.transform(Xs)
km0 = KMeans(n_clusters=2, n_init=25, random_state=RNG).fit(pc0[:, :2])
lab0 = km0.labels_
# name the smaller cluster "minor"
minor0 = 1 if (lab0 == 1).sum() < (lab0 == 0).sum() else 0
is_minor = (lab0 == minor0)
print(f"[baseline] cluster sizes: minor={is_minor.sum()} main={(~is_minor).sum()}  "
      f"PC1 var={pca0.explained_variance_ratio_[0]*100:.1f}%")

# ---- (i) haploidised calls: het -> major-allele homozygote ----
alt_freq = Xi.mean(axis=0) / 2.0            # per marker, on imputed matrix
alt_major = (alt_freq > 0.5).values         # True -> het maps to 2, else to 0
H = X.copy().values
het_mask = (H == 1)
# assign het to major-allele homozygote
col_major_val = np.where(alt_major, 2.0, 0.0)
H = np.where(het_mask, np.broadcast_to(col_major_val, H.shape), H)
Hdf = pd.DataFrame(H, index=X.index, columns=X.columns)
Hi = Hdf.fillna(Hdf.mean(axis=0))
resid_het = float((Hi == 1).values.mean())  # should be ~0 (imputation can create fractional)
Hs = StandardScaler().fit_transform(Hi.values)
# drop zero-variance columns produced by haploidisation
nzv = Hs.std(axis=0) > 0
Hs = Hs[:, nzv]
pcaH = PCA(n_components=5, random_state=RNG).fit(Hs)
pcH = pcaH.transform(Hs)
kmH = KMeans(n_clusters=2, n_init=25, random_state=RNG).fit(pcH[:, :2])
ari_H = adjusted_rand_score(lab0, kmH.labels_)
print(f"[haploid] markers used={nzv.sum()}  PC1 var={pcaH.explained_variance_ratio_[0]*100:.1f}%  "
      f"ARI vs baseline={ari_H:.3f}")

# ---- (ii) IBS distance -> PCoA (classical MDS) ----
D = Xi.values
n = D.shape[0]
ibs_dist = np.zeros((n, n))
for i in range(n):
    diff = np.abs(D[i] - D).mean(axis=1) / 2.0   # mean allele-count difference / 2
    ibs_dist[i] = diff
# classical MDS on squared distances
D2 = ibs_dist ** 2
J = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * J @ D2 @ J
w, V = np.linalg.eigh(B)
order = np.argsort(w)[::-1]
w = w[order]; V = V[:, order]
pos = w > 1e-9
coords = V[:, pos] * np.sqrt(w[pos])
pcoa1_var = w[pos][0] / w[pos].sum() * 100
kmI = KMeans(n_clusters=2, n_init=25, random_state=RNG).fit(coords[:, :2])
ari_I = adjusted_rand_score(lab0, kmI.labels_)
print(f"[ibs] PCoA1 var={pcoa1_var:.1f}%  ARI vs baseline={ari_I:.3f}")

# ---- orient axes so minor cluster is on the same side for readability ----
def orient(x):
    return x if np.median(x[is_minor]) >= np.median(x[~is_minor]) else -x
pc0_1 = orient(pc0[:, 0]); pcH_1 = orient(pcH[:, 0]); pcI_1 = orient(coords[:, 0])

# ---- correlations of PC1 with heterozygosity ----
r_pear, p_pear = stats.pearsonr(pc0_1, obs_het.values)
r_spear, p_spear = stats.spearmanr(pc0_1, obs_het.values)
# variance in PC1 explained by het (R^2) and cross-method PC1 concordance
r2_het = r_pear ** 2
r_h0, _ = stats.pearsonr(pc0_1, pcH_1)
r_i0, _ = stats.pearsonr(pc0_1, pcI_1)
# het difference between clusters
het_minor = obs_het[is_minor].mean(); het_main = obs_het[~is_minor].mean()

# ---- tables ----
summ = pd.DataFrame({
    "sample": X.index,
    "obs_het_postQC": obs_het.values,
    "cluster": np.where(is_minor, "minor", "main"),
    "PC1_baseline": pc0_1, "PC1_haploid": pcH_1, "PCoA1_ibs": pcI_1,
})
summ.to_csv(TAB / "pc1_het_summary.csv", index=False)

conc = pd.DataFrame([
    {"metric": "PC1_baseline_vs_obs_het_pearson_r", "value": r_pear, "p": p_pear},
    {"metric": "PC1_baseline_vs_obs_het_spearman_rho", "value": r_spear, "p": p_spear},
    {"metric": "PC1_variance_explained_by_het_R2", "value": r2_het, "p": np.nan},
    {"metric": "haploid_PCA_ARI_vs_baseline_kmeans", "value": ari_H, "p": np.nan},
    {"metric": "IBS_PCoA_ARI_vs_baseline_kmeans", "value": ari_I, "p": np.nan},
    {"metric": "PC1_baseline_vs_haploid_pearson_r", "value": r_h0, "p": np.nan},
    {"metric": "PC1_baseline_vs_IBS_pearson_r", "value": r_i0, "p": np.nan},
    {"metric": "residual_het_fraction_after_haploidisation", "value": resid_het, "p": np.nan},
    {"metric": "mean_obs_het_minor", "value": het_minor, "p": np.nan},
    {"metric": "mean_obs_het_main", "value": het_main, "p": np.nan},
    {"metric": "n_markers_baseline", "value": len(keep_m), "p": np.nan},
    {"metric": "n_markers_haploid_nonzerovar", "value": int(nzv.sum()), "p": np.nan},
])
conc.to_csv(TAB / "concordance_summary.csv", index=False)
print(conc.to_string(index=False))

# ---- figure (Nature style: sans-serif, lowercase bold panel letters, colourblind-safe) ----
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})
CB_MINOR = "#D55E00"   # vermillion
CB_MAIN = "#0072B2"    # blue
cols = np.where(is_minor, CB_MINOR, CB_MAIN)

fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 150/25.4))
(axA, axB), (axC, axD) = axes

# (a) baseline PCA coloured by observed het
sc = axA.scatter(pc0_1, orient(pc0[:, 1]) if False else pc0[:, 1], c=obs_het.values,
                 cmap="viridis", s=26, edgecolor="k", linewidth=0.3)
axA.set_xlabel(f"PC1 ({pca0.explained_variance_ratio_[0]*100:.1f}%)")
axA.set_ylabel(f"PC2 ({pca0.explained_variance_ratio_[1]*100:.1f}%)")
cb = fig.colorbar(sc, ax=axA, fraction=0.046, pad=0.04); cb.set_label("observed heterozygosity", fontsize=7)

# (b) PC1 vs observed het
axB.scatter(obs_het.values, pc0_1, c=cols, s=26, edgecolor="k", linewidth=0.3)
axB.set_xlabel("per-sample observed heterozygosity")
axB.set_ylabel("PC1 (baseline)")
axB.text(0.04, 0.94, f"Pearson r = {r_pear:.2f}\n$R^2$ = {r2_het:.2f}",
         transform=axB.transAxes, va="top", ha="left", fontsize=7.5)

# (c) haploidised-call PCA
axC.scatter(pcH_1, pcH[:, 1], c=cols, s=26, edgecolor="k", linewidth=0.3)
axC.set_xlabel(f"PC1 haploidised ({pcaH.explained_variance_ratio_[0]*100:.1f}%)")
axC.set_ylabel(f"PC2 haploidised ({pcaH.explained_variance_ratio_[1]*100:.1f}%)")
axC.text(0.04, 0.94, f"ARI vs baseline = {ari_H:.2f}", transform=axC.transAxes,
         va="top", ha="left", fontsize=7.5)

# (d) IBS PCoA
axD.scatter(pcI_1, coords[:, 1], c=cols, s=26, edgecolor="k", linewidth=0.3)
axD.set_xlabel(f"PCoA1 IBS ({pcoa1_var:.1f}%)")
axD.set_ylabel(f"PCoA2 IBS ({w[pos][1]/w[pos].sum()*100:.1f}%)")
axD.text(0.04, 0.94, f"ARI vs baseline = {ari_I:.2f}", transform=axD.transAxes,
         va="top", ha="left", fontsize=7.5)

# legend for cluster colours
from matplotlib.lines import Line2D
handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=CB_MINOR, markeredgecolor='k',
                  markersize=6, label='minor cluster (n=11)'),
           Line2D([0],[0], marker='o', color='w', markerfacecolor=CB_MAIN, markeredgecolor='k',
                  markersize=6, label='main cluster (n=84)')]
axB.legend(handles=handles, fontsize=6.5, loc="lower right", frameon=True)

for ax, L in zip([axA, axB, axC, axD], "abcd"):
    ax.text(-0.12, 1.06, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

fig.tight_layout()
fig.savefig(FIG / "fig_pc1_heterozygosity.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "fig_pc1_heterozygosity.pdf", bbox_inches="tight")
print("[fig] wrote fig_pc1_heterozygosity.{png,pdf}")
