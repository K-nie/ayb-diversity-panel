#!/usr/bin/env python3
"""
Genomic relationship matrix (G) and breeding-cross shortlist for the AYB panel.

1. Build G via VanRaden (2008) method 1:
       Z = M - 2P    (P_jk = allele freq of marker j)
       G = Z Z^T / (2 sum_j p_j (1 - p_j))
2. Reorder by hierarchical clustering of (1 - G); plot a heatmap with PCA-
   cluster colour-bars on the side.
3. From the off-diagonal of G, rank cross pairs in two directions:
     a. Least-related (lowest G_ij)   -> maximal Mendelian sampling
        variance, useful for outcrossing to expand allelic diversity in F1.
     b. Most-related (highest G_ij)   -> useful for line-maintenance /
        seed-rejuvenation crosses where genetic identity is the goal.
4. Cross-join each top pair to the full 13-trait phenotype frame so the
   breeder reads the genetic and phenotypic complement together. For each
   pair, every trait gets midparent value and |A - B| absolute gap.

Outputs (results/13_grm_crosspairs/)
------------------------------------
tables/
    grm_vanraden.csv                 105 x 105 G
    cross_pairs_least_related.csv    top 50 least-related pairs + 13 traits
    cross_pairs_most_related.csv     top 50 most-related pairs + 13 traits
    cross_pair_diversity_audit.csv   G summary (mean / sd / range)
figures/
    fig31_grm_heatmap.png/.pdf
    fig32_grm_offdiag_distribution.png/.pdf
    fig33_crosspair_trait_scatter.png/.pdf   gap vs G for the top traits

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from _plotstyle import apply, WONG, CLUSTER_PAL
from _pheno import load_phenotypes, TRAITS_ALL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "13_grm_crosspairs"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
TOP_N_PAIRS = 50
# Traits to highlight in fig33 (the four with strongest cross-pair signal —
# any pair with a wide gap on these is breeder-interesting).
HIGHLIGHT = ["Crude_Protein", "Total_Oxalate", "Antioxidant", "Mass_of_Seeds"]


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage -> G (VanRaden method 1)
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
M = dosage.loc[keep_m, keep_s].T
M = M.fillna(M.mean(axis=0)).values
samples = dosage.columns[keep_s].tolist()
print(f"[load] M shape = ({len(samples)}, {M.shape[1]})")

p = M.mean(axis=0) / 2.0
P = 2.0 * p
Z = M - P[None, :]
denom = (2.0 * (p * (1.0 - p))).sum()
G = Z @ Z.T / denom
G_df = pd.DataFrame(G, index=samples, columns=samples)
G_df.to_csv(TAB / "grm_vanraden.csv")
print(f"[grm] G diag mean = {np.diag(G).mean():.3f}, "
      f"off-diag mean = {G[np.triu_indices(len(samples), 1)].mean():.3f}")


# ---------------------------------------------------------------------------
# Hierarchical reorder
# ---------------------------------------------------------------------------
D_grm = 1.0 - G
np.fill_diagonal(D_grm, 0.0)
D_sym = (D_grm + D_grm.T) / 2.0
D_sym = np.maximum(D_sym, 0.0)
condensed = squareform(D_sym, checks=False)
Z_link = linkage(condensed, method="average")
order = leaves_list(Z_link)
G_ord = G[np.ix_(order, order)]

pca = pd.read_csv(PCA_TAB).set_index("sample").loc[samples]
clu = pca["cluster"].astype(int).values
clu_ord = clu[order]


# ---------------------------------------------------------------------------
# Fig 31. GRM heatmap with PCA-cluster colour bars
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(8.5, 7.4), constrained_layout=False)
gs = fig.add_gridspec(2, 2, width_ratios=[0.03, 1], height_ratios=[0.03, 1],
                      wspace=0.02, hspace=0.02)
ax_top = fig.add_subplot(gs[0, 1])
ax_left = fig.add_subplot(gs[1, 0])
ax_main = fig.add_subplot(gs[1, 1])
upper = G_ord[np.triu_indices(len(samples), 1)]
vmin, vmax = np.percentile(upper, [2, 98])
im = ax_main.imshow(G_ord, cmap="RdBu_r", aspect="auto", vmin=vmin, vmax=vmax)
ax_main.set_xticks([]); ax_main.set_yticks([])
ax_main.set_xlabel("samples (UPGMA leaf order)")
ax_main.set_ylabel("samples (UPGMA leaf order)")
clu_colors = [CLUSTER_PAL[c] for c in clu_ord]
ax_top.imshow([list(range(len(clu_ord)))], aspect="auto",
              cmap=plt.matplotlib.colors.ListedColormap(clu_colors))
ax_left.imshow(np.arange(len(clu_ord))[:, None], aspect="auto",
               cmap=plt.matplotlib.colors.ListedColormap(clu_colors))
for ax in (ax_top, ax_left):
    ax.set_xticks([]); ax.set_yticks([])
cb_ax = fig.add_axes([0.93, 0.13, 0.018, 0.55])
fig.colorbar(im, cax=cb_ax, label=r"G$_{ij}$ (VanRaden method 1)")
fig.suptitle("Genomic relationship matrix — reordered by UPGMA on (1 - G)",
             y=0.94, fontsize=12)
fig.savefig(FIG / "fig31_grm_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "fig31_grm_heatmap.pdf", bbox_inches="tight")
plt.close(fig)
print("[fig] fig31_grm_heatmap")


# ---------------------------------------------------------------------------
# Fig 32. Off-diagonal G distribution
# ---------------------------------------------------------------------------
ti, tj = np.triu_indices(len(samples), 1)
g_off = G[ti, tj]
same = clu[ti] == clu[tj]
fig, ax = plt.subplots(figsize=(7.5, 4.4), constrained_layout=True)
bins = np.linspace(g_off.min(), g_off.max(), 50)
ax.hist(g_off[same], bins=bins, alpha=0.7, color=WONG["blue"],
        edgecolor="white", linewidth=0.4, label="within-cluster pair")
ax.hist(g_off[~same], bins=bins, alpha=0.7, color=WONG["vermillion"],
        edgecolor="white", linewidth=0.4, label="between-cluster pair")
ax.axvline(g_off.mean(), color="black", linewidth=1.0,
           label=f"overall mean = {g_off.mean():.3f}")
ax.set_xlabel(r"G$_{ij}$ (off-diagonal entries)")
ax.set_ylabel("count of pairs")
ax.set_title("Distribution of genomic relationships across the panel")
ax.legend(loc="upper right")
ax.grid(True, axis="y")
save(fig, "fig32_grm_offdiag_distribution")
print("[fig] fig32_grm_offdiag_distribution")


# ---------------------------------------------------------------------------
# Cross-pair tables — full 13-trait context per pair
# ---------------------------------------------------------------------------
pheno = load_phenotypes()


def build_pair_table(idx_pairs):
    rows = []
    for i, j in idx_pairs:
        a, b = samples[i], samples[j]
        row = {
            "parent_A": a, "parent_B": b,
            "G_ij": float(G[i, j]),
            "cluster_A": int(clu[i]) + 1,
            "cluster_B": int(clu[j]) + 1,
        }
        for t in TRAITS_ALL:
            va = pheno.at[a, t] if a in pheno.index else np.nan
            vb = pheno.at[b, t] if b in pheno.index else np.nan
            row[f"{t}_A"] = va
            row[f"{t}_B"] = vb
            if pd.notna(va) and pd.notna(vb):
                row[f"{t}_midparent"] = (va + vb) / 2.0
                row[f"{t}_absgap"] = abs(va - vb)
            else:
                row[f"{t}_midparent"] = np.nan
                row[f"{t}_absgap"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# Least-related (lowest G_ij)
order_low = np.argsort(g_off)[:TOP_N_PAIRS]
least = [(ti[k], tj[k]) for k in order_low]
least_df = build_pair_table(least)
least_df.to_csv(TAB / "cross_pairs_least_related.csv", index=False)
print(f"[pair] least-related top {TOP_N_PAIRS}: "
      f"G range [{least_df['G_ij'].min():.3f}, {least_df['G_ij'].max():.3f}]")

# Most related (highest off-diag G_ij — usually within-cluster duplicates or
# near-sibs; useful for line maintenance, not breeding gain)
order_high = np.argsort(-g_off)[:TOP_N_PAIRS]
most = [(ti[k], tj[k]) for k in order_high]
most_df = build_pair_table(most)
most_df.to_csv(TAB / "cross_pairs_most_related.csv", index=False)
print(f"[pair] most-related top {TOP_N_PAIRS}: "
      f"G range [{most_df['G_ij'].min():.3f}, {most_df['G_ij'].max():.3f}]")

# Diversity audit
audit = pd.DataFrame([{
    "metric": "all_pairs",
    "n": len(g_off),
    "mean_G": g_off.mean(),
    "sd_G": g_off.std(ddof=1),
    "min_G": g_off.min(),
    "max_G": g_off.max(),
    "p05_G": np.percentile(g_off, 5),
    "p95_G": np.percentile(g_off, 95),
}, {
    "metric": "within_cluster",
    "n": int(same.sum()),
    "mean_G": g_off[same].mean(),
    "sd_G": g_off[same].std(ddof=1),
    "min_G": g_off[same].min(),
    "max_G": g_off[same].max(),
    "p05_G": np.percentile(g_off[same], 5),
    "p95_G": np.percentile(g_off[same], 95),
}, {
    "metric": "between_cluster",
    "n": int((~same).sum()),
    "mean_G": g_off[~same].mean(),
    "sd_G": g_off[~same].std(ddof=1),
    "min_G": g_off[~same].min(),
    "max_G": g_off[~same].max(),
    "p05_G": np.percentile(g_off[~same], 5),
    "p95_G": np.percentile(g_off[~same], 95),
}])
audit.to_csv(TAB / "cross_pair_diversity_audit.csv", index=False)
print("\n[audit]")
print(audit.to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 33. Trait gap vs G for the top 50 least-related pairs, for 4 traits
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.5), constrained_layout=True)
PALETTE = {
    "Crude_Protein":   "#7570b3",
    "Total_Oxalate":   "#1b9e77",
    "Antioxidant":     "#CC79A7",
    "Mass_of_Seeds":   "#6a3d9a",
}
for ax, t in zip(axes.ravel(), HIGHLIGHT):
    col = f"{t}_absgap"
    sub = least_df.dropna(subset=[col])
    ax.scatter(sub["G_ij"], sub[col], s=35,
               color=PALETTE.get(t, WONG["blue"]), edgecolor="black",
               linewidth=0.4, alpha=0.85)
    ax.set_xlabel(r"G$_{ij}$ (least-related top 50)")
    ax.set_ylabel(f"|{t}$_A$ - {t}$_B$|")
    ax.set_title(f"{t}: phenotype gap vs genetic distance\n"
                 f"({len(sub)} of {TOP_N_PAIRS} pairs with both phenotypes)",
                 fontsize=9)
    ax.grid(True)
fig.suptitle(f"Top {TOP_N_PAIRS} least-related cross pairs — "
             f"breeder-relevant trait gaps", y=1.02, fontsize=12)
save(fig, "fig33_crosspair_trait_scatter")
print("[fig] fig33_crosspair_trait_scatter")

print(f"\nOutputs in: {OUT}")
