#!/usr/bin/env python3
"""
Pairwise genetic distance and UPGMA / NJ dendrograms for the 105-line AYB panel.

Distance: Rogers' modified genetic distance for biallelic SNP dosages
    D_R(i, j) = sqrt( (1 / (8 L)) * sum_l (d_il - d_jl)^2 )
where d is 0/1/2 dosage and L is the number of QC-filtered markers.

Tree: average-linkage UPGMA from the upper-triangular Rogers' distance vector.
Leaves are coloured by the k=2 PCA cluster from script 01. Output includes a
sample-by-sample distance heatmap with rows / columns reordered by the same
linkage so the two genetic subgroups are visually contiguous.

Outputs (results/05_distance_dendrogram/)
-----------------------------------------
tables/
    rogers_distance_matrix.csv     105 x 105 pairwise Rogers' distance
    sample_dendro_order.csv         leaf order (rebreed-ready)
figures/
    fig15_dendrogram_upgma.png/.pdf
    fig16_distance_heatmap.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform

from _plotstyle import apply, WONG, CLUSTER_PAL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "05_distance_dendrogram"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Rebuild filtered dosage (same filter as scripts 01-03)
# ---------------------------------------------------------------------------
print("[load] reading HapMap...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META]
calls = hm[sample_cols].astype(str)


def dosage_row(row, alleles):
    ref, alt = alleles.split("/")
    out = np.full(len(row), np.nan)
    arr = row.values
    out[arr == ref + ref] = 0.0
    out[(arr == ref + alt) | (arr == alt + ref)] = 1.0
    out[arr == alt + alt] = 2.0
    return out


dosage = np.vstack([dosage_row(calls.iloc[i], hm["alleles"].iloc[i])
                    for i in range(len(hm))])
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)

call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]

keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values

M = dosage.loc[keep_m, keep_s].T  # samples x markers
M = M.fillna(M.mean(axis=0))
samples = M.index.tolist()
L = M.shape[1]
print(f"[load] M shape = {M.shape}  (samples x markers)")


# ---------------------------------------------------------------------------
# Rogers' modified distance
# ---------------------------------------------------------------------------
print("[dist] computing Rogers' modified distance...")
A = M.values
sq = (A ** 2).sum(axis=1)
D2 = sq[:, None] + sq[None, :] - 2.0 * (A @ A.T)
D2 = np.maximum(D2, 0.0)
np.fill_diagonal(D2, 0.0)
D = np.sqrt(D2 / (8.0 * L))
np.fill_diagonal(D, 0.0)

D_df = pd.DataFrame(D, index=samples, columns=samples)
D_df.to_csv(TAB / "rogers_distance_matrix.csv")
print(f"[dist] median = {D[np.triu_indices(len(samples), 1)].mean():.3f}, "
      f"max = {D.max():.3f}")


# ---------------------------------------------------------------------------
# UPGMA linkage
# ---------------------------------------------------------------------------
print("[tree] UPGMA linkage...")
condensed = squareform(D, checks=False)
Z = linkage(condensed, method="average")

# Group labels from PCA k=2 clustering
pca_tab = pd.read_csv(PCA_TAB).set_index("sample").loc[samples]
groups = pca_tab["cluster"].astype(int).values
g_color = [CLUSTER_PAL[g] for g in groups]


# ---------------------------------------------------------------------------
# Fig 15. UPGMA dendrogram coloured by PCA cluster
# ---------------------------------------------------------------------------
def color_func(_):
    return "#555555"


fig, ax = plt.subplots(figsize=(13, 5.5), constrained_layout=True)
dd = dendrogram(
    Z, labels=samples, ax=ax, leaf_font_size=6,
    leaf_rotation=90, color_threshold=0, above_threshold_color="#555555",
    link_color_func=color_func,
)
# colour leaf-tick labels by group
xt = ax.get_xticklabels()
sample_to_group = dict(zip(samples, g_color))
for tick in xt:
    tick.set_color(sample_to_group.get(tick.get_text(), "black"))
ax.set_ylabel("Rogers' modified distance")
ax.set_title(f"UPGMA dendrogram — 105 Nigerian AYB lines, {L:,} QC SNPs\n"
             f"leaf colour = PCA k-means cluster")
# legend
import matplotlib.patches as mpatches
handles = [mpatches.Patch(color=CLUSTER_PAL[k],
                          label=f"Cluster {k + 1}  (n = {(groups == k).sum()})")
           for k in sorted(set(groups))]
ax.legend(handles=handles, loc="upper right", frameon=False)
save(fig, "fig15_dendrogram_upgma")
print("[fig] fig15_dendrogram_upgma")


# Save dendrogram leaf order so other scripts can reorder samples consistently
leaf_order = pd.DataFrame({"leaf_index": np.arange(len(dd["ivl"])),
                           "sample": dd["ivl"]})
leaf_order["cluster"] = leaf_order["sample"].map(
    dict(zip(samples, groups)))
leaf_order.to_csv(TAB / "sample_dendro_order.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 15b. Circular (radial) UPGMA dendrogram
# Walks the same scipy linkage Z but renders in polar coordinates: each leaf
# is placed at theta = 2 * pi * i / N, and parent nodes connect by an arc at
# their distance r. Matches the standard radial-tree convention used in
# population-genetics papers (e.g. fineSTRUCTURE, GerpFold).
# ---------------------------------------------------------------------------
print("[fig] circular UPGMA dendrogram...")
from scipy.cluster.hierarchy import dendrogram as _dendro

# Get the icoord / dcoord drawing arrays from a non-plotted dendrogram so we
# can re-render them in polar coordinates. scipy returns coordinates on a
# "logical" x-axis 5, 15, 25 ... (one per leaf). We rescale to [0, 2*pi).
dd_polar = _dendro(Z, labels=samples, no_plot=True,
                   color_threshold=0)
leaf_order_polar = dd_polar["ivl"]
N = len(leaf_order_polar)

# convert each link's icoords/dcoords to (theta, r); each icoord vector is
# length-4 [x1, x1, x2, x2] horizontals + verticals.
def x_to_theta(x):
    # scipy puts leaves at 5, 15, 25... so leaf index = (x - 5) / 10
    return 2.0 * np.pi * ((x - 5.0) / 10.0) / N


fig = plt.figure(figsize=(8.5, 8.5), constrained_layout=False)
ax_r = fig.add_subplot(111, projection="polar")
for icoord, dcoord in zip(dd_polar["icoord"], dd_polar["dcoord"]):
    # icoord is [x1, x1, x2, x2]; convert each to theta
    thetas = [x_to_theta(x) for x in icoord]
    rs = list(dcoord)
    # render the U: vertical from (theta1, r=0_link) to apex, arc at apex,
    # vertical down from apex to (theta2, r=0_link).
    # scipy gives points: (x1, y1)=(left_low) (x1, y2)=(left_high)
    # (x2, y2)=(right_high) (x2, y1)=(right_low)
    # Map: left vertical = constant theta=theta1, r varies from rs[0] to rs[1].
    ax_r.plot([thetas[0], thetas[0]], [rs[0], rs[1]],
              color="#555555", linewidth=0.7)
    ax_r.plot([thetas[3], thetas[3]], [rs[3], rs[2]],
              color="#555555", linewidth=0.7)
    # arc at the top of the U: draw a smooth arc between theta1 and theta2
    # at r = rs[1] = rs[2]
    n_arc = 20
    if thetas[1] != thetas[2]:
        arc_thetas = np.linspace(thetas[1], thetas[2], n_arc)
        ax_r.plot(arc_thetas, np.full(n_arc, rs[1]),
                  color="#555555", linewidth=0.7)

# leaf positions + cluster-coloured tip points
for i, lab in enumerate(leaf_order_polar):
    theta_i = 2.0 * np.pi * i / N
    g = groups[samples.index(lab)]
    ax_r.scatter(theta_i, 0, s=22, color=CLUSTER_PAL[g],
                 edgecolor="white", linewidth=0.4, zorder=5)
    # tip labels — small, rotated to read radially
    rot_deg = np.degrees(theta_i)
    if 90 < rot_deg < 270:
        rot_deg -= 180
        ha = "right"
    else:
        ha = "left"
    ax_r.text(theta_i, max(D.max() * 1.04, 0.05), lab,
              rotation=rot_deg, rotation_mode="anchor",
              fontsize=5, ha=ha, va="center",
              color=CLUSTER_PAL[g])

ax_r.set_theta_zero_location("N")
ax_r.set_theta_direction(-1)
ax_r.set_rlabel_position(0)
ax_r.set_rlim(0, D.max() * 1.12)
ax_r.set_yticklabels([])
ax_r.set_xticklabels([])
ax_r.grid(True, alpha=0.25)
ax_r.set_title(f"Circular UPGMA dendrogram — {N} AYB lines, {L:,} QC SNPs",
               fontsize=11, pad=18)
# legend
handles = [mpatches.Patch(color=CLUSTER_PAL[k],
                          label=f"Cluster {k + 1}  (n = {(groups == k).sum()})")
           for k in sorted(set(groups))]
ax_r.legend(handles=handles, loc="lower right",
            bbox_to_anchor=(1.05, -0.05), frameon=False, fontsize=9)
fig.savefig(FIG / "fig15b_dendrogram_circular.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "fig15b_dendrogram_circular.pdf", bbox_inches="tight")
plt.close(fig)
print("[fig] fig15b_dendrogram_circular")


# ---------------------------------------------------------------------------
# Fig 16. Distance heatmap reordered by linkage, group bar on the side
# ---------------------------------------------------------------------------
print("[fig] distance heatmap...")
order_idx = [samples.index(s) for s in dd["ivl"]]
D_ord = D[np.ix_(order_idx, order_idx)]
g_ord = groups[order_idx]

fig = plt.figure(figsize=(8.5, 7.5), constrained_layout=False)
gs = fig.add_gridspec(2, 2, width_ratios=[0.03, 1], height_ratios=[0.03, 1],
                      wspace=0.02, hspace=0.02)
ax_top = fig.add_subplot(gs[0, 1])
ax_left = fig.add_subplot(gs[1, 0])
ax_main = fig.add_subplot(gs[1, 1])

# heatmap
im = ax_main.imshow(D_ord, cmap="magma", aspect="auto",
                    vmin=0, vmax=np.percentile(D_ord, 99))
ax_main.set_xticks([]); ax_main.set_yticks([])
ax_main.set_xlabel("Lines (UPGMA leaf order)")
ax_main.set_ylabel("Lines (UPGMA leaf order)")

# group color bars
for ax, axis in [(ax_top, 0), (ax_left, 1)]:
    cb = np.array([[CLUSTER_PAL[g] for g in g_ord]])
    if axis == 0:
        ax.imshow([list(range(len(g_ord)))], aspect="auto",
                  cmap=plt.matplotlib.colors.ListedColormap(
                      [CLUSTER_PAL[g] for g in g_ord]))
    else:
        ax.imshow(np.arange(len(g_ord))[:, None], aspect="auto",
                  cmap=plt.matplotlib.colors.ListedColormap(
                      [CLUSTER_PAL[g] for g in g_ord]))
    ax.set_xticks([]); ax.set_yticks([])

# colorbar
cb_ax = fig.add_axes([0.92, 0.13, 0.018, 0.55])
fig.colorbar(im, cax=cb_ax, label="Rogers' distance")

fig.suptitle("Pairwise genetic distance — reordered by UPGMA linkage",
             y=0.94, fontsize=12)
fig.savefig(FIG / "fig16_distance_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "fig16_distance_heatmap.pdf", bbox_inches="tight")
plt.close(fig)
print("[fig] fig16_distance_heatmap")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print(f"\nOutputs in: {OUT}")
print(f"  rogers_distance_matrix.csv  ({len(samples)} x {len(samples)})")
print(f"  sample_dendro_order.csv     ({len(samples)} rows)")
print(f"  fig15_dendrogram_upgma.png  + .pdf")
print(f"  fig16_distance_heatmap.png  + .pdf")
