#!/usr/bin/env python3
"""
Breeding-decision outputs from the AYB diversity panel.

1. Trait-extreme accession shortlists
   * Top-10 and bottom-10 lines per trait (mean basis).
   * Bivariate pairs of interest for crossing:
       low Tannin  x high Antioxidant
       low Phenol  x high Flavonoid
       high Antioxidant x high Flavonoid

2. Core-collection construction
   Greedy MaxMin selection on Rogers' distance — at each step, add the line
   that maximises the minimum distance to the already-selected core. Seeded
   with the most distant pair. Target size 20 lines (~ 19 % of the panel),
   which is the conventional 10-20 % core size in plant-genetic-resources
   work (Frankel & Brown 1984; van Hintum 2000).

Outputs (results/06_breeding/)
------------------------------
tables/
    trait_extremes_<TRAIT>.csv          top 10 and bottom 10 per trait
    bivariate_extreme_pairs.csv         joint trait extremes
    core_collection_top20.csv           ordered list of selected lines
    core_collection_diversity_audit.csv pairwise-distance summary core vs full
figures/
    fig17_trait_extreme_lines.png/.pdf  4-panel ranked bar chart
    fig18_core_collection_pca.png/.pdf  PCA with core lines highlighted

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG, CLUSTER_PAL, TRAIT_PAL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PHENO_XLSX = ROOT / "Wet_Chemistry_Data.xlsx"
DIST_CSV = ROOT / "results" / "05_distance_dendrogram" / "tables" / "rogers_distance_matrix.csv"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "06_breeding"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

CORE_SIZE = 20
TRAITS = ["Tannin", "Phenol", "Flavonoid", "Antioxidant"]


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
pheno = (pd.read_excel(PHENO_XLSX, sheet_name="Means")
         .rename(columns={"Genotypes": "sample", "Flavinoid": "Flavonoid"}))
D = pd.read_csv(DIST_CSV, index_col=0)
samples = D.index.tolist()
pca = pd.read_csv(PCA_TAB).set_index("sample").loc[samples]


# ---------------------------------------------------------------------------
# 1. Trait extremes
# ---------------------------------------------------------------------------
print("[extremes] per-trait top/bottom 10...")
for t in TRAITS:
    df = pheno[["sample", t]].dropna().copy()
    df = df.sort_values(t)
    bottom = df.head(10).copy()
    top = df.tail(10).copy()
    bottom["rank"] = "bottom_10"
    top["rank"] = "top_10"
    out = pd.concat([bottom, top.iloc[::-1]], ignore_index=True)
    out["rank_within_panel_pct"] = out[t].rank(method="min", pct=True).round(3)
    out.to_csv(TAB / f"trait_extremes_{t}.csv", index=False)
    print(f"  {t}: bottom range {bottom[t].min():.1f}-{bottom[t].max():.1f}; "
          f"top range {top[t].min():.1f}-{top[t].max():.1f}")


# Bivariate pairs (low / high quantiles)
def quant_set(df, col, low_q, high_q):
    lo = df[df[col] <= df[col].quantile(low_q)]["sample"].tolist()
    hi = df[df[col] >= df[col].quantile(high_q)]["sample"].tolist()
    return lo, hi


print("[extremes] bivariate pairs for crossing...")
piv = pheno.dropna(subset=TRAITS)
lo_tan, hi_tan = quant_set(piv, "Tannin", 0.25, 0.75)
lo_phe, hi_phe = quant_set(piv, "Phenol", 0.25, 0.75)
hi_fla = set(quant_set(piv, "Flavonoid", 0.25, 0.75)[1])
hi_ant = set(quant_set(piv, "Antioxidant", 0.25, 0.75)[1])

biv_rows = []
for s in piv["sample"]:
    val = piv[piv["sample"] == s].iloc[0]
    if (s in lo_tan) and (s in hi_ant):
        biv_rows.append({"sample": s, "category": "low Tannin x high Antioxidant",
                         "Tannin": val["Tannin"], "Antioxidant": val["Antioxidant"]})
    if (s in lo_phe) and (s in hi_fla):
        biv_rows.append({"sample": s, "category": "low Phenol x high Flavonoid",
                         "Phenol": val["Phenol"], "Flavonoid": val["Flavonoid"]})
    if (s in hi_ant) and (s in hi_fla):
        biv_rows.append({"sample": s, "category": "high Antioxidant x high Flavonoid",
                         "Antioxidant": val["Antioxidant"], "Flavonoid": val["Flavonoid"]})
biv = pd.DataFrame(biv_rows)
biv.to_csv(TAB / "bivariate_extreme_pairs.csv", index=False)
print(f"  {len(biv)} bivariate-extreme records written.")


# ---------------------------------------------------------------------------
# 2. Core collection — greedy MaxMin
# ---------------------------------------------------------------------------
print(f"[core] greedy MaxMin selection of {CORE_SIZE} lines...")
N = len(samples)
Dm = D.values.copy()

# initial pair = most-distant pair
iu, ju = np.unravel_index(np.argmax(Dm + np.tril(np.full_like(Dm, -1.0))), Dm.shape)
selected = [samples[iu], samples[ju]]
remaining = [s for s in samples if s not in selected]

while len(selected) < CORE_SIZE:
    sel_idx = [samples.index(s) for s in selected]
    best_s, best_min = None, -np.inf
    for s in remaining:
        i = samples.index(s)
        d_to_core = Dm[i, sel_idx].min()
        if d_to_core > best_min:
            best_min = d_to_core
            best_s = s
    selected.append(best_s)
    remaining.remove(best_s)

core = pd.DataFrame({"order_added": range(1, len(selected) + 1),
                     "sample": selected})
core = core.merge(pca[["PC1", "PC2", "cluster"]], left_on="sample", right_index=True, how="left")
core = core.merge(pheno, on="sample", how="left")
core.to_csv(TAB / "core_collection_top20.csv", index=False)
print(f"[core] core sample IDs (order of addition):")
for r in core.itertuples():
    print(f"  {r.order_added:>2d}  {r.sample}  cluster {r.cluster}")

# Diversity audit
all_pairs = Dm[np.triu_indices(N, 1)]
core_idx = [samples.index(s) for s in selected]
core_pairs = Dm[np.ix_(core_idx, core_idx)][np.triu_indices(len(core_idx), 1)]
audit = pd.DataFrame({
    "metric": ["mean_pairwise_distance",
               "median_pairwise_distance",
               "min_pairwise_distance",
               "max_pairwise_distance",
               "n_pairs"],
    "full_panel": [all_pairs.mean(), np.median(all_pairs),
                   all_pairs.min(), all_pairs.max(), len(all_pairs)],
    "core_20":    [core_pairs.mean(), np.median(core_pairs),
                   core_pairs.min(), core_pairs.max(), len(core_pairs)],
})
audit.to_csv(TAB / "core_collection_diversity_audit.csv", index=False)
print("\n[core] diversity audit:")
print(audit.to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 17. Trait-extreme bar chart
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.5), constrained_layout=True)
for ax, t in zip(axes.ravel(), TRAITS):
    df = pheno[["sample", t]].dropna().sort_values(t)
    bottom = df.head(10)
    top = df.tail(10).iloc[::-1]
    show = pd.concat([bottom, top])
    colors = ([WONG["skyblue"]] * 10) + ([TRAIT_PAL[t]] * 10)
    y = np.arange(len(show))
    ax.barh(y, show[t], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(show["sample"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(t)
    ax.set_title(f"{t} — bottom-10 (blue) and top-10 (colour)")
    ax.grid(True, axis="x")
fig.suptitle("Trait-extreme accessions per seed-biochemistry trait",
             y=1.02, fontsize=12)
save(fig, "fig17_trait_extreme_lines")
print("[fig] fig17_trait_extreme_lines")


# ---------------------------------------------------------------------------
# Fig 18. PCA with core collection highlighted
# ---------------------------------------------------------------------------
pca_full = pd.read_csv(PCA_TAB)
core_set = set(selected)
fig, ax = plt.subplots(figsize=(6.4, 5.6), constrained_layout=True)
for k in sorted(pca_full["cluster"].unique()):
    m = (pca_full["cluster"] == k) & (~pca_full["sample"].isin(core_set))
    ax.scatter(pca_full.loc[m, "PC1"], pca_full.loc[m, "PC2"],
               s=22, color=CLUSTER_PAL[k], edgecolor="white",
               linewidth=0.4, alpha=0.55,
               label=f"Cluster {k + 1} ({m.sum()})")
core_mask = pca_full["sample"].isin(core_set)
ax.scatter(pca_full.loc[core_mask, "PC1"], pca_full.loc[core_mask, "PC2"],
           s=85, facecolor="none", edgecolor="black", linewidth=1.5,
           marker="o", label=f"Core ({CORE_SIZE})", zorder=5)
# label core sample IDs
for _, r in pca_full[core_mask].iterrows():
    ax.annotate(r["sample"], (r["PC1"], r["PC2"]),
                xytext=(4, 4), textcoords="offset points",
                fontsize=7, color="black")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title(f"Core collection — {CORE_SIZE}-line MaxMin set on Rogers' distance")
ax.legend(loc="upper right")
ax.grid(True)
save(fig, "fig18_core_collection_pca")
print("[fig] fig18_core_collection_pca")


# ---------------------------------------------------------------------------
print(f"\nOutputs in: {OUT}")
