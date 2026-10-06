#!/usr/bin/env python3
"""Regenerate the Ss05 outlier-attribution strip on the DE-DUPLICATED
unique-SNP definition (Reviewer R1-M7 / R2-C).

Background
----------
The submitted Figure 7C plotted the PCAdapt and DAPC top-50 Ss05 outliers as
separate method rows (script 55, plot_ss05_attribution_strip). A SNP flagged
by both scans therefore appeared twice, and the "Cluster-2-driven" legend
counted method x SNP rows (6 PCAdapt + 7 DAPC = 13) rather than unique loci.
The manuscript now reports the de-duplicated figure: of the 16 unique Ss05
outlier SNPs, 7 (44 %) are Cluster-2-driven. This script re-draws the strip on
the unique-SNP set so the figure matches the corrected text and the permutation
null (script 94).

Attribution (unchanged thresholds, from script 55):
    Cluster-2-driven  if MAF_C2 >= 0.40 AND MAF_C1 <= 0.10
    Cluster-1-driven  if MAF_C1 >= 0.40 AND MAF_C2 <= 0.10
    panel-wide        otherwise
MAF_C1 = main cluster (n = 84, "Cluster 1"); MAF_C2 = minor cluster
(n = 11, "Cluster 2"). This matches the R2-C standardised numbering.

Input
-----
results/55_private_alleles_per_cluster/tables/Ss05_outlier_attribution.csv
    22 method x SNP rows, 16 unique SNPs.

Outputs
-------
results/95_ss05_attribution_unique_figure/figures/fig_Ss05_attribution_strip.png/.pdf
Also overwrites the manuscript copy share/figures/main/A1_NEW_3a_fig_Ss05_attribution_strip.*

Author: Benjamin Narh-Madey
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
ATTR_CSV = ROOT / "results" / "55_private_alleles_per_cluster" / "tables" / "Ss05_outlier_attribution.csv"
OUT = ROOT / "results" / "95_ss05_attribution_unique_figure"
FIG = OUT / "figures"
TAB = OUT / "tables"
SHARE = ROOT / "manuscript" / "paper_A1" / "share" / "figures" / "main"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

ATTR_C_HIGH = 0.40
ATTR_C_LOW = 0.10


def main() -> None:
    attr = pd.read_csv(ATTR_CSV)
    # De-duplicate to unique SNPs. A SNP's attribution is identical across the
    # methods that flagged it (it is a function of MAF_C1/MAF_C2 only), so the
    # first row per rs is representative.
    uniq = attr.drop_duplicates("rs").reset_index(drop=True)
    n_total = len(uniq)
    counts = uniq["attribution"].value_counts().to_dict()
    n_c2 = counts.get("Cluster-2-driven", 0)
    n_c1 = counts.get("Cluster-1-driven", 0)
    n_pw = counts.get("panel-wide", 0)
    print(f"[attr] unique SNPs={n_total}: Cluster-2-driven={n_c2}, "
          f"Cluster-1-driven={n_c1}, panel-wide={n_pw}")
    uniq.to_csv(TAB / "Ss05_outlier_attribution_unique.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 6.2))

    ax.add_patch(plt.Rectangle((0, ATTR_C_HIGH), ATTR_C_LOW, 1 - ATTR_C_HIGH,
                               facecolor=WONG["vermillion"], alpha=0.10, zorder=0))
    ax.text(ATTR_C_LOW / 2, (ATTR_C_HIGH + 1.0) / 2 - 0.05,
            "Cluster-2-driven\nregion", ha="center", va="center",
            fontsize=9, color=WONG["vermillion"], alpha=0.75)
    ax.add_patch(plt.Rectangle((ATTR_C_HIGH, 0), 1 - ATTR_C_HIGH, ATTR_C_LOW,
                               facecolor=WONG["blue"], alpha=0.10, zorder=0))
    ax.text((ATTR_C_HIGH + 1.0) / 2 - 0.05, ATTR_C_LOW / 2,
            "Cluster-1-driven\nregion", ha="center", va="center",
            fontsize=9, color=WONG["blue"], alpha=0.75)

    colour_map = {"Cluster-2-driven": WONG["vermillion"],
                  "Cluster-1-driven": WONG["blue"],
                  "panel-wide": WONG["yellow"]}
    for tag, colour in colour_map.items():
        sub = uniq.loc[uniq["attribution"] == tag]
        if len(sub) == 0:
            continue
        ax.scatter(sub["MAF_C1"], sub["MAF_C2"], color=colour, marker="o",
                   s=90, edgecolor="black", linewidth=0.6,
                   label=f"{tag} (n = {len(sub)})", alpha=0.9, zorder=3)

    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=0.7,
            label="equal between clusters", zorder=1)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Freq of panel-wide minor allele in Cluster 1 (main, n = 84)")
    ax.set_ylabel("Freq of panel-wide minor allele in Cluster 2 (minor, n = 11)")
    ax.set_title(f"Ss05 outlier attribution: {n_c2} of {n_total} unique outliers "
                 f"Cluster-2-driven")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    for stem in (FIG / "fig_Ss05_attribution_strip",
                 SHARE / "A1_NEW_3a_fig_Ss05_attribution_strip"):
        fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[done] wrote {FIG/'fig_Ss05_attribution_strip.png'} and share copy")


if __name__ == "__main__":
    main()
