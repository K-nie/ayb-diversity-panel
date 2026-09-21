#!/usr/bin/env python3
"""
Cross-panel diversity placement: this panel against published IITA AYB
DArTseq panels (Shitta 2022, n = 169; Olomitutu 2022, n = 195).

This is the LIGHTWEIGHT version of the cross-panel comparison originally
planned in analysis_plans_A1.md. The original plan was to re-anchor the
Shitta and Olomitutu marker tables to the S. stenocarpa reference and
recompute every diversity statistic under matched QC. That plan blocked
on data access: neither paper publishes the raw marker table in its
supplementary (Shitta's data is "available from the corresponding author
on reasonable request"; Olomitutu's supplementary only carries the 12
seed-size-associated SNPs).

The pragmatic alternative is to pull the published headline diversity
statistics from each paper and place our panel inside the published AYB
DArTseq diversity range using the numbers we can directly cite. The
comparison is weaker than a matched-pipeline re-analysis -- each panel's
QC thresholds differ -- but is honest, immediately tractable, and
delivers the "reviewer-disarming context" the manuscript needs without
waiting on the author-email round trip.

Published numbers and the section / table they come from are encoded in
PUB_STATS below; the original PMC URLs are documented in the README.

Outputs (results/53_cross_panel_diversity/)
-------------------------------------------
tables/
    per_panel_diversity_stats.csv   wide table: stat x (this, Shitta, Olomitutu)
    cross_panel_qc_comparison.csv    QC threshold and n_marker side-by-side
figures/
    fig_cross_panel_forest.png/.pdf  forest of published values + our panel

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
MARKER_QC_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "marker_qc.csv"
OUT = ROOT / "results" / "53_cross_panel_diversity"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# QC thresholds applied here so the panel-side numbers exactly match the
# manuscript's headline stats (call rate >= 0.90, MAF >= 0.05).
CALL_RATE_MIN = 0.90
MAF_MIN = 0.05


# ---------------------------------------------------------------------------
# Published numbers from the two comparator panels.
# Numbers are quoted verbatim from the cited PMC version of each paper.
# Provenance:
#   Shitta 2022   -> Sci Rep 12:4437; PMC8924269
#                    https://pmc.ncbi.nlm.nih.gov/articles/PMC8924269/
#   Olomitutu 2022 -> Genes 13:2350; PMC9777823
#                     https://pmc.ncbi.nlm.nih.gov/articles/PMC9777823/
# "n.r." = not reported in the paper.
# ---------------------------------------------------------------------------
PUB_STATS = {
    "Shitta_2022": {
        "n_accessions": 169,
        "n_raw_snps": 7930,
        "n_postqc_snps": 1789,
        "qc_call_rate_min": 0.80,
        "qc_maf_min": 0.05,
        "qc_maf_max": 0.95,
        "mean_MAF": 0.22,
        "mean_Ho": 0.15,
        "Ho_sem": 0.002,
        "mean_PIC": None,
        "mean_n_effective_alleles": 1.61,
        "F_IS": None,
        "F_ST_pairwise_low": 0.14,
        "F_ST_pairwise_high": 0.39,
        "amova_among_pops_pct": 13.0,
        "amova_within_pops_pct": 87.0,
        "admixture_optimal_K": 3,
        "K_criterion": "Evanno delta K",
        "LD_halfdecay_kb": None,
        "LDNe": None,
        "source": "Sci Rep 12:4437 (PMC8924269); Tables 2-4 + Results",
    },
    "Olomitutu_2022": {
        "n_accessions": 195,
        "n_raw_snps": 5416,
        "n_postqc_snps": 2491,
        "qc_call_rate_min": 0.70,
        "qc_maf_min": 0.01,
        "qc_maf_max": None,
        "mean_MAF": 0.16,
        "mean_Ho": 0.15,
        "Ho_sem": None,
        "mean_PIC": None,
        "mean_n_effective_alleles": None,
        "F_IS": None,
        "F_ST_pairwise_low": None,
        "F_ST_pairwise_high": None,
        "amova_among_pops_pct": None,
        "amova_within_pops_pct": None,
        "admixture_optimal_K": None,
        "K_criterion": None,
        "LD_halfdecay_kb": None,
        "LDNe": None,
        "source": "Genes 13:2350 (PMC9777823); Results 3.2",
    },
}


# ---------------------------------------------------------------------------
# Recompute our panel-side stats on the same QC subset reported in the
# manuscript headline. Drawing from results/01_qc_pca_power/tables/marker_qc.csv
# (which carries the full 3,204-marker pre-QC table); applying call
# rate >= 0.90 AND MAF >= 0.05 reproduces the 1,625-marker post-QC stats
# the rest of the pipeline uses.
# ---------------------------------------------------------------------------
def this_panel_stats() -> dict:
    mqc = pd.read_csv(MARKER_QC_CSV)
    qc = mqc.loc[(mqc["call_rate"] >= CALL_RATE_MIN) & (mqc["maf"] >= MAF_MIN)]
    p = qc["maf"]
    ne = 1.0 / (p ** 2 + (1 - p) ** 2)
    return {
        "n_accessions": 95,
        "n_raw_snps": int(len(mqc)),
        "n_postqc_snps": int(len(qc)),
        "qc_call_rate_min": CALL_RATE_MIN,
        "qc_maf_min": MAF_MIN,
        "qc_maf_max": None,
        "mean_MAF": float(p.mean()),
        "mean_Ho": float(qc["het"].mean()),
        "Ho_sem": None,
        "mean_PIC": float(qc["pic"].mean()),
        "mean_n_effective_alleles": float(ne.mean()),
        # Panel-wide F_IS from script 22 (Weir-Cockerham); locked in
        # because the per-locus table is in results/22_fis_inbreeding.
        "F_IS": 0.011,
        # Between PCA-derived clusters 1 and 2 from script 10.
        "F_ST_pairwise_low": 0.165,
        "F_ST_pairwise_high": 0.165,
        # From script 02 AMOVA: Phi_ST = 0.273; among-population fraction =
        # 100 * Phi_ST = 27.3.
        "amova_among_pops_pct": 27.3,
        "amova_within_pops_pct": 72.7,
        # From script 09 ADMIXTURE: CV minimum at K=3; K=2 = cleanest
        # visual split. Reported in the manuscript as K=2 / K=3 jointly.
        "admixture_optimal_K": 3,
        "K_criterion": "ADMIXTURE 5-fold CV (CV-error minimum)",
        # From script 11 LD decay (Hill-Weir fit).
        "LD_halfdecay_kb": 73.85,
        "LDNe": 659.0,
        "source": "this panel (Narh-Madey, scripts 01-13)",
    }


# ---------------------------------------------------------------------------
# Wide comparison table
# ---------------------------------------------------------------------------
LABEL_ORDER = [
    ("n_accessions", "n accessions"),
    ("n_raw_snps", "Raw SNPs"),
    ("n_postqc_snps", "Post-QC SNPs"),
    ("qc_call_rate_min", "QC: call rate floor"),
    ("qc_maf_min", "QC: MAF floor"),
    ("mean_MAF", "Mean MAF"),
    ("mean_Ho", "Mean observed heterozygosity"),
    ("mean_PIC", "Mean PIC"),
    ("mean_n_effective_alleles", "Mean n_e (effective alleles per locus, Kimura-Crow)"),
    ("F_IS", "Wright F_IS"),
    ("F_ST_pairwise_low", "Pairwise F_ST (between subpops) -- low"),
    ("F_ST_pairwise_high", "Pairwise F_ST (between subpops) -- high"),
    ("amova_among_pops_pct", "AMOVA among-populations %"),
    ("admixture_optimal_K", "ADMIXTURE optimal K"),
    ("LD_halfdecay_kb", "LD half-decay (kb)"),
    ("LDNe", "LD-based effective population size (Ne)"),
]


def build_wide_table(this: dict) -> pd.DataFrame:
    rows = []
    for key, label in LABEL_ORDER:
        rows.append({
            "statistic": label,
            "this_panel_n95": this[key],
            "Shitta_2022_n169": PUB_STATS["Shitta_2022"][key],
            "Olomitutu_2022_n195": PUB_STATS["Olomitutu_2022"][key],
        })
    df = pd.DataFrame(rows)
    return df


def build_qc_comparison(this: dict) -> pd.DataFrame:
    rows = []
    for panel, source_key, stats in (
        ("This panel (n=95)", None, this),
        ("Shitta 2022 (n=169)", "Shitta_2022", PUB_STATS["Shitta_2022"]),
        ("Olomitutu 2022 (n=195)", "Olomitutu_2022", PUB_STATS["Olomitutu_2022"]),
    ):
        rows.append({
            "panel": panel,
            "n_accessions": stats["n_accessions"],
            "n_raw_snps": stats["n_raw_snps"],
            "n_postqc_snps": stats["n_postqc_snps"],
            "qc_call_rate_min": stats["qc_call_rate_min"],
            "qc_maf_min": stats["qc_maf_min"],
            "qc_maf_max": stats["qc_maf_max"],
            "post_qc_retention_pct": 100.0 * stats["n_postqc_snps"] / stats["n_raw_snps"],
            "source": stats["source"],
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------
FOREST_STATS = [
    ("mean_MAF",            "Mean MAF"),
    ("mean_Ho",             "Mean H_o"),
    ("mean_n_effective_alleles", "Mean n_e (alleles/locus)"),
    ("amova_among_pops_pct", "AMOVA among-populations %"),
    ("F_ST_pairwise_low",    "Pairwise F_ST (low end)"),
]


def plot_forest(this: dict, out_path: Path) -> None:
    fig, axes = plt.subplots(1, len(FOREST_STATS), figsize=(13, 4.0),
                             sharey=False)

    panel_labels = ["This panel\n(n=95)", "Shitta 2022\n(n=169)", "Olomitutu 2022\n(n=195)"]
    panel_colours = [WONG["vermillion"], WONG["blue"], WONG["green"]]

    for ax, (key, label) in zip(axes, FOREST_STATS):
        vals = [this[key],
                PUB_STATS["Shitta_2022"][key],
                PUB_STATS["Olomitutu_2022"][key]]
        for y, (val, colour, lab) in enumerate(zip(vals, panel_colours, panel_labels)):
            if val is None:
                ax.text(0.5, y, "n.r.", ha="center", va="center",
                        fontsize=10, color="grey",
                        transform=ax.get_yaxis_transform())
                continue
            ax.scatter([val], [y], color=colour, s=85, edgecolor="black",
                       linewidth=0.6, zorder=3)
            ax.text(val, y + 0.18, f"{val:.3g}", ha="center", va="bottom",
                    fontsize=9)

        ax.set_yticks(range(len(panel_labels)))
        ax.set_yticklabels(panel_labels, fontsize=9)
        ax.set_title(label, fontsize=10)
        ax.set_xlim(left=0)
        ax.invert_yaxis()
        ax.grid(axis="x", linewidth=0.4, alpha=0.4)

    fig.suptitle(
        "Cross-panel placement of headline diversity statistics\n"
        "(this n=95 panel vs Shitta 2022 n=169 vs Olomitutu 2022 n=195;"
        " each panel under its own QC)",
        fontsize=11,
    )
    fig.text(0.5, 0.01,
             "n.r. = not reported in the source publication. "
             "Olomitutu et al. 2022 (Genes 13:2350) report only mean MAF "
             "and mean Ho; PIC, n_e, AMOVA and pairwise F_ST are absent "
             "from the published tables and the supplementary marker matrix "
             "is not available for re-derivation. "
             "Shitta et al. 2022 (Sci Rep 12:4437) does not report n_e "
             "or Ho SEM. Cells marked n.r. reflect this data-availability "
             "asymmetry, not missing analysis on our side.",
             ha="center", va="bottom", fontsize=7.5, color="grey",
             wrap=True)
    fig.tight_layout(rect=[0, 0.10, 1, 0.94])
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("[load] computing this-panel post-QC stats")
    this = this_panel_stats()
    print(f"  this panel post-QC: n = {this['n_accessions']}, markers = {this['n_postqc_snps']:,}")
    print(f"  mean MAF = {this['mean_MAF']:.3f}, mean Ho = {this['mean_Ho']:.3f}, mean PIC = {this['mean_PIC']:.3f}")
    print(f"  mean n_e (effective alleles per locus) = {this['mean_n_effective_alleles']:.3f}")

    print("\n[table] wide cross-panel comparison")
    wide = build_wide_table(this)
    wide.to_csv(TAB / "per_panel_diversity_stats.csv", index=False)
    print(wide.to_string(index=False))

    print("\n[table] QC and marker-retention comparison")
    qc = build_qc_comparison(this)
    qc.to_csv(TAB / "cross_panel_qc_comparison.csv", index=False)
    print(qc.to_string(index=False))

    print("\n[plot] cross-panel forest")
    plot_forest(this, FIG / "fig_cross_panel_forest")

    print("\n=== manuscript-text summary ===")
    print(f"  this panel mean MAF                       : {this['mean_MAF']:.3f}")
    print(f"   vs Shitta 2022 mean MAF                  : {PUB_STATS['Shitta_2022']['mean_MAF']}")
    print(f"   vs Olomitutu 2022 mean MAF               : {PUB_STATS['Olomitutu_2022']['mean_MAF']}")
    print(f"  this panel mean Ho                        : {this['mean_Ho']:.3f}")
    print(f"   vs Shitta 2022 mean Ho                   : {PUB_STATS['Shitta_2022']['mean_Ho']}")
    print(f"   vs Olomitutu 2022 mean Ho                : {PUB_STATS['Olomitutu_2022']['mean_Ho']}")
    print(f"  this panel mean n_e                       : {this['mean_n_effective_alleles']:.3f}")
    print(f"   vs Shitta 2022 mean n_e                  : {PUB_STATS['Shitta_2022']['mean_n_effective_alleles']}")
    print(f"  this panel AMOVA among-pops %             : {this['amova_among_pops_pct']:.1f}")
    print(f"   vs Shitta 2022 AMOVA among-pops %        : {PUB_STATS['Shitta_2022']['amova_among_pops_pct']}")
    print(f"  this panel pairwise F_ST                  : {this['F_ST_pairwise_low']:.3f}")
    print(f"   vs Shitta 2022 F_ST range                : "
          f"[{PUB_STATS['Shitta_2022']['F_ST_pairwise_low']}, "
          f"{PUB_STATS['Shitta_2022']['F_ST_pairwise_high']}]")


if __name__ == "__main__":
    main()
