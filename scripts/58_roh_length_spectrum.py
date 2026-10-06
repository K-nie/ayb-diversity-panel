#!/usr/bin/env python3
"""
F_ROH segment-length spectrum and recent-selfing-generation inference.

Builds on script 34 (ROH calling). The panel-level F_ROH = 0.256 reported
in script 34 resolves the per-locus F_IS = 0.011 paradox at the genome-wide
scale; what it does not say is whether the autozygosity load on each
accession is ancestral (background relatedness, indexed by short ROH
segments) or recent (within-pedigree inbreeding from within the last few
selfing generations, indexed by long segments).

This script stratifies each ROH segment by length into four classes
following the Kardos et al. (2018) framework as adapted for the panel-
specific length distribution (chromosomal-scale ROH on many accessions):

  short        0.5 to 1 Mb     ancestral / drift autozygosity
  intermediate 1 to 5 Mb       mid-generation pedigree autozygosity
  long         5 to 20 Mb      recent autozygosity (within ~5-10 gen)
  very_long    > 20 Mb         very recent, near-chromosome-scale ROH

Per-accession F_ROH within each class is reported. The recent effective
selfing generations per accession are estimated from the mean length of
long-class ROH under the standard ROH-coalescent expectation
L_expected ~ 100 / (2g) cM, converted at the cowpea proxy map rate of
1.06 cM/Mb (Lonardi 2019), giving g_eff = 47.17 / mean_long_length_Mb.

Outputs (results/58_roh_length_spectrum/)
-----------------------------------------
tables/
    per_accession_roh_class_spectrum.csv   per-accession F_ROH per class + g_eff
    panel_summary.csv                       panel-level means per class
figures/
    fig_roh_length_class_spectrum.png/.pdf  stacked-bar per accession
    fig_geff_distribution.png/.pdf          g_eff histogram across the panel

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
ROH_CALLS = ROOT / "results" / "34_roh" / "tables" / "roh_calls.csv"
F_ROH_SUMM = ROOT / "results" / "34_roh" / "tables" / "f_roh_per_sample.csv"
# Cluster labels from the PCA / K-means analysis (script 01) so the
# stacked-bar figure can annotate the Cluster-2 outlier subgroup.
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "58_roh_length_spectrum"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# Reference-genome span used as F_ROH denominator. 11 AYB pseudo-chromosomes
# total 649.8 Mb (Shorinola et al. 2024 ENA OY731398-OY731408).
GENOME_BP = 649_800_000

# Segment-length class cutoffs, in bp. The four classes follow Kardos 2018
# but with an added "very long" bin above 20 Mb because many AYB accessions
# in this panel carry chromosome-scale ROH (TSs156A, TSs361, TSs358 etc.)
# that would otherwise saturate a single "long" class.
CLASS_EDGES_BP = {
    "short":        (500_000,   1_000_000),
    "intermediate": (1_000_000, 5_000_000),
    "long":         (5_000_000, 20_000_000),
    "very_long":    (20_000_000, np.inf),
}
CLASS_ORDER = ["short", "intermediate", "long", "very_long"]
CLASS_LABEL = {
    "short":        "0.5-1 Mb (ancestral)",
    "intermediate": "1-5 Mb (mid-generation)",
    "long":         "5-20 Mb (recent)",
    "very_long":    ">20 Mb (very recent)",
}

# Cowpea genetic map rate from Lonardi et al. 2019, used as proxy for the
# AYB recombination scale. Under L_expected ~ 100 / (2g) cM and 1.06 cM/Mb,
# the constant relating mean long-segment length to g_eff is 100/(2*1.06) =
# 47.17 Mb-generations.
CM_PER_MB = 1.06
KARDOS_C = 100.0 / (2.0 * CM_PER_MB)  # = 47.17


def classify_segment(length_bp: float) -> str:
    """Return the class label for a single ROH-segment length in bp."""
    for cls in CLASS_ORDER:
        lo, hi = CLASS_EDGES_BP[cls]
        if lo <= length_bp < hi:
            return cls
    # Should not reach here given the inf upper bound on very_long, but
    # the explicit return guards against any unforeseen length.
    return "very_long"


def per_accession_spectrum(roh: pd.DataFrame, all_samples: list[str]) -> pd.DataFrame:
    """Build the per-accession F_ROH-by-class table.

    roh         : columns sample, length_bp (and others, ignored here)
    all_samples : full panel sample list, so accessions with zero ROH
                  segments are still represented as all-zero rows
    """
    roh = roh.copy()
    roh["class"] = roh["length_bp"].apply(classify_segment)

    # Pre-seed an empty record per panel accession so zero-ROH accessions
    # (the outbred tail) carry through with explicit zeros, rather than
    # silently dropping out of the panel-level statistics.
    records: dict[str, dict] = {}
    for sample in all_samples:
        rec = {"sample": sample}
        for cls in CLASS_ORDER:
            rec[f"n_seg_{cls}"] = 0
            rec[f"sum_bp_{cls}"] = 0.0
            rec[f"F_ROH_{cls}"] = 0.0
        rec["mean_long_segment_Mb"] = np.nan
        rec["g_eff_from_long"] = np.nan
        rec["mean_verylong_segment_Mb"] = np.nan
        rec["g_eff_from_verylong"] = np.nan
        records[sample] = rec

    for sample, sub in roh.groupby("sample", sort=False):
        rec = records[sample]
        for cls in CLASS_ORDER:
            class_sub = sub.loc[sub["class"] == cls]
            n_seg = len(class_sub)
            sum_bp = float(class_sub["length_bp"].sum())
            rec[f"n_seg_{cls}"] = n_seg
            rec[f"sum_bp_{cls}"] = sum_bp
            rec[f"F_ROH_{cls}"] = sum_bp / GENOME_BP

        # Recent-selfing-generation estimate from the long-class mean length.
        # If the accession carries no long-class segment, leave g_eff undefined.
        long_sub = sub.loc[sub["class"] == "long"]
        if len(long_sub) > 0:
            mean_long_mb = float(long_sub["length_bp"].mean()) / 1e6
            rec["mean_long_segment_Mb"] = mean_long_mb
            rec["g_eff_from_long"] = KARDOS_C / mean_long_mb

        # Very-long-class mean and g_eff for accessions whose autozygosity
        # is dominated by chromosome-scale ROH. Reported alongside the long-
        # class estimate so the breeder can see both regimes.
        vl_sub = sub.loc[sub["class"] == "very_long"]
        if len(vl_sub) > 0:
            mean_vl_mb = float(vl_sub["length_bp"].mean()) / 1e6
            rec["mean_verylong_segment_Mb"] = mean_vl_mb
            rec["g_eff_from_verylong"] = KARDOS_C / mean_vl_mb

    rows = []
    for sample in all_samples:
        rec = records[sample]
        rec["F_ROH_total"] = sum(rec[f"F_ROH_{c}"] for c in CLASS_ORDER)
        rows.append(rec)

    return pd.DataFrame(rows)


def panel_summary(spectrum: pd.DataFrame) -> pd.DataFrame:
    """Compute panel-level summary numbers per class."""
    rows = []
    for cls in CLASS_ORDER:
        col = f"F_ROH_{cls}"
        n_col = f"n_seg_{cls}"
        rows.append({
            "class": cls,
            "label": CLASS_LABEL[cls],
            "mean_F_ROH": float(spectrum[col].mean()),
            "median_F_ROH": float(spectrum[col].median()),
            "fraction_panel_with_class": float((spectrum[n_col] > 0).mean()),
            "mean_segments_per_accession": float(spectrum[n_col].mean()),
            "total_segments_in_panel": int(spectrum[n_col].sum()),
        })
    return pd.DataFrame(rows)


def plot_stacked_bar(spectrum: pd.DataFrame, out_path: Path) -> None:
    """Stacked bar per accession, ordered by F_ROH_total descending."""
    s = spectrum.sort_values("F_ROH_total", ascending=False).reset_index(drop=True)
    n = len(s)

    # Class colours: short = pale blue (background), long = vermillion (signal).
    colours = {
        "short":        WONG["skyblue"],
        "intermediate": WONG["blue"],
        "long":         WONG["orange"],
        "very_long":    WONG["vermillion"],
    }

    fig, ax = plt.subplots(figsize=(13, 5.0))
    bottom = np.zeros(n)
    x = np.arange(n)

    for cls in CLASS_ORDER:
        vals = s[f"F_ROH_{cls}"].values
        ax.bar(x, vals, bottom=bottom, color=colours[cls],
               edgecolor="white", linewidth=0.2,
               label=CLASS_LABEL[cls])
        bottom = bottom + vals

    ax.set_xticks(x)
    ax.set_xticklabels(s["sample"], rotation=90, fontsize=6)
    ax.set_ylabel("F_ROH (fraction of AYB-anchored genome)")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, max(1.0, bottom.max() * 1.05))
    ax.set_title("Per-accession F_ROH segment-length spectrum (95 AYB lines)")
    ax.legend(loc="upper right", title="ROH length class")
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def plot_geff_distribution(spectrum: pd.DataFrame, out_path: Path) -> None:
    """g_eff histogram across the panel, with the panel median annotated."""
    geff_long = spectrum["g_eff_from_long"].dropna()
    geff_vl = spectrum["g_eff_from_verylong"].dropna()

    fig, ax = plt.subplots(figsize=(7.5, 4.2))

    # Stacked log-axis histogram so the very-recent tail is visible alongside
    # the longer-history mid-generation accessions.
    bins = np.logspace(np.log10(0.3), np.log10(100), 30)
    ax.hist(geff_long, bins=bins, color=WONG["orange"], alpha=0.75,
            edgecolor="white", linewidth=0.4, label="from long (5-20 Mb) segments")
    ax.hist(geff_vl, bins=bins, color=WONG["vermillion"], alpha=0.55,
            edgecolor="white", linewidth=0.4, label="from very-long (>20 Mb) segments")

    ax.set_xscale("log")
    ax.set_xlabel("Estimated recent selfing generations (g_eff)")
    ax.set_ylabel("Number of accessions")
    ax.set_title("Recent effective selfing generations per accession")

    # Reference markers: 1, 5, 10 generations.
    for g_ref, lab in [(1.0, "1 gen"), (5.0, "5 gen"), (10.0, "10 gen")]:
        ax.axvline(g_ref, color="grey", linestyle="--", linewidth=0.7, alpha=0.5)
        ax.text(g_ref, ax.get_ylim()[1] * 0.95, lab, rotation=90,
                fontsize=8, color="grey", va="top", ha="right")

    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    print(f"reading ROH calls from {ROH_CALLS}")
    roh = pd.read_csv(ROH_CALLS)
    print(f"  loaded {len(roh)} segments across {roh['sample'].nunique()} accessions with at least 1 ROH")

    # The script 34 summary table is the authoritative panel sample list;
    # 6 accessions carried zero ROH segments (TSs15, TSs23, TSs294, TSs325,
    # TSs63A, TSs82) and would otherwise drop silently here.
    summ = pd.read_csv(F_ROH_SUMM)
    panel_samples = summ["sample"].tolist()
    print(f"  full panel n = {len(panel_samples)} accessions (6 carry zero ROH segments)")

    # Sanity: the input from script 34 already applies the 500 kb floor, but
    # we re-enforce it here so the class boundaries hold regardless of
    # upstream changes.
    n_before = len(roh)
    roh = roh.loc[roh["length_bp"] >= 500_000].copy()
    n_after = len(roh)
    if n_after < n_before:
        print(f"  dropped {n_before - n_after} segments < 500 kb (below class floor)")

    print("building per-accession class spectrum")
    spectrum = per_accession_spectrum(roh, panel_samples)
    spectrum_path = TAB / "per_accession_roh_class_spectrum.csv"
    spectrum.to_csv(spectrum_path, index=False)
    print(f"  wrote {spectrum_path}")

    print("computing panel-level summary")
    summary = panel_summary(spectrum)
    summary_path = TAB / "panel_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"  wrote {summary_path}")
    print(summary.to_string(index=False))

    print("rendering stacked-bar figure")
    plot_stacked_bar(spectrum, FIG / "fig_roh_length_class_spectrum")
    print("rendering g_eff distribution figure")
    plot_geff_distribution(spectrum, FIG / "fig_geff_distribution")

    # Quick text summary for the manuscript.
    n_panel = len(spectrum)
    fr_recent = (spectrum["F_ROH_long"] + spectrum["F_ROH_very_long"]).mean()
    fr_ancestral = spectrum["F_ROH_short"].mean()
    n_with_vl = int((spectrum["n_seg_very_long"] > 0).sum())
    g_long_med = float(spectrum["g_eff_from_long"].median(skipna=True))
    g_vl_med = float(spectrum["g_eff_from_verylong"].median(skipna=True))

    print("\n=== manuscript-text summary ===")
    print(f"panel n                       : {n_panel}")
    print(f"mean F_ROH ancestral (<1 Mb)  : {fr_ancestral:.3f}")
    print(f"mean F_ROH recent (>5 Mb)     : {fr_recent:.3f}")
    print(f"accessions with >20 Mb ROH    : {n_with_vl} ({100 * n_with_vl / n_panel:.1f}%)")
    print(f"median g_eff from long class  : {g_long_med:.2f} generations")
    print(f"median g_eff from very-long   : {g_vl_med:.2f} generations")


if __name__ == "__main__":
    main()
