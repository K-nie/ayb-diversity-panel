"""
87 -- Panel-wide triangular LD r^2 heatmap across all 11 S. stenocarpa
pseudo-chromosomes (Paper A1 supplementary).

Companion to script 81 which focuses on Ss05 / Ss10 / Ss04 (the three
chromosomes carrying the load-bearing A1 + A2 signals). Script 87 extends
the same draw_triangle pipeline to the full 11-chromosome panel so the
reader can place the per-chromosome LD architecture against the panel's
genome-wide background.

Two outputs:
  fig_panel_wide_ld_11chrs.png/.pdf  -- 11-row stack, one triangle per
                                          chromosome, shared colourbar.
  fig_genome_wide_ld_concat.png/.pdf -- single triangle across the
                                          concatenated genome, with per-
                                          chromosome bands shaded along
                                          the x axis so the reader can see
                                          which off-diagonal blocks sit on
                                          which pseudo-chromosome.

Author: Benjamin Narh-Madey
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon, Rectangle
from matplotlib.collections import PatchCollection

import sys
sys.path.insert(0, str(Path(__file__).parent))
from _figstyle import apply, WONG, publishable_axes
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = PROJ / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = PROJ / "refs/ayb_genome/ayb_marker_anchoring.csv"
OUT = PROJ / "results/87_panel_wide_ld_heatmap"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05

CHRS = [f"Ss{i:02d}" for i in range(1, 12)]
CALLOUTS = {
    "Ss04": [("MATE cluster", 67_792_286)],
    "Ss05": [("PCAdapt + DAPC outliers", 10_900_000)],
    "Ss10": [("ALMT4_2", 15_394_673)],
}

LD_CMAP = LinearSegmentedColormap.from_list(
    # Haploview-style white-to-dark-red, but with a faint cream base so
    # r^2 values near the panel median (~0.012) still render visibly.
    "ld",
    ["#FFF3E0", "#FFD8A8", "#F4A261", "#E76F51", "#C0392B", "#7B241C"],
    N=256,
)


def load_dosage() -> pd.DataFrame:
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in META]
    ref_alt = hm["alleles"].str.split("/", expand=True)
    ref_alt.columns = ["ref", "alt"]
    dosage = np.full((len(hm), len(sample_cols)), np.nan)
    calls = hm[sample_cols].astype(str)
    for i in range(len(hm)):
        r_, a_ = ref_alt.iloc[i]
        arr = calls.iloc[i].values
        dosage[i, arr == r_ + r_] = 0.0
        dosage[i, (arr == r_ + a_) | (arr == a_ + r_)] = 1.0
        dosage[i, arr == a_ + a_] = 2.0
    dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
    cr_m = dosage.notna().sum(axis=1) / dosage.shape[1]
    maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2,
                     1 - dosage.mean(axis=1, skipna=True) / 2)
    cr_s = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_m = ((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values
    keep_s = (cr_s >= SAMP_CR).values
    return dosage.loc[keep_m, keep_s]


def pairwise_r2(X: np.ndarray) -> np.ndarray:
    Xs = X - np.nanmean(X, axis=0)
    Xs = np.where(np.isnan(Xs), 0.0, Xs)
    norms = np.linalg.norm(Xs, axis=0)
    norms[norms == 0] = 1.0
    Xn = Xs / norms
    r = Xn.T @ Xn
    return r ** 2


def draw_triangle(r2: np.ndarray, positions: np.ndarray, ax,
                    vmax: float = 0.5, x_offset_mb: float = 0.0,
                    half_width_floor_mb: float = 0.02,
                    gamma: float = 0.40) -> None:
    """Render the upper triangle of `r2` as rotated diamonds on the
    horizontal axis. Pair (i, j) sits at x = midpoint(pos_i, pos_j) +
    x_offset, y = half-distance. Uses PowerNorm with gamma < 1 so the
    low end of the r^2 range (panel median ~ 0.012) renders visibly."""
    import matplotlib.colors as mcolors
    n = r2.shape[0]
    pos_mb = positions / 1e6 + x_offset_mb
    patches: list[Polygon] = []
    vals: list[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            x = (pos_mb[i] + pos_mb[j]) * 0.5
            y = (pos_mb[j] - pos_mb[i]) * 0.5
            hw = max((pos_mb[j] - pos_mb[i]) * 0.25, half_width_floor_mb)
            patches.append(Polygon([
                (x - hw, y), (x, y - hw), (x + hw, y), (x, y + hw),
            ]))
            vals.append(r2[i, j])
    norm = mcolors.PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)
    pc = PatchCollection(patches, cmap=LD_CMAP, norm=norm, alpha=0.95,
                          edgecolor="none")
    pc.set_array(np.array(vals))
    ax.add_collection(pc)


def build_per_chr_panel(dosage: pd.DataFrame, anc: pd.DataFrame) -> None:
    """11-row stack: one triangle per Ss chromosome."""
    rows = []
    fig, axes = plt.subplots(len(CHRS), 1,
                              figsize=(13.5, 2.8 * len(CHRS)),
                              gridspec_kw={"hspace": 0.55})
    for ax, chrom in zip(axes, CHRS):
        sub = anc[anc["chr_ayb"] == chrom].sort_values("snp_pos_ayb").copy()
        if sub.empty:
            ax.text(0.5, 0.5, f"No anchored markers on {chrom}",
                    ha="center", transform=ax.transAxes, fontsize=10)
            publishable_axes(ax, grid=None)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        markers = sub["rs"].tolist()
        positions = sub["snp_pos_ayb"].values.astype(float)
        X = dosage.loc[markers].T.values.astype(float)
        col_mu = np.nanmean(X, axis=0)
        X = np.where(np.isnan(X), col_mu, X)
        r2 = pairwise_r2(X)
        triu_i, triu_j = np.triu_indices_from(r2, k=1)
        long_df = pd.DataFrame({
            "rs_i": np.array(markers)[triu_i],
            "rs_j": np.array(markers)[triu_j],
            "pos_i": positions[triu_i],
            "pos_j": positions[triu_j],
            "r2": r2[triu_i, triu_j],
        })
        long_df.to_csv(OUT / f"tables/ld_pairs_{chrom}.csv", index=False)
        med_r2 = float(np.median(long_df["r2"]))
        p95 = float(np.percentile(long_df["r2"], 95))
        rows.append({
            "chr": chrom, "n_markers": len(markers),
            "n_pairs": len(long_df), "r2_median": med_r2, "r2_p95": p95,
            "span_mb": float((positions[-1] - positions[0]) / 1e6),
        })

        draw_triangle(r2, positions, ax, vmax=0.5)
        x_pad = max((positions.max() - positions.min()) / 2e6 * 0.03, 0.5)
        ax.set_xlim(positions.min() / 1e6 - x_pad,
                    positions.max() / 1e6 + x_pad)
        ax.set_ylim(0, (positions.max() - positions.min()) / 2e6 * 1.05)
        ax.set_aspect("equal")
        ax.set_ylabel("Pair half-\ndistance (Mb)", fontsize=9)
        ax.set_xlabel(f"{chrom} position (Mb)", fontsize=9.5)
        ax.set_title(f"{chrom}: n = {len(markers)} anchored markers; "
                     f"median r$^2$ = {med_r2:.3f}; "
                     f"95th pct = {p95:.3f}",
                     fontsize=10, pad=6)
        publishable_axes(ax, grid=None)
        for nm, pos in CALLOUTS.get(chrom, []):
            ax.axvline(pos / 1e6, color=WONG["vermillion"], lw=0.9,
                        ls=":", alpha=0.85)
            ax.text(pos / 1e6, ax.get_ylim()[1] * 0.92, nm,
                     ha="center", fontsize=7.5,
                     color=WONG["vermillion"],
                     bbox=dict(facecolor="white", edgecolor="none",
                                 alpha=0.9, pad=1.5))

    pd.DataFrame(rows).to_csv(OUT / "tables/per_chr_summary.csv",
                                index=False)
    print(pd.DataFrame(rows).to_string(index=False))

    import matplotlib.colors as mcolors
    cb_norm = mcolors.PowerNorm(gamma=0.40, vmin=0.0, vmax=0.5)
    sm = plt.cm.ScalarMappable(cmap=LD_CMAP, norm=cb_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.018, pad=0.012)
    cbar.set_label(r"pairwise $r^2$")

    fig.suptitle("Panel-wide triangular LD heatmaps across all 11 "
                  "S. stenocarpa pseudo-chromosomes\n"
                  "AYB-anchored DArTseq markers, n = 95",
                  fontsize=12, y=1.005)
    fig.savefig(OUT / "figures/fig_panel_wide_ld_11chrs.png", dpi=300,
                 bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "figures/fig_panel_wide_ld_11chrs.pdf",
                 bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("\nWrote", OUT / "figures/fig_panel_wide_ld_11chrs.png")


def build_concatenated_panel(dosage: pd.DataFrame, anc: pd.DataFrame) -> None:
    """Single triangle across all 11 chromosomes concatenated head-to-tail,
    with per-chromosome bands shaded along the x axis. Limits LD pairs to
    within-chromosome only (cross-chromosome pairs are biologically
    uninformative for LD and would crowd the diamond pile)."""
    fig, ax = plt.subplots(figsize=(15.5, 6.5))

    cursor_mb = 0.0
    chrom_ranges: list[tuple[str, float, float]] = []
    band_palette = ["#F0F0F0", "#FFFFFF"]
    GAP_MB = 8.0
    for i, chrom in enumerate(CHRS):
        sub = anc[anc["chr_ayb"] == chrom].sort_values("snp_pos_ayb").copy()
        if sub.empty:
            continue
        markers = sub["rs"].tolist()
        positions = sub["snp_pos_ayb"].values.astype(float)
        X = dosage.loc[markers].T.values.astype(float)
        col_mu = np.nanmean(X, axis=0)
        X = np.where(np.isnan(X), col_mu, X)
        r2 = pairwise_r2(X)
        chrom_start_mb = cursor_mb
        draw_triangle(r2, positions, ax, vmax=0.5,
                        x_offset_mb=cursor_mb - positions.min() / 1e6)
        chrom_end_mb = cursor_mb + (positions.max() - positions.min()) / 1e6
        chrom_ranges.append((chrom, chrom_start_mb, chrom_end_mb))
        cursor_mb = chrom_end_mb + GAP_MB

    # Background bands per chromosome so the reader can see boundaries
    y_max = 30.0
    for i, (chrom, x0, x1) in enumerate(chrom_ranges):
        ax.add_patch(Rectangle(
            (x0, 0), x1 - x0, y_max,
            facecolor=band_palette[i % 2], edgecolor="none",
            alpha=0.55, zorder=0,
        ))

    # Per-chromosome callouts
    for chrom, x0, x1 in chrom_ranges:
        for nm, pos in CALLOUTS.get(chrom, []):
            x = x0 + (pos - anc[anc["chr_ayb"] == chrom]["snp_pos_ayb"].min()) / 1e6
            ax.axvline(x, color=WONG["vermillion"], lw=0.9,
                        ls=":", alpha=0.85, ymax=0.7)
            ax.text(x, y_max * 0.62, nm, ha="center",
                     fontsize=7.5, color=WONG["vermillion"], rotation=0,
                     bbox=dict(facecolor="white", edgecolor="none",
                                 alpha=0.95, pad=1.5))

    # x-axis: per-chromosome ticks at centre of each band
    centres = [(x0 + x1) / 2.0 for _, x0, x1 in chrom_ranges]
    labels = [c for c, _, _ in chrom_ranges]
    ax.set_xticks(centres)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlim(0, cursor_mb)
    ax.set_ylim(0, y_max)
    ax.set_aspect("auto")
    ax.set_xlabel("Concatenated S. stenocarpa pseudo-chromosomes "
                   "(AYB-anchored markers)")
    ax.set_ylabel("Pair half-distance (Mb)")
    ax.set_title("Genome-wide triangular LD heatmap across the "
                  "11 S. stenocarpa pseudo-chromosomes "
                  "(within-chromosome pairs only)",
                  fontsize=11)
    publishable_axes(ax, grid=None)

    import matplotlib.colors as mcolors
    cb_norm = mcolors.PowerNorm(gamma=0.40, vmin=0.0, vmax=0.5)
    sm = plt.cm.ScalarMappable(cmap=LD_CMAP, norm=cb_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.020, pad=0.015)
    cbar.set_label(r"pairwise $r^2$")
    fig.tight_layout()
    fig.savefig(OUT / "figures/fig_genome_wide_ld_concat.png", dpi=300,
                 bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "figures/fig_genome_wide_ld_concat.pdf",
                 bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Wrote", OUT / "figures/fig_genome_wide_ld_concat.png")


def main() -> None:
    print("[load] dosage matrix...")
    dosage = load_dosage()
    print(f"[load] dosage shape: {dosage.shape}")
    anc = pd.read_csv(ANCHOR)
    anc = anc[anc["rs"].isin(dosage.index)].copy()
    anc = anc[["rs", "chr_ayb", "snp_pos_ayb"]].dropna()
    print(f"[load] AYB-anchored markers in QC-passed dosage: {len(anc)}")

    print("\n[build] per-chromosome 11-row stack...")
    build_per_chr_panel(dosage, anc)
    print("\n[build] genome-wide concatenated single panel...")
    build_concatenated_panel(dosage, anc)


if __name__ == "__main__":
    main()
