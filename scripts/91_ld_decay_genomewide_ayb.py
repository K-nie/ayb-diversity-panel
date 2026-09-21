#!/usr/bin/env python3
"""
91 - Genome-wide LD decay and Ne on the AYB-anchored marker set.

Reviewer 1 major point 6: the submitted panel-wide Ne (approx 659) and the
74 kb global half-decay came from a global Hill-Weir fit on the *cowpea*-
anchored SNP subset (176 SNPs with paired Vu positions in the original DArT
report). That is precisely the legacy anchoring the rest of the revision
abandons in favour of the *S. stenocarpa* assembly. The reviewer asked us to
refit the global decay on the AYB-anchored set or drop the Ne claim.

This script pools the per-chromosome AYB-anchored r^2 pair tables already
computed for the panel-wide LD heatmap (script 87) into a single genome-wide
distance-decay curve, fits the Hill & Weir (1988) drift-recombination
expectation, and derives the half-decay distance, the distance at which
r^2 < 0.1, and an Ne estimate under the same 1.06 cM/Mb cowpea map-rate proxy
used in the submission (so the AYB-anchored global value is directly
comparable to the submitted per-chromosome AYB fits and to the legacy cowpea
global value it replaces).

Inputs
  - results/87_panel_wide_ld_heatmap/tables/ld_pairs_Ss01..Ss11.csv
    (columns: rs_i, rs_j, pos_i, pos_j, r2; all within-chromosome pairs)

Outputs (results/91_ld_decay_genomewide_ayb/)
  tables/ld_genomewide_ayb_summary.csv   - rho, Ne, half-decay, r2<0.1, n
  tables/ld_genomewide_ayb_bins.csv      - the binned decay curve
  figures/fig_ld_decay_genomewide_ayb.png/.pdf  - replacement for Figure 3
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _plotstyle import apply, WONG  # noqa: E402
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PAIR_DIR = ROOT / "results" / "87_panel_wide_ld_heatmap" / "tables"
OUT = ROOT / "results" / "91_ld_decay_genomewide_ayb"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

MAX_DIST_BP = 5_000_000
N_BINS = 40
CM_PER_MB = 1.06                       # cowpea map-rate proxy, as in the paper
MORGANS_PER_BP = (CM_PER_MB / 100.0) / 1e6
CHRS = [f"Ss{str(i).zfill(2)}" for i in range(1, 12)]


def hill_weir(d_bp, rho):
    C = rho * MORGANS_PER_BP * d_bp
    return (10.0 + C) / ((2.0 + C) * (11.0 + C))


def main() -> None:
    frames = []
    markers = set()
    for c in CHRS:
        path = PAIR_DIR / f"ld_pairs_{c}.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path)
        df["dist_bp"] = (df["pos_j"] - df["pos_i"]).abs().astype(float)
        df["chr"] = c
        markers.update(df["rs_i"])
        markers.update(df["rs_j"])
        frames.append(df[["chr", "dist_bp", "r2"]])
    pairs = pd.concat(frames, ignore_index=True)
    pairs = pairs[pairs["dist_bp"] <= MAX_DIST_BP].copy()
    n_pairs = len(pairs)
    n_markers = len(markers)
    print(f"[load] {n_pairs:,} within-chromosome pairs <= {MAX_DIST_BP/1e6:.0f} Mb "
          f"across {len(CHRS)} pseudo-chromosomes; {n_markers:,} anchored markers")

    # Bin by physical distance; mean r^2 + 95% CI per bin.
    edges = np.linspace(0, MAX_DIST_BP, N_BINS + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    pairs["bin"] = np.digitize(pairs["dist_bp"], edges) - 1
    binned = (pairs.groupby("bin")["r2"]
              .agg(["mean", "std", "count"])
              .reindex(range(N_BINS)))
    binned["mid_bp"] = mids
    binned["se"] = binned["std"] / np.sqrt(binned["count"].replace(0, np.nan))
    binned.to_csv(TAB / "ld_genomewide_ayb_bins.csv")

    valid = binned["mean"].notna() & (binned["count"] >= 5)
    xfit = binned.loc[valid, "mid_bp"].values
    yfit = binned.loc[valid, "mean"].values
    popt, _ = curve_fit(hill_weir, xfit, yfit, p0=[1500.0],
                        bounds=(1e-3, 1e6), maxfev=20000)
    rho_hat = float(popt[0])
    Ne_hat = rho_hat / 4.0
    d_grid = np.linspace(1, MAX_DIST_BP, 500_000)
    y_grid = hill_weir(d_grid, rho_hat)
    half_bp = float(d_grid[np.argmin(np.abs(y_grid - y_grid[0] / 2.0))])
    r2_01 = float(d_grid[np.argmax(y_grid < 0.1)]) if (y_grid < 0.1).any() else np.nan
    print(f"[fit] genome-wide AYB rho_hat = {rho_hat:.1f}  ->  Ne = {Ne_hat:.0f}")
    print(f"[fit] half-decay = {half_bp/1e3:.1f} kb; r^2<0.1 by {r2_01/1e3:.1f} kb")

    pd.DataFrame([{
        "anchoring": "AYB (S. stenocarpa Ss01-Ss11)",
        "n_pairs": n_pairs,
        "n_markers": n_markers,
        "max_dist_bp": MAX_DIST_BP,
        "cm_per_mb_assumed": CM_PER_MB,
        "rho_hat": round(rho_hat, 1),
        "Ne_estimate": round(Ne_hat, 1),
        "half_decay_bp": round(half_bp, 0),
        "half_decay_kb": round(half_bp / 1e3, 1),
        "dist_r2_below_0.1_bp": round(r2_01, 0),
        "dist_r2_below_0.1_kb": round(r2_01 / 1e3, 1),
    }]).to_csv(TAB / "ld_genomewide_ayb_summary.csv", index=False)

    # Figure (replacement for the cowpea-anchored Figure 3).
    fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
    samp = pairs.sample(min(40_000, n_pairs), random_state=20260920)
    ax.scatter(samp["dist_bp"] / 1e3, samp["r2"], s=3, alpha=0.05,
               color=WONG["skyblue"], edgecolor="none", rasterized=True,
               label=r"per-pair r$^2$")
    ax.errorbar(binned["mid_bp"] / 1e3, binned["mean"],
                yerr=1.96 * binned["se"], fmt="o", color=WONG["blue"],
                markersize=4, elinewidth=0.8, capsize=2,
                label="bin mean ± 95 % CI")
    d_plot = np.linspace(1, MAX_DIST_BP, 1000)
    ax.plot(d_plot / 1e3, hill_weir(d_plot, rho_hat),
            color=WONG["vermillion"], linewidth=1.6,
            label=(f"Hill–Weir fit, $\\hat\\rho$ = {rho_hat:.0f}; "
                   f"half-decay {half_bp/1e3:.0f} kb; N$_e$ ≈ {Ne_hat:.0f}"))
    ax.axvline(half_bp / 1e3, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xlabel("physical distance (kb)")
    ax.set_ylabel(r"r$^2$")
    ax.set_title(f"Genome-wide LD decay across AYB-anchored markers "
                 f"({n_pairs:,} pairs, Ss01–Ss11)")
    ax.set_xlim(0, MAX_DIST_BP / 1e3)
    ax.set_ylim(0, max(0.3, float(binned["mean"].max()) * 1.15))
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True)
    fig.savefig(FIG / "fig_ld_decay_genomewide_ayb.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "fig_ld_decay_genomewide_ayb.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[fig] fig_ld_decay_genomewide_ayb")


if __name__ == "__main__":
    main()
