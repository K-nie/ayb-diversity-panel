#!/usr/bin/env python3
"""
Per-chromosome LD (r^2) decay across the AYB-anchored 11 Ss pseudo-chromosomes.

Reads the AYB-anchored per-Ss pair tables already computed in
`results/87_panel_wide_ld_heatmap/tables/ld_pairs_Ss*.csv`, bins by physical
distance within each chromosome, fits the Hill & Weir (1988) curve per
chromosome, and renders an 11-panel small-multiples figure with consistent
axes. Replaces the legacy `fig28_ld_decay_per_chrom` panel that used cowpea
Vu01-Vu11 anchoring (mostly empty after re-anchoring).

Output: results/11_ld_decay/figures/fig28_ld_decay_per_chrom_ayb.png/.pdf
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PAIR_DIR = ROOT / "results/87_panel_wide_ld_heatmap/tables"
OUT_FIG = ROOT / "results/11_ld_decay/figures"
OUT_TAB = ROOT / "results/11_ld_decay/tables"
for d in (OUT_FIG, OUT_TAB):
    d.mkdir(parents=True, exist_ok=True)

MAX_DIST_BP = 5_000_000
N_BINS = 20
CM_PER_MB = 1.06
MORGANS_PER_BP = (CM_PER_MB / 100.0) / 1e6

CHRS = [f"Ss{str(i).zfill(2)}" for i in range(1, 12)]


def hill_weir(d_bp, rho):
    C = rho * MORGANS_PER_BP * d_bp
    return (10.0 + C) / ((2.0 + C) * (11.0 + C))


def per_chr_summary(pairs: pd.DataFrame, chr_name: str) -> dict:
    sub = pairs[pairs["dist_bp"] <= MAX_DIST_BP].copy()
    n_pairs = len(sub)
    if n_pairs < 30:
        return {
            "chr": chr_name, "n_pairs": n_pairs, "rho_hat": np.nan,
            "half_decay_bp": np.nan,
        }
    edges = np.linspace(0, MAX_DIST_BP, N_BINS + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    sub["bin"] = np.digitize(sub["dist_bp"], edges) - 1
    binned = sub.groupby("bin")["r2"].agg(["mean", "count"]).reindex(
        range(N_BINS))
    binned["mid_bp"] = mids
    valid = binned["mean"].notna() & (binned["count"] >= 5)
    if valid.sum() < 5:
        return {
            "chr": chr_name, "n_pairs": n_pairs, "rho_hat": np.nan,
            "half_decay_bp": np.nan, "binned": binned,
        }
    try:
        popt, _ = curve_fit(
            hill_weir,
            binned.loc[valid, "mid_bp"].values,
            binned.loc[valid, "mean"].values,
            p0=[100.0], bounds=(1e-3, 1e6), maxfev=10000,
        )
        rho_hat = float(popt[0])
        d_grid = np.linspace(1, MAX_DIST_BP, 100_000)
        y_grid = hill_weir(d_grid, rho_hat)
        half = float(d_grid[np.argmin(np.abs(y_grid - y_grid[0] / 2.0))])
    except Exception:
        rho_hat = np.nan
        half = np.nan
    return {
        "chr": chr_name, "n_pairs": n_pairs,
        "rho_hat": rho_hat, "half_decay_bp": half,
        "binned": binned,
    }


def main() -> None:
    summaries = []
    fig, axes = plt.subplots(
        3, 4, figsize=(12.0, 8.4),
        constrained_layout=True, sharex=True, sharey=True,
    )
    axes_flat = axes.ravel()

    for ax, chr_name in zip(axes_flat[: len(CHRS)], CHRS):
        path = PAIR_DIR / f"ld_pairs_{chr_name}.csv"
        if not path.exists():
            ax.set_title(f"{chr_name} (no pairs)")
            ax.set_visible(False)
            continue
        df = pd.read_csv(path)
        df["dist_bp"] = (df["pos_j"] - df["pos_i"]).abs().astype(float)
        df = df[df["dist_bp"] <= MAX_DIST_BP]
        s = per_chr_summary(df, chr_name)
        summaries.append({k: v for k, v in s.items() if k != "binned"})

        ax.scatter(
            df["dist_bp"] / 1e3, df["r2"], s=5, alpha=0.10,
            color=WONG["skyblue"], edgecolor="none", rasterized=True,
        )
        if "binned" in s and s["binned"] is not None:
            b = s["binned"]
            ax.plot(
                b["mid_bp"] / 1e3, b["mean"], "o-",
                color=WONG["blue"], markersize=3, linewidth=1.0,
                label="bin mean",
            )
        if np.isfinite(s["rho_hat"]):
            d_plot = np.linspace(1, MAX_DIST_BP, 600)
            ax.plot(
                d_plot / 1e3, hill_weir(d_plot, s["rho_hat"]),
                color=WONG["vermillion"], linewidth=1.4,
                label=f"H-W fit, half = {s['half_decay_bp']/1e3:.0f} kb",
            )
            title = (f"{chr_name}  (n = {s['n_pairs']:,};  "
                     f"half-decay {s['half_decay_bp']/1e3:.0f} kb)")
        else:
            title = f"{chr_name}  (n = {s['n_pairs']:,})"
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, MAX_DIST_BP / 1e3)
        ax.set_ylim(0, 0.6)
        ax.grid(True)
        ax.legend(fontsize=7, loc="upper right")

    for ax in axes_flat[len(CHRS):]:
        ax.set_visible(False)
    for ax in axes[-1, :]:
        ax.set_xlabel("physical distance (kb)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"r$^2$")

    fig.suptitle(
        "LD decay per AYB pseudo-chromosome (Ss01-Ss11, "
        "Hill-Weir fit per chromosome)",
        y=1.02, fontsize=12,
    )

    out_png = OUT_FIG / "fig28_ld_decay_per_chrom_ayb.png"
    out_pdf = OUT_FIG / "fig28_ld_decay_per_chrom_ayb.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(summaries).to_csv(
        OUT_TAB / "ld_decay_per_chrom_ayb_summary.csv", index=False)

    print(f"[fig] wrote {out_png}")
    print(f"[tab] wrote {OUT_TAB / 'ld_decay_per_chrom_ayb_summary.csv'}")
    for s in summaries:
        if np.isfinite(s["half_decay_bp"]):
            print(f"  {s['chr']}: n = {s['n_pairs']:,}, "
                  f"rho-hat = {s['rho_hat']:.1f}, "
                  f"half-decay = {s['half_decay_bp']/1e3:.0f} kb")
        else:
            print(f"  {s['chr']}: n = {s['n_pairs']:,}, fit failed")


if __name__ == "__main__":
    main()
