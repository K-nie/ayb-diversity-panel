#!/usr/bin/env python3
"""
96 - Genome-wide and per-chromosome LD decay / Ne on the AYB-anchored marker
set, recomputed on the DE-DUPLICATED panel with the SAMPLE-SIZE-CORRECTED
Hill-Weir expectation. Reviewer 1 major point 6 (second pass).

Background
----------
Script 91 refit the genome-wide decay on the AYB-anchored markers (replacing the
legacy cowpea-anchored fit) but (a) used the uncorrected Hill & Weir (1988)
drift-recombination expectation, which omits the finite-sample term, and
(b) kept all n = 95 lines, including the near-clonal Cluster-2 trio
(TSs151B / TSs358 / TSs361, pairwise Gij approximately 1.60, F_ROH approximately
0.99; manuscript Table 3). Near-identical lines inflate LD and the sample-size
term is non-negligible at n = 95 (expected r^2 under no LD approximately 1/n).

This script addresses both:
  (i)  fits the sample-size-corrected expectation
       E[r^2] = base * (1 + ((3+C)(12+12C+C^2)) / (n (2+C)(11+C))),
       base = (10+C)/((2+C)(11+C)), C = rho * d(Morgans), n = sample size
       (Hill & Weir 1988; Remington et al. 2001 PNAS 98:11479).
  (ii) collapses the near-clonal trio to a single representative (drop TSs358,
       TSs361; keep TSs151B), giving n = 93. The near-clonal set is the panel's
       own duplicate call at Gij >= 0.85; the five close-pedigree pairs at
       Gij 0.81-0.82 are NOT duplicates and are retained.

Decisions (documented for the README / reviewer response)
---------------------------------------------------------
- The marker set is held fixed at the 1,473 AYB-anchored QC markers used by
  scripts 87/91 (call-rate >= 0.90, MAF >= 0.05 computed on the full panel).
  Removing two lines is not allowed to change the marker panel, so the fit
  isolates the effect of the duplicate lines and the sample-size term from any
  change in the marker set. (Re-applying the MAF filter on n = 93 drops ~29
  markers and moves Ne by < 3 %; reported as a sensitivity note only.)
- Recombination-rate proxy unchanged: 1.06 cM/Mb from the cowpea genetic map,
  so absolute Ne stays comparable to scripts 91/11b. The proxy is an
  acknowledged absolute-scale caveat; it does not affect the relative decay
  profile or the physical half-decay distance.

Inputs
------
- AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv  (genotypes)
- refs/ayb_genome/ayb_marker_anchoring.csv                            (Ss positions)

Outputs (results/96_ld_decay_dedup_sscorr/)
  tables/ld_genomewide_dedup_sscorr_summary.csv
  tables/ld_genomewide_dedup_sscorr_bins.csv
  tables/ld_per_chrom_dedup_sscorr_summary.csv
  figures/fig_ld_decay_genomewide_ayb.{png,pdf}   (replaces script-91 Figure 3)
  figures/fig28_ld_decay_per_chrom_ayb.{png,pdf}  (replaces script-11b Supp S4A)

Author: Benjamin Narh-Madey
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
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
OUT = ROOT / "results" / "96_ld_decay_dedup_sscorr"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = MARK_CR = 0.90
MIN_MAF = 0.05
MAX_DIST_BP = 5_000_000
N_BINS = 40
CM_PER_MB = 1.06
MORGANS_PER_BP = (CM_PER_MB / 100.0) / 1e6
CHRS = [f"Ss{i:02d}" for i in range(1, 12)]
DROP_NEAR_CLONAL = ["TSs358", "TSs361"]   # keep TSs151B as the trio representative


def load_dosage() -> pd.DataFrame:
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sc = [c for c in hm.columns if c not in META]
    ra = hm["alleles"].str.split("/", expand=True); ra.columns = ["ref", "alt"]
    dos = np.full((len(hm), len(sc)), np.nan); calls = hm[sc].astype(str)
    for i in range(len(hm)):
        r_, a_ = ra.iloc[i]; arr = calls.iloc[i].values
        dos[i, arr == r_ + r_] = 0.0
        dos[i, (arr == r_ + a_) | (arr == a_ + r_)] = 1.0
        dos[i, arr == a_ + a_] = 2.0
    dos = pd.DataFrame(dos, index=hm["rs#"].values, columns=sc)
    cr_m = dos.notna().sum(1) / dos.shape[1]
    maf = np.minimum(dos.mean(1, skipna=True) / 2, 1 - dos.mean(1, skipna=True) / 2)
    cr_s = dos.notna().sum(0) / dos.shape[0]
    return dos.loc[((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values, (cr_s >= SAMP_CR).values]


def r2_matrix(X: np.ndarray) -> np.ndarray:
    Xs = X - np.nanmean(X, axis=0)
    Xs = np.where(np.isnan(Xs), 0.0, Xs)
    norms = np.linalg.norm(Xs, axis=0); norms[norms == 0] = 1.0
    Xn = Xs / norms
    return (Xn.T @ Xn) ** 2


def chr_pairs(d: pd.DataFrame, anc: pd.DataFrame, chrom: str):
    sub = anc[anc["chr_ayb"] == chrom].sort_values("snp_pos_ayb")
    mk = [m for m in sub["rs"].tolist() if m in d.index]
    if len(mk) < 2:
        return None
    pos = sub.set_index("rs").loc[mk, "snp_pos_ayb"].values.astype(float)
    X = d.loc[mk].T.values.astype(float)
    mu = np.nanmean(X, axis=0); X = np.where(np.isnan(X), mu, X)
    R = r2_matrix(X); ii, jj = np.triu_indices_from(R, 1)
    dist = np.abs(pos[ii] - pos[jj]); keep = dist <= MAX_DIST_BP
    return dist[keep], R[ii, jj][keep], len(mk)


def hw_corr(d_bp, rho, n):
    C = rho * MORGANS_PER_BP * d_bp
    base = (10.0 + C) / ((2.0 + C) * (11.0 + C))
    corr = 1.0 + ((3.0 + C) * (12.0 + 12.0 * C + C ** 2)) / (n * (2.0 + C) * (11.0 + C))
    return base * corr


def fit_curve(dist, r2, n):
    edges = np.linspace(0, MAX_DIST_BP, N_BINS + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(dist, edges) - 1, 0, N_BINS - 1)
    means = np.full(N_BINS, np.nan); stds = np.full(N_BINS, np.nan); cnts = np.zeros(N_BINS, int)
    for b in range(N_BINS):
        m = idx == b
        cnts[b] = int(m.sum())
        if m.sum() > 0:
            means[b] = r2[m].mean(); stds[b] = r2[m].std()
    valid = ~np.isnan(means) & (cnts >= 5)
    popt, _ = curve_fit(lambda d, rho: hw_corr(d, rho, n), mids[valid], means[valid],
                        p0=[1500.0], bounds=(1e-3, 1e6), maxfev=40000)
    rho = float(popt[0]); Ne = rho / 4.0
    grid = np.linspace(1, MAX_DIST_BP, 500_000)
    y = hw_corr(grid, rho, n)
    half = float(grid[np.argmin(np.abs(y - y[0] / 2.0))])
    d01 = float(grid[np.argmax(y < 0.1)]) if (y < 0.1).any() else np.nan
    return dict(rho=rho, Ne=Ne, half=half, d01=d01,
                mids=mids, means=means, stds=stds, cnts=cnts, valid=valid)


def main() -> None:
    anc = pd.read_csv(ANCHOR)[["rs", "chr_ayb", "snp_pos_ayb"]].dropna()
    dos = load_dosage()
    full_n = dos.shape[1]
    ded = dos[[s for s in dos.columns if s not in DROP_NEAR_CLONAL]]
    n = ded.shape[1]
    print(f"[panel] full n={full_n}; dedup n={n} (dropped {DROP_NEAR_CLONAL}, kept TSs151B)")

    # ---- genome-wide ----
    dists, r2s, markers = [], [], set()
    for c in CHRS:
        res = chr_pairs(ded, anc, c)
        if res is None:
            continue
        dd, rr, nmk = res
        dists.append(dd); r2s.append(rr)
        sub = anc[anc["chr_ayb"] == c]
        markers.update(m for m in sub["rs"].tolist() if m in ded.index)
    dist = np.concatenate(dists); r2 = np.concatenate(r2s)
    n_pairs = len(dist); n_markers = len(markers)
    gw = fit_curve(dist, r2, n)
    print(f"[gw] rho={gw['rho']:.1f} Ne={gw['Ne']:.1f} half={gw['half']/1e3:.1f}kb "
          f"r2<0.1@{gw['d01']/1e3:.1f}kb pairs={n_pairs} markers={n_markers}")

    pd.DataFrame([{
        "anchoring": "AYB (S. stenocarpa Ss01-Ss11)", "panel": f"dedup n={n}",
        "n_pairs": n_pairs, "n_markers": n_markers,
        "sample_size_correction": "Hill-Weir 1988 / Remington 2001",
        "cm_per_mb_assumed": CM_PER_MB,
        "rho_hat": round(gw["rho"], 1), "Ne_estimate": round(gw["Ne"], 1),
        "half_decay_bp": round(gw["half"], 0), "half_decay_kb": round(gw["half"] / 1e3, 1),
        "dist_r2_below_0.1_bp": round(gw["d01"], 0),
        "dist_r2_below_0.1_kb": round(gw["d01"] / 1e3, 1),
    }]).to_csv(TAB / "ld_genomewide_dedup_sscorr_summary.csv", index=False)

    pd.DataFrame({"bin": range(N_BINS), "mid_bp": gw["mids"], "mean": gw["means"],
                  "std": gw["stds"], "count": gw["cnts"]}
                 ).to_csv(TAB / "ld_genomewide_dedup_sscorr_bins.csv", index=False)

    # ---- per-chromosome ----
    rows = []
    per_chr_fits = {}
    for c in CHRS:
        res = chr_pairs(ded, anc, c)
        if res is None:
            continue
        dd, rr, nmk = res
        f = fit_curve(dd, rr, n)
        per_chr_fits[c] = (dd, rr, f)
        rows.append({"chr": c, "n_markers": nmk, "n_pairs": len(dd),
                     "rho_hat": round(f["rho"], 1), "Ne_estimate": round(f["Ne"], 1),
                     "half_decay_kb": round(f["half"] / 1e3, 1),
                     "dist_r2_below_0.1_kb": round(f["d01"] / 1e3, 1)})
    per = pd.DataFrame(rows)
    per.to_csv(TAB / "ld_per_chrom_dedup_sscorr_summary.csv", index=False)
    print("[per-chr] Ne range %.0f-%.0f, half %.0f-%.0f kb, rho %.0f-%.0f, median Ne %.0f" % (
        per["Ne_estimate"].min(), per["Ne_estimate"].max(),
        per["half_decay_kb"].min(), per["half_decay_kb"].max(),
        per["rho_hat"].min(), per["rho_hat"].max(), per["Ne_estimate"].median()))

    # ---- Figure 3 (genome-wide) ----
    fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
    rng = np.random.default_rng(20260929)
    if n_pairs > 40000:
        sel = rng.choice(n_pairs, 40000, replace=False)
    else:
        sel = np.arange(n_pairs)
    ax.scatter(dist[sel] / 1e3, r2[sel], s=3, alpha=0.05, color=WONG["skyblue"],
               edgecolor="none", rasterized=True, label=r"per-pair r$^2$")
    se = gw["stds"] / np.sqrt(np.where(gw["cnts"] > 0, gw["cnts"], np.nan))
    ax.errorbar(gw["mids"] / 1e3, gw["means"], yerr=1.96 * se, fmt="o",
                color=WONG["blue"], markersize=4, elinewidth=0.8, capsize=2,
                label="bin mean ± 95 % CI")
    dpl = np.linspace(1, MAX_DIST_BP, 1000)
    ax.plot(dpl / 1e3, hw_corr(dpl, gw["rho"], n), color=WONG["vermillion"], linewidth=1.6,
            label=(f"Hill–Weir fit (sample-size corrected), $\\hat\\rho$ = {gw['rho']:.0f}; "
                   f"half-decay {gw['half']/1e3:.0f} kb; N$_e$ ≈ {gw['Ne']:.0f}"))
    ax.axvline(gw["half"] / 1e3, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xlabel("physical distance (kb)")
    ax.set_ylabel(r"r$^2$")
    ax.set_xlim(0, MAX_DIST_BP / 1e3)
    ax.set_ylim(0, max(0.3, float(np.nanmax(gw["means"])) * 1.15))
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True)
    fig.savefig(FIG / "fig_ld_decay_genomewide_ayb.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "fig_ld_decay_genomewide_ayb.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---- per-chromosome figure (3x4 grid) ----
    fig, axes = plt.subplots(3, 4, figsize=(13.5, 8.5), constrained_layout=True)
    axf = axes.flatten()
    for k, c in enumerate(CHRS):
        ax = axf[k]
        dd, rr, f = per_chr_fits[c]
        ax.scatter(dd / 1e3, rr, s=3, alpha=0.06, color=WONG["skyblue"],
                   edgecolor="none", rasterized=True)
        ax.errorbar(f["mids"] / 1e3, f["means"], fmt="o", color=WONG["blue"],
                    markersize=2.5, elinewidth=0.6)
        dpl = np.linspace(1, MAX_DIST_BP, 800)
        ax.plot(dpl / 1e3, hw_corr(dpl, f["rho"], n), color=WONG["vermillion"], linewidth=1.3)
        ax.axvline(f["half"] / 1e3, color="0.5", linestyle="--", linewidth=0.7)
        ax.text(0.0, 1.02, f"{c}: N$_e$≈{f['Ne']:.0f}, {f['half']/1e3:.0f} kb",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=9)
        ax.set_xlim(0, MAX_DIST_BP / 1e3)
        ax.set_ylim(0, 0.4)
        if k % 4 == 0:
            ax.set_ylabel(r"r$^2$")
        if k >= 7:
            ax.set_xlabel("distance (kb)")
        ax.grid(True)
    axf[-1].axis("off")
    fig.savefig(FIG / "fig28_ld_decay_per_chrom_ayb.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "fig28_ld_decay_per_chrom_ayb.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[fig] wrote genome-wide + per-chromosome figures")


if __name__ == "__main__":
    main()
