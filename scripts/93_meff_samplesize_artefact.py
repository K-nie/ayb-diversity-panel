#!/usr/bin/env python3
"""
93 - m_eff is bounded by sample size, not marker independence.

Reviewer 2 point F: the submitted effective number of independent tests
m_eff = 93 (Li & Ji 2005) was used to set a multiple-testing threshold
(alpha = 0.05 / m_eff = 5.4e-4). But m_eff was derived from the eigenvalue
spectrum of a marker correlation matrix estimated from n = 95 individuals.
Such a matrix has rank <= n - 1, so *at most n - 1 = 94 non-zero eigenvalues*
exist regardless of how many markers there are. The observed m_eff = 93 is
therefore pinned to the sample-size ceiling and does not measure the number
of independent markers - it measures (nearly) the number of individuals.

This script demonstrates the artefact directly:
  1. Reproduces the full-panel Li & Ji m_eff (should be ~93).
  2. Recomputes m_eff on random subsamples of individuals at n' in
     {20, 30, 40, 60, 80, 95}. If m_eff tracked marker LD structure it would
     be roughly constant; if it is a sample-size artefact it will track
     n' - 1. We show the latter.
  3. Recomputes m_eff after permuting genotypes within each marker (destroying
     all LD). A genuine marker-independence measure would then rise toward the
     marker count m; a rank-bounded one stays at ~n - 1. We show the latter.

Conclusion for the revision: report the GWAS/outlier significance threshold at
the true marker-count Bonferroni (alpha = 0.05 / m) and/or an FDR / permutation
threshold, and either drop the m_eff-scale threshold or report it only with an
explicit statement that at n = 95 it is bounded above by n - 1 and so cannot
be interpreted as a count of independent markers.

Inputs
  - AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv

Outputs (results/93_meff_samplesize_artefact/)
  tables/meff_vs_samplesize.csv
  tables/meff_artefact_summary.csv
  figures/fig_meff_vs_samplesize.png/.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _plotstyle import apply, WONG  # noqa: E402
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT = ROOT / "results" / "93_meff_samplesize_artefact"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
SEED = 20260920
N_REPS = 20                      # subsample replicates per n'


def li_ji_meff_from_X(X: np.ndarray) -> float:
    """Li & Ji (2005) m_eff from a samples x markers matrix X.
    Uses the n x n Gram of standardised columns, whose non-zero eigenvalues
    equal those of the m x m correlation matrix; the remaining (m - n) are 0
    and contribute 0 to the Li-Ji sum. This is exactly the submitted method."""
    n, m = X.shape
    Xs = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
    Xs = np.nan_to_num(Xs, nan=0.0)
    # drop constant columns (std 0 -> all zeros after nan_to_num)
    good = Xs.std(axis=0) > 0
    Xs = Xs[:, good]
    m_used = Xs.shape[1]
    G = Xs @ Xs.T / (n - 1)
    eig = np.linalg.eigvalsh(G)
    eig = np.clip(eig, 0.0, None)
    # full m x m spectrum: these n eigenvalues plus (m_used - n) zeros
    contrib = (eig >= 1).astype(float) + np.where(eig < 1, eig - np.floor(eig), 0.0)
    # zeros contribute 0 (indicator 0, frac 0)
    return float(contrib.sum()), n, m_used


def main() -> None:
    print("[load] HapMap + QC filter")
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

    cr_m = dosage.notna().sum(axis=1) / dosage.shape[1]
    maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                     1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
    cr_s = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_m = ((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values
    keep_s = (cr_s >= SAMP_CR).values
    sub = dosage.loc[keep_m, keep_s]
    M = sub.fillna(sub.mean(axis=1)).T.values   # samples x markers
    n_full, m = M.shape
    print(f"[qc] post-QC panel: n = {n_full} individuals, m = {m} markers")

    # 1. full-panel m_eff
    meff_full, _, _ = li_ji_meff_from_X(M)
    print(f"[full] m_eff = {meff_full:.1f}  (ceiling n-1 = {n_full-1})")

    # 2. m_eff vs subsample size
    rng = np.random.default_rng(SEED)
    n_grid = [20, 30, 40, 60, 80, n_full]
    rows = []
    for npr in n_grid:
        vals = []
        reps = 1 if npr == n_full else N_REPS
        for _ in range(reps):
            idx = rng.choice(n_full, size=npr, replace=False)
            meff, _, _ = li_ji_meff_from_X(M[idx])
            vals.append(meff)
        rows.append({"n_prime": npr, "n_minus_1": npr - 1,
                     "meff_mean": float(np.mean(vals)),
                     "meff_sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                     "meff_over_ceiling": float(np.mean(vals)) / (npr - 1)})
    vs = pd.DataFrame(rows)
    vs.to_csv(TAB / "meff_vs_samplesize.csv", index=False)
    print("[subsample]")
    print(vs.to_string(index=False))

    # 3. LD-destroyed control: permute each marker's genotypes across individuals
    Mp = M.copy()
    for j in range(m):
        Mp[:, j] = rng.permutation(Mp[:, j])
    meff_perm, _, _ = li_ji_meff_from_X(Mp)
    print(f"[LD-destroyed] m_eff = {meff_perm:.1f}  "
          f"(would approach m = {m} if m_eff measured independence; "
          f"stays ~n-1 = {n_full-1} instead)")

    bonf_meff = 0.05 / meff_full
    bonf_m = 0.05 / m
    summ = pd.DataFrame([
        {"metric": "post-QC individuals n", "value": n_full},
        {"metric": "post-QC markers m", "value": m},
        {"metric": "sample-size ceiling n-1", "value": n_full - 1},
        {"metric": "m_eff full panel (Li & Ji)", "value": round(meff_full, 1)},
        {"metric": "m_eff / (n-1) ceiling ratio", "value": round(meff_full / (n_full - 1), 3)},
        {"metric": "m_eff after destroying all LD (permuted)", "value": round(meff_perm, 1)},
        {"metric": "Bonferroni at m_eff (submitted)", "value": bonf_meff},
        {"metric": "Bonferroni at true marker count m", "value": bonf_m},
        {"metric": "-log10 Bonferroni at m_eff", "value": round(-np.log10(bonf_meff), 2)},
        {"metric": "-log10 Bonferroni at m", "value": round(-np.log10(bonf_m), 2)},
    ])
    summ.to_csv(TAB / "meff_artefact_summary.csv", index=False)
    print("\n[summary]")
    print(summ.to_string(index=False))

    # figure
    fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
    ax.plot(vs["n_prime"], vs["n_minus_1"], color="0.4", linestyle="--",
            linewidth=1.2, label="sample-size ceiling (n − 1)")
    ax.errorbar(vs["n_prime"], vs["meff_mean"], yerr=vs["meff_sd"],
                fmt="o-", color=WONG["blue"], markersize=5, capsize=3,
                elinewidth=0.9, label="m$_{eff}$ (Li & Ji)")
    ax.axhline(meff_perm, color=WONG["vermillion"], linestyle=":", linewidth=1.2,
               label=f"m$_{{eff}}$ with all LD destroyed = {meff_perm:.0f}")
    ax.annotate(f"m = {m} markers\n(true independence ceiling)",
                xy=(0.02, 0.92), xycoords="axes fraction", fontsize=8,
                va="top", color=WONG["vermillion"])
    ax.set_xlabel("number of individuals subsampled (n′)")
    ax.set_ylabel(r"effective number of independent tests m$_{eff}$")
    ax.set_title("m$_{eff}$ tracks sample size, not marker count\n"
                 f"(full panel n = {n_full}: m$_{{eff}}$ = {meff_full:.0f} ≈ n−1 = {n_full-1})",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True)
    fig.savefig(FIG / "fig_meff_vs_samplesize.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "fig_meff_vs_samplesize.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[fig] fig_meff_vs_samplesize")


if __name__ == "__main__":
    main()
