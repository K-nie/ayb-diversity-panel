#!/usr/bin/env python3
"""
Permutation null for the Ss05 outlier cluster-attribution count, recomputed
on the DE-DUPLICATED unique-SNP definition (Reviewer R1-M7 / R2-C).

Background
----------
The submitted manuscript reported "13 of 20" Ss05 outliers as minor-cluster
(Cluster-2) driven and tested that count against a permutation null
("null mean 0.0 +/- 0.03, empirical p <= 0.0001"). The "13" was the
method x SNP ROW count (13 of 22 rows in script 55's attribution table),
paired with an inconsistent "20" denominator; the same SNP flagged by both
PCAdapt and DAPC was counted twice. Script 55 recommends the de-duplicated
UNIQUE-SNP count: 7 of 16 unique Ss05 outliers are Cluster-2-driven.

This script re-runs the permutation null on that corrected unique-SNP count,
using the null described in the submitted Methods: cluster sizes are held at
11 / 84, and each SNP's genotypes are re-drawn from the binomial expectation
under the observed panel-wide allele frequency (destroying any real
cluster-structure association while preserving per-SNP allele frequency and
the 11 / 84 split). The per-permutation statistic is the number of the 16
unique Ss05 outlier SNPs that would be classified Cluster-2-driven
(MAF_C2 >= 0.40 AND MAF_C1 <= 0.10) by chance.

Addresses: R1-M7 (drift attribution) and the R2-C numerical-consistency fix.

Inputs
------
- AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv  raw genotypes
- results/01_qc_pca_power/tables/pca_coords.csv                       cluster labels
- results/55_private_alleles_per_cluster/tables/Ss05_outlier_attribution.csv
      the 22-row / 16-unique Ss05 outlier set + observed attribution

Outputs (results/94_ss05_attribution_permutation_null/)
------------------------------------------------------
tables/permutation_null_summary.csv   observed count, null mean/sd, empirical p
figures/fig_ss05_attribution_null.png/.pdf   observed vs null histogram

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
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
ATTR_CSV = ROOT / "results" / "55_private_alleles_per_cluster" / "tables" / "Ss05_outlier_attribution.csv"

OUT = ROOT / "results" / "94_ss05_attribution_permutation_null"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MAF_MIN = 0.05
ATTR_C_HIGH = 0.40
ATTR_C_LOW = 0.10
N_PERM = 10_000
SEED = 20260920


def dosage_row(row_vals: np.ndarray, alleles: str) -> np.ndarray:
    ref, alt = alleles.split("/")
    het1, het2 = ref + alt, alt + ref
    out = np.full(len(row_vals), np.nan)
    out[row_vals == ref + ref] = 0.0
    out[(row_vals == het1) | (row_vals == het2)] = 1.0
    out[row_vals == alt + alt] = 2.0
    return out


def build_working_matrix() -> pd.DataFrame:
    print("[load] reading HapMap")
    hm = pd.read_csv(HAPMAP, low_memory=False)
    meta_cols = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                 "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in meta_cols]
    calls = hm[sample_cols].astype(str).values
    alleles = hm["alleles"].values
    dosage = np.vstack([dosage_row(calls[i], alleles[i]) for i in range(len(hm))])
    dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
    samp_call = dosage.notna().sum(axis=0) / dosage.shape[0]
    dosage = dosage[samp_call[samp_call >= SAMP_CR].index.tolist()]
    mark_call = dosage.notna().sum(axis=1) / dosage.shape[1]
    mean_dose = dosage.mean(axis=1, skipna=True)
    maf_all = np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)
    dosage = dosage.loc[dosage.index[(mark_call >= MARK_CR) & (maf_all >= MAF_MIN)]]
    print(f"[qc] working matrix: markers={dosage.shape[0]}, samples={dosage.shape[1]}")
    return dosage


def main() -> None:
    dosage = build_working_matrix()
    samples = dosage.columns.tolist()

    pca = pd.read_csv(PCA_CSV)
    labels = dict(zip(pca["sample"], pca["cluster"]))
    # pca_coords convention: cluster==1 large (n=84, "Cluster 1"), cluster==0 small (n=11, "Cluster 2")
    n_c1 = sum(1 for s in samples if labels.get(s) == 1)
    n_c2 = sum(1 for s in samples if labels.get(s) == 0)
    print(f"[clust] Cluster 1 (large) n={n_c1}; Cluster 2 (small) n={n_c2}")
    assert (n_c1, n_c2) == (84, 11), (n_c1, n_c2)

    # Unique Ss05 outlier SNPs and the observed Cluster-2-driven count.
    attr = pd.read_csv(ATTR_CSV)
    uniq = attr.drop_duplicates("rs")
    rs_list = uniq["rs"].tolist()
    n_snps = len(rs_list)
    observed_c2 = int((uniq["attribution"] == "Cluster-2-driven").sum())
    print(f"[obs] unique Ss05 outlier SNPs={n_snps}; observed Cluster-2-driven={observed_c2}")

    # Panel-wide alt-allele frequency for each of the outlier SNPs.
    present = [rs for rs in rs_list if rs in dosage.index]
    missing = [rs for rs in rs_list if rs not in dosage.index]
    if missing:
        print(f"[warn] {len(missing)} outlier rs not in QC matrix (dropped): {missing}")
    X = dosage.loc[present].values                      # snps x samples
    p_alt = np.nansum(X, axis=1) / (2.0 * np.sum(~np.isnan(X), axis=1))
    minor_is_alt = p_alt <= 0.5
    k = len(present)
    n_tot = n_c1 + n_c2

    print(f"[null] {N_PERM} permutations, binomial redraw under panel MAF, 11/84 split")
    rng = np.random.default_rng(SEED)
    null_counts = np.empty(N_PERM, dtype=int)
    for j in range(N_PERM):
        # Draw genotypes for all n_tot individuals per SNP ~ Binomial(2, p_alt).
        g = rng.binomial(2, p_alt[:, None], size=(k, n_tot))   # snps x individuals
        # First n_c1 columns -> Cluster 1, remaining -> Cluster 2 (iid draws, so a fixed split is valid).
        pa_c1 = g[:, :n_c1].mean(axis=1) / 2.0
        pa_c2 = g[:, n_c1:].mean(axis=1) / 2.0
        maf_c1 = np.where(minor_is_alt, pa_c1, 1 - pa_c1)
        maf_c2 = np.where(minor_is_alt, pa_c2, 1 - pa_c2)
        c2_driven = (maf_c2 >= ATTR_C_HIGH) & (maf_c1 <= ATTR_C_LOW)
        null_counts[j] = int(c2_driven.sum())

    null_mean = float(null_counts.mean())
    null_sd = float(null_counts.std(ddof=1))
    # One-sided empirical p (add-one for the observed): P(null >= observed).
    p_emp = (int((null_counts >= observed_c2).sum()) + 1) / (N_PERM + 1)
    print(f"[null] mean={null_mean:.3f} sd={null_sd:.3f} max={null_counts.max()} "
          f"-> observed {observed_c2}/{k}; empirical p={p_emp:.2e}")

    summary = pd.DataFrame([{
        "n_unique_ss05_outliers": k,
        "observed_cluster2_driven": observed_c2,
        "null_mean": round(null_mean, 4),
        "null_sd": round(null_sd, 4),
        "null_max": int(null_counts.max()),
        "n_permutations": N_PERM,
        "empirical_p": p_emp,
        "seed": SEED,
    }])
    summary.to_csv(TAB / "permutation_null_summary.csv", index=False)
    print("\n[summary]\n", summary.to_string(index=False))

    # Figure: null distribution vs observed.
    fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
    bins = np.arange(-0.5, max(null_counts.max(), observed_c2) + 1.5, 1)
    ax.hist(null_counts, bins=bins, color=WONG["skyblue"], edgecolor="black",
            linewidth=0.5, label=f"Permutation null ({N_PERM:,} perms)")
    ax.axvline(observed_c2, color=WONG["vermillion"], linewidth=2.0,
               label=f"Observed = {observed_c2} / {k}")
    ax.set_xlabel("Cluster-2-driven Ss05 outlier SNPs")
    ax.set_ylabel("Permutations")
    ax.set_title("Ss05 minor-cluster attribution vs binomial permutation null")
    ax.legend(frameon=True, framealpha=0.95)
    fig.savefig(FIG / "fig_ss05_attribution_null.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "fig_ss05_attribution_null.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[done] wrote {TAB/'permutation_null_summary.csv'}")


if __name__ == "__main__":
    main()
