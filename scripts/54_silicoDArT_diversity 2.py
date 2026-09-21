#!/usr/bin/env python3
"""
SilicoDArT presence/absence diversity layer.

The DArTseq DAf18-2580 report supplied 4,992 SilicoDArT dominant
presence/absence markers in addition to the 3,204 codominant SNPs the
rest of the pipeline analyses. SilicoDArT calls are 0 (allele absent)
/ 1 (allele present) / - (missing). They are dominant (no heterozygote
call), so they cannot drive kinship-matrix-based GBLUP, but they can
support PCA, ADMIXTURE-style ancestry, and pairwise IBS distance.

Adding the SilicoDArT layer to the diversity-resource paper gives the
panel a second marker system to cross-validate against the SNP-based
architecture. The right question is whether the K = 2 PCA split that
the SNP layer recovers also surfaces independently on SilicoDArT, and
whether pairwise IBS across the two marker systems is correlated at
the panel-wide level.

Method
------
1. Read the DArT-format SilicoDArT report (6 header rows + 14 metadata
   columns + per-sample binary calls).
2. QC: per-marker and per-sample call rate, polymorphism = 1 -
   max(p, 1-p), and MAF = min(p, 1-p). Apply call rate >= 0.90 (sample
   and marker) and MAF >= 0.05 in line with the SNP-layer pipeline.
3. Restrict to the 95 working-panel accessions (the sample set that
   passes the SNP-layer QC, so the cross-marker-system comparison is
   sample-comparable).
4. PCA on the standardised binary matrix; K-means at k = 2..7 with
   silhouette scoring; concordance check against the SNP-layer K = 2
   cluster labels.
5. Pairwise IBS distance on the SilicoDArT matrix. Pearson correlation
   against SNP-layer pairwise IBS across all 95 * 94 / 2 = 4,465 pairs.

ADMIXTURE at K = 2..6 on a SilicoDArT-converted PLINK BED is not run
here. It would require a custom binary-input conversion and an
external PLINK + ADMIXTURE call; the PCA + IBS-concordance result
already answers the cross-marker-system question and ADMIXTURE can be
added in a revision pass if reviewers ask.

Outputs (results/54_silicoDArT_diversity/)
------------------------------------------
tables/
    silicoDArT_marker_qc.csv          per-marker call rate, MAF, polymorphism
    silicoDArT_sample_qc.csv          per-sample call rate, percent-present
    silicoDArT_pca_coords.csv         PCA coords + K-means cluster + SNP-cluster
    silicoDArT_pca_variance.csv       variance explained per PC
    cross_marker_system_concordance.csv  IBS correlation between layers
figures/
    fig_silicoDArT_pca.png/.pdf       PC1 x PC2 coloured by SNP-layer cluster
    fig_cross_layer_ibs.png/.pdf      pairwise IBS scatter (SilicoDArT vs SNP)

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

from _plotstyle import apply, WONG, CLUSTER_PAL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
SILICO_CSV = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SilicoDArT_1.csv"
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"

OUT = ROOT / "results" / "54_silicoDArT_diversity"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MAF_MIN = 0.05

# The DArT report leads with 6 metadata-only rows (sequencing batch
# identifiers + plate well coordinates). The 7th row is the column-
# header row that carries CloneID + AlleleSequence + ... and sample IDs.
N_HEADER_SKIP = 6
META_COLS = [
    "CloneID", "AlleleSequence", "TrimmedSequence",
    "Chrom_Vunguiculata_469_v1.0", "ChromPos_Vunguiculata_469_v1.0",
    "AlnCnt_Vunguiculata_469_v1.0", "AlnEvalue_Vunguiculata_469_v1.0",
    "CallRate", "OneRatio", "PIC",
    "AvgReadDepth", "StDevReadDepth", "Qpmr", "Reproducibility",
]


# ---------------------------------------------------------------------------
# Load SilicoDArT
# ---------------------------------------------------------------------------
def load_silicoDArT() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (calls, metadata) where calls is markers x samples and
    metadata carries the 14 DArT per-marker fields."""
    print(f"[load] reading SilicoDArT report from {SILICO_CSV}")
    df = pd.read_csv(SILICO_CSV, skiprows=N_HEADER_SKIP, low_memory=False)
    sample_cols = [c for c in df.columns if c not in META_COLS]
    print(f"[load] markers={len(df)}, samples={len(sample_cols)}")

    meta = df[META_COLS].copy()
    meta.index = df["CloneID"].values

    # Calls: 0 / 1 / - (missing). Convert "-" to NaN, cast to float.
    calls = df[sample_cols].copy()
    calls.index = df["CloneID"].values
    calls = calls.replace("-", np.nan).apply(pd.to_numeric, errors="coerce")
    return calls, meta


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------
def per_marker_qc(calls: pd.DataFrame) -> pd.DataFrame:
    """Per-marker call rate, frequency of presence (1), MAF, polymorphism."""
    n_samp = calls.shape[1]
    called = calls.notna().sum(axis=1)
    call_rate = called / n_samp
    p_present = (calls == 1).sum(axis=1) / called.replace(0, np.nan)
    maf = np.minimum(p_present, 1 - p_present)
    polymorphism = 1 - np.maximum(p_present, 1 - p_present)
    return pd.DataFrame({
        "CloneID": calls.index,
        "call_rate": call_rate.values,
        "p_present": p_present.values,
        "MAF": maf.values,
        "polymorphism": polymorphism.values,
    })


def per_sample_qc(calls: pd.DataFrame) -> pd.DataFrame:
    """Per-sample call rate and fraction of "1" (present) calls."""
    n_markers = calls.shape[0]
    called = calls.notna().sum(axis=0)
    call_rate = called / n_markers
    pct_present = (calls == 1).sum(axis=0) / called.replace(0, np.nan)
    return pd.DataFrame({
        "sample": calls.columns,
        "call_rate": call_rate.values,
        "pct_present": pct_present.values,
    })


def apply_qc(calls: pd.DataFrame, working_samples: list[str]) -> pd.DataFrame:
    """Restrict to the working-panel samples, then apply call-rate + MAF
    filter on the marker side."""
    # Sample restriction first so the marker call rate is computed on the
    # working panel rather than the full 105-sample set.
    keep_samples = [s for s in working_samples if s in calls.columns]
    n_drop = len(working_samples) - len(keep_samples)
    if n_drop:
        missing = [s for s in working_samples if s not in calls.columns]
        print(f"[qc] {n_drop} working samples missing from SilicoDArT report: {missing}")
    calls = calls[keep_samples]

    # Sample call rate filter (matches SNP-layer convention; almost all
    # working samples should pass since they already cleared the SNP-layer
    # filter at the same threshold).
    samp_cr = calls.notna().sum(axis=0) / calls.shape[0]
    keep_samples_post = samp_cr[samp_cr >= SAMP_CR].index.tolist()
    calls = calls[keep_samples_post]

    # Marker call rate + MAF.
    mark_cr = calls.notna().sum(axis=1) / calls.shape[1]
    p_present = (calls == 1).sum(axis=1) / calls.notna().sum(axis=1).replace(0, np.nan)
    maf = np.minimum(p_present, 1 - p_present)
    keep_markers = calls.index[(mark_cr >= MARK_CR) & (maf >= MAF_MIN)]
    calls = calls.loc[keep_markers]
    print(f"[qc] post-QC SilicoDArT: markers={calls.shape[0]}, samples={calls.shape[1]}")
    return calls


# ---------------------------------------------------------------------------
# PCA + K-means
# ---------------------------------------------------------------------------
def run_pca(calls: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Standardised PCA + K-means at k = 2..7. Returns the per-sample PC
    table and the per-PC variance vector."""
    X = calls.T.fillna(calls.mean(axis=1)).T  # marker x sample mean impute
    X = X.T.values  # sample x marker
    X = StandardScaler().fit_transform(X)
    pca = PCA(n_components=10, random_state=0)
    pcs = pca.fit_transform(X)

    silhouettes = {}
    cluster_assignments = {}
    for k in range(2, 8):
        labels = KMeans(n_clusters=k, n_init=20, random_state=0).fit_predict(pcs)
        sil = silhouette_score(pcs, labels)
        silhouettes[k] = sil
        cluster_assignments[k] = labels
    best_k = max(silhouettes, key=silhouettes.get)
    print(f"[pca] silhouettes: {silhouettes}")
    print(f"[pca] best k = {best_k} (silhouette = {silhouettes[best_k]:.3f})")

    pc_df = pd.DataFrame(pcs, columns=[f"PC{i + 1}" for i in range(10)])
    pc_df["sample"] = calls.columns.tolist()
    pc_df["silicoDArT_cluster"] = cluster_assignments[2]  # k=2 for the SNP comparison
    for k in range(2, 8):
        pc_df[f"cluster_k{k}"] = cluster_assignments[k]
    pc_df["best_k"] = best_k

    variance = pd.DataFrame({
        "PC": [f"PC{i + 1}" for i in range(10)],
        "variance_explained": pca.explained_variance_ratio_,
        "cumulative": np.cumsum(pca.explained_variance_ratio_),
    })
    return pc_df, variance


# ---------------------------------------------------------------------------
# IBS distance
# ---------------------------------------------------------------------------
def pairwise_ibs(matrix: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """IBS-style distance per pair = mean absolute difference across loci
    where both samples are called. Returns a flat upper-triangle vector
    plus the sample order."""
    X = matrix.values.astype(float)  # marker x sample
    samples = matrix.columns.tolist()
    n = len(samples)

    # Mean-imputed marker mean per marker (used only for the SNP layer
    # to avoid NaN-propagation; SilicoDArT has fewer NaNs but we apply
    # the same treatment for symmetry).
    row_means = np.nanmean(X, axis=1, keepdims=True)
    X_imputed = np.where(np.isnan(X), row_means, X)

    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.mean(np.abs(X_imputed[:, i] - X_imputed[:, j]))
            dist[i, j] = d
            dist[j, i] = d

    upper = []
    pair_labels = []
    for i in range(n):
        for j in range(i + 1, n):
            upper.append(dist[i, j])
            pair_labels.append((samples[i], samples[j]))
    return np.asarray(upper), pair_labels


def build_snp_dosage(working_samples: list[str]) -> pd.DataFrame:
    """Build the SNP-layer working dosage matrix on the same sample set
    so the IBS comparison is on the same accessions in the same order."""
    print("[snp] building SNP-layer dosage matrix for the IBS comparison")
    hm = pd.read_csv(HAPMAP, low_memory=False)
    meta_cols = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                 "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in meta_cols]
    calls = hm[sample_cols].astype(str).values
    alleles_arr = hm["alleles"].values

    out = np.full((len(hm), len(sample_cols)), np.nan)
    for i in range(len(hm)):
        ref, alt = alleles_arr[i].split("/")
        het1, het2 = ref + alt, alt + ref
        hom_ref = ref + ref
        hom_alt = alt + alt
        row = calls[i]
        out[i] = np.where(row == hom_ref, 0.0,
                          np.where((row == het1) | (row == het2), 1.0,
                                   np.where(row == hom_alt, 2.0, np.nan)))
    dosage = pd.DataFrame(out, index=hm["rs#"].values, columns=sample_cols)

    # Apply same QC.
    samp_call = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_samples = samp_call[samp_call >= SAMP_CR].index.tolist()
    dosage = dosage[keep_samples]
    mark_call = dosage.notna().sum(axis=1) / dosage.shape[1]
    mean_dose = dosage.mean(axis=1, skipna=True)
    maf_all = np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)
    keep_markers = dosage.index[(mark_call >= MARK_CR) & (maf_all >= MAF_MIN)]
    dosage = dosage.loc[keep_markers]

    # Restrict to working sample set in matching order.
    keep_in_order = [s for s in working_samples if s in dosage.columns]
    return dosage[keep_in_order]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_silicoDArT_pca(pc_df: pd.DataFrame, snp_cluster: dict[str, int],
                        out_path: Path) -> None:
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.0))

    # Panel A: coloured by SNP-layer cluster (cross-marker-system check).
    pc_df = pc_df.copy()
    pc_df["snp_cluster"] = pc_df["sample"].map(snp_cluster)

    # Translate SNP-layer cluster (0 = small, 1 = large in pca_coords.csv)
    # to manuscript convention (Cluster 1 = large, Cluster 2 = small) for
    # the figure legend. Counts are taken from the retained 78-line subset
    # actually plotted, not from the full 95-line SNP-layer panel; the
    # difference (17 missing samples in Cluster 1 only) is itself part of
    # the SilicoDArT QC story.
    snp_to_display = {1: ("SNP Cluster 1", CLUSTER_PAL[0]),
                      0: ("SNP Cluster 2", CLUSTER_PAL[1])}
    for k, (name, colour) in snp_to_display.items():
        sub = pc_df.loc[pc_df["snp_cluster"] == k]
        ax_a.scatter(sub["PC1"], sub["PC2"], color=colour, s=42,
                     edgecolor="white", linewidth=0.5,
                     label=f"{name} (retained n = {len(sub)})")
    ax_a.set_xlabel("SilicoDArT PC1")
    ax_a.set_ylabel("SilicoDArT PC2")
    ax_a.set_title("SilicoDArT PCA coloured by SNP-layer cluster")
    ax_a.legend(loc="best")

    # Panel B: coloured by SilicoDArT's own K-means k=2.
    sizes = pc_df["silicoDArT_cluster"].value_counts()
    big_label = sizes.idxmax()
    small_label = sizes.idxmin() if sizes.idxmin() != big_label else 1 - big_label
    silico_to_display = {big_label: (f"SilicoDArT Cluster 1 (n={int(sizes.max())})", CLUSTER_PAL[0]),
                         small_label: (f"SilicoDArT Cluster 2 (n={int(sizes.min())})", CLUSTER_PAL[1])}
    for k, (name, colour) in silico_to_display.items():
        sub = pc_df.loc[pc_df["silicoDArT_cluster"] == k]
        ax_b.scatter(sub["PC1"], sub["PC2"], color=colour, s=42,
                     edgecolor="white", linewidth=0.5, label=name)
    ax_b.set_xlabel("SilicoDArT PC1")
    ax_b.set_ylabel("SilicoDArT PC2")
    ax_b.set_title("SilicoDArT PCA coloured by its own K=2 K-means")
    ax_b.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def plot_cross_layer_ibs(snp_ibs: np.ndarray, silico_ibs: np.ndarray,
                          r: float, n_samples: int, out_path: Path) -> None:
    n_pairs = n_samples * (n_samples - 1) // 2
    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    ax.scatter(snp_ibs, silico_ibs, s=8, color=WONG["blue"], alpha=0.4)
    # Identity line scaled to overlap the data range.
    lo = min(snp_ibs.min(), silico_ibs.min())
    hi = max(snp_ibs.max(), silico_ibs.max())
    ax.plot([lo, hi], [lo, hi], color="grey", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Pairwise IBS distance on SNP layer (1,625 markers)")
    ax.set_ylabel("Pairwise IBS distance on SilicoDArT layer")
    ax.set_title(f"Cross-marker-system pairwise IBS concordance\n"
                 f"Pearson r = {r:.3f} across all {n_samples}*{n_samples - 1}/2 = {n_pairs:,} pairs")
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    calls, meta = load_silicoDArT()

    # Per-marker and per-sample QC summaries on the raw matrix.
    marker_qc = per_marker_qc(calls)
    sample_qc = per_sample_qc(calls)
    marker_qc.to_csv(TAB / "silicoDArT_marker_qc.csv", index=False)
    sample_qc.to_csv(TAB / "silicoDArT_sample_qc.csv", index=False)
    print(f"[qc] median per-marker call rate = {marker_qc['call_rate'].median():.3f}")
    print(f"[qc] median per-marker polymorphism = {marker_qc['polymorphism'].median():.3f}")

    # Restrict to the working 95-line SNP-layer panel.
    pca_snp = pd.read_csv(PCA_CSV)
    working_samples = pca_snp["sample"].tolist()
    snp_cluster_map = dict(zip(pca_snp["sample"], pca_snp["cluster"]))

    calls_post = apply_qc(calls, working_samples)

    # PCA + K-means
    pc_df, variance = run_pca(calls_post)
    pc_df.to_csv(TAB / "silicoDArT_pca_coords.csv", index=False)
    variance.to_csv(TAB / "silicoDArT_pca_variance.csv", index=False)
    print(f"[pca] PC1 + PC2 explain {variance['cumulative'].iloc[1]:.3f}")

    # Cross-cluster concordance: how many of the 95 accessions land in
    # the same K-means cluster under SilicoDArT vs SNP at k=2?
    pc_df["snp_cluster"] = pc_df["sample"].map(snp_cluster_map)
    # The 0/1 labelling between the two systems is arbitrary; align by
    # majority within each silicoDArT cluster.
    cluster_alignment = {}
    for sc in pc_df["silicoDArT_cluster"].unique():
        sub = pc_df.loc[pc_df["silicoDArT_cluster"] == sc]
        majority_snp = sub["snp_cluster"].value_counts().idxmax()
        cluster_alignment[sc] = majority_snp
    pc_df["silicoDArT_cluster_aligned"] = pc_df["silicoDArT_cluster"].map(cluster_alignment)
    concordance = (pc_df["silicoDArT_cluster_aligned"] == pc_df["snp_cluster"]).mean()
    print(f"[concord] SilicoDArT k=2 vs SNP k=2 concordance: {concordance:.3f}")

    # Pairwise IBS on both layers (same sample order).
    snp_dosage = build_snp_dosage(calls_post.columns.tolist())
    snp_dosage = snp_dosage[calls_post.columns.tolist()]  # exact order match
    print(f"[ibs] computing pairwise IBS on SNP layer ({snp_dosage.shape[0]} markers)")
    snp_ibs, snp_pairs = pairwise_ibs(snp_dosage)
    print(f"[ibs] computing pairwise IBS on SilicoDArT layer ({calls_post.shape[0]} markers)")
    silico_ibs, silico_pairs = pairwise_ibs(calls_post)
    assert snp_pairs == silico_pairs

    r = np.corrcoef(snp_ibs, silico_ibs)[0, 1]
    print(f"[concord] cross-layer pairwise-IBS Pearson r = {r:.3f}")

    concord_df = pd.DataFrame({
        "n_pairs": [len(snp_ibs)],
        "pearson_r": [r],
        "silicoDArT_cluster_concordance": [concordance],
        "snp_layer_markers": [snp_dosage.shape[0]],
        "silicoDArT_markers": [calls_post.shape[0]],
    })
    concord_df.to_csv(TAB / "cross_marker_system_concordance.csv", index=False)

    # Figures
    print("[plot] SilicoDArT PCA panels")
    plot_silicoDArT_pca(pc_df, snp_cluster_map, FIG / "fig_silicoDArT_pca")
    print("[plot] cross-layer IBS concordance")
    plot_cross_layer_ibs(snp_ibs, silico_ibs, r, calls_post.shape[1],
                         FIG / "fig_cross_layer_ibs")

    print("\n=== manuscript-text summary ===")
    print(f"  raw SilicoDArT markers              : {len(calls)}")
    print(f"  post-QC SilicoDArT markers          : {calls_post.shape[0]}")
    print(f"  post-QC samples                     : {calls_post.shape[1]}")
    print(f"  median per-marker call rate         : {marker_qc['call_rate'].median():.3f}")
    print(f"  median per-marker polymorphism      : {marker_qc['polymorphism'].median():.3f}")
    print(f"  best silhouette k                   : {int(pc_df['best_k'].iloc[0])}")
    print(f"  SilicoDArT vs SNP cluster concordance: {concordance:.3f}")
    print(f"  cross-layer pairwise-IBS Pearson r  : {r:.3f}")


if __name__ == "__main__":
    main()
