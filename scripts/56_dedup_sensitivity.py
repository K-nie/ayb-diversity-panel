#!/usr/bin/env python3
"""
De-duplication sensitivity for the 95-line panel.

Four sets of near-clonal duplicate accessions were surfaced by script 13
in the IITA TSs panel (G_ij >= 0.85). Removing one representative per
duplicate component reduces the panel from 95 to some n' < 95. The
operational question for the breeder and for downstream re-users of the
panel is: by how much does the de-duplication shift the panel-level
statistics we report?

This script does the comparison directly. It re-uses the same dosage
construction + QC cascade as script 01, then identifies duplicate
components from the GRM (script 13) via union-find on G_ij >= 0.85
pairs, picks one representative per component, and recomputes four
panel-level statistics on both the 95-line and the de-duplicated subset:

    PCA + K-means at k=2 silhouette       : does the cluster split hold?
    Per-marker MAF distribution            : mean / median shift
    Li & Ji eigenvalue m_eff               : effective independent tests
    Weir-Cockerham F_ST between clusters   : structural signal magnitude

LD-based Ne (Hill-Weir) and ADMIXTURE Q-matrix concordance are NOT
recomputed here -- they require pairwise r-squared and an external PLINK
+ ADMIXTURE call. The four statistics above are the cheapest informative
subset; the LD-Ne and ADMIXTURE checks are flagged as supplementary in
the manuscript text.

Outputs (results/56_dedup_sensitivity/)
---------------------------------------
tables/
    duplicate_components.csv       per-component members, G_ij values, representative
    dedup_delta_summary.csv         four panel statistics x (95-line, dedup, delta)
figures/
    fig_dedup_delta_forest.png/.pdf forest of normalised deltas
    fig_pca_overlay.png/.pdf        PC1 x PC2 with dropped accessions marked

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from collections import defaultdict
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
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
GRM_CSV = ROOT / "results" / "13_grm_crosspairs" / "tables" / "grm_vanraden.csv"
SAMPLE_QC_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "sample_qc.csv"
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"

OUT = ROOT / "results" / "56_dedup_sensitivity"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MAF_MIN = 0.05
DUP_G_THRESH = 0.85


# ---------------------------------------------------------------------------
# Build the working 95 x 1625 dosage matrix (mirrors script 01 logic)
# ---------------------------------------------------------------------------
def dosage_row(row_vals: np.ndarray, alleles: str) -> np.ndarray:
    ref, alt = alleles.split("/")
    het1, het2 = ref + alt, alt + ref
    hom_ref = ref + ref
    hom_alt = alt + alt
    out = np.full(len(row_vals), np.nan)
    out[row_vals == hom_ref] = 0.0
    out[(row_vals == het1) | (row_vals == het2)] = 1.0
    out[row_vals == hom_alt] = 2.0
    return out


def build_working_matrix() -> pd.DataFrame:
    """Return the 95 x 1625 post-QC dosage matrix with marker IDs as the
    index and sample IDs as columns."""
    print("[load] reading HapMap")
    hm = pd.read_csv(HAPMAP, low_memory=False)
    meta_cols = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                 "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in meta_cols]
    calls = hm[sample_cols].astype(str).values
    alleles = hm["alleles"].values

    print("[load] converting calls to dosage")
    dosage = np.vstack([dosage_row(calls[i], alleles[i]) for i in range(len(hm))])
    dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
    print(f"[load] raw: markers={dosage.shape[0]}, samples={dosage.shape[1]}")

    # Sample QC: call rate >= 0.90
    samp_call = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_samples = samp_call[samp_call >= SAMP_CR].index.tolist()
    dosage = dosage[keep_samples]

    # Marker QC: call rate >= 0.90 AND MAF >= 0.05
    mark_call = dosage.notna().sum(axis=1) / dosage.shape[1]
    mean_dose = dosage.mean(axis=1, skipna=True)
    maf = np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)
    keep_markers = dosage.index[(mark_call >= MARK_CR) & (maf >= MAF_MIN)]
    dosage = dosage.loc[keep_markers]

    print(f"[qc] post-QC: markers={dosage.shape[0]}, samples={dosage.shape[1]}")
    return dosage


# ---------------------------------------------------------------------------
# Duplicate-component detection via union-find on G_ij >= threshold pairs
# ---------------------------------------------------------------------------
class UnionFind:
    def __init__(self, elements):
        self.parent = {e: e for e in elements}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.parent[rx] = ry

    def components(self):
        out = defaultdict(list)
        for e in self.parent:
            out[self.find(e)].append(e)
        return {root: sorted(members) for root, members in out.items()
                if len(members) > 1}


def find_duplicate_components(grm: pd.DataFrame, threshold: float) -> list[dict]:
    """Return one record per duplicate component (pair or trio): members,
    pairwise G_ij values, and the picked representative."""
    samples = grm.index.tolist()
    uf = UnionFind(samples)
    pairs = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            g = float(grm.iloc[i, j])
            if g >= threshold:
                uf.union(samples[i], samples[j])
                pairs.append((samples[i], samples[j], g))

    components = uf.components()

    sample_qc = pd.read_csv(SAMPLE_QC_CSV).set_index("sample")
    pca = pd.read_csv(PCA_CSV).set_index("sample")

    out = []
    for members in components.values():
        # Pick representative by highest sample call rate; tie-break by
        # alphabetical TSs ID so the choice is reproducible.
        call_rates = {m: sample_qc.loc[m, "call_rate"] for m in members
                      if m in sample_qc.index}
        rep = max(members, key=lambda m: (call_rates.get(m, 0.0), m))

        # Pairwise G_ij values for the methods table.
        pair_vals = []
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                pair_vals.append((a, b, float(grm.loc[a, b])))

        cluster_labels = [int(pca.loc[m, "cluster"]) for m in members
                          if m in pca.index]

        out.append({
            "n_members": len(members),
            "members": ";".join(members),
            "representative": rep,
            "dropped": ";".join(m for m in members if m != rep),
            "G_ij_pairs": ";".join(f"{a}-{b}={g:.3f}" for a, b, g in pair_vals),
            "G_ij_min": min(g for _, _, g in pair_vals),
            "G_ij_max": max(g for _, _, g in pair_vals),
            "clusters": ",".join(map(str, sorted(set(cluster_labels)))),
        })
    return out


# ---------------------------------------------------------------------------
# Recompute the four headline statistics on a sample subset
# ---------------------------------------------------------------------------
def mean_impute(dosage: pd.DataFrame) -> pd.DataFrame:
    """Per-marker mean imputation across the supplied samples."""
    return dosage.T.fillna(dosage.mean(axis=1)).T


def pca_silhouette(dosage: pd.DataFrame, k: int = 2) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """PCA on standardised dosage; K-means at k clusters; return silhouette
    score, PC1, PC2, cluster labels (per-sample order = dosage.columns).

    K-means and the silhouette are computed on the leading five PCs, matching
    the canonical PCA pipeline (script 01 `01_qc_pca_power.py`) so the 95-line
    baseline reproduces the manuscript headline silhouette exactly rather than
    the inflated PC1-2-only value.
    """
    X = mean_impute(dosage).T.values  # samples x markers
    X = StandardScaler().fit_transform(X)
    pcs = PCA(n_components=min(10, X.shape[0] - 1, X.shape[1]),
              random_state=0).fit_transform(X)
    labels = KMeans(n_clusters=k, n_init=25, random_state=0).fit_predict(pcs[:, :5])
    sil = silhouette_score(pcs[:, :5], labels)
    return float(sil), pcs[:, 0], pcs[:, 1], labels


def per_marker_maf(dosage: pd.DataFrame) -> np.ndarray:
    mean_dose = dosage.mean(axis=1, skipna=True).values
    return np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)


def li_ji_meff(dosage: pd.DataFrame) -> int:
    """Li & Ji 2005 effective number of independent tests.

    Eigenvalues of the m x m correlation matrix are obtained from the
    smaller n x n Gram matrix (which shares non-zero eigenvalues). Per
    Li & Ji the m_eff is the sum of (1 if lambda_i >= 1 else 0) + the
    fractional part of each eigenvalue.
    """
    X = mean_impute(dosage).T.values  # samples x markers
    X = StandardScaler().fit_transform(X)
    # X X^T / (m - 1) shares non-zero eigenvalues with X^T X / (m - 1).
    m = X.shape[1]
    gram = X @ X.T / (m - 1)
    eigvals = np.linalg.eigvalsh(gram)
    # Drop near-zero numerical-noise eigenvalues.
    eigvals = eigvals[eigvals > 1e-10]
    meff = float(np.sum(np.where(eigvals >= 1, 1.0, eigvals - np.floor(eigvals))))
    return int(round(meff))


def weir_cockerham_fst(dosage: pd.DataFrame, labels: np.ndarray) -> float:
    """Global Weir & Cockerham (1984) F_ST (ratio of sums) between clusters.

    Uses the full three-component estimator including the observed-heterozygote
    term, computed on the raw (non-imputed) dosage with native missing-data
    handling, identical to the canonical F_ST script (`10_fst.py`). This makes
    the 95-line baseline reproduce the manuscript headline W&C F_ST rather than
    the inbred-line-simplified value the earlier version returned.

    dosage : marker x sample (may contain NaN)
    labels : per-sample cluster assignment in dosage.columns order
    """
    M_ = dosage.values  # marker x sample, NaN preserved
    grp_ = np.asarray(labels)
    pops = sorted(set(grp_.tolist()))
    r = len(pops)
    if r < 2:
        return float("nan")
    L = M_.shape[0]
    n_i = np.zeros((r, L)); p_i = np.zeros((r, L)); h_i = np.zeros((r, L))
    for ki, k in enumerate(pops):
        sub = M_[:, np.where(grp_ == k)[0]]
        nm = ~np.isnan(sub)
        n_i[ki] = nm.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            p_i[ki] = np.nansum(sub, axis=1) / (2.0 * n_i[ki])
            h_i[ki] = (sub == 1).sum(axis=1) / n_i[ki]
    n_tot = n_i.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        n_bar = n_tot / r
        p_bar = (n_i * p_i).sum(axis=0) / n_tot
        S2 = (1.0 / ((r - 1) * n_bar)) * (n_i * (p_i - p_bar) ** 2).sum(axis=0)
        h_bar = (n_i * h_i).sum(axis=0) / n_tot
        n_c = (n_tot - (n_i ** 2).sum(axis=0) / n_tot) / (r - 1)
        common = p_bar * (1.0 - p_bar) - ((r - 1) / r) * S2 - h_bar / 4.0
        a = (n_bar / n_c) * (S2 - (1.0 / (n_bar - 1)) * common)
        b = (n_bar / (n_bar - 1)) * (
            p_bar * (1.0 - p_bar) - ((r - 1) / r) * S2
            - ((2 * n_bar - 1) / (4 * n_bar)) * h_bar
        )
        c = h_bar / 2.0
    denom = a + b + c
    valid = np.isfinite(denom) & (denom > 0)
    if not valid.any():
        return float("nan")
    return float(np.nansum(a[valid]) / np.nansum(denom[valid]))


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def summarise_subset(dosage: pd.DataFrame, label: str) -> dict:
    print(f"\n[{label}] recomputing panel statistics on {dosage.shape[1]} accessions")

    sil, pc1, pc2, kmeans_labels = pca_silhouette(dosage, k=2)
    maf = per_marker_maf(dosage)
    meff = li_ji_meff(dosage)
    fst = weir_cockerham_fst(dosage, kmeans_labels)

    # Cluster sizes from the K-means split (smaller cluster always reported
    # as Cluster 2 by convention so the size compares apples-to-apples).
    sizes = pd.Series(kmeans_labels).value_counts().sort_values(ascending=False)
    n_c1 = int(sizes.iloc[0])
    n_c2 = int(sizes.iloc[1]) if len(sizes) > 1 else 0

    return {
        "label": label,
        "n_accessions": dosage.shape[1],
        "silhouette_k2": sil,
        "cluster1_size": n_c1,
        "cluster2_size": n_c2,
        "mean_MAF": float(np.mean(maf)),
        "median_MAF": float(np.median(maf)),
        "m_eff": meff,
        "F_ST_between_clusters": fst,
        "PC1": pc1,
        "PC2": pc2,
        "kmeans_labels": kmeans_labels,
    }


def plot_delta_forest(deltas: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    y = np.arange(len(deltas))
    ax.scatter(deltas["pct_change"], y, color=WONG["vermillion"], s=70, zorder=3)
    ax.axvline(0, color="grey", linewidth=0.8)
    for i, row in deltas.iterrows():
        ax.text(row["pct_change"] + 0.5, i,
                f"{row['value_95line']:.3g} -> {row['value_dedup']:.3g}",
                va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(deltas["stat"])
    ax.set_xlabel("Percent change after dropping duplicate representatives")
    ax.text(-0.14, 1.04, "c", transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="right")
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def plot_pca_overlay(stats_95: dict, dropped: list[str], samples_95: list[str],
                     out_path: Path) -> None:
    """PC1 x PC2 with the dropped duplicates circled.

    The K-means label that appears as 0 vs 1 depends on initialisation, so
    we re-map by size: the larger cluster is plotted as Cluster 1 and the
    smaller as Cluster 2 so the figure legend stays consistent with the
    manuscript convention.
    """
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    is_dropped = np.array([s in dropped for s in samples_95])
    labels = np.asarray(stats_95["kmeans_labels"])

    sizes = pd.Series(labels).value_counts()
    big_label = sizes.idxmax()
    small_label = sizes.idxmin() if sizes.idxmin() != big_label else 1 - big_label

    label_to_display = {big_label: ("Cluster 1", CLUSTER_PAL[0]),
                        small_label: ("Cluster 2", CLUSTER_PAL[1])}

    for k, (name, colour) in label_to_display.items():
        mask = (labels == k) & (~is_dropped)
        ax.scatter(stats_95["PC1"][mask], stats_95["PC2"][mask],
                   color=colour, s=42, edgecolor="white",
                   linewidth=0.5, label=f"{name} (retained, n = {mask.sum()})")

    ax.scatter(stats_95["PC1"][is_dropped], stats_95["PC2"][is_dropped],
               facecolor="none", edgecolor="black", s=160,
               linewidth=1.6, label=f"dropped duplicate (n = {is_dropped.sum()})")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.text(-0.10, 1.03, "d", transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="right")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    dosage_95 = build_working_matrix()
    samples_95 = dosage_95.columns.tolist()

    print("\n[grm] reading VanRaden GRM")
    grm = pd.read_csv(GRM_CSV, index_col=0)
    print(f"[grm] {grm.shape[0]} x {grm.shape[1]} matrix")

    # Restrict the GRM to samples present in the working dosage matrix so
    # any historical drift between the QC step and the script-13 GRM does
    # not bleed into the duplicate detection.
    common = [s for s in samples_95 if s in grm.index]
    grm = grm.loc[common, common]

    print(f"\n[dedup] finding components at G_ij >= {DUP_G_THRESH}")
    components = find_duplicate_components(grm, DUP_G_THRESH)
    comp_df = pd.DataFrame(components)
    comp_path = TAB / "duplicate_components.csv"
    comp_df.to_csv(comp_path, index=False)
    print(f"[dedup] {len(components)} components flagged")
    print(comp_df[["n_members", "members", "representative", "G_ij_min",
                   "G_ij_max", "clusters"]].to_string(index=False))

    dropped = []
    for c in components:
        dropped.extend(c["dropped"].split(";"))
    print(f"[dedup] dropping {len(dropped)} accessions: {dropped}")

    # 95-line baseline statistics
    stats_95 = summarise_subset(dosage_95, "95-line")

    # De-duplicated subset
    keep = [s for s in samples_95 if s not in dropped]
    dosage_dd = dosage_95[keep]
    stats_dd = summarise_subset(dosage_dd, f"{len(keep)}-line")

    # Delta summary table
    deltas = []
    for stat in ("silhouette_k2", "cluster1_size", "cluster2_size",
                 "mean_MAF", "median_MAF", "m_eff", "F_ST_between_clusters"):
        v95 = stats_95[stat]
        vdd = stats_dd[stat]
        try:
            pct = 100.0 * (vdd - v95) / abs(v95) if v95 != 0 else float("nan")
        except ZeroDivisionError:
            pct = float("nan")
        deltas.append({
            "stat": stat,
            "value_95line": v95,
            "value_dedup": vdd,
            "delta": vdd - v95,
            "pct_change": pct,
        })
    delta_df = pd.DataFrame(deltas)
    delta_path = TAB / "dedup_delta_summary.csv"
    delta_df.to_csv(delta_path, index=False)
    print(f"\n[delta] wrote {delta_path}")
    print(delta_df.to_string(index=False))

    print("\n[plot] rendering delta forest")
    # m_eff here is the Li & Ji count on the n x n sample Gram of this
    # subset; it is a GWAS multiple-testing quantity, not a population-
    # structure statistic, and it differs from the full-panel GWAS m_eff
    # reported in Methods (script 12, m_eff = 93). It stays in the backing
    # CSV as a record but is omitted from this structure-robustness forest
    # so the figure cannot be misread as a second, conflicting m_eff.
    forest_df = delta_df[delta_df["stat"] != "m_eff"].reset_index(drop=True)
    plot_delta_forest(forest_df, FIG / "fig_dedup_delta_forest")

    print("[plot] rendering PCA overlay")
    plot_pca_overlay(stats_95, dropped, samples_95, FIG / "fig_pca_overlay")

    print("\n=== manuscript-text summary ===")
    print(f"  duplicate components flagged at G_ij >= {DUP_G_THRESH}: {len(components)}")
    print(f"  accessions dropped                                    : {len(dropped)}")
    print(f"  panel size after de-duplication                       : {len(keep)}")
    print(f"  silhouette at k=2  : 95-line {stats_95['silhouette_k2']:.3f}"
          f"  ->  {len(keep)}-line {stats_dd['silhouette_k2']:.3f}")
    print(f"  m_eff              : {stats_95['m_eff']}  ->  {stats_dd['m_eff']}"
          f"  (Bonferroni alpha {0.05/stats_95['m_eff']:.2e} -> {0.05/stats_dd['m_eff']:.2e})")
    print(f"  F_ST btw clusters  : {stats_95['F_ST_between_clusters']:.3f}"
          f"  ->  {stats_dd['F_ST_between_clusters']:.3f}")
    print(f"  Cluster-2 size     : {stats_95['cluster2_size']}"
          f"  ->  {stats_dd['cluster2_size']}")


if __name__ == "__main__":
    main()
