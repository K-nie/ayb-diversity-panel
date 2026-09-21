#!/usr/bin/env python3
"""
Per-cluster private alleles, fixed-difference SNPs, and Ss05 outlier
attribution.

Script 27 (PCAdapt) and script 29 (DAPC) independently flagged
chromosome Ss05 as carrying the panel's strongest selection-genomics
signal, with 8 of the top-10 PCAdapt outliers and the top DAPC
LD1-loading SNPs sitting on Ss05. Neither scan returned an FDR-
significant hit, and the convergence on a single chromosome could
either reflect (a) a Cluster-2-driven allele-frequency shift in the
small n = 11 outlier subgroup, or (b) a panel-wide drift signature
that happens to concentrate on Ss05.

This script discriminates between the two by partitioning per-locus
allele frequencies into the PCA-derived clusters (Cluster 1, n = 84;
Cluster 2, n = 11). Private alleles are defined as SNPs whose minor
allele segregates at MAF >= 0.05 in one cluster and <= 0.01 in the
other (the standard population-genetics private-allele definition,
adapted for small-sample panels by relaxing the < 0.01 lower bound
from the strict zero-frequency convention -- at n = 11 a true zero
sits behind one sampling-variance step from the next-likely value).

Fixed-difference SNPs follow the Hahn 2018 definition (alternate-
allele frequency at 0 in one cluster and at 1 in the other). Per-
locus Weir-Cockerham F_ST is computed via the same inbred-line
simplification used in script 56 (script 56 derivation reused).

Ss05 outlier attribution: the top-50 PCAdapt and top-50 DAPC
outliers are subset to Ss05 and tagged as:

    Cluster-2-driven  if Cluster-2 MAF >= 0.40 AND Cluster-1 MAF <= 0.10
    Cluster-1-driven  if the reverse pattern holds
    panel-wide        otherwise (signal not localised to either cluster)

Outputs (results/55_private_alleles_per_cluster/)
-------------------------------------------------
tables/
    private_alleles_per_cluster.csv     per-cluster private allele lists
    fixed_difference_snps.csv            cluster-fixed-difference SNPs
    Ss05_outlier_attribution.csv         PCAdapt+DAPC Ss05 outliers tagged
    per_cluster_summary.csv              counts per category
figures/
    fig_per_locus_F_ST_manhattan.png/.pdf   per-SNP F_ST with Ss05 highlighted
    fig_Ss05_attribution_strip.png/.pdf      Ss05-strip showing C1 vs C2 MAF

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG, CLUSTER_PAL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
ANCHOR_CSV = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
PCADAPT_TOP_CSV = ROOT / "results" / "27_pcadapt_outliers" / "tables" / "pcadapt_outlier_top50.csv"
DAPC_LOADINGS_CSV = ROOT / "results" / "29_dapc" / "tables" / "dapc_marker_loadings.csv"

OUT = ROOT / "results" / "55_private_alleles_per_cluster"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MAF_MIN = 0.05

# Private-allele thresholds. The 0.05 lower bound on the "carrier"
# cluster filters out singletons; the 0.01 upper bound on the
# "near-absent" cluster allows for a single chance call without
# losing the private-allele inference at n = 11.
PRIVATE_MAF_PRESENT = 0.05
PRIVATE_MAF_ABSENT = 0.01

# Ss05 outlier-attribution thresholds. The 0.40 / 0.10 split is
# robust to the Cluster-2 sampling SD of ~ 0.07 at any single SNP
# (1/sqrt(2*11) per-allele binomial SE under p = 0.5 ~ 0.21 / 2 ~ 0.11
# at the population level), so accessions that pass these thresholds
# read as cluster-driven beyond reasonable sampling variation.
ATTR_C_HIGH = 0.40
ATTR_C_LOW = 0.10


# ---------------------------------------------------------------------------
# Working dosage matrix (mirrors scripts 01 and 56)
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

    samp_call = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_samples = samp_call[samp_call >= SAMP_CR].index.tolist()
    dosage = dosage[keep_samples]

    mark_call = dosage.notna().sum(axis=1) / dosage.shape[1]
    mean_dose = dosage.mean(axis=1, skipna=True)
    maf_all = np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)
    keep_markers = dosage.index[(mark_call >= MARK_CR) & (maf_all >= MAF_MIN)]
    dosage = dosage.loc[keep_markers]
    print(f"[qc] working matrix: markers={dosage.shape[0]}, samples={dosage.shape[1]}")
    return dosage


# ---------------------------------------------------------------------------
# Per-cluster allele frequencies and private-allele calls
# ---------------------------------------------------------------------------
def per_cluster_freqs(dosage: pd.DataFrame, cluster_labels: dict[str, int]) -> pd.DataFrame:
    """For each SNP, compute the alternate-allele frequency in Cluster 1
    (large, n = 84) and Cluster 2 (small, n = 11), plus the panel-wide MAF
    and a panel-wide minor-allele direction so per-cluster MAFs use a
    consistent allele convention."""
    samples = dosage.columns.tolist()
    # Cluster index: in the pca_coords.csv convention, cluster==1 is the
    # large 84-line group (Cluster 1) and cluster==0 is the small 11-line
    # group (Cluster 2). Translate to manuscript naming.
    idx_c1 = [i for i, s in enumerate(samples) if cluster_labels.get(s) == 1]
    idx_c2 = [i for i, s in enumerate(samples) if cluster_labels.get(s) == 0]
    print(f"[clust] Cluster 1 (large) n = {len(idx_c1)}; Cluster 2 (small) n = {len(idx_c2)}")

    X = dosage.values  # markers x samples
    # Per-cluster alt-allele frequency = (mean dosage in cluster) / 2,
    # ignoring NaN per cluster.
    def alt_freq(idx):
        block = X[:, idx]
        return np.nanmean(block, axis=1) / 2.0

    p_c1 = alt_freq(idx_c1)
    p_c2 = alt_freq(idx_c2)
    p_all = (np.nansum(X, axis=1) / (2.0 * np.sum(~np.isnan(X), axis=1)))

    # Per-cluster number of sampled alleles (chromosomes) at each SNP:
    # 2 x the count of non-missing genotypes in that cluster. Needed for
    # Hudson's F_ST sample-size correction (Bhatia et al. 2013).
    n_alleles_c1 = 2.0 * np.sum(~np.isnan(X[:, idx_c1]), axis=1)
    n_alleles_c2 = 2.0 * np.sum(~np.isnan(X[:, idx_c2]), axis=1)

    # Choose the panel-wide minor allele as the rarer one. Per-cluster
    # MAFs are taken under this same direction so "MAF in cluster c"
    # always refers to the panel-wide minor allele.
    minor_is_alt = p_all <= 0.5
    maf_c1 = np.where(minor_is_alt, p_c1, 1 - p_c1)
    maf_c2 = np.where(minor_is_alt, p_c2, 1 - p_c2)

    return pd.DataFrame({
        "rs": dosage.index,
        "p_panel": p_all,
        "minor_is_alt": minor_is_alt,
        "MAF_panel": np.minimum(p_all, 1 - p_all),
        "MAF_C1": maf_c1,
        "MAF_C2": maf_c2,
        # Raw alternate-allele frequencies and per-cluster allele counts,
        # carried for Hudson's F_ST (invariant to allele labelling).
        "p_alt_C1": p_c1,
        "p_alt_C2": p_c2,
        "n_alleles_C1": n_alleles_c1,
        "n_alleles_C2": n_alleles_c2,
    })


def call_private_alleles(freqs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (private to Cluster 1, private to Cluster 2)."""
    priv_c1 = freqs.loc[(freqs["MAF_C1"] >= PRIVATE_MAF_PRESENT)
                       & (freqs["MAF_C2"] <= PRIVATE_MAF_ABSENT)].copy()
    priv_c2 = freqs.loc[(freqs["MAF_C2"] >= PRIVATE_MAF_PRESENT)
                       & (freqs["MAF_C1"] <= PRIVATE_MAF_ABSENT)].copy()
    return priv_c1, priv_c2


def call_fixed_differences(freqs: pd.DataFrame) -> pd.DataFrame:
    """Hahn 2018 fixed-difference: alt-allele freq at 0 in one cluster and
    1 in the other (or 1 and 0 respectively). Use raw alt-allele freqs
    rather than MAFs because direction matters here."""
    samples = freqs.copy()
    samples["p_C1"] = np.where(samples["minor_is_alt"], samples["MAF_C1"],
                               1 - samples["MAF_C1"])
    samples["p_C2"] = np.where(samples["minor_is_alt"], samples["MAF_C2"],
                               1 - samples["MAF_C2"])
    fix = samples.loc[
        ((samples["p_C1"] <= 0.001) & (samples["p_C2"] >= 0.999))
        | ((samples["p_C1"] >= 0.999) & (samples["p_C2"] <= 0.001))
    ]
    return fix


def per_locus_fst(freqs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-SNP Hudson's F_ST (Bhatia et al. 2013, Genome Res. 23:1514, eq. 10),
    the estimator recommended for highly unequal sample sizes such as the
    n = 84 vs n = 11 cluster split here.

    For each SNP with alternate-allele frequencies p1, p2 in the two clusters
    and sampled-allele counts n1, n2:

        numerator   = (p1 - p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)
        denominator = p1(1-p2) + p2(1-p1)
        F_ST(locus) = numerator / denominator

    The estimator is invariant to which allele is labelled "alternate" and is
    bounded at 1, unlike the previous mis-specified Weir-Cockerham variant.
    The genome-wide estimate is the ratio of summed numerators to summed
    denominators (Bhatia's "ratio of averages"), returned alongside.

    Returns (per_locus_fst, numerator, denominator).
    """
    p1 = freqs["p_alt_C1"].to_numpy(dtype=float)
    p2 = freqs["p_alt_C2"].to_numpy(dtype=float)
    n1 = freqs["n_alleles_C1"].to_numpy(dtype=float)
    n2 = freqs["n_alleles_C2"].to_numpy(dtype=float)

    with np.errstate(invalid="ignore", divide="ignore"):
        num = ((p1 - p2) ** 2
               - p1 * (1 - p1) / (n1 - 1)
               - p2 * (1 - p2) / (n2 - 1))
        den = p1 * (1 - p2) + p2 * (1 - p1)
        fst = np.where(den > 0, num / den, np.nan)
    return fst, num, den


# ---------------------------------------------------------------------------
# Ss05 outlier attribution
# ---------------------------------------------------------------------------
def attribute_outlier(maf_c1: float, maf_c2: float) -> str:
    if maf_c2 >= ATTR_C_HIGH and maf_c1 <= ATTR_C_LOW:
        return "Cluster-2-driven"
    if maf_c1 >= ATTR_C_HIGH and maf_c2 <= ATTR_C_LOW:
        return "Cluster-1-driven"
    return "panel-wide"


def build_ss05_attribution(freqs: pd.DataFrame, fst: np.ndarray,
                            anchor: pd.DataFrame) -> pd.DataFrame:
    print("[attr] loading PCAdapt top-50 and DAPC top-50 marker loadings")
    pcadapt = pd.read_csv(PCADAPT_TOP_CSV)
    dapc = pd.read_csv(DAPC_LOADINGS_CSV)
    # DAPC loadings file is already sorted by |loading| descending in
    # script 29; take the top 50 markers.
    dapc_top50 = dapc.head(50).copy()

    # Restrict each to Ss05 outliers.
    pcadapt_ss05 = pcadapt.loc[pcadapt["chr_ayb"] == "Ss05"].copy()
    dapc_ss05 = dapc_top50.loc[dapc_top50["chr_ayb"] == "Ss05"].copy()
    pcadapt_ss05["method"] = "PCAdapt"
    dapc_ss05["method"] = "DAPC"

    freq_lookup = freqs.set_index("rs")[["MAF_C1", "MAF_C2"]].to_dict("index")
    fst_lookup = dict(zip(freqs["rs"].values, fst))

    rows = []
    for src in (pcadapt_ss05, dapc_ss05):
        for _, hit in src.iterrows():
            rs = hit["rs"]
            f = freq_lookup.get(rs)
            if f is None:
                continue
            attribution = attribute_outlier(f["MAF_C1"], f["MAF_C2"])
            rows.append({
                "rs": rs,
                "method": hit["method"],
                "chr_ayb": hit["chr_ayb"],
                "snp_pos_ayb": hit["snp_pos_ayb"],
                "MAF_C1": f["MAF_C1"],
                "MAF_C2": f["MAF_C2"],
                "F_ST": fst_lookup.get(rs, np.nan),
                "attribution": attribution,
            })
    out = pd.DataFrame(rows).sort_values(["method", "snp_pos_ayb"])
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_fst_manhattan(freqs: pd.DataFrame, fst: np.ndarray,
                       anchor: pd.DataFrame, out_path: Path) -> None:
    """Per-SNP F_ST Manhattan with Ss05 highlighted in vermillion and the
    top Ss05 outliers labelled with their nearest Funannotate gene symbol
    (read across all candidate_genes_*.csv tables in 21_candidate_genes_ayb)."""
    from _figstyle import adjust_labels, publishable_axes

    df = freqs[["rs"]].copy()
    df["F_ST"] = fst
    df = df.merge(anchor[["rs", "chr_ayb", "snp_pos_ayb"]], on="rs", how="left")
    df = df.dropna(subset=["chr_ayb", "snp_pos_ayb"])
    df["chr_ayb"] = df["chr_ayb"].astype(str)

    chrom_order = [f"Ss{i:02d}" for i in range(1, 12)]
    df = df.loc[df["chr_ayb"].isin(chrom_order)]
    df = df.sort_values(["chr_ayb", "snp_pos_ayb"]).reset_index(drop=True)

    chrom_offsets = {}
    cur = 0
    for c in chrom_order:
        chrom_offsets[c] = cur
        cur += df.loc[df["chr_ayb"] == c, "snp_pos_ayb"].max() if (df["chr_ayb"] == c).any() else 0
        cur += 5_000_000  # gap between chromosomes
    df["x"] = df.apply(lambda r: chrom_offsets[r["chr_ayb"]] + r["snp_pos_ayb"], axis=1)

    # Pull gene symbols for the top-10 Ss05 outliers. First try the existing
    # candidate-gene tables (paper A1 / A2); fall back to a direct nearest-
    # neighbour lookup against the Funannotate GFF for SNPs that are
    # structure-driven and therefore absent from GWAS-thresholded tables.
    cand_dir = ROOT / "results" / "21_candidate_genes_ayb" / "tables"
    gene_lookup: dict[str, str] = {}
    if cand_dir.exists():
        for cand_csv in sorted(cand_dir.glob("candidate_genes_*.csv")):
            try:
                cg = pd.read_csv(cand_csv)
            except Exception:
                continue
            cg = cg.loc[cg["snp_chr"] == "Ss05"].copy()
            if cg.empty:
                continue
            cg = cg.sort_values(["snp_rs", "abs_dist_bp"])
            for rs, sub in cg.groupby("snp_rs"):
                if rs in gene_lookup:
                    continue
                pick = sub.iloc[0]
                nm = str(pick.get("name", "") or "").strip()
                if nm and not nm.startswith("Spste."):
                    gene_lookup[rs] = nm
                else:
                    prod = str(pick.get("product", "") or "").strip()
                    if prod and prod.lower() != "hypothetical protein":
                        gene_lookup[rs] = prod[:24]

    # Funannotate GFF fallback for SNPs not picked up by candidate-gene
    # tables (high-F_ST outliers whose GWAS p sits above the threshold).
    gff_path = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
    if gff_path.exists():
        import re
        attr_re = re.compile(r"(\w+)=([^;]+)")
        ss05_genes = []
        mrna_product: dict[str, str] = {}
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[0] != "Ss05":
                    continue
                attrs = dict(attr_re.findall(f[8]))
                if f[2] == "gene":
                    ss05_genes.append({"start": int(f[3]), "end": int(f[4]),
                                        "gene_id": attrs.get("ID", ""),
                                        "name": attrs.get("Name", "")})
                elif f[2] == "mRNA":
                    pid = attrs.get("Parent", "")
                    if pid and pid not in mrna_product:
                        mrna_product[pid] = attrs.get("product", "")
        gene_df = pd.DataFrame(ss05_genes)
        gene_df["product"] = gene_df["gene_id"].map(mrna_product).fillna("")

        def nearest_gene(pos: float) -> str | None:
            if gene_df.empty:
                return None
            mid = (gene_df["start"] + gene_df["end"]) / 2.0
            dist = (mid - pos).abs()
            # Prefer the nearest *named* (non-Spste.*) gene within +/-500 kb.
            window = gene_df.loc[dist <= 500_000].copy()
            window["dist"] = dist.loc[window.index]
            named = window.loc[~window["name"].str.startswith("Spste.")]
            if not named.empty:
                pick = named.sort_values("dist").iloc[0]
                return str(pick["name"]).strip()
            # Fall back to the closest gene's compact Spste tag
            # (Spste.TSs11.05G175210.1 -> Ss05g175210) so the label is at
            # least informative.
            i = int(dist.idxmin())
            nm = str(gene_df.iloc[i]["name"]).strip()
            m = re.match(r"Spste\.\w+\.(\d+)G(\d+)", nm)
            if m:
                return f"Ss{int(m.group(1)):02d}g{m.group(2)}"
            return nm or None
    else:
        nearest_gene = lambda _p: None

    fig, ax = plt.subplots(figsize=(12, 4.8), constrained_layout=True)
    for c in chrom_order:
        sub = df.loc[df["chr_ayb"] == c]
        if len(sub) == 0:
            continue
        colour = WONG["vermillion"] if c == "Ss05" else (WONG["blue"] if (int(c[2:]) % 2) else WONG["skyblue"])
        ax.scatter(sub["x"], sub["F_ST"], s=10, color=colour, alpha=0.85,
                   zorder=2)

    # Gene-label callouts on the top-10 Ss05 outliers
    ss05_top = (df.loc[df["chr_ayb"] == "Ss05"]
                  .sort_values("F_ST", ascending=False)
                  .drop_duplicates(subset=["rs"])
                  .head(10))
    label_texts = []
    seen_genes: set[str] = set()
    for _, r in ss05_top.iterrows():
        gene = gene_lookup.get(r["rs"]) or nearest_gene(r["snp_pos_ayb"])
        if not gene or gene in seen_genes:
            continue
        seen_genes.add(gene)
        label_texts.append(ax.text(r["x"], r["F_ST"], gene,
                                    fontsize=7.5, color=WONG["vermillion"],
                                    fontstyle="italic", zorder=10))
    adjust_labels(label_texts, ax=ax,
                   expand=(1.2, 1.5), force_text=(0.8, 1.2))

    # Chromosome midpoint labels.
    ticks = []
    for c in chrom_order:
        sub = df.loc[df["chr_ayb"] == c]
        if len(sub) > 0:
            ticks.append((c, (sub["x"].min() + sub["x"].max()) / 2))
    ax.set_xticks([t[1] for t in ticks])
    ax.set_xticklabels([t[0] for t in ticks])
    ax.set_ylabel(r"Per-locus Hudson F$_{ST}$ (Cluster 1 vs Cluster 2)")
    ax.set_xlabel("S. stenocarpa pseudo-chromosome (AYB-anchored markers)")
    ax.set_title(r"Per-locus Hudson F$_{ST}$ across AYB-anchored markers; "
                 "Ss05 highlighted with top-10 nearest-gene callouts")
    publishable_axes(ax, grid="y")
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_ss05_attribution_strip(ss05: pd.DataFrame, out_path: Path) -> None:
    """Per-cluster frequency of the panel-wide minor allele for each Ss05
    outlier. Cluster-2-driven SNPs cluster in the upper-left quadrant
    (allele rare in Cluster 1, common in Cluster 2); Cluster-1-driven SNPs
    in the lower-right; panel-wide signals along the diagonal."""
    fig, ax = plt.subplots(figsize=(8.0, 6.2))

    # Shaded attribution regions across the full [0, 1] square.
    ax.add_patch(plt.Rectangle((0, ATTR_C_HIGH), ATTR_C_LOW, 1 - ATTR_C_HIGH,
                               facecolor=WONG["vermillion"], alpha=0.10,
                               zorder=0))
    ax.text(ATTR_C_LOW / 2, (ATTR_C_HIGH + 1.0) / 2 - 0.05,
            "Cluster-2-driven\nregion", ha="center", va="center",
            fontsize=9, color=WONG["vermillion"], alpha=0.75)

    ax.add_patch(plt.Rectangle((ATTR_C_HIGH, 0), 1 - ATTR_C_HIGH, ATTR_C_LOW,
                               facecolor=WONG["blue"], alpha=0.10, zorder=0))
    ax.text((ATTR_C_HIGH + 1.0) / 2 - 0.05, ATTR_C_LOW / 2,
            "Cluster-1-driven\nregion", ha="center", va="center",
            fontsize=9, color=WONG["blue"], alpha=0.75)

    colour_map = {"Cluster-2-driven": WONG["vermillion"],
                  "Cluster-1-driven": WONG["blue"],
                  "panel-wide": WONG["yellow"]}
    marker_map = {"PCAdapt": "o", "DAPC": "s"}

    for method, marker in marker_map.items():
        sub = ss05.loc[ss05["method"] == method]
        for tag, colour in colour_map.items():
            ssub = sub.loc[sub["attribution"] == tag]
            if len(ssub) == 0:
                continue
            ax.scatter(ssub["MAF_C1"], ssub["MAF_C2"],
                       color=colour, marker=marker, s=80,
                       edgecolor="black", linewidth=0.6,
                       label=f"{method} - {tag} (n={len(ssub)})",
                       alpha=0.9, zorder=3)

    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=0.7,
            label="equal between clusters", zorder=1)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Freq of panel-wide minor allele in Cluster 1 (n = 84)")
    ax.set_ylabel("Freq of panel-wide minor allele in Cluster 2 (n = 11)")
    ax.set_title("Ss05 outlier attribution: PCAdapt + DAPC top-50 hits on Ss05")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    dosage = build_working_matrix()

    print("[clust] loading PCA cluster labels")
    pca = pd.read_csv(PCA_CSV)
    cluster_labels = dict(zip(pca["sample"], pca["cluster"]))

    print("[freq] computing per-cluster allele frequencies")
    freqs = per_cluster_freqs(dosage, cluster_labels)
    print(f"  computed for {len(freqs)} SNPs")

    print("[call] private allele detection")
    priv_c1, priv_c2 = call_private_alleles(freqs)
    print(f"  private to Cluster 1 (large): {len(priv_c1)}")
    print(f"  private to Cluster 2 (small): {len(priv_c2)}")

    print("[call] fixed-difference SNPs (Hahn 2018)")
    fix = call_fixed_differences(freqs)
    print(f"  fixed-difference SNPs: {len(fix)}")

    print("[fst] per-locus Hudson F_ST (Bhatia et al. 2013)")
    fst, fst_num, fst_den = per_locus_fst(freqs)
    genome_fst = np.nansum(fst_num) / np.nansum(fst_den)
    print(f"  per-locus mean F_ST = {np.nanmean(fst):.3f}, "
          f"median = {np.nanmedian(fst):.3f}, max = {np.nanmax(fst):.3f}")
    print(f"  genome-wide Hudson F_ST (ratio of averages) = {genome_fst:.3f}")

    print("[anchor] joining to AYB anchoring table")
    anchor = pd.read_csv(ANCHOR_CSV)
    anchor = anchor.rename(columns={"rs": "rs"})

    # Write per-cluster private allele tables, joined to anchoring info.
    priv_c1 = priv_c1.merge(anchor[["rs", "chr_ayb", "snp_pos_ayb"]], on="rs", how="left")
    priv_c2 = priv_c2.merge(anchor[["rs", "chr_ayb", "snp_pos_ayb"]], on="rs", how="left")
    fix = fix.merge(anchor[["rs", "chr_ayb", "snp_pos_ayb"]], on="rs", how="left")

    priv_c1.to_csv(TAB / "private_alleles_C1.csv", index=False)
    priv_c2.to_csv(TAB / "private_alleles_C2.csv", index=False)
    fix.to_csv(TAB / "fixed_difference_snps.csv", index=False)

    # Per-cluster summary counts.
    summary_rows = [
        {"category": "Private to Cluster 1 (large, n=84)", "count": len(priv_c1),
         "fraction_of_panel": len(priv_c1) / len(freqs)},
        {"category": "Private to Cluster 2 (small, n=11)", "count": len(priv_c2),
         "fraction_of_panel": len(priv_c2) / len(freqs)},
        {"category": "Fixed-difference SNPs (Hahn 2018)", "count": len(fix),
         "fraction_of_panel": len(fix) / len(freqs)},
    ]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TAB / "per_cluster_summary.csv", index=False)
    print("\n[summary]")
    print(summary.to_string(index=False))

    # Ss05 attribution
    print("\n[ss05] outlier attribution")
    ss05 = build_ss05_attribution(freqs, fst, anchor)
    ss05.to_csv(TAB / "Ss05_outlier_attribution.csv", index=False)
    print(f"  Ss05 outliers tagged: {len(ss05)}")
    print(ss05.groupby(["method", "attribution"]).size().rename("n").reset_index().to_string(index=False))

    # Headline: how many of the PCAdapt + DAPC top-50 Ss05 hits are
    # Cluster-2-driven?
    n_c2_driven = int((ss05["attribution"] == "Cluster-2-driven").sum())
    n_c1_driven = int((ss05["attribution"] == "Cluster-1-driven").sum())
    n_panel_wide = int((ss05["attribution"] == "panel-wide").sum())
    n_total = len(ss05)
    pct_c2 = 100.0 * n_c2_driven / n_total if n_total else float("nan")

    # Figures
    print("\n[plot] per-locus F_ST Manhattan")
    plot_fst_manhattan(freqs, fst, anchor, FIG / "fig_per_locus_F_ST_manhattan")
    print("[plot] Ss05 attribution strip")
    plot_ss05_attribution_strip(ss05, FIG / "fig_Ss05_attribution_strip")

    print("\n=== manuscript-text summary ===")
    print(f"  total SNPs analysed                 : {len(freqs)}")
    print(f"  private alleles to Cluster 1 (n=84) : {len(priv_c1)}")
    print(f"  private alleles to Cluster 2 (n=11) : {len(priv_c2)}")
    print(f"  fixed-difference SNPs (Hahn 2018)   : {len(fix)}")
    print(f"  Ss05 outliers (PCAdapt + DAPC)      : {n_total}")
    print(f"    Cluster-2-driven                   : {n_c2_driven} ({pct_c2:.1f}%)")
    print(f"    Cluster-1-driven                   : {n_c1_driven}")
    print(f"    panel-wide                         : {n_panel_wide}")


if __name__ == "__main__":
    main()
