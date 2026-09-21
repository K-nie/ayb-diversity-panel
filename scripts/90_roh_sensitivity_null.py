#!/usr/bin/env python3
"""
90 - ROH robustness: length-class degeneracy, zero-het sensitivity, and a
within-locus permutation null.

Reviewer 1 major point 2. At the panel's anchored marker density the ROH
parameters do not do what the submitted Methods claim:

  * The retained ROH have a *median* length of ~8.7 Mb and a minimum of
    ~0.58 Mb; the 20-SNP sliding window forces every call to span >= 20
    markers, so the stated ">= 500 kb AND >= 8 SNP" retention filter is
    inoperative (the 20-SNP window already dominates both bounds).
  * F_ROH at the 0.5, 1.0 and 2.0 Mb cutoffs are numerically identical
    because essentially no segments fall between 0.5 and 2 Mb, so the
    four-cutoff "ancient vs recent autozygosity" length-class inference is
    not resolvable and must be dropped / demoted.

This script quantifies all three points and adds two robustness checks the
reviewer asked for:

  1. Zero-het sensitivity: re-call ROH with max_het_per_window = 0 and show
     the small-cluster vs large-cluster F_ROH contrast survives.
  2. Within-locus permutation null: permute each SNP's genotypes across
     accessions (preserving per-locus MAF and heterozygote count exactly,
     destroying along-chromosome autozygosity correlation) and re-call ROH.
     If observed F_ROH >> null, the ROH reflect genuine contiguous
     autozygosity, not chance homozygosity clustering at sparse density.

Outputs -> results/90_roh_sensitivity_null/.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"

OUT = ROOT / "results" / "90_roh_sensitivity_null"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

# QC + ROH parameters identical to script 34.
SAMP_CR, MARK_CR, MIN_MAF = 0.90, 0.90, 0.05
WINDOW, MAX_MISS, HIT_THR = 20, 2, 0.05
MIN_LEN_KB, MIN_SNPS = 500, 8
CUTOFFS_MB = [0.5, 1.0, 2.0, 5.0]
N_PERM = 100
RNG = np.random.default_rng(20260920)


def build_anchored_matrix():
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    scols = [c for c in hm.columns if c not in META]
    calls = hm[scols].astype(str)
    ref_alt = hm["alleles"].str.split("/", expand=True)
    dosage = np.full((len(hm), len(scols)), np.nan)
    for i in range(len(hm)):
        ref, alt = ref_alt.iloc[i]
        arr = calls.iloc[i].values
        dosage[i, arr == ref + ref] = 0.0
        dosage[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
        dosage[i, arr == alt + alt] = 2.0
    dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=scols)
    cr_m = dosage.notna().sum(axis=1) / dosage.shape[1]
    maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                     1 - dosage.mean(axis=1, skipna=True) / 2.0)
    cr_s = dosage.notna().sum(axis=0) / dosage.shape[0]
    dosage = dosage.loc[((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values,
                        (cr_s >= SAMP_CR).values]
    anc = pd.read_csv(ANCHOR).dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
    anc["snp_pos_ayb"] = anc["snp_pos_ayb"].astype(int)
    dosage = dosage.loc[dosage.index.intersection(anc["rs"])]
    anc = anc[anc["rs"].isin(dosage.index)]
    return dosage, anc


def froh_sum_per_sample(D_by_chrom, pos_by_chrom, max_het):
    """Vectorised PLINK-style ROH caller. D_by_chrom: list of (n_snp, n_samp)
    arrays; returns (sum_roh_bp[n_samp], list of all segment lengths_bp)."""
    n_samp = D_by_chrom[0].shape[1]
    sum_bp = np.zeros(n_samp)
    all_len = []
    for D, pos in zip(D_by_chrom, pos_by_chrom):
        n = D.shape[0]
        if n < WINDOW:
            continue
        het = (D == 1)
        miss = np.isnan(D)
        # window sums via cumsum along SNP axis
        het_c = np.vstack([np.zeros((1, n_samp)), np.cumsum(het, axis=0)])
        mis_c = np.vstack([np.zeros((1, n_samp)), np.cumsum(miss, axis=0)])
        n_win = n - WINDOW + 1
        het_win = het_c[WINDOW:WINDOW + n_win] - het_c[:n_win]
        mis_win = mis_c[WINDOW:WINDOW + n_win] - mis_c[:n_win]
        ok = (het_win <= max_het) & (mis_win <= MAX_MISS)  # (n_win, n_samp)
        # coverage counts per SNP (same for all samples): windows covering j
        cov = np.zeros(n + 1)
        cov[:n_win] += 1
        cov[WINDOW:WINDOW + n_win] -= 1
        cov = np.cumsum(cov)[:n]
        cov[cov == 0] = 1
        for s in range(n_samp):
            diff = np.zeros(n + 1)
            oks = ok[:, s]
            idx = np.nonzero(oks)[0]
            if idx.size == 0:
                continue
            diff[idx] += 1
            diff[idx + WINDOW] -= 1
            hits = np.cumsum(diff)[:n]
            in_roh = (hits / cov) >= HIT_THR
            # contiguous runs
            padded = np.concatenate([[False], in_roh, [False]])
            edges = np.diff(padded.astype(int))
            starts = np.nonzero(edges == 1)[0]
            ends = np.nonzero(edges == -1)[0] - 1
            for a, b in zip(starts, ends):
                length_bp = pos[b] - pos[a]
                if (length_bp / 1000.0) >= MIN_LEN_KB and (b - a + 1) >= MIN_SNPS:
                    sum_bp[s] += length_bp
                    all_len.append(length_bp)
    return sum_bp, all_len


def main():
    print("[load] building anchored matrix")
    dosage, anc = build_anchored_matrix()
    samples = dosage.columns.tolist()
    print(f"  anchored markers = {len(anc)}, samples = {len(samples)}")

    chrs = sorted(anc["chr_ayb"].unique(), key=lambda s: int(re.search(r"\d+", s).group()))
    D_by_chrom, pos_by_chrom = [], []
    for c in chrs:
        sub = anc[anc["chr_ayb"] == c].sort_values("snp_pos_ayb")
        D_by_chrom.append(dosage.loc[sub["rs"].values, samples].to_numpy(dtype=float))
        pos_by_chrom.append(sub["snp_pos_ayb"].to_numpy(dtype=float))

    anc_len_bp = float((anc.groupby("chr_ayb")["snp_pos_ayb"].max()
                        - anc.groupby("chr_ayb")["snp_pos_ayb"].min()).sum())

    pca = pd.read_csv(PCA_CSV).set_index("sample").reindex(samples)
    # pca cluster==1 = large (n=84); ==0 = small (n=11). Label plainly.
    is_small = (pca["cluster"].astype(int) == 0).to_numpy()

    # ---- spacing summary ----
    gaps = []
    for pos in pos_by_chrom:
        if len(pos) > 1:
            gaps.extend(np.diff(np.sort(pos)))
    gaps = np.array(gaps)
    print(f"\n[spacing] inter-SNP: mean={gaps.mean()/1e3:.0f} kb, "
          f"median={np.median(gaps)/1e3:.0f} kb; anchored span={anc_len_bp/1e6:.1f} Mb")

    # ---- observed (het<=1) ----
    print("\n[observed] het<=1 ROH call")
    sum_obs, len_obs = froh_sum_per_sample(D_by_chrom, pos_by_chrom, max_het=1)
    froh_obs = sum_obs / anc_len_bp
    len_obs = np.array(len_obs) / 1e6
    print(f"  segments={len(len_obs)}, length Mb: min={len_obs.min():.2f} "
          f"median={np.median(len_obs):.2f} max={len_obs.max():.2f}")
    print(f"  segments <1 Mb: {(len_obs<1).sum()}, <2 Mb: {(len_obs<2).sum()}")
    print(f"  panel mean F_ROH={froh_obs.mean():.4f}; "
          f"small-cluster={froh_obs[is_small].mean():.4f}, "
          f"large-cluster={froh_obs[~is_small].mean():.4f}")

    # ---- F_ROH cutoff degeneracy ----
    print("\n[degeneracy] F_ROH at length cutoffs (recomputed from segments)")
    # re-run to capture per-sample lengths for cutoff-based F_ROH
    def froh_at_cutoffs(max_het):
        per = {c: np.zeros(len(samples)) for c in CUTOFFS_MB}
        for D, pos in zip(D_by_chrom, pos_by_chrom):
            n = D.shape[0]
            if n < WINDOW:
                continue
            het = (D == 1); miss = np.isnan(D)
            het_c = np.vstack([np.zeros((1, len(samples))), np.cumsum(het, 0)])
            mis_c = np.vstack([np.zeros((1, len(samples))), np.cumsum(miss, 0)])
            n_win = n - WINDOW + 1
            ok = ((het_c[WINDOW:WINDOW+n_win]-het_c[:n_win]) <= max_het) & \
                 ((mis_c[WINDOW:WINDOW+n_win]-mis_c[:n_win]) <= MAX_MISS)
            cov = np.zeros(n+1); cov[:n_win]+=1; cov[WINDOW:WINDOW+n_win]-=1
            cov = np.cumsum(cov)[:n]; cov[cov==0]=1
            for s in range(len(samples)):
                idx = np.nonzero(ok[:, s])[0]
                if idx.size == 0: continue
                diff = np.zeros(n+1); diff[idx]+=1; diff[idx+WINDOW]-=1
                in_roh = (np.cumsum(diff)[:n]/cov) >= HIT_THR
                padded = np.concatenate([[False], in_roh, [False]])
                edges = np.diff(padded.astype(int))
                for a, b in zip(np.nonzero(edges==1)[0], np.nonzero(edges==-1)[0]-1):
                    L = pos[b]-pos[a]
                    if (L/1000.0)>=MIN_LEN_KB and (b-a+1)>=MIN_SNPS:
                        for cut in CUTOFFS_MB:
                            if L >= cut*1e6:
                                per[cut][s] += L
        return {c: per[c]/anc_len_bp for c in CUTOFFS_MB}
    froh_cut = froh_at_cutoffs(1)
    for cut in CUTOFFS_MB:
        print(f"  >= {cut} Mb: panel mean F_ROH = {froh_cut[cut].mean():.4f}")

    # ---- zero-het sensitivity ----
    print("\n[zero-het] het==0 ROH call")
    sum_zh, len_zh = froh_sum_per_sample(D_by_chrom, pos_by_chrom, max_het=0)
    froh_zh = sum_zh / anc_len_bp
    print(f"  panel mean F_ROH={froh_zh.mean():.4f}; "
          f"small-cluster={froh_zh[is_small].mean():.4f}, "
          f"large-cluster={froh_zh[~is_small].mean():.4f}")

    # ---- permutation null ----
    print(f"\n[null] within-locus permutation, {N_PERM} perms")
    null_panel = np.zeros(N_PERM)
    null_small = np.zeros(N_PERM)
    null_large = np.zeros(N_PERM)
    for p in range(N_PERM):
        perm_D = []
        for D in D_by_chrom:
            Dp = D.copy()
            for r in range(Dp.shape[0]):
                Dp[r] = RNG.permutation(Dp[r])
            perm_D.append(Dp)
        s_bp, _ = froh_sum_per_sample(perm_D, pos_by_chrom, max_het=1)
        fr = s_bp / anc_len_bp
        null_panel[p] = fr.mean()
        null_small[p] = fr[is_small].mean()
        null_large[p] = fr[~is_small].mean()
        if (p+1) % 25 == 0:
            print(f"  perm {p+1}/{N_PERM} done")
    p_panel = (np.sum(null_panel >= froh_obs.mean()) + 1) / (N_PERM + 1)
    p_small = (np.sum(null_small >= froh_obs[is_small].mean()) + 1) / (N_PERM + 1)
    print(f"  null panel F_ROH: mean={null_panel.mean():.4f} "
          f"(95th pct {np.percentile(null_panel,95):.4f}); "
          f"observed={froh_obs.mean():.4f}; empirical p={p_panel:.4f}")
    print(f"  null small-cluster F_ROH: mean={null_small.mean():.4f}; "
          f"observed={froh_obs[is_small].mean():.4f}; empirical p={p_small:.4f}")

    # ---- write tables ----
    pd.DataFrame({
        "metric": ["anchored_markers", "samples", "mean_spacing_kb",
                   "median_spacing_kb", "anchored_span_Mb",
                   "n_segments", "min_seg_Mb", "median_seg_Mb",
                   "seg_lt_1Mb", "seg_lt_2Mb",
                   "froh_panel_obs", "froh_small_obs", "froh_large_obs",
                   "froh_ge0.5", "froh_ge1.0", "froh_ge2.0", "froh_ge5.0",
                   "froh_panel_zerohet", "froh_small_zerohet", "froh_large_zerohet",
                   "null_panel_mean", "null_panel_95pct", "p_panel",
                   "null_small_mean", "p_small"],
        "value": [len(anc), len(samples), round(gaps.mean()/1e3,1),
                  round(np.median(gaps)/1e3,1), round(anc_len_bp/1e6,1),
                  len(len_obs), round(len_obs.min(),2), round(np.median(len_obs),2),
                  int((len_obs<1).sum()), int((len_obs<2).sum()),
                  round(froh_obs.mean(),4), round(froh_obs[is_small].mean(),4),
                  round(froh_obs[~is_small].mean(),4),
                  round(froh_cut[0.5].mean(),4), round(froh_cut[1.0].mean(),4),
                  round(froh_cut[2.0].mean(),4), round(froh_cut[5.0].mean(),4),
                  round(froh_zh.mean(),4), round(froh_zh[is_small].mean(),4),
                  round(froh_zh[~is_small].mean(),4),
                  round(null_panel.mean(),4), round(np.percentile(null_panel,95),4),
                  round(p_panel,4), round(null_small.mean(),4), round(p_small,4)],
    }).to_csv(TAB / "roh_sensitivity_summary.csv", index=False)

    pd.DataFrame({"sample": samples, "is_small_cluster": is_small,
                  "F_ROH_obs": froh_obs, "F_ROH_zerohet": froh_zh}).to_csv(
        TAB / "froh_per_sample_obs_vs_zerohet.csv", index=False)

    # ---- figure: null vs observed ----
    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    ax.hist(null_panel, bins=25, color="#999999", edgecolor="white",
            label=f"permutation null (n={N_PERM})")
    ax.axvline(froh_obs.mean(), color="#D55E00", linewidth=2,
               label=f"observed panel F_ROH = {froh_obs.mean():.3f}")
    ax.set_xlabel("Panel mean F$_{ROH}$")
    ax.set_ylabel("permutations")
    ax.set_title("Observed autozygosity vs within-locus permutation null")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "fig_roh_permutation_null.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "fig_roh_permutation_null.pdf", bbox_inches="tight")
    plt.close(fig)
    print("\n[write] tables + fig_roh_permutation_null")


if __name__ == "__main__":
    main()
