#!/usr/bin/env python3
"""
Runs of Homozygosity (ROH) detection and F_ROH per accession.

ROH = contiguous tracts of homozygous genotypes; commonly used as an
individual-level inbreeding measure (Howrigan et al. 2011, BMC Genomics).
For selfers, ROH should be extensive; the panel's unexpectedly low F_IS
(0.011, script 22) is best resolved by reporting F_ROH (fraction of the
panel-anchored genome covered by long ROH segments) per accession.

Method (PLINK-style sliding window):
  - window of W consecutive SNPs slides across each AYB chromosome
  - a window is "homozygous" if it contains no more than `MAX_HET_PER_WIN`
    heterozygous sites and no more than `MAX_MISSING_PER_WIN` missing
    genotypes
  - SNPs flagged homozygous by at least `HIT_THRESHOLD` windows form a ROH
  - ROHs are merged across consecutive SNPs and required to span >=
    `MIN_LENGTH_KB` kb and >= `MIN_SNPS` SNPs

The marker density on our AYB-anchored subset (1,473 markers / 650 Mb genome
= ~1 SNP per 440 kb) is low for a fine-grained ROH scan; we report ROH
length distributions rather than precise breakpoint coordinates.

Outputs (results/34_roh/)
-------------------------
tables/
    roh_calls.csv               every ROH segment per sample
    f_roh_per_sample.csv        F_ROH per accession (4 length cutoffs)
figures/
    fig73_roh_length_distribution.png/.pdf
    fig74_f_roh_per_sample.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG, CLUSTER_PAL
from _figstyle import label_top_n, publishable_axes
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "34_roh"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
WINDOW_SNPS = 20
MAX_HET_PER_WIN = 1
MAX_MISSING_PER_WIN = 2
HIT_THRESHOLD = 0.05      # at least this fraction of containing windows
                          # must call the SNP homozygous
MIN_LENGTH_KB = 500       # minimum ROH segment length
MIN_SNPS = 8              # minimum SNPs in an ROH segment
LENGTH_CUTOFFS_MB = [0.5, 1.0, 2.0, 5.0]


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage matrix (preserve NaN; we need it for missingness in ROH)
# ---------------------------------------------------------------------------
print("[load] HapMap...")
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

call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values
dosage = dosage.loc[keep_m, keep_s]
samples = dosage.columns.tolist()
print(f"[load] dosage = {dosage.shape}")

anc = pd.read_csv(ANCHOR)
anc = anc.dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc["snp_pos_ayb"] = anc["snp_pos_ayb"].astype(int)
# restrict to anchored markers
dosage = dosage.loc[dosage.index.intersection(anc["rs"])]
anc = anc[anc["rs"].isin(dosage.index)]
print(f"[anchor] AYB-anchored markers in panel: {len(anc)}")


# ---------------------------------------------------------------------------
# Sliding-window ROH per sample
# ---------------------------------------------------------------------------
print(f"[roh] sliding-window scan (W = {WINDOW_SNPS} SNPs, "
      f"<= {MAX_HET_PER_WIN} het, <= {MAX_MISSING_PER_WIN} missing/window)...")

# per chromosome, work with the SNPs sorted by position
chrs = sorted(anc["chr_ayb"].unique(),
              key=lambda s: int(re.search(r"\d+", s).group()))
chr_to_rs = {c: anc[anc["chr_ayb"] == c].sort_values("snp_pos_ayb")
             for c in chrs}

roh_segments = []
for samp in samples:
    for c in chrs:
        sub = chr_to_rs[c]
        rs_list = sub["rs"].values
        pos_list = sub["snp_pos_ayb"].values
        D = dosage.loc[rs_list, samp].values  # 0/1/2 or NaN

        if len(D) < WINDOW_SNPS:
            continue
        # for each window starting at i, count het + missing
        homo_hits = np.zeros(len(D), dtype=int)
        windows_seen = np.zeros(len(D), dtype=int)
        for i in range(len(D) - WINDOW_SNPS + 1):
            w = D[i:i + WINDOW_SNPS]
            hets = int((w == 1).sum())
            miss = int(np.isnan(w).sum())
            ok = (hets <= MAX_HET_PER_WIN) and (miss <= MAX_MISSING_PER_WIN)
            windows_seen[i:i + WINDOW_SNPS] += 1
            if ok:
                homo_hits[i:i + WINDOW_SNPS] += 1
        # per SNP fraction of containing windows that called it homozygous
        with np.errstate(invalid="ignore", divide="ignore"):
            hit_frac = np.where(windows_seen > 0, homo_hits / windows_seen, 0)
        in_roh = hit_frac >= HIT_THRESHOLD
        # find contiguous runs
        runs = []
        in_run = False
        for k in range(len(in_roh)):
            if in_roh[k] and not in_run:
                start_idx = k; in_run = True
            elif not in_roh[k] and in_run:
                end_idx = k - 1; in_run = False
                runs.append((start_idx, end_idx))
        if in_run:
            runs.append((start_idx, len(in_roh) - 1))
        # filter by length + n_snps
        for (s_, e_) in runs:
            length_kb = (pos_list[e_] - pos_list[s_]) / 1000.0
            n_snps = e_ - s_ + 1
            if length_kb >= MIN_LENGTH_KB and n_snps >= MIN_SNPS:
                roh_segments.append({
                    "sample": samp, "chr": c,
                    "start_bp": int(pos_list[s_]),
                    "end_bp": int(pos_list[e_]),
                    "length_bp": int(pos_list[e_] - pos_list[s_]),
                    "length_kb": float(length_kb),
                    "n_snps": int(n_snps),
                })
roh_df = pd.DataFrame(roh_segments)
roh_df.to_csv(TAB / "roh_calls.csv", index=False)
print(f"[roh] {len(roh_df):,} ROH segments called, "
      f"median length {roh_df['length_kb'].median():.0f} kb")


# ---------------------------------------------------------------------------
# F_ROH per sample at multiple length cutoffs
# ---------------------------------------------------------------------------
# AYB anchored-genome length (sum of max-position per chrom is a fair proxy)
anc_len_bp = (anc.groupby("chr_ayb")["snp_pos_ayb"].max()
              - anc.groupby("chr_ayb")["snp_pos_ayb"].min()).sum()
print(f"[anchor] AYB-anchored genome span: {anc_len_bp/1e6:.1f} Mb")

f_rows = []
for samp in samples:
    samp_df = roh_df[roh_df["sample"] == samp]
    row = {"sample": samp,
           "n_roh_total": int(len(samp_df)),
           "sum_roh_bp": int(samp_df["length_bp"].sum())}
    for cut in LENGTH_CUTOFFS_MB:
        cutoff_bp = cut * 1_000_000
        sub = samp_df[samp_df["length_bp"] >= cutoff_bp]
        row[f"n_roh_ge_{cut}Mb"] = int(len(sub))
        row[f"sum_roh_ge_{cut}Mb"] = int(sub["length_bp"].sum())
        row[f"F_ROH_ge_{cut}Mb"] = float(sub["length_bp"].sum() / anc_len_bp)
    f_rows.append(row)
f_df = pd.DataFrame(f_rows)
pca = pd.read_csv(PCA_CSV).set_index("sample").reindex(samples)
f_df["cluster"] = (pca["cluster"].astype(int) + 1).values
f_df = f_df.sort_values("F_ROH_ge_1.0Mb", ascending=False)
f_df.to_csv(TAB / "f_roh_per_sample.csv", index=False)

print("\n[F_ROH] summary at length cutoffs (mean over panel):")
for cut in LENGTH_CUTOFFS_MB:
    col = f"F_ROH_ge_{cut}Mb"
    print(f"  >= {cut} Mb: mean F_ROH = {f_df[col].mean():.4f}  "
          f"range [{f_df[col].min():.4f}, {f_df[col].max():.4f}]")


# ---------------------------------------------------------------------------
# Fig 73. ROH length distribution
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
ax.hist(roh_df["length_kb"], bins=40, color=WONG["blue"],
        edgecolor="white", linewidth=0.4)
ax.set_xlabel("ROH segment length (kb)")
ax.set_ylabel("number of segments")
ax.set_title(f"ROH segment length distribution "
             f"({len(roh_df):,} segments across {len(samples)} accessions)")
ax.axvline(1000, color="black", linewidth=0.6, linestyle="--",
           label="1 Mb cutoff")
ax.axvline(5000, color="grey", linewidth=0.6, linestyle=":",
           label="5 Mb cutoff")
ax.legend()
ax.grid(True, axis="y")
save(fig, "fig73_roh_length_distribution")
print("[fig] fig73_roh_length_distribution")


# ---------------------------------------------------------------------------
# Fig 74. F_ROH per accession (sorted bar chart, coloured by cluster)
# ---------------------------------------------------------------------------
sorted_f = (f_df.sort_values("F_ROH_ge_1.0Mb", ascending=False)
              .reset_index(drop=True))
fig, ax = plt.subplots(figsize=(12, 4.6), constrained_layout=True)
x = np.arange(len(sorted_f))
colors = [CLUSTER_PAL[c - 1] for c in sorted_f["cluster"]]
ax.bar(x, sorted_f["F_ROH_ge_1.0Mb"], color=colors,
       edgecolor="black", linewidth=0.25, zorder=3)

mean_froh = float(sorted_f["F_ROH_ge_1.0Mb"].mean())
ax.axhline(mean_froh, color="black", linewidth=0.8, linestyle="--", zorder=2)

# Rank-based x-axis so the 95 accession IDs do not overlap. Per-sample IDs
# are surfaced only for the 12 most extreme accessions via adjustText.
ax.set_xticks([0, len(sorted_f) // 4, len(sorted_f) // 2,
               3 * len(sorted_f) // 4, len(sorted_f) - 1])
ax.set_xticklabels(["1", str(len(sorted_f) // 4 + 1),
                    str(len(sorted_f) // 2 + 1),
                    str(3 * len(sorted_f) // 4 + 1),
                    str(len(sorted_f))], fontsize=9)
ax.set_xlabel("Accession rank (descending F$_{ROH}$)")
ax.set_ylabel(r"F$_{ROH}$ (>= 1 Mb ROH / AYB-anchored genome span)")
ax.set_title(f"Per-accession F$_{{ROH}}$ (>= 1 Mb segments) -- selfer-paper "
             f"individual-inbreeding measure ({len(sorted_f)} accessions)")

extreme_idx = np.concatenate([np.arange(6),
                              np.arange(len(sorted_f) - 6, len(sorted_f))])
label_top_n(ax, x[extreme_idx],
            sorted_f["F_ROH_ge_1.0Mb"].values[extreme_idx],
            sorted_f["sample"].values[extreme_idx],
            n=len(extreme_idx), fontsize=7.0, color="#222222")

publishable_axes(ax, grid="y")

handles = [mpatches.Patch(color=CLUSTER_PAL[k],
                          label=f"Cluster {k+1}") for k in (0, 1)]
handles.append(plt.Line2D([0], [0], color="black", linewidth=0.8,
                          linestyle="--",
                          label=f"mean = {mean_froh:.3f}"))
ax.legend(handles=handles, loc="upper right", fontsize=8.5,
          framealpha=0.9, edgecolor="none")
save(fig, "fig74_f_roh_per_sample")
print("[fig] fig74_f_roh_per_sample")

print(f"\nOutputs in: {OUT}")
