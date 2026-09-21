#!/usr/bin/env python3
"""
ADMIXTURE Q-matrix for the 105-line AYB panel, K = 2..6.

1. Rebuild the QC-filtered dosage matrix (same cascade as scripts 01-03).
2. Write a PLINK PED + MAP and convert to BED via PLINK 1.9.
3. Run unsupervised ADMIXTURE with 5-fold CV for K = 2..6.
4. Parse the CV error from each log, pick the best-K by minimum CV.
5. Plot stacked-bar Q matrices ordered by PCA-cluster, then within cluster
   by max-component assignment.

Outputs (results/09_admixture/)
-------------------------------
plink/
    ayb.bed / ayb.bim / ayb.fam            PLINK binary input
admixture/
    ayb.K.Q  ayb.K.P  ayb.K.log            ADMIXTURE outputs per K
tables/
    admixture_cv_error.csv                 CV error vs K
    admixture_Q_K<best>.csv                Q matrix at best-K
figures/
    fig22_admixture_cv_error.png/.pdf
    fig23_admixture_Q_bars_K<best>.png/.pdf
    fig24_admixture_Q_bars_all_K.png/.pdf  multi-panel K=2..6

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG, CLUSTER_PAL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "09_admixture"
PLINK_DIR = OUT / "plink"
ADM_DIR = OUT / "admixture"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (PLINK_DIR, ADM_DIR, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

PLINK = "/Users/black_einstein/miniconda3/envs/barley_gwas/bin/plink"
ADMIXTURE = "/Users/black_einstein/miniconda3/bin/admixture"

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
K_LIST = [2, 3, 4, 5, 6]
CV_FOLDS = 5

PALETTE_K = [
    WONG["blue"], WONG["vermillion"], WONG["green"],
    WONG["orange"], WONG["purple"], WONG["yellow"],
]


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 1. Filtered dosage
# ---------------------------------------------------------------------------
print("[load] HapMap...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META]
calls = hm[sample_cols].astype(str)

# Per-marker dosage + per-marker ref/alt
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

dos_f = dosage.loc[keep_m, keep_s]
samples = dos_f.columns.tolist()
markers = dos_f.index.tolist()
print(f"[load] retained {len(markers)} markers x {len(samples)} samples")

ref_alt_f = ref_alt.loc[keep_m].copy()
ref_alt_f.index = dosage.index[keep_m]
chrom_f = hm.loc[keep_m, "chrom"].fillna("0").values
pos_f = hm.loc[keep_m, "pos"].fillna(0).astype(int).values


# ---------------------------------------------------------------------------
# 2. Write PLINK PED + MAP
# ---------------------------------------------------------------------------
print("[plink] writing PED / MAP...")
# MAP: chr rs cM bp; use numeric chrom (1..11) when available, 0 otherwise
def chrom_to_int(c):
    m = re.search(r"Vu(\d+)", str(c))
    return int(m.group(1)) if m else 0


chrom_int = np.array([chrom_to_int(c) for c in chrom_f])
map_df = pd.DataFrame({
    "chr": chrom_int,
    "rs": markers,
    "cM": 0,
    "bp": pos_f,
})
# PLINK requires monotonically increasing BP per chrom; for chr=0 (unmapped)
# assign synthetic positions so PLINK doesn't complain
unmapped_mask = chrom_int == 0
map_df.loc[unmapped_mask, "bp"] = np.arange(1, unmapped_mask.sum() + 1)
map_df.to_csv(PLINK_DIR / "ayb.map", sep="\t", header=False, index=False)


# PED: FID IID PID MID SEX PHENO geno1_A1 geno1_A2 geno2_A1 ...
def dosage_to_ped_row(sample_idx: int) -> list[str]:
    row = ["FAM", samples[sample_idx], "0", "0", "0", "-9"]
    col = dos_f.values[:, sample_idx]
    for j in range(len(markers)):
        ref = ref_alt_f.iloc[j, 0]
        alt = ref_alt_f.iloc[j, 1]
        d = col[j]
        if np.isnan(d):
            row.extend(["0", "0"])
        elif d == 0:
            row.extend([ref, ref])
        elif d == 1:
            row.extend([ref, alt])
        else:
            row.extend([alt, alt])
    return row


with open(PLINK_DIR / "ayb.ped", "w") as fh:
    for i in range(len(samples)):
        fh.write(" ".join(dosage_to_ped_row(i)) + "\n")
print(f"[plink] PED + MAP -> {PLINK_DIR}")


# ---------------------------------------------------------------------------
# 3. PED/MAP -> BED via PLINK
# ---------------------------------------------------------------------------
print("[plink] converting to BED...")
cmd = [
    PLINK, "--file", str(PLINK_DIR / "ayb"),
    "--make-bed", "--out", str(PLINK_DIR / "ayb"),
    "--allow-no-sex", "--allow-extra-chr",
]
r = subprocess.run(cmd, capture_output=True, text=True)
if r.returncode != 0:
    print("PLINK STDERR:\n", r.stderr)
    raise RuntimeError("PLINK conversion failed")
print(r.stdout.strip().splitlines()[-3:])


# ---------------------------------------------------------------------------
# 4. Run ADMIXTURE for each K (5-fold CV)
# ---------------------------------------------------------------------------
print(f"[admix] running ADMIXTURE for K in {K_LIST} (cv = {CV_FOLDS})...")
cv_rows = []
bed = PLINK_DIR / "ayb.bed"
for K in K_LIST:
    log_path = ADM_DIR / f"ayb.{K}.log"
    cmd = [ADMIXTURE, f"--cv={CV_FOLDS}", str(bed), str(K)]
    with open(log_path, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                           cwd=ADM_DIR, text=True)
    if r.returncode != 0:
        print(f"  K={K}: ADMIXTURE failed (see {log_path})")
        continue
    cv = None
    with open(log_path) as fh:
        for line in fh:
            m = re.match(r"CV error \(K=(\d+)\): ([\d.]+)", line)
            if m:
                cv = float(m.group(2))
                break
    print(f"  K={K}: CV error = {cv}")
    cv_rows.append({"K": K, "cv_error": cv})

cv_df = pd.DataFrame(cv_rows).sort_values("K")
cv_df.to_csv(TAB / "admixture_cv_error.csv", index=False)
best_K = int(cv_df.loc[cv_df["cv_error"].idxmin(), "K"])
print(f"[admix] best K by CV minimum = {best_K}")


# ---------------------------------------------------------------------------
# 5. Load Q matrices and plot
# ---------------------------------------------------------------------------
pca = pd.read_csv(PCA_TAB).set_index("sample").loc[samples]
pca_cluster = pca["cluster"].astype(int).values

def load_Q(K: int) -> np.ndarray:
    q = pd.read_csv(ADM_DIR / f"ayb.{K}.Q", sep=r"\s+", header=None).values
    return q


# Fig 22. CV error vs K
fig, ax = plt.subplots(figsize=(5.5, 3.6), constrained_layout=True)
ax.plot(cv_df["K"], cv_df["cv_error"], marker="o", linewidth=1.5,
        color=WONG["blue"])
ax.axvline(best_K, ls="--", color=WONG["vermillion"],
           label=f"best K = {best_K}")
ax.set_xlabel("K (number of ancestral populations)")
ax.set_ylabel(f"{CV_FOLDS}-fold CV error")
ax.set_title("ADMIXTURE cross-validation error vs K")
ax.set_xticks(K_LIST)
ax.grid(True)
ax.legend()
save(fig, "fig22_admixture_cv_error")
print("[fig] fig22_admixture_cv_error")


def plot_Q(ax, Q, sample_order, title=None):
    n, K = Q.shape
    Q = Q[sample_order]
    bottoms = np.zeros(n)
    for k in range(K):
        ax.bar(np.arange(n), Q[:, k], bottom=bottoms,
               width=1.0, color=PALETTE_K[k],
               edgecolor="none", linewidth=0)
        bottoms += Q[:, k]
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel("ancestry")
    if title:
        ax.set_title(title)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# Sort samples: PCA cluster first, then by Q column with highest mean within
# each cluster, so blocks are visually contiguous.
def sample_order_for_K(Q):
    order = []
    for k in sorted(set(pca_cluster)):
        idx_in_cluster = np.where(pca_cluster == k)[0]
        # within cluster, pick the Q column that's dominant and sort by it
        col = Q[idx_in_cluster].mean(axis=0).argmax()
        inner = idx_in_cluster[np.argsort(-Q[idx_in_cluster, col])]
        order.extend(inner.tolist())
    return np.array(order)


# Fig 23. Q bars at best K
Q_best = load_Q(best_K)
Q_df = pd.DataFrame(Q_best, index=samples,
                    columns=[f"Q{i+1}" for i in range(best_K)])
Q_df["pca_cluster"] = pca_cluster + 1
Q_df.to_csv(TAB / f"admixture_Q_K{best_K}.csv")

fig, ax = plt.subplots(figsize=(11, 3.0), constrained_layout=True)
order = sample_order_for_K(Q_best)
plot_Q(ax, Q_best, order,
       title=f"ADMIXTURE Q matrix, K = {best_K} (best by {CV_FOLDS}-fold CV)")
# annotate cluster-boundary tick
labels = pca_cluster[order]
boundary = np.where(np.diff(labels) != 0)[0]
for b in boundary:
    ax.axvline(b + 0.5, color="black", linewidth=1.2)
ax.set_xlabel(f"105 AYB lines, ordered by PCA cluster then by Q dominance")
save(fig, f"fig23_admixture_Q_bars_K{best_K}")
print(f"[fig] fig23_admixture_Q_bars_K{best_K}")


# Fig 24. Multi-panel K=2..max
fig, axes = plt.subplots(len(K_LIST), 1, figsize=(11, 1.7 * len(K_LIST)),
                         constrained_layout=True)
for ax, K in zip(axes, K_LIST):
    Q = load_Q(K)
    order = sample_order_for_K(Q)
    plot_Q(ax, Q, order, title=f"K = {K}")
    labels = pca_cluster[order]
    boundary = np.where(np.diff(labels) != 0)[0]
    for b in boundary:
        ax.axvline(b + 0.5, color="black", linewidth=1.0)
axes[-1].set_xlabel("105 AYB lines (per-K ordering shown)")
fig.suptitle(f"ADMIXTURE Q matrices, K = {K_LIST[0]}..{K_LIST[-1]}",
             y=1.02, fontsize=12)
save(fig, "fig24_admixture_Q_bars_all_K")
print("[fig] fig24_admixture_Q_bars_all_K")

print(f"\nOutputs in: {OUT}")
