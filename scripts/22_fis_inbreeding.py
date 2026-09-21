#!/usr/bin/env python3
"""
Wright's F_IS (inbreeding coefficient) per locus, globally, and per accession
for the 95-line AYB panel.

For a selfer such as Sphenostylis stenocarpa, heterozygosity deficit is the
expected diagnostic signature. Three measures are reported:

1. **Per-locus F_IS**: 1 - H_obs / H_exp where H_exp = 2pq (Hardy-Weinberg)
   and H_obs is the observed heterozygote frequency at each locus.
2. **Global F_IS**: weighted ratio (sum H_exp - sum H_obs) / sum H_exp,
   the Weir 1996 ratio-of-sums form preferred over per-locus averaging.
   Bootstrap 95 % CI by resampling loci.
3. **Per-accession F**: 1 - H_obs_individual / E[H_exp_individual], the
   individual inbreeding coefficient (Yang et al. 2010 method I; standard
   GREML inbreeding). Useful to flag low-heterozygosity lines (highly
   selfed) vs high-heterozygosity lines (recently crossed).

Outputs (results/22_fis_inbreeding/)
------------------------------------
tables/
    fis_per_locus.csv
    fis_global_summary.csv
    f_per_accession.csv
figures/
    fig56_fis_distribution.png/.pdf
    fig57_f_per_accession.png/.pdf
    fig58_fis_manhattan.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG, CLUSTER_PAL
from _figstyle import label_top_n, publishable_axes
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "22_fis_inbreeding"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_BOOT = 1000
RNG = np.random.default_rng(20250517)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage with NaN preserved (don't impute — F_IS depends on missing)
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
maf_pre = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                     1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf_pre >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values

M = dosage.loc[keep_m, keep_s].values   # markers x samples, NaN preserved
samples = list(dosage.columns[keep_s])
markers = list(dosage.index[keep_m])
chrom_raw = hm.loc[keep_m, "chrom"].values
pos_raw = pd.to_numeric(hm.loc[keep_m, "pos"], errors="coerce").values
print(f"[load] {len(markers)} markers x {len(samples)} samples")


# ---------------------------------------------------------------------------
# Per-locus H_obs, H_exp, F_IS
# ---------------------------------------------------------------------------
print("[fis] computing per-locus H_obs, H_exp, F_IS...")
nm = ~np.isnan(M)
n_called = nm.sum(axis=1)
p_alt = np.nansum(M, axis=1) / (2.0 * np.maximum(n_called, 1))
p_alt = np.clip(p_alt, 0.0, 1.0)
H_exp = 2.0 * p_alt * (1.0 - p_alt)
H_obs = (M == 1).sum(axis=1) / np.maximum(n_called, 1)
with np.errstate(invalid="ignore", divide="ignore"):
    F_IS = np.where(H_exp > 0, 1.0 - H_obs / H_exp, np.nan)

per_locus = pd.DataFrame({
    "rs": markers, "chrom": chrom_raw, "pos": pos_raw,
    "n_called": n_called, "maf": np.minimum(p_alt, 1 - p_alt),
    "H_exp": H_exp, "H_obs": H_obs, "F_IS": F_IS,
})
per_locus.to_csv(TAB / "fis_per_locus.csv", index=False)
ok = np.isfinite(F_IS)
print(f"[fis] {ok.sum():,} loci with finite F_IS; median F_IS = "
      f"{np.nanmedian(F_IS):.3f}, mean = {np.nanmean(F_IS):.3f}")

# Global F_IS (ratio of sums, Weir 1996)
global_fis = (H_exp[ok].sum() - H_obs[ok].sum()) / H_exp[ok].sum()
print(f"[fis] global F_IS = {global_fis:.4f}")

# Bootstrap 95 % CI by resampling loci
boot = np.empty(N_BOOT)
idx_ok = np.where(ok)[0]
for b in range(N_BOOT):
    s = RNG.choice(idx_ok, size=len(idx_ok), replace=True)
    boot[b] = (H_exp[s].sum() - H_obs[s].sum()) / H_exp[s].sum()
lo, hi = np.percentile(boot, [2.5, 97.5])
print(f"[fis] bootstrap 95 % CI = [{lo:.4f}, {hi:.4f}]")

pd.DataFrame([{
    "estimator": "Wright 1969 F_IS (ratio of sums)",
    "global_FIS": global_fis,
    "boot_ci95_low": lo, "boot_ci95_high": hi,
    "n_loci_used": int(ok.sum()),
    "n_samples": int(len(samples)),
    "interpretation": ("Sphenostylis stenocarpa is highly autogamous (selfer); "
                       "F_IS > 0.5 is expected and indicates strong "
                       "heterozygosity deficit relative to Hardy-Weinberg."),
}]).to_csv(TAB / "fis_global_summary.csv", index=False)


# ---------------------------------------------------------------------------
# Per-accession inbreeding F (proportion of expected het that's missing)
# ---------------------------------------------------------------------------
print("[F-ind] computing per-accession F...")
# For each sample: H_obs_i = mean(het site count) / n_loci_called
#                  H_exp_i = mean(2 p (1-p)) over loci where sample is called
het_counts = np.where(M == 1, 1.0, 0.0)
het_counts[~nm] = np.nan
H_obs_i = np.nanmean(het_counts, axis=0)            # per sample
# expected het per sample = mean of 2 p (1-p) over loci where called
exp_per_locus = 2.0 * p_alt * (1.0 - p_alt)
H_exp_i = np.array([np.nanmean(exp_per_locus[nm[:, i]]) for i in range(len(samples))])
F_ind = 1.0 - H_obs_i / H_exp_i

pca = pd.read_csv(PCA_TAB).set_index("sample").reindex(samples)
clu = pca["cluster"].astype(int).values

f_df = pd.DataFrame({
    "sample": samples,
    "n_called": nm.sum(axis=0),
    "H_obs": H_obs_i,
    "H_exp": H_exp_i,
    "F": F_ind,
    "cluster": clu + 1,
}).sort_values("F", ascending=False)
f_df.to_csv(TAB / "f_per_accession.csv", index=False)
print(f"[F-ind] median individual F = {np.nanmedian(F_ind):.3f}, "
      f"range [{F_ind.min():.3f}, {F_ind.max():.3f}]")


# ---------------------------------------------------------------------------
# Fig 56. Per-locus F_IS distribution
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.4), constrained_layout=True)
ax.hist(F_IS[ok], bins=60, color=WONG["blue"], edgecolor="white",
        linewidth=0.4)
ax.axvline(0, color="black", linewidth=0.6, linestyle=":",
           label="HWE expectation")
ax.axvline(global_fis, color=WONG["vermillion"], linewidth=1.5,
           label=f"global F$_{{IS}}$ = {global_fis:.3f}\n"
                 f"95 % CI [{lo:.3f}, {hi:.3f}]")
ax.set_xlabel(r"per-locus F$_{IS}$ = 1 - H$_{obs}$ / H$_{exp}$")
ax.set_ylabel("count of loci")
ax.set_title(f"Per-locus F$_{{IS}}$ across {ok.sum():,} AYB loci "
             f"(selfing-rate signature)")
ax.legend(loc="upper left")
ax.grid(True, axis="y")
save(fig, "fig56_fis_distribution")
print("[fig] fig56_fis_distribution")


# ---------------------------------------------------------------------------
# Fig 57. Per-accession F, sorted, coloured by cluster
# ---------------------------------------------------------------------------
sorted_f = f_df.copy().reset_index(drop=True)  # already sorted desc
fig, ax = plt.subplots(figsize=(12, 4.6), constrained_layout=True)
colors = [CLUSTER_PAL[c - 1] for c in sorted_f["cluster"]]
x = np.arange(len(sorted_f))
ax.bar(x, sorted_f["F"], color=colors, edgecolor="black", linewidth=0.25,
       zorder=3)

median_F = float(np.nanmedian(F_ind))
ax.axhline(median_F, color="black", linewidth=0.8, linestyle="--", zorder=2)

# Rank-based x-axis so 95 IDs do not collide -- per-accession identities
# are labelled only for the 12 most extreme samples via adjustText.
ax.set_xticks([0, len(sorted_f) // 4, len(sorted_f) // 2,
               3 * len(sorted_f) // 4, len(sorted_f) - 1])
ax.set_xticklabels(["1", str(len(sorted_f) // 4 + 1),
                    str(len(sorted_f) // 2 + 1),
                    str(3 * len(sorted_f) // 4 + 1),
                    str(len(sorted_f))], fontsize=9)
ax.set_xlabel("Accession rank (descending F)")
ax.set_ylabel(r"individual inbreeding F")
ax.set_title(f"Per-accession F (Yang 2010 method I) -- "
             f"{len(sorted_f)} lines, ranked by F descending")

# Top-6 highest-F and top-6 lowest-F accessions get adjustText labels
extreme_idx = np.concatenate([np.arange(6),
                              np.arange(len(sorted_f) - 6, len(sorted_f))])
label_top_n(ax, x[extreme_idx], sorted_f["F"].values[extreme_idx],
            sorted_f["sample"].values[extreme_idx],
            n=len(extreme_idx), fontsize=7.0, color="#222222")

publishable_axes(ax, grid="y")

import matplotlib.patches as mpatches
handles = [mpatches.Patch(color=CLUSTER_PAL[k],
                          label=f"Cluster {k+1}") for k in (0, 1)]
handles.append(plt.Line2D([0], [0], color="black", linewidth=0.8,
                          linestyle="--",
                          label=f"median F = {median_F:.3f}"))
ax.legend(handles=handles, loc="upper right", fontsize=8.5,
          framealpha=0.9, edgecolor="none")
save(fig, "fig57_f_per_accession")
print("[fig] fig57_f_per_accession")


# ---------------------------------------------------------------------------
# Fig 58. F_IS per-locus Manhattan (cowpea-anchored subset, for visual)
# ---------------------------------------------------------------------------
per_locus["chr"] = per_locus["chrom"].astype(str).str.extract(
    r"(Vu\d+)", expand=False)
anc = per_locus.dropna(subset=["chr", "pos", "F_IS"]).copy()
anc["pos"] = anc["pos"].astype(int)
anc = anc.sort_values(["chr", "pos"])
chrs = sorted(anc["chr"].unique(),
              key=lambda s: int(re.search(r"\d+", s).group()))
fig, ax = plt.subplots(figsize=(11, 3.6), constrained_layout=True)
offsets, mids, x_cursor = {}, [], 0
xs = []
for c in chrs:
    sub = anc[anc["chr"] == c]
    offsets[c] = x_cursor
    xs.extend((sub["pos"].values + x_cursor).tolist())
    mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
    x_cursor += sub["pos"].max() + 5_000_000
anc["x"] = xs
palette = [WONG["blue"], "#7c7c7c"]
for i, c in enumerate(chrs):
    sub = anc[anc["chr"] == c]
    ax.scatter(sub["x"], sub["F_IS"], s=14,
               color=palette[i % 2], alpha=0.85, edgecolor="none")
ax.axhline(global_fis, color=WONG["vermillion"], linewidth=1.0,
           label=f"global F$_{{IS}}$ = {global_fis:.3f}")
ax.axhline(0, color="black", linewidth=0.6, linestyle=":")
ax.set_xticks(mids); ax.set_xticklabels(chrs, rotation=0, fontsize=8)
ax.set_ylabel(r"F$_{IS}$")
ax.set_title(f"Per-locus F$_{{IS}}$ across cowpea-anchored {len(anc)} SNPs")
ax.set_ylim(-0.1, 1.05)
ax.legend(loc="lower right")
ax.grid(True, axis="y")
save(fig, "fig58_fis_manhattan")
print("[fig] fig58_fis_manhattan")

print(f"\nOutputs in: {OUT}")
print(f"  global F_IS = {global_fis:.4f}  [{lo:.4f}, {hi:.4f}]")
print(f"  median individual F = {np.nanmedian(F_ind):.3f}")
