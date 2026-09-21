#!/usr/bin/env python3
"""
AMOVA on AYB DArTseq dosage data — Excoffier, Smouse & Quattro (1992).

Variance partition
------------------
SS_T = (1/N) sum_{i<j} d^2_ij                (total)
SS_W = sum_g (1/n_g) sum_{i<j in g} d^2_ij   (within groups)
SS_A = SS_T - SS_W                           (among groups)

df_A = G - 1; df_W = N - G
MS_A = SS_A / df_A; MS_W = SS_W / df_W
n0 = (N - sum_g n_g^2 / N) / (G - 1)         (avg group size correction)
sigma_W^2 = MS_W
sigma_A^2 = (MS_A - MS_W) / n0
Phi_ST = sigma_A^2 / (sigma_A^2 + sigma_W^2)

Distance: squared Euclidean on QC-filtered 0/1/2 dosage matrix.

Groups: K-means k=2 clusters from PCA of the same dosage matrix.
Caveat (explicitly stated in paper): grouping is SNP-derived, so AMOVA
quantifies the structure already visible in PCA — it is not independent
evidence. Replace with provenance / geographic origin if Prof. Adewale
supplies that metadata.

Significance: 999 permutations of group labels; empirical
p = (1 + #(Phi_ST_perm >= Phi_ST_obs)) / (1 + B).

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_OUT = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "02_amova"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

sns.set_context("talk")
sns.set_style("whitegrid")

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_PERM = 999
RNG = np.random.default_rng(0)


# ---------------------------------------------------------------------------
# Load HapMap and rebuild filtered dosage (identical filter to script 01)
# ---------------------------------------------------------------------------
print("[load] reading HapMap...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META_COLS = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
             "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META_COLS]
calls = hm[sample_cols].astype(str)


def dosage_row(row, alleles):
    ref, alt = alleles.split("/")
    out = np.full(len(row), np.nan)
    arr = row.values
    out[arr == ref + ref] = 0.0
    out[(arr == ref + alt) | (arr == alt + ref)] = 1.0
    out[arr == alt + alt] = 2.0
    return out


print("[load] converting calls to dosage...")
dosage = np.vstack([dosage_row(calls.iloc[i], hm["alleles"].iloc[i])
                    for i in range(len(hm))])
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)

call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]

keep_mark = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF))
keep_samp = (call_rate_s >= SAMP_CR)

X = dosage.loc[keep_mark.values, keep_samp.values].T  # samples x markers
X = X.fillna(X.mean(axis=0)).values
samples = dosage.columns[keep_samp.values].tolist()
print(f"[load] X shape = {X.shape}  (samples x markers)")


# ---------------------------------------------------------------------------
# Group labels from PCA cluster (k = 2 from script 01)
# ---------------------------------------------------------------------------
pca_tab = pd.read_csv(PCA_OUT)
groups = pca_tab.set_index("sample").loc[samples, "cluster"].astype(int).values
G = len(np.unique(groups))
N = len(samples)
n_g = np.array([np.sum(groups == g) for g in range(G)])
print(f"[load] N = {N}, G = {G}, group sizes = {n_g.tolist()}")


# ---------------------------------------------------------------------------
# Pairwise squared Euclidean distance (samples x samples)
# ---------------------------------------------------------------------------
print("[dist] computing pairwise squared distances...")
# D^2_ij = ||x_i - x_j||^2 = ||x_i||^2 + ||x_j||^2 - 2 x_i . x_j
sq = (X ** 2).sum(axis=1)
D2 = sq[:, None] + sq[None, :] - 2.0 * (X @ X.T)
D2 = np.maximum(D2, 0.0)  # clean numerical negatives
np.fill_diagonal(D2, 0.0)


def amova_phi_st(D2: np.ndarray, groups: np.ndarray) -> tuple[float, dict]:
    """Compute Phi_ST and the AMOVA partition for given group labels."""
    N_ = D2.shape[0]
    G_ = len(np.unique(groups))
    n_g_ = np.array([np.sum(groups == g) for g in range(G_)])
    # SS_T = (1/N) * sum_{i<j} D^2_ij
    triu = np.triu_indices(N_, k=1)
    ss_t = D2[triu].sum() / N_
    # SS_W: within-group sums
    ss_w = 0.0
    for g in range(G_):
        idx = np.where(groups == g)[0]
        if len(idx) < 2:
            continue
        sub = D2[np.ix_(idx, idx)]
        triu_g = np.triu_indices(len(idx), k=1)
        ss_w += sub[triu_g].sum() / len(idx)
    ss_a = ss_t - ss_w
    df_a = G_ - 1
    df_w = N_ - G_
    ms_a = ss_a / df_a
    ms_w = ss_w / df_w
    n0 = (N_ - (n_g_ ** 2).sum() / N_) / (G_ - 1)
    sigma_w = ms_w
    sigma_a = max(0.0, (ms_a - ms_w) / n0)
    denom = sigma_a + sigma_w
    phi_st = sigma_a / denom if denom > 0 else 0.0
    return phi_st, {
        "N": N_, "G": G_, "group_sizes": n_g_.tolist(),
        "SS_among": ss_a, "SS_within": ss_w, "SS_total": ss_t,
        "df_among": df_a, "df_within": df_w,
        "MS_among": ms_a, "MS_within": ms_w,
        "sigma2_among": sigma_a, "sigma2_within": sigma_w,
        "pct_var_among": 100 * sigma_a / denom if denom > 0 else 0.0,
        "pct_var_within": 100 * sigma_w / denom if denom > 0 else 0.0,
        "n0": n0,
    }


print("[amova] observed Phi_ST...")
phi_obs, parts = amova_phi_st(D2, groups)
print(f"[amova] Phi_ST_obs = {phi_obs:.4f}")
for k, v in parts.items():
    print(f"  {k}: {v}")


# ---------------------------------------------------------------------------
# Permutation test
# Track Phi_ST and the pseudo-F ratio (MS_among / MS_within) under each
# permuted label so the AMOVA table can carry a permutation p-value for both.
# For a two-group AMOVA the two p-values are essentially identical because
# Phi_ST and pseudo-F are monotone transforms of each other, but reviewers
# expect to see the F-row p-value alongside Phi_ST.
# ---------------------------------------------------------------------------
print(f"[perm] running {N_PERM} permutations...")
phi_null = np.empty(N_PERM)
F_null = np.empty(N_PERM)
g_arr = groups.copy()
for b in range(N_PERM):
    RNG.shuffle(g_arr)
    phi_b, parts_b = amova_phi_st(D2, g_arr)
    phi_null[b] = phi_b
    F_null[b] = (parts_b["MS_among"] / parts_b["MS_within"]
                 if parts_b["MS_within"] > 0 else np.nan)
F_obs = parts["MS_among"] / parts["MS_within"] if parts["MS_within"] > 0 else np.nan
p_phi = (1 + np.sum(phi_null >= phi_obs)) / (1 + N_PERM)
p_F = (1 + np.sum(F_null >= F_obs)) / (1 + N_PERM)
print(f"[perm] Phi_ST p-value = {p_phi:.4f}; pseudo-F p-value = {p_F:.4f}")
print(f"[perm] pseudo-F (MS_A / MS_W) observed = {F_obs:.4f}")
p_value = p_phi   # primary AMOVA p-value (Phi_ST permutation, Excoffier 1992)


# ---------------------------------------------------------------------------
# Save tables
# ---------------------------------------------------------------------------
# Classical AMOVA summary table — df / SS / MS / variance / pct / F / p
# F-ratio (pseudo-F = MS_among / MS_within) and its permutation p-value are
# placed on the Among-groups row only, since that is the row being tested.
amova_tab = pd.DataFrame({
    "source":             ["Among groups", "Within groups", "Total"],
    "df":                 [parts["df_among"], parts["df_within"],
                           parts["df_among"] + parts["df_within"]],
    "SS":                 [parts["SS_among"], parts["SS_within"],
                           parts["SS_total"]],
    "MS":                 [parts["MS_among"], parts["MS_within"], np.nan],
    "variance_component": [parts["sigma2_among"], parts["sigma2_within"],
                           parts["sigma2_among"] + parts["sigma2_within"]],
    "pct_variance":       [parts["pct_var_among"], parts["pct_var_within"],
                           100.0],
    "F_ratio":            [F_obs, np.nan, np.nan],
    "p_perm_F":           [p_F, np.nan, np.nan],
    "Phi_ST":             [phi_obs, np.nan, np.nan],
    "p_perm_Phi_ST":      [p_phi, np.nan, np.nan],
})
# Round numerics for readability before write — keep full precision in
# the in-memory frame, the CSV always carries full floats.
amova_tab.to_csv(TAB / "amova_table.csv", index=False)
print("\n[table] AMOVA classical table:")
print(amova_tab.to_string(index=False))

with (TAB / "amova_summary.txt").open("w") as f:
    f.write("AMOVA on AYB DArTseq dosage matrix\n")
    f.write("==================================\n")
    f.write(f"N samples: {parts['N']}\n")
    f.write(f"Groups: {parts['G']}\n")
    f.write(f"Group sizes: {parts['group_sizes']}\n")
    f.write(f"Filtered markers: {X.shape[1]}\n")
    f.write(f"Group source: K-means k=2 from PCA (script 01)\n\n")
    f.write(f"Phi_ST (observed): {phi_obs:.4f}\n")
    f.write(f"% variance among groups: {parts['pct_var_among']:.2f}\n")
    f.write(f"% variance within groups: {parts['pct_var_within']:.2f}\n\n")
    f.write(f"Pseudo-F (MS_among / MS_within): {F_obs:.4f}\n\n")
    f.write(f"Permutations: {N_PERM}\n")
    f.write(f"Null mean Phi_ST: {phi_null.mean():.4f}\n")
    f.write(f"Null max  Phi_ST: {phi_null.max():.4f}\n")
    f.write(f"Empirical p-value (Phi_ST):    {p_phi:.4f}\n")
    f.write(f"Empirical p-value (pseudo-F):  {p_F:.4f}\n\n")
    f.write("Caveat: groups are derived from PCA of the same dosage matrix, so\n")
    f.write("AMOVA here quantifies the structure visible in PCA rather than\n")
    f.write("offering independent evidence. Re-run with provenance / geographic\n")
    f.write("metadata when supplied by Prof. Adewale.\n")

pd.DataFrame({"phi_st_null": phi_null}).to_csv(TAB / "amova_null_distribution.csv",
                                                index=False)


# ---------------------------------------------------------------------------
# Figures
# Pie chart removed — variance components are already in the printed AMOVA
# table and the permutation null histogram below tells the inferential story.
# ---------------------------------------------------------------------------
print("[fig] permutation null histogram...")
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(phi_null, bins=40, color="#7c7c7c", edgecolor="white")
ax.axvline(phi_obs, color="#c14a4a", lw=2,
           label=f"Observed $\\Phi_{{ST}}$ = {phi_obs:.3f}")
ax.set_xlabel(r"$\Phi_{ST}$ under permuted labels")
ax.set_ylabel("Count")
ax.set_title(f"Permutation test (B = {N_PERM}); p = {p_value:.3f}")
ax.legend()
fig.tight_layout()
fig.savefig(FIG / "amova_permutation_null.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "amova_permutation_null.pdf", bbox_inches="tight")
plt.close(fig)

print("\nWrote:")
print(f"  {TAB/'amova_table.csv'}")
print(f"  {TAB/'amova_summary.txt'}")
print(f"  {TAB/'amova_null_distribution.csv'}")
print(f"  {FIG/'amova_permutation_null.png'}  + .pdf")
