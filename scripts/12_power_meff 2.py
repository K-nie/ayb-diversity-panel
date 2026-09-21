#!/usr/bin/env python3
"""
Genomic-prediction accuracy ceiling (Daetwyler 2008) and effective number of
independent tests (Li & Ji 2005), for the AYB panel.

A. Daetwyler 2008 accuracy ceiling
   For prediction with a balanced GBLUP, the upper bound on accuracy r =
   correlation(GEBV, true breeding value) at training size n_p, trait
   heritability h^2, and effective number of independent chromosome segments
   M_e is:
       r^2 = (n_p h^2) / (n_p h^2 + M_e)
   M_e itself follows the Goddard 2009 / Hayes-Visscher reparameterisation:
       M_e = 2 N_e L / log(4 N_e L)
   where L is genome length in Morgans. With Ne from script 11 (LD decay)
   and cowpea genome 620.4 Mb at 1.06 cM/Mb (Lonardi et al. 2019), L ~ 6.58
   Morgans. We also compute sample size needed to reach r = 0.5 for each
   trait.

B. Effective number of independent tests m_eff
   Li & Ji (2005) Heredity 95:221-227 eigenvalue method on the SNP
   correlation matrix. m_eff = sum_i I(lambda_i >= 1) + (lambda_i - floor(lambda_i))
   for the eigenvalues of the K x K correlation matrix. With n = 105 < m
   = 1672, the rank is bounded by n - 1, so m_eff <= 104.

Outputs (results/12_power_meff/)
--------------------------------
tables/
    daetwyler_ceiling.csv    r ceiling at h^2 in {0.05, 0.10, 0.20, 0.50, 0.80}
                             with n in {50, 105, 250, 500, 1000}
    sample_size_for_r0.5.csv n needed for r = 0.5 at each h^2
    meff_summary.csv         m_eff via Li & Ji, plus Bonferroni cutoffs
                             at m vs m_eff
figures/
    fig29_daetwyler_curves.png/.pdf   r vs n curves at each h^2
    fig30_meff_eigenspectrum.png/.pdf eigenvalue spectrum + m_eff cumulative

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
LD_SUMMARY = ROOT / "results" / "11_ld_decay" / "tables" / "ld_summary.csv"
OUT = ROOT / "results" / "12_power_meff"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
GENOME_MB = 620.4         # IT97K-499-35 v1.0 assembly length (Lonardi 2019)
CM_PER_MB = 1.06
H2_LIST = [0.05, 0.10, 0.20, 0.50, 0.80]
N_LIST = [50, 105, 250, 500, 1000]
# N_PANEL and N_MARKERS are computed from the filtered dosage matrix below
# (rather than hardcoded) so the ceiling reflects the actual post-QC panel.
N_PANEL = None
N_MARKERS = None


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Inputs from prior scripts
# ---------------------------------------------------------------------------
ld = pd.read_csv(LD_SUMMARY).iloc[0]
Ne = float(ld["Ne_estimate"])
print(f"[in] Ne (LD-derived) = {Ne:.1f}")

L_morgans = GENOME_MB * CM_PER_MB / 100.0   # convert cM -> Morgans
Me = 2.0 * Ne * L_morgans / np.log(4.0 * Ne * L_morgans)
print(f"[Me] genome = {GENOME_MB:.1f} Mb -> {L_morgans:.2f} Morgans")
print(f"[Me] M_e (Goddard 2009) = {Me:.1f}")


# ---------------------------------------------------------------------------
# A. Daetwyler 2008 ceiling
# ---------------------------------------------------------------------------
def r2_ceiling(n, h2, Me_):
    return (n * h2) / (n * h2 + Me_)


def n_for_r(r_target, h2, Me_):
    """invert r^2 = n h^2 / (n h^2 + M_e) -> n = M_e r^2 / (h^2 (1 - r^2))"""
    r2 = r_target ** 2
    return Me_ * r2 / (h2 * (1.0 - r2))


rows = []
for h2 in H2_LIST:
    for n in N_LIST:
        r2 = r2_ceiling(n, h2, Me)
        rows.append({"h2": h2, "n": n, "r2_ceiling": r2,
                     "r_ceiling": np.sqrt(r2)})
daet = pd.DataFrame(rows)
daet.to_csv(TAB / "daetwyler_ceiling.csv", index=False)


ss = pd.DataFrame([
    {"h2": h2,
     "n_for_r_0.3": n_for_r(0.30, h2, Me),
     "n_for_r_0.5": n_for_r(0.50, h2, Me),
     "n_for_r_0.7": n_for_r(0.70, h2, Me)}
    for h2 in H2_LIST
])
ss.to_csv(TAB / "sample_size_for_r0.5.csv", index=False)
print("[A] daetwyler ceiling sample-size table:")
print(ss.to_string(index=False))


# Figure 29 is drawn later, after N_PANEL / N_MARKERS are determined from
# the post-QC dosage matrix.


# ---------------------------------------------------------------------------
# B. Li & Ji 2005 effective number of tests
# ---------------------------------------------------------------------------
print("[B] computing m_eff (Li & Ji 2005)...")
# Rebuild filtered dosage (same cascade as script 03)
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

M = dosage.loc[keep_m, keep_s].fillna(
    dosage.loc[keep_m, keep_s].mean(axis=1)
).T.values  # samples x markers
m = M.shape[1]
n = M.shape[0]
# Use post-QC counts for the Daetwyler ceiling + Bonferroni summary so the
# table reflects the actual analysed panel, not a stale hardcoded value.
N_PANEL = n
N_MARKERS = m
print(f"[B] M shape = ({n}, {m})  (post-QC panel)")

# Correlation matrix is m x m; rank is min(n-1, m). For n=105 we have ~104
# non-trivial eigenvalues. Compute via the n x n Gram of standardised columns
# which is far cheaper and shares non-zero eigenvalues.
Xs = (M - M.mean(axis=0)) / M.std(axis=0, ddof=1)
Xs = np.nan_to_num(Xs, nan=0.0)
# eigenvalues of m x m correlation matrix equal eigenvalues of
# (1/(n-1)) Xs Xs.T (n x n) -- scaled so they sum to m
G = Xs @ Xs.T / (n - 1)        # n x n
eig = np.linalg.eigvalsh(G)
eig = np.sort(eig)[::-1]
# Eigenvalues of the m x m correlation matrix: the non-zero ones equal these
# n eigenvalues; the remaining (m - n) eigenvalues are zero. Li & Ji's
# m_eff is the sum over all m eigenvalues of indicator(lambda >= 1) +
# (lambda - floor(lambda)). Below 1 -> contribute fractional, >= 1 -> 1.
def li_ji_meff(eigs):
    return float(np.sum((eigs >= 1).astype(float)
                        + np.where(eigs < 1, eigs - np.floor(eigs), 0.0)))


# eigenvalues of the full m x m correlation matrix (concatenate zeros)
eig_full = np.concatenate([eig, np.zeros(max(0, m - len(eig)))])
m_eff = li_ji_meff(eig_full)
print(f"[B] m_eff (Li & Ji) = {m_eff:.1f}")

# 99.5 % variance threshold (alternative count used in simpleM)
cumvar = np.cumsum(eig) / eig.sum()
k995 = int(np.searchsorted(cumvar, 0.995) + 1)
print(f"[B] # PCs to explain 99.5 % variance = {k995}")

bonf_full = 0.05 / N_MARKERS
bonf_meff = 0.05 / m_eff
print(f"[B] Bonferroni at m={N_MARKERS}: alpha = {bonf_full:.2e}; "
      f"-log10 = {-np.log10(bonf_full):.2f}")
print(f"[B] Bonferroni at m_eff={m_eff:.1f}: alpha = {bonf_meff:.2e}; "
      f"-log10 = {-np.log10(bonf_meff):.2f}")

pd.DataFrame([{
    "method": "Li & Ji 2005 eigenvalue m_eff",
    "n_markers_raw": N_MARKERS,
    "m_eff": m_eff,
    "pcs_for_99.5_pct_var": k995,
    "bonferroni_raw": bonf_full,
    "neglog10_bonf_raw": -np.log10(bonf_full),
    "bonferroni_meff": bonf_meff,
    "neglog10_bonf_meff": -np.log10(bonf_meff),
}]).to_csv(TAB / "meff_summary.csv", index=False)


# Figure 30: eigenvalue spectrum + cumulative variance
fig, axes = plt.subplots(1, 2, figsize=(10, 4.0), constrained_layout=True)
ax = axes[0]
ax.bar(np.arange(1, len(eig) + 1), eig, color=WONG["blue"],
       edgecolor="none")
ax.axhline(1.0, color=WONG["vermillion"], linestyle="--", linewidth=1.0,
           label="lambda = 1")
ax.set_xlabel("eigenvalue rank")
ax.set_ylabel("eigenvalue")
ax.set_yscale("log")
ax.set_title(f"Eigenvalue spectrum of {N_MARKERS}-SNP correlation matrix")
ax.legend()
ax.grid(True, which="both", axis="y")

ax = axes[1]
ax.plot(np.arange(1, len(cumvar) + 1), cumvar, color=WONG["blue"],
        linewidth=1.4)
ax.axhline(0.995, color="black", linestyle=":", linewidth=0.8,
           label="99.5 %")
ax.axvline(k995, color=WONG["vermillion"], linestyle="--", linewidth=1.0,
           label=f"k = {k995}")
ax.set_xlabel("number of eigenvalues retained")
ax.set_ylabel("cumulative variance fraction")
ax.set_title(f"m$_{{eff}}$ (Li & Ji) = {m_eff:.1f}")
ax.legend()
ax.grid(True)
save(fig, "fig30_meff_eigenspectrum")
print("[fig] fig30_meff_eigenspectrum")


# ---------------------------------------------------------------------------
# Figure 29 (deferred): Daetwyler r vs n curves; vline at the actual post-QC n
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
n_grid = np.linspace(20, 2000, 600)
palette = [WONG["blue"], WONG["green"], WONG["orange"],
           WONG["vermillion"], WONG["purple"]]
for h2, c in zip(H2_LIST, palette):
    r = np.sqrt(r2_ceiling(n_grid, h2, Me))
    ax.plot(n_grid, r, linewidth=1.6, color=c, label=f"h$^2$ = {h2}")
ax.axvline(N_PANEL, color="black", linestyle="--", linewidth=0.8,
           label=f"current n = {N_PANEL}")
ax.axhline(0.5, color="grey", linestyle=":", linewidth=0.8,
           label="r = 0.5 target")
ax.set_xlabel("training set size n")
ax.set_ylabel("prediction accuracy r (Daetwyler 2008 ceiling)")
ax.set_title(f"Genomic-prediction ceiling for AYB at M$_e$ = {Me:.0f}")
ax.legend()
ax.grid(True)
save(fig, "fig29_daetwyler_curves")
print("[fig] fig29_daetwyler_curves")


print(f"\nOutputs in: {OUT}")
print(f"  M_e = {Me:.1f}, m_eff = {m_eff:.1f}")
print(f"  At n = {N_PANEL}: r^2 ceiling per h^2 =")
for h2 in H2_LIST:
    print(f"    h^2 = {h2}: r = {np.sqrt(r2_ceiling(N_PANEL, h2, Me)):.3f}")
