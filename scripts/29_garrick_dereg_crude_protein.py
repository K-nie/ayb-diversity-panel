#!/usr/bin/env python3
"""
Garrick de-regression sensitivity check for Crude_Protein.

Crude_Protein is the one replicated trait where the BLUP shrinkage is
non-trivial. With H^2 = 0.71 and n_rep = 2 for most genotypes, the
per-sample reliability lambda_i ≈ 0.83, so a shrunk BLUP used as the
response in a downstream GBLUP / GWAS step inherits a ~17 % reduction
in the genetic effect (Garrick, Taylor & Fernando 2009).

This sensitivity check de-regresses the Crude_Protein BLUPs back to the
original scale, then re-runs the full-panel leave-one-out (LOO) GBLUP
and the 5-fold CV GBLUP that scripts 18 / 15 use elsewhere, and reports
the change in predictive ability. If the change is below 0.02 in r, the
shrunk BLUPs are safe to keep in the main text with a footnote; if it
exceeds 0.02, the supplement should carry both numbers.

De-regression formula (Garrick et al. 2009):
    BLUP_dereg_i = mean + (BLUP_i - mean) / reliability_i
    weight_i = (1 - h^2) / ((c + (1 - reliability_i) / reliability_i) * h^2)
where the heritability h^2 here is the trait Holland H^2 (0.71 for
Crude_Protein) and c is the proportion of genetic variance not explained
by markers (Garrick default c = 0.5).

Outputs (results/29_garrick_dereg_crude_protein/)
-------------------------------------------------
tables/
    cp_blups_shrunk_vs_dereg.csv     per-sample BLUP, reliability, dereg-BLUP
    cp_dereg_sensitivity.csv         LOO r and CV r for shrunk vs dereg
figures/
    fig67_garrick_dereg_cp.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import KFold

from _plotstyle import apply, WONG
from _pheno import load_blups
apply()
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
REL_CSV = ROOT / "results" / "00_blup_pipeline" / "tables" / "reliability_per_sample.csv"
OUT = ROOT / "results" / "29_garrick_dereg_crude_protein"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
ALPHAS = [0.5, 1.0, 2.0, 5.0, 10.0, 25.0]
N_FOLDS = 5
N_CV_REPS = 50
H2_CP = 0.71              # Holland H^2 for Crude_Protein from BLUP pipeline
C_GENETIC = 0.5           # Garrick default proportion of unexplained genetic variance


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fisher_z_ci(r, n, alpha=0.05):
    if n < 4 or not np.isfinite(r) or abs(r) >= 1:
        return float("nan"), float("nan")
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - z_crit * se)), float(np.tanh(z + z_crit * se))


# ---------------------------------------------------------------------------
# Load
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
dose_filt = dosage.loc[keep_m, keep_s]
samples_geno = dose_filt.columns.tolist()


# Shrunk BLUPs (canonical) + per-sample reliability
print("[load] shrunk BLUPs + reliability per sample for Crude_Protein...")
blups = load_blups()
cp_shrunk = blups["Crude_Protein"].dropna()
rel = pd.read_csv(REL_CSV)
rel_cp = (rel[rel["trait"] == "Crude_Protein"]
          .set_index("sample")["reliability"])

shared = sorted(set(cp_shrunk.index) & set(rel_cp.index) & set(samples_geno))
print(f"[load] {len(shared)} samples with CP BLUP + reliability + genotype")

cp_s = cp_shrunk.loc[shared]
lam = rel_cp.loc[shared]

# Garrick de-regression. The mean is added back after dividing because
# what is shrunk is the deviation from the population mean, not the mean
# itself.
mu = float(cp_s.mean())
cp_dereg = mu + (cp_s - mu) / lam.clip(lower=0.05)
weight = ((1 - H2_CP)
          / ((C_GENETIC + (1 - lam) / lam.clip(lower=0.05)) * H2_CP))

per_sample = pd.DataFrame({
    "sample": shared,
    "CP_shrunk_BLUP": cp_s.values,
    "reliability": lam.values,
    "CP_dereg_BLUP": cp_dereg.values,
    "WLS_weight": weight.values,
})
per_sample.to_csv(TAB / "cp_blups_shrunk_vs_dereg.csv", index=False)
print(f"[blups] mean shrinkage factor = {lam.mean():.3f}  "
      f"(range {lam.min():.3f} - {lam.max():.3f})")
print(f"[blups] shrunk-vs-dereg Pearson r = {cp_s.corr(cp_dereg):.4f}")


# ---------------------------------------------------------------------------
# Build K on the shared sample set
# ---------------------------------------------------------------------------
dose_sub = dose_filt[shared].T.values.astype(float)
col_mean = np.nanmean(dose_sub, axis=0)
inds = np.where(np.isnan(dose_sub))
dose_sub[inds] = np.take(col_mean, inds[1])
maf_sub = np.minimum(dose_sub.mean(axis=0) / 2.0,
                     1.0 - dose_sub.mean(axis=0) / 2.0)
var_sub = dose_sub.var(axis=0)
keep_marker = (maf_sub >= MIN_MAF) & (var_sub > 1e-8)
M_sub = dose_sub[:, keep_marker]
p = M_sub.mean(axis=0) / 2.0
p = np.clip(p, 1e-6, 1 - 1e-6)
Z = M_sub - 2.0 * p[None, :]
K = (Z @ Z.T) / (2.0 * (p * (1.0 - p)).sum())
print(f"[K] {K.shape[0]} samples x {M_sub.shape[1]} markers")


# ---------------------------------------------------------------------------
# LOO GBLUP for shrunk and de-regressed BLUPs
# ---------------------------------------------------------------------------
def fit_loo(y, K, alphas=ALPHAS):
    n = len(y)
    best_alpha, best_score = None, -np.inf
    for alpha in alphas:
        preds = np.zeros(n)
        for i in range(n):
            tr = np.arange(n) != i
            kr = KernelRidge(kernel="precomputed", alpha=alpha)
            kr.fit(K[np.ix_(tr, tr)], y[tr])
            preds[i] = kr.predict(K[i:i+1, tr])[0]
        if np.std(preds) > 0:
            r = stats.pearsonr(preds, y)[0]
            if r > best_score:
                best_score = r; best_alpha = alpha
    return best_score, best_alpha


print("\n[LOO] shrunk-BLUP Crude_Protein GBLUP...")
r_shrunk_loo, a_shrunk = fit_loo(cp_s.values, K)
r_shrunk_lo, r_shrunk_hi = fisher_z_ci(r_shrunk_loo, len(cp_s))
print(f"  alpha = {a_shrunk:.2f}, LOO r = {r_shrunk_loo:+.3f}  "
      f"95 % CI [{r_shrunk_lo:+.3f}, {r_shrunk_hi:+.3f}]")

print("[LOO] de-regressed Crude_Protein GBLUP...")
r_dereg_loo, a_dereg = fit_loo(cp_dereg.values, K)
r_dereg_lo, r_dereg_hi = fisher_z_ci(r_dereg_loo, len(cp_dereg))
print(f"  alpha = {a_dereg:.2f}, LOO r = {r_dereg_loo:+.3f}  "
      f"95 % CI [{r_dereg_lo:+.3f}, {r_dereg_hi:+.3f}]")


# ---------------------------------------------------------------------------
# 5-fold CV x 50 reps for shrunk and de-regressed BLUPs
# ---------------------------------------------------------------------------
def cv_50(y, K, n_folds=5, n_reps=50, seed=42):
    rs = np.random.default_rng(seed)
    accs = []
    for r in range(n_reps):
        kf = KFold(n_splits=n_folds, shuffle=True,
                   random_state=int(rs.integers(0, 2**31 - 1)))
        preds = np.full(len(y), np.nan)
        for tr, te in kf.split(np.arange(len(y))):
            best_a, best_s = None, -np.inf
            for a in ALPHAS:
                kr = KernelRidge(kernel="precomputed", alpha=a)
                kr.fit(K[np.ix_(tr, tr)], y[tr])
                yhat = kr.predict(K[np.ix_(te, tr)])
                if np.std(yhat) > 0:
                    sc = stats.pearsonr(yhat, y[te])[0]
                    if sc > best_s:
                        best_s = sc; best_a = a
            kr = KernelRidge(kernel="precomputed", alpha=best_a or 5.0)
            kr.fit(K[np.ix_(tr, tr)], y[tr])
            preds[te] = kr.predict(K[np.ix_(te, tr)])
        ok = ~np.isnan(preds)
        if ok.sum() > 5 and np.std(preds[ok]) > 0:
            accs.append(stats.pearsonr(preds[ok], y[ok])[0])
    return np.array(accs)


print("\n[CV] 5-fold x 50-rep CV — shrunk BLUP...")
cv_shrunk = cv_50(cp_s.values, K)
r_shrunk_cv = float(cv_shrunk.mean())
r_shrunk_cv_lo, r_shrunk_cv_hi = fisher_z_ci(r_shrunk_cv, len(cp_s))
print(f"  mean r = {r_shrunk_cv:+.3f}  sd = {cv_shrunk.std():.3f}  "
      f"Fisher-z 95 % CI [{r_shrunk_cv_lo:+.3f}, {r_shrunk_cv_hi:+.3f}]")

print("[CV] 5-fold x 50-rep CV — de-regressed BLUP...")
cv_dereg = cv_50(cp_dereg.values, K)
r_dereg_cv = float(cv_dereg.mean())
r_dereg_cv_lo, r_dereg_cv_hi = fisher_z_ci(r_dereg_cv, len(cp_dereg))
print(f"  mean r = {r_dereg_cv:+.3f}  sd = {cv_dereg.std():.3f}  "
      f"Fisher-z 95 % CI [{r_dereg_cv_lo:+.3f}, {r_dereg_cv_hi:+.3f}]")


# ---------------------------------------------------------------------------
# Verdict + write the sensitivity table
# ---------------------------------------------------------------------------
delta_loo = r_dereg_loo - r_shrunk_loo
delta_cv = r_dereg_cv - r_shrunk_cv

sensitivity = pd.DataFrame([
    {"metric": "LOO r",
     "shrunk_r": r_shrunk_loo,
     "shrunk_ci_low": r_shrunk_lo, "shrunk_ci_high": r_shrunk_hi,
     "dereg_r": r_dereg_loo,
     "dereg_ci_low": r_dereg_lo, "dereg_ci_high": r_dereg_hi,
     "delta_r_dereg_minus_shrunk": delta_loo,
     "verdict": ("material" if abs(delta_loo) > 0.02 else "negligible")},
    {"metric": "5-fold x 50-rep CV mean r",
     "shrunk_r": r_shrunk_cv,
     "shrunk_ci_low": r_shrunk_cv_lo, "shrunk_ci_high": r_shrunk_cv_hi,
     "dereg_r": r_dereg_cv,
     "dereg_ci_low": r_dereg_cv_lo, "dereg_ci_high": r_dereg_cv_hi,
     "delta_r_dereg_minus_shrunk": delta_cv,
     "verdict": ("material" if abs(delta_cv) > 0.02 else "negligible")},
])
sensitivity.to_csv(TAB / "cp_dereg_sensitivity.csv", index=False)

print("\n[verdict] LOO  delta-r = {:+.3f}  -> {}".format(
    delta_loo, "MATERIAL" if abs(delta_loo) > 0.02 else "NEGLIGIBLE (< 0.02)"))
print("[verdict]   CV  delta-r = {:+.3f}  -> {}".format(
    delta_cv, "MATERIAL" if abs(delta_cv) > 0.02 else "NEGLIGIBLE (< 0.02)"))


# ---------------------------------------------------------------------------
# Fig 67. Scatter of shrunk vs de-regressed BLUPs + sensitivity bars
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)

ax = axes[0]
ax.scatter(per_sample["CP_shrunk_BLUP"], per_sample["CP_dereg_BLUP"],
            s=24, color=WONG["blue"], edgecolor="black", linewidth=0.4,
            alpha=0.8)
lo = min(per_sample["CP_shrunk_BLUP"].min(),
          per_sample["CP_dereg_BLUP"].min())
hi = max(per_sample["CP_shrunk_BLUP"].max(),
          per_sample["CP_dereg_BLUP"].max())
ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="y = x")
ax.set_xlabel("Shrunk BLUP (Crude_Protein, %)")
ax.set_ylabel("De-regressed BLUP (Garrick 2009, Crude_Protein, %)")
r_xy = per_sample["CP_shrunk_BLUP"].corr(per_sample["CP_dereg_BLUP"])
ax.set_title(f"(A) Shrunk vs de-regressed BLUPs\n"
              f"Pearson r = {r_xy:.4f}, n = {len(per_sample)}",
              fontsize=10)
ax.legend(loc="lower right")
ax.grid(True)

ax = axes[1]
x = np.arange(2)
ax.bar(x - 0.18, [r_shrunk_loo, r_shrunk_cv],
        yerr=[[r_shrunk_loo - r_shrunk_lo, r_shrunk_cv - r_shrunk_cv_lo],
               [r_shrunk_hi - r_shrunk_loo, r_shrunk_cv_hi - r_shrunk_cv]],
        width=0.34, color="#7c7c7c", edgecolor="black", capsize=4,
        label="Shrunk BLUP")
ax.bar(x + 0.18, [r_dereg_loo, r_dereg_cv],
        yerr=[[r_dereg_loo - r_dereg_lo, r_dereg_cv - r_dereg_cv_lo],
               [r_dereg_hi - r_dereg_loo, r_dereg_cv_hi - r_dereg_cv]],
        width=0.34, color=WONG["vermillion"], edgecolor="black", capsize=4,
        label="De-regressed BLUP")
ax.axhline(0, color="black", lw=0.5)
ax.set_xticks(x)
ax.set_xticklabels(["LOO r", "5-fold x 50-rep CV r"])
ax.set_ylabel("Crude_Protein predictive ability (Pearson r)")
ax.set_title(f"(B) Predictive ability — Crude_Protein only\n"
              f"$\\Delta r_{{LOO}}$ = {delta_loo:+.3f}, "
              f"$\\Delta r_{{CV}}$ = {delta_cv:+.3f}",
              fontsize=10)
ax.legend(loc="upper right")
ax.grid(True, axis="y")

save(fig, "fig67_garrick_dereg_cp")
print("[fig] fig67_garrick_dereg_cp")
print(f"\nOutputs in: {OUT}")
