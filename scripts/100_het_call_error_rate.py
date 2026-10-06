#!/usr/bin/env python3
"""100 — Empirical heterozygous-call error rate (Reviewer B.8).

The reviewer asks us to estimate an empirical het-call error rate and discuss
its effect on the main-cluster F_ROH, given that a zero-het ROH re-call lowers
F_ROH (§3.6; results/90). Two independent, internal estimators are available,
both resting on the fact that an autozygous genome should carry no true
heterozygotes, so any het call in such a genome is an error (or, at most,
rare residual true heterozygosity — an upper bound on the error rate):

  (A) the near-clonal trio (TSs151B / TSs358 / TSs361), whose members are
      genome-wide autozygous (F_ROH ~ 1.0). At any marker the trio's true
      genotype is homozygous, so a het call in a member is a miscall.
  (B) het genotypes that fall inside called ROH segments (autozygous tracts),
      pooled across the minor cluster (n = 11), where autozygosity is strongest.

We then relate the estimate to the observed vs zero-het F_ROH bracket
(results/90) to judge whether the 1-het-per-window tolerance is warranted.

Inputs (all read from disk):
  - AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv
  - refs/ayb_genome/ayb_marker_anchoring.csv
  - results/34_roh/tables/roh_calls.csv
  - results/90_roh_sensitivity_null/tables/froh_per_sample_obs_vs_zerohet.csv
Outputs:
  - results/100_het_call_error_rate/tables/het_error_summary.csv
  - results/100_het_call_error_rate/tables/trio_pair_hetvhom.csv
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
ROH = ROOT / "results" / "34_roh" / "tables" / "roh_calls.csv"
CLUST = ROOT / "results" / "90_roh_sensitivity_null" / "tables" / "froh_per_sample_obs_vs_zerohet.csv"
OUT = ROOT / "results" / "100_het_call_error_rate"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
TRIO = ["TSs151B", "TSs358", "TSs361"]
WINDOW_SNPS = 20  # matches 34_roh

# --- build post-QC, AYB-anchored dosage (0/1/2, NaN preserved) ----------------
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

anc = pd.read_csv(ANCHOR).dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc["snp_pos_ayb"] = anc["snp_pos_ayb"].astype(int)
dosage = dosage.loc[dosage.index.intersection(anc["rs"])]
anc = anc[anc["rs"].isin(dosage.index)].set_index("rs")
print(f"[load] anchored post-QC dosage = {dosage.shape[0]} markers x {dosage.shape[1]} samples")

# --- (A) trio estimator -------------------------------------------------------
trio = dosage[TRIO].copy()
# pairwise het-vs-hom discordances (reproduce results/98)
pair_rows = []
import itertools
for a, b in itertools.combinations(TRIO, 2):
    both = trio[[a, b]].dropna()
    ident = (both[a] == both[b]).sum()
    disc = both[both[a] != both[b]]
    hvh = (((disc[a] == 1) & (disc[b] != 1)) | ((disc[b] == 1) & (disc[a] != 1))).sum()
    opp = (((disc[a] == 0) & (disc[b] == 2)) | ((disc[a] == 2) & (disc[b] == 0))).sum()
    pair_rows.append(dict(pair=f"{a}-{b}", n_called=len(both), n_identical=int(ident),
                          n_discordant=len(disc), n_het_vs_hom=int(hvh), n_opp_hom=int(opp)))
pair_df = pd.DataFrame(pair_rows)
pair_df.to_csv(TAB / "trio_pair_hetvhom.csv", index=False)
print(pair_df.to_string(index=False))

# per-genotype het-miscall rate: at markers where >=2 trio members are called
# and the majority non-het call is homozygous, any het (==1) is a miscall.
n_geno = 0          # trio genotypes examined at truly-homozygous markers
n_het_err = 0       # het calls among them
for rs, row in trio.iterrows():
    vals = row.dropna().values
    if len(vals) < 2:
        continue
    homs = vals[vals != 1]
    if len(homs) < 2:
        continue  # need >=2 homozygous calls to declare true state homozygous
    # majority homozygous value
    if len(np.unique(homs)) == 1 or (homs == np.median(homs)).sum() >= 2:
        n_geno += len(vals)
        n_het_err += int((vals == 1).sum())
trio_err = n_het_err / n_geno if n_geno else float("nan")
print(f"[trio] het-miscalls {n_het_err} / {n_geno} genotypes = {trio_err*100:.3f}%")

# --- (B) het-in-ROH estimator (minor cluster) --------------------------------
clust = pd.read_csv(CLUST)
small = set(clust.loc[clust["is_small_cluster"] == True, "sample"])  # noqa: E712
roh = pd.read_csv(ROH)
# map each marker to (chr_ayb, pos)
mk_chr = anc["chr_ayb"]
mk_pos = anc["snp_pos_ayb"]

def het_in_roh(sample_set, label):
    n_in = 0; n_het = 0
    sub_roh = roh[roh["sample"].isin(sample_set)]
    # precompute per-sample per-chr marker arrays
    for samp in sample_set:
        s_roh = sub_roh[sub_roh["sample"] == samp]
        if s_roh.empty:
            continue
        dcol = dosage[samp]
        for _, seg in s_roh.iterrows():
            in_seg = (mk_chr == seg["chr"]) & (mk_pos >= seg["start_bp"]) & (mk_pos <= seg["end_bp"])
            g = dcol[in_seg.index[in_seg]].dropna()
            n_in += len(g)
            n_het += int((g == 1).sum())
    rate = n_het / n_in if n_in else float("nan")
    print(f"[het-in-ROH:{label}] {n_het} / {n_in} = {rate*100:.3f}%")
    return n_het, n_in, rate

minor_het, minor_tot, minor_rate = het_in_roh(small & set(dosage.columns), "minor")
all_het, all_tot, all_rate = het_in_roh(set(dosage.columns), "panel")

# --- effect on F_ROH: expected false-het load per 20-SNP window ---------------
froh = pd.read_csv(CLUST)
main_obs = froh.loc[~froh["is_small_cluster"], "F_ROH_obs"].mean()
main_zero = froh.loc[~froh["is_small_cluster"], "F_ROH_zerohet"].mean()
exp_false_het_per_win = WINDOW_SNPS * trio_err
print(f"[froh] main obs {main_obs:.4f} vs zero-het {main_zero:.4f}")
print(f"[froh] expected false hets per {WINDOW_SNPS}-SNP window at trio rate = {exp_false_het_per_win:.3f}")

summary = pd.DataFrame([
    dict(estimator="trio_per_genotype", n_het=n_het_err, n_total=n_geno, rate=trio_err),
    dict(estimator="het_in_ROH_minor_cluster", n_het=minor_het, n_total=minor_tot, rate=minor_rate),
    dict(estimator="het_in_ROH_panel", n_het=all_het, n_total=all_tot, rate=all_rate),
    dict(estimator="main_cluster_FROH_obs_1het", n_het=np.nan, n_total=np.nan, rate=main_obs),
    dict(estimator="main_cluster_FROH_zerohet", n_het=np.nan, n_total=np.nan, rate=main_zero),
    dict(estimator="exp_false_hets_per_20SNP_window", n_het=np.nan, n_total=np.nan, rate=exp_false_het_per_win),
])
summary.to_csv(TAB / "het_error_summary.csv", index=False)
print("\n", summary.to_string(index=False))
print("\n[done] wrote", TAB / "het_error_summary.csv")
