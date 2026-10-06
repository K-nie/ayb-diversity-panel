#!/usr/bin/env python3
"""98 - Direct genotype concordance for the near-clonal trio (Reviewer Item 6 / duplicate stats)

The reviewer notes that the three "orthogonal" layers behind the trio call are not
mutually independent (GRM + ML phylogeny share the same 1,625 SNPs; SilicoDArT shares
the same DNA extractions and libraries), so the alternative hypothesis at stake when
calling accessions duplicates is a sample-level artefact - contamination, a plate swap,
or a mixed seed lot (the reviewer's TSs357 example). The statistic that discriminates a
true genetic duplicate from a sample artefact is the DIRECT genotype concordance:

  - a true clonal/duplicate pair -> near-100% identical calls, and the few discordances
    are het-vs-hom (residual heterozygosity / calling noise), NOT opposite homozygotes;
  - a mixed seed lot or a plate swap of a distinct genotype -> a block of opposite-
    homozygote (0<->2) differences.

Also reports, per the reviewer's "Duplicate detection criteria" ask:
  G_ii (=1+F under VanRaden method 1), G_ij, and inbreeding-normalized G_ij/sqrt(G_ii G_jj).

Inputs
------
- AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv  (same QC as script 01)
- results/13_grm_crosspairs/tables/grm_vanraden.csv                  (VanRaden method-1 GRM)

Outputs (results/98_trio_duplicate_concordance/)
------------------------------------------------
tables/trio_concordance.csv     per-pair n_called, %identical, n_discordant, het/hom & opp-hom split
tables/trio_relatedness.csv      G_ii, G_ij, normalized relatedness for trio + TSs357
tables/close_pedigree_concordance.csv  same concordance stats for the five 0.78-0.82 G_ij pairs
tables/close_pedigree_relatedness.csv  G_ii, G_ij, normalized relatedness for the five pairs

Author: Benjamin Narh-Madey
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from itertools import combinations

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
GRM = ROOT / "results" / "13_grm_crosspairs" / "tables" / "grm_vanraden.csv"
OUT = ROOT / "results" / "98_trio_duplicate_concordance"
TAB = OUT / "tables"; TAB.mkdir(parents=True, exist_ok=True)
SAMP_CR = MARK_CR = 0.90; MIN_MAF = 0.05
TRIO = ["TSs151B", "TSs358", "TSs361"]; REF357 = "TSs357"
# five "close-pedigree" pairs flagged in results/56_dedup_sensitivity (0.78-0.82 G_ij range)
CLOSE_PEDIGREE = [("TSs361","TSs357"), ("TSs156A","TSs363"), ("TSs357","TSs358"),
                  ("TSs151B","TSs357"), ("TSs60","TSs282")]

hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#","alleles","chrom","pos","strand","assembly#","center",
        "protLSID","assayLSID","panelLSID","QCcode"]
scols = [c for c in hm.columns if c not in META]
calls = hm[scols].astype(str)
def dosage_row(row, alleles):
    ref, alt = alleles.split("/")
    out = np.full(len(row), np.nan); arr = row.values
    out[arr == ref+ref] = 0.0; out[(arr == ref+alt) | (arr == alt+ref)] = 1.0; out[arr == alt+alt] = 2.0
    return out
dos = pd.DataFrame(np.vstack([dosage_row(calls.iloc[i], hm["alleles"].iloc[i]) for i in range(len(hm))]),
                   index=hm["rs#"].values, columns=scols)
# QC identical to script 01
ns = dos.shape[1]
cr_m = dos.notna().sum(axis=1)/ns
mean_dose = dos.mean(axis=1, skipna=True); maf = np.minimum(mean_dose/2, 1-mean_dose/2)
cr_s = dos.notna().sum(axis=0)/dos.shape[0]
keep_s = cr_s[cr_s >= SAMP_CR].index.tolist()
keep_m = maf[(cr_m >= MARK_CR) & (maf >= MIN_MAF)].index.tolist()
X = dos.loc[keep_m, keep_s]          # markers x samples, NaN = missing
print(f"[qc] {len(keep_m)} markers x {len(keep_s)} samples")

def concordance(a, b):
    va, vb = X[a].values, X[b].values
    both = ~np.isnan(va) & ~np.isnan(vb)
    va, vb = va[both], vb[both]
    n = both.sum()
    ident = (va == vb).sum()
    disc = n - ident
    # discordance classes
    het_hom = ((va == 1) ^ (vb == 1)).sum() - ((va == 1) & (vb == 1) & (va != vb)).sum()  # exactly one is het
    het_hom = (((va == 1) & (vb != 1)) | ((vb == 1) & (va != 1))).sum()
    opp_hom = (((va == 0) & (vb == 2)) | ((va == 2) & (vb == 0))).sum()
    return dict(pair=f"{a}-{b}", n_called=int(n), n_identical=int(ident),
                pct_identical=round(100*ident/n, 3), n_discordant=int(disc),
                n_disc_het_vs_hom=int(het_hom), n_disc_opposite_hom=int(opp_hom))

rows = [concordance(a, b) for a, b in combinations(TRIO, 2)]
rows += [concordance(REF357, t) for t in TRIO]        # heterozygous line vs trio (mixed-seed-lot control)
conc = pd.DataFrame(rows)
conc.to_csv(TAB / "trio_concordance.csv", index=False)
print(conc.to_string(index=False))

# relatedness from GRM
g = pd.read_csv(GRM, index_col=0)
rel_rows = []
allids = TRIO + [REF357]
for a, b in combinations(TRIO, 2):
    rel_rows.append(dict(pair=f"{a}-{b}", G_ij=round(g.loc[a,b],4),
                         G_ii=round(g.loc[a,a],4), G_jj=round(g.loc[b,b],4),
                         norm_rel=round(g.loc[a,b]/np.sqrt(g.loc[a,a]*g.loc[b,b]),4)))
for t in TRIO:
    rel_rows.append(dict(pair=f"{REF357}-{t}", G_ij=round(g.loc[REF357,t],4),
                         G_ii=round(g.loc[REF357,REF357],4), G_jj=round(g.loc[t,t],4),
                         norm_rel=round(g.loc[REF357,t]/np.sqrt(g.loc[REF357,REF357]*g.loc[t,t]),4)))
rel = pd.DataFrame(rel_rows)
rel.to_csv(TAB / "trio_relatedness.csv", index=False)
print("\nG_ii:", {t: round(g.loc[t,t],4) for t in allids})
print(rel.to_string(index=False))

# ---- five close-pedigree pairs (0.78-0.82 G_ij) -------------------------------
cp_conc = pd.DataFrame([concordance(a, b) for a, b in CLOSE_PEDIGREE])
cp_conc.to_csv(TAB / "close_pedigree_concordance.csv", index=False)
print("\n[close-pedigree] concordance")
print(cp_conc.to_string(index=False))

cp_rel_rows = []
for a, b in CLOSE_PEDIGREE:
    cp_rel_rows.append(dict(pair=f"{a}-{b}", G_ij=round(g.loc[a,b],4),
                            G_ii=round(g.loc[a,a],4), G_jj=round(g.loc[b,b],4),
                            norm_rel=round(g.loc[a,b]/np.sqrt(g.loc[a,a]*g.loc[b,b]),4)))
cp_rel = pd.DataFrame(cp_rel_rows)
cp_rel.to_csv(TAB / "close_pedigree_relatedness.csv", index=False)
print("\n[close-pedigree] relatedness")
print(cp_rel.to_string(index=False))

# ---- GRM off-diagonal mean (reviewer cites -0.011 as the unrelated expectation) ----
gv = g.values.astype(float)
off = ~np.eye(gv.shape[0], dtype=bool)
off_mean = float(gv[off].mean())
diag_mean = float(np.diag(gv).mean())
print(f"\n[grm] off-diagonal mean G_ij = {off_mean:.4f} (n={gv.shape[0]} accessions); "
      f"diagonal mean G_ii = {diag_mean:.4f}")
pd.DataFrame([dict(grm_offdiag_mean=round(off_mean,4), grm_diag_mean=round(diag_mean,4),
                   n_accessions=gv.shape[0])]).to_csv(TAB / "grm_offdiag_summary.csv", index=False)
