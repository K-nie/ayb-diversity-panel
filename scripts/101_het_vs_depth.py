#!/usr/bin/env python3
"""
101 - Heterozygosity versus read depth and genotyping reproducibility (Reviewer 2, point B).

Reviewer 2 asked whether the panel's unexpectedly high observed heterozygosity in a
cleistogamous selfer is a technical artefact of low read depth or poor marker
reproducibility, and specifically requested that the DArT "reproducibility"
(RepAvg) and depth (AvgCountRef / AvgCountSnp) columns be used to test
heterozygosity against depth. This script does three things:

1. Marker level: correlate per-marker observed heterozygosity against mean read
   depth (AvgCountRef + AvgCountSnp) and against RepAvg over the post-QC marker
   set, so a reader can see whether low-quality markers carry the heterozygosity.
2. Outbred-tail lines: for each of the five most heterozygosity-rich main-cluster
   accessions (Wright F ~ -0.5 to -0.8; TSs325, TSs82, TSs63A, TSs104, TSs23),
   compare the depth, reproducibility and multi-mapping fraction of the markers at
   which that line is called HET against the markers at which it is called HOM. If
   its heterozygous calls do not sit on lower-depth / lower-reproducibility /
   multi-mapping markers, the excess heterozygosity is not a genotyping artefact.
3. Ties the depth/reproducibility verdict to the paralogy classification already
   produced in results/89 (unique / multi-mapping / unplaced).

Source of truth for depth + reproducibility is the DArT singlerow report, whose
native coding is 0 = homozygous-reference, 1 = homozygous-SNP, 2 = heterozygous
(verified here by matching each code's frequency to FreqHomRef / FreqHomSnp /
FreqHets at corr = 1.000). QC matches the main pipeline (script 01): marker
call-rate >= 0.90, MAF >= 0.05, sample call-rate >= 0.90 -> 1,625 markers x 95
lines (asserted against results/01_qc_pca_power/tables/marker_qc.csv).

Author: Benjamin Narh-Madey
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
SINGLEROW = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_singlerow_2.csv"
MARKER_QC = ROOT / "results" / "01_qc_pca_power" / "tables" / "marker_qc.csv"
PLACEMENT = ROOT / "results" / "89_blast_uniqueness_reanchor" / "tables" / "marker_placement_classification.csv"
OUT = ROOT / "results" / "101_het_vs_depth"
(OUT / "tables").mkdir(parents=True, exist_ok=True)

OUTBRED = ["TSs325", "TSs82", "TSs63A", "TSs104", "TSs23"]

# ---- load DArT singlerow report (depth + reproducibility live here) ----------
df = pd.read_csv(SINGLEROW, skiprows=6, low_memory=False)
samples = [c for c in df.columns if c.startswith("TSs")]
geno = df[samples].apply(pd.to_numeric, errors="coerce")  # 0=HomRef 1=HomSnp 2=Het

# verify the native coding before trusting it
for code, col in [(0, "FreqHomRef"), (1, "FreqHomSnp"), (2, "FreqHets")]:
    fr = (geno == code).sum(axis=1) / geno.notna().sum(axis=1)
    assert np.corrcoef(fr.fillna(0), df[col].fillna(0))[0, 1] > 0.999, f"coding check failed for {code}"

meta = df[["AlleleID", "AvgCountRef", "AvgCountSnp", "RepAvg"]].copy()
meta["depth"] = meta["AvgCountRef"] + meta["AvgCountSnp"]
meta = meta.rename(columns={"AlleleID": "rs"}).set_index("rs")

# ---- QC identical to script 01 -> 1,625 markers x 95 lines --------------------
# standard alt-dosage: HomRef(0)->0, Het(2)->1, HomSnp(1)->2
dos = geno.replace({0: 0, 2: 1, 1: 2})
dos.index = df["AlleleID"].values

sample_cr = dos.notna().mean(axis=0)
keep_samples = sample_cr[sample_cr >= 0.90].index.tolist()
dos = dos[keep_samples]

marker_cr = dos.notna().mean(axis=1)
altfreq = dos.mean(axis=1) / 2.0
maf = np.minimum(altfreq, 1 - altfreq)
keep_markers = (marker_cr >= 0.90) & (maf >= 0.05)
dos = dos[keep_markers]
het_code = geno.replace({0: 0, 1: 0, 2: 1})  # 1 if heterozygous
het_code.index = df["AlleleID"].values
het_code = het_code.loc[dos.index, keep_samples]

print(f"[101] post-QC: {dos.shape[0]} markers x {dos.shape[1]} samples")

# validate against the canonical QC table
qc = pd.read_csv(MARKER_QC).set_index("rs")
shared = dos.index.intersection(qc.index)
assert dos.shape[0] == 1625, f"expected 1625 markers, got {dos.shape[0]}"
assert dos.shape[1] == 95, f"expected 95 samples, got {dos.shape[1]}"
assert len(shared) == 1625, f"marker set disagrees with marker_qc ({len(shared)} shared)"

# canonical per-marker observed het (DArT FreqHets over all 105 samples, as in
# marker_qc / the manuscript); cross-checked below against a recomputation.
het_105 = (geno == 2).sum(axis=1) / geno.notna().sum(axis=1)
het_105.index = df["AlleleID"].values
assert (het_105.loc[dos.index] - qc.loc[dos.index, "het"]).abs().max() < 1e-6, \
    "recomputed 105-sample het disagrees with marker_qc"
m = pd.DataFrame({"het": qc.loc[dos.index, "het"]})
m = m.join(meta[["depth", "RepAvg", "AvgCountRef", "AvgCountSnp"]], how="left")

# placement class from script 89
pl = pd.read_csv(PLACEMENT).set_index("rs")["placement"]
m["placement"] = pl.reindex(m.index)
m["multimap"] = m["placement"].eq("multi")

# ---- (1) marker-level het vs depth / reproducibility -------------------------
rows = []
for x in ["depth", "RepAvg", "AvgCountRef", "AvgCountSnp"]:
    rho, p = spearmanr(m["het"], m[x], nan_policy="omit")
    rows.append({"test": f"het ~ {x} (Spearman, all post-QC markers)",
                 "rho": round(rho, 4), "p": f"{p:.3g}", "n": int(m[[x, 'het']].dropna().shape[0])})
corr_tbl = pd.DataFrame(rows)
print("\n[101] marker-level correlations:")
print(corr_tbl.to_string(index=False))

# ---- (2) per outbred-tail line: het vs hom marker quality --------------------
line_rows = []
for s in OUTBRED + ["ALL95"]:
    if s == "ALL95":
        ishet = het_code.values.ravel() == 1
        base = pd.DataFrame({"het": het_code.values.ravel()},
                            index=np.repeat(dos.index.values, het_code.shape[1]))
        base = base.join(m[["depth", "RepAvg", "multimap"]], how="left")
        sub_het = base[base["het"] == 1]
        sub_hom = base[base["het"] == 0]
    else:
        col = het_code[s]
        sub_het = m.loc[col[col == 1].index]
        sub_hom = m.loc[col[col == 0].index]
    rec = {"line": s,
           "n_het": len(sub_het), "n_hom": len(sub_hom),
           "het_markers_depth_med": round(sub_het["depth"].median(), 2),
           "hom_markers_depth_med": round(sub_hom["depth"].median(), 2),
           "het_markers_RepAvg_med": round(sub_het["RepAvg"].median(), 4),
           "hom_markers_RepAvg_med": round(sub_hom["RepAvg"].median(), 4),
           "het_multimap_frac": round(sub_het["multimap"].mean(), 4),
           "hom_multimap_frac": round(sub_hom["multimap"].mean(), 4)}
    # one-sided test: are het-markers LOWER reproducibility than hom-markers?
    try:
        u, pu = mannwhitneyu(sub_het["RepAvg"].dropna(), sub_hom["RepAvg"].dropna(),
                             alternative="less")
        rec["RepAvg_het<hom_p"] = f"{pu:.3g}"
    except ValueError:
        rec["RepAvg_het<hom_p"] = "NA"
    line_rows.append(rec)
line_tbl = pd.DataFrame(line_rows)
print("\n[101] per-line het-marker vs hom-marker quality:")
print(line_tbl.to_string(index=False))

corr_tbl.to_csv(OUT / "tables" / "het_vs_depth_marker_correlations.csv", index=False)
line_tbl.to_csv(OUT / "tables" / "outbred_line_het_marker_quality.csv", index=False)
m.reset_index().rename(columns={"index": "rs"}).to_csv(
    OUT / "tables" / "per_marker_het_depth_rep_placement.csv", index=False)
print(f"\n[101] wrote tables to {OUT/'tables'}")
