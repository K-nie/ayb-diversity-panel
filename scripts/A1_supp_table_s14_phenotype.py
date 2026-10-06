#!/usr/bin/env python3
"""
Supplementary Table S14 -- per-accession phenotype matrix for paper A1.

Assembles the per-accession replicate-mean phenotype record that the A1
manuscript refers to in Section 2.1 (previously cited only as an unlabelled
"Data Sheet"). Source of truth is the project wet-chemistry / seed-metric
workbook `Updated_Wet_Chemistry_Data.xlsx`:

  - Means            per-accession anti-nutritional-factor means already
                     computed from three laboratory replicates
                     (Tannin, Phenol, Flavonoid, Antioxidant).
  - Crude_Protein    two replicates per accession -> mean.
  - Moisture_Content replicates per accession (co-variate) -> mean.
  - Seed_Metrics     ten replicates per accession -> mean
                     (Length, Width, Thickness, Weight/seed mass).

Per-accession values are replicate means, matching the text. TSs151B carries
seed metrics only (no wet-chemistry record), reproduced faithfully as blanks in
the ANF / protein / moisture columns.

Outputs
-------
results/A1_phenotype_supp_table/TableS14_phenotype_per_accession.csv
manuscript/paper_A1/share/tables/supp/TableS14_phenotype_per_accession.csv

Author: Benjamin Narh-Madey
"""

from pathlib import Path

import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
XLSX = ROOT / "Updated_Wet_Chemistry_Data.xlsx"
PANEL = (ROOT / "manuscript" / "paper_A1" / "share" / "tables" / "supp"
         / "TableS10_accession_list_qc_cluster.csv")
OUT_RES = ROOT / "results" / "A1_phenotype_supp_table"
OUT_SHARE = ROOT / "manuscript" / "paper_A1" / "share" / "tables" / "supp"
OUT_RES.mkdir(parents=True, exist_ok=True)
FNAME = "TableS14_phenotype_per_accession.csv"

# Within-study label drift: the wet-chemistry sheets (Means, Raw_ANF,
# Crude_Protein, Moisture_Content) file two accessions under TSs151 and TSs366,
# whereas the Seed_Metrics sheet and the genotyped panel (Supplementary Table
# S10) file the same two under TSs151B and TSs336. The labels never co-occur in
# any sheet, the Seed_Metrics accession set equals the panel exactly, and the
# pairs are near-identical lexically (dropped "B"; 366<->336 digit swap). They
# are the same accessions; harmonise the wet-chemistry labels to the panel label
# so their wet-chemistry joins instead of being orphaned.
LABEL_FIX = {"TSs151": "TSs151B", "TSs366": "TSs336"}

means = pd.read_excel(XLSX, sheet_name="Means").rename(
    columns={"Genotypes": "Accession", "Flavinoid": "Flavonoid"})
means["Accession"] = means["Accession"].replace(LABEL_FIX)

cp = (pd.read_excel(XLSX, sheet_name="Crude_Protein")
        .assign(Genotypes=lambda d: d["Genotypes"].replace(LABEL_FIX))
        .groupby("Genotypes", as_index=False)["Crude_Protein"].mean()
        .rename(columns={"Genotypes": "Accession"}))

mc = (pd.read_excel(XLSX, sheet_name="Moisture_Content")
        .assign(Genotypes=lambda d: d["Genotypes"].replace(LABEL_FIX))
        .groupby("Genotypes", as_index=False)["MC"].mean()
        .rename(columns={"Genotypes": "Accession", "MC": "Moisture_Content"}))

sm = (pd.read_excel(XLSX, sheet_name="Seed_Metrics")
        .groupby("Genotypes", as_index=False)[["Length", "Width", "Thickness", "Weight"]].mean()
        .rename(columns={"Genotypes": "Accession",
                         "Length": "Seed_Length", "Width": "Seed_Width",
                         "Thickness": "Seed_Thickness", "Weight": "Seed_Weight"}))

# Anchor on the 105 genotyped accessions (Supplementary Table S10), so the
# deposited phenotype matrix lines up one-to-one with the panel the paper
# reports. Left-join the four phenotype blocks so an accession missing any
# block is carried as blanks rather than dropped.
panel = (pd.read_csv(PANEL)[["accession"]]
           .rename(columns={"accession": "Accession"}))

tab = panel.merge(means, on="Accession", how="left") \
           .merge(cp, on="Accession", how="left") \
           .merge(mc, on="Accession", how="left") \
           .merge(sm, on="Accession", how="left")

# round for a clean deposited table without manufacturing precision
num = tab.columns.drop("Accession")
tab[num] = tab[num].round(4)
tab = tab.sort_values("Accession").reset_index(drop=True)

tab.to_csv(OUT_RES / FNAME, index=False)
tab.to_csv(OUT_SHARE / FNAME, index=False)

print(f"[S14] {len(tab)} accessions x {len(tab.columns)} columns")
print(f"[S14] columns: {list(tab.columns)}")
print(f"[S14] panel accessions with no wet-chemistry (ANF NaN): "
      f"{tab.loc[tab['Tannin'].isna(), 'Accession'].tolist()}")
print(f"[S14] panel accessions with no seed metrics: "
      f"{tab.loc[tab['Seed_Length'].isna(), 'Accession'].tolist()}")
