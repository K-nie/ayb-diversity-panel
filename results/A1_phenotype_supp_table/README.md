# Supplementary Table S14 — per-accession phenotype matrix (A1)

## What this does and why

Assembles the per-accession replicate-mean phenotype record that the A1
manuscript (*Sphenostylis stenocarpa* / African yam bean DArTseq panel) refers to
in §2.1. The reference previously read "(Data Sheet 1)" — an unlabelled pointer
with no deposited file behind it. To close reviewer flag F1 ("Supplementary Table
S# everywhere") the matrix is now deposited as **Supplementary Table S14** and the
body pointer rewritten to "(Supplementary Table S14)".

## Inputs

- `Updated_Wet_Chemistry_Data.xlsx` (repo root) — authoritative wet-chemistry /
  seed-metric workbook. Sheets used:
  - `Means` — per-accession anti-nutritional-factor means already computed from
    three laboratory replicates (Tannin, Phenol, Flavonoid [sheet spells it
    "Flavinoid"], Antioxidant).
  - `Crude_Protein` — ~2 replicates per accession → mean.
  - `Moisture_Content` — replicates per accession (covariate) → mean (`MC`).
  - `Seed_Metrics` — ten replicates per accession → mean (Length, Width,
    Thickness, Weight).
- `manuscript/paper_A1/share/tables/supp/TableS10_accession_list_qc_cluster.csv` —
  the 105 genotyped accessions (95 Pass + 10 Fail QC). Used as the anchor so the
  deposited phenotype matrix lines up one-to-one with the panel the paper reports.

## Method

`scripts/A1_supp_table_s14_phenotype.py`:
1. Load the four phenotype blocks; group-mean the replicate sheets per accession.
2. **Left-join** all four blocks onto the 105-accession panel (TableS10). Any
   accession missing a block is carried as blanks, not dropped.
3. Round to 4 dp; sort by accession.

Reproduce:
```
/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python \
  scripts/A1_supp_table_s14_phenotype.py
```

## Outputs

- `results/A1_phenotype_supp_table/TableS14_phenotype_per_accession.csv`
- `manuscript/paper_A1/share/tables/supp/TableS14_phenotype_per_accession.csv`

Both copies md5-identical. 105 rows × 11 columns: Accession, Tannin, Phenol,
Flavonoid, Antioxidant, Crude_Protein, Moisture_Content, Seed_Length, Seed_Width,
Seed_Thickness, Seed_Weight. All 105 carry seed metrics; two (TSs151B, TSs336)
carry blanks in the wet-chemistry columns — see the caveat below.

## Key decision — panel-anchored (105 rows)

Anchoring on TableS10 and left-joining keeps the deposited table in exact
correspondence with the genotyped panel the manuscript reports, and drops the one
genuinely non-panel wet-chemistry accession (`TSs708`, phenotyped but not
genotyped).

## RESOLVED — within-study label drift harmonised (2026-10-01)

Two accessions were filed under drifted labels between the wet-chemistry and
seed-metric measurement rounds:

| Panel / Seed_Metrics label | Wet-chemistry label (Means/Raw_ANF/Crude_Protein/Moisture) |
| --- | --- |
| `TSs151B` | `TSs151` |
| `TSs336`  | `TSs366` |

Confirmed the same-accession equivalence against the workbook (decisive, not a
guess):
- The drifted and canonical labels **never co-occur** in any sheet: every
  wet-chemistry sheet carries `TSs151` / `TSs366` and never `TSs151B` / `TSs336`;
  `Seed_Metrics` carries `TSs151B` / `TSs336` and never `TSs151` / `TSs366`.
- The `Seed_Metrics` accession set **equals the genotyped panel exactly** (zero
  differences), so Seed_Metrics is on the canonical panel labelling.
- The pairs are near-identical lexically: a dropped `B` (151 → 151B) and a
  `366` ↔ `336` digit swap.

If left unreconciled these labels split one accession into two half-populated
rows (wet-chemistry under one label with no genotype, seed metrics + genotype
under the other with no wet-chemistry) — which is exactly how the phenotype loader
consumed them before this fix.

**Fix applied in two places:**
1. **This deposit** — `scripts/A1_supp_table_s14_phenotype.py` maps the
   wet-chemistry labels to the panel label (`LABEL_FIX`) before the join, so
   `TSs151B` / `TSs336` now carry their full wet-chemistry. All 105 panel
   accessions carry both blocks; the only remaining blank is one moisture-content
   value for `TSs3` (the `Moisture_Content` sheet covers 104/105 accessions — a
   genuine gap, not a labelling issue).
2. **The source loader** — `scripts/_pheno.py` gains a `SAMPLE_ALIASES`
   harmonisation (`TSs151→TSs151B`, `TSs366→TSs336`) applied in `load_replicates`
   and `load_means`, so every downstream analysis across A1 and A2 joins
   phenotype to genotype correctly from now on.

**Impact on existing A1 numbers:** none at the stated-claim level. The only A1
manuscript number touching these traits is the broad-sense heritability
"H² ≈ 0.99" (¶207); re-fit after harmonisation it is Tannin 0.993, Phenol 0.995,
Flavonoid 0.989, Antioxidant 0.992 — still ≈ 0.99. The TSs358/TSs361 phenotype
divergence example (¶207) uses clean-label accessions and is unaffected.

**Downstream still to propagate (NOT rerun here):** any already-rendered output
built from `_pheno.py` — notably the A1 phenotype PCA (fig05, supplementary) and
the entire A2 trait-genetics pipeline (BLUPs, GWAS, GBLUP, ML prediction,
heritability) — now recovers two accessions' wet-chemistry into the
genotype × phenotype join and should be regenerated to pick up the correction.
A2's headline numbers may shift slightly; left for Benjamin to schedule.

## Tool versions

- Python 3.9.23 (`/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python`)
- pandas 2.3.1, openpyxl 3.1.5
