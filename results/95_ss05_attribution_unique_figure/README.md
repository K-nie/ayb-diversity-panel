# 95 — Ss05 outlier attribution on the de-duplicated unique-SNP set

**Script:** `scripts/95_ss05_attribution_unique_figure.py`

## What was done and why

Re-draws the Ss05 selection-outlier attribution figure (manuscript Fig 7C)
on **unique SNPs** instead of method × SNP rows, in response to reviewers
**R1-M7 / R2-C**.

The submitted Fig 7C plotted the PCAdapt and DAPC top-50 Ss05 outliers as
separate method rows (script 55). A SNP flagged by *both* scans therefore
appeared twice, and the "Cluster-2-driven" legend counted method × SNP rows
(6 PCAdapt + 7 DAPC = 13) against an untraceable denominator ("13 of 20").
De-duplicating to unique loci gives the corrected, defensible statement the
manuscript now reports: **of the 16 unique Ss05 outlier SNPs, 7 (44 %) are
Cluster-2-driven.** This matches the corrected Results text and the
permutation null in script 94.

## Method

- **Input:** `results/55_private_alleles_per_cluster/tables/Ss05_outlier_attribution.csv`
  (22 method × SNP rows spanning 16 unique SNPs).
- De-duplicate to unique SNPs with `drop_duplicates("rs")`. A SNP's
  attribution is a deterministic function of its MAF_C1 / MAF_C2, so it is
  identical across the methods that flagged it; the first row per `rs` is
  representative.
- **Attribution thresholds (unchanged from script 55):**
  - `Cluster-2-driven`  if MAF_C2 ≥ 0.40 **and** MAF_C1 ≤ 0.10
  - `Cluster-1-driven`  if MAF_C1 ≥ 0.40 **and** MAF_C2 ≤ 0.10
  - `panel-wide`        otherwise
- **Cluster convention (R2-C):** MAF_C1 = main cluster (n = 84, "Cluster 1",
  Wong blue); MAF_C2 = minor cluster (n = 11, "Cluster 2", Wong vermillion).
- Plot MAF_C1 (x) vs MAF_C2 (y), one circle marker per unique SNP, coloured
  by attribution, with the two attribution regions shaded.

## Outputs

- `figures/fig_Ss05_attribution_strip.png` / `.pdf`
- `tables/Ss05_outlier_attribution_unique.csv` — the 16 unique SNPs with
  their MAF_C1 / MAF_C2 and attribution.
- Also overwrites the manuscript copy
  `manuscript/paper_A1/share/figures/main/A1_NEW_3a_fig_Ss05_attribution_strip.png` / `.pdf`
  (deposited in the SI table index as **Table S3d**).

## Result

```
[attr] unique SNPs=16: Cluster-2-driven=7, Cluster-1-driven=1, panel-wide=8
```

7 of 16 (44 %) unique Ss05 outliers are Cluster-2-driven, 1 is
Cluster-1-driven, 8 are panel-wide.

## Reproduce

```bash
cd scripts && /Users/black_einstein/miniconda3/envs/yeast-viz/bin/python 95_ss05_attribution_unique_figure.py
```

## Tool versions

- Python 3.9.23 (conda env `yeast-viz`)
- pandas 2.3.1
- matplotlib 3.8.4

## Caveats

- The 16-unique count is the union of the PCAdapt and DAPC top-50 outlier
  sets restricted to Ss05; it is not a genome-wide FDR-significant set
  (neither scan crossed BH-FDR < 0.10). The figure documents *where* the
  scan convergence comes from (within-cluster drift in the 11-line minor
  cluster), not a panel-wide selection claim.
- Attribution is a hard-threshold partition of the per-cluster MAFs; SNPs
  near the 0.40 / 0.10 boundaries are sensitive to the exact cutoff, which is
  inherited unchanged from script 55 for consistency.
