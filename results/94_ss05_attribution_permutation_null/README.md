# 94 — Ss05 outlier cluster-attribution permutation null (de-duplicated unique-SNP definition)

## What this does and why
The submitted A1 manuscript reported that "13 of 20" Ss05 selection outliers were
minor-cluster (Cluster-2) driven and tested that count against a permutation null.
Reviewer R1-M7 (drift attribution) and the R2-C numerical-consistency point forced a
re-audit of that statistic. The audit against
`results/55_private_alleles_per_cluster/tables/Ss05_outlier_attribution.csv` showed the
"13" was a **method × SNP row count** (13 of 22 rows), which double-counts any SNP
flagged by both PCAdapt and DAPC, and it was paired with an inconsistent "20"
denominator that matches neither the 22 rows nor the 16 unique SNPs. Script 55's
recommended, de-duplicated figure is **7 of 16 unique Ss05 outliers Cluster-2-driven**.

This script re-runs the permutation null described in the submitted Methods on that
corrected unique-SNP count, so the manuscript can report a numerically consistent
attribution statistic with a matching null.

## The null (as described in the submitted Methods)
Cluster sizes are held fixed at 11 (Cluster 2, minor) / 84 (Cluster 1, main). For each
of the 16 unique Ss05 outlier SNPs, genotypes for all 95 individuals are re-drawn from
`Binomial(2, p_alt)` under the observed panel-wide alt-allele frequency `p_alt`. This
destroys any real cluster-structure association while preserving each SNP's panel-wide
allele frequency and the 11/84 split. The first 84 draws are assigned to Cluster 1, the
remaining 11 to Cluster 2 (valid because the draws are i.i.d.). The per-permutation
statistic is the number of the 16 SNPs that would be classified **Cluster-2-driven**
(MAF_C2 ≥ 0.40 AND MAF_C1 ≤ 0.10 — the same thresholds used in script 55) by chance.

## Inputs
- `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv` — raw genotypes.
- `results/01_qc_pca_power/tables/pca_coords.csv` — cluster labels
  (convention: `cluster==1` large/main n=84 = "Cluster 1"; `cluster==0` small/minor n=11 = "Cluster 2").
- `results/55_private_alleles_per_cluster/tables/Ss05_outlier_attribution.csv` —
  the 22-row / 16-unique Ss05 outlier set and the observed per-SNP attribution.

## How to reproduce
```
conda activate yeast-viz
cd scripts
python 94_ss05_attribution_permutation_null.py
```
QC filters applied to build the working matrix (match the panel-wide pipeline):
sample call-rate ≥ 0.90, marker call-rate ≥ 0.90, MAF ≥ 0.05 → 1,625 markers × 95 samples.
Attribution thresholds: `ATTR_C_HIGH = 0.40`, `ATTR_C_LOW = 0.10`.
Permutations `N_PERM = 10,000`; `SEED = 20260920`.

## Outputs
- `tables/permutation_null_summary.csv` — observed count, null mean/sd/max, empirical p, seed.
- `figures/fig_ss05_attribution_null.png` / `.pdf` — observed vs null histogram.

## Result
| n unique Ss05 outliers | observed Cluster-2-driven | null mean | null sd | null max | n perms | empirical p |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 7 | 0.0006 | 0.0245 | 1 | 10,000 | ≤ 1.0 × 10⁻⁴ |

The corrected unique-SNP count (7 of 16) departs from the permutation null just as
decisively as the previously reported row-level count: across 10,000 permutations the
null never produced more than a single Cluster-2-driven SNP, so 7 is far outside the
null (empirical p ≤ 0.0001, the add-one one-sided floor at 10,000 permutations).

## Software versions
- Python 3 (conda env `yeast-viz`): numpy, pandas, matplotlib.
- Panel-wide plotting style via `scripts/_plotstyle.py` (Wong colourblind-safe palette).

## Key decisions / caveats
- **Why the unique-SNP count, not the row count:** a SNP detected by both PCAdapt and
  DAPC is one locus, not two; the row-level "13" inflates the attribution count and the
  "20" denominator is not traceable to any panel. 7 of 16 is the defensible figure.
- **Why the null never exceeds 1:** the joint requirement (MAF_C2 ≥ 0.40 *and*
  MAF_C1 ≤ 0.10) is very hard to satisfy by chance in an 11-sample cluster drawn from
  the same allele frequency as the 84-sample cluster, so the null mass sits at 0.
- **Empirical p floor:** with 10,000 permutations and the add-one correction the
  smallest reportable one-sided p is 1/(10,000+1) ≈ 1.0 × 10⁻⁴; report as "≤ 0.0001".
