# 97 — Is PC1 driven by heterozygosity itself?

**Script:** `scripts/97_pc1_heterozygosity_robustness.py`

Reviewer robustness check for Item 5 (part B). The panel splits on PC1 into 11 near-homozygous minor-cluster lines and 84 more heterozygous main-cluster lines. A reviewer asked us to rule out the trivial explanation that **PC1 is a heterozygosity axis rather than a population-structure axis**, suggesting a PCA on an IBS matrix or on haploidised calls. This script does both.

## Why

If PC1 merely tracked how many heterozygous sites each line carries, the 11/84 split would be an artefact of genotype-call coding, not of ancestry. Two heterozygosity-insensitive re-derivations settle it:

1. **Haploidised-call PCA** — every heterozygous genotype is resolved to its major-allele homozygote, producing a strictly homozygous 0/2 matrix in which heterozygosity **cannot** contribute to the ordination. This is the same major-allele resolution used to build the phylogenetic SNP alignment (script 19).
2. **IBS PCoA** — principal-coordinate analysis of an identity-by-state (allele-sharing) distance matrix, an ordination of pairwise similarity rather than of per-marker dosage variance.

If the split and the PC1 ordering survive both, PC1 indexes structure.

## Method

- **Input / QC identical to script 01:** `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv` → 0/1/2 dosage; sample call-rate ≥ 0.90, marker call-rate ≥ 0.90 & MAF ≥ 0.05 → **1,625 markers × 95 lines**. Missing calls mean-imputed per marker.
- **Baseline:** column-standardised dosage, `sklearn` PCA; k-means (k = 2, `n_init` = 25) labels the 11/84 split. Reproduces manuscript Figure 2A (PC1 = 8.5 % variance, minor n = 11, main n = 84).
- **Per-sample observed heterozygosity** computed on the post-QC marker set; correlated with baseline PC1 (Pearson, Spearman) and its R² reported.
- **Haploidised PCA:** het (dosage = 1) → major-allele homozygote per marker (alt-frequency > 0.5 → 2, else 0); mean-impute; standardise; drop zero-variance columns (1,429 informative positions remain); PCA + k-means; agreement with baseline labels by adjusted Rand index (ARI) and PC1↔PC1 Pearson r.
- **IBS PCoA:** pairwise IBS distance = mean |dosage_i − dosage_j| / 2 on the imputed matrix; classical MDS (double-centred eigendecomposition); k-means on the leading two coordinates; ARI and PCoA1↔baseline-PC1 Pearson r.
- Axes oriented so the minor cluster sits on the positive side for readability.

## Findings

| Quantity | Value |
|---|---|
| Baseline PC1 variance | 8.5 % (minor n = 11, main n = 84) |
| Mean observed het, minor vs main cluster | 0.086 vs 0.341 |
| **PC1 vs observed het** | Pearson r = **−0.587** (p = 4.1 × 10⁻¹⁰); Spearman ρ = −0.469; **R² = 0.34** |
| **Haploidised-call PCA** (1,429 strictly homozygous positions) | PC1 = 8.7 % variance; **ARI = 0.865** vs baseline k-means; **PC1↔baseline-PC1 r = 0.874** |
| **IBS PCoA** | PCoA1 = 17.3 % variance; ARI = 0.571; **PCoA1↔baseline-PC1 r = 0.879** |

**Interpretation:** PC1 *is* correlated with heterozygosity (r = −0.59), but that reflects a genuine biological contrast — the minor cluster is far more autozygous (mean het 0.086 vs 0.341), which is itself part of the structure. It is **not** an artefact: a strictly heterozygosity-free representation (haploidised calls) recovers the same 11/84 split at ARI 0.87 and reproduces PC1 at r = 0.87, and an allele-sharing IBS ordination reproduces it on the leading axis at r = 0.88. PC1 therefore indexes population structure, not the heterozygosity coding.

## Outputs

- `tables/pc1_het_summary.csv` — per-sample PC1 (baseline / haploid / IBS-PCoA1), observed het, cluster label.
- `tables/concordance_summary.csv` — all correlation / ARI / recovery statistics above.
- `figures/fig_pc1_heterozygosity.{png,pdf}` — 4-panel supplementary figure (**manuscript Supplementary Figure S16**): (a) baseline PCA coloured by het; (b) PC1 vs het; (c) haploidised-call PCA; (d) IBS PCoA.

## Reproduce

```bash
cd scripts && /Users/black_einstein/miniconda3/envs/yeast-viz/bin/python 97_pc1_heterozygosity_robustness.py
```

## Caveats

- IBS k-means agreement (ARI 0.57) is lower than the leading-axis correlation (r = 0.88) because 2-means on the IBS coordinates draws the boundary slightly differently near the near-clonal trio; the leading-axis correlation is the relevant statistic for "does the split survive."
- Haploidisation removes genuine heterozygosity information; it is a stress test of PC1, not the analysis of record for structure (that remains the standardised-dosage PCA, Figure 2A).

## Tool versions

Python 3.9.23 (conda env `yeast-viz`); numpy, pandas, scipy, scikit-learn, matplotlib, Pillow.
