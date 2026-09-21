# 10 — Weir & Cockerham F_ST per locus and global

**Script:** `scripts/10_fst.py`

## What was done

Weir & Cockerham (1984) F_ST estimator per locus and global, computed between the two PCA-derived K-means clusters (n = 84 vs n = 11), with bootstrap 95 % CI on the global statistic.

## Method

Per-locus variance components a (between-population), b (between-individual within-population) and c (within-individual) computed by the Weir & Cockerham 1984 equations 5–7 (vectorised across loci). Global F_ST = Σ_l a_l / Σ_l (a_l + b_l + c_l) — the **ratio-of-sums** form (preferred over mean-of-per-locus-ratios). Bootstrap 95 % CI by resampling loci with replacement (B = 1,000).

## Quick findings

- **Global F_ST = 0.165**, bootstrap 95 % CI [0.155, 0.176].
- Per-locus F_ST 99th percentile = 0.69 — visually flagged in the Manhattan as outliers.
- Consistent with the AMOVA Φ_ST of 0.273 (different statistic, same partition).

## Caveats

- Cluster labels were derived from PCA on the same SNPs — descriptive of the partition, not independent inference of structure.
- F_ST is not a test of selection in this implementation; the per-locus 99th-pct flag is a descriptive outlier tail. For a proper selection-genomics scan see `27_pcadapt_outliers/`.
- F_ST of 0.165 is in Wright's "moderate" range (0.05–0.25). For a putatively-autogamous landrace, this is on the lower end — consistent with low effective migration / low residual outcrossing.
