# 29 — Discriminant Analysis of Principal Components

**Script:** `scripts/29_dapc.py`

## What was done

DAPC (Jombart, Devillard & Balloux 2010) on the QC-filtered dosage with the marker-PCA k = 2 cluster labels as discriminant groups. Identifies the linear combinations of markers that maximise between-cluster separation.

## Method

Step 1: PCA on the standardised dosage matrix; the first 30 PCs are retained (covers ~75 % of marker variance). Step 2: Linear Discriminant Analysis (`sklearn.discriminant_analysis.LinearDiscriminantAnalysis`) fitted on the 30 PCs using cluster labels as targets. The LD1 axis = maximum-discrimination direction. Per-marker LD1 loadings recovered as `pca.components_.T @ lda.scalings_[:, 0]` — the contribution of each marker's standardised genotype to LD1.

## Quick findings

- LD1 cleanly separates the two clusters (density histograms barely overlap).
- **Top-5 markers by |LD1 loading| are all on Ss05** (Ss05:9.7, 10.9, 20.4, 51.0 Mb plus closely-linked SNPs).
- **DAPC and PCAdapt independently converge on Ss05 as the most-differentiated chromosome** — two independent methods agree.

## Caveats

- DAPC group labels were learned from PCA (same data) — discriminant axis is by construction maximally-cluster-separating, no independent inference about *real* group structure.
- 30 retained PCs is a reasonable default (~75 % variance) but choosing too few or too many shifts the LD1 weighting. The Ss05 finding is robust across 20–40 PCs retained.
- LD1 marker loadings are signed; the magnitude is what we plot. The Ss05 convergence with PCAdapt (an independent method) is the paper-bearing finding.
