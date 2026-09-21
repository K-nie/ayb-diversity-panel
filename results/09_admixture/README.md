# 09 — ADMIXTURE K = 2..6 with 5-fold cross-validation

**Script:** `scripts/09_admixture.py`

## What was done

Unsupervised ADMIXTURE (Alexander et al. 2009) clustering of the 95-line AYB panel across K = 2 to 6 with 5-fold cross-validation, plus stacked-bar Q-matrix renderings ordered by PCA cluster and within-cluster Q dominance.

## Method

PLINK 1.9 was used to convert the QC-filtered dosage matrix (95 × 1,625) to PED → BED format. ADMIXTURE 1.3.0 with `--cv=5` ran for K = {2, 3, 4, 5, 6}. CV error was parsed from each log; "best K" was selected by minimum CV. Q matrices were loaded and rendered as stacked bars, with within-cluster sample ordering by the dominant ancestry component.

## Quick findings

- CV error monotonically falling: K=2 = 0.521, K=3 = 0.514, K=4 = 0.521, K=5 = 0.519, K=6 = 0.530.
- **Best K by CV minimum = 3** (CV = 0.514) — a single clear inflection, unlike the earlier 0.80-QC run which had no clear elbow.
- K = 2 visually matches the PCA cluster split exactly. K = 3..6 subdivide the large Cluster 1 progressively but show signs of over-fitting at K ≥ 5.

## Caveats

- The CV-error span (0.514 to 0.530) is narrow (~3 %); the inflection at K = 3 is real but not strong. K = 2 remains the cleanest interpretive choice.
- ADMIXTURE Q-matrix interpretation can shift across runs due to label switching; the Q matrices here use a fixed random seed.
- The small Cluster 2 (n = 11) is consistently identified as a discrete ancestry component at K ≥ 2.
