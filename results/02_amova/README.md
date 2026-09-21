# 02 — AMOVA between PCA clusters

**Script:** `scripts/02_amova.py`

## What was done

Analysis of Molecular Variance (Excoffier, Smouse & Quattro 1992) on the QC-filtered dosage matrix, partitioning marker variance between the two PCA-derived K-means clusters. The classical df / SS / MS / pseudo-F / p table is the standard reviewer expectation.

## Method

Pairwise squared Euclidean distance D² was computed on the standardised dosage. Variance components: SS_T = (1/N) Σ D²; SS_W (within-group); SS_A = SS_T − SS_W. df_A = G − 1; df_W = N − G. The pseudo-F ratio = MS_A / MS_W. Φ_ST = σ²_A / (σ²_A + σ²_W) with the average-group-size correction (Excoffier 1992 eq. 5). Significance from 999 permutations of group labels; empirical p reported for both Φ_ST and the pseudo-F ratio.

## Quick findings

- **Φ_ST = 0.273** between Cluster 1 (n = 84) and Cluster 2 (n = 11).
- Pseudo-F = **8.39**, permutation p = **0.001** (both Φ_ST and F).
- 27.3 % of variance among groups; 72.7 % within groups.
- AMOVA pie chart was removed at user request — variance components are in the table.

## Caveats

- Cluster labels were learned from the same marker data (PCA + k-means), so AMOVA quantifies the partition already visible in PCA — **not** independent inferential evidence of structure.
- Replace with provenance / geographic origin if such metadata becomes available; the test would then carry true inferential weight.
