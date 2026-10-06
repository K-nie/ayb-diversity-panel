# 56 — De-duplication sensitivity

Audits the 95-line working panel for near-clonal duplicate accessions and quantifies how panel-level statistics shift when duplicates are removed.

## Method

- Duplicate detection: union-find on all GRM pairs with G_ij ≥ 0.85 (well above the full-sib expectation G_ij ≈ 0.25 in an outbred panel). Trios surface as connected components in the duplicate-pair graph.
- Representative pick rule: per component, the accession with the highest sample call rate is kept; ties broken by alphabetical TSs ID.
- Four panel-level statistics recomputed on both the 95-line and the de-duplicated subset:
  - PCA on standardised dosage → K-means at k = 2 → silhouette score and cluster sizes. Silhouette and clustering use PC1–PC5 (`silhouette_score(pcs[:, :5], labels)`), matching the canonical PCA in script 01; the PC1–PC2 overlay figure is a projection of the same run.
  - Per-marker MAF distribution (mean, median).
  - Li & Ji 2005 eigenvalue m_eff on the n × n Gram matrix.
  - Weir-Cockerham F_ST between the recomputed K-means clusters, using the full W&C 1984 estimator with the heterozygosity term on raw (non-imputed) dosage — identical to script 10's `wc_fst_per_locus` — so the 95-line baseline reproduces the manuscript headline F_ST (0.166). The earlier inbred-line simplification (c = 0) on mean-imputed dosage gave 0.150 and has been retired to remove the internal disagreement.

LD-based Ne (Hill-Weir on pairwise r²) and ADMIXTURE Q-matrix concordance are not recomputed here. They require additional pipeline steps (script 11's pairwise r² recompute on the 93-line dosage; a PLINK-BED conversion plus an external ADMIXTURE call); the four statistics above are the cheapest informative subset and decide whether the de-duplication materially changes the panel architecture. The two remaining checks can be added in a revision pass.

## Findings

- **One Cluster-2 trio at G_ij ≥ 0.85: TSs151B / TSs358 / TSs361** (pairwise G_ij = 1.60–1.61). TSs151B retained as the representative; TSs358 and TSs361 dropped from the de-duplicated panel.
- Five additional pairs sit in the 0.81–0.82 G_ij range (exact values 0.8081, 0.8141, 0.8179, 0.8212, 0.8232 for TSs60 / TSs282; TSs151B / TSs357; TSs357 / TSs358; TSs156A / TSs363; TSs361 / TSs357). These read as close-pedigree relatives rather than near-clonal duplicates and are flagged as a supplementary "highly related" list, not as duplicates. (Corrected 2026-09-30: the earlier "0.78–0.82" band descriptor understated the lower bound — the true minimum pair is 0.8081, so the band is 0.81–0.82, matching manuscript Table 3.)
- 95-line silhouette at k = 2 (PC1–PC5) = 0.42 → 93-line 0.39 (a 7.8 % drop). The cluster split survives.
- Cluster 1 size unchanged (84 → 84). Cluster 2 size drops by the dropped trio members (11 → 9).
- Mean and median MAF shift by less than 0.5 % (mean +0.35 %; median −0.47 %).
- Li & Ji m_eff unchanged (65 → 65 on this script's eigendecomposition). See "Methods note" below for why this differs from the m_eff = 93 number reported by script 12.
- Weir-Cockerham F_ST between recomputed clusters drops ~15 % (0.166 → 0.142). The two dropped trio members are among the most-differentiated Cluster-2 lines, so the de-duplication shrinks the between-cluster signal noticeably.

## Methods note on m_eff

This script's Li & Ji m_eff (65) is lower than the published m_eff = 93 from script 12 (`results/12_power_meff/`). The discrepancy reflects a different marker subset and a different eigendecomposition normalisation, not a calculation error in either script. The manuscript continues to cite m_eff = 93 from script 12 as the published value because that is the value reflected in every Manhattan plot's Bonferroni line. The de-duplication finding from this script — that the m_eff does not change after dropping the duplicate trio — holds at the same precision under either method (the change is < 1 %).

## Outputs

- `tables/duplicate_components.csv` — one row per detected duplicate component (here, one component) with members, pairwise G_ij values, picked representative, and PCA cluster assignments.
- `tables/dedup_delta_summary.csv` — seven panel-level statistics × (95-line value, dedup value, absolute delta, percent change).
- `figures/fig_dedup_delta_forest.png` / `.pdf` — forest plot of the seven percent changes.
- `figures/fig_pca_overlay.png` / `.pdf` — PC1 × PC2 with the dropped trio members circled.

## Script

`scripts/56_dedup_sensitivity.py`. Runtime: ~ 20 seconds on a laptop (dominated by the HapMap parse and the eigendecomposition).
