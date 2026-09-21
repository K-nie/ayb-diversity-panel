# 05 — Pairwise genetic distance + UPGMA dendrogram

**Script:** `scripts/05_distance_dendrogram.py`

## What was done

Rogers' modified genetic distance computed pairwise across the 95 AYB accessions, UPGMA hierarchical clustering, and three dendrogram renderings: linear, circular (radial), and a distance heatmap reordered by linkage leaf order.

## Method

Rogers' modified distance: D_R(i, j) = √[(1 / (8 L)) Σ_l (d_il − d_jl)²] where d is 0/1/2 dosage and L is the number of QC-filtered markers. Computed via the squared-Euclidean shortcut (sq[:,None] + sq[None,:] − 2 A@A.T). UPGMA via scipy `linkage(method="average")`. Leaves coloured by the k = 2 PCA cluster.

## Quick findings

- **Median pairwise Rogers' distance = 0.281**, max = 0.393.
- UPGMA cleanly splits the 11 Cluster-2 outliers as a tight monophyletic subgroup at one end of the tree.
- Distance heatmap (UPGMA-reordered) shows the same partition as a dark-block in one corner.
- Circular dendrogram added at user request.

## Caveats

- Rogers' modified distance is one of several reasonable choices (alternatives: Nei's standard distance, IBS-based distance); the results are robust to the choice for this panel.
- UPGMA assumes ultrametricity; the ML phylogenies (script 19) drop this assumption and give a slightly less constrained topology (cophenetic r = 0.84 between UPGMA and raw Rogers' D vs 0.91 for the ML methods).
- The circular dendrogram is a within-species clustering, not a phylogenetic tree in the comparative-biology sense.
