# 06 — Trait extremes, bivariate pairs, and core collections

**Script:** `scripts/06_trait_extremes_corecollection.py` (+ stratified variant in `scripts/08_stratified_core.py`)

## What was done

Per-trait top-10 / bottom-10 accession shortlists, bivariate phenotype-extreme pairs for crossing decisions, and a 20-line MaxMin core collection on Rogers' distance with a stratified variant ensuring representation of the small Cluster-2 subgroup.

## Method

Trait extremes are simple sorted rankings per trait. Bivariate pairs use the 25th and 75th percentiles per trait to define low/high quartiles, then intersect to find accessions in the favourable corner. Core collection: greedy MaxMin (Frankel & Brown 1984; van Hintum 2000) — seed with the most-distant pair, iteratively add the line that maximises minimum distance to the already-selected core. Stratified variant forces 4 lines from the smaller cluster and 16 from the larger, with the size-based allocation auto-detected at run time.

## Quick findings

- Trait-extreme tables for each of the 4 biochem traits (and 13 traits in `17_trait_extremes_13/`).
- 19 bivariate-extreme records for the three breeder-relevant criteria (low Tannin × high Antioxidant; low Phenol × high Flavonoid; high Antioxidant × high Flavonoid).
- Unstratified core: 19 of 20 lines from the large cluster, 1 from the small.
- **Stratified core** ({large = 16, small = 4}): mean pairwise distance 0.36 vs 0.28 in the full panel; min pairwise distance 0.27 (smaller than the unstratified core's 0.31 because the small Cluster-2 lines are pairwise-close among themselves).

## Caveats

- The MaxMin core collection is a one-shot greedy heuristic, not a global optimum; the unstratified version is over-biased toward the dense cluster.
- Both core variants are based on the marker distance matrix alone — no phenotype is used in the selection. A phenotype-weighted core (e.g. CoreHunter, GenoCore) is a future option.
- 8_stratified_core.py auto-detects cluster sizes at runtime so the k-means label flipping between QC versions doesn't break the allocation.
