# 13 — Genomic relationship matrix + cross-pair shortlists

**Script:** `scripts/13_grm_crosspairs.py`

## What was done

VanRaden method-1 genomic relationship matrix **G** for the 95-line panel, UPGMA-reordered heatmap rendering, off-diagonal distribution by within-vs-between-cluster pair, and 50-pair shortlists for both least-related (outcrossing) and most-related (line maintenance) cross decisions, joined with full 13-trait phenotype data per pair.

## Method

G = Z Z' / (2 Σ_j p_j (1 − p_j)) with Z = M − 2P (centred dosage; VanRaden 2008). Reordering by hierarchical clustering on (1 − G). Each top-50 pair is joined to the unified 13-trait phenotype frame for breeder review; both A and B parents' trait values plus midparent and absolute trait gap are reported.

## Quick findings

- Diagonal mean = 1.008, off-diagonal mean = −0.011 (correctly centred).
- **Four sets of near-clonal duplicate accessions surfaced**: TSs377 = TSs297 (G = 1.62); Cluster-2 trio {TSs361, TSs358, TSs151B} (G ≈ 1.61); TSs89 ≈ TSs138 (G = 1.09); Cluster-1 trio {TSs84, TSs89, TSs138} (G > 0.85).
- Top-50 least-related: G range [−0.32, −0.23]. These are the maximum-Mendelian-sampling outcrosses.
- Top-50 most-related: G range [0.43, 1.62]. Useful for line maintenance / seed rejuvenation, not breeding gain.

## Caveats

- Duplicates collapse the breeder-elite recommendations. TSs361 (top high-protein × low-oxalate candidate) is a member of a near-clonal trio with TSs358 and TSs151B — it represents one selection, not three.
- Off-diagonal G > 1 is mathematically possible under method-1 (the values are *not* identity-by-descent probabilities; they are scaled covariances) and indicates near-identical genotypes.
- The cross-pair tables don't use trait midparent variance (the Allier 2019 / Lehermeier 2017 "usefulness criterion" would require it but needs replicate phenotypes — blocked).
