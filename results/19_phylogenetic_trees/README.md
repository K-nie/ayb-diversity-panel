# 19 — ML phylogenetic inference (IQ-TREE + RAxML) + cophenetic comparison + tanglegrams

**Scripts:** `scripts/19_phylogenetic_trees.py` + `scripts/30_tanglegrams.R`

## What was done

Maximum-likelihood phylogeny inferred independently by IQ-TREE 3.0.1 and RAxML-NG, cross-validated against the UPGMA dendrogram and the raw Rogers' distance via cophenetic correlation. Tanglegrams between each tree pair via R `dendextend`.

## Method

A haploid SNP alignment was built from the QC-filtered dosage: hets resolved to the major allele at each site, sites that became invariant after resolution removed (196 sites dropped, 1,429 variant sites retained). IQ-TREE: `-m GTR+ASC -B 1000` (1000 ultrafast bootstraps); RAxML-NG: `--all --model GTR+G+ASC_LEWIS --bs-trees 100`. Trees parsed via Bio.Phylo, cophenetic distance matrices extracted, pairwise Pearson r computed across 4 matrices (IQ-TREE, RAxML, UPGMA, raw Rogers'). Tanglegrams via `dendextend::tanglegram` with `dendextend::untangle(method = "step2side")`; entanglement score (Galili 2015); Robinson-Foulds via `ape::dist.topo`.

## Quick findings

- IQ-TREE: 95 tips × 1,429 sites, GTR+ASC, log-likelihood −41,927.3, runtime 3 min 47 s.
- RAxML-NG: same alignment, GTR+G+ASC_LEWIS, log-likelihood −40,481.1, runtime 20 min 50 s.
- **Cophenetic correlations**: IQ-TREE ↔ RAxML r = **0.997**; ML ↔ raw Rogers' D r = 0.91; UPGMA ↔ raw Rogers' D r = 0.84.
- Tanglegram entanglement: IQ-TREE vs RAxML 0.24; IQ-TREE vs UPGMA 0.15; RAxML vs UPGMA 0.22.
- RF distances: 124 / 184 (IQ vs RAx), 161 / 184 (IQ vs UPGMA), 155 / 184 (RAx vs UPGMA).
- Both ML trees collapse the 11 Cluster-2 outliers as a monophyletic group at one side of the tree.

## Caveats

- "Phylogenetic tree" here is within-species (95 lines of one species), not a comparative-biology phylogeny. Frame as a within-panel clustering / genealogy.
- Major-allele imputation of hets is acceptable for an autogamous selfer but introduces a small bias toward ref-allele homozygous calls.
- RAxML flagged 11 near-zero branches in the best ML tree — consistent with the near-clonal pairs already identified via GRM (`13_grm_crosspairs`).
- Robinson-Foulds distances are inflated at n = 95 with low sequence diversity (fine-grain branch swaps drive most of the RF count); not a reliable measure of "tree disagreement" in this regime. Cophenetic r and entanglement are the right summaries.
