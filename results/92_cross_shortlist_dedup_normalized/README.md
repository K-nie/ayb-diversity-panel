# 92 — De-duplicated, inbreeding-normalized cross shortlist

Rebuilds the breeding cross shortlist after (i) collapsing the cryptic near-clonal duplicates and (ii) replacing raw VanRaden G_ij with an inbreeding-normalized relatedness. Addresses Reviewer 2 point D.

## Why

The submitted 50-pair "least-related" cross shortlist (Data Sheet from script 13) was ranked on **raw VanRaden (2008) G_ij** off-diagonals over the **full 95-line panel**, including the near-clonal trio (TSs151B, TSs358, TSs361). Two defects the reviewer flagged:

1. **Duplicates were not removed before shortlisting.** The trio are cryptic near-clones (pairwise G_ij ≈ 1.60; see script 13 / Table 3). Leaving all three in the candidate set lets a nominally "diverse" cross pair one line against a near-identical copy of another candidate, and overstates the number of genetically independent lines available.
2. **Raw G_ij is inbreeding-scaled.** Because G_ii = 1 + F_i, a highly inbred line has an inflated self-relationship and is pushed to extreme off-diagonal values against everyone. The trio's G_ii ≈ 1.6 (implied F ≈ 0.6) versus the panel-mean diagonal of 1.008 (F ≈ 0.008). Ranking on raw G_ij therefore confounds genuine complementarity with each line's own inbreeding level — so the "least-related" tail is dominated by the most-inbred lines rather than the least co-ancestral pairs.

The quantitative signature of the defect: **32 of the submitted 50 least-related pairs involve a trio member, and 21 of the 50 involve one of the two duplicate copies that should have been dropped.** The shortlist was largely an artefact of the trio's inbreeding, not a ranking of complementary parents.

## Inputs

- `results/13_grm_crosspairs/tables/grm_vanraden.csv` — the 95 × 95 VanRaden G (with diagonal), unchanged from the submission.
- `results/13_grm_crosspairs/tables/cross_pairs_least_related.csv` — the submitted 50-pair shortlist (for reconciliation).
- `results/01_qc_pca_power/tables/sample_qc.csv` — per-sample marker call rate (chooses the duplicate representative).
- `results/01_qc_pca_power/tables/pca_coords.csv` — two-cluster labels.
- Phenotypes via `scripts/_pheno.load_phenotypes` (13-trait BLUP table) for the trait-gap context columns.

## Commands (reproduce)

```bash
python scripts/92_cross_shortlist_dedup_normalized.py   # env: yeast-viz
```

## Method

- **De-duplication.** Duplicate pairs are those with G_ij ≥ 0.85 (the same threshold used for duplicate detection in the paper). Transitive closure on the duplicate-pair graph yields one component — the trio {TSs151B, TSs358, TSs361}. The component is collapsed to a **single representative, the member with the highest marker call rate** (TSs151B, call rate 0.966); the other two (TSs358, TSs361) are dropped. Panel: 95 → **93 lines**.
- **Inbreeding-normalized relatedness.** r_ij = G_ij / √(G_ii·G_jj), the genomic analogue of a correlation coefficient. Dividing by each line's own self-relationship removes the inbreeding scaling, so r_ij is bounded and comparable across lines of differing F; "least related" then means genuinely least co-ancestry rather than "most inbred".
- The 50 least-related pairs are re-ranked on the 93-line panel by r_ij, with raw G_ij retained alongside for reference and the 13-trait mid-parent / absolute-gap context carried through.

## Findings

| Metric | Value |
|---|---|
| Panel size (submitted → de-duplicated) | 95 → 93 |
| Duplicate lines dropped | TSs358, TSs361 |
| Representative kept | TSs151B (highest call rate, 0.966) |
| **Submitted shortlist pairs involving a trio member** | **32 / 50** |
| Submitted shortlist pairs involving a *dropped* duplicate | 21 / 50 |
| **New shortlist pairs involving a dropped duplicate** | **0 / 50** |
| Shortlist pairs shared (submitted ∩ new) | 17 / 50 |
| Shortlist pairs replaced | 33 / 50 |
| New normalized r_ij range | −0.237 to −0.182 |

**The reviewer is correct and the corrected shortlist changes materially.** Two-thirds of the submitted "most diverse" pairs (33/50) are replaced once duplicates are collapsed and relatedness is inbreeding-normalized. Every pair that had leaned on a now-dropped duplicate copy (21) is gone. The 17 pairs that survive are genuinely low-co-ancestry crosses robust to the correction. The revision should replace the submitted shortlist Data Sheet with `cross_pairs_least_related_dedup_normalized.csv` and state in Methods that the shortlist is built on the de-duplicated panel using inbreeding-normalized relatedness r_ij = G_ij/√(G_ii G_jj).

## Outputs

- `tables/cross_pairs_least_related_dedup_normalized.csv` — the corrected 50-pair shortlist (r_ij, raw G_ij, cluster labels, 13-trait mid-parent + absolute-gap context). **Replaces the submitted least-related Data Sheet.**
- `tables/shortlist_reconciliation_summary.csv` — every headline number above.
- `figures/fig_cross_shortlist_rawG_vs_normalized.png` / `.pdf` — raw G_ij vs normalized r_ij for all pairs, with trio-involving pairs highlighted to show why the submitted ranking was trio-dominated.

## Caveats

- r_ij normalizes by each line's genomic self-relationship, which is estimated from the same marker set; for the least-inbred lines (G_ii near 1) it is numerically close to raw G_ij, so the reshuffle is driven mainly by removing the inbreeding-inflated lines from the tail, not by a wholesale re-scaling.
- The shortlist is a *screening* tool: it flags low-co-ancestry pairs with wide trait gaps for a breeder to consider, and still requires denser genotyping / controlled crosses before use, exactly as the manuscript states.

## Script

`scripts/92_cross_shortlist_dedup_normalized.py`. Env: `yeast-viz` (numpy, pandas, matplotlib; `_pheno`, `_plotstyle`). Runtime: < 20 s (dominated by the phenotype BLUP load).
