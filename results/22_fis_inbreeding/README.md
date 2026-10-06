# 22 — Wright F_IS + per-accession inbreeding F

**Script:** `scripts/22_fis_inbreeding.py`

## What was done

Wright's F_IS per locus and global (between-cluster-independent — purely heterozygosity deficit relative to Hardy-Weinberg), plus per-accession inbreeding coefficient F (Yang et al. 2010 method I).

## Method

Per-locus H_obs = (het count) / (called samples), H_exp = 2pq (HWE), F_IS = 1 − H_obs / H_exp. Global F_IS = (Σ H_exp − Σ H_obs) / Σ H_exp (ratio-of-sums, Weir 1996). Bootstrap 95 % CI by locus-resampling (B = 1,000). Per-accession F = 1 − H_obs_i / E[H_exp_i], with H_exp_i averaged over loci where the sample is called.

## Quick findings

- **Global F_IS = 0.011, bootstrap 95 % CI [0.003, 0.020]** — unexpectedly low for a putatively-autogamous species.
- Per-locus F_IS distribution centred near zero with a long right tail (some loci with F_IS > 0.5).
- Per-accession F median = −0.033 (slightly negative; "more het than expected"), range [−0.799, 0.976]. Wide range with no clear cluster-by-F pattern.

## Caveats

- The low global F_IS conflicts with the AYB selfing-rate expectation (typically F_IS > 0.5 for selfers). Three possible explanations:
  1. DArTseq het-calling bias toward over-reporting heterozygosity at low coverage (Sansaloni et al. 2011 documented this).
  2. Genebank-maintenance / regeneration practices at IITA that introduce residual heterozygosity.
  3. Genuine low residual outcrossing in the panel — possible but contradicts the literature.
- The result matches Shitta et al. 2022 IITA AYB collection (F_IS ≈ 0.01–0.02), suggesting it is a property of the IITA TSs panel rather than our QC choices.
- Per-accession F has limited interpretive value at this density (the wide range is mostly noise at n = 1,625 markers per accession).

## Revision note — R2-C cluster-label standardisation (2026-09-21)

Reviewer R2-C asked for a single, consistent cluster numbering across the
paper. The panel PCA (`pca_coords.csv`) originally coded the **main**
cluster (n = 84) as `cluster == 1` and the **minor** cluster (n = 11) as
`cluster == 0`. The manuscript standard is now **Cluster 1 = main (n = 84,
Wong blue)** and **Cluster 2 = minor (n = 11, Wong vermillion)**.

- Change: the plotting cluster field was remapped from `clu + 1` to
  `2 - clu` (script line ~174), so `cluster == 0` (minor) → 2 and
  `cluster == 1` (main) → 1.
- Effect: **`fig57_f_per_accession` only** — the near-clonal trio and the
  other high-F minor-cluster lines now render **vermillion (Cluster 2)** and
  the main-cluster bulk **blue (Cluster 1)**, matching Fig 7C and every other
  cluster-coloured panel in the paper.
- No statistic changed: global F_IS (0.011), the bootstrap CI, per-locus
  F_IS, and the per-accession F values are cluster-label-independent. Only
  the legend/colour assignment moved.
