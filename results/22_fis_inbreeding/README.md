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
