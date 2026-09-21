# 90 — ROH robustness: length-class degeneracy, zero-het sensitivity, permutation null

Tests whether the manuscript's ROH inferences survive the sparse-marker-density critique (Reviewer 1 major point 2). Validates the ROH signal that *does* hold and identifies the inference that does not.

## Why

The submitted Methods called ROH with a 20-SNP sliding window and retained segments "≥ 500 kb AND ≥ 8 SNP," then reported F_ROH at four length cutoffs (0.5, 1.0, 2.0, 5.0 Mb) as if they separated ancient from recent autozygosity. At the anchored marker density this is not valid:

- A 20-SNP window forces every call to span ≥ 20 markers, so the "≥ 8 SNP" and "≥ 500 kb" retention bounds never bind — they are inoperative.
- Almost no segments fall between 0.5 and 2 Mb, so the four length-cutoff F_ROH values are numerically identical and cannot resolve segment-length classes.

This script quantifies both, then runs the two robustness checks the reviewer asked for (a zero-het re-call and a permutation null) to separate the claims that survive from the one that does not.

## Inputs

- `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv` — raw genotypes.
- `refs/ayb_genome/ayb_marker_anchoring.csv` — AYB anchoring; ROH runs on the **1,473 AYB-anchored QC markers** (call-rate ≥ 0.90, MAF ≥ 0.05).
- `results/01_qc_pca_power/tables/pca_coords.csv` — the two-cluster labels (small n = 11 vs large n = 84).

## Commands (reproduce)

```bash
python scripts/90_roh_sensitivity_null.py   # env: yeast-viz
```

The ROH caller is a vectorised (cumsum-based) re-implementation of script 34's PLINK-style algorithm; it reproduces script 34 exactly (719 segments, panel F_ROH 0.2560, small-cluster 0.7816, large-cluster 0.1871), so the sensitivity/null runs are on the same call logic as the paper.

## Parameters

Identical to script 34: window = 20 SNPs, ≤ 2 missing/window, hit-fraction ≥ 0.05, retention ≥ 500 kb & ≥ 8 SNPs, F_ROH = ROH span / anchored genome span (631.1 Mb). Sensitivity run flips max-het/window from 1 → 0. Null: 100 within-locus permutations, seed 20260920.

## Findings

**1. Marker spacing is uneven and the "600 kb mean" is wrong.** Inter-SNP spacing on the anchored set is mean **432 kb**, median **85 kb** (the mean is inflated by large inter-cluster gaps). The Methods figure of "600 kb mean SNP density" should be corrected to ~432 kb mean / ~85 kb median, and reported as both statistics because the distribution is strongly right-skewed.

**2. Length classes are degenerate — drop the ancient-vs-recent inference.** Retained ROH have min 0.58 Mb, **median 8.75 Mb**, max 104 Mb; only 3 of 719 segments are < 1 Mb and 45 < 2 Mb. Panel-mean F_ROH is numerically flat across cutoffs:

| Cutoff | F_ROH |
|---|---|
| ≥ 0.5 Mb | 0.2560 |
| ≥ 1.0 Mb | 0.2559 |
| ≥ 2.0 Mb | 0.2549 |
| ≥ 5.0 Mb | 0.2452 |

The 0.5/1.0/2.0 Mb cutoffs are effectively the same measurement. **Recommendation: report a single F_ROH (≥ 1 Mb) and delete the four-cutoff length-class inference and any g_eff / effective-segment argument that depends on segment-length resolution.** The "≥ 500 kb & ≥ 8 SNP" retention text should note the 20-SNP window is the binding constraint.

**3. The autozygosity contrast is robust to het tolerance.** Re-calling with zero heterozygotes per window lowers absolute F_ROH (panel 0.174) but *preserves and sharpens* the cluster contrast: small-cluster 0.666 vs large-cluster 0.110. The bimodal autozygosity landscape is not an artefact of the 1-het tolerance.

**4. The ROH signal vastly exceeds chance homozygosity at this density.** Within-locus permutation (preserving per-SNP MAF and heterozygote count, destroying along-chromosome autozygosity correlation) gives a null panel F_ROH of mean **0.047** (95th percentile 0.056). Observed is **0.256** (empirical p = 0.0099, the floor for 100 permutations); the small cluster's observed 0.782 vs null 0.047 (p = 0.0099). Long ROH therefore reflect genuine contiguous autozygosity, not chance homozygous runs produced by sparse markers.

**Bottom line for the revision:** keep the F_ROH inbreeding story and the small-vs-large cluster autozygosity contrast (both survive sensitivity and the null); collapse to one F_ROH cutoff; **drop the length-class / segment-age inference and any g_eff derived from it**; correct the spacing figure; add a limitations paragraph on marker density.

## Outputs

- `tables/roh_sensitivity_summary.csv` — every headline number above.
- `tables/froh_per_sample_obs_vs_zerohet.csv` — per-accession F_ROH, observed vs zero-het, with cluster label.
- `figures/fig_roh_permutation_null.png` / `.pdf` — observed panel F_ROH vs the permutation null.

## Caveats

- The permutation null breaks *all* along-genome correlation, so it is a conservative test of "are these runs longer/more numerous than independence predicts" — it does not model LD, which would raise the null slightly; the ~5× gap between observed and null is large enough that LD cannot close it.
- 100 permutations set an empirical-p floor of 0.0099; the observed values sit far outside the null envelope, so more permutations would only lower p further.

## Script

`scripts/90_roh_sensitivity_null.py`. Env: `yeast-viz` (numpy, pandas, matplotlib). Runtime: ~1–2 min (100 permutations × 95 samples × 11 chromosomes, vectorised).
