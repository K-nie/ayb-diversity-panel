# 96 — LD decay / Ne on the de-duplicated panel with the sample-size-corrected Hill–Weir fit

**Script:** `scripts/96_ld_decay_dedup_sscorr.py`

Second-pass refinement of the genome-wide and per-chromosome LD-decay / effective-population-size (Ne) estimates for Reviewer 1 major point 6. Supersedes the fit in script 91 (genome-wide) and script 11b (per-chromosome).

## Why

Script 91 refit the genome-wide decay on the AYB anchoring (replacing the legacy cowpea-anchored fit) but left two things the reviewer flagged:

1. **No finite-sample correction.** It fit the plain Hill & Weir (1988) drift–recombination expectation `E[r²] = (10+C)/((2+C)(11+C))`, which omits the sampling term. At n = 95 the expected r² between unlinked sites under no LD is ≈ 1/n ≈ 0.011, a non-negligible floor that biases ρ̂ (and hence Ne) downward and lengthens the apparent half-decay.
2. **Near-clonal duplicates left in the panel.** The near-clonal Cluster-2 trio (TSs151B / TSs358 / TSs361; pairwise Gij ≈ 1.60, F_ROH ≈ 0.99; manuscript Table 3) are effectively the same genotype replicated three times. Near-identical lines inflate LD.

This script addresses both, on the same AYB-anchored marker set, so the panel-wide summary is internally consistent with the rest of the revision.

## Method

- **Genotypes:** `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv`, converted to 0/1/2 dosage; markers filtered call-rate ≥ 0.90 & MAF ≥ 0.05, samples call-rate ≥ 0.90 (identical to scripts 87/91) → 1,473 AYB-anchored markers, 95 lines.
- **De-duplication (point ii):** collapse the near-clonal trio to a single representative — **keep TSs151B, drop TSs358 and TSs361 → n = 93.** The near-clonal set is the panel's own duplicate call at Gij ≥ 0.85; the five *close-pedigree* pairs at Gij 0.81–0.82 are **not** duplicates and are retained.
- **Sample-size correction (point i):** fit the Hill & Weir (1988) expectation with Weir's finite-sample term (as in Remington et al. 2001, *PNAS* 98:11479):
  `E[r²] = [(10+C)/((2+C)(11+C))] × [1 + ((3+C)(12+12C+C²)) / (n(2+C)(11+C))]`, `C = ρ·d(Morgans)`, `n` = number of lines. For an autogamous panel of near-homozygous lines `n` is taken as the line count (93), not 2 × individuals.
- **Fit:** ≤ 5 Mb within-chromosome pairs, 40 equal-width distance bins, mean r² per bin, `scipy.optimize.curve_fit` (bounds 1e-3–1e6). Ne = ρ̂ / 4 under the **1.06 cM/Mb cowpea map-rate proxy** (unchanged, so absolute Ne stays comparable to scripts 91/11b).

## Decisions (documented for the reviewer response)

- **Marker set held fixed at the 1,473-marker panel** (MAF/call-rate computed on the full n = 95). Removing two lines is not allowed to change the marker panel, so the fit isolates the effect of the duplicate lines + sampling term from any change in markers. *Sensitivity:* re-applying the MAF filter on n = 93 drops ~29 markers and moves Ne by < 3 % (467 vs 476 uncorrected) — reported here, not adopted.
- **Structure retained (not "within-main-cluster").** The reviewer offered de-duplication *or* a within-cluster fit; we took de-duplication. Residual LD inflation from the two-cluster structure therefore remains and is caveated in the manuscript: the value is an LD-based effective size for this panel, not a species demographic Ne.
- **Map-rate proxy unchanged.** 1.06 cM/Mb (cowpea) is an acknowledged absolute-scale caveat; it does not affect the relative decay profile or the physical half-decay distance.

## Findings

**Genome-wide (dedup n = 93, sample-size corrected):**

| Quantity | This run (96) | Script 91 (full n=95, no corr) | Submitted (cowpea-anchored) |
|---|---|---|---|
| ρ̂ = 4Ne | **2,446** | 1,904 | 2,636 |
| **Ne (1.06 cM/Mb)** | **≈ 611** | ≈ 476 | ≈ 659 |
| Half-decay | **84 kb** | 102 kb | 74 kb |
| r² < 0.1 by | **330 kb** | 370 kb | 267 kb |
| pairs / markers | 32,507 / 1,473 | 32,507 / 1,473 | — |

Both adjustments push Ne **up** and shorten the half-decay: the sample-size term is the larger driver (476 → 582 on the full panel), de-duplication adds a smaller lift (582 → 611). The corrected value (611) lands close to the originally submitted cowpea figure (659) — coincidentally, because the AYB re-anchoring lowered Ne and the sampling correction raised it back.

**Per-chromosome (dedup n = 93, sample-size corrected):** ρ̂ 1,806 (Ss01) – 4,818 (Ss07); Ne ≈ **452 (Ss01) – 1,205 (Ss07)**; half-decay 42 kb (Ss07) – 113 kb (Ss01); **median Ne 585.** The genome-wide 611 sits near this per-chromosome median.

**Bin resolution (reviewer ask — pairs per bin):** the first 125 kb bin holds 2,014 pairs (mean r² 0.339), the second 1,463 pairs (0.092); counts decline to ≥ 505 pairs per 125 kb bin across the 5 Mb window. All 40 bins in `ld_genomewide_dedup_sscorr_bins.csv`.

## Outputs

- `tables/ld_genomewide_dedup_sscorr_summary.csv` — ρ̂, Ne, half-decay, r²<0.1, n pairs/markers.
- `tables/ld_genomewide_dedup_sscorr_bins.csv` — 40-bin decay curve (mean/std/count/mid_bp).
- `tables/ld_per_chrom_dedup_sscorr_summary.csv` — per-chromosome ρ̂/Ne/half-decay/r²<0.1.
- `figures/fig_ld_decay_genomewide_ayb.{png,pdf}` — **replaces manuscript Figure 3** (share `main/8_fig_ld_decay_genomewide_ayb.*`).
- `figures/fig28_ld_decay_per_chrom_ayb.{png,pdf}` — **replaces Supplementary Figure S4A** (share `supp/9_fig28_ld_decay_per_chrom.*`).

## Reproduce

```bash
cd scripts && /Users/black_einstein/miniconda3/envs/yeast-viz/bin/python 96_ld_decay_dedup_sscorr.py
```

## Caveats

- Absolute Ne is proportional to the assumed 1.06 cM/Mb rate (no AYB genetic map exists); the relative decay profile and physical half-decay are unaffected.
- De-duplication removes the near-clonal trio only; broader population structure (the n = 84 vs n = 11 cluster split) still inflates LD, so Ne is a panel-specific LD-based effective size, not a demographic parameter.
- A citation for the finite-sample correction (Weir 1979; Remington et al. 2001, *PNAS*) should be added to the manuscript reference list (tracked as an Item-15 reference follow-up).

## Tool versions

Python 3.9.23 (conda env `yeast-viz`); numpy, pandas 2.3.1, scipy, matplotlib 3.8.4.
