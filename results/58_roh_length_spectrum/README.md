# 58 — F_ROH segment-length spectrum

Stratifies ROH segments from script 34 into four length classes and computes per-accession F_ROH within each class plus the inferred recent effective selfing generations (g_eff). Closes the F_IS / F_ROH paradox into a per-accession reproductive-biology readout for Paper A1.

## Inputs

- `results/34_roh/tables/roh_calls.csv` — 719 ROH segments across 89 accessions (6 outbred accessions carry zero ROH).
- `results/34_roh/tables/f_roh_per_sample.csv` — full 95-accession panel sample list.

## Method

- Segments stratified into: short (0.5–1 Mb), intermediate (1–5 Mb), long (5–20 Mb), very long (> 20 Mb).
- Per-accession F_ROH per class as the within-class summed segment length divided by the 649.8 Mb AYB reference span.
- Recent effective selfing generations per accession from the Kardos et al. 2018 expectation, converted at 1.06 cM/Mb (Lonardi 2019 cowpea proxy map): `g_eff = 47.17 / mean_long_segment_length_in_Mb`.

## Outputs

- `tables/per_accession_roh_class_spectrum.csv` — 95 rows; per-class F_ROH, segment counts, and per-class g_eff.
- `tables/panel_summary.csv` — 4 rows; panel-level mean / median F_ROH per class, plus segment counts.
- `figures/fig_roh_length_class_spectrum.png` / `.pdf` — stacked-bar per accession, ordered by total F_ROH descending.
- `figures/fig_geff_distribution.png` / `.pdf` — log-axis histogram of g_eff across the panel, with 1-, 5-, and 10-generation reference lines.

## Headline numbers (for the manuscript)

- Mean F_ROH ancestral (< 1 Mb)         : 0.000 (3 segments across the panel)
- Mean F_ROH intermediate (1–5 Mb)      : 0.010
- Mean F_ROH long (5–20 Mb)             : 0.043
- Mean F_ROH very long (> 20 Mb)        : 0.195 (dominant contribution)
- Accessions with at least one > 20 Mb ROH : 58 / 95 (61.1 %)
- Median g_eff from long (5–20 Mb)      : 5.3 generations
- Median g_eff from very long (> 20 Mb) : 1.1 generations
- Outbred accessions (zero ROH)         : TSs15, TSs23, TSs294, TSs325, TSs63A, TSs82

## Caveats

- The 0.5 Mb floor follows script 34's PLINK-style minimum; ROH below that length is not resolved at the panel's ~440 kb median marker spacing.
- The 1.06 cM/Mb map rate is a cowpea proxy; the g_eff estimates are ordinal rather than absolute generation counts.
- The Kardos formula assumes ROH lengths are sampled from the gamma expectation; segment-level noise at small sample size makes individual g_eff estimates uncertain, but the panel-level median is robust.

## Script

`scripts/58_roh_length_spectrum.py`. Runtime: ~ 2 seconds on a laptop.
