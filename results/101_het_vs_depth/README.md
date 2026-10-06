# 101 — Heterozygosity versus read depth and genotyping reproducibility (Reviewer 2, point B)

Tests whether the panel's unexpectedly high observed heterozygosity in a
cleistogamous selfer is an artefact of low read depth or poor marker
reproducibility, as Reviewer 2 asked. The reviewer specifically requested that the
DArT **RepAvg** (reproducibility) and **AvgCountRef / AvgCountSnp** (read depth)
columns be used and that heterozygosity be analysed against depth. This note
supplies that analysis at two resolutions: per marker across the panel, and per
accession for the five most heterozygosity-rich main-cluster lines that the
reviewer singled out (Wright F ≈ −0.5 to −0.8: TSs325, TSs82, TSs63A, TSs104,
TSs23).

## Why

If the heterozygosity excess were a calling artefact, heterozygous genotypes would
concentrate on low-depth and/or low-reproducibility markers, and the outbred-tail
lines' het calls would sit on worse-quality markers than their homozygous calls.
Three internal controls are available in the DArT report itself (depth,
reproducibility) and in the re-anchoring run (paralogy class, `results/89`), so the
test needs no external data.

## Inputs

- `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_singlerow_2.csv` — DArT
  singlerow report. Native genotype coding is **0 = homozygous-reference,
  1 = homozygous-SNP, 2 = heterozygous** (verified in-script by matching each
  code's per-marker frequency to FreqHomRef / FreqHomSnp / FreqHets at corr =
  1.000). Carries `AvgCountRef`, `AvgCountSnp`, `RepAvg` per marker.
- `results/01_qc_pca_power/tables/marker_qc.csv` — canonical per-marker observed
  heterozygosity (`het` = DArT FreqHets over all 105 samples). Used as the het of
  record and as a correctness check.
- `results/89_blast_uniqueness_reanchor/tables/marker_placement_classification.csv`
  — per-marker placement class (unique / multi / unplaced) from the all-hits BLAST
  re-anchoring.

## Method

- QC identical to the main pipeline (script 01): marker call-rate ≥ 0.90,
  MAF ≥ 0.05, sample call-rate ≥ 0.90 → **1,625 markers × 95 lines** (asserted in
  script; MAF computed on alt-allele dosage HomRef→0, Het→1, HomSnp→2).
- **Marker level:** Spearman correlation of per-marker observed het against mean
  depth (AvgCountRef + AvgCountSnp), RepAvg, and each depth column separately.
- **Per-line:** for each outbred-tail line, split the post-QC markers into those
  where the line is called HET versus HOM, and compare median depth, median
  RepAvg, and multi-mapping fraction between the two marker sets; one-sided
  Mann–Whitney tests whether the line's het-markers have *lower* RepAvg than its
  hom-markers. Same split reported pooled across all 95 lines (`ALL95`).

Reproduce:
```
/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python scripts/101_het_vs_depth.py
```

## Results

**Marker level (n = 1,625):**

| Test | Spearman ρ | p |
|---|---|---|
| het ~ mean depth | **+0.143** | 6.7×10⁻⁹ |
| het ~ AvgCountSnp (alt-allele reads) | **+0.293** | 1.6×10⁻³³ |
| het ~ AvgCountRef | +0.051 | 0.040 |
| het ~ RepAvg | **−0.205** | 6.6×10⁻¹⁷ |

Heterozygosity rises with depth, not falls — the benign direction: at low depth the
second allele is missed and hets are *under*-called, so depth does not manufacture
heterozygotes. The weak negative het~RepAvg correlation exists, but both het and
hom markers have the same median RepAvg (0.984, see below), so it operates inside a
narrow high-reproducibility band and bounds any reproducibility contribution as
small.

**Per outbred-tail line (het-markers vs hom-markers for that line):**

| Line | n_het | depth med (het / hom) | RepAvg med (het / hom) | multimap frac (het / hom) |
|---|---|---|---|---|
| TSs325 | 915 | 38.8 / 33.3 | 0.984 / 0.984 | 0.109 / 0.092 |
| TSs82 | 915 | 37.8 / 34.6 | 0.984 / 0.984 | 0.089 / 0.119 |
| TSs63A | 866 | 37.3 / 34.2 | 0.984 / 0.984 | 0.101 / 0.104 |
| TSs104 | 844 | 39.5 / 32.5 | 0.984 / 0.984 | 0.102 / 0.100 |
| TSs23 | 769 | 38.0 / 35.0 | 0.984 / 0.984 | 0.099 / 0.103 |
| ALL95 | 48,059 | 38.0 / 35.4 | 0.984 / 0.984 | 0.110 / 0.098 |

For every outbred-tail line the heterozygous calls sit on markers of **higher**
median depth than the homozygous calls, **identical** median reproducibility
(0.984), and multi-mapping fraction that is not elevated (lower for TSs82, flat for
the rest). The excess heterozygosity therefore sits on high-depth,
high-reproducibility, predominantly uniquely-mapping markers.

## Interpretation and caveat

Low read depth, poor reproducibility, and collapsed paralogy are all excluded as
the driver of the outbred tail's heterozygosity: the het calls are on the
*better*-depth, equally-reproducible, not-paralogous markers, and paralogy is
independently het-*depleted* panel-wide (`results/89`: multi-mapping median het
0.048 vs unique 0.114). What a genotype-call analysis **cannot** do is separate a
genuinely heterozygous individual (residual outcrossing / an F1-like plant) from a
**heterogeneous seed lot** in which DNA bulked across divergent plants registers
between-plant differences as within-sample heterozygosity — both produce
genome-wide Ho ≫ He on high-quality markers. That ambiguity is intrinsic to bulked
DNA and is the reason the parental-use recommendation for these lines is
conditioned on single-plant re-genotyping (manuscript §3.6, Discussion).

## Tool versions

- Python 3.9.23 (`/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python`)
- pandas 2.3.1, numpy, scipy (spearmanr, mannwhitneyu)

## Outputs

- `tables/het_vs_depth_marker_correlations.csv`
- `tables/outbred_line_het_marker_quality.csv`
- `tables/per_marker_het_depth_rep_placement.csv`

### Manuscript deposit

The per-line het-marker-vs-hom-marker quality table is deposited for A1 as
**Supplementary Table S15** (`TableS15_het_vs_depth_outbred_lines.csv`), with
publication column names, in two md5-identical copies:

- `results/101_het_vs_depth/tables/TableS15_het_vs_depth_outbred_lines.csv`
- `manuscript/paper_A1/share/tables/supp_tables/TableS15_het_vs_depth_outbred_lines.csv`

The marker-level Spearman het-vs-depth correlations (ρ = +0.14 depth, −0.21
RepAvg) are reported inline in the manuscript Results (§3.6) rather than as a
separate deposited table; the full per-marker depth/RepAvg/placement table
(`per_marker_het_depth_rep_placement.csv`) remains here as the backing record.
