# 93 — m_eff is bounded by sample size, not marker independence

Demonstrates that the submitted effective-number-of-tests m_eff ≈ 93 is a sample-size artefact (rank-bounded by n − 1), not a count of independent markers, and derives the correct multiple-testing threshold. Addresses Reviewer 2 point F.

## Why

The submission reported m_eff = 93 (Li & Ji 2005 eigenvalue method) and used it to set the multiple-testing threshold α = 0.05 / m_eff = 5.4 × 10⁻⁴. But m_eff was computed from the eigenvalue spectrum of a marker correlation matrix estimated from **n = 95 individuals**. Any such matrix has rank ≤ n − 1, so it has at most **n − 1 = 94 non-zero eigenvalues** no matter how many markers enter it. m_eff = 93 therefore sits at the sample-size ceiling and cannot be read as "the number of independent markers" — it is essentially the number of individuals. The reviewer is correct.

## Inputs

- `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv` — raw genotypes; the same QC cascade as the paper (sample & marker call-rate ≥ 0.90, MAF ≥ 0.05) yields **n = 95 individuals × m = 1,625 markers**.

## Commands (reproduce)

```bash
python scripts/93_meff_samplesize_artefact.py   # env: yeast-viz
```

## Method

Li & Ji (2005) m_eff = Σ_i [ I(λ_i ≥ 1) + (λ_i − ⌊λ_i⌋) ] over the eigenvalues λ_i of the marker correlation matrix, computed (as in the submission) via the n × n Gram of standardised genotype columns — whose non-zero eigenvalues equal those of the m × m correlation matrix; the remaining m − n eigenvalues are exactly zero and contribute nothing. Three tests:

1. **Full-panel m_eff** — reproduce the submitted value.
2. **m_eff vs subsample size** — recompute m_eff on random subsamples of n′ ∈ {20, 30, 40, 60, 80, 95} individuals (20 replicates each). If m_eff measured marker LD it would be roughly flat; if it is a sample-size artefact it tracks n′ − 1.
3. **LD-destroyed control** — permute each marker's genotypes across individuals (destroying all LD) and recompute. A genuine independence measure would rise toward m = 1,625; a rank-bounded one stays at ≈ n − 1.

## Findings

| Quantity | Value |
|---|---|
| Post-QC panel | n = 95 individuals, m = 1,625 markers |
| Sample-size ceiling (n − 1) | 94 |
| **Full-panel m_eff (Li & Ji)** | **92.4** (ratio to ceiling 0.983) |
| m_eff with all LD destroyed (permuted) | **94.0** (not → 1,625) |
| Submitted Bonferroni at m_eff | α = 5.4 × 10⁻⁴ (−log₁₀ = 3.27) |
| **Correct Bonferroni at true m** | **α = 3.1 × 10⁻⁵ (−log₁₀ = 4.51)** |

m_eff vs subsample size (mean over 20 reps):

| n′ | n′ − 1 | m_eff | m_eff / (n′ − 1) |
|---|---|---|---|
| 20 | 19 | 19.0 | 0.999 |
| 30 | 29 | 28.9 | 0.997 |
| 40 | 39 | 38.9 | 0.997 |
| 60 | 59 | 58.6 | 0.993 |
| 80 | 79 | 78.0 | 0.987 |
| 95 | 94 | 92.4 | 0.983 |

**m_eff tracks the number of individuals almost exactly, at every sample size, and is completely insensitive to LD.** Even after destroying all marker-to-marker correlation by permutation, m_eff stays at 94 (the n − 1 ceiling) rather than climbing toward the 1,625 markers — the opposite of what a true independence measure would do. The submitted m_eff = 93 reproduces here as 92.4 (the small difference is the current post-QC marker count 1,625 vs the earlier snapshot); either way it is pinned at ~n − 1.

**Consequence for the threshold.** The submitted m_eff-scale threshold (−log₁₀ = 3.27) is ~1.2 log-units more lenient than a Bonferroni correction on the true marker count (−log₁₀ = 4.51). Because m_eff is not a valid independence count at this n, the revision should:

- report GWAS / outlier significance at the **true marker-count Bonferroni** (α = 0.05 / m = 3.1 × 10⁻⁵) and/or an FDR / permutation threshold, and
- either **drop the m_eff-scale threshold** or report it only with an explicit statement that at n = 95 it is bounded above by n − 1 and therefore cannot be interpreted as the number of independent markers.

The Methods sentence "m_eff was derived from … the n × n Gram" should be amended to note the n − 1 rank bound.

## Outputs

- `tables/meff_vs_samplesize.csv` — m_eff (mean, sd) vs n′ with the n′ − 1 ceiling and ratio.
- `tables/meff_artefact_summary.csv` — every headline number, both Bonferroni thresholds.
- `figures/fig_meff_vs_samplesize.png` / `.pdf` — m_eff vs n′ against the n − 1 line, with the LD-destroyed control and the true marker-count ceiling annotated.

## Caveats

- This does not invalidate any GWAS/outlier hit that already clears the stricter marker-count Bonferroni; it only removes the over-lenient m_eff-scale threshold. Hits between the two thresholds move from "significant" to "suggestive".
- The Daetwyler M_e / prediction-ceiling parts of script 12 use the LD-derived Ne, which was itself refit on the AYB anchoring (script 91, Ne ≈ 476 vs the legacy 659); that is a separate correction (R1-M6) and only rescales M_e, not the m_eff argument here.

## Script

`scripts/93_meff_samplesize_artefact.py`. Env: `yeast-viz` (numpy, pandas, matplotlib). Runtime: ~1 min (dosage rebuild + 100 eigendecompositions).
