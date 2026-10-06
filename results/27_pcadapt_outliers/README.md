# 27 — PCAdapt F_ST outlier scan (selection-genomics)

**Script:** `scripts/27_pcadapt_outliers.py`

## What was done

PCAdapt-style selection-genomics scan (Luu, Bazin & Blum 2017; pcadapt R package) using the leading two PCs of the standardised dosage as the latent structure and computing per-marker Mahalanobis distance.

## Method

PCA on standardised dosage (per-marker centring, scale by √(2pq)). For each marker j, the K-vector of regression slopes z_j against the leading K = 2 PCs is computed. Mahalanobis distance from the multivariate mean using the inverse covariance of z follows χ²_K under neutral drift. Per-marker p-values → Benjamini-Hochberg q-values. Genomic-control λ_GC reported.

## Quick findings

- **No outlier reaches BH q < 0.10**. Top q = 0.25.
- The top-10 outlier candidates by Mahalanobis distance are concentrated on **Ss05** (6 of 10; the other four are Ss08×2, Ss09×1, Ss03×1). The distinct Ss05 peaks are Ss05:2.2 Mb, 9.7 Mb, 10.9 Mb, 20.4 Mb (two closely-linked SNPs), 51.0 Mb. Most extreme p = 6.8 × 10⁻⁴ (Ss05:20.4 Mb), BH q = 0.25.
- Convergence with DAPC LD1-loadings (`29_dapc/`) — same Ss05 positions are top-ranked by an independent method.

**Correction (2026-09-30):** the earlier "8 of 10 on Ss05" figure was wrong — the exact top-10 chromosome tally from `tables/pcadapt_outlier_top50.csv` is Ss05×6, Ss08×2, Ss09×1, Ss03×1. Manuscript body and Fig 7 caption corrected to "6 of the top-10".

## Fig 7A x-positioning bugfix (2026-09-30)

The Manhattan panel (`figures/fig64_pcadapt_manhattan.png/.pdf`) previously mis-placed points and gene call-outs onto the wrong chromosomes. The genome-wide x-coordinates were built in chromosome-grouped order but assigned back to the DataFrame in original row order (`anc_ok["x"] = xs`), misaligning every point to a row it did not belong to. Fixed by mapping each SNP's x from its own chromosome offset (`anc_ok["x"] = anc_ok["pos"] + anc_ok["chr_ayb"].map(offsets)`) and centring the x-tick labels on `(min+max)/2`. After the fix the strongest outliers (y ≈ 3.1) correctly cluster on Ss05, with Ss09 (ACR3_2) and Ss08 (HPGT2/FREE1) as the next tier — matching `tables/pcadapt_outlier_top50.csv`. Figure regenerated; `word/media/image15.png` in `A1_manuscript_final_submitted.docx` re-embedded with the corrected panel.

## Caveats

- PCAdapt assumes that drift / migration follow the PC structure; if structure is more complex (admixture, recent selection), it can miss signals.
- At K = 2 PCs and n = 95, the per-marker Mahalanobis stat has limited power; q ≈ 0.25 means the signal could be drift-only.
- The Ss05 convergence is suggestive of a coherent diversity / selection signal but not statistically significant. Validation requires either a larger panel or a direct test (e.g. iHS / nSL on haplotype data — beyond the scope of DArTseq).
