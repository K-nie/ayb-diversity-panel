# 27 — PCAdapt F_ST outlier scan (selection-genomics)

**Script:** `scripts/27_pcadapt_outliers.py`

## What was done

PCAdapt-style selection-genomics scan (Luu, Bazin & Blum 2017; pcadapt R package) using the leading two PCs of the standardised dosage as the latent structure and computing per-marker Mahalanobis distance.

## Method

PCA on standardised dosage (per-marker centring, scale by √(2pq)). For each marker j, the K-vector of regression slopes z_j against the leading K = 2 PCs is computed. Mahalanobis distance from the multivariate mean using the inverse covariance of z follows χ²_K under neutral drift. Per-marker p-values → Benjamini-Hochberg q-values. Genomic-control λ_GC reported.

## Quick findings

- **No outlier reaches BH q < 0.10**. Top q = 0.25.
- The top-10 outlier candidates by p are concentrated on **Ss05** (8 of 10): Ss05:9.7 Mb, 10.9 Mb, 20.4 Mb, 51.0 Mb, plus several closely-linked SNPs.
- Convergence with DAPC LD1-loadings (`29_dapc/`) — same Ss05 positions are top-ranked by an independent method.

## Caveats

- PCAdapt assumes that drift / migration follow the PC structure; if structure is more complex (admixture, recent selection), it can miss signals.
- At K = 2 PCs and n = 95, the per-marker Mahalanobis stat has limited power; q ≈ 0.25 means the signal could be drift-only.
- The Ss05 convergence is suggestive of a coherent diversity / selection signal but not statistically significant. Validation requires either a larger panel or a direct test (e.g. iHS / nSL on haplotype data — beyond the scope of DArTseq).
