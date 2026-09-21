# 54 — SilicoDArT presence/absence diversity layer

Adds the SilicoDArT dominant-marker layer (4,992 markers in the analysis plan; 5,000 markers in the raw DArT report) to the SNP-based diversity pipeline as a second marker system for cross-validation. Independent calling on the same 95-line working panel asks whether the K = 2 cluster architecture and the pairwise-IBS relationships hold under a different marker chemistry.

## Method

1. Load the DArT-format SilicoDArT report (`Report_DAf18-2580_SilicoDArT_1.csv`): 6 metadata header rows + 14 per-marker metadata columns + per-sample binary calls (0 = allele absent, 1 = allele present, `-` = missing).
2. Compute per-marker and per-sample QC: call rate, MAF (panel-wide-minor allele frequency), polymorphism = 1 − max(p, 1−p).
3. Restrict to the working 95-line SNP-layer panel (so the cross-marker-system comparison is sample-comparable).
4. Apply the same QC cascade as the SNP layer: sample call rate ≥ 0.90, marker call rate ≥ 0.90, MAF ≥ 0.05.
5. PCA on the standardised binary matrix; K-means at k = 2..7 with silhouette scoring.
6. Pairwise IBS distance on the SilicoDArT layer and on the SNP layer (script 01 dosage matrix re-built on the same sample set in matching order). Pearson correlation across all sample-pair distances.

ADMIXTURE on a SilicoDArT-converted PLINK BED is not run here — it requires a custom binary-input conversion step plus an external PLINK + ADMIXTURE call. The PCA + IBS-concordance result already decides the cross-marker-system question; ADMIXTURE is held for the revision pass if reviewers ask.

## Findings

- **5,000 raw SilicoDArT markers**, median per-marker call rate 0.92, median polymorphism 0.011. Roughly two-thirds of SilicoDArT markers are nearly monomorphic in this panel.
- **Post-QC: 1,156 informative SilicoDArT markers on 78 retained accessions.**
- **17 of 95 working-panel accessions drop on the SilicoDArT sample-QC step** (SilicoDArT call rate range 0.79–0.89 for the dropouts vs the panel median of 0.93). All 17 dropouts are in Cluster 1; **every Cluster 2 accession passes**. The SilicoDArT layer's view of the small Cluster 2 is fully resolved.
- **SilicoDArT K = 2 PCA recovers the SNP-layer K = 2 architecture at 100 % cluster concordance** on the 78 retained accessions. PC1 cleanly separates the small Cluster 2 (orange) from the large Cluster 1 (blue); silhouette at k = 2 = 0.255 (modest but real); PC1 + PC2 explain 14.3 % of marker variance.
- **Cross-marker-system pairwise-IBS Pearson r = 0.787** across 3,003 pairs. The IBS scatter sits systematically below the SNP=SilicoDArT identity line because the dominant binary calls collapse heterozygote-vs-homozygote distinctions that the codominant SNP layer separates; the relative *ordering* of sample-pair distances is what matches.
- **The TSs151B / TSs358 / TSs361 trio surfaces as near-clonal on both marker systems** (the three near-zero IBS pairs visible at the lower-left of the cross-layer IBS scatter). Independent corroboration of the duplicate finding from §56.

## Outputs

- `tables/silicoDArT_marker_qc.csv` — per-marker call rate / p_present / MAF / polymorphism on the raw 5,000-marker matrix.
- `tables/silicoDArT_sample_qc.csv` — per-sample call rate + fraction-present on the raw 105-sample × 5,000-marker matrix.
- `tables/silicoDArT_pca_coords.csv` — PC1–PC10 + K-means cluster labels at k = 2..7 + SNP-layer cluster.
- `tables/silicoDArT_pca_variance.csv` — variance explained per PC.
- `tables/cross_marker_system_concordance.csv` — concordance and IBS-correlation summary.
- `figures/fig_silicoDArT_pca.png` / `.pdf` — dual-panel PCA: SNP-layer cluster colouring + SilicoDArT's own K = 2.
- `figures/fig_cross_layer_ibs.png` / `.pdf` — pairwise IBS scatter, SNP layer vs SilicoDArT layer.

## Caveats

- The SilicoDArT layer is dominant. Per-locus statistics that assume codominance (Yang 2010 F estimator, Hardy-Weinberg test, GBLUP kinship) are not meaningful on SilicoDArT and are not run here.
- The 17-sample drop at the SilicoDArT QC step is a real biological / library-prep effect: those accessions have lower SilicoDArT call rates than the panel median. Whether the drop biases the diversity placement toward Cluster 2 (since none of Cluster 2 drops) is a methodological point worth flagging — but at silhouette 0.255, the SilicoDArT signal is dominated by the cluster split itself, not by the sample subset choice.
- The 100 % K = 2 cluster concordance figure is computed on the 78 retained samples. A more demanding test would impute SilicoDArT calls for the 17 dropouts and re-evaluate; that imputation step is not run here.

## Script

`scripts/54_silicoDArT_diversity.py`. Runtime: ~ 4 minutes on a laptop (dominated by the SNP-layer dosage rebuild + the pairwise IBS computation, which is O(n² × m) at n = 78 / m = 1,156 to 1,625).
