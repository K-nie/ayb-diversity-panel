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

## Asymmetric-dropout bias on r — direction of the subset effect (Reviewer 3 point E, 2026-10-05)

Reviewer 3 asked whether the all-main-cluster dropout (17 of the 84 main-cluster accessions drop at the SilicoDArT sample-call-rate filter; all 11 minor-cluster accessions pass) biases the two-layer comparison, and whether the retained 78-accession subset is representative. The manuscript (Results §3.9) had described r = 0.787 as a "lower bound on cross-layer agreement" on the strength of the main-cluster depletion. A direct sensitivity test shows that framing is **directionally wrong** — the depletion *inflates* r, it does not deflate it.

Decomposing the 3,003 retained pairs by cluster-pair type (cross-layer Pearson r, SNP-layer IBS vs SilicoDArT-layer IBS):

| pair type | n pairs | cross-layer r | mean SNP IBS | mean SilicoDArT IBS |
|---|---|---|---|---|
| within-main  | 2,211 | 0.709 | 0.504 | 0.341 |
| between-cluster | 737 | 0.763 | 0.578 | 0.371 |
| within-minor | 55 | 0.974 | 0.340 | 0.208 |
| **pooled (retained)** | **3,003** | **0.787** | — | — |

The fully-retained minor cluster supplies the tightest cross-layer agreement (within-minor r = 0.974), while the 20 %-depleted main cluster supplies the loosest (within-main r = 0.709). The retained subset therefore over-represents the most concordant pair class and under-represents the least concordant one. Reweighting the retained pairs to the full-panel cluster mix (84 main + 11 minor) — which handles both the per-class r and the between-cluster range-extension effect on the pooled correlation — moves the pooled estimate **down from 0.787 to 0.773**. This is the *optimistic* bound: it assumes the 17 dropouts would behave like the retained main-cluster samples, whereas they are precisely the low-call-rate samples (SilicoDArT call rate 0.79–0.89 vs panel median 0.93), so their pairs would be noisier and pull r lower still.

An assumption-free subsampling test confirms the same direction without any per-class-r or reweighting assumption. Fixing the 11 minor-cluster accessions and randomly subsampling the main cluster (200 draws each), pooled r rises monotonically as the main cluster shrinks: n_main = 30 → r = 0.858; 40 → 0.831; 50 → 0.810; 60 → 0.795; 67 (all retained) → 0.787. The trend is driven by combinatorics: between-cluster pairs scale as n_main × 11 (linear) while within-main pairs scale as C(n_main, 2) (quadratic), so a smaller main cluster raises the *fraction* of high-IBS between-cluster pairs that extend the range and inflate the pooled correlation. Extrapolating toward the true n_main = 84 lands near r ≈ 0.767, below the retained 0.787 and consistent with the reweighting estimate (0.773) — and optimistic, since the real dropouts are the low-call-rate samples.

**Conclusion for the manuscript.** Two biases act in opposite directions and must not be conflated: (i) the dominant SilicoDArT chemistry collapses heterozygote/homozygote contrasts, compressing IBS and pushing r *below* the true structural agreement (a genuine downward pressure, from the marker system); (ii) the all-main-cluster dropout removes the lowest-concordance within-main and noisiest low-call-rate pairs and over-weights the tightly-concordant minor cluster, pushing r *upward* (an inflation, from the sample subset). The subset depletion therefore does **not** make r = 0.787 a conservative lower bound — if anything the subset inflates it. The honest statement is that r = 0.787 is an estimate carrying a modest upward bias from the subset composition and a downward bias from the dominant-marker chemistry, best read as corroboration of structure (which is robust: 100 % K = 2 concordance, trio recovered on both layers) rather than as a bound on panel-wide concordance.

Reproduce: load `scripts/54_silicoDArT_diversity.py` as a module (`sys.path.insert(0,"scripts")`), rebuild `calls_post` (78 × 1,156) and the matched SNP dosage, compute `pairwise_ibs` on both layers, tag each pair by the `results/01_qc_pca_power/tables/pca_coords.csv` cluster code (1 = main, 0 = minor), and (a) correlate within each pair-type class and (b) compute a full-panel-mix-weighted Pearson r. No new output files were written; this is an interpretation-level sensitivity check feeding the Results §3.9 wording.

## Script

`scripts/54_silicoDArT_diversity.py`. Runtime: ~ 4 minutes on a laptop (dominated by the SNP-layer dosage rebuild + the pairwise IBS computation, which is O(n² × m) at n = 78 / m = 1,156 to 1,625).
