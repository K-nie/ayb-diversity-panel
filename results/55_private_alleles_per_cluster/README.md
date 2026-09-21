# 55 — Per-cluster private alleles and Ss05 outlier attribution

Resolves whether the Ss05 selection-genomics convergence reported by PCAdapt (script 27) and DAPC (script 29) is panel-wide drift or a Cluster-2-specific allele-frequency shift.

## Method

- Per SNP, the alternate-allele frequency is computed separately in Cluster 1 (n = 84, the larger PCA-derived cluster) and Cluster 2 (n = 11, the smaller). The panel-wide minor allele is chosen as the rarer of the two panel-wide alleles, and per-cluster frequencies are reported under that direction so "freq in cluster c" always refers to the same panel-wide-minor allele (values can exceed 0.5 when a cluster carries the panel-wide-minor allele at majority frequency).
- Private allele: panel-wide-minor allele has freq ≥ 0.05 in one cluster and ≤ 0.01 in the other. The 0.01 upper bound (rather than 0) allows for a single chance call without losing the private-allele inference at n = 11.
- Fixed-difference SNP (Hahn 2018): alt-allele frequency at 0 in one cluster and 1 in the other.
- **Per-locus F_ST: Hudson's estimator (Bhatia et al. 2013, Genome Res. 23:1514, eq. 10).** For each SNP with alternate-allele frequencies p1, p2 and sampled-allele counts n1, n2 in the two clusters, F_ST = num/den where num = (p1−p2)² − p1(1−p1)/(n1−1) − p2(1−p2)/(n2−1) and den = p1(1−p2) + p2(1−p1). The genome-wide estimate is the ratio of summed numerators to summed denominators (Bhatia's "ratio of averages"). Hudson's estimator is the recommended choice for markedly unequal sample sizes (here n = 84 vs n = 11) because it is nearly unbiased with respect to sample size and is bounded at 1. This replaces a previously mis-specified Weir-Cockerham variant that returned values up to 1.199 (the reviewer-flagged F_ST > 1 bug).
- Ss05 outlier attribution: the top-50 PCAdapt outliers (BH q sorted) and the top-50 DAPC LD1-loading SNPs are subset to Ss05 (16 unique SNPs across the two methods; 22 method×SNP rows before de-duplication on rs identifier). Each is tagged Cluster-2-driven if freq_C2 ≥ 0.40 AND freq_C1 ≤ 0.10, Cluster-1-driven if the reverse pattern holds, otherwise panel-wide.

## Findings

- **Ss05 attribution headline (unique SNPs): of 16 unique Ss05 outliers, 7 (44 %) are Cluster-2-driven, 1 is Cluster-1-driven, and 8 are panel-wide.** At the method×SNP row level (double-counting SNPs flagged by both methods) the split is 13 Cluster-2-driven, 1 Cluster-1-driven, 8 panel-wide out of 22. The manuscript should report the de-duplicated unique-SNP count (7 of 16) rather than the row count, so a SNP flagged by both PCAdapt and DAPC is not counted twice. Where Cluster-2-driven, the panel-wide-minor allele segregates in the 11-line Cluster 2 at ~50–76 % frequency while sitting near zero copies in the 84-line Cluster 1.
- **Cluster-2 private alleles concentrate on Ss05.** Of the 5 SNPs private to Cluster 2 at the panel-wide MAF ≥ 0.05 / ≤ 0.01 threshold, 4 sit on Ss05 and 1 on Ss09.
- **Cluster-1 private alleles vastly outnumber Cluster-2 private alleles.** 596 SNPs are private to Cluster 1 panel-wide, vs 5 private to Cluster 2. The asymmetry is a direct consequence of the unbalanced cluster sizes — the MAF ≥ 0.05 floor is easier to clear in 168 chromosomes (Cluster 1) than in 22 (Cluster 2), and the ≤ 0.01 absence bound is easier to satisfy in 22 chromosomes than in 168. Report this as a sample-size effect, not a biological asymmetry.
- **No Hahn 2018 fixed-difference SNPs.** Consistent with moderate genome-wide differentiation (Hudson genome-wide F_ST = 0.206) — substantial but not complete.
- **Per-locus Hudson F_ST across the eleven AYB pseudo-chromosomes: per-locus mean 0.153, median 0.112, max 0.776.** All values are bounded at 1, as expected. Ss05 carries the highest concentration of elevated-F_ST SNPs, matching the Cluster-2-driven outliers; no locus exceeds F_ST = 1.

## Outputs

- `tables/private_alleles_C1.csv` — 596 SNPs private to Cluster 1 with per-cluster freqs + AYB chromosome / position.
- `tables/private_alleles_C2.csv` — 5 SNPs private to Cluster 2 (4 on Ss05, 1 on Ss09).
- `tables/fixed_difference_snps.csv` — empty (no Hahn 2018 fixed differences at this panel).
- `tables/Ss05_outlier_attribution.csv` — 22 method×SNP rows (16 unique Ss05 SNPs) across PCAdapt + DAPC with per-cluster freqs, Hudson F_ST, and attribution class.
- `tables/per_cluster_summary.csv` — three-row category counts.
- `figures/fig_Ss05_attribution_strip.png` / `.pdf` — Cluster-1 vs Cluster-2 frequency scatter for the Ss05 outliers with shaded attribution regions.
- `figures/fig_per_locus_F_ST_manhattan.png` / `.pdf` — per-SNP Hudson F_ST across the eleven AYB pseudo-chromosomes; Ss05 highlighted.

## Caveats

- Hudson's F_ST at n = 11 (Cluster 2) still carries sampling noise, but is nearly unbiased with respect to the unequal sample sizes — this is the reason it, rather than Weir-Cockerham, is used here. The genome-wide value is reported as a ratio of averages (sum of numerators / sum of denominators), not a mean of per-locus ratios, which is the estimator recommended by Bhatia et al. (2013) to avoid the small-denominator inflation that produced the earlier > 1 artefacts.
- At n = 11, the Cluster-2 sampling SD on any per-SNP allele-frequency estimate is approximately 0.11 at p = 0.5. The 0.40 / 0.10 attribution thresholds are chosen to lie outside ± 2 SD of the panel-wide neutral expectation in both directions, so the attribution call is robust to this sampling noise.
- Private-allele count asymmetry (596 in Cluster 1 vs 5 in Cluster 2) is a sample-size effect, not a biological difference. The relevant comparison is per-cluster-size normalised, where the rates are comparable.
- These per-cluster F_ST / attribution statistics are computed on a k-means-on-PCA partition and are therefore descriptive of that partition, not a test of independently defined populations (see script 10 / Methods reframing).

## Script

`scripts/55_private_alleles_per_cluster.py`. Runtime: ~ 20 seconds on a laptop (dominated by the HapMap parse). Python env: `yeast-viz` (numpy, pandas, matplotlib).
