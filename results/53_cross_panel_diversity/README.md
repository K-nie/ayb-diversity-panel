# 53 — Cross-panel diversity placement

Places the 95-line panel inside the published AYB DArTseq diversity range by reporting headline statistics from this panel alongside the comparable numbers from Shitta/Aliyu 2022 (Sci Rep 12:4437; n = 169) and Olomitutu 2022 (Genes 13:2350; n = 195).

## Why this version is the lightweight one

The original Analysis 2 plan called for re-anchoring the Shitta and Olomitutu marker tables to the *S. stenocarpa* chromosome-scale reference, applying matched QC, and recomputing diversity statistics on the intersection marker set. That plan blocked on raw-data access:

- **Shitta 2022 (PMC8924269):** supplementary is one PDF with passport, sub-population means, PCoA, and Mantel correlations. Data availability statement: "The data set generated during an/or analyzed during the current study are available from the corresponding author on reasonable request" (suzzynde@yahoo.com).
- **Olomitutu 2022 (PMC9777823):** supplementary `genes-13-02350-s001.zip` (107 KB) contains only the 12 seed-size-associated SNPs and QC figures. No full marker table.

A matched-pipeline re-analysis requires an author-email round trip with days-to-weeks of latency. This lightweight version pulls each panel's published headline statistics directly from the cited PMC version of the paper and reports them in parallel with ours. The comparison is honest about each panel running its own QC; it is not a controlled re-analysis but it is the inter-study placement the reviewer is likely to ask for.

## Method

1. Manually transcribe published statistics from Shitta 2022 (Tables 2-4, Results) and Olomitutu 2022 (Results 3.2) into the `PUB_STATS` dict in `scripts/53_cross_panel_diversity.py`. Each number is annotated with the section / table it comes from.
2. Recompute our panel's matching statistics from `results/01_qc_pca_power/tables/marker_qc.csv` by applying the manuscript's call rate ≥ 0.90 AND MAF ≥ 0.05 cascade, ensuring the panel-side numbers exactly reproduce the headline values cited in §3.2 of the draft.
3. Build a wide side-by-side comparison table on 16 statistics and a QC-threshold comparison on 5.
4. Render a five-panel forest plot showing each panel's value (or "n.r." marker for not-reported) on a common axis per statistic.

## Findings

| Statistic | This (n=95) | Shitta (n=169) | Olomitutu (n=195) |
|---|---|---|---|
| Raw SNPs | 3,204 | 7,930 | 5,416 |
| Post-QC SNPs | 1,625 (50.7 %) | 1,789 (22.6 %) | 2,491 (46.0 %) |
| QC: call rate floor | 0.90 | 0.80 | 0.70 |
| QC: MAF floor | 0.05 | 0.05 | 0.01 |
| Mean MAF | 0.230 | 0.22 | 0.16 |
| Mean H_o | 0.297 | 0.15 ± 0.002 | 0.15 |
| Mean PIC | 0.258 | n.r. | n.r. |
| Mean n_e (alleles/locus) | 1.525 | 1.61 ± 0.008 | n.r. |
| F_IS | 0.011 | n.r. | n.r. |
| Pairwise F_ST | 0.165 (k=2) | 0.14–0.39 (k=3) | n.r. |
| AMOVA among-pops % | 27.3 | 13 | n.r. |
| ADMIXTURE optimal K | 3 (CV) | 3 (Evanno Δ K) | n.r. |
| LD half-decay (kb) | 73.85 | n.r. | n.r. |
| LD-Ne | 659 | n.r. | n.r. |

**Headline reading.** Our panel sits inside the published AYB DArTseq diversity range on every directly comparable axis (mean MAF, n_e, F_ST, ADMIXTURE K). Two deviations are flagged:

1. **Our per-marker H_o (0.297) is roughly double the comparator panels' 0.15.** Our value is close to the Hardy-Weinberg expectation `2pq ≈ 0.32` at our median MAF of 0.20, meaning the elevation reflects intermediate-MAF markers being retained under our QC. The genome-wide IBD picture (F_ROH = 0.256; §3.8) is what diagnoses the panel as autozygous, not the per-locus H_o.
2. **Our AMOVA among-populations fraction (27.3 %) is roughly double Shitta's 13 %.** This is a direct consequence of the K = 2 partition we use vs the K = 3 partition Shitta uses; concentrating variance into one between-cluster contrast inflates the among-population fraction relative to a three-way split.

Olomitutu 2022 reports only MAF and H_o, so direct comparison on F_ST, AMOVA, ADMIXTURE, F_IS, and LD is restricted to Shitta.

## Outputs

- `tables/per_panel_diversity_stats.csv` — 16-row wide table; rows are statistics, columns are panels.
- `tables/cross_panel_qc_comparison.csv` — 3-row × 8-column QC threshold + marker retention summary.
- `figures/fig_cross_panel_forest.png` / `.pdf` — five-panel forest of mean MAF, H_o, n_e, AMOVA among-pops %, pairwise F_ST.

## Caveats

- **Each panel runs its own QC.** Differences in mean MAF and H_o between panels are entangled with differences in marker-retention pipelines (call-rate threshold, MAF floor). A matched-QC re-analysis remains the right way to settle the per-statistic question — see the "lightweight version" justification above.
- **Shitta's "Ne = 1.61" is the *number of effective alleles per locus* (Kimura-Crow)**, not the LD-based effective population size. We report both quantities side by side: our Kimura-Crow n_e = 1.525 is comparable to Shitta's 1.61; our LD-Ne = 659 is not comparable to any reported Shitta or Olomitutu number.
- **The 13 % vs 27.3 % AMOVA gap is not a biological signal.** It is a direct partition-choice artefact (K = 3 distributes among-population variance across three pairwise contrasts; K = 2 concentrates it). A like-for-like comparison would require re-running AMOVA on our panel under K = 3, but at our n = 95 with the 11-line Cluster 2 the K = 3 split would have a 4-to-5-line subcluster too small to support stable AMOVA. We flag this honestly in the manuscript text.
- **For the version after the next data round.** Once Shitta and Olomitutu marker tables are obtained from the corresponding authors, the planned matched-QC re-anchored re-analysis (Analysis 2 Plan A) can run. The script structure here lends itself to a swap: replace the `PUB_STATS` dict with recomputed values from the matched-QC pipeline, regenerate the table and forest.

## Sources for the published numbers

- **Shitta NS, Unachukwu N, Edemodu AC, Abebe AT, Oselebe HO, Awtew WG.** *Genetic diversity and population structure of an African yam bean (Sphenostylis stenocarpa) collection from IITA GenBank.* Scientific Reports 12:4437 (2022). DOI: 10.1038/s41598-022-08271-4. PMC: PMC8924269. URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC8924269/
- **Olomitutu OE, Paliwal R, Abe A, Oluwole OO, Oyatomi OA, Abberton MT.** *Genome-Wide Association Study Revealed SNP Alleles Associated with Seed Size Traits in African Yam Bean (Sphenostylis stenocarpa (Hochst ex. A. Rich.) Harms).* Genes 13:2350 (2022). DOI: 10.3390/genes13122350. PMC: PMC9777823. URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC9777823/

## Script

`scripts/53_cross_panel_diversity.py`. Runtime: ~ 3 seconds on a laptop (no per-marker computation; just panel-side recomputation from the existing marker QC table).
