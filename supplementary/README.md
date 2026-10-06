# Supplementary tables and figures

Manuscript supplementary items for Narh-Madey et al., 2026, numbered as they
appear in the Supplementary Information (SI) document. The SI PDF is the
authoritative source for table and figure legends; this directory provides the
underlying CSV tables and the figure image files.

```
tables/     Supplementary Tables S1–S15 (CSV)
figures/    Supplementary Figures (PDF vector + PNG raster)
```

## Tables (`tables/`)

| File | Content |
|---|---|
| `TableS1_marker_qc.csv` | Per-marker QC (call rate, MAF, heterozygosity) |
| `TableS2_silicoDArT_marker_qc.csv` | SilicoDArT per-marker QC |
| `TableS2b_cross_marker_system_concordance.csv` | SNP vs SilicoDArT cross-system concordance |
| `TableS3_Ss05_outlier_attribution.csv` | Ss05 outlier-locus attribution |
| `TableS3b_private_alleles_C2.csv` | Private alleles, minor cluster (C2) |
| `TableS3c_private_alleles_C1.csv` | Private alleles, main cluster (C1) |
| `TableS3d_Ss05_outlier_attribution_unique.csv` | Ss05 attribution, unique-placement subset |
| `TableS3e_top_outliers_pcadapt_dapc.csv` | Top pcadapt / DAPC outlier loci |
| `TableS4_f_roh_per_sample.csv` | Per-sample F_ROH |
| `TableS5_core_collection_top20_stratified.csv` | Stratified top-20 core collection |
| `TableS6_Table_S6_reanchoring_per_marker_summary.csv` | Per-marker re-anchoring summary |
| `TableS7_per_accession_roh_class_spectrum.csv` | Per-accession ROH length-class spectrum |
| `TableS8_panel_summary.csv` | Panel-level summary statistics |
| `TableS9_dedup_delta_summary.csv` | Pre/post deduplication deltas |
| `TableS10_accession_list_qc_cluster.csv` | Accession list with QC + cluster label |
| `TableS11_marker_placement_uniqueness.csv` | BLAST placement-uniqueness classification |
| `TableS12_cross_shortlist_dedup_normalized.csv` | Inbreeding-normalized cross shortlist |
| `TableS13_duplicate_detection_criteria.csv` | Cryptic-duplicate detection criteria |
| `TableS14_phenotype_per_accession.csv` | Per-accession phenotype values |
| `TableS15_het_vs_depth_outbred_lines.csv` | Heterozygosity vs read depth, outbred lines |

Most of these tables are re-derivable from the per-stage outputs under
`results/`; they are collected here under manuscript numbering for direct
reference from the SI. The pipeline-named outputs in `results/<NN>_*/tables/`
retain their stage-local names and may use a different numbering.

## Figures (`figures/`)

Supplementary figures are provided as vector PDF and 300 dpi PNG. File names
carry the internal pipeline figure IDs (e.g. `fig28`, `fig67`); the SI document
maps each to its S-number and legend.
