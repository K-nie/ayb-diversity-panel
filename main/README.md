# Main-text tables and figures

Manuscript main-text tables and figures for Narh-Madey et al., 2026, numbered as
they appear in the paper body. The manuscript is the authoritative source for
table and figure legends; this directory provides the underlying CSV tables and
the figure image files.

```
tables/     Main-text tables (CSV)
figures/    Main-text figures (PDF vector + PNG raster)
```

## Tables (`tables/`)

| File | Content |
|---|---|
| `Table1_amova_table.csv` | AMOVA variance partition |
| `Table2_cophenetic_correlation_matrix.csv` | Cophenetic correlation across tree methods |
| `Table3_duplicate_components.csv` | Cryptic-duplicate connected components |
| `Table5_per_panel_diversity_stats.csv` | Per-panel diversity statistics |
| `Table6_cross_panel_qc_comparison.csv` | Cross-panel QC comparison |

The pipeline-named outputs under `results/<NN>_*/tables/` retain their
stage-local names and may use a different numbering; the files here follow the
manuscript's main-text table numbering.

## Figures (`figures/`)

Main-text figures are provided as vector PDF and 300 dpi PNG. File names carry
the internal pipeline figure IDs (e.g. `fig01`, `fig68`); the manuscript maps
each to its figure number and legend.
