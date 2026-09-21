# Chromosome-anchored DArTseq genotyping reveals genetic architecture and cryptic duplicates in an African yam bean diversity panel

Analysis code and result tables/figures for a 95-line *Sphenostylis stenocarpa* (African yam bean) DArTseq diversity panel re-anchored to the chromosome-scale *S. stenocarpa* reference genome. Covers QC and PCA, ADMIXTURE, AMOVA, F_ST, genome-wide and per-chromosome LD decay, maximum-likelihood phylogenies (IQ-TREE 2 and RAxML-NG), pcadapt and DAPC selection scans, the VanRaden genomic relationship matrix with cryptic-duplicate detection, F_IS and runs of homozygosity, SilicoDArT cross-marker-system validation, and an inbreeding-normalized cross shortlist.

This repository contains the **analysis code and derived result tables/figures** for
the study. Raw genotype and phenotype data are archived separately (see Data below);
manuscript drafts are not included.

## Repository layout

```
scripts/    numbered Python and R analysis scripts (shared helpers: _plotstyle, _pheno, _figstyle)
results/    one directory per analysis stage, each with figures/, tables/, and a README.md
            documenting method, inputs, outputs, findings, and caveats
refs/       machine-learning best-practice reference material
```

Each `results/<NN>_*/README.md` is the reproducibility and methods record for that
stage: what it does and why, its inputs, the exact commands and thresholds, the outputs,
and the caveats.

## Reproducing

Scripts run from a single conda environment (`yeast-viz`: numpy, pandas, scikit-learn,
matplotlib, plus the R toolchain for the `.R` steps). Stages are numbered by pipeline order;
run the lower-numbered QC/PCA stage first, then the analysis of interest. Per-stage READMEs
give the exact command and parameters.

## Data

- **Raw DArTseq genotypes** (order DAf18-2580, 105 accessions; post-QC 95 accessions x 1,625
  markers at sample and marker call rate >= 0.90, MAF >= 0.05): Zenodo [10.5281/zenodo.20348832](https://doi.org/10.5281/zenodo.20348832).
- **Reference genome**: *S. stenocarpa* chromosome-scale assembly, ENA [PRJEB57813](https://www.ebi.ac.uk/ena/browser/view/PRJEB57813)
  (Shorinola et al. 2024); Funannotate annotation at Zenodo [10.5281/zenodo.13853757](https://doi.org/10.5281/zenodo.13853757).
- Phenotype data are held by the breeding program and available from the authors on
  reasonable request.

## Headline findings

- 89.3% of DArTseq markers re-anchor directly to the *S. stenocarpa* reference (11 pseudo-chromosomes, 649.8 Mb), up from 11.9% under the earlier cowpea proxy.
- A near-clonal trio (TSs151B, TSs358, TSs361; VanRaden G_ij approximately 1.6) is detected and collapsed before downstream diversity and breeding-shortlist steps.
- Genome-wide LD half-decay approximately 102 kb with an LD-based effective population size N_e approximately 476, refit on the AYB-anchored marker set.
- F_ROH approximately 0.256 reconciles the otherwise low F_IS, consistent with a predominantly selfing mating system.
- pcadapt and DAPC scans converge on pseudo-chromosome Ss05; per-locus Hudson F_ST is reported across the 11 pseudo-chromosomes.

## Citation

If you use this code or its outputs, please cite this repository (see `CITATION.cff`)
and the accompanying manuscript.

## License

Code is released under the MIT License (see `LICENSE`).
