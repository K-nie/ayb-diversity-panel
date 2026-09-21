# 34 — Runs of homozygosity (ROH) + F_ROH

**Script:** `scripts/34_roh.py`

## What was done

Per-accession runs of homozygosity from the AYB-anchored dosage matrix
and a per-accession inbreeding coefficient F_ROH computed at four
ROH-length cutoffs.

## Method

ROH segments were called with a PLINK-style sliding 20-SNP window
(Howrigan et al. 2011) allowing at most one heterozygous and two
missing genotypes per window. SNPs flagged homozygous in ≥ 5 % of
containing windows were merged into ROH segments. Segments ≥ 500 kb
and ≥ 8 SNPs were retained. F_ROH was reported per accession at four
length cutoffs (≥ 0.5, 1.0, 2.0, 5.0 Mb) as the fraction of the
AYB-anchored genome span covered by ROH at each threshold.

## Inputs

- HapMap dosage matrix (post-QC, MAF ≥ 0.05)
- AYB marker anchoring table

## Outputs

- `tables/roh_segments.csv` — sample, chrom, start_bp, end_bp, n_snps,
  length_bp
- `tables/f_roh_per_sample.csv` — sample × four F_ROH cutoffs
- `figures/fig73_roh_length_distribution.png/.pdf` — log-binned
  ROH-length histogram across all segments
- `figures/fig74_f_roh_per_sample.png/.pdf` — per-accession F_ROH at the
  ≥ 0.5 Mb cutoff, sorted

## Quick findings

- **Mean F_ROH (≥ 0.5 Mb) = 0.256** — substantial autozygosity, consistent
  with the AYB selfing-rate literature.
- Resolves the global F_IS = 0.011 anomaly from script 22: the
  heterozygosity-deficit estimator misses long-stretch autozygosity at
  this marker density (~500 kb mean spacing), while ROH catches it.
- Per-accession F_ROH ranges from ~0.05 (least inbred) to ~0.60 (most
  inbred); the spread should be a routine QC flag in genebank
  regeneration cycles.

## Caveats

- DArTseq density limits ROH-call sensitivity to segments ≥ 500 kb;
  shorter autozygous tracts (typical of recent crosses) are not detected.
- Het-call inflation at low DArTseq coverage works in the opposite
  direction from autozygosity-call inflation — both biases cancel only
  approximately. Treat F_ROH as a relative-rank ordering, not an
  absolute fraction.
