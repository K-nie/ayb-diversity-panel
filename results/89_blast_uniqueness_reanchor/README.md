# 89 — BLAST placement uniqueness and paralogy → heterozygosity test

Re-derives the *S. stenocarpa* anchoring uniqueness statistics from a BLAST run that reports all qualifying hits per tag, replacing the original `-max_target_seqs 1` run that cannot count placement multiplicity. Directly addresses Reviewer 1 major point 3 (and feeds Reviewer 2 point B on the heterozygosity anomaly).

## Why

The submitted Methods reported megablast with `-max_target_seqs 1` and claimed "2,862 of 3,204 markers (89.3 %) mapped to a single position." `-max_target_seqs 1` returns *one* hit per query regardless of how many equally good hits exist and is explicitly not a "best hit" or "unique hit" filter (Shah et al. 2019, Bioinformatics 35:1613; NCBI BLAST user guidance). A "single position" claim therefore cannot be made from that run. This analysis re-runs the search reporting up to 50 hits per tag under the identical significance thresholds and counts genomic placement multiplicity per marker.

## Inputs

- `refs/ayb_genome/dartseq_tags.fasta` — 6,200 DArT tag sequences (reference + SNP-allele tag for each marker).
- `refs/ayb_genome/Sphenostylis_stenocarpa_chrom.fasta` — chromosome-scale assembly (Waweru et al. 2024, *Sci. Data* 11:1384 = manuscript ref [7]; ENA OY731398–OY731408; 11 pseudo-chromosomes; **649,801,261 bp** confirmed from the rebuilt BLAST DB, matching the 649.8 Mb reported in Methods). The 11 pseudo-molecules sum to 649.80 Mb; adding the 8,422 unscaffolded contigs (`AYB_PhaseGenomics_8422_unscaffolded_contigs.fasta`, 51.76 Mb) gives a full deposited assembly of 701.56 Mb (≈ the 701.3 Mb genome size reported by Waweru et al.). Only the 11 pseudo-molecules were used for positional anchoring.
- `refs/ayb_genome/ayb_marker_anchoring.csv` — the 3,204 study SNP-markers (`rs`); these collapse to **3,100 unique DArT clones** (104 clones carry >1 SNP).
- `results/01_qc_pca_power/tables/marker_qc.csv` — per-marker observed heterozygosity (`het`).

## Commands (reproduce)

BLAST DB and tag paths must be free of spaces (BLAST+ splits paths on whitespace), so the DB was built under `/tmp/ayb_blast/` via symlinks:

```bash
# BLAST+ 2.17.0 (submitted Methods cited 2.13; version updated in revision)
makeblastdb -in ayb_chrom.fasta -dbtype nucl -out ayb_chrom -title ayb_chrom

blastn -task megablast -query dartseq_tags.fasta -db ayb_chrom \
  -perc_identity 95 -evalue 1e-10 -max_target_seqs 50 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads 4 -out dartseq_to_ayb.multihit.tsv

python scripts/89_blast_uniqueness_reanchor.py
```

Thresholds `-perc_identity 95 -evalue 1e-10` are identical to the submitted anchoring run; only `-max_target_seqs 1` is replaced by `-max_target_seqs 50` so multiplicity can be counted.

## Method

- A tag's qualifying hits are those with percent identity ≥ 95 and e-value ≤ 1×10⁻¹⁰.
- Placement is a property of the **clone/tag**, so each clone is classified once and the class is propagated to every SNP-marker on it (keeps the manuscript's 3,204 SNP-marker denominator while counting loci at clone resolution).
- Hits on the same chromosome within **1 kb** collapse to one genomic locus (a 67 bp tag cannot resolve finer, and a marker's reference and SNP-allele tags land at the same spot).
- Classification per marker: **uniquely placed** = exactly one locus; **multi-mapping** = ≥ 2 loci; **unplaced** = 0 qualifying hits.
- For multi-mappers, the best-vs-second-locus bitscore gap distinguishes genuinely ambiguous placements (gap < 2 bits) from those with an unambiguous best locus despite a weaker secondary hit (gap ≥ 2 bits).
- Paralogy → heterozygosity test: observed heterozygosity of multi-mapping vs uniquely placed markers, one-sided Mann-Whitney U (H₁: multi > unique). Collapsed paralogues read as constitutive heterozygotes, so if the panel's H_obs excess were a mapping artefact, multi-mappers should carry inflated het.

## Findings

**Placement uniqueness (n = 3,204 SNP-markers):**

| Class | Markers | % of all |
|---|---|---|
| Uniquely placed (1 locus) | 2,603 | 81.2 % |
| Multi-mapping (≥ 2 loci) | 488 | 15.2 % |
| — of which genuinely ambiguous (gap < 2 bits) | 200 | 6.2 % |
| — of which unambiguous best locus (gap ≥ 2 bits) | 288 | 9.0 % |
| Unplaced (0 hits) | 113 | 3.5 % |
| **Confidently anchored (unique OR clear best locus)** | **2,891** | **90.2 %** |

- **The "89.3 % mapped to a single position" claim is not supported as worded.** Only **81.2 %** have a strictly unique qualifying hit. The ~89–90 % figure survives only under the weaker, correct claim that **90.2 % are confidently anchored to a unique *best* locus** (single locus, or a clear best locus with a ≥ 2-bit margin). **6.2 %** are genuinely ambiguous (≥ 2 near-equal loci) and 3.5 % are unplaced. The revision should replace "mapped to a single position" with "confidently anchored to a unique best locus (90.2 %; 81.2 % with a strictly unique hit, 6.2 % ambiguous, 3.5 % unplaced)."
- Clone-level view: 2,533 unique / 460 multi / 107 unplaced of 3,100 clones — same picture.

**Paralogy does NOT explain the heterozygosity anomaly:**

- Multi-mapping markers are **not** heterozygosity-enriched. Median het: uniquely placed 0.114 vs multi-mapping 0.048; means 0.179 vs 0.132. One-sided Mann-Whitney U (multi > unique) p = 1.000 — multi-mappers are, if anything, het-*depleted* (rank-biserial −0.23).
- Collapsed paralogy is therefore rejected as the source of the panel's observed heterozygosity; the H_obs pattern (Reviewer 2 point B) must be explained by residual outcrossing / seed-lot admixture / genotyping, not by mapping artefact. This is a clean negative result that removes one candidate artefact.

## Outputs

- `tables/marker_placement_classification.csv` — per SNP-marker: clone id, n_loci, placement class, top locus (chrom/pos/bitscore), best-second bitscore gap, observed het. **This is the per-marker anchoring table to deposit as a Data Sheet.**
- `tables/placement_summary.csv` — headline counts and the het-test p-value.

## Caveats

- Percent-identity ≥ 95 over a 67 bp tag tolerates up to ~3 mismatches; a small fraction of "multi-mapping" calls will be recent segmental duplications rather than true paralogues, which is why the bitscore-gap split (ambiguous vs clear-best) is reported alongside the raw multi count.
- The 1 kb locus-merge window is deliberately generous for 67 bp tags; tightening it would raise the multi count slightly (more within-region HSPs counted separately), so 81.2 % unique is a mild upper bound on ambiguity, not a lower bound.
- BLAST+ 2.17.0 was used here vs 2.13 in the submitted Methods; megablast placements are stable across these minor versions, but the version bump is noted for the revision.

## Script

`scripts/89_blast_uniqueness_reanchor.py`. Env: `yeast-viz` (numpy, pandas, scipy). Runtime: megablast ~5 s (4 threads); classification < 5 s.
