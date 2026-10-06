# 98 — Direct genotype concordance for the near-clonal trio

**Script:** `scripts/98_trio_duplicate_concordance.py`

Reviewer Item 6 (trio independence framing + duplicate statistics). The reviewer's substantive objection is that the three lines of support for the near-clonal trio (TSs151B / TSs358 / TSs361) are **not mutually independent**: the genomic relationship matrix (GRM) and the maximum-likelihood phylogeny are both computed from the same 1,625 SNPs, and the SilicoDArT layer, although a different marker chemistry, comes from the same DNA extractions, restriction digest, and sequencing libraries. The alternative hypothesis when calling accessions duplicates is therefore a **sample-level artefact** — contamination, a plate swap, or a mixed seed lot (the reviewer's TSs357 example).

## Why

Marker-layer agreement cannot distinguish a true genetic duplicate from a sample artefact, because a shared artefact propagates to every layer built from the same DNA. The statistic that *does* discriminate is the **direct genotype concordance** between the accessions:

- a true clonal/duplicate pair → near-100 % identical calls, and the few discordances are heterozygous-versus-homozygous (residual heterozygosity / calling noise), **not** opposite homozygotes;
- a mixed seed lot or a swap of a genuinely different genotype → a block of opposite-homozygote (0 ↔ 2) differences.

The reviewer also asked (duplicate-detection criteria) for G_ii (= 1 + F under VanRaden method 1), G_ij, and the inbreeding-normalized relatedness G_ij / √(G_ii·G_jj), because raw G_ij is inflated by (1 + F) in inbred lines. The reviewer further asked that the **same statistics be reported for the five "close-pedigree" pairs** (0.81–0.82 G_ij) so that readers can judge where the 0.85 duplicate threshold sits relative to the highly related but non-duplicate tail, and noted that the correct unrelated-pair expectation is **0 (off-diagonal mean −0.011)**, not the ≈ 0.25 the Methods text originally stated.

## Method

- Genotypes from `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv`, converted to 0/1/2 dosage, QC identical to script 01 (sample & marker call-rate ≥ 0.90, MAF ≥ 0.05) → **1,625 markers × 95 lines**.
- **Concordance** per pair: over markers called in both members, count identical calls (% identical), discordant calls, and split discordances into heterozygous-vs-homozygous vs opposite-homozygote (0 ↔ 2).
- **Relatedness** read from `results/13_grm_crosspairs/tables/grm_vanraden.csv`: G_ii, G_ij, and G_ij / √(G_ii·G_jj).
- Pairs: the three trio pairs, plus TSs357 (heterozygous, H_o ≈ 0.16, F_ROH ≈ 0.07) against each trio member as the reviewer's mixed-seed-lot control, plus the five **close-pedigree pairs** flagged in `results/56_dedup_sensitivity` (TSs361/TSs357, TSs156A/TSs363, TSs357/TSs358, TSs151B/TSs357, TSs60/TSs282) at 0.81–0.82 raw G_ij.
- **GRM off-diagonal mean** computed across all 95 × 95 off-diagonal entries as the empirical unrelated-pair baseline the reviewer cites.

## Findings

**Trio internal concordance (true-duplicate signature):**

| Pair | markers co-called | % identical | discordant | het-vs-hom | opposite-hom |
|---|---|---|---|---|---|
| TSs151B–TSs358 | 1,613 | **99.07** | 15 | 15 | 0 |
| TSs151B–TSs361 | 1,613 | **98.95** | 17 | 17 | 0 |
| TSs358–TSs361 | 1,611 | **99.57** | 7 | 6 | 1 |

The trio pairs share **98.9–99.6 %** of calls, and the 7–17 discordances per pair are almost entirely het-vs-homozygous with **at most one opposite-homozygote call** — the signature of a true genetic duplicate, not a mixed seed lot or plate swap.

**TSs357 control (reviewer's mixed-seed-lot hypothesis):** only **73.1–73.3 %** identical to the trio (430–434 discordant), normalized relatedness 0.68 — a close-pedigree relative, an order of magnitude less concordant than the trio members are with each other. TSs357 is therefore **not** a trio member.

**Relatedness (duplicate-detection criteria):**

| | G_ii | trio-pair G_ij | G_ij / √(G_ii·G_jj) |
|---|---|---|---|
| TSs151B | 1.614 | 1.600 (151B–358), 1.600 (151B–361) | 0.990, 0.989 |
| TSs358 | 1.618 | 1.610 (358–361) | 0.993 |
| TSs361 | 1.623 | — | — |
| TSs357 | 0.893 | 0.814–0.823 vs trio | 0.678–0.684 |

Even after normalizing for inbreeding, trio-pair relatedness is **0.99**, so the near-identity is not an artefact of the (1 + F) inflation of raw G_ij. TSs357's normalized relatedness to the trio (0.68) is far below the trio's internal 0.99.

**Five close-pedigree pairs (0.81–0.82 raw G_ij) — where the 0.85 threshold sits:**

| Pair | markers co-called | % identical | discordant | het-vs-hom | opposite-hom | G_ij | G_ii / G_jj | G_ij / √(G_ii·G_jj) |
|---|---|---|---|---|---|---|---|---|
| TSs361–TSs357 | 1,613 | 73.34 | 430 | 426 | 4 | 0.823 | 1.623 / 0.893 | 0.684 |
| TSs156A–TSs363 | 1,613 | 87.17 | 207 | 27 | **180** | 0.821 | 1.528 / 1.589 | 0.527 |
| TSs357–TSs358 | 1,613 | 73.34 | 430 | 425 | 5 | 0.818 | 0.893 / 1.618 | 0.681 |
| TSs151B–TSs357 | 1,615 | 73.13 | 434 | 430 | 4 | 0.814 | 1.614 / 0.893 | 0.678 |
| TSs60–TSs282 | 1,601 | 73.08 | 431 | 430 | 1 | 0.808 | 0.870 / 1.609 | 0.683 |

None of the five close-pedigree pairs approaches the trio's duplicate signature. Four of the five sit at **73 %** direct concordance and normalized relatedness **0.68**, an order of magnitude less concordant than the trio's 99 %. The fifth, TSs156A–TSs363, is more identical at the call level (87 %) but carries **180 opposite-homozygote discordances** (versus ≤ 1 for every trio pair) — the signature of two genuinely different genotypes sharing ancestry, not a duplicate; its normalized relatedness (0.53) is the lowest of the five. The gap between the highest close-pedigree normalized relatedness (0.68) and the lowest trio value (0.99) is wide and empty, so the 0.85 raw-G_ij threshold cleanly separates the near-clonal trio from the close-pedigree tail regardless of where in the 0.68–0.99 gap the cut is placed.

**GRM off-diagonal mean = −0.011** (diagonal mean G_ii = 1.008, n = 95). The unrelated-pair expectation is therefore ≈ 0, confirming the reviewer's correction of the original "≈ 0.25" Methods statement.

## Outputs

- `tables/trio_concordance.csv` — per-pair co-called markers, % identical, discordant, het/hom & opposite-hom split.
- `tables/trio_relatedness.csv` — G_ii, G_ij, and G_ij / √(G_ii·G_jj) for trio pairs + TSs357.
- `tables/close_pedigree_concordance.csv` — same concordance statistics for the five close-pedigree pairs.
- `tables/close_pedigree_relatedness.csv` — G_ii, G_ij, and G_ij / √(G_ii·G_jj) for the five close-pedigree pairs.
- `tables/grm_offdiag_summary.csv` — GRM off-diagonal mean, diagonal mean, n accessions.

## Reproduce

```bash
cd scripts && /Users/black_einstein/miniconda3/envs/yeast-viz/bin/python 98_trio_duplicate_concordance.py
```

## Caveats

- Concordance is on the post-QC 1,625-marker set; discordance counts scale with co-called markers (~1,613 per trio pair).
- TSs357's 73 % concordance is consistent with either a close-pedigree relative or a *partial* mixed seed lot; the trio-internal concordance (99 %+) is not.

## Tool versions

Python 3.9.23 (conda env `yeast-viz`); numpy, pandas.
