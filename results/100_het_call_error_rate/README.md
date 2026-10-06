# 100 — Empirical heterozygous-call error rate (Reviewer B.8)

Reviewer B.8 asked us to estimate an **empirical heterozygous-call (het-call) error rate** and to discuss its effect on the main-cluster F_ROH, given that a zero-heterozygosity ROH re-call lowers F_ROH relative to the 1-het-tolerant call (§3.6; `results/90_roh_sensitivity_null`). This note supplies that estimate from two internal, independent estimators and interprets the main-cluster F_ROH accordingly.

## Why

An autozygous genome should carry no true heterozygous sites, so any het call in a genuinely autozygous context is either a genotyping error or, at most, rare residual true heterozygosity — i.e. an **upper bound** on the error rate. Two such contexts exist in this panel:

1. **The near-clonal trio** (TSs151B / TSs358 / TSs361), whose members are genome-wide autozygous (F_ROH → 1.0) and mutually 99 % identical. At any marker the trio's true genotype is homozygous, so a het call in a member is a miscall. This is the cleanest, purest estimator (error-only).
2. **Heterozygous calls inside called ROH segments**, pooled across the minor cluster (n = 11), where autozygosity is strongest. This includes hets tolerated by the ROH window rule and therefore mixes genotyping error with any genuine low-level heterozygosity.

## Inputs (provenance)

- `AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv` — raw DArTseq HapMap calls.
- `refs/ayb_genome/ayb_marker_anchoring.csv` — AYB (Ss) chromosome/position per marker.
- `results/34_roh/tables/roh_calls.csv` — called ROH segments per accession (1-het-per-20-SNP-window tolerance).
- `results/90_roh_sensitivity_null/tables/froh_per_sample_obs_vs_zerohet.csv` — per-accession F_ROH (observed vs zero-het) and `is_small_cluster` label.

QC matches the ROH pipeline exactly: marker call rate ≥ 0.90, MAF ≥ 0.05, sample call rate ≥ 0.90, restricted to AYB-anchored markers → **1,473 markers × 95 samples**.

## Method

- **Trio per-genotype rate.** For each anchored post-QC marker, take the three trio calls (drop missing). Where ≥ 2 members are homozygous and agree, the true state is homozygous; count het calls (dosage = 1) among the trio at those markers as miscalls. Rate = het miscalls / trio genotypes examined. Pairwise het-vs-homozygous and opposite-homozygote discordances are also tabulated to reproduce `results/98_trio_duplicate_concordance`.
- **Het-in-ROH rate.** For each minor-cluster accession, over the markers inside each of its called ROH segments (same chromosome, position within `[start_bp, end_bp]`), count het genotypes / total non-missing genotypes. Pool. The panel-wide figure is computed the same way over all 95 samples.
- **Effect on F_ROH.** Compare main-cluster mean F_ROH under the 1-het tolerance (observed) vs zero-het re-call, and compute the expected number of error-driven false hets per 20-SNP window at the trio rate.

## Results

| Estimator | count | rate |
|---|---|---|
| Trio per-genotype het-miscall | 16 / 4,368 | **0.37 %** |
| Het inside called ROH, minor cluster | 191 / 12,587 | **1.52 %** |
| Het inside called ROH, whole panel | 1,177 / 39,166 | **3.01 %** |
| Main-cluster F_ROH (1-het tolerance, observed) | — | **0.187** |
| Main-cluster F_ROH (zero-het re-call) | — | **0.109** |
| Expected error-hets per 20-SNP window (at 0.37 %) | — | **0.073** |

Pairwise trio het-vs-hom discordances: TSs151B–TSs358 = 15, TSs151B–TSs361 = 17, TSs358–TSs361 = 6 (0/0/1 opposite-homozygote) — identical to `results/98`.

## Interpretation

1. **The genotyping het-call error floor is ≈ 0.37 %** per genotype (trio clones) — low, and consistent with the 15/17/6 het-vs-hom trio discordances already reported.
2. **Inside called ROH the het density is higher** (1.5 % minor cluster, 3.0 % panel) than the 0.37 % error floor. Because these tracts are autozygous, the excess above the floor reflects the 20-SNP window's single-heterozygote tolerance admitting **genuine low-level heterozygosity**, not error alone.
3. **The main-cluster F_ROH is bracket-dependent, not a robust point estimate.** At the 0.37 % error rate only ~0.073 false hets are expected per 20-SNP window (~7 % of windows carry one), which is far too few to explain the 0.187 → 0.109 drop under the zero-het re-call. The tolerance is therefore capturing real heterozygosity as well as error, so the honest statement is a **0.11–0.19 range** bracketed by the zero- and one-het calls.
4. **The bimodal contrast — the paper's central finding — is unaffected.** Under both calls the minor cluster (0.78 / 0.67) stands far apart from the main cluster (0.19 / 0.11).
5. **Bearing on the "bulking" question (§4).** DNA extraction is now confirmed to have been **bulked across multiple plants per accession** (Prof. Adewale, 2026-09-29). Bulking of *divergent* genotypes would inflate heterozygosity uniformly and erase ROH — but the strong per-accession autozygosity and the 99 %-identical trio show the pooled plants were themselves near-identical, as expected of a selfing landrace, so the bulk did not pool divergent genotypes. Residual heterozygosity is therefore attributed to low-rate residual outcrossing plus the ≈ 0.37 % genotyping error floor, not to between-genotype admixture in the bulk.

## Outputs

- `tables/het_error_summary.csv` — all estimators above.
- `tables/trio_pair_hetvhom.csv` — per-pair trio concordance (reproduces results/98).

## Reproduce

```bash
cd "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python scripts/100_het_call_error_rate.py
```

## Caveats

- Both estimators are **upper bounds** on genotyping error: some trio/minor-cluster hets may be genuine residual heterozygosity rather than miscalls, which would push the true error rate below 0.37 %.
- The het-in-ROH rate is not a pure error rate — it deliberately includes tolerated heterozygosity to show the window rule admits real hets.
- Extraction DNA was **bulked across multiple plants per accession** (confirmed by Prof. Adewale, 2026-09-29; now stated in Methods §2.2). This is not recorded in any data file — it is a protocol fact. The autozygosity pattern shows the bulked plants were near-identical, so the bulk does not pool divergent genotypes.

## Tool versions

Python 3.9.23 (conda env `yeast-viz`); pandas, numpy.
