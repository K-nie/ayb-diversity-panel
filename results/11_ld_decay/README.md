# 11 — Linkage-disequilibrium decay + N_e estimate

**Script:** `scripts/11_ld_decay.py`

## What was done

Pairwise r² LD computed for all SNP pairs within chromosomes up to 5 Mb physical separation, binned by distance, then fit to the Hill & Weir (1988) drift-recombination expectation to estimate effective population size.

## Method

Standardised dosage → pairwise correlation r → r². The Hill & Weir (1988) E[r²] = (10 + C) / ((2 + C)(11 + C)) form was fit by non-linear least squares against the per-bin means, with C = ρ × c (recombination distance in Morgans). ρ = 4 N_e c assumed 1.06 cM/Mb (Lonardi et al. 2019 cowpea). Half-decay distance = where E[r²] = max/2.

## Quick findings

- 478 SNP pairs computed across 175 cowpea-anchored SNPs within 5 Mb.
- **LD half-decay = 74 kb**.
- r² < 0.1 by **267 kb**.
- ρ̂ = 4 N_e = 2,636 → **N_e ≈ 659**.

## Caveats

- LD decay is computed on the **cowpea-anchored subset only** (175 SNPs, ~11 % of the panel). The AYB-anchored re-run would give similar half-decay but is not the bottleneck for this stat.
- Hill-Weir is the canonical inbred-population expectation; the assumption of constant N_e through time is a known simplification.
- Half-decay of 74 kb is on the longer end for legumes (cowpea ~50–100 kb in Muñoz-Amatriaín et al. 2017), consistent with the autogamous reproduction mode.
