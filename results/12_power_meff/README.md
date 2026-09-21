# 12 — Daetwyler ceiling + Li & Ji effective number of tests

**Script:** `scripts/12_power_meff.py`

## What was done

Two power / multiple-testing summaries: (i) the Daetwyler (2008) genomic-prediction accuracy ceiling for the panel across a h² × n grid, and (ii) the Li & Ji (2005) eigenvalue-based effective number of independent markers (m_eff) for the Bonferroni threshold.

## Method

**Daetwyler ceiling**: r² = n h² / (n h² + M_e), with M_e = 2 N_e L / ln(4 N_e L) (Goddard 2009). L = genome length in Morgans, derived from N_e = 659 (from script 11) and cowpea genome length 620.4 Mb × 1.06 cM/Mb → L = 6.58 M, giving M_e = 888.

**m_eff**: eigenvalues of the m × m marker correlation matrix (computed via the n × n Gram of standardised columns, since rank ≤ n − 1). Li & Ji 2005 formula: m_eff = Σ I(λ_i ≥ 1) + (λ_i − ⌊λ_i⌋) for λ_i < 1.

## Quick findings

- **M_e (Goddard 2009) = 888**.
- **m_eff (Li & Ji 2005) = 93** (of 1,625 raw markers).
- Bonferroni at raw m = 1,625: α = 3.08 × 10⁻⁵ (−log₁₀ = 4.51).
- Bonferroni at m_eff = 93: α = 5.35 × 10⁻⁴ (−log₁₀ = 3.27). *This is the threshold used in the GWAS Manhattans.*
- Daetwyler ceiling at n = 95: r ≤ 0.07 (h² = 0.05); **0.15 (h² = 0.20)**; 0.22 (h² = 0.50); 0.28 (h² = 0.80).
- To reach r = 0.5 at h² = 0.50, n must reach ~671 — a 7× expansion of the current panel.

## Caveats

- Daetwyler ceiling assumes the trait is genuinely polygenic with the given h²; large-effect QTL trump this expectation.
- m_eff is bounded above by n − 1 = 94 at this sample size, and our value 93 is close to that bound. The bound itself reflects rank limitations, not a biological feature.
- Observed Antioxidant GBLUP r = 0.15 (script 03) matches the h² = 0.20 ceiling exactly — the panel is operating at its theoretical limit.
