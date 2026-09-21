# 91 — Genome-wide LD decay and Ne on the AYB-anchored marker set

Refits the panel-wide LD-decay curve and the LD-based effective-population-size (Ne) estimate on the *Sphenostylis stenocarpa* (AYB) anchoring, replacing the submitted global fit that rested on the legacy cowpea-anchored SNP subset. Addresses Reviewer 1 major point 6.

## Why

The submitted manuscript reported a panel-wide **Ne ≈ 659** and a **74 kb** global half-decay from a Hill–Weir fit run on the **cowpea-anchored** SNP subset (176 SNPs carrying paired *Vigna unguiculata* positions in the original DArT report). That is the same legacy anchoring the revision otherwise abandons in favour of the AYB assembly, so the headline demographic number was inconsistent with the rest of the paper. The reviewer asked us to **refit the global decay on the AYB-anchored set or drop the Ne claim**.

The submitted per-chromosome fits (Supplementary Figure S4A) were already AYB-anchored, but there was no *genome-wide* AYB fit — the main-text Figure 3 and the headline Ne came only from the cowpea subset. This analysis pools the AYB-anchored per-chromosome r² pair tables into a single genome-wide decay curve and refits Hill–Weir, so the panel-wide summary is on the same anchoring as everything else.

## Inputs

- `results/87_panel_wide_ld_heatmap/tables/ld_pairs_Ss01..Ss11.csv` — all within-chromosome SNP pairs on the AYB-anchored QC marker set, with pre-computed r² (columns `rs_i, rs_j, pos_i, pos_j, r2`). These are the same pair tables used for the panel-wide LD heatmap (script 87), computed on the **1,473 AYB-anchored QC markers** (call-rate ≥ 0.90, MAF ≥ 0.05).

## Commands (reproduce)

```bash
python scripts/91_ld_decay_genomewide_ayb.py   # env: yeast-viz
```

## Parameters

Identical modelling choices to the submitted fit, only the anchoring changes:

- Pairs restricted to ≤ 5 Mb within each chromosome (same window as the submission).
- 40 equal-width distance bins; mean r² per bin; Hill & Weir (1988) drift–recombination expectation E[r²] = (10 + C) / ((2 + C)(11 + C)), C = ρ·d(Morgans).
- ρ̂ fit by non-linear least squares on the binned means (`scipy.optimize.curve_fit`, bounds 1e-3–1e6).
- Ne = ρ̂ / 4 under the **same 1.06 cM/Mb cowpea map-rate proxy** used in the submission — so the AYB global value is directly comparable to the submitted per-chromosome AYB fits and to the cowpea global value it replaces. The map-rate proxy remains an acknowledged source of absolute-scale error (see Caveats); it affects the cowpea and AYB fits identically, so the *comparison* between them is unaffected.

## Findings

**Genome-wide AYB-anchored fit** (32,507 within-chromosome pairs ≤ 5 Mb, 1,473 anchored markers):

| Quantity | AYB-anchored (this run) | Submitted cowpea-anchored |
|---|---|---|
| ρ̂ = 4Ne | 1,904 | 2,636 |
| **Ne (1.06 cM/Mb proxy)** | **≈ 476** | ≈ 659 |
| Half-decay distance | **102 kb** | 74 kb |
| r² < 0.1 by | 370 kb | 267 kb |

**1. The headline Ne moves from ≈ 659 to ≈ 476 on the correct anchoring, and the half-decay lengthens from 74 kb to 102 kb.** The AYB-anchored genome-wide Ne (476) sits squarely inside the submitted per-chromosome AYB-anchored Ne range (295–759; ρ̂ 1,179–3,034) and near its centre, whereas the cowpea global value (659) sat at the upper edge. Refitting on the AYB assembly therefore lowers the panel-wide Ne modestly and lengthens the half-decay, but does not change the qualitative story: a small, highly autogamous panel with LD that decays over ~100 kb.

**2. The autogamy interpretation is unchanged.** A ~102 kb half-decay remains at the long end of legume LD-decay values for selfing species (cowpea ≈ 50–100 kb; soybean ≈ 100–200 kb; chickpea ≈ 100 kb), consistent with the panel's cleistogamous/autogamous reproductive biology. The revision should report the **AYB-anchored genome-wide fit** (Ne ≈ 476, half-decay 102 kb, r² < 0.1 by 370 kb) as the headline and can retain the cowpea value only as a cross-check, explicitly labelled legacy.

## Outputs

- `tables/ld_genomewide_ayb_summary.csv` — ρ̂, Ne, half-decay, r²<0.1 distance, n pairs / markers.
- `tables/ld_genomewide_ayb_bins.csv` — the binned decay curve (mean/sd/count/SE per distance bin).
- `figures/fig_ld_decay_genomewide_ayb.png` / `.pdf` — genome-wide AYB decay + Hill–Weir fit; **replacement for the submitted Figure 3**.

## Caveats

- **Absolute Ne depends on the assumed recombination rate.** No AYB genetic map exists, so the 1.06 cM/Mb cowpea rate is a proxy; the true AYB rate could shift absolute Ne up or down proportionally. This limitation is shared with the submitted fit and does not affect the cowpea-vs-AYB *comparison* (both use the same rate) or the relative decay profile / half-decay in physical distance.
- The pair tables pool all within-chromosome pairs across the 11 pseudo-chromosomes; per-chromosome heterogeneity in decay rate (ρ̂ 1,179–3,034) is retained in the per-chromosome fits (Supplementary Figure S4A, script 11b) and is summarised, not erased, by the genome-wide fit.

## Script

`scripts/91_ld_decay_genomewide_ayb.py`. Env: `yeast-viz` (numpy, pandas, scipy, matplotlib). Runtime: < 10 s.
