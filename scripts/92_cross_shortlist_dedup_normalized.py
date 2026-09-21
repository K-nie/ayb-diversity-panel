#!/usr/bin/env python3
"""
92 - De-duplicated, inbreeding-normalized cross shortlist.

Reviewer 2 point D: the submitted 50-pair "least-related" cross shortlist was
ranked on raw VanRaden (2008) G_ij off-diagonals over the *full* 95-line panel,
including the near-clonal trio (TSs151B, TSs358, TSs361). Two problems:

  1. The trio are cryptic duplicates -- keeping all three in the candidate set
     lets a "diverse" cross pair a line against a near-identical copy of another
     candidate, and inflates the apparent number of independent lines.
  2. Raw G_ij is inbreeding-scaled: G_ii = 1 + F_i, so a highly inbred line
     (the trio have G_ii ~ 1.6, i.e. F ~ 0.6, vs panel mean 1.008) is pushed to
     extreme off-diagonal values against everyone. Ranking on raw G_ij therefore
     confounds genuine complementarity with each line's own inbreeding level,
     and the trio dominates the "least-related" tail.

Fix (this script):
  - De-duplicate: collapse the Gij>=0.85 duplicate set (the trio) to a single
    representative -- the member with the highest marker call rate -- and drop
    the redundant copies before shortlisting.
  - Inbreeding-normalize relatedness: r_ij = G_ij / sqrt(G_ii * G_jj), the
    genomic analogue of a correlation. This is bounded and comparable across
    lines of differing inbreeding, so "least related" means genuinely least
    co-ancestry rather than "most inbred".
  - Re-rank the 50 least-related pairs on the de-duplicated panel using r_ij,
    carry the 13-trait phenotype context, and quantify how much the shortlist
    changes vs the submitted raw-G/full-panel version.

Inputs
  - results/13_grm_crosspairs/tables/grm_vanraden.csv        (95 x 95 G, incl. diag)
  - results/13_grm_crosspairs/tables/cross_pairs_least_related.csv (submitted top 50)
  - results/01_qc_pca_power/tables/{sample_qc.csv, pca_coords.csv}
  - phenotypes via scripts/_pheno.load_phenotypes

Outputs (results/92_cross_shortlist_dedup_normalized/)
  tables/cross_pairs_least_related_dedup_normalized.csv
  tables/shortlist_reconciliation_summary.csv
  figures/fig_cross_shortlist_rawG_vs_normalized.png/.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _plotstyle import apply, WONG  # noqa: E402
from _pheno import load_phenotypes, TRAITS_ALL  # noqa: E402
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
GRM = ROOT / "results" / "13_grm_crosspairs" / "tables" / "grm_vanraden.csv"
SUBMITTED = ROOT / "results" / "13_grm_crosspairs" / "tables" / "cross_pairs_least_related.csv"
SAMPLE_QC = ROOT / "results" / "01_qc_pca_power" / "tables" / "sample_qc.csv"
PCA = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"

OUT = ROOT / "results" / "92_cross_shortlist_dedup_normalized"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

DUP_GIJ_THRESHOLD = 0.85     # same threshold used for duplicate detection in the paper
TOP_N_PAIRS = 50


def main() -> None:
    G = pd.read_csv(GRM, index_col=0)
    samples = list(G.index)
    Gv = G.values
    n = len(samples)
    diag = np.diag(Gv)
    print(f"[grm] {n} lines; diag mean {diag.mean():.3f} "
          f"(min {diag.min():.3f}, max {diag.max():.3f}); implied F = diag - 1")

    # ---- identify duplicate set(s) at Gij >= threshold, collapse by transitive closure
    ti, tj = np.triu_indices(n, 1)
    gij = Gv[ti, tj]
    dup_pairs = [(samples[ti[k]], samples[tj[k]]) for k in np.where(gij >= DUP_GIJ_THRESHOLD)[0]]
    # transitive closure -> connected components
    from collections import defaultdict
    adj = defaultdict(set)
    for a, b in dup_pairs:
        adj[a].add(b)
        adj[b].add(a)
    seen, comps = set(), []
    for s in samples:
        if s in adj and s not in seen:
            stack, comp = [s], []
            while stack:
                u = stack.pop()
                if u in seen:
                    continue
                seen.add(u)
                comp.append(u)
                stack.extend(adj[u] - seen)
            comps.append(sorted(comp))
    print(f"[dedup] duplicate groups at Gij>={DUP_GIJ_THRESHOLD}: {comps}")

    # ---- pick representative per group = highest marker call rate
    qc = pd.read_csv(SAMPLE_QC).set_index("sample")["call_rate"]
    drop = set()
    reps = {}
    for comp in comps:
        crs = {m: float(qc.get(m, np.nan)) for m in comp}
        rep = max(comp, key=lambda m: (crs[m] if np.isfinite(crs[m]) else -1))
        reps[rep] = comp
        for m in comp:
            if m != rep:
                drop.add(m)
        print(f"[dedup] group {comp}: keep {rep} "
              f"(call rates {{'{rep}': {crs[rep]:.3f}}}), drop {sorted(set(comp)-{rep})}")

    keep = [s for s in samples if s not in drop]
    keep_idx = [samples.index(s) for s in keep]
    print(f"[dedup] panel {n} -> {len(keep)} lines after collapsing duplicates")

    # ---- inbreeding-normalized relatedness on the FULL matrix, then subset
    d_sqrt = np.sqrt(np.clip(diag, 1e-9, None))
    R = Gv / np.outer(d_sqrt, d_sqrt)      # r_ij = G_ij / sqrt(G_ii G_jj)
    Gk = Gv[np.ix_(keep_idx, keep_idx)]
    Rk = R[np.ix_(keep_idx, keep_idx)]

    pca = pd.read_csv(PCA).set_index("sample")["cluster"].astype(int)
    pheno = load_phenotypes()

    kti, ktj = np.triu_indices(len(keep), 1)
    r_off = Rk[kti, ktj]
    g_off = Gk[kti, ktj]
    order_low = np.argsort(r_off)[:TOP_N_PAIRS]     # least related by normalized r
    rows = []
    for k in order_low:
        i, j = kti[k], ktj[k]
        a, b = keep[i], keep[j]
        row = {
            "parent_A": a, "parent_B": b,
            "r_ij_normalized": float(Rk[i, j]),
            "G_ij_raw": float(Gk[i, j]),
            "cluster_A": int(pca.get(a, -2)) + 1,
            "cluster_B": int(pca.get(b, -2)) + 1,
        }
        for t in TRAITS_ALL:
            va = pheno.at[a, t] if a in pheno.index else np.nan
            vb = pheno.at[b, t] if b in pheno.index else np.nan
            row[f"{t}_midparent"] = (va + vb) / 2.0 if pd.notna(va) and pd.notna(vb) else np.nan
            row[f"{t}_absgap"] = abs(va - vb) if pd.notna(va) and pd.notna(vb) else np.nan
        rows.append(row)
    new_df = pd.DataFrame(rows)
    new_df.to_csv(TAB / "cross_pairs_least_related_dedup_normalized.csv", index=False)
    print(f"[shortlist] new top {TOP_N_PAIRS} least-related (normalized r_ij): "
          f"r range [{new_df['r_ij_normalized'].min():.3f}, {new_df['r_ij_normalized'].max():.3f}]")

    # ---- reconciliation vs the submitted raw-G / full-panel shortlist
    sub = pd.read_csv(SUBMITTED)
    trio = sorted({m for comp in comps for m in comp})
    sub_pairset = {frozenset((r.parent_A, r.parent_B)) for r in sub.itertuples()}
    new_pairset = {frozenset((r.parent_A, r.parent_B)) for r in new_df.itertuples()}
    overlap = sub_pairset & new_pairset
    sub_with_trio = sum(1 for p in sub_pairset if p & set(trio))
    sub_with_dropped = sum(1 for p in sub_pairset if p & drop)
    new_with_dropped = sum(1 for p in new_pairset if p & drop)

    summ = pd.DataFrame([
        {"metric": "panel size submitted", "value": n},
        {"metric": "panel size after dedup", "value": len(keep)},
        {"metric": "duplicate lines dropped", "value": ";".join(sorted(drop))},
        {"metric": "duplicate representative kept", "value": ";".join(sorted(reps))},
        {"metric": "submitted shortlist pairs involving a trio member", "value": sub_with_trio},
        {"metric": "submitted shortlist pairs involving a DROPPED duplicate", "value": sub_with_dropped},
        {"metric": "new shortlist pairs involving a DROPPED duplicate", "value": new_with_dropped},
        {"metric": "shortlist pairs shared (submitted ∩ new)", "value": len(overlap)},
        {"metric": "shortlist pairs replaced", "value": TOP_N_PAIRS - len(overlap)},
        {"metric": "new normalized r_ij min", "value": round(float(new_df['r_ij_normalized'].min()), 4)},
        {"metric": "new normalized r_ij max", "value": round(float(new_df['r_ij_normalized'].max()), 4)},
    ])
    summ.to_csv(TAB / "shortlist_reconciliation_summary.csv", index=False)
    print("\n[reconcile]")
    print(summ.to_string(index=False))

    # ---- figure: raw G_ij vs normalized r_ij for every retained pair,
    #      with dropped-duplicate-involving pairs (submitted artefact) marked
    fig, ax = plt.subplots(figsize=(7.0, 5.2), constrained_layout=True)
    # all pairs on the FULL panel (to show where the trio sat)
    full_r = R[np.triu_indices(n, 1)]
    full_g = Gv[np.triu_indices(n, 1)]
    involves_drop = np.array([
        (samples[ti[k]] in drop) or (samples[tj[k]] in drop)
        for k in range(len(ti))
    ])
    ax.scatter(full_g[~involves_drop], full_r[~involves_drop], s=8, alpha=0.25,
               color=WONG["skyblue"], edgecolor="none",
               label="pair (both retained lines)")
    ax.scatter(full_g[involves_drop], full_r[involves_drop], s=22, alpha=0.9,
               color=WONG["vermillion"], edgecolor="black", linewidth=0.3,
               label="pair involving a dropped duplicate (trio)")
    # mark the submitted least-related tail
    ax.axvline(sub["G_ij"].max(), color="0.4", linestyle="--", linewidth=0.8,
               label=f"submitted least-related cutoff (G ≤ {sub['G_ij'].max():.2f})")
    ax.set_xlabel(r"raw VanRaden G$_{ij}$ (submitted ranking axis)")
    ax.set_ylabel(r"inbreeding-normalized relatedness r$_{ij}$ = G$_{ij}$/$\sqrt{G_{ii}G_{jj}}$")
    ax.set_title("Cross-pair relatedness: raw G$_{ij}$ vs inbreeding-normalized r$_{ij}$\n"
                 f"{sub_with_trio}/{TOP_N_PAIRS} submitted least-related pairs "
                 f"involve the inbred trio", fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True)
    fig.savefig(FIG / "fig_cross_shortlist_rawG_vs_normalized.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "fig_cross_shortlist_rawG_vs_normalized.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[fig] fig_cross_shortlist_rawG_vs_normalized")


if __name__ == "__main__":
    main()
