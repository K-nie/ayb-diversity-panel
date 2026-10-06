# A1 cluster-label harmonization (Reviewer 2, point C, item 11)

## What this does and why

Reviewer 2 flagged that the two population clusters are numbered inconsistently
across the A1 manuscript text, tables, figures, and deposited supplementary
files. The two clusters are a **main bulk of 84 accessions** and a **minor /
PC1-outlier group of 11 accessions**. The manuscript body (para ~77 and Table 3)
refers to them by the words *main* and *minor*, but the deposited per-accession
and per-pair supplementary tables carried an **integer** cluster code whose
ordering disagreed with the body convention and, internally, between the two
generator scripts. An integer label whose meaning depends on which script wrote
it is exactly the ambiguity the reviewer caught.

The fix removes the ordering question entirely: the deposited cluster columns now
carry the **descriptive words `main` / `minor`**, matching Table 3's wording, so
the label is convention-independent and cannot be misread.

## Authoritative membership (source of truth)

`results/01_qc_pca_power/tables/pca_coords.csv` assigns the raw k-means/PCA
cluster code per accession:

| raw PCA code | size | identity |
|---|---|---|
| `1` | 84 | main bulk |
| `0` | 11 | minor / PC1 outliers |

## The integer-convention drift (the root cause)

The two generator scripts map the raw PCA code to a published integer
differently:

- `scripts/34_roh.py` L211: `cluster = 2 - pca_cluster` → main(1)→**1**, minor(0)→**2**
  (intends *Cluster 1 = main*; comment L208–209 says so).
- `scripts/92_cross_shortlist_dedup_normalized.py` L143–144: `cluster = pca_cluster + 1`
  → minor(0)→**1**, main(1)→**2** (*Cluster 1 = minor*).

The two scripts therefore emit opposite integer orderings. The deposited
`TableS4` additionally predates the `34_roh.py` `2 - pca` fix, so it was written
with the `pca+1` ordering (Cluster 1 = minor), the same ordering as the current
`TableS12`. Both deposited files thus arrived at integer **1 = minor, 2 = main** —
the opposite of the manuscript body's *Cluster 1 = main*.

## What was changed

The integer cluster codes in the deposited supplementary tables were mapped
to descriptive words, touching only the cluster column(s) — all other
columns (relatedness, F_ROH, phenotype mid-parent/gap values) are byte-for-byte
unchanged (minimal token replacement, no float reformatting). **Two distinct
integer conventions were present**, so the word map differs per file; each was
verified against the authoritative PCA membership before relabeling:

- `manuscript/paper_A1/share/tables/supp_tables/TableS4_f_roh_per_sample.csv`
  — column `cluster`, map `1 → minor`, `2 → main`: 95 cells relabeled (84 `main`, 11 `minor`).
- `manuscript/paper_A1/share/tables/supp_tables/TableS12_cross_shortlist_dedup_normalized.csv`
  — columns `cluster_A`, `cluster_B`, same `1 → minor`, `2 → main` map: 100 cells relabeled.
- `manuscript/paper_A1/share/tables/supp_tables/TableS5_core_collection_top20_stratified.csv`
  — column `cluster`, map `0 → minor`, `1 → main` (the **raw PCA code**, not the
  shifted published integer used in S4/S12): 20 cells relabeled (16 `main`, 4 `minor`).
  Added 2026-10-05: S5 was missed in the first sweep and is the table Reviewer 2
  explicitly asks for under R2-D ("provide the 20 core members as a table"), so an
  integer code here would have re-triggered R2-C point 11 on that very artifact.

## Verification

- `TableS4` cluster 1 had exactly 11 members, cluster 2 had 84 — matches the PCA
  membership, so `1 = minor`, `2 = main`.
- Every `TableS12` parent-cluster assignment was cross-checked against `TableS4`
  by accession: **100/100 agreed, 0 mismatches** — the two files share one integer
  convention, so the same word map applies to both.
- `TableS5` uses the raw PCA cluster code directly: each of its 20 members was
  cross-checked against `results/01_qc_pca_power/tables/pca_coords.csv`
  (`0` = minor bulk of 11, `1` = main bulk of 84) — **20/20 agreed, 0 mismatches**,
  confirming `0 → minor`, `1 → main`. The relabel touched only the `cluster`
  column; **0 non-cluster cells changed**.
- Post-relabel counts: `TableS4` {main: 84, minor: 11};
  `TableS12` cluster_A {minor: 29, main: 21}, cluster_B {main: 45, minor: 5};
  `TableS5` {main: 16, minor: 4}.

## Caveat / remaining source cleanup

The deposited copies are now unambiguous and final for this submission. The
upstream integer conventions in `scripts/34_roh.py` (`2 - pca`) and
`scripts/92_cross_shortlist_dedup_normalized.py` (`pca + 1`) still disagree with
each other; a future full re-run would re-emit integers. Harmonizing both
scripts to emit the descriptive `main` / `minor` label (or a single agreed
integer ordering) is the recommended source-level follow-up. It was **not**
re-run here — the relabel was applied to the deposited artifacts directly, which
is why this note records the exact mapping and verification rather than a command
to reproduce.

The figure-level cluster key (Cluster 1 = main = blue, Cluster 2 = minor =
vermillion) is handled in the separate A1 figure-convention sweep.

## Tool versions

- Python 3.9.23 (`/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python`),
  pandas 2.3.1.
