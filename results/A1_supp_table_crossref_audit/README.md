# A1 supplementary-table cross-reference audit (Reviewer 2, points C & D follow-through)

## What this does and why

While resolving Reviewer 2's cross-shortlist / de-duplication point (R2-D) the
deposited supplementary-table set was audited end-to-end: every supplementary-table
cross-reference in the A1 manuscript, the Supplementary Information (SI), and the
response letter was traced to the deposited CSV it names, and every deposited
CSV in `manuscript/paper_A1/share/tables/supp_tables/` was traced back to a
citation. The panel now carries **15 numbered tables (S1–S15)** plus the
sub-lettered variants (S2b, S3b–S3e). The audit found four internal
inconsistencies, all editorial (no analysis re-run); each is recorded below with
the verification that justified the fix.

## Defects found and fixes applied (2026-10-05)

### 1. Supplementary Table S15 had no SI caption
`TableS15_het_vs_depth_outbred_lines.csv` (the Reviewer-2-point-B het-vs-depth
outbred-tail table, analysis in `results/101_het_vs_depth/`) was deposited and
referenced in the manuscript Results (§3.6) but the SI caption list stopped at
S14. A three-paragraph caption block (title / filename / description) was inserted
after the S14 block in `SI_A1.docx` as a tracked insertion, cloning the S14 block's
formatting and matching the SI tracked-change convention (author
`Benjamin Narh-Madey`, date `2026-09-30T00:00:00Z`, ins ids 103–109).

### 2. Stale "S1–S13" supplementary-table range
Both the manuscript Data-Availability paragraph and the response letter stated the
deposited range as "Supplementary Tables S1–S13", predating S14 (phenotype matrix)
and S15 (het-vs-depth). Corrected to **S1–S15** in:
- `A1_manuscript_final_submitted.docx` — the paragraph already carried a tracked
  edit `del "9" / ins "13"` (an in-round revision of the original submitted
  "S1–S9"); the insertion text was advanced `13 → 15`, yielding a clean
  `del "9" / ins "15"`. Tracked-change counts unchanged (ins 221, del 121).
- `RESPONSE_LETTER_A1.docx` — plain text edit (the letter carries no tracked
  changes; 0 ins / 0 del), "(S1–S13)" → "(S1–S15)".

### 3. Table S10 caption mislabelled its own row count
The S14 phenotype-matrix caption cross-references S10 as the accession list
("…for the 105 genotyped accessions (Supplementary Table S10)"), but S10's own
caption read "Full **95**-accession manifest…". Verification of the deposited file
settled which number is right:
- `TableS10_accession_list_qc_cluster.csv` has **105 data rows** and carries
  `qc_status` + `qc_fail_reason` columns; the QC split is **95 Pass, 10 Fail**
  (all 10 fail on `call rate < 0.90`). A manifest that encodes QC pass/fail must
  list all 105 genotyped accessions, so "95-accession manifest" was the error —
  not the S14 cross-reference.
- The S10 caption was corrected (tracked del+ins, same SI convention) to
  "Full **105**-accession manifest of all genotyped accessions (95 pass QC;
  10 fail on call rate < 0.90): …", making S10 and the S14 cross-reference
  internally consistent. The S14 caption was left unchanged (its pointer is now
  correct).

### 4. Nine orphan `TableS_NEW_*` provisional-name files in the deposit
The deposited `supp_tables/` folder still held nine `TableS_NEW_*.csv` files —
provisional names from the 2026-09-20 reviewer-response runs, superseded by the
numbered S1–S15 scheme. None is referenced by name in the manuscript, SI, or
response letter (the two apparent substring hits were false positives:
`TableS11_marker_placement_uniqueness.csv` in the SI caption, and the path
`results/90_roh_sensitivity_null/` in the letter). Each was confirmed
**byte-identical (md5)** to a canonical `results/` source before removal:

| removed file | byte-identical canonical source |
|---|---|
| `TableS_NEW_cross_shortlist_dedup_normalized.csv` | `results/92_cross_shortlist_dedup_normalized/tables/cross_pairs_least_related_dedup_normalized.csv` (final deposit is the harmonised `TableS12`) |
| `TableS_NEW_froh_obs_vs_zerohet.csv` | `results/90_roh_sensitivity_null/tables/froh_per_sample_obs_vs_zerohet.csv` |
| `TableS_NEW_ld_genomewide_ayb_summary.csv` | `results/91_ld_decay_genomewide_ayb/tables/ld_genomewide_ayb_summary.csv` |
| `TableS_NEW_marker_placement_uniqueness.csv` | `results/89_blast_uniqueness_reanchor/tables/marker_placement_classification.csv` (= deposited `TableS11`) |
| `TableS_NEW_meff_artefact_summary.csv` | `results/93_meff_samplesize_artefact/tables/meff_artefact_summary.csv` |
| `TableS_NEW_meff_vs_samplesize.csv` | `results/93_meff_samplesize_artefact/tables/meff_vs_samplesize.csv` |
| `TableS_NEW_placement_uniqueness_summary.csv` | `results/89_blast_uniqueness_reanchor/tables/placement_summary.csv` |
| `TableS_NEW_roh_sensitivity_summary.csv` | `results/90_roh_sensitivity_null/tables/roh_sensitivity_summary.csv` |
| `TableS_NEW_shortlist_reconciliation.csv` | `results/92_cross_shortlist_dedup_normalized/tables/shortlist_reconciliation_summary.csv` |

Content is preserved in `results/` (and in the GitHub deposit cited in the Data
Availability Statement). Safety copies of the removed files are at
`/tmp/A1_TableS_NEW_removed_20261005/`.

## Verification

- SI docx after edits: paragraphs 109 → 112 (+3 S15 caption paras), tables 1
  (unchanged), ins 32 → 38, del 21 → 22; S15 accepted-view shows title/filename/
  description; S10 accepted-view reads the corrected 105-accession wording.
- Manuscript docx after edit: paragraphs 426, tables 5, ins 221, del 121 — all
  identical to the pre-edit baseline (the range change edited inserted text only);
  accepted-view reads "Supplementary Figures S1–S16 and Supplementary Tables
  S1–S15 …".
- Response letter after edit: 63 paragraphs, 0 ins / 0 del (plain letter); reads
  "Supplementary Tables (S1–S15)".
- `supp_tables/` now holds 20 numbered deposit files (S1–S15 + S2b, S3b–S3e); no
  `TableS_NEW_*` remain.

## Caveat / remaining cleanup (not applied here)

- `TableS5_core_collection_top20_stratified.csv.PRErelabel` (a working backup from
  the cluster-label relabel) still sits in the deposited `supp_tables/` folder, as
  do several `*.PRE*` / `*.BACKUP*` working copies at `share/` top level. These are
  local working artifacts, not deposit tables; they were left in place because this
  pass was scoped to the nine `TableS_NEW_*` files. Recommend removing them from the
  share deposit before the package is zipped for submission.
- The share package (`share/`) does not itself contain a `results/` tree; the
  `results/` artifacts reach reviewers through the GitHub repository named in the
  Data Availability Statement, not through the submission bundle.

## Tool versions

- Python 3.9.23 (`/Users/black_einstein/miniconda3/envs/yeast-viz/bin/python`),
  lxml; docx edited as raw `word/document.xml` (zip round-trip), md5 via `md5 -q`.
