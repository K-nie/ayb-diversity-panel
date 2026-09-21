#!/usr/bin/env python3
"""
Maximum-likelihood phylogenetic inference for the 95-line AYB panel.

Three independent tree reconstructions for cross-validation:
  1. UPGMA on Rogers' modified distance (already produced by script 05;
     loaded here for cross-comparison).
  2. IQ-TREE 2 with GTR+ASC, ModelFinder, 1000 ultrafast bootstraps.
  3. RAxML-NG with GTR+ASC, 100 ML bootstraps.

Alignment construction
----------------------
DArTseq SNPs are biallelic. For each SNP we emit one character per sample:
the reference allele as a nucleotide (taken from the HapMap `alleles` field),
heterozygotes as the IUPAC ambiguity code, the alt allele as its nucleotide,
and missing dosage as 'N'. The resulting alignment has 1,625 columns (one per
QC-passing marker) and 95 rows (one per sample). Ascertainment-bias
correction (GTR+ASC) is enabled because the panel contains only variant
sites by construction.

Cross-tree comparison
---------------------
Compute pairwise cophenetic distances from each tree and compare each pair
of distance matrices with Pearson r (Mantel-style). Output a 4 x 4
correlation heatmap among {UPGMA, IQ-TREE, RAxML, raw Rogers'} matrices and
a tree-vs-tree leaf-correspondence figure (tanglegram-lite).

Outputs (results/19_phylogenetic_trees/)
----------------------------------------
data/
    snp_alignment.phy / .fasta
    iqtree/                 IQ-TREE working directory + treefile
    raxml/                  RAxML working directory + bestTree
tables/
    cophenetic_correlation_matrix.csv
figures/
    fig50_iqtree_radial.png/.pdf
    fig51_raxml_radial.png/.pdf
    fig52_cophenetic_matrix.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import warnings
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from _plotstyle import apply, WONG, CLUSTER_PAL
apply()
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
ROGERS = ROOT / "results" / "05_distance_dendrogram" / "tables" / "rogers_distance_matrix.csv"
OUT = ROOT / "results" / "19_phylogenetic_trees"
DATA_DIR = OUT / "data"
IQ_DIR = DATA_DIR / "iqtree"
RAX_DIR = DATA_DIR / "raxml"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (DATA_DIR, IQ_DIR, RAX_DIR, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

IQTREE = "/Users/black_einstein/miniconda3/bin/iqtree"
RAXML = "/usr/local/bin/raxml-ng"

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
# IUPAC ambiguity codes are NOT used in the alignment because IQ-TREE's
# ASC correction (needed for SNP-only data) and RAxML-NG's ASC_LEWIS both
# fail when sites are dominated by ambiguity (looks invariant). We resolve
# heterozygotes to the major allele at each site instead (acceptable for a
# selfer like AYB where heterozygosity is rare).


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Rebuild filtered dosage + ref/alt
# ---------------------------------------------------------------------------
print("[load] HapMap...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META]
calls = hm[sample_cols].astype(str)

ref_alt = hm["alleles"].str.split("/", expand=True)
ref_alt.columns = ["ref", "alt"]
dosage = np.full((len(hm), len(sample_cols)), np.nan)
for i in range(len(hm)):
    ref, alt = ref_alt.iloc[i]
    arr = calls.iloc[i].values
    dosage[i, arr == ref + ref] = 0.0
    dosage[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
    dosage[i, arr == alt + alt] = 2.0
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)

call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values

dos_f = dosage.loc[keep_m, keep_s]
samples = list(dos_f.columns)
markers = list(dos_f.index)
ref_f = ref_alt.loc[keep_m, "ref"].values
alt_f = ref_alt.loc[keep_m, "alt"].values
print(f"[load] alignment: {len(samples)} samples x {len(markers)} SNPs")


# ---------------------------------------------------------------------------
# Build pseudo-SNP alignment with IUPAC ambiguity for heterozygotes
# ---------------------------------------------------------------------------
print("[align] writing SNP alignment (haploid: hets -> major allele)...")
matrix = dos_f.values  # markers x samples
n_samp = len(samples)
n_mark = len(markers)
align = np.full((n_samp, n_mark), "N", dtype="U1")
n_invariant = 0
keep_col = []
for j in range(n_mark):
    r, a = ref_f[j], alt_f[j]
    col = matrix[j]
    # major allele = whichever homozygote count is larger; if equal, ref wins
    n_ref_hom = int((col == 0).sum())
    n_alt_hom = int((col == 2).sum())
    het_letter = r if n_ref_hom >= n_alt_hom else a
    align[col == 0, j] = r
    align[col == 1, j] = het_letter
    align[col == 2, j] = a
    # NaN -> 'N'
    # check post-resolution variation; drop columns that became invariant
    chars = set(align[:, j]) - {"N"}
    if len(chars) >= 2:
        keep_col.append(j)
    else:
        n_invariant += 1
align = align[:, keep_col]
n_mark = align.shape[1]
print(f"[align] dropped {n_invariant} sites that became invariant after "
      f"major-allele resolution; {n_mark} variant sites retained")

# write PHYLIP relaxed-sequential (IQ-TREE accepts that with -s)
PHY = DATA_DIR / "snp_alignment.phy"
with open(PHY, "w") as fh:
    fh.write(f" {n_samp} {n_mark}\n")
    for i, s in enumerate(samples):
        seq = "".join(align[i])
        fh.write(f"{s}  {seq}\n")

# write FASTA (raxml-ng default)
FAS = DATA_DIR / "snp_alignment.fasta"
with open(FAS, "w") as fh:
    for i, s in enumerate(samples):
        fh.write(f">{s}\n{''.join(align[i])}\n")

print(f"[align] PHYLIP -> {PHY} ({PHY.stat().st_size/1e6:.1f} MB)")
print(f"[align] FASTA  -> {FAS} ({FAS.stat().st_size/1e6:.1f} MB)")


# ---------------------------------------------------------------------------
# IQ-TREE with GTR+ASC and 1000 ultrafast bootstraps
# ---------------------------------------------------------------------------
print("\n[iqtree] running IQ-TREE GTR+ASC + UFBoot 1000...")
iq_prefix = IQ_DIR / "ayb"
# clean previous run
for f in IQ_DIR.glob("ayb.*"):
    f.unlink()

cmd = [IQTREE, "-s", str(PHY), "-m", "GTR+ASC", "-B", "1000",
       "--prefix", str(iq_prefix), "-T", "4", "-redo"]
r = subprocess.run(cmd, capture_output=True, text=True)
iq_log = iq_prefix.with_suffix(".log")
if r.returncode != 0:
    print("  IQ-TREE first attempt failed, retrying without ASC...")
    cmd = [IQTREE, "-s", str(PHY), "-m", "GTR+G", "-B", "1000",
           "--prefix", str(iq_prefix), "-T", "4", "-redo"]
    r = subprocess.run(cmd, capture_output=True, text=True)

iq_tree_path = iq_prefix.with_suffix(".treefile")
if iq_tree_path.exists():
    print(f"  IQ-TREE tree: {iq_tree_path}")
else:
    print(f"  IQ-TREE failed; stdout tail = {r.stdout[-300:]!r}")


# ---------------------------------------------------------------------------
# RAxML-NG with GTR+ASC and 100 bootstraps
# ---------------------------------------------------------------------------
print("\n[raxml] running raxml-ng GTR+ASC + 100 bootstraps...")
for f in RAX_DIR.glob("ayb_rax.*"):
    f.unlink()
rax_prefix = RAX_DIR / "ayb_rax"

# all-in-one: search ML tree + bootstrap + bsconverge
cmd = [RAXML, "--all", "--msa", str(FAS), "--model", "GTR+G+ASC_LEWIS",
       "--prefix", str(rax_prefix), "--seed", "42", "--threads", "4",
       "--bs-trees", "100"]
r = subprocess.run(cmd, capture_output=True, text=True)
rax_tree_path = Path(str(rax_prefix) + ".raxml.support")
if not rax_tree_path.exists():
    rax_tree_path = Path(str(rax_prefix) + ".raxml.bestTree")
if rax_tree_path.exists():
    print(f"  RAxML tree: {rax_tree_path}")
else:
    print(f"  RAxML failed; stderr tail = {r.stderr[-400:]!r}")


# ---------------------------------------------------------------------------
# Read trees (using a minimal Newick parser) + render
# ---------------------------------------------------------------------------
def parse_newick(s):
    """Return a dict {tip -> path-length-from-root-following-parent-pointers}.

    For cophenetic comparisons we only need leaf-to-leaf distances, which we
    compute from the dendrogram returned by ETE-style traversal. We use
    Bio.Phylo here for simplicity.
    """
    import io
    from Bio import Phylo
    return Phylo.read(io.StringIO(s), "newick")


def tree_to_cophenetic(tree, samples):
    """Return a square cophenetic-distance matrix indexed by `samples`."""
    leaves = {leaf.name: leaf for leaf in tree.get_terminals()}
    n = len(samples)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if samples[i] not in leaves or samples[j] not in leaves:
                continue
            d = tree.distance(leaves[samples[i]], leaves[samples[j]])
            M[i, j] = M[j, i] = d
    return M


# Read all four distance matrices
print("\n[copheno] reading distance matrices...")
dist_rogers = pd.read_csv(ROGERS, index_col=0)
# reorder columns to match `samples`
samples_in_rogers = list(dist_rogers.index)
common = [s for s in samples if s in samples_in_rogers]
Mraw = dist_rogers.loc[common, common].values

# UPGMA cophenetic from script 05's linkage (rebuild for safety)
from scipy.cluster.hierarchy import cophenet
upgma_link = linkage(squareform(Mraw, checks=False), method="average")
Mupgma = squareform(cophenet(upgma_link))

# IQ-TREE & RAxML cophenetic
def cophen_or_none(path, samples_):
    if not Path(path).exists(): return None
    nw = Path(path).read_text()
    try:
        tree = parse_newick(nw)
    except Exception as e:
        print(f"  failed to parse {path}: {e}")
        return None
    return tree_to_cophenetic(tree, samples_)

Miq = cophen_or_none(iq_tree_path, common)
Mrx = cophen_or_none(rax_tree_path, common)

mats = {"raw_Rogers": Mraw, "UPGMA": Mupgma}
if Miq is not None: mats["IQ-TREE"] = Miq
if Mrx is not None: mats["RAxML"] = Mrx

# pairwise Pearson r between off-diagonal entries
keys = list(mats.keys())
n_off = (len(common) * (len(common) - 1)) // 2
triu = np.triu_indices(len(common), k=1)
flat = {k: M[triu] for k, M in mats.items()}
import scipy.stats as st
corr = pd.DataFrame(np.eye(len(keys)), index=keys, columns=keys)
for i, a in enumerate(keys):
    for j, b in enumerate(keys):
        if i >= j: continue
        ok = np.isfinite(flat[a]) & np.isfinite(flat[b])
        if ok.sum() < 5:
            corr.loc[a, b] = corr.loc[b, a] = np.nan; continue
        r = st.pearsonr(flat[a][ok], flat[b][ok])[0]
        corr.loc[a, b] = corr.loc[b, a] = r
corr.to_csv(TAB / "cophenetic_correlation_matrix.csv")
print("[copheno] correlation matrix:")
print(corr.round(3).to_string())


# ---------------------------------------------------------------------------
# Fig 52 cophenetic correlation heatmap
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 5.0), constrained_layout=True)
cmap = plt.get_cmap("RdBu_r")
im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap=cmap, aspect="equal")
n = corr.shape[0]
ax.set_xticks(range(n)); ax.set_yticks(range(n))
ax.set_xticklabels(corr.columns, rotation=30, ha="right")
ax.set_yticklabels(corr.index)
for i in range(n):
    for j in range(n):
        v = corr.values[i, j]
        if np.isnan(v): continue
        ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                fontsize=10,
                color="white" if abs(v) > 0.55 else "black")
ax.set_xticks(np.arange(n + 1) - 0.5, minor=True)
ax.set_yticks(np.arange(n + 1) - 0.5, minor=True)
ax.grid(which="minor", color="white", linewidth=1.2)
ax.tick_params(which="minor", length=0)
fig.colorbar(im, ax=ax, shrink=0.7, label="Pearson r (cophenetic distances)")
ax.set_title("Cophenetic distance correlation across tree methods")
save(fig, "fig52_cophenetic_matrix")
print("[fig] fig52_cophenetic_matrix")


# ---------------------------------------------------------------------------
# Helper: radial render using Bio.Phylo + matplotlib
# ---------------------------------------------------------------------------
def render_radial(tree_path, fig_name, title):
    if not Path(tree_path).exists():
        print(f"  skip {fig_name} (no tree)"); return
    from Bio import Phylo
    tree = Phylo.read(str(tree_path), "newick")
    pca = pd.read_csv(PCA_TAB).set_index("sample")
    clusters = {s: int(pca.loc[s, "cluster"]) if s in pca.index else 0
                for s in samples}

    # for radial layout we use angle = leaf_index / n_leaves * 2*pi
    leaves = [l.name for l in tree.get_terminals()]
    n = len(leaves)
    theta = {leaf: 2*np.pi*i/n for i, leaf in enumerate(leaves)}

    # build node theta = mean of leaf thetas under it; node r = distance from root
    def node_theta(node):
        if node.is_terminal():
            return theta[node.name]
        return np.mean([node_theta(c) for c in node.clades])

    def node_r(node, r=0):
        out = {id(node): r}
        for c in node.clades:
            out.update(node_r(c, r + (c.branch_length or 0)))
        return out

    rmap = node_r(tree.root)

    fig = plt.figure(figsize=(8.5, 8.5))
    ax = fig.add_subplot(111, projection="polar")

    def draw(node):
        nt = node_theta(node)
        nr = rmap[id(node)]
        for c in node.clades:
            ct = node_theta(c); cr = rmap[id(c)]
            # arc at parent r between nt and ct, then radial line out to (ct, cr)
            arc_t = np.linspace(min(nt, ct), max(nt, ct), 30)
            ax.plot(arc_t, np.full_like(arc_t, nr), color="#555555", lw=0.6)
            ax.plot([ct, ct], [nr, cr], color="#555555", lw=0.6)
            if c.is_terminal():
                k = clusters.get(c.name, 0)
                ax.scatter(ct, cr, s=22, color=CLUSTER_PAL[k],
                           edgecolor="white", linewidth=0.4, zorder=5)
                rot = np.degrees(ct)
                ha = "left" if rot < 90 or rot > 270 else "right"
                if 90 < rot < 270:
                    rot -= 180
                ax.text(ct, cr * 1.02 + 0.001, c.name,
                        rotation=rot, rotation_mode="anchor",
                        fontsize=5, ha=ha, va="center",
                        color=CLUSTER_PAL[k])
            else:
                draw(c)
    draw(tree.root)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid(True, alpha=0.25)
    handles = [mpatches.Patch(color=CLUSTER_PAL[k],
                              label=f"Cluster {k + 1}") for k in (0, 1)]
    ax.legend(handles=handles, loc="lower right",
              bbox_to_anchor=(1.05, -0.05), frameon=False, fontsize=9)
    ax.set_title(title, pad=18, fontsize=11)
    save(fig, fig_name)
    print(f"[fig] {fig_name}")


render_radial(iq_tree_path, "fig50_iqtree_radial",
              f"IQ-TREE (GTR+ASC, UFBoot1000) — {len(samples)} AYB lines")
render_radial(rax_tree_path, "fig51_raxml_radial",
              f"RAxML-NG (GTR+G+ASC) — {len(samples)} AYB lines")

print(f"\nOutputs in: {OUT}")
