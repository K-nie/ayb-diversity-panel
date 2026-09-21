#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# IQ-TREE 3.0.1 ML phylogeny with ADMIXTURE Q (K = 4) overlaid as per-tip
# stacked bars. The combination shows topology and ancestry side-by-side --
# tips that the tree groups together should share their dominant Q component
# at K = 4, and the figure also exposes admixed individuals that span clades.
#
# Two layouts are written:
#   fig83_iqtree_admixture_Q4_rectangular.png/.pdf  -- standard left-right
#   fig83_iqtree_admixture_Q4_circular.png/.pdf     -- fan layout for the
#                                                       supplementary panel
#
# Author: Benjamin Narh-Madey
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(tidyr)
})

ROOT <- "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
IQ_TREE <- file.path(ROOT, "results/19_phylogenetic_trees/data/iqtree/ayb.treefile")
Q_FILE <- file.path(ROOT, "results/09_admixture/admixture/ayb.4.Q")
FAM_FILE <- file.path(ROOT, "results/09_admixture/plink/ayb.fam")
FIG <- file.path(ROOT, "results/19_phylogenetic_trees/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# Wong palette for the four ADMIXTURE clusters
WONG_K4 <- c(K1 = "#0072B2", K2 = "#D55E00",
              K3 = "#009E73", K4 = "#E69F00")

cat("[load] IQ-TREE...\n")
tree <- read.tree(IQ_TREE)
cat(sprintf("  %d tips\n", length(tree$tip.label)))

cat("[load] ADMIXTURE Q (K=4)...\n")
Q <- read.table(Q_FILE)
fam <- read.table(FAM_FILE, stringsAsFactors = FALSE)
sample_ids <- fam$V2
stopifnot(nrow(Q) == length(sample_ids))
colnames(Q) <- paste0("K", seq_len(ncol(Q)))
Q$sample <- sample_ids

# Restrict to tips present in BOTH IQ-TREE and Q
common <- intersect(tree$tip.label, Q$sample)
cat(sprintf("  IQ-TREE tips         : %d\n", length(tree$tip.label)))
cat(sprintf("  Q samples            : %d\n", length(sample_ids)))
cat(sprintf("  intersection         : %d\n", length(common)))
tree <- keep.tip(tree, common)
Q <- Q[match(common, Q$sample), ]

# Per-tip "dominant cluster" for tip-label colouring
Q$dominant <- apply(Q[, paste0("K", 1:4)], 1, function(r) {
  paste0("K", which.max(r))
})

# ---------------------------------------------------------------------------
# Rectangular layout
# ---------------------------------------------------------------------------
p_rect <- ggtree(tree, layout = "rectangular", size = 0.4, colour = "grey30") %<+%
  Q[, c("sample", "dominant")] +
  geom_tiplab(aes(colour = dominant), size = 2.6, offset = 0.001) +
  scale_colour_manual(values = WONG_K4,
                       name = "Dominant K=4 cluster") +
  ggtitle("IQ-TREE 3.0.1 ML phylogeny with tips coloured by dominant ADMIXTURE K=4 cluster") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

# Tip-point overlay also coloured by dominant cluster, plus an emphasised
# coloured-square at each tip so the figure works without the Q-bar strip.
final_rect <- p_rect +
  geom_tippoint(aes(colour = dominant), size = 2.2)

ggsave(file.path(FIG, "fig83_iqtree_admixture_Q4_rectangular.png"),
       final_rect, width = 11, height = 14, dpi = 300, bg = "white")
ggsave(file.path(FIG, "fig83_iqtree_admixture_Q4_rectangular.pdf"),
       final_rect, width = 11, height = 14, bg = "white")
cat("[fig] fig83_iqtree_admixture_Q4_rectangular\n")

# ---------------------------------------------------------------------------
# Circular layout
# ---------------------------------------------------------------------------
p_circ <- ggtree(tree, layout = "fan", open.angle = 20,
                  size = 0.55, colour = "grey30") %<+%
  Q[, c("sample", "dominant")] +
  geom_tippoint(aes(colour = dominant), size = 3.4,
                  show.legend = TRUE) +
  geom_tiplab(aes(colour = dominant), size = 3.4, offset = 0.001,
               fontface = "bold", family = "sans") +
  scale_colour_manual(values = WONG_K4,
                       name = "Dominant K=4 cluster",
                       guide = guide_legend(override.aes =
                                              list(size = 6))) +
  ggtitle(paste0("IQ-TREE 3.0.1 ML phylogeny of ", length(common),
                   " accessions\n",
                   "tip colour = dominant ADMIXTURE K=4 cluster")) +
  theme(plot.title = element_text(size = 14, hjust = 0.5,
                                    face = "bold"),
        legend.position = c(0.93, 0.13),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))

# Larger canvas so 95 accession labels read clearly around the perimeter.
ggsave(file.path(FIG, "fig83_iqtree_admixture_Q4_circular.png"),
       p_circ, width = 16, height = 16, dpi = 300, bg = "white")
ggsave(file.path(FIG, "fig83_iqtree_admixture_Q4_circular.pdf"),
       p_circ, width = 16, height = 16, bg = "white")
cat("[fig] fig83_iqtree_admixture_Q4_circular\n")

cat("\n[done] both layouts written to ", FIG, "\n")
