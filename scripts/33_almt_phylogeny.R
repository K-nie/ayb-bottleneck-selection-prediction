#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Render the ALMT-family phylogeny with AYB clades highlighted.
#
# Input  : results/33_almt_phylogeny/data/ALMT_tree.contree (IQ-TREE consensus
#          with ultrafast-bootstrap support values)
# Output : results/33_almt_phylogeny/figures/fig72_almt_family_tree.png/.pdf
#
# AYB tips are pre-named AYB_ALMT*_Ss* (script 33_almt_phylogeny.py step).
# Reference tips are UniProt headers; we annotate species + AtALMT number.
#
# Author: Benjamin Narh-Madey
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
})

ROOT <- "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
TREE_FILE <- file.path(ROOT,
                       "results/33_almt_phylogeny/data/ALMT_tree.contree")
FIG <- file.path(ROOT, "results/33_almt_phylogeny/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

cat("[load] tree...\n")
tree <- read.tree(TREE_FILE)
tree <- midpoint.root(tree)
cat(sprintf("  %d tips\n", length(tree$tip.label)))

# Classify each tip
classify_tip <- function(label) {
  if (grepl("^AYB_", label)) return("AYB")
  if (grepl("ARATH", label)) return("Arabidopsis")
  if (grepl("SOYBN|MAX$", label, ignore.case = TRUE)) return("Soybean")
  if (grepl("PHAVU", label, ignore.case = TRUE)) return("Common bean")
  if (grepl("VIGUN", label, ignore.case = TRUE)) return("Cowpea")
  if (grepl("LOTUS|LOTJA", label, ignore.case = TRUE)) return("Other")
  return("Other")
}
species <- sapply(tree$tip.label, classify_tip)

# Wong-style palette
COLORS <- c(AYB = "#D55E00",  Arabidopsis = "#0072B2",
            Soybean = "#009E73",  `Common bean` = "#CC79A7",
            Cowpea = "#E69F00", Other = "#999999")
tip_col <- COLORS[species]

# Tidy display labels (drop UniProt prefixes for the species ones; AYB labels
# are already display-ready)
display_label <- function(lbl) {
  if (grepl("^AYB_", lbl)) return(lbl)
  # UniProt header is like "sp|Q9C6L8|ALMT4_ARATH ..."; keep the protein code
  m <- regmatches(lbl, regexec("^[a-z]+\\|[A-Z0-9]+\\|([A-Z0-9_]+)", lbl))[[1]]
  if (length(m) >= 2) return(m[2])
  return(lbl)
}
tree$tip.label_display <- sapply(tree$tip.label, display_label)

# Render as a ggtree circular (fan) tree with species-coloured tip labels,
# extra emphasis on the AYB tips, and bootstrap support shown only on
# well-supported internal nodes. Per-tip text fontsize is tuned to a 14"
# square canvas so all 142 tips stay legible.

tip_df <- data.frame(label = tree$tip.label,
                      display = tree$tip.label_display,
                      species = species,
                      is_ayb = species == "AYB",
                      stringsAsFactors = FALSE)

# Bootstrap support frame for ggtree nodes; ufboot < 70 is dropped to avoid
# cluttering the figure with low-confidence labels.
bs_vals <- suppressWarnings(as.numeric(tree$node.label))
n_tips <- length(tree$tip.label)
bs_df <- data.frame(
  node = (n_tips + 1):(n_tips + tree$Nnode),
  ufboot = bs_vals,
  ufboot_label = ifelse(!is.na(bs_vals) & bs_vals >= 70,
                          as.character(round(bs_vals)),
                          NA_character_)
)

p <- ggtree(tree, layout = "fan", open.angle = 18,
             size = 0.55, colour = "grey30") %<+% tip_df %<+% bs_df +
  # Faint species-coloured radial strip at each tip so the species
  # identity reads even when the tip label itself is hard to track
  # in the dense outer ring.
  geom_tippoint(aes(colour = species, size = is_ayb),
                  show.legend = TRUE) +
  scale_size_manual(values = c(`TRUE` = 4.6, `FALSE` = 2.0),
                      guide = "none") +
  geom_tiplab(aes(label = display, colour = species,
                    fontface = ifelse(is_ayb, "bold", "plain")),
               size = 3.5, offset = 0.025, show.legend = FALSE,
               family = "sans") +
  scale_colour_manual(values = COLORS, name = "Species",
                       breaks = names(COLORS),
                       guide = guide_legend(override.aes =
                                              list(size = 5.5))) +
  geom_text2(aes(label = ufboot_label, subset = !is.na(ufboot_label)),
              size = 2.5, colour = "grey30",
              hjust = 1.15, vjust = -0.35) +
  ggtitle(paste0("ALMT-family protein phylogeny -- 142 sequences,\n",
                   "LG+G4 with 1000 UFBoot; AYB tips in vermillion bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14,
                                    face = "bold"),
        legend.position = c(0.94, 0.12),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))

# Larger canvas so the 142 tip labels actually have radial real estate.
ggsave(file.path(FIG, "fig72_almt_family_tree.png"), p,
       width = 18, height = 18, dpi = 300, bg = "white")
ggsave(file.path(FIG, "fig72_almt_family_tree.pdf"), p,
       width = 18, height = 18, bg = "white")
cat(sprintf("[fig] fig72_almt_family_tree (ggtree, PNG + PDF) written to %s\n", FIG))

# also pull out the AYB clade context: which Arabidopsis/legume ALMT is each
# AYB tip's nearest neighbour?
cat("\n[neighbour] AYB ALMT nearest-neighbour in the tree:\n")
d_mat <- cophenetic(tree)
for (ayb_tip in tree$tip.label[species == "AYB"]) {
  others <- setdiff(tree$tip.label, ayb_tip)
  others_d <- d_mat[ayb_tip, others]
  top3 <- sort(others_d)[1:3]
  cat(sprintf("  %s -> closest 3:\n", ayb_tip))
  for (n in names(top3)) {
    cat(sprintf("     %s   d = %.4f\n", n, top3[n]))
  }
}
