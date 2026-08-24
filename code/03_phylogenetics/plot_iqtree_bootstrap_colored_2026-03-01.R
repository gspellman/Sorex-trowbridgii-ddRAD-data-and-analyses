#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(ape))

analysis_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23/iqtree_analysis_2026-03-01"
contree_file <- file.path(analysis_dir, "trow_individuals_biallelic_iqtree.contree")
popmap_file <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/popmap.tsv"

tree <- read.tree(contree_file)
pop <- read.table(popmap_file, sep = "\t", header = FALSE, col.names = c("Sample","Population"), quote = "", comment.char = "")

# population palette requested by user
pop_cols <- c(
  North = "#4169E1",        # royal blue
  North_Coast = "#4CBB17",  # Kelly green
  Sierra_1 = "#800080",     # purple
  Sierra_2 = "#00FFFF",     # cyan
  Sierra_3 = "#FF0000",     # red
  South_Coast = "#000000"   # black
)

# tip populations/colors
pop_map <- setNames(pop$Population, pop$Sample)
tip_pop <- pop_map[tree$tip.label]
missing_pop <- is.na(tip_pop)
if (any(missing_pop)) tip_pop[missing_pop] <- "Unknown"

tip_col <- pop_cols[tip_pop]
tip_col[is.na(tip_col)] <- "#666666"

Ntip <- length(tree$tip.label)
Nedge <- nrow(tree$edge)

# Support values from internal node labels in contree
support_by_node <- rep(NA_real_, Ntip + tree$Nnode)
if (!is.null(tree$node.label)) {
  vals <- suppressWarnings(as.numeric(tree$node.label))
  support_by_node[(Ntip + 1):(Ntip + tree$Nnode)] <- vals
}

edge_support <- rep(NA_real_, Nedge)
for (i in seq_len(Nedge)) {
  child <- tree$edge[i, 2]
  if (child > Ntip) edge_support[i] <- support_by_node[child]
}

# branch width by support bins
edge_lwd <- rep(0.7, Nedge)
edge_lwd[!is.na(edge_support) & edge_support >= 70 & edge_support < 80] <- 1.6
edge_lwd[!is.na(edge_support) & edge_support >= 80 & edge_support < 90] <- 2.6
edge_lwd[!is.na(edge_support) & edge_support >= 90] <- 3.8

# save support table
support_tbl <- data.frame(
  parent_node = tree$edge[,1],
  child_node = tree$edge[,2],
  is_internal_child = tree$edge[,2] > Ntip,
  bootstrap_support = edge_support,
  lwd = edge_lwd
)
write.table(support_tbl, file.path(analysis_dir, "trow_individuals_biallelic_iqtree_branch_support.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# publication figure
png(file.path(analysis_dir, "trow_individuals_biallelic_iqtree_bootstrap_weighted_colored.png"),
    width = 3200, height = 3200, res = 400)
par(mar = c(1.2, 1.2, 1.8, 1.2), xpd = NA)
plot(tree,
     type = "unrooted",
     show.tip.label = FALSE,
     edge.color = "#4A4A4A",
     edge.width = edge_lwd,
     cex = 0.6,
     no.margin = TRUE)

# add colored tips
tiplabels(pch = 21, bg = tip_col, col = "white", cex = 1.15, lwd = 0.3)

title(main = "IQ-TREE ML phylogeny (biallelic SNPs) with bootstrap-weighted branches", cex.main = 1.2)

# legends
legend("topleft", legend = names(pop_cols), pt.bg = unname(pop_cols), pch = 21,
       pt.cex = 1.2, cex = 0.9, bty = "n", title = "Population")
legend("bottomleft", legend = c("70-80", "80-90", "90-100"),
       lwd = c(1.6, 2.6, 3.8), col = "#4A4A4A", cex = 0.9, bty = "n",
       title = "Bootstrap support")
dev.off()

pdf(file.path(analysis_dir, "trow_individuals_biallelic_iqtree_bootstrap_weighted_colored.pdf"),
    width = 8.2, height = 8.2)
par(mar = c(1.2, 1.2, 1.8, 1.2), xpd = NA)
plot(tree,
     type = "unrooted",
     show.tip.label = FALSE,
     edge.color = "#4A4A4A",
     edge.width = edge_lwd,
     cex = 0.6,
     no.margin = TRUE)
tiplabels(pch = 21, bg = tip_col, col = "white", cex = 1.1, lwd = 0.3)
title(main = "IQ-TREE ML phylogeny (biallelic SNPs) with bootstrap-weighted branches", cex.main = 1.05)
legend("topleft", legend = names(pop_cols), pt.bg = unname(pop_cols), pch = 21,
       pt.cex = 1.0, cex = 0.82, bty = "n", title = "Population")
legend("bottomleft", legend = c("70-80", "80-90", "90-100"),
       lwd = c(1.6, 2.6, 3.8), col = "#4A4A4A", cex = 0.82, bty = "n",
       title = "Bootstrap support")
dev.off()

cat("Wrote figure and support table to:\n", analysis_dir, "\n", sep = "")
