#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
})

base <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24"
work <- file.path(base, "work")
fig <- file.path(base, "figures")
tab <- file.path(base, "tables")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
dir.create(tab, recursive = TRUE, showWarnings = FALSE)

tree_file <- file.path(work, "ml_biallelic_individuals.treefile")
contree_file <- file.path(work, "ml_biallelic_individuals.contree")
meta_file <- file.path(tab, "ml_biallelic_sample_metadata.tsv")

tr <- read.tree(contree_file)
meta <- read.table(meta_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pop_map <- setNames(meta$Population, meta$Sample)

pop_cols <- c(
  North = "#4169E1",
  North_Coast = "#4CBB17",
  Sierra_1 = "#800080",
  Sierra_2 = "#00FFFF",
  Sierra_3 = "#FF0000",
  South_Coast = "#000000",
  Unknown = "#7f7f7f"
)

tip_pop <- pop_map[tr$tip.label]
tip_pop[is.na(tip_pop)] <- "Unknown"
tip_col <- unname(pop_cols[tip_pop])
tip_col[is.na(tip_col)] <- "#7f7f7f"

n_tip <- length(tr$tip.label)
edge_child <- tr$edge[, 2]
node_support <- suppressWarnings(as.numeric(tr$node.label))

support_by_edge <- rep(NA_real_, nrow(tr$edge))
for (i in seq_along(edge_child)) {
  ch <- edge_child[i]
  if (ch > n_tip) {
    idx <- ch - n_tip
    if (idx >= 1 && idx <= length(node_support)) {
      support_by_edge[i] <- node_support[idx]
    }
  }
}

edge_width <- rep(0.8, nrow(tr$edge))
bin_70_80 <- !is.na(support_by_edge) & support_by_edge >= 70 & support_by_edge < 80
bin_80_90 <- !is.na(support_by_edge) & support_by_edge >= 80 & support_by_edge < 90
bin_90_100 <- !is.na(support_by_edge) & support_by_edge >= 90 & support_by_edge <= 100
edge_width[bin_70_80] <- 1.6
edge_width[bin_80_90] <- 2.8
edge_width[bin_90_100] <- 4.0

plot_one <- function(outfile, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(outfile, width = 3200, height = 3200, res = 360)
  } else {
    pdf(outfile, width = 10, height = 10, useDingbats = FALSE)
  }

  par(mar = c(0.8, 0.8, 2.6, 0.8), xpd = NA)
  plot.phylo(
    tr,
    type = "unrooted",
    use.edge.length = TRUE,
    no.margin = TRUE,
    show.tip.label = FALSE,
    edge.width = edge_width,
    edge.color = "#555555",
    cex = 0.50
  )

  tiplabels(
    pch = 21,
    bg = tip_col,
    col = "white",
    cex = 1.25
  )

  title("Unrooted ML phylogeny (biallelic SNPs); branch width scaled by bootstrap bins", cex.main = 0.95, line = 0.2)

  legend(
    "topleft",
    legend = names(pop_cols)[names(pop_cols) != "Unknown"],
    pch = 21,
    pt.bg = unname(pop_cols[names(pop_cols) != "Unknown"]),
    col = "white",
    bty = "n",
    pt.cex = 1.2,
    cex = 0.9,
    y.intersp = 1.0,
    x.intersp = 0.7
  )

  legend(
    "topright",
    legend = c("70-80", "80-90", "90-100"),
    lwd = c(1.6, 2.8, 4.0),
    col = "#555555",
    bty = "n",
    cex = 0.9,
    title = "Bootstrap support"
  )

  dev.off()
}

plot_one(file.path(fig, "ML_unrooted_biallelic_bootstrap_weighted_publication.png"), "png")
plot_one(file.path(fig, "ML_unrooted_biallelic_bootstrap_weighted_publication.pdf"), "pdf")

write.tree(tr, file = file.path(tab, "ml_biallelic_individuals_bootstrap_contree.newick"))
write.tree(read.tree(tree_file), file = file.path(tab, "ml_biallelic_individuals_mltree.newick"))

edge_tab <- data.frame(
  parent = tr$edge[,1],
  child = tr$edge[,2],
  bootstrap_support = support_by_edge,
  edge_width = edge_width
)
write.table(edge_tab, file = file.path(tab, "ml_biallelic_edge_support_weights.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_txt <- file.path(tab, "ml_biallelic_phylogeny_summary.txt")
cat(
  "Maximum likelihood phylogeny summary\n",
  sprintf("Tips (samples): %d\n", n_tip),
  sprintf("Internal nodes with support labels: %d\n", sum(!is.na(node_support))),
  sprintf("Edges weighted by support 70-80: %d\n", sum(bin_70_80)),
  sprintf("Edges weighted by support 80-90: %d\n", sum(bin_80_90)),
  sprintf("Edges weighted by support 90-100: %d\n", sum(bin_90_100)),
  "Figures:\n",
  "- ML_unrooted_biallelic_bootstrap_weighted_publication.png\n",
  "- ML_unrooted_biallelic_bootstrap_weighted_publication.pdf\n",
  file = summary_txt,
  sep = ""
)
