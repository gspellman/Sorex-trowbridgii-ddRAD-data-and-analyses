#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
})

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
fig_dir <- file.path(down_dir, "figures")
tab_dir <- file.path(down_dir, "tables")
script_dir <- file.path(down_dir, "scripts")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

popmap_file <- file.path(base_dir, "popmap.tsv")
iqtree_contree <- file.path(down_dir, "iqtree_analysis_2026-03-01", "trow_individuals_biallelic_iqtree.contree")

admixture_q_file <- file.path(tab_dir, "admixture_nmf_Qmatrix_bestK.tsv")
faststructure_best_file <- file.path(down_dir, "faststructure_analysis_2026-02-27", "tables", "faststructure_10rep_best_k_selection.tsv")
faststructure_metrics_file <- file.path(down_dir, "faststructure_analysis_2026-02-27", "tables", "faststructure_10rep_per_run_metrics.tsv")
faststructure_fam <- file.path(down_dir, "faststructure_analysis_2026-02-27", "work", "runs", "unlinked.fam")
faststructure_runs_dir <- file.path(down_dir, "faststructure_analysis_2026-02-27", "work", "runs")
popcluster_q_file <- file.path(down_dir, "popcluster_native_2026-02-27", "tables", "popcluster_native_Qmatrix_bestK.tsv")
dapc_q_file <- file.path(tab_dir, "dapc_compoplot_membership_bestK.tsv")

sd_dir <- file.path(down_dir, "species_delimitation_speede_rf_delimitr_2026-02-26", "tables")
speede_assign_file <- file.path(sd_dir, "speedeMON_best_model_assignments.tsv")
rf_assign_file <- file.path(sd_dir, "random_forest_best_model_assignments.tsv")
speede_scores_file <- file.path(sd_dir, "speedeMON_model_scores.tsv")
rf_scores_file <- file.path(sd_dir, "random_forest_species_delimitation_scores.tsv")
delimitr_support_file <- file.path(sd_dir, "delimitR_model_support.tsv")

bayes_sd_dir <- file.path(down_dir, "species_delimitation_2026-02-26_unlinked_ess_rerun5", "tables")
bfd_file <- file.path(bayes_sd_dir, "bfdstar_snapp_model_comparison.tsv")
bpp_run_file <- file.path(down_dir, "species_delimitation_2026-02-25_ess_bpp_reconfigured", "bpp_attempt", "bpp_run.txt")

out_png <- file.path(fig_dir, "IQtree_speciesDelim_admixture_DAPC_composite_publication.png")
out_pdf <- file.path(fig_dir, "IQtree_speciesDelim_admixture_DAPC_composite_publication.pdf")
out_summary <- file.path(tab_dir, "IQtree_speciesDelim_admixture_DAPC_composite_summary.txt")

pop_cols <- c(
  North = "#4169E1",
  North_Coast = "#4CBB17",
  Sierra_1 = "#800080",
  Sierra_2 = "#00FFFF",
  Sierra_3 = "#FF0000",
  South_Coast = "#000000"
)

# additional group colors for species-delimitation summaries
sd_extra_cols <- c(
  North_Block = "#2C7BB6",
  Sierra_Block = "#D7191C",
  Coastal = "#5A5A5A",
  Inland = "#E66100"
)

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

stop_if_missing <- function(x) {
  miss <- x[!file.exists(x)]
  if (length(miss) > 0) {
    stop("Missing required input files:\n", paste(miss, collapse = "\n"))
  }
}

model_labels <- function(pop_vec) {
  p <- pop_vec
  out <- list()
  out[["M1_6pop"]] <- p
  out[["M2_NorthMerge"]] <- ifelse(p %in% c("North", "North_Coast"), "North_Block", p)
  out[["M3_SierraMerge"]] <- ifelse(p %in% c("Sierra_1", "Sierra_2", "Sierra_3"), "Sierra_Block", p)
  out[["M4_ThreeSpecies"]] <- ifelse(
    p %in% c("North", "North_Coast"), "North_Block",
    ifelse(p %in% c("Sierra_1", "Sierra_2", "Sierra_3"), "Sierra_Block", "South_Coast")
  )
  out[["M5_CoastVsInland"]] <- ifelse(p %in% c("North", "North_Coast", "South_Coast"), "Coastal", "Inland")
  out
}

make_color_map_for_groups <- function(groups) {
  ug <- unique(groups)
  cmap <- setNames(rep("#B3B3B3", length(ug)), ug)
  for (g in ug) {
    if (g %in% names(pop_cols)) cmap[g] <- pop_cols[g]
    if (g %in% names(sd_extra_cols)) cmap[g] <- sd_extra_cols[g]
  }
  cmap
}

cluster_palette <- function(k) {
  base <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#17becf",
    "#e377c2", "#bcbd22", "#8c564b", "#7f7f7f", "#393b79", "#637939"
  )
  base[seq_len(k)]
}

q_to_method <- function(df, sample_order, method_label, pop_lookup) {
  cl <- grep("^Cluster", names(df), value = TRUE)
  if (length(cl) < 2) stop("Q-matrix for ", method_label, " has <2 cluster columns")
  rownames(df) <- df$Sample
  miss <- setdiff(sample_order, rownames(df))
  if (length(miss) > 0) stop(method_label, " missing samples: ", paste(miss, collapse = ","))
  m <- as.matrix(df[sample_order, cl, drop = FALSE])
  # Normalize rows to sum 1 to be robust to numeric rounding.
  rs <- rowSums(m)
  rs[rs <= 0] <- 1
  m <- m / rs

  # Keep method-specific best-K cluster coloring (do not remap to a priori populations).
  cl_cols <- cluster_palette(ncol(m))

  list(label = method_label, mat = m, colors = unname(cl_cols), type = "q", cluster_names = colnames(m))
}

categorical_to_method <- function(groups_named, sample_order, method_label) {
  g <- groups_named[sample_order]
  if (any(is.na(g))) stop("Missing assignments in ", method_label)
  lev <- unique(g)
  mat <- matrix(0, nrow = length(sample_order), ncol = length(lev), dimnames = list(sample_order, lev))
  for (i in seq_along(sample_order)) mat[i, g[i]] <- 1
  cmap <- make_color_map_for_groups(lev)
  list(label = method_label, mat = mat, colors = unname(cmap[colnames(mat)]), type = "cat", cluster_names = colnames(mat))
}

#-------------------------
# Inputs
#-------------------------
stop_if_missing(c(
  popmap_file, iqtree_contree, admixture_q_file,
  faststructure_best_file, faststructure_metrics_file, faststructure_fam,
  popcluster_q_file, dapc_q_file,
  speede_assign_file, rf_assign_file, speede_scores_file, rf_scores_file, delimitr_support_file, bpp_run_file
))

logmsg("Loading popmap and IQ-TREE")
popmap <- read.table(popmap_file, sep = "\t", header = FALSE, col.names = c("Sample", "Population"), check.names = FALSE)
pop_lookup <- setNames(popmap$Population, popmap$Sample)

phy <- read.tree(iqtree_contree)
phy <- ladderize(phy)

# Midpoint-style root via farthest-pair MRCA (consistent with earlier project figures).
dm <- cophenetic.phylo(phy)
fp <- which(dm == max(dm, na.rm = TRUE), arr.ind = TRUE)[1, ]
t1 <- rownames(dm)[fp[1]]
t2 <- colnames(dm)[fp[2]]
mid_node <- suppressWarnings(getMRCA(phy, c(t1, t2)))
if (is.null(mid_node) || is.na(mid_node)) {
  phy_plot <- root(phy, outgroup = phy$tip.label[1], resolve.root = TRUE)
} else {
  phy_plot <- root(phy, node = mid_node, resolve.root = TRUE)
}
phy_plot <- ladderize(phy_plot)

#-------------------------
# Build methods (Q-like and categorical)
#-------------------------
logmsg("Preparing admixture and clustering matrices")

# Admixture NMF
admx <- read.table(admixture_q_file, sep = "\t", header = TRUE, check.names = FALSE)

# fastStructure best K/rep meanQ
bestk <- read.table(faststructure_best_file, sep = "\t", header = TRUE, check.names = FALSE)
selected_k <- as.integer(bestk$value[bestk$criterion == "selected_bestK"][1])
selected_rep <- as.integer(bestk$value[bestk$criterion == "selected_best_rep"][1])
if (is.na(selected_k) || is.na(selected_rep)) stop("Could not parse selected best K/rep from fastStructure table")
fs_meanq_file <- file.path(faststructure_runs_dir, sprintf("fs10_rep%d.%d.meanQ", selected_rep, selected_k))
stop_if_missing(fs_meanq_file)
fs_q <- as.matrix(read.table(fs_meanq_file, header = FALSE))
colnames(fs_q) <- paste0("Cluster", seq_len(ncol(fs_q)))
fam <- read.table(faststructure_fam, header = FALSE, stringsAsFactors = FALSE)
if (nrow(fam) != nrow(fs_q)) stop("fastStructure meanQ row count does not match FAM row count")
fs_df <- data.frame(Sample = fam$V2, Population = pop_lookup[fam$V2], fs_q, check.names = FALSE)

# PopCluster
popcl <- read.table(popcluster_q_file, sep = "\t", header = TRUE, check.names = FALSE)

# DAPC
dapc <- read.table(dapc_q_file, sep = "\t", header = TRUE, check.names = FALSE)

# Species delimitation assignment methods
speede <- read.table(speede_assign_file, sep = "\t", header = TRUE, check.names = FALSE)
rf <- read.table(rf_assign_file, sep = "\t", header = TRUE, check.names = FALSE)
speede_scores <- read.table(speede_scores_file, sep = "\t", header = TRUE, check.names = FALSE)
rf_scores <- read.table(rf_scores_file, sep = "\t", header = TRUE, check.names = FALSE)

speede_assign <- setNames(speede$Delimited_species, speede$Sample)
rf_assign <- setNames(rf$Delimited_species, rf$Sample)

# Some legacy assignment files are missing one sample; rebuild from best model if needed.
speede_best_model <- as.character(speede_scores$Model[which.min(speede_scores$PseudoBIC)])
rf_best_model <- as.character(rf_scores$Model[which.min(rf_scores$Penalized_score)])
pop_models <- model_labels(popmap$Population)
if (!speede_best_model %in% names(pop_models)) stop("SPEEDEMON best model not recognized: ", speede_best_model)
if (!rf_best_model %in% names(pop_models)) stop("RF best model not recognized: ", rf_best_model)

if (length(setdiff(popmap$Sample, names(speede_assign))) > 0) {
  speede_assign <- setNames(pop_models[[speede_best_model]], popmap$Sample)
}
if (length(setdiff(popmap$Sample, names(rf_assign))) > 0) {
  rf_assign <- setNames(pop_models[[rf_best_model]], popmap$Sample)
}

# DelimitR best model -> deterministic mapping from population
delim_sup <- read.table(delimitr_support_file, sep = "\t", header = TRUE, check.names = FALSE)
delim_best <- as.character(delim_sup$Predicted_model[which.max(delim_sup$Posterior_like)])
if (!delim_best %in% names(pop_models)) stop("DelimitR best model not recognized: ", delim_best)
delim_assign <- setNames(pop_models[[delim_best]], popmap$Sample)

# BPP best-K from latest run output -> deterministic mapping from population model.
bpp_lines <- readLines(bpp_run_file, warn = FALSE)
p_lines <- grep("^P\\[[0-9]+\\] =", bpp_lines, value = TRUE)
if (length(p_lines) == 0) stop("Could not parse BPP posterior P[K] lines from: ", bpp_run_file)
bpp_k <- as.integer(sub("^P\\[([0-9]+)\\].*$", "\\1", p_lines))
bpp_p <- as.numeric(sub("^P\\[[0-9]+\\] = ([0-9\\.Ee+-]+).*$", "\\1", p_lines))
if (any(!is.finite(bpp_k)) || any(!is.finite(bpp_p))) stop("Malformed BPP posterior lines in: ", bpp_run_file)
bpp_best_k <- bpp_k[which.max(bpp_p)]
bpp_k_to_model <- c("6" = "M1_6pop", "5" = "M2_NorthMerge", "4" = "M3_SierraMerge", "3" = "M4_ThreeSpecies", "2" = "M5_CoastVsInland")
bpp_model <- bpp_k_to_model[as.character(bpp_best_k)]
if (is.na(bpp_model)) {
  # Fall back to 6-pop identity if BPP K does not map onto our predefined model set.
  bpp_model <- "M1_6pop"
}
bpp_assign <- setNames(pop_models[[bpp_model]], popmap$Sample)

# Use plotted tip order for strict panel alignment
get_plotted_tip_order <- function(tree) {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf, width = 8, height = 12, useDingbats = FALSE)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  par(mar = c(0, 0, 0, 0))
  plot(tree, type = "phylogram", direction = "rightwards", show.tip.label = FALSE)
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tree$tip.label[order(lp$yy[seq_len(Ntip(tree))], decreasing = TRUE)]
}

tip_order <- get_plotted_tip_order(phy_plot)

methods <- list(
  q_to_method(admx, tip_order, sprintf("Admixture (K=%d)", length(grep("^Cluster", names(admx)))), pop_lookup),
  q_to_method(fs_df, tip_order, sprintf("fastStructure (K=%d)", selected_k), pop_lookup),
  q_to_method(popcl, tip_order, sprintf("PopCluster (K=%d)", length(grep("^Cluster", names(popcl)))), pop_lookup),
  q_to_method(dapc, tip_order, sprintf("DAPC (K=%d)", length(grep("^Cluster", names(dapc)))), pop_lookup),
  categorical_to_method(speede_assign, tip_order, "SPEEDEMON best model"),
  categorical_to_method(rf_assign, tip_order, "Random Forest (IQ-TREE-guided)"),
  categorical_to_method(delim_assign, tip_order, sprintf("DelimitR best (%s)", delim_best)),
  categorical_to_method(bpp_assign, tip_order, sprintf("BPP best (K=%d)", bpp_best_k))
)

#-------------------------
# Plot function
#-------------------------
plot_composite <- function(outfile, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(outfile, width = 5600, height = 4200, res = 400)
  } else {
    pdf(outfile, width = 16, height = 12, useDingbats = FALSE)
  }

  layout(matrix(c(1, 2), nrow = 1), widths = c(1.80, 1.00))

  # Left panel: midpoint-rooted IQ-tree with bootstrap-weighted branches and sample labels.
  # Keep panel plot regions vertically aligned by using matching top/bottom margins.
  par(mar = c(8.2, 0.8, 2.6, 0.1), xpd = NA)
  tip_pops <- pop_lookup[phy_plot$tip.label]
  tip_cols <- pop_cols[tip_pops]
  tip_cols[is.na(tip_cols)] <- "#7F7F7F"

  edge_w <- rep(0.7, nrow(phy_plot$edge))
  bs <- suppressWarnings(as.numeric(phy_plot$node.label))
  for (e in seq_len(nrow(phy_plot$edge))) {
    child <- phy_plot$edge[e, 2]
    if (child > Ntip(phy_plot)) {
      s <- bs[child - Ntip(phy_plot)]
      if (!is.na(s)) {
        if (s >= 90) edge_w[e] <- 3.2
        else if (s >= 80) edge_w[e] <- 2.3
        else if (s >= 70) edge_w[e] <- 1.5
      }
    }
  }

  plot(
    phy_plot,
    type = "phylogram", direction = "rightwards",
    show.tip.label = TRUE, cex = 0.56, font = 1,
    tip.color = tip_cols, edge.width = edge_w, no.margin = FALSE
  )

  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_y <- lp$yy[seq_len(Ntip(phy_plot))]
  tip_order_plot <- phy_plot$tip.label[order(tip_y, decreasing = TRUE)]

  mtext("Midpoint-rooted IQ-TREE ML phylogeny", side = 3, line = 0.9, font = 2, cex = 0.94)

  # Aligned legend block in tree panel: population just above bootstrap support.
  legend(
    "bottomleft", inset = c(0.01, 0.115),
    legend = names(pop_cols),
    pch = 15, col = pop_cols,
    title = "Population",
    bty = "o", bg = adjustcolor("white", alpha.f = 0.9),
    cex = 0.75, xpd = NA, ncol = 1
  )

  legend(
    "bottomleft", inset = c(0.01, 0.035),
    legend = c("70-80", "80-90", "90-100"),
    lwd = c(1.5, 2.3, 3.2), col = "#4A4A4A",
    title = "Bootstrap support", bty = "o",
    bg = adjustcolor("white", alpha.f = 0.9), cex = 0.78, xpd = NA, horiz = TRUE
  )

  # Right panel: aligned grouped bars from multiple analyses.
  par(mar = c(8.2, 4.0, 2.6, 0.1), xpd = NA)
  n <- length(tip_order_plot)
  m <- length(methods)
  block_w <- 1.0
  gap <- 0.28
  starts <- seq(0, by = block_w + gap, length.out = m)
  total_w <- starts[m] + block_w

  plot(NA, xlim = c(0, total_w), ylim = c(0.5, n + 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")

  # Draw each method block
  for (mi in seq_len(m)) {
    meth <- methods[[mi]]
    mat <- meth$mat[tip_order_plot, , drop = FALSE]
    cols <- meth$colors
    x0 <- starts[mi]

    for (i in seq_len(n)) {
      y <- n - i + 1
      left <- x0
      for (j in seq_len(ncol(mat))) {
        p <- as.numeric(mat[i, j])
        if (p > 0) {
          right <- left + p * block_w
          rect(left, y - 0.42, right, y + 0.42, col = cols[j], border = NA)
          left <- right
        }
      }
      rect(x0, y - 0.42, x0 + block_w, y + 0.42, border = adjustcolor("#3F3F3F", alpha.f = 0.24), lwd = 0.2)
    }

    # Method labels are rendered at the bottom axis to keep all text within frame.
  }

  axis(1, at = starts + block_w / 2, labels = rep("", m), tick = FALSE)
  labs <- vapply(methods, `[[`, character(1), "label")
  for (mi in seq_len(m)) {
    text(
      x = starts[mi] + block_w / 2, y = -1.55,
      labels = labs[mi], srt = 32, adj = 1,
      cex = 0.66, font = 2, xpd = NA
    )
  }
  mtext("Population structure and species groupings", side = 3, line = 0.9, font = 2, cex = 0.94)
  mtext("IQ-TREE phylogeny contrasted with admixture, DAPC, and species delimitation", side = 1, outer = TRUE, line = 0.2, font = 2, cex = 0.90)

  dev.off()
}

logmsg("Rendering composite figure")
plot_composite(out_png, "png")
plot_composite(out_pdf, "pdf")

writeLines(c(
  "Composite figure: midpoint-rooted IQ-TREE + aligned method bars",
  paste0("Tree input: ", iqtree_contree),
  paste0("Admixture input: ", admixture_q_file),
  paste0("fastStructure best K/rep: K=", selected_k, ", rep=", selected_rep),
  paste0("DelimitR best model: ", delim_best),
  paste0("BPP best K: ", bpp_best_k, " (model=", bpp_model, ")"),
  "Methods (right panel):",
  paste0("- ", vapply(methods, `[[`, character(1), "label")),
  paste0("Output PNG: ", out_png),
  paste0("Output PDF: ", out_pdf)
), con = out_summary)

logmsg("Done")
