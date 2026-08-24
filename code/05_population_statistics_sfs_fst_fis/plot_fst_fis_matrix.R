#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

out_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

fst <- as.matrix(read.table(file.path(tab_dir, "pairwise_fst_weir_cockerham.tsv"), header = TRUE, row.names = 1, check.names = FALSE))
p_fst <- as.matrix(read.table(file.path(tab_dir, "pairwise_fst_permutation_pvalues.tsv"), header = TRUE, row.names = 1, check.names = FALSE))
fis_df <- read.table(file.path(tab_dir, "population_fis_permutation.tsv"), header = TRUE, sep = "\t", check.names = FALSE)

pops <- rownames(fst)
np <- length(pops)

stars_from_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

plot_matrix <- function(file, width = 9, height = 8.2) {
  if (grepl("\\.pdf$", file)) {
    pdf(file, width = width, height = height, useDingbats = FALSE)
  } else {
    png(file, width = width * 320, height = height * 320, res = 320)
  }

  par(mar = c(7.5, 8.5, 2.0, 4.6), oma = c(0.5, 0.5, 3.4, 0.5), xpd = NA)
  plot(NA, xlim = c(0.5, np + 0.5), ylim = c(0.5, np + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", asp = 1)

  fst_vals <- fst[upper.tri(fst)]
  fst_vals <- fst_vals[is.finite(fst_vals)]
  fst_min <- if (length(fst_vals)) min(fst_vals) else 0
  fst_max <- if (length(fst_vals)) max(fst_vals) else 0.2
  if (fst_max <= fst_min) fst_max <- fst_min + 1e-6

  pal_fst <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(120)
  get_fst_col <- function(v) {
    if (!is.finite(v)) return("#f0f0f0")
    t <- (v - fst_min) / (fst_max - fst_min)
    t <- min(max(t, 0), 1)
    pal_fst[1 + floor(t * 119)]
  }

  pal_p <- colorRampPalette(c("#ffffff", "#fdd0a2", "#e6550d"))(120)
  get_p_col <- function(p) {
    if (!is.finite(p)) return("#f0f0f0")
    v <- -log10(p + 1e-12)
    t <- min(v / 3, 1)
    pal_p[1 + floor(t * 119)]
  }

  for (i in seq_len(np)) {
    for (j in seq_len(np)) {
      xleft <- j - 0.5
      xright <- j + 0.5
      ybottom <- np - i + 0.5
      ytop <- np - i + 1.5

      if (i < j) {
        v <- fst[i, j]
        p <- p_fst[i, j]
        rect(xleft, ybottom, xright, ytop, col = get_fst_col(v), border = "grey70", lwd = 0.8)
        lab <- if (is.finite(v)) sprintf("%.3f\n%s", v, stars_from_p(p)) else "NA"
        text(j, np - i + 1, lab, cex = 0.85, font = 2)
      } else if (i > j) {
        p <- p_fst[i, j]
        rect(xleft, ybottom, xright, ytop, col = get_p_col(p), border = "grey70", lwd = 0.8)
        lab <- if (is.finite(p)) sprintf("p=%s\n%s", fmt_p(p), stars_from_p(p)) else "NA"
        text(j, np - i + 1, lab, cex = 0.72)
      } else {
        fp <- fis_df[fis_df$Population == pops[i], , drop = FALSE]
        fv <- fp$Fis[1]
        pv <- fp$p_value[1]
        fis_col <- if (!is.finite(fv)) "#f0f0f0" else if (fv >= 0) "#fee0d2" else "#deebf7"
        rect(xleft, ybottom, xright, ytop, col = fis_col, border = "grey50", lwd = 1.2)
        lab <- if (is.finite(fv)) sprintf("Fis=%.3f\n%s", fv, stars_from_p(pv)) else "Fis=NA"
        text(j, np - i + 1, lab, cex = 0.76, font = 2)
      }
    }
  }

  axis(1, at = seq_len(np), labels = pops, las = 2, tick = FALSE, cex.axis = 0.9)
  axis(2, at = np:1, labels = pops, las = 2, tick = FALSE, cex.axis = 0.9)

  # Titles inside full figure frame via outer margin
  mtext("Pairwise Fst and Fis Matrix", side = 3, outer = TRUE, line = 1.3, font = 2, cex = 1.0)
  mtext("Upper: Fst | Lower: p-value | Diagonal: Fis", side = 3, outer = TRUE, line = 0.4, cex = 0.85)

  legend("right", inset = c(-0.08, 0), bty = "n", cex = 0.82, y.intersp = 1.1,
         title = "p-value",
         legend = c("*** p<0.001", "** p<0.01", "* p<0.05", "ns p>=0.05"))

  dev.off()
}

plot_matrix(file.path(fig_dir, "Fst_Fis_two_sided_matrix_publication.pdf"))
plot_matrix(file.path(fig_dir, "Fst_Fis_two_sided_matrix_publication.png"))

cat("Redrew Fst/Fis matrix figure with title inside frame\n")
