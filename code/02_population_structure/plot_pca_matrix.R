#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

out_root <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23"
fig_dir <- file.path(out_root, "figures")
tab_dir <- file.path(out_root, "tables")
log_dir <- file.path(out_root, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stacks_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
vcf_file <- file.path(stacks_dir, "populations.snps.vcf")
pop_file <- file.path(stacks_dir, "popmap.tsv")

pop_cols <- c(
  North = "#4169E1",
  North_Coast = "#4CBB17",
  Sierra_1 = "#800080",
  Sierra_2 = "#00FFFF",
  Sierra_3 = "#FF0000",
  South_Coast = "#000000"
)
pop_order <- names(pop_cols)

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

logmsg("Loading popmap")
popmap <- read.table(pop_file, header = FALSE, sep = "\t", col.names = c("Sample", "Pop"))
popmap$Pop <- factor(popmap$Pop, levels = pop_order)

logmsg("Reading VCF header/body")
header_lines <- readLines(vcf_file, n = 5000)
chrom_line <- header_lines[grepl("^#CHROM", header_lines)]
if (length(chrom_line) != 1) stop("Could not find #CHROM header in VCF.")
samples <- strsplit(chrom_line, "\t")[[1]][10:length(strsplit(chrom_line, "\t")[[1]])]

vcf <- read.table(vcf_file, sep = "\t", comment.char = "#", header = FALSE,
                  quote = "", fill = TRUE)
colnames(vcf)[1:9] <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
colnames(vcf)[10:ncol(vcf)] <- samples

popmap <- popmap[match(samples, popmap$Sample), ]
if (any(is.na(popmap$Sample))) stop("Some VCF samples are missing from popmap.")

as_gt <- function(x) {
  x <- sub(":.*", "", x)
  x[x %in% c("./.", ".|.", ".")] <- NA
  gsub("|", "/", x, fixed = TRUE)
}

gt_to_nonref <- function(x) {
  x <- as_gt(x)
  a1 <- sub("/.*", "", x)
  a2 <- sub(".*/", "", x)
  good <- !(is.na(x) | a1 == "." | a2 == ".")
  out <- rep(NA_real_, length(x))
  out[good] <- as.numeric(a1[good] != "0") + as.numeric(a2[good] != "0")
  out
}

gt_to_bial <- function(x) {
  x <- as_gt(x)
  a1 <- sub("/.*", "", x)
  a2 <- sub(".*/", "", x)
  good <- !(is.na(x) | a1 == "." | a2 == ".")
  valid <- good & a1 %in% c("0", "1") & a2 %in% c("0", "1")
  out <- rep(NA_real_, length(x))
  out[valid] <- as.numeric(a1[valid] == "1") + as.numeric(a2[valid] == "1")
  out
}

prep_matrix <- function(x) {
  keep <- colSums(!is.na(x)) > 0
  x <- x[, keep, drop = FALSE]
  cm <- colMeans(x, na.rm = TRUE)
  na_idx <- which(is.na(x), arr.ind = TRUE)
  if (nrow(na_idx) > 0) x[na_idx] <- cm[na_idx[,2]]
  v <- apply(x, 2, var)
  x[, v > 0, drop = FALSE]
}

logmsg("Converting all informative and biallelic genotypes")
gt_mat <- as.matrix(vcf[, samples])
all_dos <- apply(gt_mat, 2, gt_to_nonref)
rownames(all_dos) <- vcf$ID
colnames(all_dos) <- samples
X_all <- prep_matrix(t(all_dos))

bial_idx <- !grepl(",", vcf$ALT, fixed = TRUE)
gt_bial <- gt_mat[bial_idx, , drop = FALSE]
ids_bial <- vcf$ID[bial_idx]
bial_dos <- apply(gt_bial, 2, gt_to_bial)
rownames(bial_dos) <- ids_bial
colnames(bial_dos) <- samples
X_bial <- prep_matrix(t(bial_dos))

logmsg("Running PCA")
pca_all <- prcomp(X_all, center = TRUE, scale. = TRUE)
var_all <- 100 * (pca_all$sdev^2 / sum(pca_all$sdev^2))
scores_all <- pca_all$x[, 1:4, drop = FALSE]

pca_bial <- prcomp(X_bial, center = TRUE, scale. = TRUE)
var_bial <- 100 * (pca_bial$sdev^2 / sum(pca_bial$sdev^2))
scores_bial <- pca_bial$x[, 1:4, drop = FALSE]

write.table(
  data.frame(PC = paste0("PC", 1:10),
             all_informative_percent = round(var_all[1:10], 4),
             biallelic_percent = round(var_bial[1:10], 4)),
  file.path(tab_dir, "pca_variance_explained.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

point_cex <- 1.45
legend_cex <- 0.95
legend_pt_cex <- 1.7

plot_pca_matrix <- function(file, width, height) {
  if (grepl("\\.pdf$", file)) {
    pdf(file, width = width, height = height, useDingbats = FALSE)
  } else {
    png(file, width = width * 320, height = height * 320, res = 320)
  }

  par(mfrow = c(4, 4), mar = c(2.9, 2.9, 1.2, 0.7), oma = c(1.2, 1.2, 3.4, 0.8))

  for (i in 1:4) {
    for (j in 1:4) {
      if (i == j) {
        plot.new()
        text(0.5, 0.64, paste0("PC", i), cex = 1.25, font = 2)
        text(0.5, 0.45, paste0("Bial: ", sprintf("%.2f", var_bial[i]), "%"), cex = 0.95)
        text(0.5, 0.30, paste0("All:  ", sprintf("%.2f", var_all[i]), "%"), cex = 0.95)
      } else {
        if (i < j) {
          sc <- scores_bial
          vv <- var_bial
          subttl <- "Biallelic"
        } else {
          sc <- scores_all
          vv <- var_all
          subttl <- "All informative"
        }

        plot(sc[, j], sc[, i],
             col = pop_cols[as.character(popmap$Pop)],
             pch = 16,
             cex = point_cex,
             xlab = paste0("PC", j, " (", sprintf("%.2f", vv[j]), "%)"),
             ylab = paste0("PC", i, " (", sprintf("%.2f", vv[i]), "%)"),
             main = subttl)
      }
    }
  }

  mtext("PCA Matrix: Upper Triangle = Biallelic Sites, Lower Triangle = All Informative Sites",
        side = 3, outer = TRUE, line = 1.8, cex = 1.1, font = 2)

  # Draw figure-level legend inside the device panel (top-right corner).
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE, xpd = FALSE)
  plot.new()
  legend(x = 0.985, y = 0.985,
         legend = pop_order,
         col = pop_cols[pop_order],
         pch = 16,
         pt.cex = legend_pt_cex,
         bty = "o",
         box.col = "grey50",
         bg = grDevices::adjustcolor("white", alpha.f = 0.8),
         cex = legend_cex,
         title = "Population",
         xjust = 1, yjust = 1)

  dev.off()
}

logmsg("Writing updated PCA matrix figures")
plot_pca_matrix(file.path(fig_dir, "PCA_matrix_first4PCs_publication.pdf"), 14, 14)
plot_pca_matrix(file.path(fig_dir, "PCA_matrix_first4PCs_publication.png"), 14, 14)

logmsg("Done")
