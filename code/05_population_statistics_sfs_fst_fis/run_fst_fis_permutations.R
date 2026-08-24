#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

stacks_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
out_dir <- file.path(stacks_dir, "downstream_analysis_2026-02-23")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(stacks_dir, "populations.snps.vcf")
pop_file <- file.path(stacks_dir, "popmap.tsv")

# Publication palette provided by user
pop_cols <- c(
  North = "#4169E1",
  North_Coast = "#4CBB17",
  Sierra_1 = "#800080",
  Sierra_2 = "#00FFFF",
  Sierra_3 = "#FF0000",
  South_Coast = "#000000"
)
pop_order <- names(pop_cols)

# Permutation settings
n_perm_fst <- 199
n_perm_fis <- 199
set.seed(20260223)

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

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

# Weir & Cockerham pairwise theta estimator (biallelic)
calc_pair_wc_fst <- function(geno_mat, idx1, idx2) {
  g1 <- geno_mat[, idx1, drop = FALSE]
  g2 <- geno_mat[, idx2, drop = FALSE]

  n1 <- rowSums(!is.na(g1))
  n2 <- rowSums(!is.na(g2))
  valid <- (n1 >= 2) & (n2 >= 2)
  if (!any(valid)) return(NA_real_)

  p1 <- rowSums(g1, na.rm = TRUE) / (2 * n1)
  p2 <- rowSums(g2, na.rm = TRUE) / (2 * n2)
  h1 <- rowSums(g1 == 1, na.rm = TRUE) / n1
  h2 <- rowSums(g2 == 1, na.rm = TRUE) / n2

  nbar <- (n1 + n2) / 2
  nc <- (2 * nbar - (n1^2 + n2^2) / (2 * nbar))
  pbar <- (n1 * p1 + n2 * p2) / (2 * nbar)
  hbar <- (n1 * h1 + n2 * h2) / (2 * nbar)
  s2 <- (n1 * (p1 - pbar)^2 + n2 * (p2 - pbar)^2) / nbar

  ok <- valid & is.finite(pbar) & is.finite(hbar) & is.finite(s2) & (nbar > 1) & (nc > 0)
  if (!any(ok)) return(NA_real_)

  a <- (nbar / nc) * (s2 - (1 / (nbar - 1)) * (pbar * (1 - pbar) - 0.5 * s2 - hbar / 4))
  b <- (nbar / (nbar - 1)) * (pbar * (1 - pbar) - 0.5 * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
  c <- hbar / 2

  num <- sum(a[ok], na.rm = TRUE)
  den <- sum((a + b + c)[ok], na.rm = TRUE)
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

# Population-level Fis (heterozygote deficit estimator)
calc_pop_fis <- function(geno_mat_pop) {
  n <- rowSums(!is.na(geno_mat_pop))
  k <- rowSums(geno_mat_pop, na.rm = TRUE)
  ho_count <- rowSums(geno_mat_pop == 1, na.rm = TRUE)

  p <- k / (2 * n)
  he <- 2 * p * (1 - p)
  valid <- (n >= 2) & is.finite(he) & (he > 0)
  if (!any(valid)) return(NA_real_)

  hobs_total <- sum(ho_count[valid], na.rm = TRUE)
  hexp_total <- sum(he[valid] * n[valid], na.rm = TRUE)
  if (!is.finite(hexp_total) || hexp_total <= 0) return(NA_real_)

  1 - (hobs_total / hexp_total)
}

# Allele permutation test for Fis within a population
perm_test_fis <- function(geno_mat_pop, n_perm = 199) {
  n <- rowSums(!is.na(geno_mat_pop))
  k <- rowSums(geno_mat_pop, na.rm = TRUE)
  ho_count <- rowSums(geno_mat_pop == 1, na.rm = TRUE)
  p <- k / (2 * n)
  he <- 2 * p * (1 - p)
  valid <- (n >= 2) & is.finite(he) & (he > 0)

  if (!any(valid)) return(list(obs = NA_real_, p = NA_real_))

  n <- n[valid]
  k <- k[valid]
  ho_count <- ho_count[valid]
  he <- he[valid]

  obs_fis <- 1 - sum(ho_count) / sum(he * n)

  perm_fis <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    het_total <- 0
    for (i in seq_along(n)) {
      ni <- n[i]
      ki <- k[i]
      alleles <- c(rep.int(1L, ki), rep.int(0L, 2L * ni - ki))
      shuf <- sample(alleles, length(alleles), replace = FALSE)
      geno2 <- matrix(shuf, ncol = 2, byrow = TRUE)
      het_total <- het_total + sum(rowSums(geno2) == 1L)
    }
    perm_fis[b] <- 1 - het_total / sum(he * n)
  }

  p_two <- (1 + sum(abs(perm_fis) >= abs(obs_fis))) / (n_perm + 1)
  list(obs = obs_fis, p = p_two)
}

#-------------------------
# Read metadata and VCF
#-------------------------
logmsg("Loading population map")
popmap <- read.table(pop_file, header = FALSE, sep = "\t", col.names = c("Sample", "Pop"))
popmap$Pop <- factor(popmap$Pop, levels = pop_order)

logmsg("Reading VCF header")
header_lines <- readLines(vcf_file, n = 5000)
chrom_line <- header_lines[grepl("^#CHROM", header_lines)]
if (length(chrom_line) != 1) stop("Could not find #CHROM line in VCF")
samples <- strsplit(chrom_line, "\t")[[1]][10:length(strsplit(chrom_line, "\t")[[1]])]

logmsg("Reading VCF body")
vcf <- read.table(vcf_file, sep = "\t", comment.char = "#", header = FALSE,
                  quote = "", fill = TRUE)
colnames(vcf)[1:9] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
colnames(vcf)[10:ncol(vcf)] <- samples

# Reorder popmap to VCF sample order
popmap <- popmap[match(samples, popmap$Sample), ]
if (any(is.na(popmap$Sample))) stop("Some VCF samples missing in popmap")

#-------------------------
# ddRAD-appropriate filtering
#-------------------------
logmsg("Applying standard ddRAD SNP filters")
# Keep SNPs that are biallelic and single nucleotide
is_snp <- nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1
is_bial <- !grepl(",", vcf$ALT, fixed = TRUE)
keep_sites <- is_snp & is_bial

vcf2 <- vcf[keep_sites, c("ID", "REF", "ALT", samples)]

# Convert GT to 0/1/2 dosage
as_dosage <- function(x) {
  gt <- sub(":.*", "", x)
  gt <- gsub("\\|", "/", gt)
  out <- rep(NA_real_, length(gt))
  out[gt == "0/0"] <- 0
  out[gt %in% c("0/1", "1/0")] <- 1
  out[gt == "1/1"] <- 2
  out
}

logmsg("Converting genotypes to dosage matrix")
geno <- sapply(samples, function(s) as_dosage(vcf2[[s]]))
if (!is.matrix(geno)) geno <- as.matrix(geno)
rownames(geno) <- vcf2$ID
colnames(geno) <- samples

# Missingness and MAF filters
site_missing <- rowMeans(is.na(geno))
p_all <- rowMeans(geno, na.rm = TRUE) / 2
maf <- pmin(p_all, 1 - p_all)

keep <- (site_missing <= 0.20) & is.finite(maf) & (maf >= 0.05)
geno_f <- geno[keep, , drop = FALSE]

logmsg("Retained loci after filters: ", nrow(geno_f))
if (nrow(geno_f) < 10) stop("Too few loci after filtering")

#-------------------------
# Pairwise Fst + permutation
#-------------------------
pops <- pop_order[pop_order %in% unique(as.character(popmap$Pop))]
np <- length(pops)

fst_mat <- matrix(NA_real_, np, np, dimnames = list(pops, pops))
p_fst_mat <- matrix(NA_real_, np, np, dimnames = list(pops, pops))

logmsg("Computing pairwise Weir-Cockerham Fst + permutation tests")
for (i in seq_len(np - 1)) {
  for (j in (i + 1):np) {
    p1 <- pops[i]; p2 <- pops[j]
    idx1 <- which(popmap$Pop == p1)
    idx2 <- which(popmap$Pop == p2)

    sub_idx <- c(idx1, idx2)
    sub_geno <- geno_f[, sub_idx, drop = FALSE]

    obs <- calc_pair_wc_fst(sub_geno, seq_along(idx1), length(idx1) + seq_along(idx2))

    perm <- rep(NA_real_, n_perm_fst)
    n1 <- length(idx1)
    n_sub <- ncol(sub_geno)

    for (b in seq_len(n_perm_fst)) {
      ord <- sample.int(n_sub, n_sub, replace = FALSE)
      g1 <- ord[1:n1]
      g2 <- ord[(n1 + 1):n_sub]
      perm[b] <- calc_pair_wc_fst(sub_geno, g1, g2)
    }

    p_two <- (1 + sum(abs(perm) >= abs(obs), na.rm = TRUE)) / (n_perm_fst + 1)

    fst_mat[i, j] <- fst_mat[j, i] <- obs
    p_fst_mat[i, j] <- p_fst_mat[j, i] <- p_two

    logmsg("  ", p1, " vs ", p2, ": Fst=", sprintf("%.4f", obs), ", p=", fmt_p(p_two))
  }
}

diag(fst_mat) <- NA_real_
diag(p_fst_mat) <- NA_real_

#-------------------------
# Population Fis + allele permutation
#-------------------------
logmsg("Computing population Fis + allele permutation tests")
fis_df <- data.frame(Population = pops, Fis = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)

for (i in seq_along(pops)) {
  p <- pops[i]
  idx <- which(popmap$Pop == p)
  g <- geno_f[, idx, drop = FALSE]
  res <- perm_test_fis(g, n_perm = n_perm_fis)
  fis_df$Fis[i] <- res$obs
  fis_df$p_value[i] <- res$p
  logmsg("  ", p, ": Fis=", sprintf("%.4f", res$obs), ", p=", fmt_p(res$p))
}

#-------------------------
# Save tables
#-------------------------
write.table(round(fst_mat, 6), file.path(tab_dir, "pairwise_fst_weir_cockerham.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(round(p_fst_mat, 6), file.path(tab_dir, "pairwise_fst_permutation_pvalues.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(fis_df, file.path(tab_dir, "population_fis_permutation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------
# Publication-style two-sided matrix plot
# Upper triangle: Fst estimates (+ significance)
# Lower triangle: Fst permutation p-values
# Diagonal: Fis per population (+ significance)
#-------------------------
plot_matrix <- function(file, width = 9, height = 8.2) {
  if (grepl("\\.pdf$", file)) {
    pdf(file, width = width, height = height, useDingbats = FALSE)
  } else {
    png(file, width = width * 320, height = height * 320, res = 320)
  }

  par(mar = c(7.5, 8.5, 3.5, 2), xpd = NA)
  plot(NA, xlim = c(0.5, np + 0.5), ylim = c(0.5, np + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", asp = 1)

  fst_vals <- fst_mat[upper.tri(fst_mat)]
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

  # Draw cells with annotations
  for (i in seq_len(np)) {
    for (j in seq_len(np)) {
      xleft <- j - 0.5
      xright <- j + 0.5
      ybottom <- np - i + 0.5
      ytop <- np - i + 1.5

      if (i < j) {
        v <- fst_mat[i, j]
        p <- p_fst_mat[i, j]
        rect(xleft, ybottom, xright, ytop, col = get_fst_col(v), border = "grey70", lwd = 0.8)
        lab <- if (is.finite(v)) sprintf("%.3f\n%s", v, stars_from_p(p)) else "NA"
        text(j, np - i + 1, lab, cex = 0.85, font = 2)
      } else if (i > j) {
        p <- p_fst_mat[i, j]
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

  mtext("Two-Sided Matrix: Pairwise Fst (Upper) and Significance p-values (Lower); Diagonal = Population Fis",
        side = 3, line = 1.2, font = 2, cex = 1.0)
  mtext("Weir-Cockerham estimator on filtered biallelic ddRAD SNPs; permutation p-values", side = 3, line = 0.2, cex = 0.85)

  legend("topright", inset = c(0, -0.18), horiz = TRUE, bty = "n", cex = 0.82,
         legend = c("*** p<0.001", "** p<0.01", "* p<0.05", "ns p>=0.05"))

  dev.off()
}

plot_matrix(file.path(fig_dir, "Fst_Fis_two_sided_matrix_publication.pdf"))
plot_matrix(file.path(fig_dir, "Fst_Fis_two_sided_matrix_publication.png"))

summary_txt <- c(
  "Estimator:",
  "  Pairwise Fst: Weir & Cockerham theta (1984), ratio-of-components across loci.",
  "  Fis: population-level heterozygote deficit estimator (1 - Hobs/Hexp).",
  "Filters:",
  "  Biallelic SNPs only (single nucleotide REF/ALT)",
  "  Site missingness <= 20%",
  "  Minor allele frequency >= 0.05",
  paste0("Permutation tests:"),
  paste0("  Pairwise Fst permutations per comparison: ", n_perm_fst),
  paste0("  Fis allele-pairing permutations per population: ", n_perm_fis),
  "Outputs:",
  paste0("  ", file.path(tab_dir, "pairwise_fst_weir_cockerham.tsv")),
  paste0("  ", file.path(tab_dir, "pairwise_fst_permutation_pvalues.tsv")),
  paste0("  ", file.path(tab_dir, "population_fis_permutation.tsv")),
  paste0("  ", file.path(fig_dir, "Fst_Fis_two_sided_matrix_publication.pdf")),
  paste0("  ", file.path(fig_dir, "Fst_Fis_two_sided_matrix_publication.png"))
)
writeLines(summary_txt, con = file.path(tab_dir, "fst_fis_analysis_summary.txt"))

logmsg("Finished Fst/Fis analysis")
