#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
vcf_file <- file.path(base_dir, "populations.snps.vcf")
pop_file <- file.path(base_dir, "popmap.tsv")
tab_dir <- file.path(down_dir, "tables")
fig_dir <- file.path(down_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

is_missing_gt <- function(x) {
  gt <- sub(":.*", "", x)
  gt <- gsub("\\|", "/", gt)
  gt %in% c("./.", ".")
}

con <- file(vcf_file, "r")
on.exit(close(con), add = TRUE)

samples <- NULL
missing_counts <- NULL
total_sites <- 0L

repeat {
  ln <- readLines(con, n = 1)
  if (!length(ln)) break
  if (startsWith(ln, "#CHROM")) {
    hdr <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    samples <- hdr[10:length(hdr)]
    missing_counts <- integer(length(samples))
    next
  }
  if (startsWith(ln, "#")) next
  p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  if (length(p) < 10) next
  g <- p[10:length(p)]
  total_sites <- total_sites + 1L
  miss <- vapply(g, is_missing_gt, logical(1))
  missing_counts <- missing_counts + as.integer(miss)
}

if (is.null(samples) || total_sites == 0) stop("Failed to parse VCF for samples/sites")

miss_prop <- missing_counts / total_sites
callrate <- 1 - miss_prop

popmap <- read.table(pop_file, header = TRUE, sep = "\t", check.names = FALSE)
colnames(popmap) <- c("Sample", "Population")

df <- data.frame(
  Sample = samples,
  N_sites = total_sites,
  N_missing = missing_counts,
  Missingness = miss_prop,
  Callrate = callrate,
  stringsAsFactors = FALSE
)
df <- merge(df, popmap, by = "Sample", all.x = TRUE)

q1 <- as.numeric(stats::quantile(df$Missingness, 0.25, na.rm = TRUE))
q3 <- as.numeric(stats::quantile(df$Missingness, 0.75, na.rm = TRUE))
iqr <- q3 - q1
lower <- max(0, q1 - 1.5 * iqr)
upper <- q3 + 1.5 * iqr

df$Outlier_IQR <- (df$Missingness < lower) | (df$Missingness > upper)
mu <- mean(df$Missingness, na.rm = TRUE)
sd_m <- stats::sd(df$Missingness, na.rm = TRUE)
if (is.finite(sd_m) && sd_m > 0) {
  df$Missingness_Z <- (df$Missingness - mu) / sd_m
} else {
  df$Missingness_Z <- 0
}
df$Outlier_Z3 <- abs(df$Missingness_Z) > 3

df <- df[order(df$Missingness, decreasing = TRUE), , drop = FALSE]

full_tab <- file.path(tab_dir, "sample_missingness_full_dataset.tsv")
out_tab <- file.path(tab_dir, "sample_missingness_outliers.tsv")
sum_txt <- file.path(tab_dir, "sample_missingness_summary.txt")
write.table(df, full_tab, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df[df$Outlier_IQR | df$Outlier_Z3, , drop = FALSE], out_tab, sep = "\t", quote = FALSE, row.names = FALSE)

cat(
  "Per-sample missingness summary (full dataset)\n",
  sprintf("VCF: %s\n", vcf_file),
  sprintf("Samples: %d\n", nrow(df)),
  sprintf("Sites: %d\n", total_sites),
  sprintf("Mean missingness: %.6f\n", mean(df$Missingness)),
  sprintf("Median missingness: %.6f\n", median(df$Missingness)),
  sprintf("Q1: %.6f\n", q1),
  sprintf("Q3: %.6f\n", q3),
  sprintf("IQR: %.6f\n", iqr),
  sprintf("IQR lower bound: %.6f\n", lower),
  sprintf("IQR upper bound: %.6f\n", upper),
  sprintf("IQR outliers: %d\n", sum(df$Outlier_IQR)),
  sprintf("|Z| > 3 outliers: %d\n", sum(df$Outlier_Z3)),
  file = sum_txt
)

hist_png <- file.path(fig_dir, "sample_missingness_histogram_full_dataset.png")
hist_pdf <- file.path(fig_dir, "sample_missingness_histogram_full_dataset.pdf")

plot_hist <- function(dev_fun, file) {
  dev_fun(file)
  par(mar = c(5, 5, 3, 1))
  vals <- 100 * df$Missingness
  h <- hist(vals, breaks = 20, col = "#4C78A8", border = "white",
            main = "Per-sample missingness distribution (full SNP dataset)",
            xlab = "Per-sample missingness (%)", ylab = "Number of samples")
  abline(v = 100 * upper, col = "#D62728", lwd = 2, lty = 2)
  abline(v = 100 * median(df$Missingness), col = "#2F2F2F", lwd = 1.5)
  if (sum(df$Outlier_IQR) > 0) {
    for (x in 100 * df$Missingness[df$Outlier_IQR]) {
      abline(v = x, col = "#E45756", lwd = 1, lty = 1)
    }
  }
  legend("topright",
         legend = c(sprintf("Median = %.2f%%", 100 * median(df$Missingness)),
                    sprintf("IQR upper = %.2f%%", 100 * upper)),
         lty = c(1, 2), lwd = c(1.5, 2), col = c("#2F2F2F", "#D62728"), bty = "n")
  dev.off()
}

plot_hist(function(f) png(f, width = 3000, height = 2000, res = 320), hist_png)
plot_hist(function(f) pdf(f, width = 9.5, height = 6.2, useDingbats = FALSE), hist_pdf)

cat(sprintf("Wrote:\n%s\n%s\n%s\n%s\n%s\n", full_tab, out_tab, sum_txt, hist_png, hist_pdf))
