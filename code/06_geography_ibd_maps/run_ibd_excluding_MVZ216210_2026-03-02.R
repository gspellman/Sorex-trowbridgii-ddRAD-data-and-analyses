#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
fig_dir <- file.path(down_dir, "figures")
tab_dir <- file.path(down_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

coord_file <- "/Users/gspellman/Trowbridgii_analyses/Sample_geographic_coordinants.txt"
popmap_file <- file.path(base_dir, "popmap.tsv")
vcf_file <- file.path(base_dir, "populations.snps.vcf")

out_stats <- file.path(tab_dir, "IBD_statistics_summary_excl_MVZ216210.tsv")
out_pairs <- file.path(tab_dir, "IBD_pairwise_distances_excl_MVZ216210.tsv")
out_bins <- file.path(tab_dir, "IBD_distance_bin_summary_excl_MVZ216210.tsv")
out_summary <- file.path(tab_dir, "IBD_analysis_summary_excl_MVZ216210.txt")
out_png <- file.path(fig_dir, "IBD_scatter_regression_excl_MVZ216210_publication.png")
out_pdf <- file.path(fig_dir, "IBD_scatter_regression_excl_MVZ216210_publication.pdf")
out_bin_png <- file.path(fig_dir, "IBD_distance_bin_trend_excl_MVZ216210_publication.png")
out_bin_pdf <- file.path(fig_dir, "IBD_distance_bin_trend_excl_MVZ216210_publication.pdf")

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
stop_if_missing <- function(x) {
  m <- x[!file.exists(x)]
  if (length(m) > 0) stop("Missing required files:\n", paste(m, collapse = "\n"))
}

parse_gt <- function(x) {
  gt <- sub(":.*$", "", x)
  gt <- gsub("\\|", "/", gt)
  out <- rep(NA_real_, length(gt))
  out[gt == "0/0"] <- 0
  out[gt %in% c("0/1", "1/0")] <- 1
  out[gt == "1/1"] <- 2
  out
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371.0088
  p1 <- lat1 * pi / 180
  p2 <- lat2 * pi / 180
  dphi <- (lat2 - lat1) * pi / 180
  dl <- (lon2 - lon1) * pi / 180
  a <- sin(dphi / 2)^2 + cos(p1) * cos(p2) * sin(dl / 2)^2
  2 * r * asin(sqrt(a))
}

mantel_test <- function(d1, d2, perms = 4999, seed = 20260302) {
  set.seed(seed)
  n <- nrow(d1)
  idx <- upper.tri(d1)
  x <- d1[idx]
  y <- d2[idx]
  r_obs <- suppressWarnings(cor(x, y, method = "pearson"))
  more <- 0L
  for (i in seq_len(perms)) {
    p <- sample.int(n)
    yp <- d2[p, p][idx]
    r <- suppressWarnings(cor(x, yp, method = "pearson"))
    if (is.finite(r) && abs(r) >= abs(r_obs)) more <- more + 1L
  }
  pval <- (more + 1) / (perms + 1)
  c(r = r_obs, p = pval)
}

stop_if_missing(c(coord_file, popmap_file, vcf_file))

logmsg("Loading coordinates and popmap")
coord <- read.delim(coord_file, header = TRUE, sep = "\t", check.names = FALSE)
pop <- read.delim(popmap_file, header = FALSE, sep = "\t", col.names = c("Sample", "Population"), check.names = FALSE)
coord$Sample <- as.character(coord$Sample)
pop$Sample <- as.character(pop$Sample)

if ("MVZ216210" %in% pop$Sample) stop("MVZ216210 found in excluded dataset popmap")

meta <- merge(pop, coord[, c("Sample", "Latitude", "Longitude")], by = "Sample", all = FALSE)
meta <- meta[order(meta$Population, meta$Sample), ]
meta$Latitude <- as.numeric(as.character(meta$Latitude))
meta$Longitude <- as.numeric(as.character(meta$Longitude))
if (any(!is.finite(meta$Latitude)) || any(!is.finite(meta$Longitude))) {
  stop("Non-finite coordinates detected after numeric conversion.")
}
if (min(meta$Latitude) < -90 || max(meta$Latitude) > 90 || min(meta$Longitude) < -180 || max(meta$Longitude) > 180) {
  stop("Coordinate ranges are invalid; check latitude/longitude parsing.")
}

sample_set <- unique(meta$Sample)

logmsg("Parsing VCF and extracting genotypes for matched samples")
lns <- readLines(vcf_file, warn = FALSE)
header_line <- lns[grep("^#CHROM", lns)]
if (length(header_line) != 1) stop("VCF #CHROM header not found")
h <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
all_samples <- h[10:length(h)]
keep_samples <- all_samples[all_samples %in% sample_set]
if (length(keep_samples) < 3) stop("Too few matched samples in VCF")
keep_idx <- match(keep_samples, all_samples)

variant_lines <- lns[!grepl("^#", lns)]
G_list <- vector("list", length(variant_lines))
keep_row <- logical(length(variant_lines))
ri <- 0L
for (k in seq_along(variant_lines)) {
  p <- strsplit(variant_lines[k], "\t", fixed = TRUE)[[1]]
  ref <- p[4]
  alt <- p[5]
  if (nchar(ref) != 1 || nchar(alt) != 1 || grepl(",", alt, fixed = TRUE)) next
  gfields <- p[10:length(p)]
  gt <- parse_gt(gfields[keep_idx])
  ri <- ri + 1L
  G_list[[ri]] <- gt
  keep_row[ri] <- TRUE
}
G_list <- G_list[seq_len(ri)]
G <- do.call(rbind, G_list)
colnames(G) <- keep_samples

# Align meta order to genotype order
meta2 <- meta[match(keep_samples, meta$Sample), ]
if (any(!is.finite(meta2$Latitude)) || any(!is.finite(meta2$Longitude))) {
  stop("Failed to align finite coordinates to genotype sample order.")
}

# Preflight geographic sanity check to prevent coordinate-parsing artifacts.
geo_check <- matrix(0, nrow(meta2), nrow(meta2))
for (ii in seq_len(nrow(meta2) - 1)) {
  for (jj in (ii + 1):nrow(meta2)) {
    dchk <- haversine_km(meta2$Latitude[ii], meta2$Longitude[ii], meta2$Latitude[jj], meta2$Longitude[jj])
    geo_check[ii, jj] <- dchk
    geo_check[jj, ii] <- dchk
  }
}
max_geo_check <- max(geo_check, na.rm = TRUE)
if (max_geo_check > 2000) {
  stop(sprintf("Geographic distance sanity check failed (max=%.3f km > 2000 km).", max_geo_check))
}

logmsg("Applying SNP filters (missing<=20%, MAF>=0.05)")
site_missing <- rowMeans(is.na(G))
p <- rowMeans(G, na.rm = TRUE) / 2
maf <- pmin(p, 1 - p)
keep <- is.finite(site_missing) & is.finite(maf) & site_missing <= 0.20 & maf >= 0.05
Gf <- G[keep, , drop = FALSE]
if (nrow(Gf) < 10) stop("Too few SNPs remain after filtering")

logmsg("Computing pairwise geographic and genetic distances")
n <- ncol(Gf)
geo <- matrix(0, n, n)
gen <- matrix(0, n, n)
pairs <- vector("list", n * (n - 1) / 2)
pair_i <- 0L
for (i in seq_len(n - 1)) {
  gi <- Gf[, i]
  for (j in (i + 1):n) {
    gj <- Gf[, j]
    ok <- is.finite(gi) & is.finite(gj)
    dgen <- if (any(ok)) mean(abs(gi[ok] - gj[ok]) / 2) else NA_real_
    dgeo <- haversine_km(meta2$Latitude[i], meta2$Longitude[i], meta2$Latitude[j], meta2$Longitude[j])
    gen[i, j] <- dgen
    gen[j, i] <- dgen
    geo[i, j] <- dgeo
    geo[j, i] <- dgeo

    pair_i <- pair_i + 1L
    pairs[[pair_i]] <- data.frame(
      Sample1 = keep_samples[i],
      Sample2 = keep_samples[j],
      Population1 = meta2$Population[i],
      Population2 = meta2$Population[j],
      Geographic_distance_km = dgeo,
      Genetic_distance = dgen,
      stringsAsFactors = FALSE
    )
  }
}

mgen <- median(gen[is.finite(gen)], na.rm = TRUE)
gen[!is.finite(gen)] <- mgen

idx <- upper.tri(gen)
x <- geo[idx]
y <- gen[idx]

logmsg("Running significance tests")
pear <- cor.test(x, y, method = "pearson")
spea <- cor.test(x, y, method = "spearman", exact = FALSE)
man <- mantel_test(gen, geo, perms = 4999, seed = 20260302)
lmfit <- lm(y ~ x)
lms <- summary(lmfit)

stats_df <- data.frame(
  Dataset = "stacks_refmap_sorex_excl_MVZ216210",
  Excluded_sample = "MVZ216210",
  N_samples = n,
  N_pairwise_comparisons = length(x),
  N_SNPs_after_filtering = nrow(Gf),
  Pearson_r = unname(pear$estimate),
  Pearson_p = pear$p.value,
  Spearman_rho = unname(spea$estimate),
  Spearman_p = spea$p.value,
  Mantel_r = unname(man["r"]),
  Mantel_p_permutation = unname(man["p"]),
  Mantel_permutations = 4999,
  Linear_slope = coef(lmfit)[2],
  Linear_intercept = coef(lmfit)[1],
  Linear_r2 = lms$r.squared,
  Linear_p = coef(lms)[2, 4],
  Linear_slope_stderr = coef(lms)[2, 2],
  check.names = FALSE
)
write.table(stats_df, out_stats, sep = "\t", row.names = FALSE, quote = FALSE)

pairs_df <- do.call(rbind, pairs)
write.table(pairs_df, out_pairs, sep = "\t", row.names = FALSE, quote = FALSE)

bins <- as.numeric(quantile(x, probs = seq(0, 1, 0.1), na.rm = TRUE))
bins <- unique(bins)
bin_id <- cut(x, breaks = bins, include.lowest = TRUE, labels = FALSE)
bin_rows <- list()
for (b in seq_len(length(bins) - 1)) {
  m <- which(bin_id == b)
  if (length(m) == 0) next
  yy <- y[m]
  bin_rows[[length(bin_rows) + 1L]] <- data.frame(
    Bin = b,
    Geo_km_min = bins[b],
    Geo_km_max = bins[b + 1],
    N_pairs = length(m),
    Mean_genetic_distance = mean(yy, na.rm = TRUE),
    SE_genetic_distance = if (length(yy) > 1) sd(yy, na.rm = TRUE) / sqrt(length(yy)) else NA_real_,
    check.names = FALSE
  )
}
bin_df <- do.call(rbind, bin_rows)
write.table(bin_df, out_bins, sep = "\t", row.names = FALSE, quote = FALSE)

logmsg("Rendering publication IBD figures")
plot_scatter <- function(outf, is_png = TRUE) {
  if (is_png) png(outf, width = 3200, height = 2500, res = 400) else pdf(outf, width = 8.0, height = 6.25, useDingbats = FALSE)
  par(mar = c(5.0, 5.2, 3.2, 1.2))
  plot(x, y, pch = 16, col = rgb(76/255, 120/255, 168/255, 0.40), cex = 0.75,
       xlab = "Geographic distance (km)",
       ylab = "Genetic distance (mean allele-dosage difference / 2)",
       main = "Isolation by distance (MVZ216210 excluded)")
  abline(lmfit, col = "#D62728", lwd = 2.3)
  grid(col = "#E0E0E0", lty = 3)
  legend("topleft", bty = "n", cex = 0.90, text.col = "#202020",
         legend = c(
           sprintf("Pearson r = %.3f (p = %.2e)", unname(pear$estimate), pear$p.value),
           sprintf("Spearman rho = %.3f (p = %.2e)", unname(spea$estimate), spea$p.value),
           sprintf("Mantel r = %.3f (p = %.4f, 4,999 perms)", unname(man["r"]), unname(man["p"]))
         ))
  dev.off()
}

plot_scatter(out_png, TRUE)
plot_scatter(out_pdf, FALSE)

plot_bins <- function(outf, is_png = TRUE) {
  if (is_png) png(outf, width = 3000, height = 2200, res = 400) else pdf(outf, width = 7.5, height = 5.5, useDingbats = FALSE)
  par(mar = c(5.0, 5.2, 3.2, 1.2))
  xc <- (bin_df$Geo_km_min + bin_df$Geo_km_max) / 2
  plot(xc, bin_df$Mean_genetic_distance, type = "b", pch = 16, cex = 0.9, lwd = 1.8,
       col = "#1F77B4",
       xlab = "Geographic distance bin midpoint (km)",
       ylab = "Mean genetic distance",
       main = "IBD distance-class trend (MVZ216210 excluded)")
  arrows(x0 = xc, y0 = bin_df$Mean_genetic_distance - bin_df$SE_genetic_distance,
         x1 = xc, y1 = bin_df$Mean_genetic_distance + bin_df$SE_genetic_distance,
         angle = 90, code = 3, length = 0.03, col = "#1F77B4", lwd = 1.3)
  grid(col = "#E0E0E0", lty = 3)
  dev.off()
}

plot_bins(out_bin_png, TRUE)
plot_bins(out_bin_pdf, FALSE)

writeLines(c(
  "Isolation-by-distance analysis summary (MVZ216210 excluded)",
  paste0("Samples analyzed: ", n),
  paste0("SNPs after filters: ", nrow(Gf)),
  paste0("Pairwise comparisons: ", length(x)),
  sprintf("Max pairwise geographic distance (sanity check): %.3f km", max_geo_check),
  sprintf("Pearson r = %.6f, p = %.6e", unname(pear$estimate), pear$p.value),
  sprintf("Spearman rho = %.6f, p = %.6e", unname(spea$estimate), spea$p.value),
  sprintf("Mantel r = %.6f, p = %.6f, permutations = 4999", unname(man["r"]), unname(man["p"])),
  paste0("Scatter figure: ", out_png),
  paste0("Distance-bin figure: ", out_bin_png),
  paste0("Stats table: ", out_stats),
  paste0("Pairwise table: ", out_pairs),
  paste0("Bin table: ", out_bins)
), con = out_summary)

logmsg("IBD analysis complete")
