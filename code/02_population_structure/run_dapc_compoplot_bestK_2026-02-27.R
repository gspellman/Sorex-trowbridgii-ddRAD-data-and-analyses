#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

base <- '/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean'
out_root <- file.path(base, 'downstream_analysis_2026-02-23')
out_fig <- file.path(out_root, 'figures')
out_tab <- file.path(out_root, 'tables')
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(base, 'populations.snps.vcf')
pop_file <- file.path(base, 'popmap.tsv')

if (!file.exists(vcf_file)) stop('Missing VCF: ', vcf_file)
if (!file.exists(pop_file)) stop('Missing popmap: ', pop_file)

pop_cols <- c(
  North = '#4169E1',
  North_Coast = '#4CBB17',
  Sierra_1 = '#800080',
  Sierra_2 = '#00FFFF',
  Sierra_3 = '#FF0000',
  South_Coast = '#000000'
)
pop_order <- names(pop_cols)

as_gt <- function(x) {
  x <- sub(':.*', '', x)
  x[x %in% c('./.', '.|.', '.')] <- NA
  gsub('|', '/', x, fixed = TRUE)
}

gt_to_bial <- function(x) {
  x <- as_gt(x)
  a1 <- sub('/.*', '', x)
  a2 <- sub('.*/', '', x)
  good <- !(is.na(x) | a1 == '.' | a2 == '.')
  valid <- good & a1 %in% c('0', '1') & a2 %in% c('0', '1')
  out <- rep(NA_real_, length(x))
  out[valid] <- as.numeric(a1[valid] == '1') + as.numeric(a2[valid] == '1')
  out
}

prep_matrix <- function(x) {
  keep <- colSums(!is.na(x)) > 0
  x <- x[, keep, drop = FALSE]
  cm <- colMeans(x, na.rm = TRUE)
  na_idx <- which(is.na(x), arr.ind = TRUE)
  if (nrow(na_idx) > 0) x[na_idx] <- cm[na_idx[, 2]]
  v <- apply(x, 2, var)
  x <- x[, v > 0, drop = FALSE]
  x
}

# Read VCF sample header
hdr <- readLines(vcf_file, n = 5000)
chrom_line <- hdr[grepl('^#CHROM', hdr)]
if (length(chrom_line) != 1) stop('Could not parse #CHROM header in VCF.')
samples <- strsplit(chrom_line, '\t')[[1]][10:length(strsplit(chrom_line, '\t')[[1]])]

# Read VCF body
vcf <- read.table(vcf_file, sep = '\t', comment.char = '#', header = FALSE, quote = '', fill = TRUE)
colnames(vcf)[1:9] <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
colnames(vcf)[10:ncol(vcf)] <- samples

# Biallelic matrix
bial_idx <- !grepl(',', vcf$ALT, fixed = TRUE)
gt_bial <- as.matrix(vcf[bial_idx, samples, drop = FALSE])
bial_dos <- apply(gt_bial, 2, gt_to_bial)
colnames(bial_dos) <- samples
X <- prep_matrix(t(bial_dos))   # samples x loci

# Load population map for ordering only
popmap <- read.table(pop_file, header = FALSE, sep = '\t', col.names = c('Sample', 'Pop'))
popmap$Pop <- factor(popmap$Pop, levels = pop_order)
pop_for_samples <- popmap$Pop[match(rownames(X), popmap$Sample)]

# PCA and retain best_npc from prior DAPC-like CV (or fallback)
cv_file <- file.path(out_tab, 'dapc_like_cv_accuracy.tsv')
best_npc <- 20
if (file.exists(cv_file)) {
  cv <- read.table(cv_file, header = TRUE, sep = '\t')
  if (all(c('n_pcs','cv_accuracy') %in% colnames(cv))) {
    best_npc <- cv$n_pcs[which.max(cv$cv_accuracy)][1]
  }
}

pca <- prcomp(X, center = TRUE, scale. = TRUE)
npc <- min(best_npc, ncol(pca$x), nrow(pca$x) - 1)
Z <- pca$x[, seq_len(npc), drop = FALSE]

# DAPC-style K selection by BIC over kmeans on retained PCs
bic_from_kmeans <- function(km, Z) {
  n <- nrow(Z)
  k <- nrow(km$centers)
  wss <- sum(km$withinss)
  if (wss <= 0) wss <- .Machine$double.eps
  # Pseudo-BIC used in our clustering workflows; less over-penalizing than full parameter-count BIC
  n * log((wss / n) + 1e-12) + (k - 1) * log(max(n, 2))
}

K_grid <- 1:10
set.seed(20260227)
km_list <- lapply(K_grid, function(k) kmeans(Z, centers = k, nstart = 50, iter.max = 200))
bic_vals <- sapply(seq_along(K_grid), function(i) bic_from_kmeans(km_list[[i]], Z))

# Average silhouette (base-R implementation, no extra packages)
dist_mat <- as.matrix(dist(Z))
avg_silhouette <- function(cluster) {
  cluster <- as.integer(cluster)
  n <- length(cluster)
  u <- sort(unique(cluster))
  if (length(u) < 2) return(NA_real_)
  s <- numeric(n)
  for (i in seq_len(n)) {
    ci <- cluster[i]
    same <- which(cluster == ci)
    same <- same[same != i]
    a_i <- if (length(same) > 0) mean(dist_mat[i, same]) else 0
    b_i <- Inf
    for (cj in u[u != ci]) {
      other <- which(cluster == cj)
      if (length(other) > 0) {
        d <- mean(dist_mat[i, other])
        if (d < b_i) b_i <- d
      }
    }
    den <- max(a_i, b_i)
    s[i] <- ifelse(is.finite(den) && den > 0, (b_i - a_i) / den, 0)
  }
  mean(s)
}
sil_vals <- sapply(km_list, function(km) avg_silhouette(km$cluster))

# Use silhouette as primary best-K criterion; retain BIC as secondary diagnostic
if (all(is.na(sil_vals))) {
  best_k <- K_grid[which.min(bic_vals)][1]
} else {
  k2 <- K_grid[K_grid >= 2]
  s2 <- sil_vals[K_grid >= 2]
  best_k <- k2[which.max(s2)][1]
}
best_km <- km_list[[which(K_grid == best_k)[1]]]

# DAPC-like soft membership from centroid distances in retained PCA space
centers <- best_km$centers
if (best_k == 1) centers <- matrix(centers, nrow = 1)
cl <- best_km$cluster

# Cluster-specific scale from within-cluster mean squared distance
sigma2 <- rep(1, best_k)
for (k in seq_len(best_k)) {
  idxk <- which(cl == k)
  if (length(idxk) > 1) {
    d2k <- rowSums((Z[idxk, , drop = FALSE] - matrix(centers[k, ], nrow = length(idxk), ncol = ncol(Z), byrow = TRUE))^2)
    s2 <- mean(d2k)
    sigma2[k] <- ifelse(is.finite(s2) && s2 > 1e-8, s2, 1)
  }
}

d2 <- matrix(0, nrow = nrow(Z), ncol = best_k)
for (k in seq_len(best_k)) {
  d2[, k] <- rowSums((Z - matrix(centers[k, ], nrow = nrow(Z), ncol = ncol(Z), byrow = TRUE))^2)
}
logw <- sweep(-0.5 * d2, 2, sigma2, FUN = "/")
logw <- sweep(logw, 1, apply(logw, 1, max), FUN = "-")
w <- exp(logw)
post <- w / rowSums(w)
colnames(post) <- paste0('Cluster', seq_len(ncol(post)))

# Save tables
k_sel <- data.frame(K = K_grid, BIC = bic_vals, AvgSilhouette = sil_vals)
write.table(k_sel, file.path(out_tab, 'dapc_compoplot_k_bic.tsv'), sep='\t', row.names = FALSE, quote = FALSE)
write.table(data.frame(best_k = best_k, n_pcs = npc, selection = 'max_average_silhouette'),
            file.path(out_tab, 'dapc_compoplot_best_k.tsv'), sep='\t', row.names = FALSE, quote = FALSE)

qtab <- data.frame(Sample = rownames(Z), Population = as.character(pop_for_samples), post, check.names = FALSE)
write.table(qtab, file.path(out_tab, 'dapc_compoplot_membership_bestK.tsv'), sep='\t', row.names = FALSE, quote = FALSE)

# Order samples by population then max-cluster assignment
pop_chr <- as.character(pop_for_samples)
max_cl <- apply(post, 1, which.max)
ord <- order(factor(pop_chr, levels = pop_order), max_cl, rownames(Z))
post_o <- post[ord, , drop = FALSE]
pop_o <- pop_chr[ord]

# Population axis helpers
idx <- seq_len(length(pop_o))
split_idx <- split(idx, factor(pop_o, levels = pop_order))
split_idx <- split_idx[sapply(split_idx, length) > 0]
centers <- sapply(split_idx, function(v) mean(range(v)))
labels <- names(centers)
seps <- cumsum(sapply(split_idx, length)) + 0.5
seps <- seps[seps < length(pop_o)]

# Distinct cluster palette
pal <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#1f78b4','#b15928','#17becf')[seq_len(ncol(post_o))]

plot_comp <- function(file, w, h, png_dev = FALSE) {
  if (png_dev) png(file, width = w*300, height = h*300, res = 300) else pdf(file, width = w, height = h, useDingbats = FALSE)
  layout(matrix(c(1,2), nrow = 1), widths = c(5.0, 1.2))
  par(mar = c(5.5, 4.5, 2.5, 0.5), xaxs = 'i', yaxs = 'i')
  barplot(t(post_o), col = pal, border = NA, space = 0, axes = FALSE)
  axis(2, las = 1)
  axis(1, at = centers, labels = labels, cex.axis = 0.85)
  abline(v = seps, col = 'white', lwd = 1)
  box()
  mtext(sprintf('DAPC compoplot (best K=%d)', best_k), side = 3, line = 1, font = 2)
  mtext('Membership probability', side = 2, line = 3)
  mtext('Population', side = 1, line = 4.2)

  par(mar = c(5.5, 0.2, 2.5, 1.0))
  plot.new()
  legend('center', legend = colnames(post_o), fill = pal, bty = 'n', cex = 0.9)
  dev.off()
}

plot_comp(file.path(out_fig, 'DAPC_compoplot_bestK_publication.pdf'), 12, 4.8, png_dev = FALSE)
plot_comp(file.path(out_fig, 'DAPC_compoplot_bestK_publication.png'), 12, 4.8, png_dev = TRUE)

# BIC diagnostic figure
png(file.path(out_fig, 'DAPC_compoplot_K_BIC_diagnostic.png'), width = 1800, height = 1200, res = 300)
par(mar = c(4.5, 4.5, 2.2, 1))
plot(K_grid, bic_vals, type = 'b', pch = 19, xlab = 'K', ylab = 'BIC', main = 'DAPC-style K selection (k-means on retained PCs)')
abline(v = best_k, col = '#CC0000', lty = 2, lwd = 1.5)
dev.off()

pdf(file.path(out_fig, 'DAPC_compoplot_K_BIC_diagnostic.pdf'), width = 6, height = 4)
par(mar = c(4.5, 4.5, 2.2, 1))
plot(K_grid, bic_vals, type = 'b', pch = 19, xlab = 'K', ylab = 'BIC', main = 'DAPC-style K selection (k-means on retained PCs)')
abline(v = best_k, col = '#CC0000', lty = 2, lwd = 1.5)
dev.off()

png(file.path(out_fig, 'DAPC_compoplot_K_silhouette_diagnostic.png'), width = 1800, height = 1200, res = 300)
par(mar = c(4.5, 4.5, 2.2, 1))
plot(K_grid, sil_vals, type = 'b', pch = 19, xlab = 'K', ylab = 'Average silhouette', main = 'DAPC-style K selection (silhouette)')
abline(v = best_k, col = '#CC0000', lty = 2, lwd = 1.5)
dev.off()

pdf(file.path(out_fig, 'DAPC_compoplot_K_silhouette_diagnostic.pdf'), width = 6, height = 4)
par(mar = c(4.5, 4.5, 2.2, 1))
plot(K_grid, sil_vals, type = 'b', pch = 19, xlab = 'K', ylab = 'Average silhouette', main = 'DAPC-style K selection (silhouette)')
abline(v = best_k, col = '#CC0000', lty = 2, lwd = 1.5)
dev.off()

cat('Done\n')
