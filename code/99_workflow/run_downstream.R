#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
script_dir <- dirname(script_path)
default_out_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
default_stacks_dir <- normalizePath(file.path(default_out_root, ".."), mustWork = FALSE)

out_root <- Sys.getenv("DOWNSTREAM_ROOT", unset = default_out_root)
fig_dir <- file.path(out_root, "figures")
tab_dir <- file.path(out_root, "tables")
log_dir <- file.path(out_root, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stacks_dir <- Sys.getenv("STACKS_BASE_DIR", unset = default_stacks_dir)
vcf_file <- file.path(stacks_dir, "populations.snps.vcf")
pop_file <- file.path(stacks_dir, "popmap.tsv")
prune_in <- file.path(out_root, "work/admixture/trow_pruned.prune.in")
treemix_plotting_funcs <- Sys.getenv("TREEMIX_PLOTTING_FUNCS",
                                     unset = file.path(dirname(stacks_dir), ".envs", "treemix_x86", "bin", "plotting_funcs.R"))

pop_cols <- c(
  North = "#4169E1",        # royal blue
  North_Coast = "#4CBB17",  # Kelly green
  Sierra_1 = "#800080",     # purple
  Sierra_2 = "#00FFFF",     # cyan
  Sierra_3 = "#FF0000",     # red
  South_Coast = "#000000"   # black
)
pop_order <- names(pop_cols)

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

stopifnot_file <- function(f, label) {
  if (!file.exists(f)) stop(sprintf("Missing %s: %s", label, f))
}

stopifnot_file(vcf_file, "VCF input")
stopifnot_file(pop_file, "population map")

#-------------------------
# Load metadata
#-------------------------
logmsg("Loading popmap")
popmap <- read.table(pop_file, header = FALSE, sep = "\t", col.names = c("Sample", "Pop"))
popmap$Pop <- factor(popmap$Pop, levels = pop_order)

#-------------------------
# Parse VCF into genotype matrices
#-------------------------
logmsg("Reading VCF header")
header_lines <- readLines(vcf_file, n = 5000)
chrom_line <- header_lines[grepl("^#CHROM", header_lines)]
if (length(chrom_line) != 1) stop("Could not find #CHROM header in VCF.")
samples <- strsplit(chrom_line, "\t")[[1]][10:length(strsplit(chrom_line, "\t")[[1]])]

logmsg("Reading VCF body")
vcf <- read.table(vcf_file, sep = "\t", comment.char = "#", header = FALSE,
                  quote = "", fill = TRUE)
colnames(vcf)[1:9] <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
colnames(vcf)[10:ncol(vcf)] <- samples

# Reorder popmap to VCF sample order
popmap <- popmap[match(samples, popmap$Sample), ]
if (any(is.na(popmap$Sample))) stop("Some VCF samples are missing from popmap.")

# Genotype conversion helpers
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

logmsg("Converting all informative genotypes")
gt_mat <- as.matrix(vcf[, samples])
all_dos <- apply(gt_mat, 2, gt_to_nonref)   # variants x samples
rownames(all_dos) <- vcf$ID
colnames(all_dos) <- samples
X_all <- prep_matrix(t(all_dos))             # samples x variants

logmsg("Converting biallelic genotypes")
bial_idx <- !grepl(",", vcf$ALT, fixed = TRUE)
gt_bial <- gt_mat[bial_idx, , drop = FALSE]
ids_bial <- vcf$ID[bial_idx]
bial_dos <- apply(gt_bial, 2, gt_to_bial)
rownames(bial_dos) <- ids_bial
colnames(bial_dos) <- samples
X_bial <- prep_matrix(t(bial_dos))

write.table(data.frame(Sample = samples, Population = popmap$Pop),
            file.path(tab_dir, "sample_population_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

#-------------------------
# PCA analysis
#-------------------------
logmsg("Running PCA (all informative sites)")
pca_all <- prcomp(X_all, center = TRUE, scale. = TRUE)
var_all <- 100 * (pca_all$sdev^2 / sum(pca_all$sdev^2))
scores_all <- pca_all$x[, 1:4, drop = FALSE]

logmsg("Running PCA (biallelic sites)")
pca_bial <- prcomp(X_bial, center = TRUE, scale. = TRUE)
var_bial <- 100 * (pca_bial$sdev^2 / sum(pca_bial$sdev^2))
scores_bial <- pca_bial$x[, 1:4, drop = FALSE]

pca_var_tab <- data.frame(
  PC = paste0("PC", 1:10),
  all_informative_percent = round(var_all[1:10], 4),
  biallelic_percent = round(var_bial[1:10], 4)
)
write.table(pca_var_tab, file.path(tab_dir, "pca_variance_explained.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

plot_pca_matrix <- function(file, width, height) {
  if (grepl("\\.pdf$", file)) pdf(file, width = width, height = height, useDingbats = FALSE)
  else png(file, width = width * 300, height = height * 300, res = 300)

  par(mfrow = c(4, 4), mar = c(2.7, 2.7, 1.0, 0.6), oma = c(1.0, 1.0, 3.2, 0.6))

  for (i in 1:4) {
    for (j in 1:4) {
      if (i == j) {
        plot.new()
        text(0.5, 0.62, paste0("PC", i), cex = 1.2, font = 2)
        text(0.5, 0.45, paste0("Bial: ", sprintf("%.2f", var_bial[i]), "%"), cex = 0.9)
        text(0.5, 0.30, paste0("All:  ", sprintf("%.2f", var_all[i]), "%"), cex = 0.9)
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
             col = pop_cols[as.character(popmap$Pop)], pch = 16, cex = 0.75,
             xlab = paste0("PC", j, " (", sprintf("%.2f", vv[j]), "%)"),
             ylab = paste0("PC", i, " (", sprintf("%.2f", vv[i]), "%)"),
             main = subttl)
      }
    }
  }

  mtext("PCA Matrix: Upper Triangle = Biallelic Sites, Lower Triangle = All Informative Sites",
        side = 3, outer = TRUE, line = 0.5, cex = 1.1, font = 2)
  par(xpd = NA)
  legend("top", inset = c(0, -0.02), horiz = TRUE,
         legend = pop_order, col = pop_cols[pop_order], pch = 16, bty = "n", cex = 0.8)
  dev.off()
}

logmsg("Plotting PCA matrix figure")
plot_pca_matrix(file.path(fig_dir, "PCA_matrix_first4PCs_publication.pdf"), 14, 14)
plot_pca_matrix(file.path(fig_dir, "PCA_matrix_first4PCs_publication.png"), 14, 14)

#-------------------------
# Admixture-like analysis via cross-validated NMF
#-------------------------
logmsg("Preparing matrix for admixture-like NMF")
X_admix <- t(bial_dos)  # samples x biallelic variants (with NAs)
canon_id <- function(x) {
  y <- gsub(":", "_", x, fixed = TRUE)
  y <- sub("_[+-]$", "", y)
  y
}
col_can <- canon_id(colnames(X_admix))
if (file.exists(prune_in)) {
  prune_ids <- scan(prune_in, what = "character", quiet = TRUE)
} else {
  logmsg("LD-pruned loci list not found; using all biallelic loci for admixture matrix.")
  prune_ids <- col_can
}
keep <- col_can %in% prune_ids
X_admix <- X_admix[, keep, drop = FALSE]
colnames(X_admix) <- col_can[keep]
if (ncol(X_admix) == 0) stop("No overlap between pruned loci and biallelic loci for admixture matrix.")
X_admix <- prep_matrix(X_admix) / 2  # scale to [0,1]

nmf_fit_weighted <- function(X, K, W = NULL, maxiter = 400, tol = 1e-5, seed = 1) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  if (is.null(W)) W <- matrix(1, n, p)

  Q <- matrix(runif(n * K, 0.1, 1), n, K)
  Q <- Q / rowSums(Q)
  F <- matrix(runif(K * p, 0.1, 1), K, p)

  eps <- 1e-9
  prev <- Inf
  for (it in 1:maxiter) {
    XF <- W * X
    WH <- W * (Q %*% F)

    F <- F * ((t(Q) %*% XF) / (t(Q) %*% WH + eps))
    WH <- W * (Q %*% F)
    Q <- Q * ((XF %*% t(F)) / (WH %*% t(F) + eps))

    F[!is.finite(F)] <- eps
    F[F < eps] <- eps
    Q[!is.finite(Q)] <- eps
    Q[Q < eps] <- eps

    rs <- rowSums(Q)
    rs[!is.finite(rs) | rs < eps] <- 1
    Q <- Q / rs

    err <- sum((W * (X - Q %*% F))^2) / sum(W)
    if (!is.finite(err)) {
      err <- prev
      break
    }
    if (abs(prev - err) < tol) break
    prev <- err
  }
  list(Q = Q, F = F, error = prev)
}

cv_nmf <- function(X, K, folds = 5, reps = 3, seed = 1) {
  set.seed(seed + K)
  errs <- numeric(reps)
  n <- nrow(X); p <- ncol(X)
  for (r in 1:reps) {
    mask <- matrix(runif(n * p) > (1/folds), n, p)
    Xmask <- X
    Xmask[!mask] <- 0
    fit <- nmf_fit_weighted(Xmask, K, W = mask, maxiter = 300, seed = seed + K * 100 + r)
    pred <- fit$Q %*% fit$F
    errs[r] <- mean((X[!mask] - pred[!mask])^2)
  }
  c(mean = mean(errs), sd = sd(errs))
}

logmsg("Running cross-validated NMF for K=1..8")
Kvals <- 1:8
cv_tab <- data.frame(K = Kvals, cv_mse_mean = NA_real_, cv_mse_sd = NA_real_)
for (k in Kvals) {
  cv <- cv_nmf(X_admix, k, folds = 5, reps = 3, seed = 42)
  cv_tab$cv_mse_mean[cv_tab$K == k] <- cv[1]
  cv_tab$cv_mse_sd[cv_tab$K == k] <- cv[2]
  logmsg("K=", k, " CV-MSE=", sprintf("%.6f", cv[1]))
}
write.table(cv_tab, file.path(tab_dir, "admixture_nmf_cv_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

bestK <- cv_tab$K[which.min(cv_tab$cv_mse_mean)]
logmsg("Best K from CV NMF: ", bestK)

# Final NMF fit at best K (multi-start)
fits <- lapply(1:5, function(s) nmf_fit_weighted(X_admix, bestK, maxiter = 600, seed = 100 + s))
errs <- sapply(fits, `[[`, "error")
best_fit <- fits[[which.min(errs)]]
Qbest <- best_fit$Q
colnames(Qbest) <- paste0("Cluster", 1:ncol(Qbest))
rownames(Qbest) <- rownames(X_admix)

q_out <- data.frame(Sample = rownames(Qbest), Population = popmap$Pop[match(rownames(Qbest), popmap$Sample)], Qbest)
write.table(q_out, file.path(tab_dir, "admixture_nmf_Qmatrix_bestK.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot CV + compoplot
plot_admix <- function(file, width, height) {
  if (grepl("\\.pdf$", file)) pdf(file, width = width, height = height, useDingbats = FALSE)
  else png(file, width = width * 300, height = height * 300, res = 300)

  par(mfrow = c(2,1), mar = c(4,4,2,1))

  # CV curve
  plot(cv_tab$K, cv_tab$cv_mse_mean, type = "b", pch = 16, lwd = 2,
       xlab = "K (number of ancestral clusters)", ylab = "Cross-validated MSE",
       main = "Admixture-like NMF model selection")
  arrows(cv_tab$K, cv_tab$cv_mse_mean - cv_tab$cv_mse_sd,
         cv_tab$K, cv_tab$cv_mse_mean + cv_tab$cv_mse_sd,
         angle = 90, code = 3, length = 0.05)
  abline(v = bestK, col = "red", lty = 2, lwd = 2)
  mtext(paste("Best K =", bestK), side = 3, line = -1.2, adj = 0.95, col = "red", cex = 0.9)

  # Compoplot
  ord <- order(popmap$Pop[match(rownames(Qbest), popmap$Sample)], rownames(Qbest))
  Qord <- t(Qbest[ord, , drop = FALSE])
  pops_ord <- popmap$Pop[match(rownames(Qbest)[ord], popmap$Sample)]
  clust_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")[1:nrow(Qord)]

  bp <- barplot(Qord, col = clust_cols, border = NA, space = 0,
                xlab = "Samples (ordered by population)", ylab = "Ancestry proportion",
                main = paste0("Compoplot (Best K=", bestK, ")"))

  # population separators and labels
  run_lengths <- rle(as.character(pops_ord))
  ends <- cumsum(run_lengths$lengths)
  starts <- c(1, head(ends, -1) + 1)
  for (e in ends[-length(ends)]) abline(v = e + 0.5, lty = 3)
  mids <- (starts + ends) / 2
  text(mids, -0.07, labels = run_lengths$values, xpd = NA, srt = 45, adj = 1, cex = 0.7)

  legend("topright", legend = colnames(Qbest), fill = clust_cols, bty = "n", cex = 0.8)
  dev.off()
}

logmsg("Plotting admixture figure")
plot_admix(file.path(fig_dir, "Admixture_bestK_compoplot_publication.pdf"), 14, 10)
plot_admix(file.path(fig_dir, "Admixture_bestK_compoplot_publication.png"), 14, 10)

#-------------------------
# DAPC-like analysis (PCA + Fisher discriminant)
#-------------------------
logmsg("Running DAPC-like analysis (PCA + discriminant)")
X_dapc <- X_bial
y <- as.character(popmap$Pop[match(rownames(X_dapc), popmap$Sample)])

lda_fit <- function(X, y) {
  cls <- unique(y)
  p <- ncol(X)
  mu <- colMeans(X)
  Sw <- matrix(0, p, p)
  Sb <- matrix(0, p, p)
  for (cl in cls) {
    Xi <- X[y == cl, , drop = FALSE]
    mui <- colMeans(Xi)
    if (nrow(Xi) > 1) Sw <- Sw + (nrow(Xi) - 1) * cov(Xi)
    d <- matrix(mui - mu, ncol = 1)
    Sb <- Sb + nrow(Xi) * (d %*% t(d))
  }
  eig <- eigen(solve(Sw + diag(1e-6, p), Sb), symmetric = FALSE)
  list(W = Re(eig$vectors), values = Re(eig$values), classes = cls)
}

predict_lda <- function(train_scores, train_y, test_scores, W, n_ld) {
  A <- W[, 1:n_ld, drop = FALSE]
  tr <- train_scores %*% A
  te <- test_scores %*% A
  cent <- sapply(unique(train_y), function(cl) colMeans(tr[train_y == cl, , drop = FALSE]))
  if (is.null(dim(cent))) cent <- matrix(cent, ncol = 1)
  cl_names <- colnames(cent)
  pred <- character(nrow(te))
  for (i in 1:nrow(te)) {
    d <- colSums((t(cent) - matrix(te[i,], nrow = ncol(cent), ncol = ncol(te), byrow = TRUE))^2)
    pred[i] <- cl_names[which.min(d)]
  }
  pred
}

# choose number of PCs by 5-fold CV
set.seed(1)
n <- nrow(X_dapc)
fold <- sample(rep(1:5, length.out = n))
pc_grid <- seq(5, min(40, n - 2), by = 5)
acc <- numeric(length(pc_grid))

for (ii in seq_along(pc_grid)) {
  npc <- pc_grid[ii]
  preds <- character(n)
  for (f in 1:5) {
    tr_idx <- which(fold != f)
    te_idx <- which(fold == f)
    pca <- prcomp(X_dapc[tr_idx, , drop = FALSE], center = TRUE, scale. = TRUE)
    tr_scores <- pca$x[, 1:npc, drop = FALSE]
    te_scores <- scale(X_dapc[te_idx, , drop = FALSE], center = pca$center, scale = pca$scale) %*% pca$rotation[, 1:npc, drop = FALSE]
    fit <- lda_fit(tr_scores, y[tr_idx])
    n_ld <- min(length(unique(y)) - 1, ncol(tr_scores))
    preds[te_idx] <- predict_lda(tr_scores, y[tr_idx], te_scores, fit$W, n_ld)
  }
  acc[ii] <- mean(preds == y)
  logmsg("DAPC CV npc=", npc, " accuracy=", sprintf("%.4f", acc[ii]))
}

best_npc <- pc_grid[which.max(acc)]
logmsg("Selected PCs for DAPC-like model: ", best_npc)

# fit final model
pca_final <- prcomp(X_dapc, center = TRUE, scale. = TRUE)
pcs <- pca_final$x[, 1:best_npc, drop = FALSE]
lda_final <- lda_fit(pcs, y)
ndim <- min(length(unique(y)) - 1, ncol(pcs), 2)
ld_scores <- pcs %*% lda_final$W[, 1:ndim, drop = FALSE]

dapc_tab <- data.frame(Sample = rownames(X_dapc), Population = y,
                       LD1 = ld_scores[,1],
                       LD2 = if (ndim >= 2) ld_scores[,2] else 0)
write.table(dapc_tab, file.path(tab_dir, "dapc_like_scores.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(n_pcs = pc_grid, cv_accuracy = acc),
            file.path(tab_dir, "dapc_like_cv_accuracy.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

plot_dapc <- function(file, width, height) {
  if (grepl("\\.pdf$", file)) pdf(file, width = width, height = height, useDingbats = FALSE)
  else png(file, width = width * 300, height = height * 300, res = 300)

  par(mfrow = c(1,2), mar = c(4,4,2,1))
  plot(pc_grid, acc, type = "b", pch = 16, lwd = 2,
       xlab = "Number of retained PCs", ylab = "5-fold CV accuracy",
       main = "DAPC-like model tuning")
  abline(v = best_npc, col = "red", lty = 2)

  plot(dapc_tab$LD1, dapc_tab$LD2,
       col = pop_cols[dapc_tab$Population], pch = 16, cex = 0.9,
       xlab = "Discriminant axis 1", ylab = "Discriminant axis 2",
       main = paste0("DAPC-like scatter (", best_npc, " PCs retained)"))
  for (pp in pop_order) {
    idx <- which(dapc_tab$Population == pp)
    if (length(idx) >= 3) {
      h <- chull(dapc_tab$LD1[idx], dapc_tab$LD2[idx])
      polygon(dapc_tab$LD1[idx][h], dapc_tab$LD2[idx][h], border = pop_cols[pp], lwd = 1.2)
    }
  }
  legend("topright", legend = pop_order, col = pop_cols[pop_order], pch = 16, bty = "n", cex = 0.8)
  dev.off()
}

logmsg("Plotting DAPC-like figure")
plot_dapc(file.path(fig_dir, "DAPC_like_publication.pdf"), 14, 6)
plot_dapc(file.path(fig_dir, "DAPC_like_publication.png"), 14, 6)

#-------------------------
# TreeMix figures
#-------------------------
logmsg("Preparing TreeMix figures")
llik_files <- file.path(out_root, "work/treemix", paste0("tm_m", 0:8, ".llik"))
missing_llik <- llik_files[!file.exists(llik_files)]
best_m <- NA_integer_
if (length(missing_llik)) {
  logmsg("TreeMix files not present in work/treemix; skipping TreeMix plotting in run_downstream.")
  logmsg("Run scripts/run_treemix_models_2026-02-27.R then scripts/run_treemix_figs.R to generate TreeMix outputs.")
} else {
  get_llik <- function(f) {
    ln <- readLines(f)
    x <- sub(".*Exiting ln\\(likelihood\\) with [0-9]+ migration events: ([0-9.-]+).*", "\\1", ln[2])
    as.numeric(x)
  }
  llik <- sapply(llik_files, get_llik)
  km <- 0:8
  max_ll <- max(llik, na.rm = TRUE)
  best_m <- min(km[which(llik >= max_ll - 1e-6)])

  write.table(data.frame(migration_edges = km, log_likelihood = llik),
              file.path(tab_dir, "treemix_loglikelihood_by_m.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Use treemix plotting helpers
  stopifnot_file(treemix_plotting_funcs, "TreeMix plotting helper")
  source(treemix_plotting_funcs)
  label_cols <- cbind(pop_order, pop_cols[pop_order])
  pop_order_vec <- pop_order

  plot_treemix <- function(file, width, height) {
    if (grepl("\\.pdf$", file)) pdf(file, width = width, height = height, useDingbats = FALSE)
    else png(file, width = width * 300, height = height * 300, res = 300)

    par(mfrow = c(1,2), mar = c(4,4,2,1))
    plot(km, llik, type = "b", pch = 16, lwd = 2,
         xlab = "Migration edges (m)", ylab = "TreeMix log-likelihood",
         main = "TreeMix model fit")
    abline(v = best_m, col = "red", lty = 2, lwd = 2)
    mtext(paste("Chosen m =", best_m), side = 3, line = -1.2, adj = 0.95, col = "red")

    stem <- file.path(out_root, "work/treemix", paste0("tm_m", best_m))
    plot_tree(stem, o = NA, cex = 0.9, disp = 0.01, plus = 0.02,
              mbar = TRUE, scale = TRUE, plotmig = TRUE, lwd = 1.4)
    title(main = paste0("TreeMix graph (m=", best_m, ")"))

    dev.off()
  }

  logmsg("Plotting TreeMix figure")
  plot_treemix(file.path(fig_dir, "TreeMix_publication.pdf"), 14, 6)
  plot_treemix(file.path(fig_dir, "TreeMix_publication.png"), 14, 6)

  # residual plot (extra)
  stem <- file.path(out_root, "work/treemix", paste0("tm_m", best_m))
  for (sx in c(".cov.gz", ".treeout.gz", ".llik")) {
    stopifnot_file(paste0(stem, sx), paste0("TreeMix residual input ", sx))
  }

  pdf(file.path(fig_dir, "TreeMix_residuals_publication.pdf"), width = 7, height = 6, useDingbats = FALSE)
  plot_resid(stem, pop_order_vec, cex = 0.8)
  title(main = paste0("TreeMix residuals (m=", best_m, ")"))
  dev.off()

  png(file.path(fig_dir, "TreeMix_residuals_publication.png"), width = 2100, height = 1800, res = 300)
  plot_resid(stem, pop_order_vec, cex = 0.8)
  title(main = paste0("TreeMix residuals (m=", best_m, ")"))
  dev.off()
}

# summary report
summary_txt <- file.path(tab_dir, "analysis_summary.txt")
cat(
  "Stacks downstream analysis summary\n",
  "================================\n",
  "Samples:", length(samples), "\n",
  "All informative loci used in PCA:", ncol(X_all), "\n",
  "Biallelic loci used in PCA/DAPC:", ncol(X_bial), "\n",
  "NMF-admixture loci (LD-pruned):", ncol(X_admix), "\n",
  "Best K (NMF CV):", bestK, "\n",
  "Best TreeMix m:", ifelse(is.na(best_m), "not_run", best_m), "\n",
  sep = "",
  file = summary_txt
)

logmsg("All downstream analyses completed.")
