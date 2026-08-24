#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
native_dir <- file.path(down_dir, "popcluster_native_2026-02-27")
fig_dir <- file.path(native_dir, "figures")
tab_dir <- file.path(native_dir, "tables")
work_dir <- file.path(native_dir, "work")

pop_file <- file.path(base_dir, "popmap.tsv")
if (!file.exists(pop_file)) stop("Missing popmap")
if (!file.exists(file.path(tab_dir, "popcluster_native_k_summary.tsv"))) stop("Missing popcluster_native_k_summary.tsv")

pop_cols <- c(
  North = "#4169E1",
  North_Coast = "#4CBB17",
  Sierra_1 = "#800080",
  Sierra_2 = "#00FFFF",
  Sierra_3 = "#FF0000",
  South_Coast = "#000000"
)
pop_order <- names(pop_cols)

pop <- read.table(pop_file, header = FALSE, sep = "\t", fill = TRUE, comment.char = "")
pop <- pop[,1:2, drop = FALSE]
colnames(pop) <- c("Sample","Population")
if (tolower(pop$Sample[1]) %in% c("sample","samples","id") && tolower(pop$Population[1]) %in% c("population","pop","group")) pop <- pop[-1,,drop=FALSE]
pop$Population <- factor(pop$Population, levels = pop_order)

ks <- read.delim(file.path(tab_dir, "popcluster_native_k_summary.tsv"), sep = "\t", stringsAsFactors = FALSE)

kset <- c(6,7,8)
get_best_run <- function(K) {
  row <- ks[ks$K == K, , drop = FALSE]
  if (!nrow(row)) stop("No K row in summary for K=", K)
  row$BestRun[1]
}

load_q <- function(best_run) {
  q_name <- gsub("_K_", ".K.", best_run, fixed = TRUE)
  q_name <- gsub("_R_", ".R.", q_name, fixed = TRUE)
  q_file <- file.path(work_dir, paste0(q_name, "_Q"))
  if (!file.exists(q_file)) stop("Missing Q file: ", q_file)
  q <- read.csv(q_file, header = FALSE)
  colnames(q) <- c("Index", paste0("Cluster", seq_len(ncol(q)-1)))
  q$Index <- as.integer(q$Index)
  q$Sample <- pop$Sample[q$Index]
  q$Population <- as.character(pop$Population[q$Index])
  q
}

q_list <- lapply(kset, function(K) {
  br <- get_best_run(K)
  q <- load_q(br)
  q
})
names(q_list) <- paste0("K", kset)

# consistent ordering across all panels
ord <- order(pop$Population, pop$Sample)
ord_samples <- pop$Sample[ord]

plot_multik <- function(file, width = 14, height = 10) {
  if (grepl("\\.pdf$", file)) pdf(file, width = width, height = height, useDingbats = FALSE)
  else png(file, width = width * 320, height = height * 320, res = 320)

  par(mfrow = c(3,1), mar = c(3.0, 4.0, 1.9, 1.0), oma = c(6.0, 0.2, 1.2, 0.2), xpd = NA)

  for (i in seq_along(kset)) {
    K <- kset[i]
    q <- q_list[[paste0("K", K)]]
    q <- q[match(ord_samples, q$Sample), , drop = FALSE]

    qmat <- t(as.matrix(q[, grep("^Cluster", colnames(q)), drop = FALSE]))
    cols <- grDevices::hcl.colors(nrow(qmat), palette = "Dark 3")

    barplot(qmat, col = cols, border = NA, space = 0, axes = FALSE, axisnames = FALSE, yaxs = "i")
    axis(2, at = c(0,0.5,1), las = 1)
    if (i == 1) mtext("Ancestry proportion", side = 2, line = 2.7, cex = 0.95)

    rl <- rle(q$Population)
    ends <- cumsum(rl$lengths)
    for (e in ends[-length(ends)]) abline(v = e + 0.5, lty = 3, col = "grey45")
    starts <- c(1, head(ends, -1) + 1)
    mids <- (starts + ends) / 2
    if (i < length(kset)) {
      axis(1, at = numeric(0), labels = FALSE, tck = 0)
    } else {
      axis(1, at = mids, labels = rl$values, las = 2, cex.axis = 0.85)
      mtext("Population", side = 1, line = 4.6, cex = 0.95)
    }
    mtext(paste0("K = ", K), side = 3, line = 0.3, font = 2, cex = 1.0)
  }

  mtext("Native PopCluster compoplots", side = 3, outer = TRUE, line = 0.2, cex = 1.1, font = 2)
  dev.off()
}

plot_multik(file.path(fig_dir, "PopCluster_native_compoplot_K6_K7_K8_publication.pdf"))
plot_multik(file.path(fig_dir, "PopCluster_native_compoplot_K6_K7_K8_publication.png"))

writeLines(c(
  "Native PopCluster side-by-side compoplots",
  "Panels: K=6, K=7, K=8",
  "No cluster legend included",
  paste0("Figure PDF: ", file.path(fig_dir, "PopCluster_native_compoplot_K6_K7_K8_publication.pdf")),
  paste0("Figure PNG: ", file.path(fig_dir, "PopCluster_native_compoplot_K6_K7_K8_publication.png"))
), con = file.path(tab_dir, "popcluster_native_compoplot_K6_K7_K8_summary.txt"))

cat("Done\n")
