#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
out_dir <- file.path(down_dir, "popcluster_native_2026-02-27")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
work_copy <- file.path(out_dir, "work")
for (d in c(out_dir, fig_dir, tab_dir, log_dir, work_copy)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_dir <- "/tmp/popcluster_native_excl_MVZ216210_2026-02-27"
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(base_dir, "populations.snps.vcf")
pop_file <- file.path(base_dir, "popmap.tsv")
popcluster_bin <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/PopCluster_2025/Bin/PopClusterMac"

if (!file.exists(vcf_file)) stop("Missing VCF: ", vcf_file)
if (!file.exists(pop_file)) stop("Missing popmap: ", pop_file)
if (!file.exists(popcluster_bin)) stop("Missing PopCluster binary: ", popcluster_bin)

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
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(line, "\n", sep = "")
  write(line, file = file.path(log_dir, "run_popcluster_native.log"), append = TRUE)
  flush.console()
}

parse_popmap <- function(f) {
  p <- read.table(f, header = FALSE, sep = "\t", fill = TRUE, comment.char = "")
  if (ncol(p) < 2) stop("Invalid popmap")
  p <- p[, 1:2, drop = FALSE]
  colnames(p) <- c("Sample", "Population")
  if (tolower(p$Sample[1]) %in% c("sample", "samples", "id") && tolower(p$Population[1]) %in% c("population", "pop", "group")) {
    p <- p[-1, , drop = FALSE]
  }
  p$Population <- factor(p$Population, levels = pop_order)
  p
}

as_gt <- function(x) {
  x <- sub(":.*", "", x)
  x <- gsub("\\|", "/", x)
  x[x %in% c("./.", ".|.", ".")] <- NA
  x
}

gt_to_bial_dos <- function(x) {
  x <- as_gt(x)
  a1 <- sub("/.*", "", x)
  a2 <- sub(".*/", "", x)
  good <- !(is.na(x) | a1 == "." | a2 == ".")
  valid <- good & a1 %in% c("0", "1") & a2 %in% c("0", "1")
  out <- rep(NA_integer_, length(x))
  out[valid] <- as.integer(a1[valid] == "1") + as.integer(a2[valid] == "1")
  out
}

logmsg("Reading VCF and building PopCluster input matrix")
head_lines <- readLines(vcf_file, n = 5000)
chrom_line <- head_lines[grepl("^#CHROM", head_lines)]
if (length(chrom_line) != 1) stop("No #CHROM line in VCF")
samples <- strsplit(chrom_line, "\t")[[1]][10:length(strsplit(chrom_line, "\t")[[1]])]

pop <- parse_popmap(pop_file)
pop <- pop[match(samples, pop$Sample), , drop = FALSE]
if (any(is.na(pop$Sample))) stop("Some VCF samples missing from popmap")

vcf <- read.table(vcf_file, sep = "\t", comment.char = "#", header = FALSE, quote = "", fill = TRUE)
colnames(vcf)[1:9] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
colnames(vcf)[10:ncol(vcf)] <- samples

is_snp <- nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1
is_bial <- !grepl(",", vcf$ALT, fixed = TRUE)
keep <- is_snp & is_bial
vcf <- vcf[keep, , drop = FALSE]

mat <- as.matrix(vcf[, samples, drop = FALSE])
dos <- apply(mat, 2, gt_to_bial_dos) # variants x samples

# Standard filters for ddRAD SNPs
callrate <- rowMeans(!is.na(dos))
af <- rowMeans(dos, na.rm = TRUE) / 2
maf <- pmin(af, 1 - af)
keep2 <- is.finite(maf) & callrate >= 0.80 & maf >= 0.01
dos <- dos[keep2, , drop = FALSE]

# Remove invariant after filtering
v <- apply(dos, 1, function(x) var(x, na.rm = TRUE))
dos <- dos[is.finite(v) & v > 0, , drop = FALSE]

n_ind <- ncol(dos)
n_loci <- nrow(dos)
if (n_ind < 10 || n_loci < 100) stop("Insufficient data after filters: inds=", n_ind, " loci=", n_loci)

logmsg("PopCluster native input dimensions: individuals=", n_ind, " loci=", n_loci)

# Write .dat in PopCluster 1-row format: ind_index a1_l1 a2_l1 a1_l2 a2_l2 ...
# Alleles encoded as 1/2 for ref/alt, missing as 0/0
make_alleles <- function(g) {
  out <- matrix(0L, nrow = length(g), ncol = 2)
  idx0 <- which(g == 0L)
  idx1 <- which(g == 1L)
  idx2 <- which(g == 2L)
  if (length(idx0)) out[idx0, ] <- cbind(1L, 1L)
  if (length(idx1)) out[idx1, ] <- cbind(1L, 2L)
  if (length(idx2)) out[idx2, ] <- cbind(2L, 2L)
  out
}

dat_file <- file.path(run_dir, "trow_exclMVZ216210.dat")
con <- file(dat_file, open = "wt")
on.exit(close(con), add = TRUE)

for (i in seq_len(n_ind)) {
  g <- dos[, i]
  aa <- make_alleles(g)
  vec <- c(i, as.vector(t(aa)))
  writeLines(paste(vec, collapse = " "), con = con)
}
close(con)

# Write .PcPjt parameter file
k_min <- 1
k_max <- 10
n_rep <- 5
seed <- 20260227

pjt_file <- file.path(run_dir, "trow_exclMVZ216210.PcPjt")
out_base <- "trow_exclMVZ216210"
pjt_lines <- c(
  sprintf("%d                                                   !Integer, #Individuals", n_ind),
  sprintf("%d                                                     !Integer, #Loci", n_loci),
  "1                                                     !Boolean, All loci SNP (1/0=Y/N)",
  "0                                                     !String, Missing allele",
  sprintf("%d                                                   !Integer, Random number seed", seed),
  "trow_exclMVZ216210.dat                                !String, Genotypefilename",
  "trow_exclMVZ216210                                    !String, Outputfilename",
  "0                                                     !integer, 3/2/1/0 = strong/medium/weak/no scaling",
  sprintf("%d                                                    !Integer, Minimum K", k_min),
  sprintf("%d                                                    !Integer, Maximum K", k_max),
  sprintf("%d                                                     !Integer, Num replicate runs per K", n_rep),
  "0                                                     !Integer, 0/1=Search using assignment_prob/relatedness",
  "0                                                     !Boolean, 1/0=Output allele frequency:YES/NO",
  "0                                                     !Boolean, 1/0=PopData available:YES/NO",
  "0                                                     !Boolean, 1/0=PopFlag available:YES/NO",
  "2                                                     !Integer, 1/2/3/4=Mixture/Admixture/Hybridyzation/Migration model",
  "1                                                     !Boolean, 1/0=Estimate locus-specific F-Statistics=Y/N",
  "1                                                     !Boolean, 1/0=Use K-Means clustering method=Y/N",
  "0                                                     !Boolean, 1/0=Individual location data available=Y/N",
  "0                                                     !Integer, 0/1/2=Individual data in 1-row/2-rows/1-column",
  "1                                                     !Boolean, 1/0=Infer locus admixture=Y/N",
  "0                                                     !Boolean, 0/1/2: Compute relatedness = No/Wang/LynchRitland",
  "0                                                     !Boolean, 0/1 estimate kinship = No/Yes",
  "2                                                     !Integer, 0/1/2=Undefined/equal/unequal prior allele freq"
)
writeLines(pjt_lines, con = pjt_file)

logmsg("Running native PopCluster K=", k_min, "..", k_max, " reps=", n_rep)
cmd_args <- c(paste0("INP:", basename(pjt_file)), "LON:0", "OMP:4", "KIN:0", "COP:0")
k_file <- file.path(run_dir, paste0(out_base, ".K"))
if (!file.exists(k_file)) {
  old_wd <- getwd()
  setwd(run_dir)
  on.exit(setwd(old_wd), add = TRUE)
  run_log <- system2(popcluster_bin, args = cmd_args, stdout = TRUE, stderr = TRUE)
  setwd(old_wd)
  writeLines(run_log, con = file.path(log_dir, "popcluster_native_stdout.log"))
} else {
  logmsg("Existing PopCluster outputs detected in run dir; skipping rerun and reusing results.")
}

if (!file.exists(k_file)) stop("PopCluster did not produce K summary file: ", k_file)

# Parse .K summary table
k_lines <- readLines(k_file, warn = FALSE)
start <- grep("^\\s*K\\s+BestRun", k_lines)
if (!length(start)) stop("Cannot parse K summary header in ", k_file)

# data block ends before Method section
meth <- grep("^\\s*Method\\s+Best_K", k_lines)
end <- if (length(meth)) meth[1] - 1 else length(k_lines)
blk <- k_lines[(start[1] + 1):end]
blk <- blk[nzchar(trimws(blk))]
k_raw <- read.table(text = blk, header = FALSE, quote = "\"", fill = TRUE, stringsAsFactors = FALSE)
if (ncol(k_raw) < 8) stop("Failed to parse K rows from ", k_file)
k_raw <- k_raw[, 1:8, drop = FALSE]
colnames(k_raw) <- c("K", "BestRun", "LogL_Mean", "LogL_Min", "LogL_Max", "DLK1", "DLK2", "FST_FIS")
num_or_na <- function(x) {
  x <- as.character(x)
  x[x == "-"] <- NA_character_
  suppressWarnings(as.numeric(x))
}
k_tab <- data.frame(
  K = as.integer(k_raw$K),
  BestRun = as.character(k_raw$BestRun),
  LogL_Mean = num_or_na(k_raw$LogL_Mean),
  LogL_Min = num_or_na(k_raw$LogL_Min),
  LogL_Max = num_or_na(k_raw$LogL_Max),
  DLK1 = num_or_na(k_raw$DLK1),
  DLK2 = num_or_na(k_raw$DLK2),
  FST_FIS = num_or_na(k_raw$FST_FIS),
  stringsAsFactors = FALSE
)
if (!nrow(k_tab)) stop("Failed to parse K rows from ", k_file)

best_k_dlk2 <- if (any(is.finite(k_tab$DLK2))) k_tab$K[which.max(k_tab$DLK2)] else NA_integer_
best_k_fst <- if (any(is.finite(k_tab$FST_FIS))) k_tab$K[which.max(k_tab$FST_FIS)] else NA_integer_

best_k <- if (!is.na(best_k_dlk2)) best_k_dlk2 else best_k_fst
if (is.na(best_k)) stop("Failed to infer best K from DLK2 or FST/FIS")

best_run <- k_tab$BestRun[k_tab$K == best_k][1]
q_file <- file.path(run_dir, paste0(best_run, "_Q"))
if (!file.exists(q_file)) {
  q_alt <- gsub("_K_", ".K.", best_run, fixed = TRUE)
  q_alt <- gsub("_R_", ".R.", q_alt, fixed = TRUE)
  q_file <- file.path(run_dir, paste0(q_alt, "_Q"))
}
if (!file.exists(q_file)) stop("Missing Q output for best run: ", q_file)

q <- read.csv(q_file, header = FALSE)
if (nrow(q) != n_ind) stop("Q matrix row count mismatch")
colnames(q) <- c("Index", paste0("Cluster", seq_len(ncol(q) - 1)))
q$Sample <- samples[match(q$Index, seq_len(n_ind))]
q$Population <- as.character(pop$Population[match(q$Sample, pop$Sample)])

# Save tables
write.table(k_tab, file.path(tab_dir, "popcluster_native_k_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(best_k = best_k, best_k_DLK2 = best_k_dlk2, best_k_FSTFIS = best_k_fst, best_run = best_run),
            file.path(tab_dir, "popcluster_native_best_k.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(q[, c("Sample", "Population", grep("^Cluster", colnames(q), value = TRUE))],
            file.path(tab_dir, "popcluster_native_Qmatrix_bestK.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Plot publication-quality figure (K diagnostics + compoplot)
ord <- order(factor(q$Population, levels = pop_order), q$Sample)
q_ord <- q[ord, , drop = FALSE]
cluster_cols <- grep("^Cluster", colnames(q_ord), value = TRUE)
Qmat <- t(as.matrix(q_ord[, cluster_cols, drop = FALSE]))

cluster_palette <- grDevices::hcl.colors(nrow(Qmat), palette = "Dark 3")

plot_figure <- function(file, width = 14, height = 9) {
  if (grepl("\\.pdf$", file)) pdf(file, width = width, height = height, useDingbats = FALSE)
  else png(file, width = width * 320, height = height * 320, res = 320)

  layout(matrix(c(1,2), nrow = 2), heights = c(1, 1.7))

  par(mar = c(4,4,2,1))
  plot(k_tab$K, k_tab$LogL_Mean, type = "b", pch = 16, lwd = 2,
       xlab = "K", ylab = "Mean log-likelihood", main = "PopCluster model fit")
  if (any(is.finite(k_tab$DLK2))) {
    par(new = TRUE)
    plot(k_tab$K, k_tab$DLK2, type = "b", pch = 17, lwd = 1.8, axes = FALSE, xlab = "", ylab = "", col = "#D55E00")
    axis(4, col.axis = "#D55E00")
    mtext("DLK2", side = 4, line = 2.5, col = "#D55E00")
    legend("bottomright", legend = c("LogL mean", "DLK2"), col = c("black", "#D55E00"), pch = c(16,17), lwd = c(2,1.8), bty = "n")
  }
  abline(v = best_k, col = "red", lty = 2, lwd = 2)
  mtext(sprintf("Best K = %d", best_k), side = 3, line = -1.2, adj = 0.95, col = "red", cex = 0.9)

  par(mar = c(6,4,2,1), xpd = NA)
  barplot(Qmat, col = cluster_palette, border = NA, space = 0, axisnames = FALSE,
          ylab = "Ancestry proportion", main = paste0("PopCluster compoplot (K=", best_k, ")"))
  rl <- rle(q_ord$Population)
  ends <- cumsum(rl$lengths)
  starts <- c(1, head(ends, -1) + 1)
  mids <- (starts + ends) / 2
  for (e in ends[-length(ends)]) abline(v = e + 0.5, lty = 3, col = "grey40")
  axis(1, at = mids, labels = rl$values, las = 2, cex.axis = 0.85)
  mtext("Population", side = 1, line = 4.5)

  # no cluster legend per recent preference
  mtext("Native PopCluster analysis (MVZ216210 excluded)", side = 3, outer = TRUE, line = -1.0, font = 2, cex = 1.05)
  dev.off()
}

plot_figure(file.path(fig_dir, "PopCluster_native_bestK_compoplot_publication.pdf"))
plot_figure(file.path(fig_dir, "PopCluster_native_bestK_compoplot_publication.png"))

# Copy run products for reproducibility
run_files <- list.files(run_dir, full.names = TRUE)
file.copy(run_files, work_copy, overwrite = TRUE)

writeLines(c(
  "PopCluster native analysis summary",
  paste0("Binary: ", popcluster_bin),
  paste0("Input VCF: ", vcf_file),
  paste0("Individuals: ", n_ind),
  paste0("Loci retained: ", n_loci),
  paste0("K range: ", k_min, "-", k_max),
  paste0("Replicates per K: ", n_rep),
  paste0("Best K by DLK2: ", ifelse(is.na(best_k_dlk2), "NA", best_k_dlk2)),
  paste0("Best K by FST/FIS: ", ifelse(is.na(best_k_fst), "NA", best_k_fst)),
  paste0("Selected best K: ", best_k),
  paste0("Best run file stem: ", best_run),
  paste0("Figure PDF: ", file.path(fig_dir, "PopCluster_native_bestK_compoplot_publication.pdf")),
  paste0("Figure PNG: ", file.path(fig_dir, "PopCluster_native_bestK_compoplot_publication.png"))
), con = file.path(tab_dir, "popcluster_native_analysis_summary.txt"))

logmsg("Native PopCluster analysis complete")
