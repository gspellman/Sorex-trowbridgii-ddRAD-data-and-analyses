#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ape))

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
run_dir <- file.path(down_dir, "species_delimitation_2026-02-25_ess_bpp_reconfigured")
bpp_dir <- file.path(run_dir, "bpp_attempt")
tab_dir <- file.path(run_dir, "tables")
log_dir <- file.path(run_dir, "logs")
for (d in c(run_dir, bpp_dir, tab_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

vcf_file <- file.path(base_dir, "populations.snps.vcf")
pop_file <- file.path(base_dir, "popmap.tsv")
bpp_bin <- "/Users/gspellman/Downloads/bpp-4.8.7-macos-aarch64/bin/bpp"
if (!file.exists(bpp_bin)) stop("BPP binary not found: ", bpp_bin)
iqtree_contree <- file.path(down_dir, "iqtree_analysis_2026-03-01", "trow_individuals_biallelic_iqtree.contree")
if (!file.exists(iqtree_contree)) stop("IQ-TREE file not found: ", iqtree_contree)

samples_per_pop <- 2
max_snps <- 100
set.seed(20260226)

parse_gt <- function(gt) {
  gt <- sub(":.*", "", gt)
  gt <- gsub("\\|", "/", gt)
  if (gt %in% c("./.", ".|.", ".")) return(NA_real_)
  if (gt == "0/0") return(0)
  if (gt %in% c("0/1", "1/0")) return(1)
  if (gt == "1/1") return(2)
  NA_real_
}

calc_ess <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 50) return(NA_real_)
  ac <- stats::acf(x, lag.max = min(1000, floor(n / 2)), plot = FALSE, demean = TRUE)$acf
  rho <- as.numeric(ac[-1])
  pos <- rho[rho > 0]
  if (!length(pos)) return(n)
  tau <- 1 + 2 * sum(pos)
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  n / tau
}

calc_min_ess_from_mcmc <- function(mcmc_file) {
  if (!file.exists(mcmc_file)) return(NA_real_)
  d <- tryCatch(read.table(mcmc_file, header = TRUE, check.names = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) < 150) return(NA_real_)
  keep <- names(d)[vapply(d, is.numeric, logical(1))]
  keep <- setdiff(keep, c("Gen", "gen", "sample", "Sample"))
  if (!length(keep)) return(NA_real_)
  burn <- floor(0.2 * nrow(d))
  p <- d[(burn + 1):nrow(d), , drop = FALSE]
  ess <- vapply(keep, function(k) calc_ess(p[[k]]), numeric(1))
  suppressWarnings(min(ess, na.rm = TRUE))
}

pop <- read.table(pop_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(pop) <- c("Sample", "Population")
pop_levels <- unique(pop$Population)
selected <- unlist(lapply(pop_levels, function(pp) {
  ids <- pop$Sample[pop$Population == pp]
  sample(ids, min(samples_per_pop, length(ids)))
}), use.names = FALSE)
selected <- unique(selected)
sel_pop <- setNames(pop$Population[match(selected, pop$Sample)], selected)

# Build a population-level BPP guide tree from the IQ-TREE individual tree.
iq <- read.tree(iqtree_contree)
iq_keep <- intersect(selected, iq$tip.label)
if (length(iq_keep) < 6) stop("Too few selected samples overlap with IQ-TREE tips")
iq_sub <- drop.tip(iq, setdiff(iq$tip.label, iq_keep))
dm <- cophenetic.phylo(iq_sub)
spp_all <- sort(unique(sel_pop[selected]))
pop_dm <- matrix(0, nrow = length(spp_all), ncol = length(spp_all), dimnames = list(spp_all, spp_all))
for (i in seq_along(spp_all)) {
  for (j in seq_along(spp_all)) {
    if (i == j) next
    ai <- names(sel_pop)[sel_pop == spp_all[i]]
    bj <- names(sel_pop)[sel_pop == spp_all[j]]
    ai <- intersect(ai, rownames(dm))
    bj <- intersect(bj, colnames(dm))
    if (!length(ai) || !length(bj)) stop("Cannot compute inter-population IQ-tree distances for ", spp_all[i], " vs ", spp_all[j])
    pop_dm[i, j] <- mean(dm[ai, bj, drop = FALSE], na.rm = TRUE)
  }
}
diag(pop_dm) <- 0
guide_tree <- nj(as.dist(pop_dm))
guide_tree$tip.label <- spp_all
# BPP requires a rooted binary species tree.
outgroup_pop <- names(which.max(rowMeans(pop_dm, na.rm = TRUE)))[1]
guide_tree <- root(guide_tree, outgroup = outgroup_pop, resolve.root = TRUE)
if (!is.null(guide_tree$edge.length)) {
  guide_tree$edge.length[!is.finite(guide_tree$edge.length) | guide_tree$edge.length <= 0] <- 1e-06
}
guide_tree_newick <- write.tree(guide_tree)

con <- file(vcf_file, "r")
on.exit(close(con), add = TRUE)
hdr <- character()
repeat {
  ln <- readLines(con, n = 1)
  if (!length(ln)) stop("VCF header not found")
  if (startsWith(ln, "#CHROM")) {
    hdr <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    break
  }
}
vcf_samples <- hdr[10:length(hdr)]
sample_idx <- match(selected, vcf_samples)
if (any(is.na(sample_idx))) stop("Selected samples missing from VCF")

rows <- list()
ri <- 0L
repeat {
  ln <- readLines(con, n = 1)
  if (!length(ln)) break
  if (!nzchar(ln) || substr(ln, 1, 1) == "#") next
  p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  if (length(p) < 10 + max(sample_idx)) next
  ref <- p[4]
  alt <- p[5]
  if (nchar(ref) != 1 || nchar(alt) != 1 || grepl(",", alt, fixed = TRUE)) next
  g <- p[10:length(p)]
  vals <- vapply(sample_idx, function(j) parse_gt(g[j]), numeric(1))
  if (any(!is.finite(vals))) next
  if (length(unique(vals)) < 2) next
  ri <- ri + 1L
  rows[[ri]] <- vals
}

if (length(rows) < 80) stop("Too few complete SNPs for BPP run")
G <- do.call(rbind, rows)
colnames(G) <- selected
if (nrow(G) > max_snps) {
  keep <- sort(sample(seq_len(nrow(G)), max_snps, replace = FALSE))
  G <- G[keep, , drop = FALSE]
}

dosage_to_char <- function(v) ifelse(v == 0, "A", ifelse(v == 1, "M", "T"))
seq_char <- apply(G, 2, function(col) paste0(dosage_to_char(as.integer(col)), collapse = ""))
seq_char <- seq_char[selected]
if (any(is.na(seq_char) | !nzchar(seq_char))) stop("Sequence construction failed")

seq_file <- file.path(bpp_dir, "bpp_input_single_locus.phy")
imap_file <- file.path(bpp_dir, "bpp_imap.txt")
ctl_file <- file.path(bpp_dir, "bpp.ctl")
run_log <- file.path(log_dir, "bpp_attempt.log")
jobname <- "bpp_run"
mcmc_file <- file.path(bpp_dir, paste0(jobname, ".mcmc.txt"))

cat(sprintf("%d %d\n", length(selected), nrow(G)), file = seq_file)
for (tx in selected) cat(sprintf("%-18s %s\n", paste0("^", tx), seq_char[tx]), file = seq_file, append = TRUE)
write.table(data.frame(Sample = selected, Species = sel_pop[selected], stringsAsFactors = FALSE),
            imap_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write_ctl <- function(burnin, nsample, sampfreq) {
  spp <- unique(sel_pop[selected])
  spp <- spp[order(spp)]
  spp_counts <- vapply(spp, function(sp) sum(sel_pop[selected] == sp), integer(1))
  spp_line <- paste(c(length(spp), spp), collapse = " ")
  counts_line <- paste(spp_counts, collapse = " ")
  # Starting tree derived from IQ-TREE individual phylogeny collapsed to populations.
  start_tree <- guide_tree_newick
  phase_line <- paste(rep("1", length(spp)), collapse = " ")

  writeLines(c(
    "seed = 20260226",
    paste0("jobname = ", jobname),
    "seqfile = bpp_input_single_locus.phy",
    "Imapfile = bpp_imap.txt",
    "speciesdelimitation = 1 0 2",
    "speciestree = 1 0.1 0.1 0.1",
    "speciesmodelprior = 1",
    paste0("species&tree = ", spp_line),
    paste0("               ", counts_line),
    paste0("               ", start_tree),
    paste0("phase = ", phase_line),
    "usedata = 1",
    "nloci = 1",
    "cleandata = 0",
    "thetaprior = gamma 2 2000",
    "tauprior = gamma 2 1000",
    "finetune = 1",
    paste0("burnin = ", burnin),
    paste0("sampfreq = ", sampfreq),
    paste0("nsample = ", nsample)
  ), con = ctl_file)
}

attempts <- data.frame(
  Attempt = integer(), burnin = integer(), nsample = integer(), sampfreq = integer(),
  min_ESS = numeric(), status = integer(), stringsAsFactors = FALSE
)

plan <- data.frame(
  burnin = c(2000L, 5000L, 10000L, 15000L, 25000L, 40000L),
  nsample = c(10000L, 25000L, 50000L, 100000L, 160000L, 240000L),
  sampfreq = c(10L, 10L, 10L, 10L, 10L, 10L)
)
ess_target <- 200
ok <- FALSE

for (i in seq_len(nrow(plan))) {
  write_ctl(plan$burnin[i], plan$nsample[i], plan$sampfreq[i])
  if (file.exists(mcmc_file)) file.remove(mcmc_file)
  if (file.exists(file.path(bpp_dir, paste0(jobname, ".txt")))) file.remove(file.path(bpp_dir, paste0(jobname, ".txt")))
  if (file.exists(file.path(bpp_dir, paste0(jobname, ".SeedUsed")))) file.remove(file.path(bpp_dir, paste0(jobname, ".SeedUsed")))
  cmd <- sprintf("cd \"%s\" && \"%s\" --cfile \"%s\" > \"%s\" 2>&1",
                 bpp_dir, bpp_bin, ctl_file, run_log)
  status <- system(cmd)
  min_ess <- if (status == 0) calc_min_ess_from_mcmc(mcmc_file) else NA_real_
  attempts <- rbind(attempts, data.frame(
    Attempt = i, burnin = plan$burnin[i], nsample = plan$nsample[i], sampfreq = plan$sampfreq[i],
    min_ESS = min_ess, status = status, stringsAsFactors = FALSE
  ))
  if (status == 0 && is.finite(min_ess) && min_ess >= ess_target) {
    ok <- TRUE
    break
  }
}

write.table(attempts, file.path(tab_dir, "bpp_ess_by_attempt.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(
  BPP_binary = bpp_bin,
  Status = if (nrow(attempts)) tail(attempts$status, 1) else 127,
  Message = if (ok) paste0("BPP completed with ESS >= ", ess_target) else "BPP did not reach ESS target",
  Control_file = ctl_file,
  Sequence_file = seq_file,
  Imap_file = imap_file,
  stringsAsFactors = FALSE
), file.path(tab_dir, "bpp_run_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "BPP reconfigured run summary",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("IQ-TREE guide source: ", iqtree_contree),
  paste0("IQ-TREE-derived population guide tree: ", guide_tree_newick),
  paste0("Samples: ", length(selected)),
  paste0("SNPs: ", nrow(G)),
  paste0("ESS target: ", ess_target),
  paste0("Reached target: ", ok)
), con = file.path(tab_dir, "bpp_summary.txt"))
