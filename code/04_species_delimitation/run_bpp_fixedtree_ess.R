#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
run_dir <- file.path(down_dir, "species_delimitation_2026-02-25_bpp_fixedtree_ess")
bpp_dir <- file.path(run_dir, "bpp_fixedtree")
tab_dir <- file.path(run_dir, "tables")
log_dir <- file.path(run_dir, "logs")
for (d in c(run_dir, bpp_dir, tab_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

bpp_bin <- "/Users/gspellman/Downloads/bpp-4.8.7-macos-aarch64/bin/bpp"
if (!file.exists(bpp_bin)) stop("Missing BPP binary: ", bpp_bin)

src_seq <- file.path(down_dir, "species_delimitation_2026-02-25_ess_bpp_reconfigured", "bpp_attempt", "bpp_input_single_locus.phy")
src_imap <- file.path(down_dir, "species_delimitation_2026-02-25_ess_bpp_reconfigured", "bpp_attempt", "bpp_imap.txt")
if (!file.exists(src_seq) || !file.exists(src_imap)) stop("Source BPP input files missing")
file.copy(src_seq, file.path(bpp_dir, "bpp_input_single_locus.phy"), overwrite = TRUE)
file.copy(src_imap, file.path(bpp_dir, "bpp_imap.txt"), overwrite = TRUE)

ess_target <- 200
jobname <- "bpp_fixedtree"

calc_ess <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 50) return(NA_real_)
  ac <- stats::acf(x, lag.max = min(1000, floor(n / 2)), plot = FALSE, demean = TRUE)$acf[-1]
  pos <- as.numeric(ac)[as.numeric(ac) > 0]
  if (!length(pos)) return(n)
  tau <- 1 + 2 * sum(pos)
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  n / tau
}

min_ess <- function(mcmc_file) {
  d <- tryCatch(read.table(mcmc_file, header = TRUE, check.names = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) < 200) return(NA_real_)
  keep <- names(d)[vapply(d, is.numeric, logical(1))]
  keep <- setdiff(keep, c("Gen", "gen", "sample", "Sample"))
  if (!length(keep)) return(NA_real_)
  burn <- floor(0.2 * nrow(d))
  p <- d[(burn + 1):nrow(d), , drop = FALSE]
  ess <- vapply(keep, function(k) calc_ess(p[[k]]), numeric(1))
  suppressWarnings(min(ess, na.rm = TRUE))
}

write_ctl <- function(path, burnin, nsample, sampfreq) {
  writeLines(c(
    "seed = 20260226",
    paste0("jobname = ", jobname),
    "seqfile = bpp_input_single_locus.phy",
    "Imapfile = bpp_imap.txt",
    "speciesdelimitation = 0",
    "speciestree = 0",
    "speciesmodelprior = 1",
    "species&tree = 6 North North_Coast Sierra_1 Sierra_2 Sierra_3 South_Coast",
    "               2 2 2 2 2 2",
    "               ((North,North_Coast),((Sierra_1,Sierra_2),(Sierra_3,South_Coast)));",
    "phase = 1 1 1 1 1 1",
    "usedata = 1",
    "nloci = 1",
    "cleandata = 0",
    "thetaprior = gamma 2 2000",
    "tauprior = gamma 2 1000",
    "finetune = 1",
    "print = 1 0 0 0",
    paste0("burnin = ", burnin),
    paste0("sampfreq = ", sampfreq),
    paste0("nsample = ", nsample)
  ), con = path)
}

plan <- data.frame(
  burnin = c(5000L, 10000L, 20000L, 40000L, 60000L),
  nsample = c(50000L, 100000L, 200000L, 350000L, 500000L),
  sampfreq = c(10L, 10L, 10L, 10L, 10L)
)

attempts <- data.frame(
  Attempt = integer(), burnin = integer(), nsample = integer(), sampfreq = integer(),
  rows = integer(), min_ESS = numeric(), status = integer(), stringsAsFactors = FALSE
)

ctl_file <- file.path(bpp_dir, "bpp_fixedtree.ctl")
log_file <- file.path(log_dir, "bpp_fixedtree.log")
mcmc_file <- file.path(bpp_dir, paste0(jobname, ".mcmc.txt"))
ok <- FALSE

for (i in seq_len(nrow(plan))) {
  write_ctl(ctl_file, plan$burnin[i], plan$nsample[i], plan$sampfreq[i])
  for (f in c(mcmc_file, file.path(bpp_dir, paste0(jobname, ".txt")), file.path(bpp_dir, paste0(jobname, ".SeedUsed")))) {
    if (file.exists(f)) file.remove(f)
  }
  cmd <- sprintf("cd \"%s\" && \"%s\" --cfile \"%s\" > \"%s\" 2>&1", bpp_dir, bpp_bin, ctl_file, log_file)
  status <- system(cmd)
  d <- if (file.exists(mcmc_file)) tryCatch(read.table(mcmc_file, header = TRUE, check.names = FALSE), error = function(e) NULL) else NULL
  nrows <- if (is.null(d)) 0L else nrow(d)
  me <- if (status == 0 && !is.null(d)) min_ess(mcmc_file) else NA_real_
  attempts <- rbind(attempts, data.frame(
    Attempt = i, burnin = plan$burnin[i], nsample = plan$nsample[i], sampfreq = plan$sampfreq[i],
    rows = nrows, min_ESS = me, status = status, stringsAsFactors = FALSE
  ))
  if (status == 0 && is.finite(me) && me >= ess_target) {
    ok <- TRUE
    break
  }
}

write.table(attempts, file.path(tab_dir, "bpp_fixedtree_ess_by_attempt.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(
  BPP_binary = bpp_bin,
  Status = if (nrow(attempts)) tail(attempts$status, 1) else 127,
  ESS_target = ess_target,
  ESS_reached = ok,
  Message = if (ok) paste0("ESS target reached (>= ", ess_target, ")") else "ESS target not reached",
  Control_file = ctl_file,
  MCMC_file = mcmc_file,
  stringsAsFactors = FALSE
), file.path(tab_dir, "bpp_fixedtree_run_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

