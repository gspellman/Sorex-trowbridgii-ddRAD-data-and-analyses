#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
source_run <- file.path(down_dir, "eems_analysis_rerun_2026-07-29")
run_dir <- file.path(down_dir, "eems_analysis_ess200_2026-07-29")
data_dir <- file.path(run_dir, "data")
grid_dir <- file.path(run_dir, "fixed_grid")
mcmc_root <- file.path(run_dir, "mcmc")
log_dir <- file.path(run_dir, "logs")
table_dir <- file.path(run_dir, "tables")

for (d in c(run_dir, data_dir, grid_dir, mcmc_root, log_dir, table_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

source_prefix <- file.path(source_run, "data", "trow_eems")
data_prefix <- file.path(data_dir, "trow_eems")
grid_prefix <- file.path(grid_dir, "fixed_250_grid")
grid_source <- file.path(source_run, "mcmc", "chain2_demes250")
runeems <- file.path(down_dir, "eems-master", "runeems_snps", "src", "runeems_snps")
boost_lib <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/.envs/treemix_x86/lib"

required <- c(
  paste0(source_prefix, c(".coord", ".diffs", ".outer", ".order")),
  file.path(grid_source, c("demes.txt", "edges.txt")),
  runeems
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required EEMS files:\n", paste(missing, collapse = "\n"))

for (ext in c(".coord", ".diffs", ".outer", ".order")) {
  file.copy(paste0(source_prefix, ext), paste0(data_prefix, ext), overwrite = TRUE)
}
file.copy(file.path(grid_source, "demes.txt"), paste0(grid_prefix, ".demes"), overwrite = TRUE)
file.copy(file.path(grid_source, "edges.txt"), paste0(grid_prefix, ".edges"), overwrite = TRUE)

samples <- readLines(paste0(data_prefix, ".order"), warn = FALSE)
samples <- samples[nzchar(samples)]
if ("MVZ216210" %in% samples) stop("Excluded sample MVZ216210 is present")
if (length(samples) != 98L) stop("Expected 98 EEMS samples; found ", length(samples))

coords <- as.matrix(read.table(paste0(data_prefix, ".coord"), header = FALSE))
diffs <- as.matrix(read.table(paste0(data_prefix, ".diffs"), header = FALSE))
if (!all(dim(diffs) == c(98L, 98L))) stop("Dissimilarity matrix dimensions are incorrect")
if (max(abs(diffs - t(diffs))) > 1e-8) stop("Dissimilarity matrix is not symmetric")
if (max(abs(diag(diffs))) > 1e-8) stop("Dissimilarity matrix diagonal is nonzero")
if (nrow(coords) != 98L) stop("Coordinate count is incorrect")

grid_demes <- read.table(paste0(grid_prefix, ".demes"), header = FALSE)
grid_edges <- read.table(paste0(grid_prefix, ".edges"), header = FALSE)
if (nrow(grid_demes) != 210L) stop("Expected 210 clipped cells in fixed grid")
if (ncol(grid_edges) != 2L) stop("Fixed-grid edge file must contain two columns")

chain_plan <- data.frame(
  Chain = 1:4,
  Seed = 20260732:20260735,
  Nominal_nDemes = 250L,
  Actual_grid_cells = nrow(grid_demes),
  Iterations = 7000000L,
  Burnin = 1000000L,
  Thin = 999L
)
write.table(
  chain_plan,
  file.path(table_dir, "eems_ess200_chain_plan.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

run_chain <- function(i) {
  chain <- chain_plan$Chain[i]
  mcmc_dir <- file.path(mcmc_root, sprintf("chain%d_fixed250", chain))
  params <- file.path(run_dir, sprintf("params_chain%d.ini", chain))
  log_file <- file.path(log_dir, sprintf("chain%d.log", chain))
  if (dir.exists(mcmc_dir)) unlink(mcmc_dir, recursive = TRUE)

  writeLines(c(
    paste0("datapath = ", data_prefix),
    paste0("mcmcpath = ", mcmc_dir),
    paste0("gridpath = ", grid_prefix),
    "nIndiv = 98",
    "nSites = 14578",
    "nDemes = 250",
    "diploid = true",
    paste0("numMCMCIter = ", chain_plan$Iterations[i]),
    paste0("numBurnIter = ", chain_plan$Burnin[i]),
    paste0("numThinIter = ", chain_plan$Thin[i]),
    "mEffctProposalS2 = 0.20",
    "mrateMuProposalS2 = 0.005",
    "qSeedsProposalS2 = 0.06",
    "mSeedsProposalS2 = 0.01",
    "qEffctProposalS2 = 0.001"
  ), params)

  started <- Sys.time()
  status <- system2(
    "arch",
    args = c(
      "-x86_64", "env",
      shQuote(paste0("DYLD_FALLBACK_LIBRARY_PATH=", boost_lib)),
      shQuote(runeems),
      "--params", shQuote(params),
      "--seed", chain_plan$Seed[i]
    ),
    stdout = log_file,
    stderr = log_file
  )

  data.frame(
    Chain = chain,
    Seed = chain_plan$Seed[i],
    Return_code = status,
    Runtime_hours = as.numeric(difftime(Sys.time(), started, units = "hours")),
    MCMC_directory = mcmc_dir,
    Log = log_file
  )
}

results <- parallel::mclapply(1:4, run_chain, mc.cores = 4L)
run_summary <- do.call(rbind, results)
write.table(
  run_summary,
  file.path(table_dir, "eems_ess200_run_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
if (any(run_summary$Return_code != 0)) stop("One or more long EEMS chains failed")

writeLines(
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  file.path(run_dir, "LONG_CHAINS_COMPLETE.txt")
)
message("Four-chain ESS-targeted EEMS run completed: ", run_dir)
