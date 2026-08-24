#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
source_run <- file.path(down_dir, "eems_analysis_ess200_extension1_2026-07-29")
input_run <- file.path(down_dir, "eems_analysis_ess200_2026-07-29")
run_dir <- file.path(down_dir, "eems_analysis_ess200_extension2_2026-07-30")
data_prefix <- file.path(input_run, "data", "trow_eems")
grid_prefix <- file.path(input_run, "fixed_grid", "fixed_250_grid")
mcmc_root <- file.path(run_dir, "mcmc")
log_dir <- file.path(run_dir, "logs")
table_dir <- file.path(run_dir, "tables")
runeems <- file.path(down_dir, "eems-master", "runeems_snps", "src", "runeems_snps")
boost_lib <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/.envs/treemix_x86/lib"

for (d in c(run_dir, mcmc_root, log_dir, table_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

previous_paths <- file.path(
  source_run,
  "mcmc",
  paste0("chain", 1:4, "_fixed250_extension1")
)
required <- c(
  paste0(data_prefix, c(".coord", ".diffs", ".outer", ".order")),
  paste0(grid_prefix, c(".demes", ".edges")),
  unlist(lapply(previous_paths, function(p) file.path(
    p,
    c(
      "lastqtiles.txt", "lastmtiles.txt", "lastthetas.txt", "lastdfpars.txt",
      "lastqhyper.txt", "lastmhyper.txt", "lastqeffct.txt", "lastmeffct.txt",
      "lastqseeds.txt", "lastmseeds.txt"
    )
  ))),
  runeems
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing continuation inputs:\n", paste(missing, collapse = "\n"))

samples <- readLines(paste0(data_prefix, ".order"), warn = FALSE)
samples <- samples[nzchar(samples)]
if (length(samples) != 98L || "MVZ216210" %in% samples) {
  stop("Sample exclusion or sample count validation failed")
}

chain_plan <- data.frame(
  Chain = 1:4,
  Seed = 20260740:20260743,
  Iterations = 8000000L,
  Burnin = 1000000L,
  Thin = 999L,
  Retained_states = 7000L,
  Previous_path = previous_paths
)
write.table(
  chain_plan,
  file.path(table_dir, "eems_extension2_chain_plan.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

run_chain <- function(i) {
  chain <- chain_plan$Chain[i]
  mcmc_dir <- file.path(mcmc_root, sprintf("chain%d_fixed250_extension2", chain))
  params <- file.path(run_dir, sprintf("params_chain%d.ini", chain))
  log_file <- file.path(log_dir, sprintf("chain%d.log", chain))
  if (dir.exists(mcmc_dir)) unlink(mcmc_dir, recursive = TRUE)

  writeLines(c(
    paste0("datapath = ", data_prefix),
    paste0("mcmcpath = ", mcmc_dir),
    paste0("prevpath = ", previous_paths[i]),
    paste0("gridpath = ", grid_prefix),
    "nIndiv = 98",
    "nSites = 14578",
    "nDemes = 250",
    "diploid = true",
    paste0("numMCMCIter = ", chain_plan$Iterations[i]),
    paste0("numBurnIter = ", chain_plan$Burnin[i]),
    paste0("numThinIter = ", chain_plan$Thin[i]),
    "qVoronoiPr = 0.45",
    "qSeedsProposalS2 = 0.025",
    "qEffctProposalS2 = 0.001",
    "mSeedsProposalS2 = 0.01",
    "mEffctProposalS2 = 0.20",
    "mrateMuProposalS2 = 0.0015"
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
    Previous_path = previous_paths[i],
    MCMC_directory = mcmc_dir,
    Log = log_file
  )
}

results <- parallel::mclapply(1:4, run_chain, mc.cores = 4L)
run_summary <- do.call(rbind, results)
write.table(
  run_summary,
  file.path(table_dir, "eems_extension2_run_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
if (any(run_summary$Return_code != 0)) stop("One or more EEMS continuation chains failed")

writeLines(
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  file.path(run_dir, "EXTENSION2_COMPLETE.txt")
)
message("Second tuned EEMS continuation completed: ", run_dir)
