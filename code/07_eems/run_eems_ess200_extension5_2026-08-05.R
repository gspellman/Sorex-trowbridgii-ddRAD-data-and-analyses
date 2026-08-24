#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
source_run <- file.path(down_dir, "eems_analysis_ess200_extension4_2026-08-04")
input_run <- file.path(down_dir, "eems_analysis_ess200_2026-07-29")
run_dir <- file.path(down_dir, "eems_analysis_ess200_extension5_2026-08-05")
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

source_paths <- file.path(
  source_run,
  "mcmc",
  paste0("chain", 1:4, "_fixed250_extension4")
)
restart_files <- c(
  "lastqtiles.txt", "lastmtiles.txt", "lastthetas.txt", "lastdfpars.txt",
  "lastqhyper.txt", "lastmhyper.txt", "lastqeffct.txt", "lastmeffct.txt",
  "lastqseeds.txt", "lastmseeds.txt"
)
required <- c(
  paste0(data_prefix, c(".coord", ".diffs", ".outer", ".order")),
  paste0(grid_prefix, c(".demes", ".edges")),
  unlist(lapply(source_paths, function(p) file.path(p, restart_files))),
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
  Tuning_seed = 20260810:20260813,
  Production_seed = 20260820:20260823,
  Tuning_iterations = 2000000L,
  Tuning_burnin = 1000000L,
  Production_iterations = 20000000L,
  Production_burnin = 2000000L,
  Thin = 999L,
  Production_retained_states = 18000L,
  Source_path = source_paths
)
write.table(
  chain_plan,
  file.path(table_dir, "eems_extension5_chain_plan.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write_params <- function(path, mcmc_path, previous_path, iterations, burnin) {
  writeLines(c(
    paste0("datapath = ", data_prefix),
    paste0("mcmcpath = ", mcmc_path),
    paste0("prevpath = ", previous_path),
    paste0("gridpath = ", grid_prefix),
    "nIndiv = 98",
    "nSites = 14578",
    "nDemes = 250",
    "diploid = true",
    paste0("numMCMCIter = ", iterations),
    paste0("numBurnIter = ", burnin),
    "numThinIter = 999",
    "qVoronoiPr = 0.30",
    "qSeedsProposalS2 = 0.025",
    "qEffctProposalS2 = 0.001",
    "mSeedsProposalS2 = 0.01",
    "mEffctProposalS2 = 0.15",
    "mrateMuProposalS2 = 0.0015"
  ), path)
}

run_eems <- function(params, seed, log_file) {
  started <- Sys.time()
  status <- system2(
    "arch",
    args = c(
      "-x86_64", "env",
      shQuote(paste0("DYLD_FALLBACK_LIBRARY_PATH=", boost_lib)),
      shQuote(runeems),
      "--params", shQuote(params),
      "--seed", seed
    ),
    stdout = log_file,
    stderr = log_file
  )
  c(
    Return_code = status,
    Runtime_hours = as.numeric(difftime(Sys.time(), started, units = "hours"))
  )
}

run_chain <- function(i) {
  chain <- chain_plan$Chain[i]
  tuning_dir <- file.path(mcmc_root, sprintf("chain%d_tuning", chain))
  production_dir <- file.path(mcmc_root, sprintf("chain%d_fixed250_extension5", chain))
  tuning_params <- file.path(run_dir, sprintf("params_chain%d_tuning.ini", chain))
  production_params <- file.path(run_dir, sprintf("params_chain%d_production.ini", chain))
  tuning_log <- file.path(log_dir, sprintf("chain%d_tuning.log", chain))
  production_log <- file.path(log_dir, sprintf("chain%d_production.log", chain))

  for (d in c(tuning_dir, production_dir)) {
    if (dir.exists(d)) unlink(d, recursive = TRUE)
  }

  write_params(
    tuning_params, tuning_dir, source_paths[i],
    chain_plan$Tuning_iterations[i], chain_plan$Tuning_burnin[i]
  )
  tuning <- run_eems(tuning_params, chain_plan$Tuning_seed[i], tuning_log)
  if (tuning[["Return_code"]] != 0) {
    return(data.frame(
      Chain = chain, Tuning_return_code = tuning[["Return_code"]],
      Tuning_runtime_hours = tuning[["Runtime_hours"]],
      Production_return_code = NA, Production_runtime_hours = NA
    ))
  }

  missing_tuning <- file.path(tuning_dir, restart_files)
  if (any(!file.exists(missing_tuning))) stop("Tuning restart files missing for chain ", chain)

  write_params(
    production_params, production_dir, tuning_dir,
    chain_plan$Production_iterations[i], chain_plan$Production_burnin[i]
  )
  production <- run_eems(
    production_params, chain_plan$Production_seed[i], production_log
  )

  data.frame(
    Chain = chain,
    Tuning_return_code = tuning[["Return_code"]],
    Tuning_runtime_hours = tuning[["Runtime_hours"]],
    Production_return_code = production[["Return_code"]],
    Production_runtime_hours = production[["Runtime_hours"]]
  )
}

results <- parallel::mclapply(1:4, run_chain, mc.cores = 4L)
run_summary <- do.call(rbind, results)
write.table(
  run_summary,
  file.path(table_dir, "eems_extension5_run_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
if (any(run_summary$Tuning_return_code != 0) ||
    any(run_summary$Production_return_code != 0)) {
  stop("One or more EEMS tuning or production chains failed")
}

writeLines(
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  file.path(run_dir, "EXTENSION5_COMPLETE.txt")
)
message("Fifth EEMS continuation completed: ", run_dir)
