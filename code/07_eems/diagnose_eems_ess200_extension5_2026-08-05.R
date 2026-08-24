#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

run_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23/eems_analysis_ess200_extension5_2026-08-05"
table_dir <- file.path(run_dir, "tables")
paths <- file.path(
  run_dir, "mcmc", paste0("chain", 1:4, "_fixed250_extension5")
)

if (!all(file.exists(file.path(paths, "eemsrun.txt")))) {
  stop("Fifth-continuation production chains are incomplete")
}

effective_size <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 4L || var(x) == 0) return(NA_real_)
  rho <- as.numeric(acf(
    x, plot = FALSE, lag.max = min(n - 1L, floor(10 * log10(n)))
  )$acf)[-1]
  pair_count <- floor(length(rho) / 2)
  pair_sums <- rho[seq_len(pair_count) * 2 - 1] + rho[seq_len(pair_count) * 2]
  first_negative <- which(pair_sums < 0)[1]
  if (!is.na(first_negative)) pair_sums <- pair_sums[seq_len(first_negative - 1)]
  tau <- 1 + 2 * sum(pair_sums)
  n / max(tau, 1)
}

split_chains <- function(chains) {
  halves <- unlist(lapply(chains, function(x) {
    n <- floor(length(x) / 2)
    list(x[seq_len(n)], tail(x, n))
  }), recursive = FALSE)
  n <- min(vapply(halves, length, integer(1)))
  lapply(halves, function(x) tail(x, n))
}

classic_rhat <- function(chains) {
  n <- min(vapply(chains, length, integer(1)))
  chains <- lapply(chains, function(x) tail(x, n))
  within <- mean(vapply(chains, var, numeric(1)))
  if (!is.finite(within) || within == 0) return(NA_real_)
  between <- n * var(vapply(chains, mean, numeric(1)))
  sqrt((((n - 1) / n) * within + between / n) / within)
}

rank_normalize <- function(x) {
  n <- length(x)
  stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4))
}

rank_split_rhat <- function(chains) {
  halves <- split_chains(chains)
  lengths <- vapply(halves, length, integer(1))
  pooled <- unlist(halves, use.names = FALSE)
  ranked <- rank_normalize(pooled)
  ranked_halves <- split(ranked, rep(seq_along(halves), lengths))
  bulk <- classic_rhat(ranked_halves)

  folded <- abs(pooled - median(pooled))
  folded_ranked <- rank_normalize(folded)
  folded_halves <- split(folded_ranked, rep(seq_along(halves), lengths))
  folded_rhat <- classic_rhat(folded_halves)

  c(Bulk_Rhat = bulk, Folded_Rhat = folded_rhat,
    Rank_normalized_Rhat = max(bulk, folded_rhat, na.rm = TRUE))
}

file_specs <- list(
  mcmcpilogl.txt = c("LogPrior", "LogLikelihood"),
  mcmcmhyper.txt = c("MigrationMean", "MigrationVariance"),
  mcmcqhyper.txt = c("DiversityMean", "DiversityVariance"),
  mcmcthetas.txt = c("Theta", "DegreesFreedom"),
  mcmcmtiles.txt = "MigrationTiles",
  mcmcqtiles.txt = "DiversityTiles"
)

all_values <- list()
ess_rows <- list()
half_rows <- list()
for (filename in names(file_specs)) {
  parameter_names <- file_specs[[filename]]
  per_chain <- lapply(paths, function(p) {
    x <- read.table(file.path(p, filename), header = FALSE)
    colnames(x) <- parameter_names
    x
  })
  for (parameter in parameter_names) {
    key <- paste(filename, parameter, sep = "::")
    all_values[[key]] <- lapply(per_chain, `[[`, parameter)
    for (i in 1:4) {
      values <- per_chain[[i]][[parameter]]
      midpoint <- floor(length(values) / 2)
      ess_rows[[length(ess_rows) + 1L]] <- data.frame(
        Chain = i, File = filename, Parameter = parameter,
        Retained_states = length(values), ESS = effective_size(values),
        Mean = mean(values), SD = sd(values)
      )
      half_rows[[length(half_rows) + 1L]] <- data.frame(
        Chain = i, File = filename, Parameter = parameter,
        First_half_mean = mean(values[seq_len(midpoint)]),
        Second_half_mean = mean(values[(midpoint + 1):length(values)])
      )
    }
  }
}

ess <- do.call(rbind, ess_rows)
halves <- do.call(rbind, half_rows)
rhat <- do.call(rbind, lapply(names(all_values), function(key) {
  pieces <- strsplit(key, "::", fixed = TRUE)[[1]]
  values <- unlist(all_values[[key]], use.names = FALSE)
  if (var(values) == 0) {
    metrics <- c(Bulk_Rhat = NA, Folded_Rhat = NA, Rank_normalized_Rhat = NA)
  } else {
    metrics <- rank_split_rhat(all_values[[key]])
  }
  data.frame(File = pieces[1], Parameter = pieces[2], t(metrics),
             check.names = FALSE)
}))

write.table(ess, file.path(table_dir, "eems_extension5_parameter_ess.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(rhat, file.path(table_dir, "eems_extension5_rank_normalized_rhat.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(halves, file.path(table_dir, "eems_extension5_half_means.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

finite_ess <- ess[is.finite(ess$ESS), ]
finite_rhat <- rhat[is.finite(rhat$Rank_normalized_Rhat), ]
summary <- data.frame(
  Minimum_per_chain_ESS = min(finite_ess$ESS),
  Maximum_rank_normalized_split_Rhat = max(finite_rhat$Rank_normalized_Rhat),
  ESS_target = 200,
  Rhat_target = 1.01
)
summary$ESS_pass <- summary$Minimum_per_chain_ESS >= summary$ESS_target
summary$Rhat_pass <- summary$Maximum_rank_normalized_split_Rhat <= summary$Rhat_target
summary$Overall_pass <- summary$ESS_pass && summary$Rhat_pass
write.table(summary, file.path(table_dir, "eems_extension5_convergence_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

status_file <- if (summary$Overall_pass) {
  "CONVERGENCE_PASS.txt"
} else {
  "CONVERGENCE_EXTENSION_REQUIRED.txt"
}
writeLines(c(
  paste("Minimum per-chain ESS:", round(summary$Minimum_per_chain_ESS, 2)),
  paste("Maximum rank-normalized split R-hat:",
        round(summary$Maximum_rank_normalized_split_Rhat, 4))
), file.path(run_dir, status_file))
print(summary)
