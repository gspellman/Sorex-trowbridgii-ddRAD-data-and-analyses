#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
out_dir <- file.path(base_dir, "downstream_analysis_2026-02-23", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pipeline_log <- file.path(base_dir, "pipeline.log")
gstacks_log <- file.path(base_dir, "gstacks.log")
pop_log <- file.path(base_dir, "populations.log")
popsum_file <- file.path(base_dir, "populations.sumstats_summary.tsv")
popmap_file <- file.path(base_dir, "popmap.tsv")
vcf_file <- file.path(base_dir, "populations.snps.vcf")

get_first <- function(lines, pattern) {
  idx <- grep(pattern, lines, perl = TRUE)
  if (length(idx) == 0) return(NA_character_)
  lines[idx[1]]
}

extract_num <- function(line, pattern) {
  if (is.na(line)) return(NA_real_)
  m <- regexec(pattern, line, perl = TRUE)
  r <- regmatches(line, m)[[1]]
  if (length(r) < 2) return(NA_real_)
  as.numeric(gsub(",", "", r[2]))
}

extract_chr <- function(line, pattern) {
  if (is.na(line)) return(NA_character_)
  m <- regexec(pattern, line, perl = TRUE)
  r <- regmatches(line, m)[[1]]
  if (length(r) < 2) return(NA_character_)
  r[2]
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}

write_md_table <- function(df, file, title, caption = NULL) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste0("# ", title), con)
  writeLines("", con)
  if (!is.null(caption)) {
    writeLines(caption, con)
    writeLines("", con)
  }
  hdr <- paste(names(df), collapse = " | ")
  sep <- paste(rep("---", ncol(df)), collapse = " | ")
  writeLines(paste0("| ", hdr, " |"), con)
  writeLines(paste0("| ", sep, " |"), con)
  for (i in seq_len(nrow(df))) {
    row_txt <- paste(as.character(df[i, ]), collapse = " | ")
    writeLines(paste0("| ", row_txt, " |"), con)
  }
}

pipeline_lines <- readLines(pipeline_log)
gstacks_lines <- readLines(gstacks_log)
pop_lines <- readLines(pop_log)

run_start <- extract_chr(get_first(pipeline_lines, "^\\[[^]]+\\]"), "^\\[([^]]+)\\]")
run_end <- extract_chr(tail(grep("Copying completed outputs", pipeline_lines, value = TRUE), 1), "^\\[([^]]+)\\]")

samples_line <- get_first(pop_lines, "The population map contained")
n_samples <- extract_num(samples_line, "contained ([0-9,]+) samples")
n_pops <- extract_num(samples_line, "samples, ([0-9,]+) population")

reads_line <- get_first(gstacks_lines, "^Read [0-9,]+ BAM records")
reads_total <- extract_num(reads_line, "Read ([0-9,]+) BAM records")
kept_line <- get_first(gstacks_lines, "kept [0-9,]+ primary alignments")
kept_primary <- extract_num(kept_line, "kept ([0-9,]+) primary alignments")
kept_pct <- extract_num(kept_line, "\\(([0-9.]+)%\\)")
loci_built_line <- get_first(gstacks_lines, "^Built [0-9,]+ loci")
loci_built <- extract_num(loci_built_line, "Built ([0-9,]+) loci")
genotyped_line <- get_first(gstacks_lines, "^Genotyped [0-9,]+ loci")
loci_genotyped <- extract_num(genotyped_line, "Genotyped ([0-9,]+) loci")
coverage_line <- get_first(gstacks_lines, "effective per-sample coverage")
cov_mean <- extract_num(coverage_line, "mean=([0-9.]+)x")
cov_min <- extract_num(coverage_line, "min=([0-9.]+)x")
cov_max <- extract_num(coverage_line, "max=([0-9.]+)x")

filter_samples_per_pop <- extract_num(get_first(pop_lines, "Percent samples limit per population"), ": ([0-9.]+)$")
filter_min_pops <- extract_num(get_first(pop_lines, "Locus Population limit"), ": ([0-9.]+)$")
filter_min_maf <- extract_num(get_first(pop_lines, "Minor allele frequency cutoff"), ": ([0-9.]+)$")

removed_line <- get_first(pop_lines, "^Removed [0-9,]+ loci")
removed_loci <- extract_num(removed_line, "Removed ([0-9,]+) loci")
from_loci <- extract_num(removed_line, "from ([0-9,]+) loci")
kept_line2 <- get_first(pop_lines, "^Kept [0-9,]+ loci")
kept_loci <- extract_num(kept_line2, "Kept ([0-9,]+) loci")
kept_sites <- extract_num(kept_line2, "composed of ([0-9,]+) sites")
filtered_sites <- extract_num(kept_line2, "; ([0-9,]+) of those sites were filtered")
variant_sites <- extract_num(kept_line2, ", ([0-9,]+) variant sites remained")

stacks_tbl <- data.frame(
  Metric = c(
    "Run start",
    "Run end",
    "Reference genome",
    "Input samples",
    "Populations",
    "Total BAM records",
    "Primary alignments retained",
    "Retained primary alignments (%)",
    "Loci built (gstacks)",
    "Loci genotyped (gstacks)",
    "Per-sample coverage mean (x)",
    "Per-sample coverage min (x)",
    "Per-sample coverage max (x)",
    "populations: min samples per pop",
    "populations: min populations per locus",
    "populations: min MAF",
    "Loci removed by filters",
    "Total loci evaluated",
    "Loci retained",
    "Total sites in retained loci",
    "Sites filtered",
    "Variant sites retained"
  ),
  Value = c(
    run_start,
    run_end,
    "Sorex_ornatus (GCA_041430635.1)",
    as.integer(n_samples),
    as.integer(n_pops),
    as.integer(reads_total),
    as.integer(kept_primary),
    fmt(kept_pct, 1),
    as.integer(loci_built),
    as.integer(loci_genotyped),
    fmt(cov_mean, 1),
    fmt(cov_min, 1),
    fmt(cov_max, 1),
    fmt(filter_samples_per_pop, 2),
    fmt(filter_min_pops, 0),
    fmt(filter_min_maf, 2),
    as.integer(removed_loci),
    as.integer(from_loci),
    as.integer(kept_loci),
    as.integer(kept_sites),
    as.integer(filtered_sites),
    as.integer(variant_sites)
  ),
  stringsAsFactors = FALSE
)

write.table(
  stacks_tbl,
  file.path(out_dir, "Table1_STACKS_run_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write_md_table(
  stacks_tbl,
  file.path(out_dir, "Table1_STACKS_run_summary.md"),
  "Table 1. STACKS Pipeline Summary",
  "Reference-based STACKS v2.68 analysis of trowdata ddRAD reads against Sorex ornatus reference assembly."
)

popmap <- read.delim(popmap_file, header = FALSE)
colnames(popmap) <- c("Sample", "Population")
pop_sizes <- as.data.frame(table(popmap$Population), stringsAsFactors = FALSE)
colnames(pop_sizes) <- c("Population", "N_samples")

popsum_lines <- readLines(popsum_file)
idx_all <- grep("^# All positions", popsum_lines)
if (length(idx_all) != 1) stop("Could not find '# All positions' block in populations.sumstats_summary.tsv")
hdr <- strsplit(sub("^# ", "", popsum_lines[idx_all + 1]), "\t")[[1]]
rows <- popsum_lines[(idx_all + 2):length(popsum_lines)]
rows <- rows[nzchar(rows) & !grepl("^#", rows)]
allpos <- read.delim(text = paste(rows, collapse = "\n"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(allpos) <- hdr
num_cols <- setdiff(colnames(allpos), "Pop ID")
for (cc in num_cols) allpos[[cc]] <- as.numeric(allpos[[cc]])
colnames(allpos)[colnames(allpos) == "Pop ID"] <- "Population"

h <- readLines(vcf_file, n = 5000)
ch <- h[grepl("^#CHROM", h)]
if (length(ch) != 1) stop("Unable to parse VCF header line.")
vcf_samples <- strsplit(ch, "\t")[[1]][10:length(strsplit(ch, "\t")[[1]])]

vcf <- read.table(vcf_file, sep = "\t", comment.char = "#", header = FALSE, quote = "", fill = TRUE)
colnames(vcf)[1:9] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
colnames(vcf)[10:ncol(vcf)] <- vcf_samples

is_bial <- nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1 & !grepl(",", vcf$ALT)
vcf <- vcf[is_bial, , drop = FALSE]
vcf_samples <- intersect(vcf_samples, popmap$Sample)
GT <- as.matrix(vcf[, vcf_samples, drop = FALSE])
GT <- sub(":.*", "", GT)

n_sites <- nrow(GT)
dosage <- matrix(NA_real_, nrow = nrow(GT), ncol = ncol(GT), dimnames = list(NULL, colnames(GT)))
dosage[GT %in% c("0/0", "0|0")] <- 0
dosage[GT %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
dosage[GT %in% c("1/1", "1|1")] <- 2

pop_levels <- c("North", "North_Coast", "Sierra_1", "Sierra_2", "Sierra_3", "South_Coast")

calc_pop <- function(pop_name) {
  s <- popmap$Sample[popmap$Population == pop_name]
  s <- intersect(s, colnames(dosage))
  d <- dosage[, s, drop = FALSE]
  called <- is.finite(d)
  site_called_n <- rowSums(called)

  polym <- rep(FALSE, nrow(d))
  ok <- site_called_n > 1
  if (any(ok)) {
    mn <- apply(d[ok, , drop = FALSE], 1, function(v) min(v[is.finite(v)]))
    mx <- apply(d[ok, , drop = FALSE], 1, function(v) max(v[is.finite(v)]))
    polym[ok] <- mn < mx
  }

  bial_sites <- site_called_n > 0

  sample_missing <- 1 - colMeans(called)
  mean_called_bial_snp <- mean(colSums(called[bial_sites, , drop = FALSE]))
  sd_called_bial_snp <- sd(colSums(called[bial_sites, , drop = FALSE]))
  mean_var_snp_per_sample <- mean(colSums(d[polym, , drop = FALSE] > 0, na.rm = TRUE))
  sd_var_snp_per_sample <- sd(colSums(d[polym, , drop = FALSE] > 0, na.rm = TRUE))

  af <- rep(NA_real_, sum(polym))
  if (sum(polym) > 0) {
    pdat <- d[polym, , drop = FALSE]
    af <- apply(pdat, 1, function(v) {
      vv <- v[is.finite(v)]
      if (length(vv) == 0) return(NA_real_)
      sum(vv) / (2 * length(vv))
    })
  }
  maf <- pmin(af, 1 - af)
  obs_het <- rep(NA_real_, sum(polym))
  if (sum(polym) > 0) {
    pdat <- d[polym, , drop = FALSE]
    obs_het <- apply(pdat, 1, function(v) {
      vv <- v[is.finite(v)]
      if (length(vv) == 0) return(NA_real_)
      mean(vv == 1)
    })
  }

  data.frame(
    Population = pop_name,
    N_samples = length(s),
    Total_biallelic_SNP_sites = sum(bial_sites),
    Polymorphic_SNP_sites = sum(polym),
    Percent_polymorphic_sites = 100 * sum(polym) / max(1, sum(bial_sites)),
    Mean_missing_genotypes_pct = 100 * (1 - mean(called)),
    Mean_sample_missingness_pct = 100 * mean(sample_missing),
    SD_sample_missingness_pct = 100 * sd(sample_missing),
    Mean_called_biallelic_SNPs_per_sample = mean_called_bial_snp,
    SD_called_biallelic_SNPs_per_sample = sd_called_bial_snp,
    Mean_variable_SNP_genotypes_per_sample = mean_var_snp_per_sample,
    SD_variable_SNP_genotypes_per_sample = sd_var_snp_per_sample,
    Mean_MAF_polymorphic_sites = mean(maf, na.rm = TRUE),
    Mean_observed_heterozygosity_polymorphic = mean(obs_het, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

pop_stats <- do.call(rbind, lapply(pop_levels, calc_pop))

tab2 <- merge(pop_stats, allpos[, c("Population", "Private", "Num_Indv", "Pi", "Fis", "%Polymorphic_Loci")], by = "Population", all.x = TRUE)
tab2 <- merge(tab2, pop_sizes, by = "Population", all.x = TRUE, suffixes = c("", "_popmap"))
if ("N_samples_popmap" %in% names(tab2)) {
  tab2$N_samples <- tab2$N_samples_popmap
  tab2$N_samples_popmap <- NULL
}
tab2$STACKS_Mean_called_individuals_per_locus <- tab2$Num_Indv
tab2$STACKS_locus_missingness_pct <- 100 * (1 - tab2$Num_Indv / tab2$N_samples)

tab2 <- tab2[, c(
  "Population",
  "N_samples",
  "Total_biallelic_SNP_sites",
  "Polymorphic_SNP_sites",
  "Percent_polymorphic_sites",
  "Mean_missing_genotypes_pct",
  "Mean_sample_missingness_pct",
  "SD_sample_missingness_pct",
  "Mean_called_biallelic_SNPs_per_sample",
  "SD_called_biallelic_SNPs_per_sample",
  "Mean_variable_SNP_genotypes_per_sample",
  "SD_variable_SNP_genotypes_per_sample",
  "Mean_MAF_polymorphic_sites",
  "Mean_observed_heterozygosity_polymorphic",
  "Private",
  "Pi",
  "Fis",
  "STACKS_Mean_called_individuals_per_locus",
  "STACKS_locus_missingness_pct"
)]

for (cc in setdiff(names(tab2), c("Population", "N_samples", "Total_biallelic_SNP_sites", "Polymorphic_SNP_sites", "Private"))) {
  tab2[[cc]] <- round(tab2[[cc]], 4)
}

write.table(
  tab2,
  file.path(out_dir, "Table2_population_SNP_variation_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

md2 <- tab2
md2$Percent_polymorphic_sites <- paste0(sprintf("%.2f", md2$Percent_polymorphic_sites), "%")
md2$Mean_missing_genotypes_pct <- paste0(sprintf("%.2f", md2$Mean_missing_genotypes_pct), "%")
md2$Mean_sample_missingness_pct <- paste0(sprintf("%.2f", md2$Mean_sample_missingness_pct), "%")
md2$SD_sample_missingness_pct <- paste0(sprintf("%.2f", md2$SD_sample_missingness_pct), "%")
md2$STACKS_locus_missingness_pct <- paste0(sprintf("%.2f", md2$STACKS_locus_missingness_pct), "%")

write_md_table(
  md2,
  file.path(out_dir, "Table2_population_SNP_variation_summary.md"),
  "Table 2. Per-Population SNP Variation and Missingness Summary",
  "Derived from populations.snps.vcf (biallelic SNPs) and populations.sumstats_summary.tsv."
)

cat("Wrote:\n")
cat(file.path(out_dir, "Table1_STACKS_run_summary.tsv"), "\n")
cat(file.path(out_dir, "Table1_STACKS_run_summary.md"), "\n")
cat(file.path(out_dir, "Table2_population_SNP_variation_summary.tsv"), "\n")
cat(file.path(out_dir, "Table2_population_SNP_variation_summary.md"), "\n")
