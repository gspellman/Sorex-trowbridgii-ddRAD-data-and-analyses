#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

BASE <- '/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean'
OUT <- file.path(BASE, 'downstream_analysis_2026-02-23')
TAB <- file.path(OUT, 'tables')

vcf_file <- file.path(BASE, 'populations.snps.vcf')
pop_file <- file.path(BASE, 'popmap.tsv')
table2_file <- file.path(TAB, 'Table2_population_SNP_variation_summary.tsv')
demog_file <- file.path(TAB, 'demography_population_summary.tsv')
fst_file <- file.path(TAB, 'pairwise_fst_weir_cockerham.tsv')
fis_file <- file.path(TAB, 'population_fis_permutation.tsv')

out_tsv <- file.path(TAB, 'Table3_population_SNP_diversity_publication.tsv')
out_md <- file.path(TAB, 'Table3_population_SNP_diversity_publication.md')
out_dxy <- file.path(TAB, 'pairwise_dxy_matrix.tsv')

pop_order <- c('North', 'North_Coast', 'Sierra_1', 'Sierra_2', 'Sierra_3', 'South_Coast')

parse_gt <- function(gt) {
  gt <- sub(':.*$', '', gt)
  gt <- gsub('\\|', '/', gt)
  if (gt %in% c('./.', '.', '.|.')) return(NA_real_)
  if (gt == '0/0') return(0)
  if (gt %in% c('0/1', '1/0')) return(1)
  if (gt == '1/1') return(2)
  NA_real_
}

cat('[1/4] Loading metadata\n')
popmap <- read.table(pop_file, sep='\t', header=FALSE, col.names=c('Sample', 'Population'))
popmap <- subset(popmap, Population %in% pop_order)
pop_lookup <- setNames(popmap$Population, popmap$Sample)

cat('[2/4] Reading VCF dosage matrix\n')
con <- file(vcf_file, 'r')
on.exit(close(con), add=TRUE)

samples <- NULL
keep_idx <- NULL
rows <- list()
ri <- 0L

repeat {
  ln <- readLines(con, n=1)
  if (length(ln) == 0) break
  if (startsWith(ln, '##')) next
  if (startsWith(ln, '#CHROM')) {
    parts <- strsplit(ln, '\t', fixed=TRUE)[[1]]
    all_samples <- parts[10:length(parts)]
    samples <- all_samples[all_samples %in% names(pop_lookup)]
    keep_idx <- match(samples, all_samples)
    next
  }
  p <- strsplit(ln, '\t', fixed=TRUE)[[1]]
  ref <- p[4]
  alt <- p[5]
  if (nchar(ref) != 1 || nchar(alt) != 1 || grepl(',', alt, fixed=TRUE)) next
  g <- p[10:length(p)]
  vals <- vapply(keep_idx, function(j) parse_gt(g[j]), numeric(1))
  ri <- ri + 1L
  rows[[ri]] <- vals
}

G <- do.call(rbind, rows)

cat('[3/4] Applying global SNP filters and computing Dxy\n')
site_missing <- rowMeans(is.na(G))
p <- rowMeans(G, na.rm=TRUE) / 2
maf <- pmin(p, 1 - p)
mask <- (site_missing <= 0.20) & is.finite(maf) & (maf >= 0.05)
G <- G[mask, , drop=FALSE]

pidx <- lapply(pop_order, function(pop) which(pop_lookup[samples] == pop))
names(pidx) <- pop_order

n_pop <- length(pop_order)
dxy_mat <- matrix(NA_real_, nrow=n_pop, ncol=n_pop, dimnames=list(pop_order, pop_order))

for (i in seq_len(n_pop)) {
  for (j in seq_len(n_pop)) {
    if (i == j) next
    if (j < i) {
      dxy_mat[i, j] <- dxy_mat[j, i]
      next
    }
    g1 <- G[, pidx[[i]], drop=FALSE]
    g2 <- G[, pidx[[j]], drop=FALSE]
    p1 <- rowMeans(g1, na.rm=TRUE) / 2
    p2 <- rowMeans(g2, na.rm=TRUE) / 2
    valid <- is.finite(p1) & is.finite(p2)
    if (!any(valid)) {
      dxy <- NA_real_
    } else {
      dxy_site <- p1[valid] * (1 - p2[valid]) + p2[valid] * (1 - p1[valid])
      dxy <- mean(dxy_site, na.rm=TRUE)
    }
    dxy_mat[i, j] <- dxy
    dxy_mat[j, i] <- dxy
  }
}
mean_dxy <- rowMeans(dxy_mat, na.rm=TRUE)

fst <- read.table(fst_file, sep='\t', header=TRUE, check.names=FALSE, row.names=1)
fst <- fst[pop_order, pop_order]
diag(fst) <- NA_real_
mean_fst <- rowMeans(as.matrix(fst), na.rm=TRUE)

cat('[4/4] Building publication-ready table\n')
t2 <- read.table(table2_file, sep='\t', header=TRUE, check.names=FALSE)
dem <- read.table(demog_file, sep='\t', header=TRUE, check.names=FALSE)
fis <- read.table(fis_file, sep='\t', header=TRUE, check.names=FALSE)

m <- data.frame(Population=pop_order)
m <- merge(m, t2[, c('Population','N_samples','Polymorphic_SNP_sites','Private','Mean_MAF_polymorphic_sites')], by='Population', all.x=TRUE)
m <- merge(m, dem[, c('Population','Segregating_sites','Mean_MAF','TajimasD_ascertained','Pi_per_site')], by='Population', all.x=TRUE)
m <- merge(m, fis[, c('Population','Fis')], by='Population', all.x=TRUE)

m$SNPs_total <- m$Segregating_sites
m$SNPs_biallelic <- m$Polymorphic_SNP_sites
m$Allelic_diversity_total <- 2 * m$Mean_MAF * (1 - m$Mean_MAF)
m$Allelic_diversity_biallelic <- 2 * m$Mean_MAF_polymorphic_sites * (1 - m$Mean_MAF_polymorphic_sites)
m$Dxy_mean_to_other_populations <- mean_dxy[m$Population]
m$Fst_mean_to_other_populations <- mean_fst[m$Population]

final <- m[, c(
  'Population','N_samples','SNPs_total','SNPs_biallelic',
  'Allelic_diversity_total','Allelic_diversity_biallelic',
  'Private','TajimasD_ascertained','Dxy_mean_to_other_populations',
  'Fst_mean_to_other_populations','Fis','Pi_per_site'
)]

colnames(final) <- c(
  'Population','Sample_size','SNPs_total','SNPs_biallelic',
  'Allelic_diversity_total','Allelic_diversity_biallelic',
  'Private_alleles','Tajimas_D','Dxy_mean_to_other_populations',
  'Fst_mean_to_other_populations','Fis','Pi_nucleotide_diversity'
)

num_cols <- c('Allelic_diversity_total','Allelic_diversity_biallelic','Tajimas_D','Dxy_mean_to_other_populations','Fst_mean_to_other_populations','Fis','Pi_nucleotide_diversity')
for (cc in num_cols) final[[cc]] <- round(final[[cc]], 6)

write.table(final, out_tsv, sep='\t', quote=FALSE, row.names=FALSE)
write.table(as.data.frame(dxy_mat), out_dxy, sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)

md_header <- c(
  '# Table 3. SNP diversity summary by population',
  '',
  "| Population | Sample size | SNPs (total) | SNPs (biallelic) | Allelic diversity (total) | Allelic diversity (biallelic) | Private alleles | Tajima's D | Dxy (mean to other pops) | Fst (mean to other pops) | Fis | Pi (nucleotide diversity) |",
  '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|'
)
md_rows <- apply(final, 1, function(r) {
  sprintf('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |',
          r['Population'], r['Sample_size'], r['SNPs_total'], r['SNPs_biallelic'],
          r['Allelic_diversity_total'], r['Allelic_diversity_biallelic'],
          r['Private_alleles'], r['Tajimas_D'], r['Dxy_mean_to_other_populations'],
          r['Fst_mean_to_other_populations'], r['Fis'], r['Pi_nucleotide_diversity'])
})
md_notes <- c(
  '',
  'Notes:',
  '- `SNPs (total)` are segregating sites from the demographic summary table.',
  '- `SNPs (biallelic)` are polymorphic biallelic SNP sites in the filtered panel.',
  '- `Allelic diversity` is expected heterozygosity (2p(1-p)) from the corresponding mean MAF.',
  '- `Dxy` and `Fst` are mean pairwise values for each focal population vs the other five populations.'
)
writeLines(c(md_header, md_rows, md_notes), con=out_md)

cat('Wrote:', out_tsv, '\n')
cat('Wrote:', out_md, '\n')
cat('Wrote:', out_dxy, '\n')
