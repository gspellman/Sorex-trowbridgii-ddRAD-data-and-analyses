#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(grid)
})

base_dir <- '/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23'
tab_file <- file.path(base_dir, 'tables', 'Table3_population_SNP_diversity_publication.tsv')
table1_file <- file.path(base_dir, 'tables', 'Table1_STACKS_run_summary.tsv')
table2_file <- file.path(base_dir, 'tables', 'Table2_population_SNP_variation_summary.tsv')
fig_dir <- file.path(base_dir, 'figures')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(fig_dir, 'Table3_population_SNP_diversity_publication.pdf')
out_png <- file.path(fig_dir, 'Table3_population_SNP_diversity_publication.png')

df <- read.table(tab_file, sep='\t', header=TRUE, check.names=FALSE)
t1 <- read.table(table1_file, sep='\t', header=TRUE, check.names=FALSE)
t2 <- read.table(table2_file, sep='\t', header=TRUE, check.names=FALSE)

total_snps <- t1$Value[t1$Metric == 'Variant sites retained'][1]
biallelic_snps <- max(t2$Total_biallelic_SNP_sites, na.rm = TRUE)
title_text <- sprintf(
  'Table 3. SNP diversity summary by population (total SNPs=%s; biallelic SNPs=%s)',
  format(as.integer(total_snps), big.mark = ','),
  format(as.integer(biallelic_snps), big.mark = ',')
)

fmt_num <- function(x, digits = 6) sprintf(paste0('%.', digits, 'f'), as.numeric(x))

disp <- data.frame(
  Population = df$Population,
  `Sample size` = df$Sample_size,
  `Allelic diversity\n(total)` = fmt_num(df$Allelic_diversity_total, 6),
  `Allelic diversity\n(biallelic)` = fmt_num(df$Allelic_diversity_biallelic, 6),
  `Private\nalleles` = format(df$Private_alleles, big.mark=','),
  `Tajima\'s D` = fmt_num(df$Tajimas_D, 6),
  `Dxy\n(mean to others)` = fmt_num(df$Dxy_mean_to_other_populations, 6),
  `Fst\n(mean to others)` = fmt_num(df$Fst_mean_to_other_populations, 6),
  `Fis` = fmt_num(df$Fis, 6),
  `Pi\n(nucleotide diversity)` = fmt_num(df$Pi_nucleotide_diversity, 6),
  check.names = FALSE
)

headers <- colnames(disp)
mat <- as.matrix(disp)
n_row <- nrow(mat)
n_col <- ncol(mat)

# Wider columns for long headers
col_w <- c(1.75, 1.05, 1.30, 1.36, 1.00, 1.00, 1.15, 1.10, 0.95, 1.40)
col_w <- col_w / sum(col_w)
col_left <- c(0, cumsum(col_w))[1:n_col]

row_h <- 1 / (n_row + 2.6)  # includes title and footnote space

draw_table <- function() {
  grid.newpage()

  # Title
  grid.text(title_text,
            x = unit(0.5, 'npc'), y = unit(0.975, 'npc'),
            gp = gpar(fontsize = 12.2, fontface = 'bold'))

  top <- 0.92

  # Header background
  grid.rect(x = unit(0.5, 'npc'), y = unit(top - row_h/2, 'npc'),
            width = unit(0.98, 'npc'), height = unit(row_h, 'npc'),
            gp = gpar(fill = '#E6E6E6', col = NA))

  # Zebra rows
  for (i in seq_len(n_row)) {
    y <- top - row_h * (i + 0.5)
    fill <- if (i %% 2 == 0) '#F7F7F7' else 'white'
    grid.rect(x = unit(0.5, 'npc'), y = unit(y, 'npc'),
              width = unit(0.98, 'npc'), height = unit(row_h, 'npc'),
              gp = gpar(fill = fill, col = NA))
  }

  # Outer border
  table_height <- row_h * (n_row + 1)
  grid.rect(x = unit(0.5, 'npc'), y = unit(top - table_height/2, 'npc'),
            width = unit(0.98, 'npc'), height = unit(table_height, 'npc'),
            gp = gpar(fill = NA, col = '#404040', lwd = 1.2))

  # Vertical grid lines + header text + body text
  left_margin <- 0.01
  table_w <- 0.98

  for (j in seq_len(n_col)) {
    x0 <- left_margin + table_w * col_left[j]
    cw <- table_w * col_w[j]

    # Vertical separator
    if (j > 1) {
      grid.lines(x = unit(c(x0, x0), 'npc'),
                 y = unit(c(top - table_height, top), 'npc'),
                 gp = gpar(col = '#B0B0B0', lwd = 0.7))
    }

    # Header
    grid.text(headers[j],
              x = unit(x0 + cw/2, 'npc'),
              y = unit(top - row_h/2, 'npc'),
              gp = gpar(fontsize = 8.5, fontface = 'bold'))

    # Column alignment
    just <- if (j == 1) c('left', 'center') else c('right', 'center')
    xpos <- if (j == 1) x0 + 0.008 else x0 + cw - 0.008

    # Cells
    for (i in seq_len(n_row)) {
      y <- top - row_h * (i + 0.5)
      grid.text(mat[i, j],
                x = unit(xpos, 'npc'), y = unit(y, 'npc'),
                just = just,
                gp = gpar(fontsize = 8.6))
    }
  }

  # Horizontal line under header
  grid.lines(x = unit(c(0.01, 0.99), 'npc'),
             y = unit(c(top - row_h, top - row_h), 'npc'),
             gp = gpar(col = '#707070', lwd = 1.0))

  # Footnote
  foot <- 'Dxy and Fst are reported as mean pairwise values for each focal population versus the other five populations.'
  grid.text(foot, x = unit(0.01, 'npc'), y = unit(0.03, 'npc'),
            just = c('left', 'bottom'), gp = gpar(fontsize = 8.2, col = '#333333'))
}

pdf(out_pdf, width = 14, height = 4.8, useDingbats = FALSE)
draw_table()
dev.off()

png(out_png, width = 4200, height = 1450, res = 350)
draw_table()
dev.off()

cat('Wrote:', out_pdf, '\n')
cat('Wrote:', out_png, '\n')
