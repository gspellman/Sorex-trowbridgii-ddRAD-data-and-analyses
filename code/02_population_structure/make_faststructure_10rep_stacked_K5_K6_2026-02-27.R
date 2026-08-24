base <- '/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23/faststructure_analysis_2026-02-27'
run_dir <- file.path(base, 'work', 'runs')
out_figs <- file.path(base, 'figures')
out_tables <- file.path(base, 'tables')

dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

res <- read.table(file.path(out_tables, 'faststructure_10rep_per_run_metrics.tsv'), sep='\t', header=TRUE, stringsAsFactors=FALSE)

pick_best_rep <- function(k) {
  x <- res[res$K == k, ]
  x <- x[order(-x$ml), ]
  x$rep[1]
}

rep5 <- pick_best_rep(5)
rep6 <- pick_best_rep(6)

fam <- read.table(file.path(run_dir, 'unlinked.fam'), header=FALSE, stringsAsFactors=FALSE)
colnames(fam) <- c('FID','IID','PID','MID','SEX','PHENO')
pop <- read.table('/Users/gspellman/Trowbridgii_analyses/Trowbridgii_popfile.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
ord <- merge(data.frame(IID=fam$IID, idx=seq_along(fam$IID), stringsAsFactors=FALSE), pop, by.x='IID', by.y='Sample', all.x=TRUE, sort=FALSE)
ord$Pop <- factor(ord$Pop, levels=c('North','North_Coast','Sierra_1','Sierra_2','Sierra_3','South_Coast'))
ord <- ord[order(ord$Pop, ord$IID), ]

split_idx <- split(seq_len(nrow(ord)), ord$Pop)
split_idx <- split_idx[names(split_idx) != 'NA']
centers <- sapply(split_idx, function(v) mean(range(v)))
labels <- names(centers)
seps <- cumsum(sapply(split_idx, length)) + 0.5
seps <- seps[seps < nrow(ord)]

getQ <- function(rep_id, k) {
  qfile <- file.path(run_dir, sprintf('fs10_rep%d.%d.meanQ', rep_id, k))
  Q <- as.matrix(read.table(qfile, header=FALSE))
  Q <- Q[ord$idx, , drop=FALSE]
  Q
}

Q5 <- getQ(rep5, 5)
Q6 <- getQ(rep6, 6)

pal5 <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e')
pal6 <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')

png(file.path(out_figs, 'faststructure_10rep_stacked_compoplot_K5_K6.png'), width=3600, height=2200, res=300)
layout(matrix(c(1,2), nrow=2), heights=c(1,1))
par(oma=c(0.2,0.2,0.2,0.2))

par(mar=c(1.2,4.5,2.4,1), xaxs='i', yaxs='i')
barplot(t(Q5), col=pal5, border=NA, space=0, axes=FALSE)
axis(2, las=1)
abline(v=seps, col='white', lwd=1)
box()
mtext(sprintf('K=5 (best rep=%d)', rep5), side=3, line=1, font=2)
mtext('Ancestry proportion', side=2, line=3)

par(mar=c(3.8,4.5,2.4,1), xaxs='i', yaxs='i')
barplot(t(Q6), col=pal6, border=NA, space=0, axes=FALSE)
axis(2, las=1)
axis(1, at=centers, labels=labels, cex.axis=0.85, line=0.2)
abline(v=seps, col='white', lwd=1)
box()
mtext(sprintf('K=6 (best rep=%d)', rep6), side=3, line=1, font=2)
mtext('Ancestry proportion', side=2, line=3)
mtext('Population', side=1, line=2.4)

dev.off()

pdf(file.path(out_figs, 'faststructure_10rep_stacked_compoplot_K5_K6.pdf'), width=12, height=7.2)
layout(matrix(c(1,2), nrow=2), heights=c(1,1))
par(oma=c(0.2,0.2,0.2,0.2))

par(mar=c(1.2,4.5,2.4,1), xaxs='i', yaxs='i')
barplot(t(Q5), col=pal5, border=NA, space=0, axes=FALSE)
axis(2, las=1)
abline(v=seps, col='white', lwd=1)
box()
mtext(sprintf('K=5 (best rep=%d)', rep5), side=3, line=1, font=2)
mtext('Ancestry proportion', side=2, line=3)

par(mar=c(3.8,4.5,2.4,1), xaxs='i', yaxs='i')
barplot(t(Q6), col=pal6, border=NA, space=0, axes=FALSE)
axis(2, las=1)
axis(1, at=centers, labels=labels, cex.axis=0.85, line=0.2)
abline(v=seps, col='white', lwd=1)
box()
mtext(sprintf('K=6 (best rep=%d)', rep6), side=3, line=1, font=2)
mtext('Ancestry proportion', side=2, line=3)
mtext('Population', side=1, line=2.4)

dev.off()

write.table(data.frame(K=c(5,6), best_rep=c(rep5,rep6)),
            file.path(out_tables, 'faststructure_10rep_stacked_K5_K6_best_reps.tsv'),
            sep='\t', row.names=FALSE, quote=FALSE)

cat('Done\n')
