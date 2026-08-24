base <- '/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23/faststructure_analysis_2026-02-27'
run_dir <- file.path(base, 'work', 'runs')
out_tables <- file.path(base, 'tables')
out_figs <- file.path(base, 'figures')

extract_metrics <- function(logfile) {
  x <- readLines(logfile, warn = FALSE)
  ml_line <- x[grepl('Marginal Likelihood\\s*=', x)]
  cv_line <- x[grepl('CV error\\s*=', x)]
  ml <- if (length(ml_line)) as.numeric(sub('.*=\\s*', '', tail(ml_line, 1))) else NA
  cv <- NA; cv_se <- NA
  if (length(cv_line)) {
    nums <- regmatches(tail(cv_line,1), gregexpr('[-]?[0-9]+\\.?[0-9]*', tail(cv_line,1)))[[1]]
    nums <- as.numeric(nums)
    if (length(nums) >= 1) cv <- nums[1]
    if (length(nums) >= 2) cv_se <- nums[2]
  }
  c(ml=ml, cv=cv, cv_se=cv_se)
}

man <- read.table(file.path(out_tables, 'faststructure_run_manifest_10rep.tsv'), sep='\t', header=TRUE, stringsAsFactors=FALSE)
man <- man[man$state == 'ok', ]

rows <- list(); idx <- 1
for (i in seq_len(nrow(man))) {
  r <- man$rep[i]; k <- man$K[i]
  lf <- file.path(run_dir, sprintf('fs10_rep%d.%d.log', r, k))
  met <- extract_metrics(lf)
  rows[[idx]] <- data.frame(rep=r, K=k, seed=man$seed[i], logfile=lf, ml=met['ml'], cv=met['cv'], cv_se=met['cv_se'])
  idx <- idx + 1
}
res <- do.call(rbind, rows)
res <- res[order(res$rep, res$K), ]
write.table(res, file.path(out_tables, 'faststructure_10rep_per_run_metrics.tsv'), sep='\t', row.names=FALSE, quote=FALSE)

Ks <- sort(unique(res$K))
sumK <- data.frame(
  K = Ks,
  n_runs = sapply(Ks, function(k) sum(res$K == k)),
  ml_mean = sapply(Ks, function(k) mean(res$ml[res$K==k], na.rm=TRUE)),
  ml_sd = sapply(Ks, function(k) sd(res$ml[res$K==k], na.rm=TRUE)),
  cv_mean = sapply(Ks, function(k) mean(res$cv[res$K==k], na.rm=TRUE)),
  cv_sd = sapply(Ks, function(k) sd(res$cv[res$K==k], na.rm=TRUE))
)
write.table(sumK, file.path(out_tables, 'faststructure_10rep_k_summary.tsv'), sep='\t', row.names=FALSE, quote=FALSE)

parse_choose <- function(f) {
  x <- readLines(f, warn = FALSE)
  ml <- as.integer(sub('.*=\\s*', '', x[grepl('maximizes marginal likelihood', x)][1]))
  comp <- as.integer(sub('.*=\\s*', '', x[grepl('used to explain structure', x)][1]))
  rp <- as.integer(sub('.*_rep([0-9]+)\\.txt$', '\\1', basename(f)))
  data.frame(rep=rp, chooseK_ml=ml, chooseK_components=comp)
}
choose_files <- Sys.glob(file.path(out_tables, 'chooseK_10rep_rep*.txt'))
choose <- do.call(rbind, lapply(choose_files, parse_choose))
choose <- choose[order(choose$rep), ]
write.table(choose, file.path(out_tables, 'faststructure_10rep_chooseK.tsv'), sep='\t', row.names=FALSE, quote=FALSE)

mode_k <- function(v) {
  tb <- sort(table(v), decreasing=TRUE)
  as.integer(names(tb)[1])
}
best_k_choose <- mode_k(choose$chooseK_components)
best_k_ml <- sumK$K[which.max(sumK$ml_mean)]
best_k_cv <- sumK$K[which.min(sumK$cv_mean)]
best_k <- best_k_choose
best_rows <- res[res$K == best_k, ]
best_row <- best_rows[order(-best_rows$ml), ][1, ]

sel <- data.frame(
  criterion = c('chooseK_mode_components','mean_marginal_likelihood','mean_cv_error','selected_bestK','selected_best_rep'),
  value = c(best_k_choose,best_k_ml,best_k_cv,best_k,best_row$rep)
)
write.table(sel, file.path(out_tables, 'faststructure_10rep_best_k_selection.tsv'), sep='\t', row.names=FALSE, quote=FALSE)

# compoplot
qfile <- file.path(run_dir, sprintf('fs10_rep%d.%d.meanQ', best_row$rep, best_k))
Q <- as.matrix(read.table(qfile, header=FALSE))
fam <- read.table(file.path(run_dir, 'unlinked.fam'), header=FALSE, stringsAsFactors=FALSE)
colnames(fam) <- c('FID','IID','PID','MID','SEX','PHENO')
pop <- read.table('/Users/gspellman/Trowbridgii_analyses/Trowbridgii_popfile.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
ord <- merge(data.frame(IID=fam$IID, idx=seq_along(fam$IID), stringsAsFactors=FALSE), pop, by.x='IID', by.y='Sample', all.x=TRUE, sort=FALSE)
ord$Pop <- factor(ord$Pop, levels=c('North','North_Coast','Sierra_1','Sierra_2','Sierra_3','South_Coast'))
ord <- ord[order(ord$Pop, ord$IID), ]
Q <- Q[ord$idx, , drop=FALSE]

split_idx <- split(seq_len(nrow(ord)), ord$Pop)
split_idx <- split_idx[names(split_idx) != 'NA']
centers <- sapply(split_idx, function(v) mean(range(v)))
labels <- names(centers)
seps <- cumsum(sapply(split_idx, length)) + 0.5
seps <- seps[seps < nrow(ord)]

k <- ncol(Q)
cluster_cols <- c(
  '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
  '#e6ab02', '#a6761d', '#1f78b4', '#b15928', '#17becf'
)[seq_len(k)]
cluster_labs <- paste0('Cluster', seq_len(k))

png(file.path(out_figs, 'faststructure_10rep_bestK_compoplot_publication.png'), width=3600, height=1400, res=300)
layout(matrix(c(1,2), nrow=1), widths=c(4.9,1.1))
par(mar=c(6,4.5,3,0.5), xaxs='i', yaxs='i')
barplot(t(Q), col=cluster_cols, border=NA, space=0, axes=FALSE)
axis(2, las=1)
axis(1, at=centers, labels=labels, cex.axis=0.9)
abline(v=seps, col='white', lwd=1)
box()
mtext(sprintf('fastStructure compoplot (10 reps/K, best K=%d; rep=%d)', best_k, best_row$rep), side=3, line=1, font=2)
mtext('Ancestry proportion', side=2, line=3)
mtext('Population', side=1, line=4.5)
par(mar=c(6,0.5,3,1))
plot.new()
legend('center', legend=cluster_labs, fill=cluster_cols, bty='n', cex=1.0)
dev.off()

pdf(file.path(out_figs, 'faststructure_10rep_bestK_compoplot_publication.pdf'), width=12, height=4.2)
layout(matrix(c(1,2), nrow=1), widths=c(4.9,1.1))
par(mar=c(6,4.5,3,0.5), xaxs='i', yaxs='i')
barplot(t(Q), col=cluster_cols, border=NA, space=0, axes=FALSE)
axis(2, las=1)
axis(1, at=centers, labels=labels, cex.axis=0.9)
abline(v=seps, col='white', lwd=1)
box()
mtext(sprintf('fastStructure compoplot (10 reps/K, best K=%d; rep=%d)', best_k, best_row$rep), side=3, line=1, font=2)
mtext('Ancestry proportion', side=2, line=3)
mtext('Population', side=1, line=4.5)
par(mar=c(6,0.5,3,1))
plot.new()
legend('center', legend=cluster_labs, fill=cluster_cols, bty='n', cex=0.9)
dev.off()

png(file.path(out_figs, 'faststructure_10rep_K_diagnostics.png'), width=2200, height=2200, res=300)
par(mfrow=c(2,1), mar=c(4,4.5,2,1))
plot(sumK$K, sumK$ml_mean, type='b', pch=19, xlab='K', ylab='Mean marginal likelihood', main='fastStructure K diagnostics (10 reps/K)')
plot(sumK$K, sumK$cv_mean, type='b', pch=19, xlab='K', ylab='Mean CV error')
dev.off()

pdf(file.path(out_figs, 'faststructure_10rep_K_diagnostics.pdf'), width=6.5, height=6.5)
par(mfrow=c(2,1), mar=c(4,4.5,2,1))
plot(sumK$K, sumK$ml_mean, type='b', pch=19, xlab='K', ylab='Mean marginal likelihood', main='fastStructure K diagnostics (10 reps/K)')
plot(sumK$K, sumK$cv_mean, type='b', pch=19, xlab='K', ylab='Mean CV error')
dev.off()

cat('Done\n')
