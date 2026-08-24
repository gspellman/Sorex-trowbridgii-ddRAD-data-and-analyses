#!/usr/bin/env python3
import os
import itertools
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
OUT = os.path.join(BASE, "downstream_analysis_2026-02-23", "additional_analyses_2026-02-24")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
for d in (FIG, TAB):
    os.makedirs(d, exist_ok=True)

VCF = os.path.join(BASE, "populations.snps.vcf")
POPMAP = os.path.join(BASE, "popmap.tsv")

pop_cols = {
    "North": "#4169E1",
    "North_Coast": "#4CBB17",
    "Sierra_1": "#800080",
    "Sierra_2": "#00FFFF",
    "Sierra_3": "#FF0000",
    "South_Coast": "#000000",
}
pop_order = ["North", "North_Coast", "Sierra_1", "Sierra_2", "Sierra_3", "South_Coast"]


def parse_gt(gt):
    g = gt.split(':', 1)[0].replace('|', '/')
    if g in ("./.", "."):
        return np.nan
    if g == "0/0":
        return 0.0
    if g in ("0/1", "1/0"):
        return 1.0
    if g == "1/1":
        return 2.0
    return np.nan


def read_vcf(vcf, sample_set):
    chrom, pos, vid, rows = [], [], [], []
    samples = None
    idx = None
    with open(vcf, 'r') as f:
        for ln in f:
            if ln.startswith('##'):
                continue
            if ln.startswith('#CHROM'):
                parts = ln.rstrip('\n').split('\t')
                all_s = parts[9:]
                samples = [s for s in all_s if s in sample_set]
                idx = [all_s.index(s) for s in samples]
                continue
            p = ln.rstrip('\n').split('\t')
            ref, alt = p[3], p[4]
            if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                continue
            g = p[9:]
            row = [parse_gt(g[j]) for j in idx]
            chrom.append(p[0])
            pos.append(int(p[1]))
            vid.append(p[2] if p[2] != '.' else f"{p[0]}:{p[1]}")
            rows.append(row)
    return samples, np.array(chrom), np.array(pos), np.array(vid), np.array(rows, dtype=float)


def hudson_fst(a1, n1, a2, n2):
    # a1/a2 = alt allele counts, n1/n2 = called chromosomes
    if n1 <= 1 or n2 <= 1:
        return np.nan
    p1 = a1 / n1
    p2 = a2 / n2
    num = (p1 - p2) ** 2 - (p1 * (1 - p1) / (n1 - 1)) - (p2 * (1 - p2) / (n2 - 1))
    den = p1 * (1 - p2) + p2 * (1 - p1)
    if den <= 0:
        return np.nan
    return num / den


# load metadata
popmap = pd.read_csv(POPMAP, sep='\t', header=None, names=['Sample', 'Population'])
popmap = popmap[popmap['Population'].isin(pop_order)].copy()
sample_set = set(popmap['Sample'])

samples, chrom, pos, vid, G = read_vcf(VCF, sample_set)
popmap = popmap.set_index('Sample').loc[samples].reset_index()
spop = dict(zip(popmap['Sample'], popmap['Population']))

# base filters
site_missing = np.mean(np.isnan(G), axis=1)
maf = np.minimum(np.nanmean(G, axis=1) / 2, 1 - np.nanmean(G, axis=1) / 2)
mask = (site_missing <= 0.20) & np.isfinite(maf) & (maf >= 0.05)
G = G[mask]
chrom = chrom[mask]
pos = pos[mask]
vid = vid[mask]

# per-pop indices
pidx = {p: np.array([i for i, s in enumerate(samples) if spop[s] == p], dtype=int) for p in pop_order}

# coordinate offsets for concatenated Manhattan x-axis
chrom_order = pd.unique(chrom)
offset = {}
cum = 0
for c in chrom_order:
    m = int(pos[chrom == c].max())
    offset[c] = cum
    cum += m + 1
pos_cum = np.array([p + offset[c] for c, p in zip(chrom, pos)], dtype=float)

pair_rows = []
summary_rows = []
window_rows = []

pairs = list(itertools.combinations(pop_order, 2))

for p1, p2 in pairs:
    i1, i2 = pidx[p1], pidx[p2]
    g1, g2 = G[:, i1], G[:, i2]

    n1 = 2 * np.sum(~np.isnan(g1), axis=1)
    n2 = 2 * np.sum(~np.isnan(g2), axis=1)
    a1 = np.nansum(g1, axis=1)
    a2 = np.nansum(g2, axis=1)

    fst = np.array([hudson_fst(x1, y1, x2, y2) for x1, y1, x2, y2 in zip(a1, n1, a2, n2)], dtype=float)

    ok = np.isfinite(fst)
    fst_use = fst[ok]
    if fst_use.size == 0:
        continue

    thr = np.quantile(fst_use, 0.99)
    outlier = ok & (fst >= thr)

    # per-SNP table
    df = pd.DataFrame({
        'Comparison': f"{p1}_vs_{p2}",
        'Pop1': p1,
        'Pop2': p2,
        'CHROM': chrom,
        'POS': pos,
        'ID': vid,
        'Fst_Hudson': fst,
        'Outlier_top1pct': outlier.astype(int)
    })
    pair_rows.append(df)

    # fixed windows for overlap comparison
    win_size = 50000
    wstart = ((pos - 1) // win_size) * win_size + 1
    wend = wstart + win_size - 1
    wdf = pd.DataFrame({
        'Comparison': f"{p1}_vs_{p2}",
        'Pop1': p1,
        'Pop2': p2,
        'CHROM': chrom,
        'WIN_START': wstart,
        'WIN_END': wend,
        'Fst': fst,
        'Outlier': outlier
    })
    wsum = wdf.groupby(['Comparison', 'Pop1', 'Pop2', 'CHROM', 'WIN_START', 'WIN_END'], as_index=False).agg(
        Mean_Fst=('Fst', 'mean'),
        N_SNP=('Fst', lambda x: np.sum(np.isfinite(x))),
        Any_outlier=('Outlier', 'max')
    )
    # classify divergent windows by top 1% of mean window Fst for this comparison
    wf = wsum[np.isfinite(wsum['Mean_Fst']) & (wsum['N_SNP'] >= 5)].copy()
    if len(wf):
        wthr = np.quantile(wf['Mean_Fst'], 0.99)
        wf['Divergent_window'] = (wf['Mean_Fst'] >= wthr).astype(int)
        window_rows.append(wf)

    summary_rows.append({
        'Comparison': f"{p1}_vs_{p2}",
        'Pop1': p1,
        'Pop2': p2,
        'N_SNP_used': int(np.sum(ok)),
        'Fst_mean': float(np.nanmean(fst_use)),
        'Fst_median': float(np.nanmedian(fst_use)),
        'Fst_p99_threshold': float(thr),
        'N_outlier_SNP': int(np.sum(outlier))
    })

# write core tables
all_pairs = pd.concat(pair_rows, ignore_index=True)
all_pairs.to_csv(os.path.join(TAB, 'pairwise_fst_by_snp_all_comparisons.tsv'), sep='\t', index=False)
pd.DataFrame(summary_rows).to_csv(os.path.join(TAB, 'pairwise_fst_comparison_summary.tsv'), sep='\t', index=False)

# Manhattan facet figure
n = len(summary_rows)
cols = 3
rows = int(np.ceil(n / cols))
fig, axes = plt.subplots(rows, cols, figsize=(cols * 5.5, rows * 3.8), squeeze=False)
axes = axes.flatten()

for k, s in enumerate(summary_rows):
    comp = s['Comparison']
    p1 = s['Pop1']
    p2 = s['Pop2']
    d = all_pairs[all_pairs['Comparison'] == comp].copy()
    x = np.array([p + offset[c] for c, p in zip(d['CHROM'], d['POS'])], dtype=float)
    y = d['Fst_Hudson'].values
    o = d['Outlier_top1pct'].values.astype(bool)

    ax = axes[k]
    # alternating scaffold colors
    ch_idx = {c: i for i, c in enumerate(chrom_order)}
    colv = np.array(['#4c78a8' if ch_idx[c] % 2 == 0 else '#9ecae9' for c in d['CHROM']])
    ax.scatter(x[~o], y[~o], c=colv[~o], s=4, alpha=0.5, linewidths=0)
    ax.scatter(x[o], y[o], c='#d62728', s=6, alpha=0.9, linewidths=0)
    ax.axhline(s['Fst_p99_threshold'], ls='--', lw=0.8, color='black')
    ax.set_title(f"{p1} vs {p2}", fontsize=9)
    ax.set_ylim(bottom=min(-0.05, np.nanmin(y) - 0.02))
    if k % cols == 0:
        ax.set_ylabel('Per-SNP Hudson FST')
    else:
        ax.set_ylabel('')
    ax.set_xlabel('Genome coordinate (concatenated scaffolds)')

for j in range(k + 1, len(axes)):
    axes[j].axis('off')

fig.suptitle('Pairwise population FST Manhattan scans', fontsize=14, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(os.path.join(FIG, 'Pairwise_Fst_manhattan_facets_publication.png'), dpi=320)
fig.savefig(os.path.join(FIG, 'Pairwise_Fst_manhattan_facets_publication.pdf'))
plt.close(fig)

# overlap analysis of divergent windows across comparisons
w_all = pd.concat(window_rows, ignore_index=True)
w_all.to_csv(os.path.join(TAB, 'pairwise_fst_window_stats_all_comparisons.tsv'), sep='\t', index=False)

# universe of windows with enough SNPs for each comparison
comps = sorted(w_all['Comparison'].unique().tolist())
win_key = w_all['CHROM'].astype(str) + ':' + w_all['WIN_START'].astype(str) + '-' + w_all['WIN_END'].astype(str)
w_all = w_all.assign(WindowKey=win_key)

# sets of divergent windows per comparison
div_sets = {}
univ_sets = {}
for c in comps:
    d = w_all[w_all['Comparison'] == c]
    univ_sets[c] = set(d['WindowKey'])
    div_sets[c] = set(d.loc[d['Divergent_window'] == 1, 'WindowKey'])

# shared divergence table
ov_rows = []
for c1, c2 in itertools.combinations(comps, 2):
    inter = div_sets[c1].intersection(div_sets[c2])
    union = div_sets[c1].union(div_sets[c2])
    j = len(inter) / len(union) if len(union) else np.nan

    # fisher exact on overlapping universe
    u = univ_sets[c1].intersection(univ_sets[c2])
    a = len((div_sets[c1] & div_sets[c2]) & u)
    b = len((div_sets[c1] - div_sets[c2]) & u)
    c = len((div_sets[c2] - div_sets[c1]) & u)
    d = len((u - div_sets[c1] - div_sets[c2]))
    table = np.array([[a, b], [c, d]], dtype=int)
    _, p = stats.fisher_exact(table, alternative='greater') if np.all(table >= 0) else (np.nan, np.nan)

    ov_rows.append({
        'Comparison1': c1,
        'Comparison2': c2,
        'N_shared_divergent_windows': a,
        'Jaccard_shared_divergence': j,
        'FisherP_enrichment': p,
        'Universe_windows_intersection': len(u)
    })

ov = pd.DataFrame(ov_rows).sort_values('Jaccard_shared_divergence', ascending=False)
# BH FDR
if len(ov):
    pvals = ov['FisherP_enrichment'].values
    order = np.argsort(np.where(np.isfinite(pvals), pvals, 1.0))
    ranked = np.empty_like(order)
    ranked[order] = np.arange(1, len(order) + 1)
    ov['FDR_BH'] = np.minimum(1.0, ov['FisherP_enrichment'] * len(ov) / ranked)
ov.to_csv(os.path.join(TAB, 'divergent_window_overlap_between_pairwise_comparisons.tsv'), sep='\t', index=False)

# matrix heatmap figure
mat = pd.DataFrame(np.nan, index=comps, columns=comps)
for c in comps:
    mat.loc[c, c] = 1.0
for _, r in ov.iterrows():
    mat.loc[r['Comparison1'], r['Comparison2']] = r['Jaccard_shared_divergence']
    mat.loc[r['Comparison2'], r['Comparison1']] = r['Jaccard_shared_divergence']

for ext in ['png', 'pdf']:
    plt.figure(figsize=(10, 8), dpi=320 if ext == 'png' else None)
    arr = mat.values.astype(float)
    im = plt.imshow(arr, cmap='magma', vmin=0, vmax=np.nanmax(arr[np.isfinite(arr)]) if np.any(np.isfinite(arr)) else 1)
    plt.xticks(range(len(comps)), comps, rotation=70, ha='right', fontsize=7)
    plt.yticks(range(len(comps)), comps, fontsize=7)
    plt.title('Shared divergent genomic windows across pairwise FST comparisons\n(Jaccard overlap)')
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.set_label('Jaccard overlap of divergent windows')
    # annotate upper triangle
    for i in range(len(comps)):
        for j in range(len(comps)):
            v = arr[i, j]
            if np.isfinite(v):
                plt.text(j, i, f"{v:.2f}", ha='center', va='center', fontsize=6, color='white' if v > 0.4 else 'black')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, f'Divergent_window_overlap_heatmap_publication.{ext}'))
    plt.close()

# also provide top shared windows table across >=2 comparisons
div_win = w_all[w_all['Divergent_window'] == 1].copy()
agg = div_win.groupby('WindowKey').agg(
    CHROM=('CHROM', 'first'),
    WIN_START=('WIN_START', 'first'),
    WIN_END=('WIN_END', 'first'),
    N_comparisons_divergent=('Comparison', 'nunique'),
    Comparisons=('Comparison', lambda x: ';'.join(sorted(set(x))))
).reset_index()
agg = agg.sort_values(['N_comparisons_divergent', 'CHROM', 'WIN_START'], ascending=[False, True, True])
agg.to_csv(os.path.join(TAB, 'shared_divergent_windows_ranked.tsv'), sep='\t', index=False)

with open(os.path.join(TAB, 'pairwise_fst_manhattan_and_overlap_summary.txt'), 'w') as out:
    out.write('Pairwise FST Manhattan and shared-divergence analysis summary\n')
    out.write(f'Comparisons evaluated: {len(summary_rows)}\n')
    out.write('Per-SNP FST estimator: Hudson FST\n')
    out.write('Divergent SNP threshold: top 1% FST within each comparison\n')
    out.write('Divergent window definition: top 1% mean-window FST (50 kb windows; >=5 SNP/window)\n')
    out.write('Overlap metrics: pairwise Jaccard and Fisher exact enrichment\n')

print('Done', flush=True)
